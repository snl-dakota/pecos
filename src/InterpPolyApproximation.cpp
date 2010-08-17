/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        InterpPolyApproximation
//- Description:  Implementation code for InterpPolyApproximation class
//-               
//- Owner:        Mike Eldred

#include "InterpPolyApproximation.hpp"
#include "TensorProductDriver.hpp"
#include "SparseGridDriver.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

//#define DEBUG


namespace Pecos {

int InterpPolyApproximation::min_coefficients() const
{
  // return the minimum number of data instances required to build the 
  // surface of multiple dimensions
  return (expansionCoeffFlag || expansionCoeffGradFlag) ? 1 : 0;
}


void InterpPolyApproximation::allocate_arrays()
{
  PolynomialApproximation::allocate_arrays();

  if (expansionCoeffFlag && expansionCoeffs.length() != numCollocPts)
    expansionCoeffs.sizeUninitialized(numCollocPts);
  if (expansionCoeffGradFlag) {
    const SurrogateDataPoint& sdp
      = (anchorPoint.is_null()) ? dataPoints[0] : anchorPoint;
    size_t num_deriv_vars = sdp.response_gradient().length();
    if (expansionCoeffGrads.numRows() != num_deriv_vars ||
	expansionCoeffGrads.numCols() != numCollocPts)
      expansionCoeffGrads.shapeUninitialized(num_deriv_vars, numCollocPts);
  }

  // checking numCollocPts is insufficient due to anisotropy --> changes in
  // anisotropic weights could move points around without changing the total.
  //bool update_exp_form =
  //  ( (expansionCoeffFlag && expansionCoeffs.length() != numCollocPts) ||
  //    (expansionCoeffGradFlag &&
  //     expansionCoeffGrads.numCols() != numCollocPts ) );

  switch (expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    const UShortArray&   quad_order = tpq_driver->quadrature_order();

    // verify total number of collocation pts (should not include anchorPoint)
    size_t i, j, num_gauss_pts = 1;
    for (i=0; i<numVars; ++i)
      num_gauss_pts *= quad_order[i];
    if (num_gauss_pts != numCollocPts) {
      PCerr << "Error: inconsistent total Gauss point count in "
	    << "InterpPolyApproximation::find_coefficients()" << std::endl;
      abort_handler(-1);
    }

    bool update_exp_form = (quad_order != quadOrderPrev);
    // *** TO DO: capture updates to parameterized/numerical polynomials
    //if (update_exp_form)
      // collocKey currently updated in TensorProductDriver::compute_grid()
      //tpq_driver->allocate_collocation_arrays();

    bool update_basis_form = update_exp_form;
    // *** TO DO: carefully evaluate interdependence between exp_form/basis_form
    if (update_basis_form) {
      // size and initialize polynomialBasis, one interpolant per variable
      if (polynomialBasis.empty())
	{ polynomialBasis.resize(1); polynomialBasis[0].resize(numVars); }
      const Real2DArray& gauss_pts_1d = driverRep->gauss_points_array()[0];
      std::vector<BasisPolynomial>& poly_basis_0 = polynomialBasis[0];
      for (i=0; i<numVars; ++i) {
	bool found = false;
	for (j=0; j<i; ++j)
	  if (gauss_pts_1d[i] == gauss_pts_1d[j]) // equality in pt vector
	    { found = true; break; }
	if (found) // reuse previous instance via shared representation
	  poly_basis_0[i] = poly_basis_0[j];
	else { // instantiate a new unique instance
	  poly_basis_0[i] = BasisPolynomial(LAGRANGE);
	  poly_basis_0[i].interpolation_points(gauss_pts_1d[i]);
	}
      }
    }

    quadOrderPrev = quad_order;
    break;
  }
  case SPARSE_GRID: {
    const RealVector& aniso_wts  = driverRep->anisotropic_weights();
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    unsigned short    ssg_level  = ssg_driver->level();

    bool update_exp_form
      = (ssg_level != ssgLevelPrev || aniso_wts != ssgAnisoWtsPrev);
    // *** TO DO: capture updates to parameterized/numerical polynomials
    if (update_exp_form) {
      ssg_driver->allocate_smolyak_arrays();
      ssg_driver->allocate_collocation_arrays();
    }

    bool update_basis_form = update_exp_form;
    // *** TO DO: carefully evaluate interdependence between exp_form/basis_form
    if (update_basis_form) {
      size_t i, j, k;
      unsigned short num_levels = ssg_level + 1;
      // size and initialize polynomialBasis, multiple interpolants per variable
      if (polynomialBasis.empty()) {
	polynomialBasis.resize(num_levels);
	for (i=0; i<num_levels; ++i)
	  polynomialBasis[i].resize(numVars);
      }
      // j range is 0:w inclusive; i range is 1:w+1 inclusive
      const Real3DArray& gauss_pts_1d = driverRep->gauss_points_array();
      for (i=0; i<num_levels; ++i) { // i->num_levels-1->ssg_level
	const Real2DArray&          gauss_pts_1d_i = gauss_pts_1d[i];
	std::vector<BasisPolynomial>& poly_basis_i = polynomialBasis[i];
	// anisotropic levels: check other dims at corresponding level
	for (j=0; j<numVars; ++j) {
	  bool found = false;
	  for (k=0; k<j; ++k)
	    if (//ssg_level[k]   >= i && // level i exists for dimension k
		gauss_pts_1d_i[j] == gauss_pts_1d_i[k]) // pt vector equality
	      { found = true; break; }
	  if (found) // reuse previous instances via shared representations
	    poly_basis_i[j] = poly_basis_i[k]; // shared rep
	  else { // instantiate new unique instances
	    poly_basis_i[j] = BasisPolynomial(LAGRANGE);
	    poly_basis_i[j].interpolation_points(gauss_pts_1d_i[j]);
	  }
	}
      }
    }

    ssgLevelPrev = ssg_level; ssgAnisoWtsPrev = aniso_wts;
    break;
  }
  }
}


void InterpPolyApproximation::find_coefficients()
{
  if (!expansionCoeffFlag && !expansionCoeffGradFlag) {
    PCerr << "Warning: neither expansion coefficients nor expansion "
	  << "coefficient gradients\n         are active in "
	  << "InterpPolyApproximation::find_coefficients().\n         "
	  << "Bypassing approximation construction." << std::endl;
    return;
  }

  // For testing of anchorPoint logic:
  //anchorPoint = dataPoints.front();
  //dataPoints.pop_front();

  size_t i, offset = 0;
  numCollocPts = dataPoints.size();
  // anchorPoint, if present, is the first expansionSample.
  if (!anchorPoint.is_null()) {
    offset        = 1;
    numCollocPts += offset;
  }
  if (!numCollocPts) {
    PCerr << "Error: nonzero number of sample points required in "
	  << "InterpPolyApproximation::find_coefficients()." << std::endl;
    abort_handler(-1);
  }

  allocate_arrays();

  if (!anchorPoint.is_null()) {
    if (expansionCoeffFlag)
      expansionCoeffs[0] = anchorPoint.response_function();
    if (expansionCoeffGradFlag)
      Teuchos::setCol(anchorPoint.response_gradient(), 0, expansionCoeffGrads);
  }

  std::vector<SurrogateDataPoint>::iterator it = dataPoints.begin();
  for (i=offset; i<numCollocPts; ++i, ++it) {
    if (expansionCoeffFlag)
      expansionCoeffs[i] = it->response_function();
    if (expansionCoeffGradFlag)
      Teuchos::setCol(it->response_gradient(), (int)i, expansionCoeffGrads);
  }
}


/** Overloaded version supporting tensor-product quadrature. */
const Real& InterpPolyApproximation::tensor_product_value(const RealVector& x)
{
  TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
  const UShort2DArray& key        = tpq_driver->collocation_key();

  tpValue = 0.;
  size_t i, j;
  std::vector<BasisPolynomial>& poly_basis_0 = polynomialBasis[0];
  for (i=0; i<numCollocPts; ++i) {
    Real L_i = 1.0;
    const UShortArray& key_i = key[i];
    for (j=0; j<numVars; ++j)
      L_i *= poly_basis_0[j].get_value(x[j], key_i[j]);
    tpValue += expansionCoeffs[i] * L_i;
  }
  return tpValue;
}


/** Overloaded version supporting Smolyak sparse grids. */
const Real& InterpPolyApproximation::
tensor_product_value(const RealVector& x, size_t tp_index)
{
  SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
  const UShort2DArray& key     = ssg_driver->collocation_key()[tp_index];
  const UShortArray& sm_index  = ssg_driver->smolyak_multi_index()[tp_index];
  const SizetArray& coeff_index
    = ssg_driver->expansion_coefficient_indices()[tp_index];

  tpValue = 0.;
  size_t i, j, num_colloc_pts = key.size();
  for (i=0; i<num_colloc_pts; ++i) {
    const UShortArray& key_i = key[i];
    Real L_i = 1.0;
    for (j=0; j<numVars; ++j)
      L_i *= polynomialBasis[sm_index[j]][j].get_value(x[j], key_i[j]);
    tpValue += expansionCoeffs[coeff_index[i]] * L_i;
  }
  return tpValue;
}


/** Overloaded version supporting tensor-product quadrature. */
const RealVector& InterpPolyApproximation::
tensor_product_gradient(const RealVector& x)
{
  TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
  const UShort2DArray& key        = tpq_driver->collocation_key();

  if (tpGradient.length() != numVars)
    tpGradient.sizeUninitialized(numVars);
  tpGradient = 0.;
  std::vector<BasisPolynomial>& poly_basis_0 = polynomialBasis[0];
  size_t i, j, k;
  for (i=0; i<numCollocPts; ++i) {
    const UShortArray& key_i   = key[i];
    const Real&        coeff_i = expansionCoeffs[i];
    for (j=0; j<numVars; ++j) {
      Real term_i_grad_j = 1.0;
      for (k=0; k<numVars; ++k)
	term_i_grad_j *= (k == j) ?
	  poly_basis_0[k].get_gradient(x[k], key_i[k]) :
	  poly_basis_0[k].get_value(x[k],    key_i[k]);
      tpGradient[j] += coeff_i * term_i_grad_j;
    }
  }
  return tpGradient;
}


/** Overloaded version supporting Smolyak sparse grids. */
const RealVector& InterpPolyApproximation::
tensor_product_gradient(const RealVector& x, size_t tp_index)
{
  SparseGridDriver*    ssg_driver = (SparseGridDriver*)driverRep;
  const UShort2DArray& key        = ssg_driver->collocation_key()[tp_index];
  const UShortArray&   sm_index   = ssg_driver->smolyak_multi_index()[tp_index];
  const SizetArray&    coeff_index
    = ssg_driver->expansion_coefficient_indices()[tp_index];

  if (tpGradient.length() != numVars)
    tpGradient.sizeUninitialized(numVars);
  tpGradient = 0.;
  size_t i, j, k, num_colloc_pts = key.size();
  for (i=0; i<num_colloc_pts; ++i) {
    const UShortArray& key_i   = key[i];
    const Real&        coeff_i = expansionCoeffs[coeff_index[i]];
    for (j=0; j<numVars; ++j) {
      Real term_i_grad_j = 1.0;
      for (k=0; k<numVars; ++k)
	term_i_grad_j *= (k == j) ?
	  polynomialBasis[sm_index[k]][k].get_gradient(x[k], key_i[k]) :
	  polynomialBasis[sm_index[k]][k].get_value(x[k],    key_i[k]);
      tpGradient[j] += coeff_i * term_i_grad_j;
    }
  }
  return tpGradient;
}


/** Overloaded version supporting tensor-product quadrature. */
const RealVector& InterpPolyApproximation::
tensor_product_gradient(const RealVector& x, const UIntArray& dvv)
{
  TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
  const UShort2DArray& key        = tpq_driver->collocation_key();

  size_t num_deriv_vars = dvv.size();
  if (tpGradient.length() != num_deriv_vars)
    tpGradient.sizeUninitialized(num_deriv_vars);
  tpGradient = 0.;
  std::vector<BasisPolynomial>& poly_basis_0 = polynomialBasis[0];
  size_t i, j, k, deriv_index;
  for (i=0; i<numCollocPts; ++i) {
    const UShortArray& key_i   = key[i];
    const Real&        coeff_i = expansionCoeffs[i];
    for (j=0; j<num_deriv_vars; ++j) {
      deriv_index = dvv[j] - 1; // requires an "All" view
      Real term_i_grad_j = 1.0;
      for (k=0; k<numVars; ++k)
	term_i_grad_j *= (k == deriv_index) ?
	  poly_basis_0[k].get_gradient(x[k], key_i[k]) :
	  poly_basis_0[k].get_value(x[k],    key_i[k]);
      tpGradient[j] += coeff_i * term_i_grad_j;
    }
  }
  return tpGradient;
}


/** Overloaded version supporting Smolyak sparse grids. */
const RealVector& InterpPolyApproximation::
tensor_product_gradient(const RealVector& x, size_t tp_index,
			const UIntArray& dvv)
{
  SparseGridDriver*    ssg_driver = (SparseGridDriver*)driverRep;
  const UShort2DArray& key        = ssg_driver->collocation_key()[tp_index];
  const UShortArray&   sm_index   = ssg_driver->smolyak_multi_index()[tp_index];
  const SizetArray&    coeff_index
    = ssg_driver->expansion_coefficient_indices()[tp_index];

  size_t num_deriv_vars = dvv.size();
  if (tpGradient.length() != num_deriv_vars)
    tpGradient.sizeUninitialized(num_deriv_vars);
  tpGradient = 0.;
  size_t i, j, k, deriv_index, num_colloc_pts = key.size();
  for (i=0; i<num_colloc_pts; ++i) {
    const UShortArray& key_i   = key[i];
    const Real&        coeff_i = expansionCoeffs[coeff_index[i]];
    for (j=0; j<num_deriv_vars; ++j) {
      deriv_index = dvv[j] - 1; // requires an "All" view
      Real term_i_grad_j = 1.0;
      for (k=0; k<numVars; ++k)
	term_i_grad_j *= (k == deriv_index) ?
	  polynomialBasis[sm_index[k]][k].get_gradient(x[k], key_i[k]) :
	  polynomialBasis[sm_index[k]][k].get_value(x[k],    key_i[k]);
      tpGradient[j] += coeff_i * term_i_grad_j;
    }
  }
  return tpGradient;
}


/** Overloaded version supporting tensor-product quadrature. */
const Real& InterpPolyApproximation::tensor_product_mean(const RealVector& x)
{
  TensorProductDriver* tpq_driver   = (TensorProductDriver*)driverRep;
  const UShort2DArray& key          = tpq_driver->collocation_key();
  const Real2DArray&   gauss_wts_1d = driverRep->gauss_weights_array()[0];

  tpMean = 0.;
  std::vector<BasisPolynomial>& poly_basis_0 = polynomialBasis[0];
  size_t i, j;
  for (i=0; i<numCollocPts; ++i) {
    const UShortArray& key_i = key[i];
    Real prod_i = 1.;
    for (j=0; j<numVars; ++j)
      prod_i *= (randomVarsKey[j]) ? gauss_wts_1d[j][key_i[j]] :
	poly_basis_0[j].get_value(x[j], key_i[j]);
    tpMean += expansionCoeffs[i] * prod_i;
  }
  return tpMean;
}


/** Overloaded version supporting Smolyak sparse grids. */
const Real& InterpPolyApproximation::
tensor_product_mean(const RealVector& x, size_t tp_index)
{
  SparseGridDriver*    ssg_driver = (SparseGridDriver*)driverRep;
  const UShort2DArray& key        = ssg_driver->collocation_key()[tp_index];
  const UShortArray&   sm_index   = ssg_driver->smolyak_multi_index()[tp_index];
  const SizetArray&    coeff_index
    = ssg_driver->expansion_coefficient_indices()[tp_index];
  const Real3DArray&   gauss_wts_1d = driverRep->gauss_weights_array();

  tpMean = 0.;
  size_t i, j, num_colloc_pts = key.size();
  for (i=0; i<num_colloc_pts; ++i) {
    const UShortArray& key_i = key[i];
    Real prod_i = 1.;
    for (j=0; j<numVars; ++j)
      prod_i *= (randomVarsKey[j]) ? gauss_wts_1d[sm_index[j]][j][key_i[j]] :
	polynomialBasis[sm_index[j]][j].get_value(x[j], key_i[j]);
    tpMean += expansionCoeffs[coeff_index[i]] * prod_i;
  }
  return tpMean;
}


/** Overloaded version supporting tensor-product quadrature. */
const RealVector& InterpPolyApproximation::
tensor_product_mean_gradient(const RealVector& x, const UIntArray& dvv)
{
  TensorProductDriver* tpq_driver   = (TensorProductDriver*)driverRep;
  const UShort2DArray& key          = tpq_driver->collocation_key();
  const Real2DArray&   gauss_wts_1d = driverRep->gauss_weights_array()[0];

  size_t i, j, k, deriv_index, num_deriv_vars = dvv.size(),
    cntr = 0; // insertions in cntr order w/i tpCoeffGrads
  if (tpMeanGrad.length() != num_deriv_vars)
    tpMeanGrad.sizeUninitialized(num_deriv_vars);
  tpMeanGrad = 0.;
  std::vector<BasisPolynomial>& poly_basis_0 = polynomialBasis[0];
  SizetList::iterator it;
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    // Error check for required data
    if (randomVarsKey[deriv_index] && !expansionCoeffGradFlag) {
      PCerr << "Error: expansion coefficient gradients not defined in "
	    << "InterpPolyApproximation::get_mean_gradient()." << std::endl;
      abort_handler(-1);
    }
    else if (!randomVarsKey[deriv_index] && !expansionCoeffFlag) {
      PCerr << "Error: expansion coefficients not defined in "
	    << "InterpPolyApproximation::get_mean_gradient()" << std::endl;
      abort_handler(-1);
    }
    for (j=0; j<numCollocPts; ++j) {
      const UShortArray& key_j = key[j];
      Real wt_prod_j = 1.;
      for (it=randomIndices.begin(); it!=randomIndices.end(); it++)
	{ k = *it; wt_prod_j *= gauss_wts_1d[k][key_j[k]]; }
      if (randomVarsKey[deriv_index]) {
	// --------------------------------------------------------------------
	// derivative of All var expansion w.r.t. random var (design insertion)
	// --------------------------------------------------------------------
	Real Lsa_j = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); it++)
	  { k = *it; Lsa_j *= poly_basis_0[k].get_value(x[k], key_j[k]); }
	tpMeanGrad[i] += wt_prod_j * Lsa_j * expansionCoeffGrads(cntr, j);
      }
      else {
	// ---------------------------------------------------------------------
	// deriv of All var expansion w.r.t. nonrandom var (design augmentation)
	// ---------------------------------------------------------------------
	Real dLsa_j_dsa_i = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); it++) {
	  k = *it;
	  dLsa_j_dsa_i *= (k == deriv_index) ?
	    poly_basis_0[k].get_gradient(x[k], key_j[k]) :
	    poly_basis_0[k].get_value(x[k],    key_j[k]);
	}
	tpMeanGrad[i] += expansionCoeffs[j] * wt_prod_j * dLsa_j_dsa_i;
      }
    }
    if (randomVarsKey[deriv_index]) // deriv w.r.t. design var insertion
      cntr++;
  }

  return tpMeanGrad;
}


/** Overloaded version supporting Smolyak sparse grids. */
const RealVector& InterpPolyApproximation::
tensor_product_mean_gradient(const RealVector& x, size_t tp_index,
			     const UIntArray& dvv)
{
  SparseGridDriver*    ssg_driver = (SparseGridDriver*)driverRep;
  const UShort2DArray& key        = ssg_driver->collocation_key()[tp_index];
  const UShortArray&   sm_index   = ssg_driver->smolyak_multi_index()[tp_index];
  const SizetArray&    coeff_index
    = ssg_driver->expansion_coefficient_indices()[tp_index];
  const Real3DArray&   gauss_wts_1d = driverRep->gauss_weights_array();

  // -------------------------------------------------------------------
  // Mixed variable key:
  //   xi = ran vars, sa = augmented des vars, si = inserted design vars
  // Active variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(sa, si) L_i(xi)
  //   mu(sa, si)    = Sum_i r_i(sa, si) wt_prod_i
  //   dmu/ds        = Sum_i dr_i/ds wt_prod_i
  // All variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(si) L_i(xi, sa)
  //   mu(sa, si)    = Sum_i r_i(si) Lsa_i wt_prod_i
  //   dmu/dsa       = Sum_i r_i(si) dLsa_i/dsa wt_prod_i
  //   dmu/dsi       = Sum_i dr_i/dsi Lsa_i wt_prod_i
  // -------------------------------------------------------------------
  size_t i, j, k, deriv_index,
    cntr = 0, // insertions in cntr order w/i tpCoeffGrads
    num_deriv_vars = dvv.size(), num_colloc_pts = key.size();
  if (tpMeanGrad.length() != num_deriv_vars)
    tpMeanGrad.sizeUninitialized(num_deriv_vars);
  tpMeanGrad = 0.;
  SizetList::iterator it;
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    // Error check for required data
    if (randomVarsKey[deriv_index] && !expansionCoeffGradFlag) {
      PCerr << "Error: expansion coefficient gradients not defined in "
	    << "InterpPolyApproximation::get_mean_gradient()." << std::endl;
      abort_handler(-1);
    }
    else if (!randomVarsKey[deriv_index] && !expansionCoeffFlag) {
      PCerr << "Error: expansion coefficients not defined in "
	    << "InterpPolyApproximation::get_mean_gradient()" << std::endl;
      abort_handler(-1);
    }
    for (j=0; j<num_colloc_pts; ++j) {
      const UShortArray& key_j = key[j];
      Real wt_prod_j = 1.;
      for (it=randomIndices.begin(); it!=randomIndices.end(); it++) {
	k = *it;
	wt_prod_j *= gauss_wts_1d[sm_index[k]][k][key_j[k]];
      }
      if (randomVarsKey[deriv_index]) {
	// --------------------------------------------------------------------
	// derivative of All var expansion w.r.t. random var (design insertion)
	// --------------------------------------------------------------------
	Real Lsa_j = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); it++) {
	  k = *it;
	  Lsa_j *= polynomialBasis[sm_index[k]][k].get_value(x[k], key_j[k]);
	}
	tpMeanGrad[i]
	  += wt_prod_j * Lsa_j * expansionCoeffGrads(cntr, coeff_index[j]);
      }
      else {
	// ---------------------------------------------------------------------
	// deriv of All var expansion w.r.t. nonrandom var (design augmentation)
	// ---------------------------------------------------------------------
	Real dLsa_j_dsa_i = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); it++) {
	  k = *it;
	  dLsa_j_dsa_i *= (k == deriv_index) ?
	    polynomialBasis[sm_index[k]][k].get_gradient(x[k], key_j[k]) :
	    polynomialBasis[sm_index[k]][k].get_value(x[k],    key_j[k]);
	}
	tpMeanGrad[i]
	  += expansionCoeffs[coeff_index[j]] * wt_prod_j * dLsa_j_dsa_i;
      }
    }
    if (randomVarsKey[deriv_index]) // deriv w.r.t. design var insertion
      cntr++;
  }

  return tpMeanGrad;
}


/** Overloaded version supporting tensor-product quadrature. */
const Real& InterpPolyApproximation::
tensor_product_variance(const RealVector& x)
{
  TensorProductDriver* tpq_driver   = (TensorProductDriver*)driverRep;
  const UShort2DArray& key          = tpq_driver->collocation_key();
  const Real2DArray&   gauss_wts_1d = driverRep->gauss_weights_array()[0];

  tpVariance = 0.;
  size_t i, j, k;
  Real mean = 0.;
  SizetList::iterator it;
  std::vector<BasisPolynomial>& poly_basis_0 = polynomialBasis[0];
  for (i=0; i<numCollocPts; ++i) {
    const UShortArray& key_i = key[i];
    Real wt_prod_i = 1., Ls_prod_i = 1.;
    for (k=0; k<numVars; ++k)
      if (randomVarsKey[k])
	wt_prod_i *= gauss_wts_1d[k][key_i[k]];
      else
	Ls_prod_i *= poly_basis_0[k].get_value(x[k], key_i[k]);
    const Real& exp_coeff_i = expansionCoeffs[i];
    mean += exp_coeff_i * wt_prod_i * Ls_prod_i;
    for (j=0; j<numCollocPts; ++j) {
      const UShortArray& key_j = key[j];
      // to include the ij-th term, colloc pts xi_i must be the same as xi_j
      // for the ran var subset.  In this case, wt_prod_i may be reused.
      bool include = true;
      for (it=randomIndices.begin(); it!=randomIndices.end(); it++)
	if (key_i[*it] != key_j[*it])
	  { include = false; break; }
      if (include) {
	Real Ls_prod_j = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); it++) {
	  k = *it;
	  Ls_prod_j *= poly_basis_0[k].get_value(x[k], key_j[k]);
	}
	tpVariance += wt_prod_i * Ls_prod_i * Ls_prod_j * exp_coeff_i
	  * expansionCoeffs[j];
      }
    }
  }
  tpVariance -= mean*mean;
  return tpVariance;
}


/** Overloaded version supporting Smolyak sparse grids. */
const Real& InterpPolyApproximation::
tensor_product_variance(const RealVector& x, size_t tp_index)
{
  SparseGridDriver*    ssg_driver = (SparseGridDriver*)driverRep;
  const UShort2DArray& key        = ssg_driver->collocation_key()[tp_index];
  const UShortArray&   sm_index   = ssg_driver->smolyak_multi_index()[tp_index];
  const SizetArray&    coeff_index
    = ssg_driver->expansion_coefficient_indices()[tp_index];
  const Real3DArray&   gauss_wts_1d = driverRep->gauss_weights_array();

  size_t i, j, k, num_colloc_pts = key.size();
  Real mean = 0.;
  SizetList::iterator it;
  tpVariance = 0.;
  for (i=0; i<num_colloc_pts; ++i) {
    const UShortArray& key_i = key[i];
    Real wt_prod_i = 1., Ls_prod_i = 1.;
    for (k=0; k<numVars; ++k)
      if (randomVarsKey[k])
	wt_prod_i *= gauss_wts_1d[sm_index[k]][k][key_i[k]];
      else
	Ls_prod_i *= polynomialBasis[sm_index[k]][k].get_value(x[k], key_i[k]);
    const Real& exp_coeff_i = expansionCoeffs[coeff_index[i]];
    mean += exp_coeff_i * wt_prod_i * Ls_prod_i;
    for (j=0; j<num_colloc_pts; ++j) {
      const UShortArray& key_j = key[j];
      // to include the ij-th term, colloc pts xi_i must be the same as xi_j
      // for the ran var subset.  In this case, wt_prod_i may be reused.
      bool include = true;
      for (it=randomIndices.begin(); it!=randomIndices.end(); it++)
	if (key_i[*it] != key_j[*it])
	  { include = false; break; }
      if (include) {
	Real Ls_prod_j = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); it++) {
	  k = *it;
	  Ls_prod_j *= polynomialBasis[sm_index[k]][k].get_value(x[k],key_j[k]);
	}
	tpVariance += wt_prod_i * Ls_prod_i * Ls_prod_j * exp_coeff_i
	  * expansionCoeffs[coeff_index[j]];
      }
    }
  }
  tpVariance -= mean*mean;
  return tpVariance;
}


/** Overloaded version supporting tensor-product quadrature. */
const RealVector& InterpPolyApproximation::
tensor_product_variance_gradient(const RealVector& x, const UIntArray& dvv)
{
  TensorProductDriver* tpq_driver   = (TensorProductDriver*)driverRep;
  const UShort2DArray& key          = tpq_driver->collocation_key();
  const Real2DArray&   gauss_wts_1d = driverRep->gauss_weights_array()[0];

  Real mean = 0.;
  size_t i, j, k, l, deriv_index, num_deriv_vars = dvv.size(),
    cntr = 0; // insertions in cntr order w/i expansionCoeffGrads
  if (tpVarianceGrad.length() != num_deriv_vars)
    tpVarianceGrad.sizeUninitialized(num_deriv_vars);
  tpVarianceGrad = 0.;
  SizetList::iterator it;
  std::vector<BasisPolynomial>& poly_basis_0 = polynomialBasis[0];
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    if (randomVarsKey[deriv_index] && !expansionCoeffGradFlag) {
      PCerr << "Error: expansion coefficient gradients not defined in "
	    << "InterpPolyApproximation::get_variance_gradient()." << std::endl;
      abort_handler(-1);
    }
    Real mean_grad_i = 0.;
    // first loop of double sum
    for (j=0; j<numCollocPts; ++j) {
      const UShortArray& key_j = key[j];
      // compute wt_prod_j and Lsa_j
      Real wt_prod_j = 1., Lsa_j = 1., dLsa_j_dsa_i = 1.;
      for (k=0; k<numVars; ++k)
	if (randomVarsKey[k])
	  wt_prod_j *= gauss_wts_1d[k][key_j[k]];
	else
	  Lsa_j *= poly_basis_0[k].get_value(x[k], key_j[k]);
      // update mean (once) and mean_grad_i
      if (i == 0)
	mean += expansionCoeffs[j] * wt_prod_j * Lsa_j;
      if (randomVarsKey[deriv_index])
	mean_grad_i += wt_prod_j * Lsa_j * expansionCoeffGrads(cntr, j);
      else {
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); it++) {
	  k = *it;
	  dLsa_j_dsa_i *= (k == deriv_index) ?
	    poly_basis_0[k].get_gradient(x[k], key_j[k]) :
	    poly_basis_0[k].get_value(x[k],    key_j[k]);
	}
	mean_grad_i += wt_prod_j * dLsa_j_dsa_i * expansionCoeffs[j];
      }
      // second loop of double sum
      for (k=0; k<numCollocPts; ++k) {
	const UShortArray& key_k = key[k];
	// to include the jk-th term, colloc pts xi_j must be the same as xi_k
	// for the random var subset.  In this case, wt_prod_j may be reused.
	bool include = true;
	for (it=randomIndices.begin(); it!=randomIndices.end(); it++)
	  if (key_j[*it] != key_k[*it])
	    { include = false; break; }
	if (include) {
	  Real Lsa_k = 1.;
	  for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); it++) {
	    l = *it;
	    Lsa_k *= poly_basis_0[l].get_value(x[l], key_k[l]);
	  }
	  if (randomVarsKey[deriv_index])
	    // ---------------------------------------------------------------
	    // deriv of All var expansion w.r.t. random var (design insertion)
	    // ---------------------------------------------------------------
	    tpVarianceGrad[i] += wt_prod_j * Lsa_j * Lsa_k *
	      (expansionCoeffs[j] * expansionCoeffGrads(cntr, k) +
	       expansionCoeffGrads(cntr, j) * expansionCoeffs[k]);
	  else {
	    // ---------------------------------------------------------------
	    // deriv of All var exp w.r.t. nonrandom var (design augmentation)
	    // ---------------------------------------------------------------
	    Real dLsa_k_dsa_i = 1.;
	    for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); it++){
	      l = *it;
	      dLsa_k_dsa_i *= (l == deriv_index) ?
		poly_basis_0[l].get_gradient(x[l], key_k[l]):
		poly_basis_0[l].get_value(x[l],    key_k[l]);
	    }
	    tpVarianceGrad[i] +=
	      wt_prod_j * expansionCoeffs[j] * expansionCoeffs[k] *
	      (Lsa_j * dLsa_k_dsa_i + dLsa_j_dsa_i * Lsa_k);
	  }
	}
      }
    }

    // subtract 2 mu dmu/ds
    tpVarianceGrad[i] -= 2. * mean * mean_grad_i;

    if (randomVarsKey[deriv_index]) // deriv w.r.t. design var insertion
      cntr++;
  }
  return tpVarianceGrad;
}


/** Overloaded version supporting Smolyak sparse grids. */
const RealVector& InterpPolyApproximation::
tensor_product_variance_gradient(const RealVector& x, size_t tp_index,
				 const UIntArray& dvv)
{
  SparseGridDriver*    ssg_driver = (SparseGridDriver*)driverRep;
  const UShort2DArray& key        = ssg_driver->collocation_key()[tp_index];
  const UShortArray&   sm_index   = ssg_driver->smolyak_multi_index()[tp_index];
  const SizetArray&    coeff_index
    = ssg_driver->expansion_coefficient_indices()[tp_index];
  const Real3DArray&   gauss_wts_1d = driverRep->gauss_weights_array();

  // -------------------------------------------------------------------
  // Mixed variable key:
  //   xi = ran vars, sa = augmented des vars, si = inserted design vars
  // Active variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(sa, si) L_i(xi)
  //   var(sa, si)   = Sum_i Sum_j r_i(sa, si) r_j(sa, si) wt_prod_ij - mu^2
  //   dvar/ds       = Sum_i Sum_j (r_i dr_j/ds + dr_i/ds r_j) wt_prod_ij -
  //                   2 mu dmu/ds
  // All variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(si) L_i(xi, sa)
  //   var(sa, si)   = Sum_i Sum_j r_i(si) r_j(si) Lsa_i Lsa_j wt_prod_ij -
  //                   mu^2
  //   dvar/dsa      = Sum_i Sum_j r_i(si) r_j(si) (Lsa_i dLsa_j/dsa +
  //                   dLsa_i/dsa Lsa_j) wt_prod_ij - 2 mu dmu/dsa
  //   dvar/dsi      = Sum_i Sum_j (r_i dr_j/dsi + dr_i/dsi r_j)
  //                   Lsa_i Lsa_j wt_prod_ij - 2 mu dmu/dsi
  // -------------------------------------------------------------------
  size_t i, j, k, l, deriv_index,
    cntr = 0, // insertions in cntr order w/i expansionCoeffGrads
    num_deriv_vars = dvv.size(), num_colloc_pts = key.size();
  if (tpVarianceGrad.length() != num_deriv_vars)
    tpVarianceGrad.sizeUninitialized(num_deriv_vars);
  tpVarianceGrad = 0.;
  Real mean = 0.;
  SizetList::iterator it;
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    if (randomVarsKey[deriv_index] && !expansionCoeffGradFlag) {
      PCerr << "Error: expansion coefficient gradients not defined in "
	    << "InterpPolyApproximation::get_variance_gradient()." << std::endl;
      abort_handler(-1);
    }
    Real mean_grad_i = 0.;
    // first loop of double sum
    for (j=0; j<num_colloc_pts; ++j) {
      const UShortArray& key_j = key[j];
      // compute wt_prod_j and Lsa_j
      Real wt_prod_j = 1., Lsa_j = 1., dLsa_j_dsa_i = 1.;
      for (k=0; k<numVars; ++k)
	if (randomVarsKey[k])
	  wt_prod_j *= gauss_wts_1d[sm_index[k]][k][key_j[k]];
	else
	  Lsa_j *= polynomialBasis[sm_index[k]][k].get_value(x[k], key_j[k]);
      // update mean (once) and mean_grad_i
      if (i == 0)
	mean += expansionCoeffs[coeff_index[j]] * wt_prod_j * Lsa_j;
      if (randomVarsKey[deriv_index])
	mean_grad_i
	  += wt_prod_j * Lsa_j * expansionCoeffGrads(cntr,coeff_index[j]);
      else {
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); it++) {
	  k = *it;
	  dLsa_j_dsa_i *= (k == deriv_index) ?
	    polynomialBasis[sm_index[k]][k].get_gradient(x[k], key_j[k]) :
	    polynomialBasis[sm_index[k]][k].get_value(x[k],    key_j[k]);
	}
	mean_grad_i
	  += wt_prod_j * dLsa_j_dsa_i * expansionCoeffs[coeff_index[j]];
      }
      // second loop of double sum
      for (k=0; k<num_colloc_pts; ++k) {
	const UShortArray& key_k = key[k];
	// to include the jk-th term, colloc pts xi_j must be the same as xi_k
	// for the random var subset.  In this case, wt_prod_j may be reused.
	bool include = true;
	for (it=randomIndices.begin(); it!=randomIndices.end(); it++)
	  if (key_j[*it] != key_k[*it])
	    { include = false; break; }
	if (include) {
	  Real Lsa_k = 1.;
	  for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); it++) {
	    l = *it;
	    Lsa_k *= polynomialBasis[sm_index[l]][l].get_value(x[l], key_k[l]);
	  }
	  if (randomVarsKey[deriv_index])
	    // ---------------------------------------------------------------
	    // deriv of All var expansion w.r.t. random var (design insertion)
	    // ---------------------------------------------------------------
	    tpVarianceGrad[i] += wt_prod_j * Lsa_j * Lsa_k *
	      (expansionCoeffs[coeff_index[j]] *
	       expansionCoeffGrads(cntr,coeff_index[k]) +
	       expansionCoeffGrads(cntr,coeff_index[j]) *
	       expansionCoeffs[coeff_index[k]]);
	  else {
	    // ---------------------------------------------------------------
	    // deriv of All var exp w.r.t. nonrandom var (design augmentation)
	    // ---------------------------------------------------------------
	    Real dLsa_k_dsa_i = 1.;
	    for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); it++){
	      l = *it;
	      dLsa_k_dsa_i *= (l == deriv_index) ?
		polynomialBasis[sm_index[l]][l].get_gradient(x[l], key_k[l]):
		polynomialBasis[sm_index[l]][l].get_value(x[l],    key_k[l]);
	    }
	    tpVarianceGrad[i] +=
	      wt_prod_j * expansionCoeffs[coeff_index[j]] *
	      expansionCoeffs[coeff_index[k]] *
	      (Lsa_j * dLsa_k_dsa_i + dLsa_j_dsa_i * Lsa_k);
	  }
	}
      }
    }

    // subtract 2 mu dmu/ds
    tpVarianceGrad[i] -= 2. * mean * mean_grad_i;

    if (randomVarsKey[deriv_index]) // deriv w.r.t. design var insertion
      cntr++;
  }
  return tpVarianceGrad;
}


const Real& InterpPolyApproximation::get_value(const RealVector& x)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "InterpPolyApproximation::get_value()" << std::endl;
    abort_handler(-1);
  }

  // sum expansion to get response prediction
  switch (expCoeffsSolnApproach) {
  case QUADRATURE:
    return tensor_product_value(x);
    break;
  case SPARSE_GRID: {
    // Smolyak recursion of anisotropic tensor products
    approxValue = 0.;
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    const IntArray&   sm_coeffs  = ssg_driver->smolyak_coefficients();
    size_t i, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i)
      approxValue += sm_coeffs[i] * tensor_product_value(x, i);
    return approxValue;
    break;
  }
  }
}


const RealVector& InterpPolyApproximation::get_gradient(const RealVector& x)
{
  // this could define a default_dvv and call get_gradient(x, dvv),
  // but we want this fn to be as fast as possible

  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "InterpPolyApproximation::get_gradient()" << std::endl;
    abort_handler(-1);
  }

  switch (expCoeffsSolnApproach) {
  case QUADRATURE:
    return tensor_product_gradient(x);
    break;
  case SPARSE_GRID: {
    if (approxGradient.length() != numVars)
      approxGradient.sizeUninitialized(numVars);
    approxGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    const IntArray&   sm_coeffs  = ssg_driver->smolyak_coefficients();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff = sm_coeffs[i];
      const RealVector& tp_grad = tensor_product_gradient(x, i);
      for (j=0; j<numVars; ++j)
	approxGradient[j] += coeff * tp_grad[j];
    }
    return approxGradient;
    break;
  }
  }
}


const RealVector& InterpPolyApproximation::
get_gradient(const RealVector& x, const UIntArray& dvv)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "InterpPolyApproximation::get_gradient()" << std::endl;
    abort_handler(-1);
  }

  switch (expCoeffsSolnApproach) {
  case QUADRATURE:
    return tensor_product_gradient(x, dvv);
    break;
  case SPARSE_GRID: {
    size_t num_deriv_vars = dvv.size();
    if (approxGradient.length() != num_deriv_vars)
      approxGradient.sizeUninitialized(num_deriv_vars);
    approxGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    const IntArray&   sm_coeffs  = ssg_driver->smolyak_coefficients();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff = sm_coeffs[i];
      const RealVector& tp_grad = tensor_product_gradient(x, i, dvv);
      for (j=0; j<num_deriv_vars; ++j)
	approxGradient[j] += coeff * tp_grad[j];
    }
    return approxGradient;
    break;
  }
  }
}


/** In this case, all expansion variables are random variables and the
    mean of the expansion is simply the sum over i of r_i w_i. */
const Real& InterpPolyApproximation::get_mean()
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "InterpPolyApproximation::get_mean()" << std::endl;
    abort_handler(-1);
  }

  expansionMean = 0.;
  const RealVector& wt_sets = driverRep->weight_sets();
  for (size_t i=0; i<numCollocPts; ++i)
    expansionMean += expansionCoeffs[i] * wt_sets[i];
  return expansionMean;
}


/** In this case, a subset of the expansion variables are random
    variables and the mean of the expansion involves evaluating the
    expectation over this subset. */
const Real& InterpPolyApproximation::get_mean(const RealVector& x)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "InterpPolyApproximation::get_mean()" << std::endl;
    abort_handler(-1);
  }

  switch (expCoeffsSolnApproach) {
  case QUADRATURE:
    return tensor_product_mean(x);
    break;
  case SPARSE_GRID: {
    expansionMean = 0.;
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    const IntArray&   sm_coeffs  = ssg_driver->smolyak_coefficients();
    size_t i, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i)
      expansionMean += sm_coeffs[i] * tensor_product_mean(x, i);
    return expansionMean;
    break;
  }
  }
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  In
    this case, the derivative of the expectation is the expectation of
    the derivative.  The mixed derivative case (some design variables
    are inserted and some are augmented) requires no special treatment. */
const RealVector& InterpPolyApproximation::get_mean_gradient()
{
  // d/ds <R> = <dR/ds>

  // Error check for required data
  if (!expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in "
	  << "InterpPolyApproximation::get_mean_gradient()." << std::endl;
    abort_handler(-1);
  }

  const RealVector& wt_sets = driverRep->weight_sets();
  size_t i, j, num_deriv_vars = expansionCoeffGrads.numRows();
  if (expansionMeanGrad.length() != num_deriv_vars)
    expansionMeanGrad.sizeUninitialized(num_deriv_vars);
  expansionMeanGrad = 0.;
  for (i=0; i<numCollocPts; ++i) {
    const Real& wt_prod_i = wt_sets[i];
    for (j=0; j<num_deriv_vars; ++j)
      expansionMeanGrad[j] += expansionCoeffGrads(j,i) * wt_prod_i;
  }
  return expansionMeanGrad;
}


/** In this function, a subset of the expansion variables are random
    variables and any augmented design/state variables (i.e., not
    inserted as random variable distribution parameters) are included
    in the expansion.  In this case, the mean of the expansion is the
    expectation over the random subset and the derivative of the mean
    is the derivative of the remaining expansion over the non-random
    subset.  This function must handle the mixed case, where some
    design/state variables are augmented (and are part of the
    expansion: derivatives are evaluated as described above) and some
    are inserted (derivatives are obtained from expansionCoeffGrads). */
const RealVector& InterpPolyApproximation::
get_mean_gradient(const RealVector& x, const UIntArray& dvv)
{
  switch (expCoeffsSolnApproach) {
  case QUADRATURE:
    return tensor_product_mean_gradient(x, dvv);
    break;
  case SPARSE_GRID:
    size_t num_deriv_vars = dvv.size();
    if (expansionMeanGrad.length() != num_deriv_vars)
      expansionMeanGrad.sizeUninitialized(num_deriv_vars);
    expansionMeanGrad = 0.;
    // Smolyak recursion of anisotropic tensor products
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    const IntArray&   sm_coeffs  = ssg_driver->smolyak_coefficients();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff = sm_coeffs[i];
      const RealVector& tpm_grad = tensor_product_mean_gradient(x, i, dvv);
      for (j=0; j<num_deriv_vars; ++j)
	expansionMeanGrad[j] += coeff * tpm_grad[j];
    }
    return expansionMeanGrad;
    break;
  }
}


/** In this case, all expansion variables are random variables and the
    variance of the expansion is the sum over all but the first term
    of the coefficients squared times the polynomial norms squared. */
const Real& InterpPolyApproximation::get_variance()
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "InterpPolyApproximation::get_variance()" << std::endl;
    abort_handler(-1);
  }

  return get_covariance(expansionCoeffs);

  /*
  expansionVariance = 0.;
  Real mean = 0.;
  const RealVector& wt_sets = driverRep->weight_sets();
  for (size_t i=0; i<numCollocPts; ++i) {
    const Real& exp_coeff_i = expansionCoeffs[i];
    const Real& wt_prod_i   = wt_sets[i];
    Real coeff_wt_i = exp_coeff_i * wt_prod_i;
    mean              += coeff_wt_i;
    expansionVariance += exp_coeff_i * coeff_wt_i;
  }
  expansionVariance -= mean*mean;
  return expansionVariance;
  */
}


const Real& InterpPolyApproximation::
get_covariance(const RealVector& exp_coeffs_2)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "InterpPolyApproximation::get_variance()" << std::endl;
    abort_handler(-1);
  }

  expansionVariance = 0.;
  Real mean_1 = 0., mean_2 = 0.;
  const RealVector& wt_sets = driverRep->weight_sets();
  for (size_t i=0; i<numCollocPts; ++i) {
    const Real& coeff_2i  = exp_coeffs_2[i];
    const Real& wt_prod_i = wt_sets[i];
    Real coeff_wt_1i = expansionCoeffs[i] * wt_prod_i;
    mean_1 += coeff_wt_1i;
    mean_2 += coeff_2i * wt_prod_i;
    expansionVariance += coeff_wt_1i * coeff_2i;
  }
  expansionVariance -= mean_1*mean_2;
  return expansionVariance;
}


/** In this case, a subset of the expansion variables are random
    variables and the variance of the expansion involves summations
    over this subset. */
const Real& InterpPolyApproximation::get_variance(const RealVector& x)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "InterpPolyApproximation::get_variance()" << std::endl;
    abort_handler(-1);
  }

  switch (expCoeffsSolnApproach) {
  case QUADRATURE:
    return tensor_product_variance(x);
    break;
  case SPARSE_GRID:
    expansionVariance = 0.;
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    const IntArray&   sm_coeffs  = ssg_driver->smolyak_coefficients();
    size_t i, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i)
      expansionVariance	+= sm_coeffs[i] * tensor_product_variance(x, i);
    return expansionVariance;
    break;
  }
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  The
    mixed derivative case (some design variables are inserted and some
    are augmented) requires no special treatment. */
const RealVector& InterpPolyApproximation::get_variance_gradient()
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "InterpPolyApproximation::get_variance_gradient()" << std::endl;
    abort_handler(-1);
  }
  if (!expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in "
	  << "InterpPolyApproximation::get_variance_gradient()." << std::endl;
    abort_handler(-1);
  }

  const RealVector& wt_sets = driverRep->weight_sets();
  size_t i, j, num_deriv_vars = expansionCoeffGrads.numRows();
  if (expansionVarianceGrad.length() != num_deriv_vars)
    expansionVarianceGrad.sizeUninitialized(num_deriv_vars);
  expansionVarianceGrad = 0.;

  Real mean = 0.;
  for (i=0; i<numCollocPts; ++i)
    mean += expansionCoeffs[i] * wt_sets[i];
  for (i=0; i<numCollocPts; ++i) {
    Real term_i = 2. * (expansionCoeffs[i] - mean) * wt_sets[i];
    for (j=0; j<num_deriv_vars; ++j)
      expansionVarianceGrad[j] += term_i * expansionCoeffGrads(j,i);
  }

  return expansionVarianceGrad;
}


/** In this function, a subset of the expansion variables are random
    variables and any augmented design/state variables (i.e., not
    inserted as random variable distribution parameters) are included
    in the expansion.  This function must handle the mixed case, where
    some design/state variables are augmented (and are part of the
    expansion) and some are inserted (derivatives are obtained from
    expansionCoeffGrads). */
const RealVector& InterpPolyApproximation::
get_variance_gradient(const RealVector& x, const UIntArray& dvv)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "InterpPolyApproximation::get_variance_gradient()" << std::endl;
    abort_handler(-1);
  }

  switch (expCoeffsSolnApproach) {
  case QUADRATURE:
    return tensor_product_variance_gradient(x, dvv);
    break;
  case SPARSE_GRID:
    size_t num_deriv_vars = dvv.size();
    if (expansionVarianceGrad.length() != num_deriv_vars)
      expansionVarianceGrad.sizeUninitialized(num_deriv_vars);
    expansionVarianceGrad = 0.;
    // Smolyak recursion of anisotropic tensor products
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    const IntArray&   sm_coeffs  = ssg_driver->smolyak_coefficients();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff = sm_coeffs[i];
      const RealVector& tpv_grad = tensor_product_variance_gradient(x, i, dvv);
      for (j=0; j<num_deriv_vars; ++j)
	expansionVarianceGrad[j] += coeff * tpv_grad[j];
    }
    return expansionVarianceGrad;
    break;
  }
}


// Compute Sobol Indices for global sensitivity analysis
// GT: need to make change here; find constiuent sets according to output verbosity
// STATUS: not complete
void InterpPolyApproximation::compute_global_sensitivity()
{
  if (outputLevel <  NORMAL_OUTPUT)
    return;
  // need the full tensor product mean and tensor product variance
  Real total_mean = get_mean();
  Real total_variance =	get_variance();

  ////////////////////// Start Sort ///////////////////////////

  // size member variables
  // GT: these are sized correctly (follows the length of sobolIndices)
  constituentSets.resize(sobolIndices.length());
  partialVariance.size(sobolIndices.length());

  // perform subset sort
  get_subsets();

  ////////////////////// Begin Computations ////////////////////
	
  // initialize 
  partialVariance[0] = std::pow(total_mean,2.0); // init with mean sq
  totalSobolIndices = 0; // init total indices

  /*// Solve for partial variances
  for (int ct=1; ct<sobolIndices.length(); ct++) {
    // GT: instead of using counter directly, convert to "set number"
    // GT: only compute partial variances of sets of interest
    partial_variance(ct);
    sobolIndices[ct] = 1/total_variance*partialVariance[ct];
    // GT: Looks like total indices simply identify the membership of the sobolIndices 
    // and adds it to the appropriate bin
    for (int k=0; k<std::numeric_limits<int>::digits; ++k)
      if (ct & (1 << k)) // if subset ct contains variable k
	totalSobolIndices[k] += sobolIndices[ct];
  }*/
  // Solve for partial variance
  for (IntIntMIter map_iter=sobolIndexMap.begin(); map_iter!=sobolIndexMap.end(); map_iter++) {
    partial_variance((*map_iter).first);
    sobolIndices[(*map_iter).second] = 1/total_variance*partialVariance[(*map_iter).second];
    // GT: Looks like total indices simply identify the membership of the sobolIndices 
    // and adds it to the appropriate bin
    for (int k=0; k<std::numeric_limits<int>::digits; ++k)
      if ((*map_iter).first & (1 << k)) // if subset ct contains variable k
	totalSobolIndices[k] += sobolIndices[(*map_iter).second];
  }
}


/** Find constituent subsets. */
// GT: need to make change here; find constiuent sets according to output verbosity
// STATUS: not complete
void InterpPolyApproximation::get_subsets()
{
  // includes the "zero" set
  // GT: Represents all possible sets of variablesi
  // GT: change this to "numVars" if quiet or silent
  int num_subsets = sobolIndices.length(); 

  // change this to not look at ONLY the subsets according to 
  // verbosity
  /*for (int i=1; i<num_subsets; ++i) {
    // find all constituent subsets of set i and store
    lower_sets(i,constituentSets[i]);
    // pull out top_level_set from the constituent set
    constituentSets[i].erase(i);	
  }*/
  // Here we want to utilize the integer representation of the subset
  // but we want to store it in a size appropriate container
  // so finding lower sets is given the argument of the integer rep (*.first) and stored in
  // constituentSets in size-appropriate-index-map (*.second)
  for (IntIntMIter map_iter=sobolIndexMap.begin(); map_iter!=sobolIndexMap.end(); map_iter++) {
    lower_sets((*map_iter).first,constituentSets[(*map_iter).second]);
    constituentSets[(*map_iter).second].erase((*map_iter).first);
  }
}


/** For input set, recursively finds constituent subsets with one
    fewer element */
void InterpPolyApproximation::
lower_sets(int plus_one_set, IntSet& top_level_set)
{
  // if this set has been stored before, stop
  if (top_level_set.count(plus_one_set))
    return;
  // otherwise store current set
  else
    top_level_set.insert(plus_one_set);
  // and find lower level sets
  for (int k=0; k<std::numeric_limits<int>::digits; ++k)
    // GT: this performs a bitwise comparison by shifting 1 by k spaces and comparing that to a binary form of plus_one_set; this allows the variable membership using integers instead of a d-array of bools
    if (plus_one_set & (1 << k)) // if subset i contains variable k, remove that variable from the set by converting the bit-form of (1<<k) to an integer and subtract from the plus_one_set
      lower_sets(plus_one_set-(int)std::pow(2.0,k),top_level_set);
}


/** Forms an interpolant over variables that are members of the given set.
    Finds the variance of the interpolant w.r.t. the variables in the set.
    Overloaded version supporting tensor-product quadrature. */
Real InterpPolyApproximation::partial_variance_integral(int set_value)
{
  TensorProductDriver* tpq_driver   = (TensorProductDriver*)driverRep;
  const UShort2DArray& key          = tpq_driver->collocation_key();
  const Real2DArray&   gauss_wts_1d = driverRep->gauss_weights_array()[0];
  const UShortArray&   quad_order   = tpq_driver->quadrature_order();;

  // Distinguish between non-members and members of the given set, set_value
  BoolDeque nonmember_vars(numVars,true); 
  int num_mem_expansion_coeffs = 1; // number of expansion coeffs in
                                    // member-variable-only expansion 
  IntVector indexing_factor; // factors indexing member variables 
  indexing_factor.sizeUninitialized(numVars);

  // create member variable key and get number of expansion coeffs in
  // member-variable-only expansion
  indexing_factor = 1;
	
  for (int k=0; k<std::numeric_limits<int>::digits; ++k) {
    // if subset contains variable k, set key for variable k to true
    if (set_value & (1 << k)) {
      nonmember_vars[k] = false;	
      // information to properly index mem_expansion_coeffs
      indexing_factor[k] = num_mem_expansion_coeffs;
      num_mem_expansion_coeffs *= quad_order[k];
    }	
  }
	
  // Create vector to store new coefficients
  RealVector mem_expansion_coeffs(num_mem_expansion_coeffs),
    mem_weights(num_mem_expansion_coeffs);

  // Perform integration over non-member variables and store indices
  // of new expansion
  size_t i, j;
  for (i=0; i<numCollocPts; ++i) {
    const UShortArray& key_i = key[i];
    size_t mem_expansion_coeffs_index = 0;	
    Real prod_i_nonmembers = 1, prod_i_members = 1;
    for (j=0; j<numVars; ++j) {
      // Convert key to corresponding index on mem_expansion_coeffs
      mem_expansion_coeffs_index += (nonmember_vars[j]) ? 
	0 : key_i[j]*indexing_factor[j];
      // Save the product of the weights of the member and non-member variables 
      if (nonmember_vars[j])
	prod_i_nonmembers *= gauss_wts_1d[j][key_i[j]];
      else
	prod_i_members *= gauss_wts_1d[j][key_i[j]];
    }

    // mem_weights is performed more time than necessary here, but it
    // seems to be the simplest place to put it
    mem_weights[mem_expansion_coeffs_index] = prod_i_members;
    // sort coefficients by the "signature" of the member variables
    // (i.e. mem_expansion_coeffs_index)
    mem_expansion_coeffs[mem_expansion_coeffs_index]
      += expansionCoeffs[i]*prod_i_nonmembers;
  }

  // Now integrate over the remaining variables	
  Real integral = 0;
  for (i=0; i<num_mem_expansion_coeffs; ++i)
    integral += std::pow(mem_expansion_coeffs[i], 2.0)*mem_weights[i];
  return integral;	
}


/** Forms an interpolant over variables that are members of the given set.
    Finds the variance of the interpolant w.r.t. the variables in the set.
    Overloaded version supporting Smolyak sparse grids. */
Real InterpPolyApproximation::
partial_variance_integral(int set_value, size_t tp_index)
{
  SparseGridDriver*    ssg_driver = (SparseGridDriver*)driverRep;
  const UShort2DArray& key        = ssg_driver->collocation_key()[tp_index];
  const UShortArray&   sm_index   = ssg_driver->smolyak_multi_index()[tp_index];
  const SizetArray&    coeff_index
    = ssg_driver->expansion_coefficient_indices()[tp_index];
  const Real3DArray&   gauss_wts_1d = driverRep->gauss_weights_array();

  // Distinguish between non-members and members of the given set, set_value
  BoolDeque nonmember_vars(numVars,true); 
  int num_mem_expansion_coeffs = 1; // number of expansion coeffs in
                                    // member-variable-only expansion 
  IntVector indexing_factor; // factors indexing member variables 
  indexing_factor.sizeUninitialized(numVars);

  // create member variable key and get number of expansion coeffs in
  // member-variable-only expansion
  indexing_factor = 1;
	
  UShortArray quad_order;
  ssg_driver->level_to_order(sm_index, quad_order);
  for (int k=0; k<std::numeric_limits<int>::digits; ++k) {
    // if subset contains variable k, set key for variable k to true
    if (set_value & (1 << k)) {
      nonmember_vars[k] = false;	
      // information to properly index mem_expansion_coeffs
      indexing_factor[k] = num_mem_expansion_coeffs;
      num_mem_expansion_coeffs *= quad_order[k];
    }	
  }
	
  // Create vector to store new coefficients
  RealVector mem_expansion_coeffs(num_mem_expansion_coeffs),
    mem_weights(num_mem_expansion_coeffs);

  // Perform integration over non-member variables and store indices
  // of new expansion
  size_t i, j, num_colloc_pts = key.size();
  for (i=0; i <num_colloc_pts; ++i) {
    const UShortArray& key_i = key[i];
    size_t mem_expansion_coeffs_index = 0;	
    Real prod_i_nonmembers = 1, prod_i_members = 1;
    for (j=0; j<numVars; ++j) {
      // Convert key to corresponding index on mem_expansion_coeffs
      mem_expansion_coeffs_index += (nonmember_vars[j]) ? 
	0 : key_i[j]*indexing_factor[j];
      // Save the product of the weights of the member and non-member variables 
      if (nonmember_vars[j])
	prod_i_nonmembers *= gauss_wts_1d[sm_index[j]][j][key_i[j]];
      else
	prod_i_members *= gauss_wts_1d[sm_index[j]][j][key_i[j]];
    }

    // mem_weights is performed more time than necessary here, but it
    // seems to be the simplest place to put it
    mem_weights[mem_expansion_coeffs_index] = prod_i_members;
    // sort coefficients by the "signature" of the member variables
    // (i.e. mem_expansion_coeffs_index)
    mem_expansion_coeffs[mem_expansion_coeffs_index]
      += expansionCoeffs[coeff_index[i]]*prod_i_nonmembers;
  }

  // Now integrate over the remaining variables	
  Real integral = 0;
  for (i=0; i<num_mem_expansion_coeffs; ++i)
    integral += std::pow(mem_expansion_coeffs[i], 2.0)*mem_weights[i];
  return integral;	
}


/** Computes the partial expection for a certain set represented in
    integer form.  Solves for lower level subsets if necessary and
    stores all computations in subsetPartialVariance. */
void InterpPolyApproximation::partial_variance(int set_value)
{
  // Computes the integral first
  switch (expCoeffsSolnApproach) {
  case QUADRATURE:
    partialVariance[sobolIndexMap[set_value]] = partial_variance_integral(set_value);
    break;
  case SPARSE_GRID: {
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    const IntArray&   sm_coeffs  = ssg_driver->smolyak_coefficients();
    // Smolyak recursion of anisotropic tensor products
    size_t num_smolyak_indices = sm_coeffs.size();
    for (size_t i=0; i<num_smolyak_indices; ++i)
      partialVariance[sobolIndexMap[set_value]] += sm_coeffs[i]
	* partial_variance_integral(set_value, i);
    break;
  }
  }
  
  // Now subtract the contributions from constituent subsets
  IntSet::iterator itr;
  for (itr  = constituentSets[set_value].begin();
       itr != constituentSets[set_value].end(); itr++) 
    partialVariance[set_value] -= partialVariance[*itr];
}

} // namespace Pecos

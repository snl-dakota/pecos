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
  return (expansionCoeffFlag || expansionGradFlag) ? 1 : 0;
}


void InterpPolyApproximation::allocate_arrays()
{
  PolynomialApproximation::allocate_arrays();

  if (expansionCoeffFlag && expansionCoeffs.length() != numCollocPts)
    expansionCoeffs.sizeUninitialized(numCollocPts);
  if (expansionGradFlag) {
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
  //    (expansionGradFlag && expansionCoeffGrads.numCols() != numCollocPts ) );

  switch (expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    const UShortArray& quad_order = tpq_driver->quadrature_order();
    bool update_exp_form = (quad_order != quadOrderPrev);
    // *** TO DO: capture updates to parameterized/numerical polynomials
    bool update_basis_form = update_exp_form;
    // *** TO DO: carefully evaluate interdependence between exp_form/basis_form

    // verify total number of collocation pts (should not include anchorPoint)
    size_t i, j, num_gauss_pts = 1;
    for (i=0; i<numVars; i++)
      num_gauss_pts *= quad_order[i];
    if (num_gauss_pts != numCollocPts) {
      PCerr << "Error: inconsistent total Gauss point count in "
	    << "InterpPolyApproximation::find_coefficients()" << std::endl;
      abort_handler(-1);
    }

    // only one set of Smolyak indices, all 0's, indicating usage of
    // polynomialBasis[variable][0], with order as specified by quadrature spec.
    if (smolyakMultiIndex.empty()) {
      smolyakMultiIndex.resize(1);
      smolyakMultiIndex[0].assign(numVars, 0);
    }

    if (update_basis_form) {
      // size and initialize polynomialBasis, one interpolant per variable
      if (polynomialBasis.empty()) {
	polynomialBasis.resize(numVars);
	for (i=0; i<numVars; i++)
	  polynomialBasis[i].resize(1);
      }
      const Real3DArray& gauss_pts_1d = driverRep->gauss_points_array();
      for (i=0; i<numVars; i++) {
	bool found = false;
	for (j=0; j<i; j++)
	  if (gauss_pts_1d[i][0] == gauss_pts_1d[j][0]) // equality in pt vector
	    { found = true; break; }
	if (found) // reuse previous instance via shared representation
	  polynomialBasis[i][0] = polynomialBasis[j][0];
	else { // instantiate a new unique instance
	  polynomialBasis[i][0] = BasisPolynomial(LAGRANGE);
	  polynomialBasis[i][0].interpolation_points(gauss_pts_1d[i][0]);
	}
      }
    }

    if (update_exp_form) {
      // define mapping from 1:numCollocPts to set of 1d interpolation indices
      if (collocKey.empty())
	collocKey.resize(1);
      collocKey[0].resize(numCollocPts);
      tensor_product_multi_index(quad_order, collocKey[0]);

      // map indices i to i for the lone tensor product 
      if (expansionCoeffIndices.empty())
	expansionCoeffIndices.resize(1);
      expansionCoeffIndices[0].resize(numCollocPts);
      for (i=0; i<numCollocPts; i++)
	expansionCoeffIndices[0][i] = i;
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
      smolyak_multi_index(smolyakMultiIndex, smolyakCoeffs);
      allocate_collocation_arrays();
    }

    bool update_basis_form = update_exp_form;
    // *** TO DO: carefully evaluate interdependence between exp_form/basis_form
    if (update_basis_form) {
      // size and initialize polynomialBasis, multiple interpolants per variable
      if (polynomialBasis.empty())
	polynomialBasis.resize(numVars);
      // j range is 0:w inclusive; i range is 1:w+1 inclusive
      unsigned short num_levels = ssg_level + 1;
      const Real3DArray& gauss_pts_1d = driverRep->gauss_points_array();
      size_t i, j, k;
      for (i=0; i<numVars; i++) {
	polynomialBasis[i].resize(num_levels);
	for (j=0; j<num_levels; j++) { // j->num_levels-1->ssg_level
	  // anisotropic levels: check other dims at corresponding level
	  bool found = false;
	  for (k=0; k<i; k++)
	    if (//ssg_level[k]   >= j && // level j exists for dimension k
		gauss_pts_1d[i][j] == gauss_pts_1d[k][j]) // pt vector equality
	      { found = true; break; }
	  if (found) // reuse previous instances via shared representations
	    polynomialBasis[i][j] = polynomialBasis[k][j]; // shared rep
	  else { // instantiate new unique instances
	    polynomialBasis[i][j] = BasisPolynomial(LAGRANGE);
	    polynomialBasis[i][j].interpolation_points(gauss_pts_1d[i][j]);
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
  if (!expansionCoeffFlag && !expansionGradFlag) {
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
    if (expansionGradFlag)
      Teuchos::setCol(anchorPoint.response_gradient(), 0, expansionCoeffGrads);
  }

  std::vector<SurrogateDataPoint>::iterator it = dataPoints.begin();
  for (i=offset; i<numCollocPts; ++i, ++it) {
    if (expansionCoeffFlag)
      expansionCoeffs[i] = it->response_function();
    if (expansionGradFlag)
      Teuchos::setCol(it->response_gradient(), (int)i, expansionCoeffGrads);
  }
}


const Real& InterpPolyApproximation::
tensor_product_value(const RealVector& x, size_t tp_index)
{
  tpValue = 0.;
  const UShort2DArray& key         = collocKey[tp_index];
  const UShortArray&   sm_index    = smolyakMultiIndex[tp_index];
  const SizetArray&    coeff_index = expansionCoeffIndices[tp_index];
  size_t i, j, num_colloc_pts = key.size();
  for (i=0; i<num_colloc_pts; i++) {
    const UShortArray& key_i = key[i];
    Real L_i = 1.0;
    for (j=0; j<numVars; j++)
      L_i *= polynomialBasis[j][sm_index[j]].get_value(x[j], key_i[j]);
    tpValue += expansionCoeffs[coeff_index[i]] * L_i;
  }
  return tpValue;
}


const RealVector& InterpPolyApproximation::
tensor_product_gradient(const RealVector& x, size_t tp_index)
{
  if (tpGradient.length() != numVars)
    tpGradient.sizeUninitialized(numVars);
  tpGradient = 0.;

  const UShort2DArray& key         = collocKey[tp_index];
  const UShortArray&   sm_index    = smolyakMultiIndex[tp_index];
  const SizetArray&    coeff_index = expansionCoeffIndices[tp_index];
  size_t i, j, k, num_colloc_pts = key.size();
  for (i=0; i<num_colloc_pts; i++) {
    const UShortArray& key_i   = key[i];
    const Real&        coeff_i = expansionCoeffs[coeff_index[i]];
    for (j=0; j<numVars; j++) {
      Real term_i_grad_j = 1.0;
      for (k=0; k<numVars; k++)
	term_i_grad_j *= (k == j) ?
	  polynomialBasis[k][sm_index[k]].get_gradient(x[k], key_i[k]) :
	  polynomialBasis[k][sm_index[k]].get_value(x[k],    key_i[k]);
      tpGradient[j] += coeff_i * term_i_grad_j;
    }
  }
  return tpGradient;
}


const RealVector& InterpPolyApproximation::
tensor_product_gradient(const RealVector& x, size_t tp_index,
			const UIntArray& dvv)
{
  size_t num_deriv_vars = dvv.size();
  if (tpGradient.length() != num_deriv_vars)
    tpGradient.sizeUninitialized(num_deriv_vars);
  tpGradient = 0.;

  const UShort2DArray& key         = collocKey[tp_index];
  const UShortArray&   sm_index    = smolyakMultiIndex[tp_index];
  const SizetArray&    coeff_index = expansionCoeffIndices[tp_index];
  size_t i, j, k, deriv_index, num_colloc_pts = key.size();
  for (i=0; i<num_colloc_pts; i++) {
    const UShortArray& key_i   = key[i];
    const Real&        coeff_i = expansionCoeffs[coeff_index[i]];
    for (j=0; j<num_deriv_vars; j++) {
      deriv_index = dvv[j] - 1; // requires an "All" view
      Real term_i_grad_j = 1.0;
      for (k=0; k<numVars; k++)
	term_i_grad_j *= (k == deriv_index) ?
	  polynomialBasis[k][sm_index[k]].get_gradient(x[k], key_i[k]) :
	  polynomialBasis[k][sm_index[k]].get_value(x[k],    key_i[k]);
      tpGradient[j] += coeff_i * term_i_grad_j;
    }
  }
  return tpGradient;
}


const Real& InterpPolyApproximation::
tensor_product_mean(const RealVector& x, size_t tp_index)
{
  tpMean = 0.;

  const UShort2DArray& key          = collocKey[tp_index];
  const UShortArray&   sm_index     = smolyakMultiIndex[tp_index];
  const SizetArray&    coeff_index  = expansionCoeffIndices[tp_index];
  const Real3DArray&   gauss_wts_1d = driverRep->gauss_weights_array();

  size_t i, j, num_colloc_pts = key.size();
  for (i=0; i<num_colloc_pts; i++) {
    const UShortArray& key_i = key[i];
    Real prod_i = 1.;
    for (j=0; j<numVars; j++)
      prod_i *= (randomVarsKey[j]) ? gauss_wts_1d[j][sm_index[j]][key_i[j]] :
	polynomialBasis[j][sm_index[j]].get_value(x[j], key_i[j]);
    tpMean += expansionCoeffs[coeff_index[i]] * prod_i;
  }
  return tpMean;
}


const RealVector& InterpPolyApproximation::
tensor_product_mean_gradient(const RealVector& x, size_t tp_index,
			     const UIntArray& dvv)
{
  size_t i, j, k, deriv_index, num_deriv_vars = dvv.size(),
    cntr = 0; // insertions carried in order within tpCoeffGrads
  if (tpMeanGrad.length() != num_deriv_vars)
    tpMeanGrad.sizeUninitialized(num_deriv_vars);
  tpMeanGrad = 0.;

  const UShort2DArray& key          = collocKey[tp_index];
  const UShortArray&   sm_index     = smolyakMultiIndex[tp_index];
  const SizetArray&    coeff_index  = expansionCoeffIndices[tp_index];
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
  size_t num_colloc_pts = key.size();
  SizetList::iterator it;
  for (i=0; i<num_deriv_vars; i++) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    // Error check for required data
    if (randomVarsKey[deriv_index] && !expansionGradFlag) {
      PCerr << "Error: expansion coefficient gradients not defined in "
	    << "InterpPolyApproximation::get_mean_gradient()." << std::endl;
      abort_handler(-1);
    }
    else if (!randomVarsKey[deriv_index] && !expansionCoeffFlag) {
      PCerr << "Error: expansion coefficients not defined in "
	    << "InterpPolyApproximation::get_mean_gradient()" << std::endl;
      abort_handler(-1);
    }
    for (j=0; j<num_colloc_pts; j++) {
      const UShortArray& key_j = key[j];
      Real wt_prod_j = 1.;
      for (it=randomIndices.begin(); it!=randomIndices.end(); it++) {
	k = *it;
	wt_prod_j *= gauss_wts_1d[k][sm_index[k]][key_j[k]];
      }
      if (randomVarsKey[deriv_index]) {
	// -----------------------------------------------------------
	// derivative of All variable expansion w.r.t. random variable
	// (design variable insertion)
	// -----------------------------------------------------------
	Real Lsa_j = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); it++) {
	  k = *it;
	  Lsa_j *= polynomialBasis[k][sm_index[k]].get_value(x[k], key_j[k]);
	}
	tpMeanGrad[i]
	  += wt_prod_j * Lsa_j * expansionCoeffGrads(cntr,coeff_index[j]);
      }
      else {
	// --------------------------------------------------------------
	// derivative of All variable expansion w.r.t. nonrandom variable
	// (design variable augmentation)
	// --------------------------------------------------------------
	Real dLsa_j_dsa_i = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); it++) {
	  k = *it;
	  dLsa_j_dsa_i *= (k == deriv_index) ?
	    polynomialBasis[k][sm_index[k]].get_gradient(x[k], key_j[k]) :
	    polynomialBasis[k][sm_index[k]].get_value(x[k],    key_j[k]);
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


const Real& InterpPolyApproximation::
tensor_product_variance(const RealVector& x, size_t tp_index)
{
  tpVariance = 0.;
  const UShort2DArray& key          = collocKey[tp_index];
  const UShortArray&   sm_index     = smolyakMultiIndex[tp_index];
  const SizetArray&    coeff_index  = expansionCoeffIndices[tp_index];
  const Real3DArray&   gauss_wts_1d = driverRep->gauss_weights_array();

  size_t i, j, k, num_colloc_pts = key.size();
  Real mean = 0.;
  SizetList::iterator it;
  for (i=0; i<num_colloc_pts; i++) {
    const UShortArray& key_i = key[i];
    Real wt_prod_i = 1., Ls_prod_i = 1.;
    for (k=0; k<numVars; k++)
      if (randomVarsKey[k])
	wt_prod_i *= gauss_wts_1d[k][sm_index[k]][key_i[k]];
      else
	Ls_prod_i *= polynomialBasis[k][sm_index[k]].get_value(x[k], key_i[k]);
    const Real& exp_coeff_i = expansionCoeffs[coeff_index[i]];
    mean += exp_coeff_i * wt_prod_i * Ls_prod_i;
    for (j=0; j<num_colloc_pts; j++) {
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
	  Ls_prod_j *= polynomialBasis[k][sm_index[k]].get_value(x[k],key_j[k]);
	}
	tpVariance += wt_prod_i * Ls_prod_i * Ls_prod_j * exp_coeff_i
	  * expansionCoeffs[coeff_index[j]];
      }
    }
  }
  tpVariance -= mean*mean;
  return tpVariance;
}


const RealVector& InterpPolyApproximation::
tensor_product_variance_gradient(const RealVector& x, size_t tp_index,
				 const UIntArray& dvv)
{
  size_t i, j, k, l, deriv_index, num_deriv_vars = dvv.size(),
    cntr = 0; // insertions carried in order within expansionCoeffGrads
  if (tpVarianceGrad.length() != num_deriv_vars)
    tpVarianceGrad.sizeUninitialized(num_deriv_vars);
  tpVarianceGrad = 0.;

  const UShort2DArray& key          = collocKey[tp_index];
  const UShortArray&   sm_index     = smolyakMultiIndex[tp_index];
  const SizetArray&    coeff_index  = expansionCoeffIndices[tp_index];
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
  Real mean = 0.;
  size_t num_colloc_pts = key.size();
  SizetList::iterator it;
  for (i=0; i<num_deriv_vars; i++) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    if (randomVarsKey[deriv_index] && !expansionGradFlag) { // Error check 
      PCerr << "Error: expansion coefficient gradients not defined in "
	    << "InterpPolyApproximation::get_variance_gradient()." << std::endl;
      abort_handler(-1);
    }
    Real mean_grad_i = 0.;
    // first loop of double sum
    for (j=0; j<num_colloc_pts; j++) {
      const UShortArray& key_j = key[j];
      // compute wt_prod_j and Lsa_j
      Real wt_prod_j = 1., Lsa_j = 1., dLsa_j_dsa_i = 1.;
      for (k=0; k<numVars; k++)
	if (randomVarsKey[k])
	  wt_prod_j *= gauss_wts_1d[k][sm_index[k]][key_j[k]];
	else
	  Lsa_j *= polynomialBasis[k][sm_index[k]].get_value(x[k], key_j[k]);
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
	    polynomialBasis[k][sm_index[k]].get_gradient(x[k], key_j[k]) :
	    polynomialBasis[k][sm_index[k]].get_value(x[k],    key_j[k]);
	}
	mean_grad_i
	  += wt_prod_j * dLsa_j_dsa_i * expansionCoeffs[coeff_index[j]];
      }
      // second loop of double sum
      for (k=0; k<num_colloc_pts; k++) {
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
	    Lsa_k *= polynomialBasis[l][sm_index[l]].get_value(x[l], key_k[l]);
	  }
	  if (randomVarsKey[deriv_index])
	    // -----------------------------------------------------------
	    // derivative of All variable expansion w.r.t. random variable
	    // (design variable insertion)
	    // -----------------------------------------------------------
	    tpVarianceGrad[i] += wt_prod_j * Lsa_j * Lsa_k *
	      (expansionCoeffs[coeff_index[j]] *
	       expansionCoeffGrads(cntr,coeff_index[k]) +
	       expansionCoeffGrads(cntr,coeff_index[j]) *
	       expansionCoeffs[coeff_index[k]]);
	  else {
	    // ---------------------------------------------------------
	    // derivative of All variable expansion w.r.t. nonrandom var
	    // (design variable augmentation)
	    // ---------------------------------------------------------
	    Real dLsa_k_dsa_i = 1.;
	    for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); it++){
	      l = *it;
	      dLsa_k_dsa_i *= (l == deriv_index) ?
		polynomialBasis[l][sm_index[l]].get_gradient(x[l], key_k[l]):
		polynomialBasis[l][sm_index[l]].get_value(x[l],    key_k[l]);
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
    return tensor_product_value(x, 0); // first and only collocKey
    break;
  case SPARSE_GRID: {
    // Smolyak recursion of anisotropic tensor products
    approxValue = 0.;
    size_t i, num_smolyak_indices = smolyakMultiIndex.size();
    for (i=0; i<num_smolyak_indices; i++)
      approxValue += smolyakCoeffs[i] * tensor_product_value(x, i);
    return approxValue;
    break;
  }
  }
}


const RealVector& InterpPolyApproximation::
get_gradient(const RealVector& x)
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
    return tensor_product_gradient(x, 0); // first and only collocKey
    break;
  case SPARSE_GRID: {
    if (approxGradient.length() != numVars)
      approxGradient.sizeUninitialized(numVars);
    approxGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    size_t i, j, num_smolyak_indices = smolyakMultiIndex.size();
    for (i=0; i<num_smolyak_indices; i++) {
      const Real&       coeff   = smolyakCoeffs[i];
      const RealVector& tp_grad = tensor_product_gradient(x, i);
      for (j=0; j<numVars; j++)
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
    return tensor_product_gradient(x, 0, dvv); // first and only collocKey
    break;
  case SPARSE_GRID: {
    size_t num_deriv_vars = dvv.size();
    if (approxGradient.length() != num_deriv_vars)
      approxGradient.sizeUninitialized(num_deriv_vars);
    approxGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    size_t i, j, num_smolyak_indices = smolyakMultiIndex.size();
    for (i=0; i<num_smolyak_indices; i++) {
      const Real&       coeff   = smolyakCoeffs[i];
      const RealVector& tp_grad = tensor_product_gradient(x, i, dvv);
      for (j=0; j<num_deriv_vars; j++)
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
  for (size_t i=0; i<numCollocPts; i++)
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
    return tensor_product_mean(x, 0);
    break;
  case SPARSE_GRID: {
    expansionMean = 0.;
    size_t i, num_smolyak_indices = smolyakMultiIndex.size();
    for (i=0; i<num_smolyak_indices; i++)
      expansionMean += smolyakCoeffs[i] * tensor_product_mean(x, i);
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
  if (!expansionGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in "
	  << "InterpPolyApproximation::get_mean_gradient()." << std::endl;
    abort_handler(-1);
  }

  const RealVector& wt_sets = driverRep->weight_sets();
  size_t i, j, num_deriv_vars = expansionCoeffGrads.numRows();
  if (expansionMeanGrad.length() != num_deriv_vars)
    expansionMeanGrad.sizeUninitialized(num_deriv_vars);
  expansionMeanGrad = 0.;
  for (i=0; i<numCollocPts; i++) {
    const Real& wt_prod_i = wt_sets[i];
    for (j=0; j<num_deriv_vars; j++)
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
    return tensor_product_mean_gradient(x, 0, dvv);
    break;
  case SPARSE_GRID:
    size_t num_deriv_vars = dvv.size();
    if (expansionMeanGrad.length() != num_deriv_vars)
      expansionMeanGrad.sizeUninitialized(num_deriv_vars);
    expansionMeanGrad = 0.;
    // Smolyak recursion of anisotropic tensor products
    size_t i, j, num_smolyak_indices = smolyakMultiIndex.size();
    for (i=0; i<num_smolyak_indices; i++) {
      const Real&       coeff    = smolyakCoeffs[i];
      const RealVector& tpm_grad = tensor_product_mean_gradient(x, i, dvv);
      for (j=0; j<num_deriv_vars; j++)
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
  for (size_t i=0; i<numCollocPts; i++) {
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
  for (size_t i=0; i<numCollocPts; i++) {
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
    return tensor_product_variance(x, 0);
    break;
  case SPARSE_GRID:
    expansionVariance = 0.;
    size_t i, num_smolyak_indices = smolyakMultiIndex.size();
    for (i=0; i<num_smolyak_indices; i++)
      expansionVariance	+= smolyakCoeffs[i] * tensor_product_variance(x, i);
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
  if (!expansionGradFlag) {
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
  for (i=0; i<numCollocPts; i++)
    mean += expansionCoeffs[i] * wt_sets[i];
  for (i=0; i<numCollocPts; i++) {
    Real term_i = 2. * (expansionCoeffs[i] - mean) * wt_sets[i];
    for (j=0; j<num_deriv_vars; j++)
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
    return tensor_product_variance_gradient(x, 0, dvv);
    break;
  case SPARSE_GRID:
    size_t num_deriv_vars = dvv.size();
    if (expansionVarianceGrad.length() != num_deriv_vars)
      expansionVarianceGrad.sizeUninitialized(num_deriv_vars);
    expansionVarianceGrad = 0.;
    // Smolyak recursion of anisotropic tensor products
    size_t i, j, num_smolyak_indices = smolyakMultiIndex.size();
    for (i=0; i<num_smolyak_indices; i++) {
      const Real&       coeff    = smolyakCoeffs[i];
      const RealVector& tpv_grad = tensor_product_variance_gradient(x, i, dvv);
      for (j=0; j<num_deriv_vars; j++)
	expansionVarianceGrad[j] += coeff * tpv_grad[j];
    }
    return expansionVarianceGrad;
    break;
  }
}


// Compute Sobol Indices for global sensitivity analysis
void InterpPolyApproximation::compute_global_sensitivity()
{
  // need the full tensor product mean and tensor product variance
  Real total_mean = get_mean();
  Real total_variance =	get_variance();

  ////////////////////// Start Sort ///////////////////////////

  // size member variables
  constituentSets.resize(sobolIndices.length());
  partialVariance.size(sobolIndices.length());

  // perform subset sort
  get_subsets();

  ////////////////////// Begin Computations ////////////////////
	
  // initialize 
  partialVariance[0] = std::pow(total_mean,2.0); // init with mean sq
  totalSobolIndices = 0; // init total indices

  // Solve for partial variances
  for (int ct=1; ct<sobolIndices.length(); ct++){
    partial_variance(ct);
    sobolIndices[ct] = 1/total_variance*partialVariance[ct];
    for (int k=0; k<std::numeric_limits<int>::digits; k++)
      if (ct & (1 << k)) // if subset ct contains variable k
	totalSobolIndices[k] += sobolIndices[ct];
  }
}


/** Find constituent subsets. */
void InterpPolyApproximation::get_subsets()
{
  // includes the "zero" set
  int num_subsets = (int)std::pow(2.,(int)numVars);

  for (int i=1; i<num_subsets; i++) {
    // find all constituent subsets of set i and store
    lower_sets(i,constituentSets[i]);
    // pull out top_level_set from the constituent set
    constituentSets[i].erase(i);	
  }
}


/** For input set, recursively finds constituent subsets with one
    fewer element */
void InterpPolyApproximation::
lower_sets(int plus_one_set, IntSet &top_level_set)
{
  // if this set has been stored before, stop
  if (top_level_set.count(plus_one_set))
    return;
  // otherwise store current set
  else
    top_level_set.insert(plus_one_set);
  // and find lower level sets
  for (int k=0; k<std::numeric_limits<int>::digits; k++)
    if (plus_one_set & (1 << k)) // if subset i contains variable k
      lower_sets(plus_one_set-(int)std::pow(2.0,k),top_level_set);
}


/** Forms an interpolant over variables that are members of the given set.
    Finds the variance of the interpolant w.r.t. the variables in the set. */
Real InterpPolyApproximation::
partial_variance_integral(const int &set_value, size_t tp_index,
			  UShortArray &quad_order)
// quad_order is passed in order to bypass SPARSE_GRID or QUADRATURE
// specification
{
  // Distinguish between non-members and members of the given set, set_value
  BoolDeque nonmember_vars(numVars,true); 
  int num_mem_expansion_coeffs = 1; // number of expansion coeffs in
                                    // member-variable-only expansion 
  IntVector indexing_factor; // factors indexing member variables 
  indexing_factor.sizeUninitialized(numVars);

  // create member variable key and get number of expansion coeffs in
  // member-variable-only expansion
  indexing_factor = 1;
	
  for (int k=0; k<std::numeric_limits<int>::digits; k++){
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

  const UShortArray&   sm_index     = smolyakMultiIndex[tp_index];
  const SizetArray&    coeff_index  = expansionCoeffIndices[tp_index];
  const UShort2DArray& key          = collocKey[tp_index];
  const Real3DArray&   gauss_wts_1d = driverRep->gauss_weights_array();

  size_t num_colloc_pts = key.size();
  size_t i,j;

  // Perform integration over non-member variables and store indices
  // of new expansion
  for (i=0; i <num_colloc_pts; i++) {
    const UShortArray& key_i = key[i];
    size_t mem_expansion_coeffs_index = 0;	
    Real prod_i_nonmembers = 1, prod_i_members = 1;
    for (j=0; j<numVars; j++) {
      // Convert key to corresponding index on mem_expansion_coeffs
      mem_expansion_coeffs_index += (nonmember_vars[j]) ? 
	0 : key_i[j]*indexing_factor[j];
      // Save the product of the weights of the member and non-member variables 
      if (nonmember_vars[j])
	prod_i_nonmembers *= gauss_wts_1d[j][sm_index[j]][key_i[j]];
      else
	prod_i_members *= gauss_wts_1d[j][sm_index[j]][key_i[j]];
    }

    // mem_weights is performed more time than necessary here, but it
    // seems to be the simplest place to put it
    mem_weights[mem_expansion_coeffs_index] = prod_i_members;
    // sort coefficients by the "signature" of the member variables
    // (i.e. mem_expansion_coeffs_index)
    mem_expansion_coeffs[mem_expansion_coeffs_index] +=
      expansionCoeffs[coeff_index[i]]*prod_i_nonmembers;
  }

  // Now integrate over the remaining variables	
  Real integral = 0;
  for (i=0; i<num_mem_expansion_coeffs; i++)
    integral += std::pow(mem_expansion_coeffs[i], 2.0)*mem_weights[i];
  return integral;	
}


/** Computes the partial expection for a certain set represented in
    integer form.  Solves for lower level subsets if necessary and
    stores all computations in subsetPartialVariance. */
void InterpPolyApproximation::partial_variance(const int &set_value)
{
  // Number of quad points are impt to computing integral
  UShortArray quad_order;
  // Computes the integral first
  switch (expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    quad_order = tpq_driver->quadrature_order();
    partialVariance[set_value]
      = partial_variance_integral(set_value,0,quad_order);
    break;
  }
  case SPARSE_GRID: {
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    // Smolyak recursion of anisotropic tensor products
    size_t num_smolyak_indices = smolyakMultiIndex.size();
    for (size_t i=0; i<num_smolyak_indices; i++) {
      ssg_driver->level_to_order(smolyakMultiIndex[i], quad_order);
      partialVariance[set_value] += smolyakCoeffs[i]
	* partial_variance_integral(set_value,i,quad_order);
    }
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

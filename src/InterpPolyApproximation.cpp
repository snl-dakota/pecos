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
//#define INTERPOLATION_TEST


namespace Pecos {

int InterpPolyApproximation::min_coefficients() const
{
  // return the minimum number of data instances required to build the 
  // surface of multiple dimensions
  return (configOptions.expansionCoeffFlag ||
	  configOptions.expansionCoeffGradFlag) ? 1 : 0;
}


void InterpPolyApproximation::allocate_arrays()
{
  allocate_component_effects();
  allocate_total_effects();

  if (configOptions.expansionCoeffFlag &&
      expansionCoeffs.length() != numCollocPts)
    expansionCoeffs.sizeUninitialized(numCollocPts);
  if (configOptions.expansionCoeffGradFlag) {
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
  //  ( (configOptions.expansionCoeffFlag &&
  //     expansionCoeffs.length()      != numCollocPts) ||
  //    (configOptions.expansionCoeffGradFlag &&
  //     expansionCoeffGrads.numCols() != numCollocPts ) );

  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    const UShortArray&   quad_order = tpq_driver->quadrature_order();

    // verify total number of collocation pts (should not include anchorPoint)
    size_t i, j, num_gauss_pts = 1;
    for (i=0; i<numVars; ++i)
      num_gauss_pts *= quad_order[i];
    if (num_gauss_pts != numCollocPts) {
      PCerr << "Error: inconsistent total Gauss point count in "
	    << "InterpPolyApproximation::allocate_arrays()" << std::endl;
      abort_handler(-1);
    }

    //bool update_exp_form = (quad_order != quadOrderPrev);
    // *** TO DO: capture updates to parameterized/numerical polynomials
    //if (update_exp_form) {
    //}

    // can't use quad_order > quadOrderPrev logic since only 1 pt set is stored
    bool update_basis_form = (quad_order != quadOrderPrev);
    // *** TO DO: carefully evaluate interdependence between exp_form/basis_form
    if (update_basis_form) {
      // size and initialize polynomialBasis, one interpolant per variable
      if (polynomialBasis.empty())
	{ polynomialBasis.resize(1); polynomialBasis[0].resize(numVars); }
      const Real2DArray& gauss_pts_1d = tpq_driver->gauss_points_array();
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

    //bool update_exp_form
    //  = (ssg_level != ssgLevelPrev || aniso_wts != ssgAnisoWtsPrev);
    // *** TO DO: capture updates to parameterized/numerical polynomials
    //if (update_exp_form) {
    //}

    // Ignore weights since they only reduce the interpolation depth from the
    // level and the basis update uses a coarse increment based on level.  This
    // matches isotropic sparse grids, but forces fewer and larger updates in
    // the case of anisotropic or generalized grids.
    bool update_basis_form
      = (ssgLevelPrev == USHRT_MAX || ssg_level > ssgLevelPrev);
    if (update_basis_form)
      update_sparse_interpolation_basis(ssg_level);

    ssgLevelPrev = ssg_level; ssgAnisoWtsPrev = aniso_wts;
    break;
  }
  }
}


void InterpPolyApproximation::compute_coefficients()
{
  if (!configOptions.expansionCoeffFlag &&
      !configOptions.expansionCoeffGradFlag) {
    PCerr << "Warning: neither expansion coefficients nor expansion "
	  << "coefficient gradients\n         are active in "
	  << "InterpPolyApproximation::compute_coefficients().\n         "
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
	  << "InterpPolyApproximation::compute_coefficients()." << std::endl;
    abort_handler(-1);
  }

  allocate_arrays();

  if (!anchorPoint.is_null()) {
    if (configOptions.expansionCoeffFlag)
      expansionCoeffs[0] = anchorPoint.response_function();
    if (configOptions.expansionCoeffGradFlag)
      Teuchos::setCol(anchorPoint.response_gradient(), 0, expansionCoeffGrads);
  }

  std::vector<SurrogateDataPoint>::iterator it = dataPoints.begin();
  for (i=offset; i<numCollocPts; ++i, ++it) {
    if (configOptions.expansionCoeffFlag)
      expansionCoeffs[i] = it->response_function();
    if (configOptions.expansionCoeffGradFlag)
      Teuchos::setCol(it->response_gradient(), (int)i, expansionCoeffGrads);
  }

#ifdef INTERPOLATION_TEST
  // SC should accurately interpolate the collocation data for TPQ and
  // SSG with fully nested rules, but will exhibit interpolation error
  // for SSG with other rules.
  it = dataPoints.begin();
  for (i=offset; i<numCollocPts; ++i, ++it)
    PCout << "Colloc pt " << i+1 << ": coeff = " << expansionCoeffs[i]
	  << " interpolation error = " << std::abs(expansionCoeffs[i] -
	     get_value(it->continuous_variables())) << '\n';
#endif // INTERPOLATION_TEST
}


void InterpPolyApproximation::increment_coefficients()
{
  bool err_flag = false;
  switch (configOptions.expCoeffsSolnApproach) {
  case SPARSE_GRID: {
    // As for allocate_arrays(), increments are performed in coarser steps
    // than may be strictly necessary: all increments are filled in for all
    // vars for a step in level (ignoring anisotropy or generalized indices).
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    const UShortArray& trial_set = ssg_driver->trial_index_set();
    unsigned short max_trial_index = 0;
    for (size_t i=0; i<numVars; ++i)
      if (trial_set[i] > max_trial_index)
	max_trial_index = trial_set[i];
    update_sparse_interpolation_basis(max_trial_index);
    break;
  }
  default:
    err_flag = true; break;
  }
  if (err_flag) {
    PCerr << "Error: unsupported grid definition in InterpPolyApproximation::"
	  << "increment_coefficients()" << std::endl;
    abort_handler(-1);
  }

  restore_expansion_coefficients();
}


void InterpPolyApproximation::decrement_coefficients()
{
  // leave polynomialBasis as is

  switch (configOptions.expCoeffsSolnApproach) {
  case SPARSE_GRID: {
    // move previous expansion data to current expansion
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    savedSmolyakMultiIndex.push_back(ssg_driver->trial_index_set());
    break;
  }
  }

  // not necessary to prune; next increment/restore/finalize takes care of this
  //if (configOptions.expansionCoeffFlag)
  //  expansionCoeffs.resize(numCollocPts);
  //if (configOptions.expansionCoeffGradFlag) {
  //  size_t num_deriv_vars = expansionCoeffGrads.numRows();
  //  expansionCoeffGrads.reshape(num_deriv_vars, numCollocPts);
  //}

  numCollocPts = dataPoints.size(); // data already decremented
  if (!anchorPoint.is_null())
    numCollocPts += 1;
}


void InterpPolyApproximation::restore_coefficients()
{
  // leave polynomialBasis as is (a previous increment is being restored)

  switch (configOptions.expCoeffsSolnApproach) {
  case SPARSE_GRID: {
    // move previous expansion data to current expansion
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    std::deque<UShortArray>::iterator sit
      = std::find(savedSmolyakMultiIndex.begin(), savedSmolyakMultiIndex.end(),
		  ssg_driver->trial_index_set());
    if (sit != savedSmolyakMultiIndex.end())
      savedSmolyakMultiIndex.erase(sit);
    break;
  }
  }

  restore_expansion_coefficients();
}


void InterpPolyApproximation::finalize_coefficients()
{
  // leave polynomialBasis as is (all previous increments are being restored)

  switch (configOptions.expCoeffsSolnApproach) {
  case SPARSE_GRID:
    // move previous expansion data to current expansion
    savedSmolyakMultiIndex.clear();
    break;
  }

  restore_expansion_coefficients();
}


void InterpPolyApproximation::restore_expansion_coefficients()
{
  size_t i, new_colloc_pts = dataPoints.size();
  if (!anchorPoint.is_null())
    new_colloc_pts += 1;

  if (configOptions.expansionCoeffFlag)
    expansionCoeffs.resize(new_colloc_pts);
  if (configOptions.expansionCoeffGradFlag) {
    size_t num_deriv_vars = expansionCoeffGrads.numRows();
    expansionCoeffGrads.reshape(num_deriv_vars, new_colloc_pts);
  }

  std::vector<SurrogateDataPoint>::iterator it = dataPoints.begin();
  std::advance(it, numCollocPts);
  for (i=numCollocPts; i<new_colloc_pts; ++i, ++it) {
    if (configOptions.expansionCoeffFlag)
      expansionCoeffs[i] = it->response_function();
    if (configOptions.expansionCoeffGradFlag)
      Teuchos::setCol(it->response_gradient(), (int)i, expansionCoeffGrads);
  }

  numCollocPts = new_colloc_pts;
}


void InterpPolyApproximation::
update_sparse_interpolation_basis(unsigned short max_level)
{
  SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
  const Real3DArray& gauss_pts_1d = ssg_driver->gauss_points_array();

  // resize if needed (leaving previous levels unmodified)
  size_t i, j, k, basis_size = polynomialBasis.size();
  // j range is 0:w inclusive; i range is 1:w+1 inclusive
  unsigned short num_levels = max_level + 1;
  if (num_levels > basis_size) {
    polynomialBasis.resize(num_levels);
    for (i=basis_size; i<num_levels; ++i)
      polynomialBasis[i].resize(numVars);
  }

  // fill gaps that may exist within any level
  for (i=0; i<num_levels; ++i) { // i -> 0:num_levels-1 -> 0:ssg_level
    for (j=0; j<numVars; ++j) {
      const RealArray& gauss_pts_1d_ij =    gauss_pts_1d[i][j];
      BasisPolynomial&   poly_basis_ij = polynomialBasis[i][j];
      if (poly_basis_ij.is_null() && !gauss_pts_1d_ij.empty()) {
	bool found = false;
	for (k=0; k<j; ++k)
	  if (gauss_pts_1d_ij == gauss_pts_1d[i][k] &&  // vector equality
	      !polynomialBasis[i][k].is_null())
	    { found = true; break; }
	if (found) // reuse previous instances via shared representations
	  poly_basis_ij = polynomialBasis[i][k]; // shared rep
	else { // instantiate new unique instances
	  poly_basis_ij = BasisPolynomial(LAGRANGE);
	  poly_basis_ij.interpolation_points(gauss_pts_1d_ij);
	}
      }
    }
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
  SparseGridDriver*   ssg_driver = (SparseGridDriver*)driverRep;
  const UShortArray&    sm_index = ssg_driver->smolyak_multi_index()[tp_index];
  const UShort2DArray&       key = ssg_driver->collocation_key()[tp_index];
  const SizetArray& colloc_index = ssg_driver->collocation_indices()[tp_index];

  tpValue = 0.;
  size_t i, j, num_colloc_pts = key.size();
  for (i=0; i<num_colloc_pts; ++i) {
    const UShortArray& key_i = key[i];
    Real L_i = 1.0;
    for (j=0; j<numVars; ++j)
      L_i *= polynomialBasis[sm_index[j]][j].get_value(x[j], key_i[j]);
    tpValue += expansionCoeffs[colloc_index[i]] * L_i;
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
  SparseGridDriver*   ssg_driver = (SparseGridDriver*)driverRep;
  const UShortArray&    sm_index = ssg_driver->smolyak_multi_index()[tp_index];
  const UShort2DArray&       key = ssg_driver->collocation_key()[tp_index];
  const SizetArray& colloc_index = ssg_driver->collocation_indices()[tp_index];

  if (tpGradient.length() != numVars)
    tpGradient.sizeUninitialized(numVars);
  tpGradient = 0.;
  size_t i, j, k, num_colloc_pts = key.size();
  for (i=0; i<num_colloc_pts; ++i) {
    const UShortArray& key_i   = key[i];
    const Real&        coeff_i = expansionCoeffs[colloc_index[i]];
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
tensor_product_gradient(const RealVector& x, const SizetArray& dvv)
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
			const SizetArray& dvv)
{
  SparseGridDriver*   ssg_driver = (SparseGridDriver*)driverRep;
  const UShortArray&    sm_index = ssg_driver->smolyak_multi_index()[tp_index];
  const UShort2DArray&       key = ssg_driver->collocation_key()[tp_index];
  const SizetArray& colloc_index = ssg_driver->collocation_indices()[tp_index];

  size_t num_deriv_vars = dvv.size();
  if (tpGradient.length() != num_deriv_vars)
    tpGradient.sizeUninitialized(num_deriv_vars);
  tpGradient = 0.;
  size_t i, j, k, deriv_index, num_colloc_pts = key.size();
  for (i=0; i<num_colloc_pts; ++i) {
    const UShortArray& key_i   = key[i];
    const Real&        coeff_i = expansionCoeffs[colloc_index[i]];
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
  const Real2DArray&   gauss_wts_1d = tpq_driver->gauss_weights_array();

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
  const UShortArray&     sm_index = ssg_driver->smolyak_multi_index()[tp_index];
  const UShort2DArray&        key = ssg_driver->collocation_key()[tp_index];
  const SizetArray&  colloc_index = ssg_driver->collocation_indices()[tp_index];
  const Real3DArray& gauss_wts_1d = ssg_driver->gauss_weights_array();

  tpMean = 0.;
  size_t i, j, num_colloc_pts = key.size();
  for (i=0; i<num_colloc_pts; ++i) {
    const UShortArray& key_i = key[i];
    Real prod_i = 1.;
    for (j=0; j<numVars; ++j)
      prod_i *= (randomVarsKey[j]) ? gauss_wts_1d[sm_index[j]][j][key_i[j]] :
	polynomialBasis[sm_index[j]][j].get_value(x[j], key_i[j]);
    tpMean += expansionCoeffs[colloc_index[i]] * prod_i;
  }
  return tpMean;
}


/** Overloaded version supporting tensor-product quadrature. */
const RealVector& InterpPolyApproximation::
tensor_product_mean_gradient(const RealVector& x, const SizetArray& dvv)
{
  TensorProductDriver* tpq_driver   = (TensorProductDriver*)driverRep;
  const UShort2DArray& key          = tpq_driver->collocation_key();
  const Real2DArray&   gauss_wts_1d = tpq_driver->gauss_weights_array();

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
    if (randomVarsKey[deriv_index] && !configOptions.expansionCoeffGradFlag) {
      PCerr << "Error: expansion coefficient gradients not defined in "
	    << "InterpPolyApproximation::get_mean_gradient()." << std::endl;
      abort_handler(-1);
    }
    else if (!randomVarsKey[deriv_index] && !configOptions.expansionCoeffFlag) {
      PCerr << "Error: expansion coefficients not defined in "
	    << "InterpPolyApproximation::get_mean_gradient()" << std::endl;
      abort_handler(-1);
    }
    for (j=0; j<numCollocPts; ++j) {
      const UShortArray& key_j = key[j];
      Real wt_prod_j = 1.;
      for (it=randomIndices.begin(); it!=randomIndices.end(); ++it)
	{ k = *it; wt_prod_j *= gauss_wts_1d[k][key_j[k]]; }
      if (randomVarsKey[deriv_index]) {
	// --------------------------------------------------------------------
	// derivative of All var expansion w.r.t. random var (design insertion)
	// --------------------------------------------------------------------
	Real Lsa_j = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it)
	  { k = *it; Lsa_j *= poly_basis_0[k].get_value(x[k], key_j[k]); }
	tpMeanGrad[i] += wt_prod_j * Lsa_j * expansionCoeffGrads(cntr, j);
      }
      else {
	// ---------------------------------------------------------------------
	// deriv of All var expansion w.r.t. nonrandom var (design augmentation)
	// ---------------------------------------------------------------------
	Real dLsa_j_dsa_i = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
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
			     const SizetArray& dvv)
{
  SparseGridDriver*    ssg_driver = (SparseGridDriver*)driverRep;
  const UShortArray&     sm_index = ssg_driver->smolyak_multi_index()[tp_index];
  const UShort2DArray&        key = ssg_driver->collocation_key()[tp_index];
  const SizetArray&  colloc_index = ssg_driver->collocation_indices()[tp_index];
  const Real3DArray& gauss_wts_1d = ssg_driver->gauss_weights_array();

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
    if (randomVarsKey[deriv_index] && !configOptions.expansionCoeffGradFlag) {
      PCerr << "Error: expansion coefficient gradients not defined in "
	    << "InterpPolyApproximation::get_mean_gradient()." << std::endl;
      abort_handler(-1);
    }
    else if (!randomVarsKey[deriv_index] && !configOptions.expansionCoeffFlag) {
      PCerr << "Error: expansion coefficients not defined in "
	    << "InterpPolyApproximation::get_mean_gradient()" << std::endl;
      abort_handler(-1);
    }
    for (j=0; j<num_colloc_pts; ++j) {
      const UShortArray& key_j = key[j];
      Real wt_prod_j = 1.;
      for (it=randomIndices.begin(); it!=randomIndices.end(); ++it) {
	k = *it;
	wt_prod_j *= gauss_wts_1d[sm_index[k]][k][key_j[k]];
      }
      if (randomVarsKey[deriv_index]) {
	// --------------------------------------------------------------------
	// derivative of All var expansion w.r.t. random var (design insertion)
	// --------------------------------------------------------------------
	Real Lsa_j = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	  k = *it;
	  Lsa_j *= polynomialBasis[sm_index[k]][k].get_value(x[k], key_j[k]);
	}
	tpMeanGrad[i]
	  += wt_prod_j * Lsa_j * expansionCoeffGrads(cntr, colloc_index[j]);
      }
      else {
	// ---------------------------------------------------------------------
	// deriv of All var expansion w.r.t. nonrandom var (design augmentation)
	// ---------------------------------------------------------------------
	Real dLsa_j_dsa_i = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	  k = *it;
	  dLsa_j_dsa_i *= (k == deriv_index) ?
	    polynomialBasis[sm_index[k]][k].get_gradient(x[k], key_j[k]) :
	    polynomialBasis[sm_index[k]][k].get_value(x[k],    key_j[k]);
	}
	tpMeanGrad[i]
	  += expansionCoeffs[colloc_index[j]] * wt_prod_j * dLsa_j_dsa_i;
      }
    }
    if (randomVarsKey[deriv_index]) // deriv w.r.t. design var insertion
      cntr++;
  }

  return tpMeanGrad;
}


/** Overloaded version supporting tensor-product quadrature. */
const Real& InterpPolyApproximation::
tensor_product_covariance(const RealVector& x, const RealVector& exp_coeffs_2)
{
  TensorProductDriver* tpq_driver   = (TensorProductDriver*)driverRep;
  const UShort2DArray& key          = tpq_driver->collocation_key();
  const Real2DArray&   gauss_wts_1d = tpq_driver->gauss_weights_array();

  tpVariance = 0.;
  size_t i, j, k;
  Real mean1 = 0., mean2 = 0.;
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
    Real wt_Ls_prod_i = wt_prod_i * Ls_prod_i;
    mean1 += exp_coeff_i     * wt_Ls_prod_i;
    mean2 += exp_coeffs_2[i] * wt_Ls_prod_i;
    for (j=0; j<numCollocPts; ++j) {
      const UShortArray& key_j = key[j];
      // to include the ij-th term, colloc pts xi_i must be the same as xi_j
      // for the ran var subset.  In this case, wt_prod_i may be reused.
      bool include = true;
      for (it=randomIndices.begin(); it!=randomIndices.end(); ++it)
	if (key_i[*it] != key_j[*it])
	  { include = false; break; }
      if (include) {
	Real Ls_prod_j = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	  k = *it;
	  Ls_prod_j *= poly_basis_0[k].get_value(x[k], key_j[k]);
	}
	tpVariance += wt_prod_i * Ls_prod_i * Ls_prod_j * exp_coeff_i
	  * exp_coeffs_2[j];
      }
    }
  }
  tpVariance -= mean1*mean2;
  return tpVariance;
}


/** Overloaded version supporting Smolyak sparse grids. */
const Real& InterpPolyApproximation::
tensor_product_covariance(const RealVector& x, const RealVector& exp_coeffs_2,
			  size_t tp_index)
{
  SparseGridDriver*    ssg_driver = (SparseGridDriver*)driverRep;
  const UShortArray&     sm_index = ssg_driver->smolyak_multi_index()[tp_index];
  const UShort2DArray&        key = ssg_driver->collocation_key()[tp_index];
  const SizetArray&  colloc_index = ssg_driver->collocation_indices()[tp_index];
  const Real3DArray& gauss_wts_1d = ssg_driver->gauss_weights_array();

  size_t i, j, k, index, num_colloc_pts = key.size();
  Real mean1 = 0., mean2 = 0.;
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
    index = colloc_index[i];
    const Real& exp_coeff_i = expansionCoeffs[index];
    Real wt_Ls_prod_i = wt_prod_i * Ls_prod_i;
    mean1 += exp_coeff_i         * wt_Ls_prod_i;
    mean2 += exp_coeffs_2[index] * wt_Ls_prod_i;
    for (j=0; j<num_colloc_pts; ++j) {
      const UShortArray& key_j = key[j];
      // to include the ij-th term, colloc pts xi_i must be the same as xi_j
      // for the ran var subset.  In this case, wt_prod_i may be reused.
      bool include = true;
      for (it=randomIndices.begin(); it!=randomIndices.end(); ++it)
	if (key_i[*it] != key_j[*it])
	  { include = false; break; }
      if (include) {
	Real Ls_prod_j = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	  k = *it;
	  Ls_prod_j *= polynomialBasis[sm_index[k]][k].get_value(x[k],key_j[k]);
	}
	tpVariance += wt_prod_i * Ls_prod_i * Ls_prod_j * exp_coeff_i
	  * exp_coeffs_2[colloc_index[j]];
      }
    }
  }
  tpVariance -= mean1*mean2;
  return tpVariance;
}


/** Overloaded version supporting tensor-product quadrature. */
const RealVector& InterpPolyApproximation::
tensor_product_variance_gradient(const RealVector& x, const SizetArray& dvv)
{
  TensorProductDriver* tpq_driver   = (TensorProductDriver*)driverRep;
  const UShort2DArray& key          = tpq_driver->collocation_key();
  const Real2DArray&   gauss_wts_1d = tpq_driver->gauss_weights_array();

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
    if (randomVarsKey[deriv_index] && !configOptions.expansionCoeffGradFlag) {
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
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
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
	for (it=randomIndices.begin(); it!=randomIndices.end(); ++it)
	  if (key_j[*it] != key_k[*it])
	    { include = false; break; }
	if (include) {
	  Real Lsa_k = 1.;
	  for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
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
	    for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it){
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
				 const SizetArray& dvv)
{
  SparseGridDriver*    ssg_driver = (SparseGridDriver*)driverRep;
  const UShortArray&     sm_index = ssg_driver->smolyak_multi_index()[tp_index];
  const UShort2DArray&        key = ssg_driver->collocation_key()[tp_index];
  const SizetArray&  colloc_index = ssg_driver->collocation_indices()[tp_index];
  const Real3DArray& gauss_wts_1d = ssg_driver->gauss_weights_array();

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
    if (randomVarsKey[deriv_index] && !configOptions.expansionCoeffGradFlag) {
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
	mean += expansionCoeffs[colloc_index[j]] * wt_prod_j * Lsa_j;
      if (randomVarsKey[deriv_index])
	mean_grad_i
	  += wt_prod_j * Lsa_j * expansionCoeffGrads(cntr,colloc_index[j]);
      else {
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	  k = *it;
	  dLsa_j_dsa_i *= (k == deriv_index) ?
	    polynomialBasis[sm_index[k]][k].get_gradient(x[k], key_j[k]) :
	    polynomialBasis[sm_index[k]][k].get_value(x[k],    key_j[k]);
	}
	mean_grad_i
	  += wt_prod_j * dLsa_j_dsa_i * expansionCoeffs[colloc_index[j]];
      }
      // second loop of double sum
      for (k=0; k<num_colloc_pts; ++k) {
	const UShortArray& key_k = key[k];
	// to include the jk-th term, colloc pts xi_j must be the same as xi_k
	// for the random var subset.  In this case, wt_prod_j may be reused.
	bool include = true;
	for (it=randomIndices.begin(); it!=randomIndices.end(); ++it)
	  if (key_j[*it] != key_k[*it])
	    { include = false; break; }
	if (include) {
	  Real Lsa_k = 1.;
	  for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	    l = *it;
	    Lsa_k *= polynomialBasis[sm_index[l]][l].get_value(x[l], key_k[l]);
	  }
	  if (randomVarsKey[deriv_index])
	    // ---------------------------------------------------------------
	    // deriv of All var expansion w.r.t. random var (design insertion)
	    // ---------------------------------------------------------------
	    tpVarianceGrad[i] += wt_prod_j * Lsa_j * Lsa_k *
	      (expansionCoeffs[colloc_index[j]] *
	       expansionCoeffGrads(cntr,colloc_index[k]) +
	       expansionCoeffGrads(cntr,colloc_index[j]) *
	       expansionCoeffs[colloc_index[k]]);
	  else {
	    // ---------------------------------------------------------------
	    // deriv of All var exp w.r.t. nonrandom var (design augmentation)
	    // ---------------------------------------------------------------
	    Real dLsa_k_dsa_i = 1.;
	    for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it){
	      l = *it;
	      dLsa_k_dsa_i *= (l == deriv_index) ?
		polynomialBasis[sm_index[l]][l].get_gradient(x[l], key_k[l]):
		polynomialBasis[sm_index[l]][l].get_value(x[l],    key_k[l]);
	    }
	    tpVarianceGrad[i] +=
	      wt_prod_j * expansionCoeffs[colloc_index[j]] *
	      expansionCoeffs[colloc_index[k]] *
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
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "InterpPolyApproximation::get_value()" << std::endl;
    abort_handler(-1);
  }

  // sum expansion to get response prediction
  switch (configOptions.expCoeffsSolnApproach) {
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
      if (sm_coeffs[i])
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
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "InterpPolyApproximation::get_gradient()" << std::endl;
    abort_handler(-1);
  }

  switch (configOptions.expCoeffsSolnApproach) {
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
      if (coeff) {
	const RealVector& tp_grad = tensor_product_gradient(x, i);
	for (j=0; j<numVars; ++j)
	  approxGradient[j] += coeff * tp_grad[j];
      }
    }
    return approxGradient;
    break;
  }
  }
}


const RealVector& InterpPolyApproximation::
get_gradient(const RealVector& x, const SizetArray& dvv)
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "InterpPolyApproximation::get_gradient()" << std::endl;
    abort_handler(-1);
  }

  switch (configOptions.expCoeffsSolnApproach) {
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
      if (coeff) {
	const RealVector& tp_grad = tensor_product_gradient(x, i, dvv);
	for (j=0; j<num_deriv_vars; ++j)
	  approxGradient[j] += coeff * tp_grad[j];
      }
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
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "InterpPolyApproximation::get_mean()" << std::endl;
    abort_handler(-1);
  }

  const RealVector& wt_sets = driverRep->weight_sets();
  Real& mean = centralNumMoments[0]; mean = 0.;
  for (size_t i=0; i<numCollocPts; ++i)
    mean += expansionCoeffs[i] * wt_sets[i];
  return mean;
}


/** In this case, a subset of the expansion variables are random
    variables and the mean of the expansion involves evaluating the
    expectation over this subset. */
const Real& InterpPolyApproximation::get_mean(const RealVector& x)
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "InterpPolyApproximation::get_mean()" << std::endl;
    abort_handler(-1);
  }

  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE:
    return tensor_product_mean(x);
    break;
  case SPARSE_GRID: {
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    const IntArray&   sm_coeffs  = ssg_driver->smolyak_coefficients();
    size_t i, num_smolyak_indices = sm_coeffs.size();
    Real& mean = centralNumMoments[0]; mean = 0.;
    for (i=0; i<num_smolyak_indices; ++i)
      if (sm_coeffs[i])
	mean += sm_coeffs[i] * tensor_product_mean(x, i);
    return mean;
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
  if (!configOptions.expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in "
	  << "InterpPolyApproximation::get_mean_gradient()." << std::endl;
    abort_handler(-1);
  }

  const RealVector& wt_sets = driverRep->weight_sets();
  size_t i, j, num_deriv_vars = expansionCoeffGrads.numRows();
  if (meanGradient.length() != num_deriv_vars)
    meanGradient.sizeUninitialized(num_deriv_vars);
  meanGradient = 0.;
  for (i=0; i<numCollocPts; ++i) {
    const Real& wt_prod_i = wt_sets[i];
    for (j=0; j<num_deriv_vars; ++j)
      meanGradient[j] += expansionCoeffGrads(j,i) * wt_prod_i;
  }
  return meanGradient;
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
get_mean_gradient(const RealVector& x, const SizetArray& dvv)
{
  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE:
    return tensor_product_mean_gradient(x, dvv);
    break;
  case SPARSE_GRID:
    size_t num_deriv_vars = dvv.size();
    if (meanGradient.length() != num_deriv_vars)
      meanGradient.sizeUninitialized(num_deriv_vars);
    meanGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    const IntArray&   sm_coeffs  = ssg_driver->smolyak_coefficients();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff = sm_coeffs[i];
      if (coeff) {
	const RealVector& tpm_grad = tensor_product_mean_gradient(x, i, dvv);
	for (j=0; j<num_deriv_vars; ++j)
	  meanGradient[j] += coeff * tpm_grad[j];
      }
    }
    return meanGradient;
    break;
  }
}


/** In this case, all expansion variables are random variables and the
    variance of the expansion is the sum over all but the first term
    of the coefficients squared times the polynomial norms squared. */
const Real& InterpPolyApproximation::get_variance()
{
  centralNumMoments[1] = get_covariance(expansionCoeffs);
  return centralNumMoments[1];
}


/** In this case, a subset of the expansion variables are random
    variables and the variance of the expansion involves summations
    over this subset. */
const Real& InterpPolyApproximation::get_variance(const RealVector& x)
{
  centralNumMoments[1] = get_covariance(x, expansionCoeffs);
  return centralNumMoments[1];
}


Real InterpPolyApproximation::get_covariance(const RealVector& exp_coeffs_2)
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "InterpPolyApproximation::get_covariance()" << std::endl;
    abort_handler(-1);
  }

  Real var = 0., mean_1 = 0., mean_2 = 0.;
  const RealVector& wt_sets = driverRep->weight_sets();
  for (size_t i=0; i<numCollocPts; ++i) {
    const Real& coeff_2i  = exp_coeffs_2[i];
    const Real& wt_prod_i = wt_sets[i];
    Real coeff_wt_1i = expansionCoeffs[i] * wt_prod_i;
    mean_1 += coeff_wt_1i;
    mean_2 += coeff_2i * wt_prod_i;
    var    += coeff_wt_1i * coeff_2i;
  }
  var -= mean_1*mean_2;
  return var;
}


/** In this case, a subset of the expansion variables are random
    variables and the variance of the expansion involves summations
    over this subset. */
Real InterpPolyApproximation::
get_covariance(const RealVector& x, const RealVector& exp_coeffs_2)
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "InterpPolyApproximation::get_covariance()" << std::endl;
    abort_handler(-1);
  }

  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE:
    return tensor_product_covariance(x, exp_coeffs_2);
    break;
  case SPARSE_GRID:
    Real var = 0.;
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    const IntArray&   sm_coeffs  = ssg_driver->smolyak_coefficients();
    size_t i, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i)
      if (sm_coeffs[i])
	var += sm_coeffs[i] * tensor_product_covariance(x, exp_coeffs_2, i);
    return var;
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
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "InterpPolyApproximation::get_variance_gradient()" << std::endl;
    abort_handler(-1);
  }
  if (!configOptions.expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in "
	  << "InterpPolyApproximation::get_variance_gradient()." << std::endl;
    abort_handler(-1);
  }

  const RealVector& wt_sets = driverRep->weight_sets();
  size_t i, j, num_deriv_vars = expansionCoeffGrads.numRows();
  if (varianceGradient.length() != num_deriv_vars)
    varianceGradient.sizeUninitialized(num_deriv_vars);
  varianceGradient = 0.;

  Real mean = 0.;
  for (i=0; i<numCollocPts; ++i)
    mean += expansionCoeffs[i] * wt_sets[i];
  for (i=0; i<numCollocPts; ++i) {
    Real term_i = 2. * (expansionCoeffs[i] - mean) * wt_sets[i];
    for (j=0; j<num_deriv_vars; ++j)
      varianceGradient[j] += term_i * expansionCoeffGrads(j,i);
  }

  return varianceGradient;
}


/** In this function, a subset of the expansion variables are random
    variables and any augmented design/state variables (i.e., not
    inserted as random variable distribution parameters) are included
    in the expansion.  This function must handle the mixed case, where
    some design/state variables are augmented (and are part of the
    expansion) and some are inserted (derivatives are obtained from
    expansionCoeffGrads). */
const RealVector& InterpPolyApproximation::
get_variance_gradient(const RealVector& x, const SizetArray& dvv)
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "InterpPolyApproximation::get_variance_gradient()" << std::endl;
    abort_handler(-1);
  }

  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE:
    return tensor_product_variance_gradient(x, dvv);
    break;
  case SPARSE_GRID:
    size_t num_deriv_vars = dvv.size();
    if (varianceGradient.length() != num_deriv_vars)
      varianceGradient.sizeUninitialized(num_deriv_vars);
    varianceGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    const IntArray&   sm_coeffs  = ssg_driver->smolyak_coefficients();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff = sm_coeffs[i];
      if (coeff) {
	const RealVector& tpv_grad
	  = tensor_product_variance_gradient(x, i, dvv);
	for (j=0; j<num_deriv_vars; ++j)
	  varianceGradient[j] += coeff * tpv_grad[j];
      }
    }
    return varianceGradient;
    break;
  }
}


void InterpPolyApproximation::compute_component_effects()
{
  // perform subset sort
  constituentSets.resize(sobolIndices.length());
  get_subsets();

  // initialize partialVariance
  if (partialVariance.empty())
    partialVariance.sizeUninitialized(sobolIndices.length());
  partialVariance = 0.;

  Real mean = get_mean(), total_variance = get_variance();
  partialVariance[0] = mean*mean; // init with mean sq

  // Solve for partial variance
  for (IntIntMIter map_iter=sobolIndexMap.begin();
       map_iter!=sobolIndexMap.end(); ++map_iter) {
    // partialVariance[0] stores the mean; it is not a component function
    // and does not follow the procedures for obtaining variance 
    if (map_iter->first) {
      partial_variance(map_iter->first);
      sobolIndices[map_iter->second]
	= partialVariance[map_iter->second]/total_variance;
      // total indices simply identify the membership of the sobolIndices 
      // and adds it to the appropriate bin
    }
  }
#ifdef DEBUG
  PCout << "In InterpPolyApproximation::compute_component_effects(), "
	<< "sobolIndices =\n"; write_data(PCout, sobolIndices);
#endif // DEBUG
}


void InterpPolyApproximation::compute_total_effects()
{
  // iterate through existing indices if all component indices are available
  totalSobolIndices = 0; // init total indices
  if (configOptions.vbdControl == ALL_VBD) {
    for (IntIntMIter itr=sobolIndexMap.begin(); itr!=sobolIndexMap.end(); ++itr)
      for (int k=0; k<numVars; ++k) {
        if (itr->first & (1 << k))
          totalSobolIndices[k] += sobolIndices[itr->second];
        totalSobolIndices[k] = std::abs(totalSobolIndices[k]);
      }
  }
  // if not available, compute total indices independently
  // approach parallels partial_variance_integral where the algorithm is 
  // separated by integration approach
  else {
    Real total_variance = get_variance();
    int j, set_value;
    switch (configOptions.expCoeffsSolnApproach) {
      case QUADRATURE: {
        for (j=0; j<numVars; ++j) {
          // define set_value that includes all but index of interest
          set_value = (int)std::pow(2.,int(numVars)) - (int)std::pow(2.,j) - 1;
          totalSobolIndices[j]
	    = std::abs(1 - (total_effects_integral(set_value)/total_variance));
        }
        break;
      }
      case SPARSE_GRID: {
        SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
        const IntArray&   sm_coeffs  = ssg_driver->smolyak_coefficients();
        // Smolyak recursion of anisotropic tensor products
        size_t i, num_smolyak_indices = sm_coeffs.size();
        // iterate each variable 
        for (j=0; j<numVars; ++j) {
          set_value = (int)std::pow(2.,int(numVars)) - (int)std::pow(2.,j) - 1; 
          for (i=0; i<num_smolyak_indices; ++i)
	    if (sm_coeffs[i])
	      totalSobolIndices[j]
		+= sm_coeffs[i]*total_effects_integral(set_value,i);
          totalSobolIndices[j]
	    = std::abs(1 - (totalSobolIndices[j]/total_variance));
        }
        break;
      }
    }
  }
#ifdef DEBUG
  PCout << "In InterpPolyApproximation::compute_total_effects(), "
	<< "totalSobolIndices =\n"; write_data(PCout, totalSobolIndices);
#endif // DEBUG
}


Real InterpPolyApproximation::total_effects_integral(int set_value)
{
  // Some other routine here
  TensorProductDriver* tpq_driver   = (TensorProductDriver*)driverRep;
  const UShort2DArray& key          = tpq_driver->collocation_key();
  const Real2DArray&   gauss_wts_1d = tpq_driver->gauss_weights_array();
  const UShortArray&   quad_order   = tpq_driver->quadrature_order();;

  // Distinguish between non-members and members of the given set, set_value
  BoolDeque nonmember_vars(numVars,true); 
  int num_mem_expansion_coeffs = 1; // number of expansion coeffs in
                                    // member-variable-only expansion 
  IntVector indexing_factor; // factors indexing member variables 
  indexing_factor.sizeUninitialized(numVars);

  indexing_factor = 1;
        
  for (int k=0; k<numVars; ++k) {
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
  Real total_mean = get_mean();
  Real integral = 0;
  for (i=0; i<num_mem_expansion_coeffs; ++i)
    integral
      += std::pow((mem_expansion_coeffs[i]-total_mean),2.0)*mem_weights[i];
  return integral;
}


Real InterpPolyApproximation::
total_effects_integral(int set_value, size_t tp_index)
{
  SparseGridDriver*    ssg_driver = (SparseGridDriver*)driverRep;
  const UShortArray&     sm_index = ssg_driver->smolyak_multi_index()[tp_index];
  const UShort2DArray&        key = ssg_driver->collocation_key()[tp_index];
  const SizetArray&  colloc_index = ssg_driver->collocation_indices()[tp_index];
  const Real3DArray& gauss_wts_1d = ssg_driver->gauss_weights_array();

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
  for (int k=0; k<numVars; ++k) {
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
      += expansionCoeffs[colloc_index[i]]*prod_i_nonmembers;
  }
 
  // Now integrate over the remaining variables	
  Real integral = 0;
  Real total_mean = get_mean();
  for (i=0; i<num_mem_expansion_coeffs; ++i)
    integral
      += std::pow(mem_expansion_coeffs[i] - total_mean, 2.0)*mem_weights[i];
  return integral;
}


/** Find constituent subsets. */
void InterpPolyApproximation::get_subsets()
{
  // includes the "zero" set
  //int num_subsets = sobolIndices.length(); 

  // Here we want to utilize the integer representation of the subset
  // but we want to store it in a size appropriate container
  // so finding lower sets is given the argument of the integer rep (->first)
  // and stored in constituentSets in size-appropriate-index-map (->second)
  for (IntIntMIter map_iter=sobolIndexMap.begin();
       map_iter!=sobolIndexMap.end(); ++map_iter) {
    lower_sets(map_iter->first, constituentSets[map_iter->second]);
    constituentSets[map_iter->second].erase(map_iter->first);
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
  for (int k=0; k<numVars; ++k)
    // this performs a bitwise comparison by shifting 1 by k spaces 
    // and comparing that to a binary form of plus_one_set; this allows 
    // the variable membership using integers instead of a d-array of bools
    if (plus_one_set & (1 << k)) 
      // if subset i contains variable k, remove that variable from the set 
      // by converting the bit-form of (1<<k) to an integer and subtract from
      // the plus_one_set
      lower_sets(plus_one_set-(int)std::pow(2.0,k),top_level_set);
}


/** Forms an interpolant over variables that are members of the given set.
    Finds the variance of the interpolant w.r.t. the variables in the set.
    Overloaded version supporting tensor-product quadrature. */
Real InterpPolyApproximation::partial_variance_integral(int set_value)
{
  TensorProductDriver* tpq_driver   = (TensorProductDriver*)driverRep;
  const UShort2DArray& key          = tpq_driver->collocation_key();
  const Real2DArray&   gauss_wts_1d = tpq_driver->gauss_weights_array();
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
	
  for (int k=0; k<numVars; ++k) {
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
  const UShortArray&     sm_index = ssg_driver->smolyak_multi_index()[tp_index];
  const UShort2DArray&        key = ssg_driver->collocation_key()[tp_index];
  const SizetArray&  colloc_index = ssg_driver->collocation_indices()[tp_index];
  const Real3DArray& gauss_wts_1d = ssg_driver->gauss_weights_array();

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
  for (int k=0; k<numVars; ++k) {
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
      += expansionCoeffs[colloc_index[i]]*prod_i_nonmembers;
  }

  // Now integrate over the remaining variables	
  Real integral = 0;
  for (i=0; i<num_mem_expansion_coeffs; ++i)
    integral += std::pow(mem_expansion_coeffs[i], 2.0)*mem_weights[i];
  return integral;	
}


/** Computes the variance of component functions. Assumes that all
    subsets of set_value have been computed in advance which will be
    true so long as the partial_variance is called following
    appropriate enumeration of set value  */
void InterpPolyApproximation::partial_variance(int set_value)
{
  // Computes the integral first
  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE:
    partialVariance[sobolIndexMap[set_value]]
      = partial_variance_integral(set_value);
    break;
  case SPARSE_GRID: {
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    const IntArray&   sm_coeffs  = ssg_driver->smolyak_coefficients();
    // Smolyak recursion of anisotropic tensor products
    size_t num_smolyak_indices = sm_coeffs.size();
    for (size_t i=0; i<num_smolyak_indices; ++i)
      if (sm_coeffs[i])
	partialVariance[sobolIndexMap[set_value]] += sm_coeffs[i]
	  * partial_variance_integral(set_value, i);
    break;
  }
  }
  
  // Now subtract the contributions from constituent subsets
  IntSet::iterator itr;
  for (itr  = constituentSets[sobolIndexMap[set_value]].begin();
       itr != constituentSets[sobolIndexMap[set_value]].end(); ++itr) 
    partialVariance[sobolIndexMap[set_value]]
      -= partialVariance[sobolIndexMap[*itr]];
}

} // namespace Pecos

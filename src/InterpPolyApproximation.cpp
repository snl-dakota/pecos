/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
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


/** This version provides the polynomial types needed to retrieve
    collocation points and weights by an integration driver.  These
    may involve orthogonal polynomials which will differ from the
    interpolation polynomial types used in the basis. */
bool InterpPolyApproximation::
initialize_integration_basis_types(const ShortArray& u_types,
				   const BasisConfigOptions& bc_options,
				   ShortArray& basis_types)
{
  bool extra_dist_params = false;

  // Initialize basis_types and extra_dist_params from u_types.
  size_t i, num_vars = u_types.size();
  if (basis_types.size() != num_vars)
    basis_types.resize(num_vars);
  for (i=0; i<num_vars; ++i) {
    switch (u_types[i]) {
    case STD_NORMAL:
      basis_types[i] = HERMITE_ORTHOG;                                break;
    case STD_UNIFORM:
      if (bc_options.piecewiseBasis)
	basis_types[i] = (bc_options.useDerivs) ? PIECEWISE_CUBIC_INTERP :
	  PIECEWISE_LINEAR_INTERP;
      else
	basis_types[i] = (bc_options.useDerivs) ? HERMITE_INTERP :
	  LEGENDRE_ORTHOG;
      break;
    case STD_EXPONENTIAL:
      basis_types[i] = LAGUERRE_ORTHOG;                               break;
    case STD_BETA:
      basis_types[i] = JACOBI_ORTHOG;       extra_dist_params = true; break;
    case STD_GAMMA:
      basis_types[i] = GEN_LAGUERRE_ORTHOG; extra_dist_params = true; break;
    default:
      basis_types[i] = NUM_GEN_ORTHOG;      extra_dist_params = true; break;
    }
  }

  return extra_dist_params;
}


void InterpPolyApproximation::
initialize_polynomial_basis_type(short& poly_type_1d, short& rule)
{
  switch (basisType) {
  case PIECEWISE_NODAL_INTERPOLATION_POLYNOMIAL:
  case PIECEWISE_HIERARCHICAL_INTERPOLATION_POLYNOMIAL:
    poly_type_1d = (basisConfigOptions.useDerivs) ?
      PIECEWISE_CUBIC_INTERP : PIECEWISE_LINEAR_INTERP;
    rule = NEWTON_COTES;                    break;
  case GLOBAL_NODAL_INTERPOLATION_POLYNOMIAL:
  case GLOBAL_HIERARCHICAL_INTERPOLATION_POLYNOMIAL:
    poly_type_1d = (basisConfigOptions.useDerivs) ?
      HERMITE_INTERP : LAGRANGE_INTERP;
    rule = NO_RULE;                         break;
  default:
    poly_type_1d = NO_POLY; rule = NO_RULE; break;
  }
}


int InterpPolyApproximation::min_coefficients() const
{
  // return the minimum number of data instances required to build the 
  // surface of multiple dimensions
  return (expConfigOptions.expansionCoeffFlag ||
	  expConfigOptions.expansionCoeffGradFlag) ? 1 : 0;
}


void InterpPolyApproximation::allocate_arrays()
{
  allocate_component_effects();
  allocate_total_effects();
  allocate_expansion_coefficients();

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    const UShortArray&   quad_order = tpq_driver->quadrature_order();

    // verify total number of collocation pts (should not include anchor pt)
    size_t i, j, num_colloc_pts = 1;
    for (i=0; i<numVars; ++i)
      num_colloc_pts *= quad_order[i];
    if (num_colloc_pts != numCollocPts) {
      PCerr << "Error: inconsistent total collocation point count in "
	    << "InterpPolyApproximation::allocate_arrays()" << std::endl;
      abort_handler(-1);
    }

    //bool update_exp_form = (quad_order != quadOrderPrev);
    // *** TO DO: capture updates to parameterized/numerical polynomials
    //if (update_exp_form) {
    //}

    // can't use quad_order > quadOrderPrev logic since only 1 pt set is stored
    bool update_basis_form = (quad_order != quadOrderPrev);
    if (update_basis_form)
      update_tensor_interpolation_basis();

    quadOrderPrev = quad_order;
    break;
  }
  case COMBINED_SPARSE_GRID: case HIERARCHICAL_SPARSE_GRID: {
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    unsigned short    ssg_level  = ssg_driver->level();
    const RealVector& aniso_wts  = ssg_driver->anisotropic_weights();

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
  if (!expConfigOptions.expansionCoeffFlag &&
      !expConfigOptions.expansionCoeffGradFlag) {
    PCerr << "Warning: neither expansion coefficients nor expansion "
	  << "coefficient gradients\n         are active in "
	  << "InterpPolyApproximation::compute_coefficients().\n         "
	  << "Bypassing approximation construction." << std::endl;
    return;
  }

  // For testing of anchor point logic:
  //size_t index = surrData.size() - 1;
  //surrData.anchor_point(surrData.variables_data()[index],
  //                      surrData.response_data()[index]);
  //surrData.pop(1);

  numCollocPts = surrData.size();
  // anchor point, if present, is the first expansionSample.
  if (surrData.anchor())
    ++numCollocPts;
  if (!numCollocPts) {
    PCerr << "Error: nonzero number of sample points required in "
	  << "InterpPolyApproximation::compute_coefficients()." << std::endl;
    abort_handler(-1);
  }

  allocate_arrays();
  compute_expansion_coefficients();

#ifdef INTERPOLATION_TEST
  // SC should accurately interpolate the collocation data for TPQ and
  // SSG with fully nested rules, but will exhibit interpolation error
  // for SSG with other rules.
  if (expConfigOptions.expansionCoeffFlag) {
    size_t i, index = 0, offset = (surrData.anchor()) ? 1 : 0,
      w7 = WRITE_PRECISION+7;
    Real interp_val, err, val_max_err = 0., grad_max_err = 0.,
      val_rmse = 0., grad_rmse = 0.;
    PCout << std::scientific << std::setprecision(WRITE_PRECISION);
    for (i=offset; i<numCollocPts; ++i, ++index) {
      const RealVector& c_vars = surrData.continuous_variables(index);
      const Real&      resp_fn = surrData.response_function(index);
      interp_val = value(c_vars);
      err = (std::abs(resp_fn) > DBL_MIN) ? std::abs(1. - interp_val/resp_fn) :
	                                    std::abs(resp_fn - interp_val);
      PCout << "Colloc pt " << std::setw(3) << i+1
	    << ": truth value  = "  << std::setw(w7) << resp_fn
	    << " interpolant = "    << std::setw(w7) << interp_val
	    << " relative error = " << std::setw(w7) << err <<'\n';
      if (err > val_max_err) val_max_err = err;
      val_rmse += err * err;
      if (basisConfigOptions.useDerivs) {
	const RealVector& resp_grad   = surrData.response_gradient(index);
	const RealVector& interp_grad = gradient_basis_variables(c_vars);
	for (size_t j=0; j<numVars; ++j) {
	  err = (std::abs(resp_grad[j]) > DBL_MIN) ?
	    std::abs(1. - interp_grad[j]/resp_grad[j]) :
	    std::abs(resp_grad[j] - interp_grad[j]);
	  PCout << "               " << "truth grad_" << j+1 << " = "
		<< std::setw(w7) << resp_grad[j]   << " interpolant = "
		<< std::setw(w7) << interp_grad[j] << " relative error = "
		<< std::setw(w7) << err << '\n';
	  if (err > grad_max_err) grad_max_err = err;
	  grad_rmse += err * err;
	}
      }
    }
    val_rmse = std::sqrt(val_rmse/(numCollocPts-offset));
    PCout << "\nValue interpolation errors:    " << std::setw(w7) << val_max_err
	  << " (max) " << std::setw(w7) << val_rmse << " (RMS)\n";
    if (basisConfigOptions.useDerivs) {
      grad_rmse = std::sqrt(grad_rmse/(numCollocPts-offset)/numVars);
      PCout << "Gradient interpolation errors: " << std::setw(w7)
	    << grad_max_err << " (max) " << std::setw(w7) << grad_rmse
	    << " (RMS)\n";
    }
  }
#endif // INTERPOLATION_TEST
}


void InterpPolyApproximation::increment_coefficients()
{
  bool err_flag = false;
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case COMBINED_SPARSE_GRID: case HIERARCHICAL_SPARSE_GRID: {
    // As for allocate_arrays(), increments are performed in coarser steps
    // than may be strictly necessary: all increments are filled in for all
    // vars for a step in level (ignoring anisotropy or generalized indices).
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    const UShortArray& trial_set = ssg_driver->trial_set();
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

  increment_expansion_coefficients();
}


void InterpPolyApproximation::decrement_coefficients()
{
  // leave polynomialBasis as is

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case COMBINED_SPARSE_GRID: case HIERARCHICAL_SPARSE_GRID: {
    // move previous expansion data to current expansion
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    savedLevMultiIndex.push_back(ssg_driver->trial_set());
    break;
  }
  }

  decrement_expansion_coefficients();

  numCollocPts = surrData.size(); // data already decremented
  if (surrData.anchor())
    ++numCollocPts;
}


void InterpPolyApproximation::restore_coefficients()
{
  // leave polynomialBasis as is (a previous increment is being restored)

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case COMBINED_SPARSE_GRID: case HIERARCHICAL_SPARSE_GRID: {
    // move previous expansion data to current expansion
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    std::deque<UShortArray>::iterator sit
      = std::find(savedLevMultiIndex.begin(), savedLevMultiIndex.end(),
		  ssg_driver->trial_set());
    if (sit != savedLevMultiIndex.end())
      savedLevMultiIndex.erase(sit);
    break;
  }
  }

  restore_expansion_coefficients();
}


void InterpPolyApproximation::finalize_coefficients()
{
  // leave polynomialBasis as is (all previous increments are being restored)

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case COMBINED_SPARSE_GRID: case HIERARCHICAL_SPARSE_GRID:
    // move previous expansion data to current expansion
    savedLevMultiIndex.clear();
    break;
  }

  finalize_expansion_coefficients();
}


void InterpPolyApproximation::update_tensor_interpolation_basis()
{
  TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
  const UShortArray&   quad_order = tpq_driver->quadrature_order();
  const UShortArray&    lev_index = tpq_driver->level_index();

  // resize if needed (leaving previous levels unmodified)
  size_t i, j, k, basis_size = polynomialBasis.size();
  unsigned short max_order = quad_order[0];
  for (i=1; i<numVars; ++i)
    if (quad_order[i] > max_order)
      max_order = quad_order[i];
  // quad_order range is 1:m; quad_index range is 0:m-1
  if (max_order > basis_size) {
    polynomialBasis.resize(max_order);
    for (i=basis_size; i<max_order; ++i)
      polynomialBasis[i].resize(numVars);
  }

  // fill any required gaps in polynomialBasis.
  const Real3DArray& colloc_pts_1d = driverRep->collocation_points_array();
  short poly_type_1d; short rule; bool found; unsigned short l_index;
  initialize_polynomial_basis_type(poly_type_1d, rule);
  for (j=0; j<numVars; ++j) {
    l_index = lev_index[j];
    const RealArray& colloc_pts_1d_ij          =   colloc_pts_1d[l_index][j];
    std::vector<BasisPolynomial>& poly_basis_i = polynomialBasis[l_index];
    BasisPolynomial& poly_basis_ij = poly_basis_i[j];
    if (poly_basis_ij.is_null()) { // does not account for parametric changes
                                   // resulting in new pts for existing orders
      found = false;
      for (k=0; k<numVars; ++k)
	if (k != j && !poly_basis_i[k].is_null() &&
	    colloc_pts_1d_ij == poly_basis_i[k].interpolation_points())
	  { found = true; break; }
      if (found) // reuse previous instance via shared representation
	poly_basis_ij = poly_basis_i[k]; // shared rep
      else { // instantiate and initialize a new unique instance
	poly_basis_ij = BasisPolynomial(poly_type_1d, rule);
	poly_basis_ij.interpolation_points(colloc_pts_1d_ij);
      }
    }
  }
}


void InterpPolyApproximation::
update_sparse_interpolation_basis(unsigned short max_level)
{
  // resize if needed (leaving previous levels unmodified)
  size_t i, j, k, basis_size = polynomialBasis.size();
  // j range is 0:w inclusive; i range is 1:w+1 inclusive
  unsigned short num_levels = max_level + 1;
  if (num_levels > basis_size) {
    polynomialBasis.resize(num_levels);
    for (i=basis_size; i<num_levels; ++i)
      polynomialBasis[i].resize(numVars);
  }

  // fill gaps that may exist within any level (SparseGridDriver::
  // update_1d_collocation_points_weights() updates in an unstructured manner)
  const Real3DArray& colloc_pts_1d = driverRep->collocation_points_array();
  short poly_type_1d; short rule; bool found;
  initialize_polynomial_basis_type(poly_type_1d, rule);
  for (i=0; i<num_levels; ++i) { // i -> 0:num_levels-1 -> 0:ssg_level
    const Real2DArray& colloc_pts_1d_i = colloc_pts_1d[i];
    std::vector<BasisPolynomial>& poly_basis_i = polynomialBasis[i];
    for (j=0; j<numVars; ++j) {
      const RealArray& colloc_pts_1d_ij = colloc_pts_1d_i[j];
      BasisPolynomial&    poly_basis_ij =    poly_basis_i[j];
      if ( poly_basis_ij.is_null() &&  // doesn't account for parametric changes
	  !colloc_pts_1d_ij.empty()) { // resulting in new pts for existing levs
	found = false;
	for (k=0; k<numVars; ++k)
	  if (k != j && !poly_basis_i[k].is_null() &&
	      colloc_pts_1d_ij == colloc_pts_1d_i[k])  // vector equality
	    { found = true; break; }
	if (found) // reuse previous instance via shared representation
	  poly_basis_ij = poly_basis_i[k]; // shared rep
	else { // instantiate and initialize a new unique instance
	  poly_basis_ij = BasisPolynomial(poly_type_1d, rule);
	  poly_basis_ij.interpolation_points(colloc_pts_1d_ij);
	}
      }
    }
  }
}


Real InterpPolyApproximation::
tensor_product_value(const RealVector& x, const RealVector& exp_t1_coeffs,
		     const RealMatrix& exp_t2_coeffs,
		     const UShortArray& basis_index, const UShort2DArray& key,
		     const SizetArray& colloc_index)
{
  Real tp_val = 0.; size_t i, num_colloc_pts = key.size();
  if (exp_t2_coeffs.empty()) {
    if (colloc_index.empty())
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += exp_t1_coeffs[i] * 
	          type1_interpolant_value(x, key[i], basis_index);
    else
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += exp_t1_coeffs[colloc_index[i]] *
	          type1_interpolant_value(x, key[i], basis_index);
  }
  else {
    size_t j, c_index;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      c_index = (colloc_index.empty()) ? i : colloc_index[i];
      tp_val += exp_t1_coeffs[c_index] *
	        type1_interpolant_value(x, key_i, basis_index);
      const Real* exp_t2_coeff_i = exp_t2_coeffs[c_index];
      for (j=0; j<numVars; ++j)
	tp_val += exp_t2_coeff_i[j] *
	          type2_interpolant_value(x, j, key_i, basis_index);
    }
  }
  return tp_val;
}


const RealVector& InterpPolyApproximation::
tensor_product_gradient_basis_variables(const RealVector& x,
					const RealVector& exp_t1_coeffs,
					const RealMatrix& exp_t2_coeffs,
					const UShortArray& basis_index,
					const UShort2DArray& key,
					const SizetArray& colloc_index)
{
  if (tpGradient.length() != numVars)
    tpGradient.sizeUninitialized(numVars);
  tpGradient = 0.;
  size_t i, j, num_colloc_pts = key.size();
  if (exp_t2_coeffs.empty()) {
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      const Real& exp_t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      for (j=0; j<numVars; ++j)
	tpGradient[j] += exp_t1_coeff_i *
	  type1_interpolant_gradient(x, j, key_i, basis_index);
    }
  }
  else {
    size_t k;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      const Real& exp_t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      const Real* exp_t2_coeff_i = (colloc_index.empty()) ?
	exp_t2_coeffs[i] : exp_t2_coeffs[colloc_index[i]];
      for (j=0; j<numVars; ++j) { // ith contribution to jth grad component
	tpGradient[j] += exp_t1_coeff_i *
	  type1_interpolant_gradient(x, j, key_i, basis_index);
	for (k=0; k<numVars; ++k) // type2 interpolant for kth grad comp
	  tpGradient[j] += exp_t2_coeff_i[k] *
	    type2_interpolant_gradient(x, j, k, key_i, basis_index);
      }
    }
  }
  return tpGradient;
}


const RealVector& InterpPolyApproximation::
tensor_product_gradient_basis_variables(const RealVector& x,
					const RealVector& exp_t1_coeffs,
					const RealMatrix& exp_t2_coeffs,
					const UShortArray& basis_index,
					const UShort2DArray& key,
					const SizetArray& colloc_index,
					const SizetArray& dvv)
{
  size_t i, j, deriv_index, num_colloc_pts = key.size(),
    num_deriv_vars = dvv.size();
  if (tpGradient.length() != num_deriv_vars)
    tpGradient.sizeUninitialized(num_deriv_vars);
  tpGradient = 0.;
  if (exp_t2_coeffs.empty()) {
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i   = key[i];
      const Real& exp_t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      for (j=0; j<num_deriv_vars; ++j) {
	deriv_index = dvv[j] - 1; // requires an "All" view
	tpGradient[j] += exp_t1_coeff_i *
	  type1_interpolant_gradient(x, deriv_index, key_i, basis_index);
      }
    }
  }
  else {
    size_t k;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      const Real& exp_t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      const Real* exp_t2_coeff_i = (colloc_index.empty()) ?
	exp_t2_coeffs[i] : exp_t2_coeffs[colloc_index[i]];
      for (j=0; j<num_deriv_vars; ++j) {
	deriv_index = dvv[j] - 1; // requires an "All" view
	tpGradient[j] += exp_t1_coeff_i *
	  type1_interpolant_gradient(x, deriv_index, key_i, basis_index);
	for (k=0; k<numVars; ++k)
	  tpGradient[j] += exp_t2_coeff_i[k] *
	    type2_interpolant_gradient(x, deriv_index, k, key_i, basis_index);
      }
    }
  }
  return tpGradient;
}


const RealVector& InterpPolyApproximation::
tensor_product_gradient_nonbasis_variables(const RealVector& x,
					   const RealMatrix& exp_t1_coeff_grads,
					   const UShortArray& basis_index,
					   const UShort2DArray& key,
					   const SizetArray& colloc_index)
{
  size_t i, j, num_colloc_pts = key.size(),
    num_deriv_vars = exp_t1_coeff_grads.numRows();
  if (tpGradient.length() != num_deriv_vars)
    tpGradient.sizeUninitialized(num_deriv_vars);
  tpGradient = 0.;
  for (i=0; i<num_colloc_pts; ++i) {
    const Real* exp_t1_coeff_grad_i = (colloc_index.empty()) ?
      exp_t1_coeff_grads[i] : exp_t1_coeff_grads[colloc_index[i]];
    Real t1_val = type1_interpolant_value(x, key[i], basis_index);
    for (j=0; j<numVars; ++j)
      tpGradient[j] += exp_t1_coeff_grad_i[j] * t1_val;
  }
  return tpGradient;
}


void InterpPolyApproximation::compute_component_effects()
{
  // perform subset sort
  size_t sobol_len = sobolIndices.length();
  constituentSets.resize(sobol_len);
  get_subsets();

  // initialize partialVariance
  if (partialVariance.empty())
    partialVariance.sizeUninitialized(sobol_len);
  partialVariance = 0.;

  // initialize with mean^2
  const Real& nm0 = numericalMoments[0]; // standardized, if not num exception
  partialVariance[0] = nm0*nm0;

  const Real& nm1 = numericalMoments[1]; // standardized, if not num exception
  Real total_variance = (nm1 > 0.) ? nm1*nm1 : nm1;

  // Solve for partial variance
  for (IntIntMIter map_iter=sobolIndexMap.begin();
       map_iter!=sobolIndexMap.end(); ++map_iter) {
    // partialVariance[0] stores the mean; it is not a component function
    // and does not follow the procedures for obtaining variance 
    if (map_iter->first) {
      compute_partial_variance(map_iter->first);
      sobolIndices[map_iter->second] = partialVariance[map_iter->second]
	                             / total_variance;
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
  totalSobolIndices = 0.; // init total indices
  if (expConfigOptions.vbdControl == ALL_VBD)
    for (IntIntMIter itr=sobolIndexMap.begin(); itr!=sobolIndexMap.end(); ++itr)
      for (int k=0; k<numVars; ++k) {
        if (itr->first & (1 << k))
          totalSobolIndices[k] += sobolIndices[itr->second];
        totalSobolIndices[k] = std::abs(totalSobolIndices[k]);
      }

  // If not available, compute total indices independently.
  // This approach parallels partial_variance_integral where the algorithm is 
  // separated by integration approach.
  else
    compute_total_sobol_indices(); // virtual

#ifdef DEBUG
  PCout << "In InterpPolyApproximation::compute_total_effects(), "
	<< "totalSobolIndices =\n"; write_data(PCout, totalSobolIndices);
#endif // DEBUG
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


/** Computes the variance of component functions. Assumes that all
    subsets of set_value have been computed in advance which will be
    true so long as the partial_variance is called following
    appropriate enumeration of set value  */
void InterpPolyApproximation::compute_partial_variance(int set_value)
{
  // derived classes override to define partialVariance and then invoke
  // base version for constituentSets post-processing

  // Now subtract the contributions from constituent subsets
  IntSet::iterator itr; int set_index = sobolIndexMap[set_value];
  for (itr  = constituentSets[set_index].begin();
       itr != constituentSets[set_index].end(); ++itr) 
    partialVariance[set_index] -= partialVariance[sobolIndexMap[*itr]];
}


/** Forms an interpolant over variables that are members of the given set.
    Finds the variance of the interpolant w.r.t. the variables in the set.
    Overloaded version supporting Smolyak sparse grids. */
Real InterpPolyApproximation::
partial_variance_integral(int set_value, const UShortArray& quad_order,
			  const UShortArray& lev_index, 
			  const UShort2DArray& key,
			  const SizetArray& colloc_index)
{
  RealVector member_coeffs, member_wts;
  member_coefficients_weights(set_value, quad_order, lev_index, key,
			      colloc_index, member_coeffs, member_wts);

  // Now integrate over the remaining variables	
  Real integral = 0.;
  size_t i, num_member_coeffs = member_coeffs.length();
  for (i=0; i<num_member_coeffs; ++i)
    integral += std::pow(member_coeffs[i], 2.) * member_wts[i];
  return integral;	
}


Real InterpPolyApproximation::
total_effects_integral(int set_value, const UShortArray& quad_order,
		       const UShortArray& lev_index, const UShort2DArray& key,
		       const SizetArray& colloc_index)
{
  RealVector member_coeffs, member_wts;
  member_coefficients_weights(set_value, quad_order, lev_index, key,
			      colloc_index, member_coeffs, member_wts);

  // Now integrate over the remaining variables	
  Real integral = 0.;
  const Real& total_mean = numericalMoments[0];
  size_t i, num_member_coeffs = member_coeffs.length();
  for (i=0; i<num_member_coeffs; ++i)
    integral += std::pow(member_coeffs[i] - total_mean, 2.) * member_wts[i];
  return integral;
}

} // namespace Pecos

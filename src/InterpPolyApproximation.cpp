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
#include "CombinedSparseGridDriver.hpp"
#include "HierarchSparseGridDriver.hpp"
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

  bool param_update = false;
  const std::vector<BasisPolynomial>& num_int_poly_basis
    = driverRep->polynomial_basis();
  for (size_t i=0; i<numVars; ++i)
    if (num_int_poly_basis[i].parametric_update())
      { param_update = true; break; }

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

    // can't use quad_order > quadOrderPrev logic since only 1 pt set is stored
    bool update_basis_form = (quad_order != quadOrderPrev);
    if (update_basis_form || param_update)
      update_tensor_interpolation_basis(tpq_driver->level_index());

    quadOrderPrev = quad_order;
    break;
  }
  case COMBINED_SPARSE_GRID: case HIERARCHICAL_SPARSE_GRID: {
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    unsigned short    ssg_level  = ssg_driver->level();
    const RealVector& aniso_wts  = ssg_driver->anisotropic_weights();

    // Ignore weights since they only reduce the interpolation depth from the
    // level and the basis update uses a coarse increment based on level.  This
    // matches isotropic sparse grids, but forces fewer and larger updates in
    // the case of anisotropic or generalized grids.
    bool update_basis_form
      = (ssgLevelPrev == USHRT_MAX || ssg_level > ssgLevelPrev);
    if (update_basis_form || param_update)
      update_sparse_interpolation_basis(ssg_level);

    ssgLevelPrev = ssg_level; ssgAnisoWtsPrev = aniso_wts;
    break;
  }
  }

  if (numericalMoments.empty()) {
    size_t num_moments = (nonRandomIndices.empty()) ? 4 : 2;
    numericalMoments.sizeUninitialized(num_moments);
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
  if (surrData.anchor()) // anchor point, if present, is first expansionSample
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
  unsigned short max_set_index = 0;
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case COMBINED_SPARSE_GRID: {
    // As for allocate_arrays(), increments are performed in coarser steps
    // than may be strictly necessary: all increments are filled in for all
    // vars for a step in level (ignoring anisotropy or generalized indices).
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    const UShortArray& trial_set = csg_driver->trial_set();
    for (size_t i=0; i<numVars; ++i)
      if (trial_set[i] > max_set_index)
	max_set_index = trial_set[i];
    break;
  }
  case HIERARCHICAL_SPARSE_GRID: {
    HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
    switch (expConfigOptions.refinementControl) {
    case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED: { // generalized sparse grids
      const UShortArray& trial_set = hsg_driver->trial_set();
      for (size_t i=0; i<numVars; ++i)
	if (trial_set[i] > max_set_index)
	  max_set_index = trial_set[i];
      break;
    }
    default: { // isotropic/anisotropic refinement
      const UShort3DArray&   sm_mi = hsg_driver->smolyak_multi_index();
      const UShortArray& incr_sets = hsg_driver->increment_sets();
      size_t lev, num_lev = sm_mi.size(), set, start_set, num_sets, v;
      for (lev=0; lev<num_lev; ++lev) {
	start_set = incr_sets[lev]; num_sets = sm_mi[lev].size();
	for (set=start_set; set<num_sets; ++set) {
	  const UShortArray& sm_set = sm_mi[lev][set];
	  for (v=0; v<numVars; ++v)
	    if (sm_set[v] > max_set_index)
	      max_set_index = sm_set[v];
	}
      }
      break;
    }
    }
    break;
  }
  default:
    PCerr << "Error: unsupported grid definition in InterpPolyApproximation::"
	  << "increment_coefficients()" << std::endl;
    abort_handler(-1);
    break;
  }
  update_sparse_interpolation_basis(max_set_index);

  increment_expansion_coefficients();
  numCollocPts = surrData.size(); if (surrData.anchor()) ++numCollocPts;
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
  numCollocPts = surrData.size(); if (surrData.anchor()) ++numCollocPts;
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
  numCollocPts = surrData.size(); if (surrData.anchor()) ++numCollocPts;
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
  numCollocPts = surrData.size(); if (surrData.anchor()) ++numCollocPts;
}


void InterpPolyApproximation::
update_tensor_interpolation_basis(const UShortArray& lev_index)
{
  // resize if needed (leaving previous levels unmodified)
  resize_polynomial_basis(lev_index);

  // fill any required gaps in polynomialBasis
  for (size_t i=0; i<numVars; ++i)
    update_interpolation_basis(lev_index[i], i);
}


void InterpPolyApproximation::
update_tensor_interpolation_basis(const UShortArray& lev_index,
				  const SizetList& subset_indices)
{
  // resize if needed (leaving previous levels unmodified)
  resize_polynomial_basis(lev_index);

  // fill any required gaps in polynomialBasis
  SizetList::const_iterator cit; size_t i;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit)
    { i = *cit; update_interpolation_basis(lev_index[i], i); }
}


void InterpPolyApproximation::
update_sparse_interpolation_basis(unsigned short max_level)
{
  // resize if needed (leaving previous levels unmodified)
  // j range is 0:w inclusive; i range is 1:w+1 inclusive
  size_t l, v, orig_size = polynomialBasis.size();
  resize_polynomial_basis(max_level);

  const std::vector<BasisPolynomial>& num_int_poly_basis
    = driverRep->polynomial_basis();
  // We must currently process all levels, even if not parameterized, since
  // update_interpolation_basis() can only update the polynomial basis if the
  // collocation points are available and the collocation points are only 
  // available if a particular index level is visited for that variable within 
  // IntegrationDriver::compute_tensor_grid() (which calls IntegrationDriver::
  // update_1d_collocation_points_weights()).  That is, there may be gaps that
  // need to be filled in levels that were previously allocated.
  for (v=0; v<numVars; ++v) {
    //if (num_int_poly_basis[v].parameterized()) // check all levels for updates
      for (l=0; l<=max_level; ++l)
	update_interpolation_basis(l, v);
    //else                                       // update only the new levels
    //  for (l=orig_size; l<=max_level; ++l)
    //    update_interpolation_basis(l, v);
  }
}


void InterpPolyApproximation::
update_interpolation_basis(unsigned short lev_index, size_t var_index)
{
  // fill gaps that may exist within any level
  const RealArray& colloc_pts_1d_lv
    = driverRep->collocation_points_1d()[lev_index][var_index];
  if (!colloc_pts_1d_lv.empty()) {
    const BasisPolynomial& num_int_poly_basis_v
      = driverRep->polynomial_basis()[var_index];
    std::vector<BasisPolynomial>& poly_basis_l  = polynomialBasis[lev_index];
    BasisPolynomial&              poly_basis_lv = poly_basis_l[var_index];
    short poly_type_1d, rule;
    // don't share reps in case of parameterized basis or barycentric interp,
    // due to need for individual updates to parameters or interpolated x value.
    if (num_int_poly_basis_v.parameterized() || barycentricFlag) {
      if (poly_basis_lv.is_null()) {
	initialize_polynomial_basis_type(poly_type_1d, rule);
	poly_basis_lv = BasisPolynomial(poly_type_1d, rule);
	poly_basis_lv.interpolation_points(colloc_pts_1d_lv);
      }
      else if (num_int_poly_basis_v.parametric_update())
	poly_basis_lv.interpolation_points(colloc_pts_1d_lv);
    }
    else if (poly_basis_lv.is_null()) { // can share reps for efficiency
      size_t var_index2;
      if (find_basis(lev_index, var_index, var_index2))
	poly_basis_lv = poly_basis_l[var_index2]; // reuse prev via shared rep
      else { // instantiate and initialize a new unique instance
	initialize_polynomial_basis_type(poly_type_1d, rule);
	poly_basis_lv = BasisPolynomial(poly_type_1d, rule);
	poly_basis_lv.interpolation_points(colloc_pts_1d_lv);
      }
    }
  }
}


bool InterpPolyApproximation::
find_basis(unsigned short level, size_t v1, size_t& v2)
{
  std::vector<BasisPolynomial>& poly_basis_l = polynomialBasis[level];
  for (v2=0; v2<numVars; ++v2)
    if (v2 != v1 && !poly_basis_l[v2].is_null() && same_basis(level, v1, v2))
      return true;
  return false; 
}


bool InterpPolyApproximation::
same_basis(unsigned short level, size_t v1, size_t v2)
{
  const ShortArray& rules = driverRep->collocation_rules();
  short rule1 = rules[v1];
  if (rules[v2] == rule1)
    switch (rule1) {
    case GAUSS_JACOBI: case GEN_GAUSS_LAGUERRE: case GOLUB_WELSCH: {
      // rule type insufficient in these cases, check collocation points
      const Real2DArray& colloc_pts_1d
	= driverRep->collocation_points_1d()[level];
      return (colloc_pts_1d[v1] == colloc_pts_1d[v2]); break;
    }
    default:
      return true;                                     break;
    }
  else
    return false;
}


/** Barycentric approach is only valid for value-based global Lagrange
    interpolation, either nodal or hierarchical.  General approach is
    valid for value-based or gradient-enhanced, local or global, and
    nodal or hierarchical. */
Real InterpPolyApproximation::
tensor_product_value(const RealVector& x, const RealVector& exp_t1_coeffs,
		     const RealMatrix& exp_t2_coeffs,
		     const UShortArray& basis_index, const UShort2DArray& key,
		     const SizetArray& colloc_index)
{
  Real tp_val = 0.; size_t i, num_colloc_pts = key.size();
  if (barycentricFlag) {

    // For barycentric interpolation: track x != newPoint within 1D basis
    set_new_point(x, basis_index, 1); // value factors needed

    // apply Horner's rule to group common factors
    // TO DO: manage inactive interpolation dimensions
    RealVector accumulator(numVars); size_t j, v_index;
    RealVectorArray bc_val_facts = barycentric_value_factors_array(basis_index);
    // TO DO: max_key not valid for delta_quad on closed rules
    UShortArray max_key(numVars);
    for (i=0; i<numVars; ++i)
      max_key[i] = bc_val_facts[i].length() - 1;
    const RealVector& bc_val_fact_0 = bc_val_facts[0];
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      accumulator[0] += (colloc_index.empty()) ?
	exp_t1_coeffs[i]               * bc_val_fact_0[key_i[0]] :
	exp_t1_coeffs[colloc_index[i]] * bc_val_fact_0[key_i[0]];
      if (key_i[0] == max_key[0]) {
	// compute number of variables (v) to accumulate
	v_index = 1;
	while (key_i[v_index] == max_key[v_index] && v_index < numVars-1)
	  ++v_index;
	// accumulate sums over v variables
	for (j=1; j<=v_index; ++j) {
	  accumulator[j]  += accumulator[j-1] * bc_val_facts[j][key_i[j]];
	  accumulator[j-1] = 0.;
	}
      }
    }
    tp_val = accumulator[v_index] * barycentric_value_scaling(basis_index);

    /*
    // delegate loops: cleaner, but some efficiency lost (and possibly precision
    // as well due to subtractive cancellation among large products)
    if (colloc_index.empty())
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += exp_t1_coeffs[i] *
	  barycentric_value_factor(key[i], basis_index);
    else
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += exp_t1_coeffs[colloc_index[i]] *
	  barycentric_value_factor(key[i], basis_index);

    // apply barycentric denominator
    tp_val *= barycentric_value_scaling(basis_index);
    */
  }
  else if (exp_t2_coeffs.empty()) {
    if (colloc_index.empty())
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += exp_t1_coeffs[i] * 
	          type1_interpolant_value(x, key[i], basis_index);
    else
      for (i=0; i<num_colloc_pts; ++i) {
	tp_val += exp_t1_coeffs[colloc_index[i]] *
	          type1_interpolant_value(x, key[i], basis_index);
	PCout << "  t1c = " << exp_t1_coeffs[colloc_index[i]] << " t1i = "
	      << type1_interpolant_value(x, key[i], basis_index)
	      << " tp_val = " << tp_val << '\n';
      }
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


/** All variables version. */
Real InterpPolyApproximation::
tensor_product_value(const RealVector& x, const RealVector& exp_t1_coeffs,
		     const RealMatrix& exp_t2_coeffs,
		     const UShortArray& basis_index, const UShort2DArray& key,
		     const SizetArray& colloc_index,
		     const SizetList& subset_indices)
{
  Real tp_val = 0.; size_t i, num_colloc_pts = key.size();
  if (barycentricFlag) {

    // For barycentric interpolation: track x != newPoint within 1D basis
    set_new_point(x, basis_index, subset_indices, 1); // value factors needed

    if (colloc_index.empty())
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += exp_t1_coeffs[i] * 
	  barycentric_value_factor(key[i], basis_index, subset_indices);
    else
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += exp_t1_coeffs[colloc_index[i]] *
	  barycentric_value_factor(key[i], basis_index, subset_indices);

    // apply barycentric denominator
    tp_val *= barycentric_value_scaling(basis_index, subset_indices);
  }
  else if (exp_t2_coeffs.empty()) {
    if (colloc_index.empty())
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += exp_t1_coeffs[i] * 
	  type1_interpolant_value(x, key[i], basis_index, subset_indices);
    else
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += exp_t1_coeffs[colloc_index[i]] *
	  type1_interpolant_value(x, key[i], basis_index, subset_indices);
  }
  else {
    size_t j, c_index;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      c_index = (colloc_index.empty()) ? i : colloc_index[i];
      tp_val += exp_t1_coeffs[c_index] *
	type1_interpolant_value(x, key_i, basis_index, subset_indices);
      const Real* exp_t2_coeff_i = exp_t2_coeffs[c_index];
      for (j=0; j<numVars; ++j)
	tp_val += exp_t2_coeff_i[j] *
	  type2_interpolant_value(x, j, key_i, basis_index, subset_indices);
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

  if (barycentricFlag) {

    // For barycentric interpolation: track x != newPoint within 1D basis
    set_new_point(x, basis_index, 3); // value+gradient factors needed

    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray&   key_i = key[i];
      const Real& exp_t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      for (j=0; j<numVars; ++j)
	tpGradient[j] += exp_t1_coeff_i *
	  barycentric_gradient_factor(j, key_i, basis_index);
    }

    tpGradient.scale(barycentric_gradient_scaling(basis_index));
  }
  else if (exp_t2_coeffs.empty()) {
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray&   key_i = key[i];
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
      const UShortArray&   key_i = key[i];
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
					const SizetList& subset_indices)
{
  if (tpGradient.length() != numVars)
    tpGradient.sizeUninitialized(numVars);
  tpGradient = 0.;
  size_t i, j, num_colloc_pts = key.size();

  if (barycentricFlag) {

    // For barycentric interpolation: track x != newPoint within 1D basis
    set_new_point(x, basis_index, subset_indices, 3);//value+grad factors needed

    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray&   key_i = key[i];
      const Real& exp_t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      for (j=0; j<numVars; ++j)
	tpGradient[j] += exp_t1_coeff_i *
	  barycentric_gradient_factor(j, key_i, basis_index, subset_indices);
    }

    tpGradient.scale(barycentric_gradient_scaling(basis_index, subset_indices));
  }
  else if (exp_t2_coeffs.empty()) {
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray&   key_i = key[i];
      const Real& exp_t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      for (j=0; j<numVars; ++j)
	tpGradient[j] += exp_t1_coeff_i *
	  type1_interpolant_gradient(x, j, key_i, basis_index, subset_indices);
    }
  }
  else {
    size_t k;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray&   key_i = key[i];
      const Real& exp_t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      const Real* exp_t2_coeff_i = (colloc_index.empty()) ?
	exp_t2_coeffs[i] : exp_t2_coeffs[colloc_index[i]];
      for (j=0; j<numVars; ++j) { // ith contribution to jth grad component
	tpGradient[j] += exp_t1_coeff_i *
	  type1_interpolant_gradient(x, j, key_i, basis_index, subset_indices);
	for (k=0; k<numVars; ++k) // type2 interpolant for kth grad comp
	  tpGradient[j] += exp_t2_coeff_i[k] *
	    type2_interpolant_gradient(x, j, k, key_i, basis_index,
				       subset_indices);
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

  if (barycentricFlag) {

    // For barycentric interpolation: track x != newPoint within 1D basis
    set_new_point(x, basis_index, 3); // value+gradient factors needed

    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray&   key_i = key[i];
      const Real& exp_t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      for (j=0; j<num_deriv_vars; ++j) {
	deriv_index = dvv[j] - 1; // requires an "All" view
	tpGradient[j] += exp_t1_coeff_i *
	  barycentric_gradient_factor(deriv_index, key_i, basis_index);
      }
    }

    tpGradient.scale(barycentric_gradient_scaling(basis_index));
  }
  else if (exp_t2_coeffs.empty()) {
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray&   key_i = key[i];
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
      const UShortArray&   key_i = key[i];
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
  if (barycentricFlag) {

    // For barycentric interpolation: track x != newPoint within 1D basis
    set_new_point(x, basis_index, 1); // value factors needed for coeff grads

    for (i=0; i<num_colloc_pts; ++i) {
      const Real* exp_t1_coeff_grad_i = (colloc_index.empty()) ?
	exp_t1_coeff_grads[i] : exp_t1_coeff_grads[colloc_index[i]];
      Real bc_fact = barycentric_value_factor(key[i], basis_index);
      for (j=0; j<num_deriv_vars; ++j)
	tpGradient[j] += exp_t1_coeff_grad_i[j] * bc_fact;
    }

    // apply barycentric denominator
    tpGradient.scale(barycentric_value_scaling(basis_index));
  }
  else {
    for (i=0; i<num_colloc_pts; ++i) {
      const Real* exp_t1_coeff_grad_i = (colloc_index.empty()) ?
	exp_t1_coeff_grads[i] : exp_t1_coeff_grads[colloc_index[i]];
      Real t1_val = type1_interpolant_value(x, key[i], basis_index);
      for (j=0; j<num_deriv_vars; ++j)
	tpGradient[j] += exp_t1_coeff_grad_i[j] * t1_val;
    }
  }
  return tpGradient;
}


void InterpPolyApproximation::compute_component_effects()
{
  // initialize partialVariance
  if (partialVariance.empty()) partialVariance.size(sobolIndices.length());
  else                         partialVariance = 0.;
  partialVariance[0] = numericalMoments[0]*numericalMoments[0];// init w/ mean^2

  // Compute the total expansion variance.  For standard mode, the full variance
  // is likely already available, as managed by computedVariance in variance().
  // For all variables mode, we use covariance(this) without passing x for the
  // nonRandomIndices (bypass computedVariance checks by not using variance()).
  Real total_variance = (nonRandomIndices.empty()) ? variance() : // std mode
                        covariance(this);                    // all vars mode

  // Solve for partial variances
  for (BAULMIter it=sobolIndexMap.begin(); it!=sobolIndexMap.end(); ++it) {
    unsigned long index = it->second;
    if (index) { // partialVariance[0] stores mean; no variance to calculate
      compute_partial_variance(it->first);
      sobolIndices[index] = partialVariance[index] / total_variance;
    }
  }
#ifdef DEBUG
  PCout << "In InterpPolyApproximation::compute_component_effects(), "
	<< "sobolIndices =\n"; write_data(PCout, sobolIndices);
#endif // DEBUG
}


void InterpPolyApproximation::compute_total_effects()
{
  totalSobolIndices = 0.; // init total indices

  // iterate through existing indices if all component indices are available.
  // totalSobolIndices simply parse the bit sets of each of the sobolIndices 
  // and add them to each matching variable bin.
  if (expConfigOptions.vbdControl == ALL_VBD) {
    for (BAULMIter it=sobolIndexMap.begin(); it!=sobolIndexMap.end(); ++it)
      for (size_t k=0; k<numVars; ++k)
        if (it->first[k]) // var k is present in this Sobol' index
          totalSobolIndices[k] += sobolIndices[it->second];
    // ensure non-negativity of indices
    //for (size_t k=0; k<numVars; ++k)
    //  totalSobolIndices[k] = std::abs(totalSobolIndices[k]);
  }

  // If not available, compute total indices independently.  This approach
  // parallels partial_variance_integral where the algorithm is separated
  // by integration approach.
  else
    compute_total_sobol_indices(); // virtual

#ifdef DEBUG
  PCout << "In InterpPolyApproximation::compute_total_effects(), "
	<< "totalSobolIndices =\n"; write_data(PCout, totalSobolIndices);
#endif // DEBUG
}


/** Computes the variance of component functions.  Assumes that partial
    variances of all subsets of set_value have been computed in advance:
    compute_component_effects() calls compute_partial_variance() using
    the ordered set_value's in sobolIndexMap. */
void InterpPolyApproximation::
compute_partial_variance(const BitArray& set_value)
{
  // derived classes override to compute partialVariance and then invoke
  // base version for post-processing of proper subsets

  // compute child subsets.  An alternate approach would be to iterate
  // over sobolIndexMap using it->is_proper_subset_of(set_value).
  BitArraySet children;
  proper_subsets(set_value, children);

  // index of parent set within sobolIndices and partialVariance
  unsigned long set_index;
  if (!children.empty())
    set_index = sobolIndexMap[set_value];

  // subtract the contributions from child subsets.  partialVariance
  // calculations are computed by ordered traversal of sobolIndexMap, which
  // uses BitArray keys that sort lexicographically (in ascending order of
  // the corresponding unsigned integer) --> all proper subsets are available.
  for (BASIter it=children.begin(); it!=children.end(); ++it) {
    unsigned long subset_index  = sobolIndexMap[*it];
    partialVariance[set_index] -= partialVariance[subset_index];
  }
}


/** For input parent set, recursively finds constituent child subsets
    with one fewer element */
void InterpPolyApproximation::
proper_subsets(const BitArray& parent_set, BitArraySet& children)
{
  for (size_t k=0; k<numVars; ++k)
    if (parent_set[k]) { // check for membership of variable k in parent set
      // remove var k from parent set to create child set
      BitArray child_set = parent_set; child_set[k].flip();
      // if child set has not been stored previously, insert it and recurse
      if (children.find(child_set) == children.end()) {
	children.insert(child_set);
	proper_subsets(child_set, children); // recurse until {0} is reached
      }
    }
}

} // namespace Pecos

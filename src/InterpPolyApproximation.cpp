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
      Real      resp_fn = surrData.response_function(index);
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
  if (barycentricFlag) {

    // For barycentric interpolation: track x != newPoint within 1D basis
    set_new_point(x, basis_index, 1); // value factors needed

    // TO DO: max key limit not valid for delta_quad on closed rules

    size_t j, num_act_v = barycentric_active_variables(basis_index);
    if (num_act_v == 0) { // convert 1-D exact indices into n-D colloc index
      size_t pt_index = barycentric_exact_index(basis_index);
      return (colloc_index.empty()) ?
	exp_t1_coeffs[pt_index] : exp_t1_coeffs[colloc_index[pt_index]];
    }
    else if (num_act_v == numVars) { // interpolation over all variables
      RealVector accumulator(numVars); // init to 0.
      BasisPolynomial&   poly_0 = polynomialBasis[basis_index[0]][0];
      const RealVector& bc_vf_0 = poly_0.barycentric_value_factors();
      size_t i, num_colloc_pts = key.size();//, ei0 = poly_0.exact_index(), eij;
      unsigned short key_i0, key_ij, max0 = poly_0.interpolation_size() - 1;
      for (i=0; i<num_colloc_pts; ++i) {
	const UShortArray& key_i = key[i]; key_i0 = key_i[0];
	//if (ei0 == _NPOS)
	accumulator[0] += (colloc_index.empty()) ?
	  exp_t1_coeffs[i]               * bc_vf_0[key_i0] :
	  exp_t1_coeffs[colloc_index[i]] * bc_vf_0[key_i0];
	//else if (ei0 == key_i0)                      // only 1 pt should match
	//  accumulator[0]  = (colloc_index.empty()) ? // (= instead of +=)
	//    exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
	if (key_i0 == max0) {
	  // accumulate sums over variables with max key value
	  for (j=1; j<numVars; ++j) {
	    BasisPolynomial& poly_j = polynomialBasis[basis_index[j]][j];
	    key_ij = key_i[j]; //eij = poly_j.exact_index();
	    //if (eij == _NPOS)
	    accumulator[j] += accumulator[j-1]
	      * poly_j.barycentric_value_factor(key_ij);
	    //else if (eij == key_ij)              // only 1 pt should match
	    //  accumulator[j] = accumulator[j-1]; // (= instead of +=)
	    accumulator[j-1] = 0.;
	    if (key_ij + 1 != poly_j.interpolation_size())
	      break;
	  }
	}
      }
      return accumulator[numVars-1] / barycentric_value_scaling(basis_index);
    }
    else { // partial interpolation over active variables
      SizetArray pt_factors(num_act_v), act_v_set(num_act_v);
      size_t pts_vj, num_act_pts = 1, num_pts = 1, ej, av_cntr, pt_index = 0;
      unsigned short bi_j;
      for (j=0, av_cntr=0; j<numVars; ++j) {
	bi_j = basis_index[j];
	if (bi_j) { // else pts_vj = 1 and ej can be taken to be 0
	  BasisPolynomial& poly_j = polynomialBasis[bi_j][j];
	  ej = poly_j.exact_index(); pts_vj = poly_j.interpolation_size();
	  if (ej == _NPOS) { // active for interpolation
	    pt_factors[av_cntr] = num_pts; act_v_set[av_cntr] = j;
	    num_act_pts *= pts_vj; ++av_cntr;
	  }
	  else             // inactive for interpolation
	    pt_index += num_pts * ej;
	  num_pts *= pts_vj;
	}
      }
      // define initial pt_index offset
      RealVector accumulator(num_act_v); // init to 0.
      size_t i, v0 = act_v_set[0], vj;
      BasisPolynomial&   poly_v0 = polynomialBasis[basis_index[v0]][v0];
      const RealVector& bc_vf_v0 = poly_v0.barycentric_value_factors();
      size_t pts_v0 = poly_v0.interpolation_size(), pf0 = pt_factors[0],
	 pts_v0_pf0 = pts_v0 * pf0, pfj, prev_pt_set;
      unsigned short key_i0, key_ij, max0 = pts_v0 - 1;
      // loop over active pts, summing contributions from active variables
      for (i=0; i<num_act_pts; ++i) {
	const UShortArray& key_i = key[pt_index]; key_i0 = key_i[v0];
	accumulator[0] += (colloc_index.empty()) ? 
	  exp_t1_coeffs[pt_index]               * bc_vf_v0[key_i0] :
	  exp_t1_coeffs[colloc_index[pt_index]] * bc_vf_v0[key_i0];
	pt_index += pf0;
	if (key_i0 == max0) {
	  // accumulate sums over variables with max key value
	  for (j=1, prev_pt_set=pts_v0_pf0; j<num_act_v; ++j) {
	    vj = act_v_set[j]; key_ij = key_i[vj]; pfj = pt_factors[j];
	    BasisPolynomial& poly_vj = polynomialBasis[basis_index[vj]][vj];
	    // update accumulators: push [j-1] entry up to [j] level
	    accumulator[j]  += accumulator[j-1]
	      * poly_vj.barycentric_value_factor(key_ij);
	    accumulator[j-1] = 0.;
	    // update pt_index: prev index rolls back to 0 and curr index +1.
	    // index increment is zero unless active vars are nonconsecutive.
	    if (pfj != prev_pt_set)
	      pt_index += pfj - prev_pt_set;
	    pts_vj = poly_vj.interpolation_size();
	    if (key_ij + 1 == pts_vj) prev_pt_set = pts_vj * pfj;
	    else                      break;
	  }
	}
      }
      return accumulator[num_act_v-1] / barycentric_value_scaling(basis_index);
    }

    /*
    // Simple bc option: delegate loops; cleaner, but some efficiency lost (and
    // precision as well due to subtractive cancellation among large products)
    size_t i, num_colloc_pts = key.size(); Real tp_val = 0.;
    if (colloc_index.empty())
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += exp_t1_coeffs[i] *
	  barycentric_value_factor(key[i], basis_index);
    else
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += exp_t1_coeffs[colloc_index[i]] *
	  barycentric_value_factor(key[i], basis_index);
    // apply barycentric denominator
    return tp_val / barycentric_value_scaling(basis_index);
    */
  }
  else if (exp_t2_coeffs.empty()) {
    /*
    // Horner's rule approach:
    RealVector accumulator(numVars); // init to 0.
    unsigned short bi_0 = basis_index[0], bi_j;
    BasisPolynomial& poly_0 = polynomialBasis[bi_0][0];
    unsigned short key_i0, key_ij, max0 = poly_0.interpolation_size() - 1;
    size_t i, j, num_colloc_pts = key.size(); Real x0 = x[0];
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i]; key_i0 = key_i[0];
      if (bi_0)
	accumulator[0] += (colloc_index.empty()) ?
	  exp_t1_coeffs[i]               * poly_0.type1_value(x0, key_i0) :
	  exp_t1_coeffs[colloc_index[i]] * poly_0.type1_value(x0, key_i0);
      else
	accumulator[0] = (colloc_index.empty()) ?
	  exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      if (key_i0 == max0) {
	// accumulate sums over variables with max key value
	for (j=1; j<numVars; ++j) {
	  key_ij = key_i[j]; bi_j = basis_index[j];
	  BasisPolynomial& poly_j = polynomialBasis[bi_j][j];
	  if (bi_j)
	    accumulator[j] += accumulator[j-1] *
	      poly_j.type1_value(x[j], key_ij);
	  else
	    accumulator[j]  = accumulator[j-1];
	  accumulator[j-1] = 0.;
	  if (key_ij + 1 != poly_j.interpolation_size())
	    break;
	}
      }
    }
    return accumulator[numVars-1];
    */

    // Simpler but more expensive approach:
    size_t i, num_colloc_pts = key.size(); Real tp_val = 0.;
    if (colloc_index.empty())
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += exp_t1_coeffs[i] * 
	          type1_interpolant_value(x, key[i], basis_index);
    else
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += exp_t1_coeffs[colloc_index[i]] *
	          type1_interpolant_value(x, key[i], basis_index);
    return tp_val;
  }
  else {
    /*
    // Horner's rule approach:
    RealVector t1_accumulator(numVars);          // init to 0.
    RealMatrix t2_accumulator(numVars, numVars); // init to 0.
    Real *t2_accum_0 = t2_accumulator[0], *t2_accum_j, *t2_accum_jm1;
    unsigned short       bi_0 = basis_index[0], bi_j;
    BasisPolynomial&   poly_0 = polynomialBasis[bi_0][0];
    size_t i, j, jm1, k, c_index, num_colloc_pts = key.size();
    unsigned short key_i0, key_ij, max0 = poly_0.interpolation_size() - 1;
    Real t1_val, x0 = x[0];
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i]; key_i0 = key_i[0];
      c_index = (colloc_index.empty()) ? i : colloc_index[i];
      const Real* t2_coeffs_i = exp_t2_coeffs[c_index];
      if (bi_0 == 0) {
	t1_accumulator[0]  = exp_t1_coeffs[c_index];  // t1 value is 1
	t2_accum_0[0]      = t2_coeffs_i[0] * poly_0.type2_value(x0, key_i0);
	for (k=1; k<numVars; ++k)
	  t2_accum_0[k]    = t2_coeffs_i[k];                 // t1 value is 1
      }
      else {
	t1_val = poly_0.type1_value(x0, key_i0);
	t1_accumulator[0] += exp_t1_coeffs[c_index] * t1_val;
	t2_accum_0[0]     += t2_coeffs_i[0] * poly_0.type2_value(x0, key_i0);
	for (k=1; k<numVars; ++k)
	  t2_accum_0[k]   += t2_coeffs_i[k] * t1_val;
      }
      if (key_i0 == max0) {
	// accumulate sums over variables with max key value
	for (j=1; j<numVars; ++j) {
	  bi_j = basis_index[j]; key_ij = key_i[j]; jm1 = j-1;
	  t2_accum_j = t2_accumulator[j]; t2_accum_jm1 = t2_accumulator[jm1];
	  BasisPolynomial& poly_j = polynomialBasis[bi_j][j];
	  if (bi_j == 0) {
	    t1_accumulator[j] = t1_accumulator[jm1];         // t1 value is 1
	    t2_accum_j[j] = t2_accum_jm1[j] * poly_j.type2_value(x[j], key_ij);
	    for (k=0; k<numVars; ++k)
	      if (k != j)
		t2_accum_j[k] = t2_accum_jm1[k];             // t1 value is 1
	  }
	  else {
	    t1_val = poly_j.type1_value(x[j], key_ij);
	    t1_accumulator[j] += t1_accumulator[jm1] * t1_val;
	    t2_accum_j[j] += t2_accum_jm1[j] * poly_j.type2_value(x[j], key_ij);
	    for (k=0; k<numVars; ++k)
	      if (k != j)
		t2_accum_j[k] += t2_accum_jm1[k] * t1_val;
	  }
	  t1_accumulator[jm1] = 0.;
	  for (k=0; k<numVars; ++k)
	    t2_accum_jm1[k] = 0.;
	  if (key_ij + 1 != poly_j.interpolation_size())
	    break;
	}
      }
    }
    Real  tp_val   = t1_accumulator[numVars-1];
    Real* t2_accum = t2_accumulator[numVars-1];
    for (j=0; j<numVars; ++j)
      tp_val += t2_accum[j];
    return tp_val;
    */

    // Simpler but more expensive approach:
    size_t i, num_colloc_pts = key.size(), j, c_index; Real tp_val = 0.;
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
    return tp_val;
  }
}


/** All variables version. */
Real InterpPolyApproximation::
tensor_product_value(const RealVector& x, const RealVector& subset_t1_coeffs,
		     const RealMatrix& subset_t2_coeffs,
		     const UShortArray& basis_index,
		     const UShort2DArray& subset_key,
		     const SizetArray& subset_colloc_index,
		     const SizetList& subset_indices)
{
  // Note: subset_* are consistent with the reduced variable subset,
  // but x and basis_index are full space.

  if (barycentricFlag) {

    // For barycentric interpolation: track x != newPoint within 1D basis
    set_new_point(x, basis_index, subset_indices, 1); // value factors needed

    size_t num_subset_v = subset_indices.size(),
      num_act_v = barycentric_active_variables(basis_index, subset_indices);
    if (num_act_v == 0) { // convert 1-D exact indices into n-D colloc index
      size_t pt_index = barycentric_exact_index(basis_index, subset_indices);
      // Note: pt_index calculation utilizes only the subset variables
      return (subset_colloc_index.empty()) ? subset_t1_coeffs[pt_index] :
	subset_t1_coeffs[subset_colloc_index[pt_index]];
    }
    else if (num_act_v == num_subset_v) { // interpolation over all of subset
      SizetList::const_iterator cit = subset_indices.begin();
      size_t i, j, num_colloc_pts = subset_key.size(), v0 = *cit, vj;
      RealVector accumulator(num_subset_v); // init to 0.
      BasisPolynomial&   poly_0 = polynomialBasis[basis_index[v0]][v0];
      const RealVector& bc_vf_0 = poly_0.barycentric_value_factors();
      unsigned short key_i0, key_ij, max0 = poly_0.interpolation_size() - 1;
      for (i=0; i<num_colloc_pts; ++i) {
	// Note: key is reduced in first dimension (pts), but not second (vars)
	const UShortArray& key_i = subset_key[i]; key_i0 = key_i[v0];
	accumulator[0] += (subset_colloc_index.empty()) ?
	  subset_t1_coeffs[i]                      * bc_vf_0[key_i0] :
	  subset_t1_coeffs[subset_colloc_index[i]] * bc_vf_0[key_i0];
	if (key_i0 == max0) {
	  // accumulate sums over variables with max key value
	  for (j=1, cit=++subset_indices.begin(); j<num_subset_v; ++j, ++cit) {
	    vj = *cit; key_ij = key_i[vj];
	    BasisPolynomial& poly_j = polynomialBasis[basis_index[vj]][vj];
	    accumulator[j] += accumulator[j-1]
	      * poly_j.barycentric_value_factor(key_ij);
	    accumulator[j-1] = 0.;
	    if (key_ij + 1 != poly_j.interpolation_size())
	      break;
	  }
	}
      }
      return accumulator[num_subset_v-1]
	/ barycentric_value_scaling(basis_index, subset_indices);
    }
    else { // partial interpolation over active variables
      unsigned short bi_j; SizetList::const_iterator cit;
      SizetArray pt_factors(num_act_v), act_v_set(num_act_v);
      size_t i, j, pts_vj, num_act_pts = 1, num_pts = 1, ej, av_cntr = 0,
	pt_index = 0;
      // define interpolation set: in subset_indices, nonzero bi, and no ei
      for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
	j = *cit; bi_j = basis_index[j];
	if (bi_j) { // else pts_vj = 1 and ej can be taken to be 0
	  BasisPolynomial& poly_j = polynomialBasis[bi_j][j];
	  ej = poly_j.exact_index(); pts_vj = poly_j.interpolation_size();
	  if (ej == _NPOS) { // active for interpolation
	    pt_factors[av_cntr] = num_pts; act_v_set[av_cntr] = j;
	    num_act_pts *= pts_vj; ++av_cntr;
	  }
	  else             // inactive for interpolation
	    pt_index += num_pts * ej;
	  num_pts *= pts_vj;
	}
      }
      // define initial pt_index offset
      RealVector accumulator(num_act_v); // init to 0.
      size_t v0 = act_v_set[0], vj;
      BasisPolynomial&   poly_v0 = polynomialBasis[basis_index[v0]][v0];
      const RealVector& bc_vf_v0 = poly_v0.barycentric_value_factors();
      size_t pts_v0 = poly_v0.interpolation_size(), pf0 = pt_factors[0],
	 pts_v0_pf0 = pts_v0 * pf0, pfj, prev_pt_set;
      unsigned short key_i0, key_ij, max0 = pts_v0 - 1;
      // loop over active pts, summing contributions from active variables
      for (i=0; i<num_act_pts; ++i) {
	const UShortArray& key_i = subset_key[pt_index]; key_i0 = key_i[v0];
	accumulator[0] += (subset_colloc_index.empty()) ? 
	  subset_t1_coeffs[pt_index]                      * bc_vf_v0[key_i0] :
	  subset_t1_coeffs[subset_colloc_index[pt_index]] * bc_vf_v0[key_i0];
	pt_index += pf0;
	if (key_i0 == max0) {
	  // accumulate sums over variables with max key value
	  for (j=1, prev_pt_set=pts_v0_pf0; j<num_act_v; ++j) {
	    vj = act_v_set[j]; key_ij = key_i[vj]; pfj = pt_factors[j];
	    BasisPolynomial& poly_vj = polynomialBasis[basis_index[vj]][vj];
	    // update accumulators: push [j-1] entry up to [j] level
	    accumulator[j]  += accumulator[j-1]
	      * poly_vj.barycentric_value_factor(key_ij);
	    accumulator[j-1] = 0.;
	    // update pt_index: prev index rolls back to 0 and curr index +1.
	    // index increment is zero unless active vars are nonconsecutive.
	    if (pfj != prev_pt_set)
	      pt_index += pfj - prev_pt_set;
	    pts_vj = poly_vj.interpolation_size();
	    if (key_ij + 1 == pts_vj) prev_pt_set = pts_vj * pfj;
	    else                      break;
	  }
	}
      }
      return accumulator[num_act_v-1]
	/ barycentric_value_scaling(basis_index, subset_indices);
    }

    /*
    // Simple bc option: delegate loops; cleaner, but less efficient/precise
    size_t i, num_colloc_pts = subset_key.size(); Real tp_val = 0.;
    if (subset_colloc_index.empty())
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += subset_t1_coeffs[i] * 
	  barycentric_value_factor(subset_key[i], basis_index, subset_indices);
    else
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += subset_t1_coeffs[subset_colloc_index[i]] *
	  barycentric_value_factor(subset_key[i], basis_index, subset_indices);
    // apply barycentric denominator
    return tp_val / barycentric_value_scaling(basis_index, subset_indices);
    */
  }
  else if (subset_t2_coeffs.empty()) {
    size_t i, num_colloc_pts = subset_key.size(); Real tp_val = 0.;
    if (subset_colloc_index.empty())
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += subset_t1_coeffs[i] * 
	  type1_interpolant_value(x, subset_key[i], basis_index,subset_indices);
    else
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += subset_t1_coeffs[subset_colloc_index[i]] *
	  type1_interpolant_value(x, subset_key[i], basis_index,subset_indices);
    return tp_val;
  }
  else {
    size_t i, j, num_colloc_pts = subset_key.size(), c_index; Real tp_val = 0.;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = subset_key[i];
      c_index = (subset_colloc_index.empty()) ? i : subset_colloc_index[i];
      tp_val += subset_t1_coeffs[c_index] *
	type1_interpolant_value(x, key_i, basis_index, subset_indices);
      const Real* subset_t2_coeff_i = subset_t2_coeffs[c_index];
      for (j=0; j<numVars; ++j)
	tp_val += subset_t2_coeff_i[j] *
	  type2_interpolant_value(x, j, key_i, basis_index, subset_indices);
    }
    return tp_val;
  }
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

    // Note: for basis gradients, one cannot eliminate exact index dimensions
    // as is managed in the barycentric case of tensor_product_value()

    unsigned short key_i0, key_ij, bi0 = basis_index[0], bij;
    BasisPolynomial&   poly_0 = polynomialBasis[bi0][0];
    const RealVector& bc_vf_0 = poly_0.barycentric_value_factors();
    const RealVector& bc_gf_0 = poly_0.barycentric_gradient_factors();
    size_t k, ei0 = poly_0.exact_index(), eij;
    unsigned short max0 = poly_0.interpolation_size() - 1;
    RealMatrix accumulator(numVars, numVars); // init to 0.
    Real *accum_0 = accumulator[0], *accum_j, *accum_jm1,
      t1_coeff_bc_vf_00, bc_vf_jj, t1_coeff;
    for (i=0; i<num_colloc_pts; ++i) {
      t1_coeff = (colloc_index.empty()) ? exp_t1_coeffs[i] :
	exp_t1_coeffs[colloc_index[i]];
      const UShortArray& key_i = key[i]; key_i0 = key_i[0];
      if (bi0) {
	// accumulate grad factor for comp 0 and value factor for other comps
	accum_0[0] += t1_coeff * bc_gf_0[key_i0];
	if (ei0 == _NPOS) {
	  t1_coeff_bc_vf_00 = t1_coeff * bc_vf_0[key_i0];
	  for (j=1; j<numVars; ++j)
	    accum_0[j] += t1_coeff_bc_vf_00;
	}
	else if (ei0 == key_i0) // value factor is 1
	  for (j=1; j<numVars; ++j)
	    accum_0[j] += t1_coeff;
	//else value factor is 0
      }
      else // grad factor is zero, value factor is omitted
	for (j=1; j<numVars; ++j)
	  accum_0[j] += t1_coeff;
      if (key_i0 == max0) {
	// accumulate sums over variables with max key value
	for (j=1; j<numVars; ++j) {
	  bij     = basis_index[j]; key_ij    = key_i[j];
	  accum_j = accumulator[j]; accum_jm1 = accumulator[j-1];
	  BasisPolynomial& poly_j = polynomialBasis[bij][j];
	  if (bij) {
	    eij = poly_j.exact_index();
	    accum_j[j] += accum_jm1[j] *
	      poly_j.barycentric_gradient_factor(key_ij);
	    if (eij == _NPOS) { // bc_vf_jj has general value
	      bc_vf_jj = poly_j.barycentric_value_factor(key_ij);
	      for (k=0; k<numVars; ++k) {
		if (k != j) accum_j[k] += accum_jm1[k] * bc_vf_jj;
		accum_jm1[k] = 0.;
	      }
	    }
	    else if (eij == key_ij) { // bc_vf_jj is 1
	      for (k=0; k<numVars; ++k) {
		if (k != j) accum_j[k] += accum_jm1[k];
		accum_jm1[k] = 0.;
	      }
	    }
	    else // bc_vf_jj is 0
	      for (k=0; k<numVars; ++k)
		accum_jm1[k] = 0.;
	  }
	  else { // grad factor is zero, value factor is omitted
	    for (k=0; k<numVars; ++k) {
	      if (k != j) accum_j[k] += accum_jm1[k];
	      accum_jm1[k] = 0.;
	    }
	  }
	  if (key_ij + 1 != poly_j.interpolation_size())
	    break;
	}
      }
    }
    Real bcg_scale = barycentric_gradient_scaling(basis_index);
    Real* accum = accumulator[numVars-1];
    for (j=0; j<numVars; ++j)
      tpGradient[j] = accum[j] * bcg_scale;

    /*
    size_t i, j, num_colloc_pts = key.size();
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray&   key_i = key[i];
      Real t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      for (j=0; j<numVars; ++j)
	tpGradient[j] += t1_coeff_i *
	  barycentric_gradient_factor(j, key_i, basis_index);
    }
    tpGradient.scale(barycentric_gradient_scaling(basis_index));
    */
  }
  else if (exp_t2_coeffs.empty()) {
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray&   key_i = key[i];
      Real t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      for (j=0; j<numVars; ++j)
	tpGradient[j] += t1_coeff_i *
	  type1_interpolant_gradient(x, j, key_i, basis_index);
    }
  }
  else {
    size_t k;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray&   key_i = key[i];
      Real t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      const Real* exp_t2_coeff_i = (colloc_index.empty()) ?
	exp_t2_coeffs[i] : exp_t2_coeffs[colloc_index[i]];
      for (j=0; j<numVars; ++j) { // ith contribution to jth grad component
	tpGradient[j] += t1_coeff_i *
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
					const RealVector& subset_t1_coeffs,
					const RealMatrix& subset_t2_coeffs,
					const UShortArray& basis_index,
					const UShort2DArray& subset_key,
					const SizetArray& subset_colloc_index,
					const SizetList& subset_indices)
{
  if (tpGradient.length() != numVars)
    tpGradient.sizeUninitialized(numVars);
  tpGradient = 0.;
  size_t i, j, num_colloc_pts = subset_key.size();

  if (barycentricFlag) {

    // For barycentric interpolation: track x != newPoint within 1D basis
    set_new_point(x, basis_index, subset_indices, 3);//value+grad factors needed

    // Note: for basis gradients, one cannot eliminate exact index dimensions
    // as is managed in the barycentric case of tensor_product_value()

    SizetList::const_iterator cit = subset_indices.begin();
    size_t k, num_subset_v = subset_indices.size(), v0 = *cit, vj;
    unsigned short key_i0, key_ij, bi0 = basis_index[v0], bij;
    BasisPolynomial&   poly_v0 = polynomialBasis[bi0][v0];
    const RealVector& bc_vf_v0 = poly_v0.barycentric_value_factors();
    const RealVector& bc_gf_v0 = poly_v0.barycentric_gradient_factors();
    size_t ei0 = poly_v0.exact_index(), eij;
    unsigned short max0 = poly_v0.interpolation_size() - 1;
    RealMatrix accumulator(numVars, num_subset_v); // init to 0.
    Real *accum_0 = accumulator[0], *accum_j, *accum_jm1, bc_vf_00, bc_vf_jj,
      t1_coeff;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = subset_key[i]; key_i0 = key_i[v0];
      t1_coeff = (subset_colloc_index.empty()) ? subset_t1_coeffs[i] :
	subset_t1_coeffs[subset_colloc_index[i]];
      if (bi0) {
	// accumulate grad factor for comp 0 and value factor for other comps
	accum_0[0] += t1_coeff * bc_gf_v0[key_i0];
	if (ei0 == _NPOS) {
	  bc_vf_00 = bc_vf_v0[key_i0];
	  for (j=1; j<numVars; ++j)
	    accum_0[j] += t1_coeff * bc_vf_00;
	}
	else if (ei0 == key_i0) // value factor is 1, else value factor is 0
	  for (j=1; j<numVars; ++j)
	    accum_0[j] += t1_coeff;
      }
      else // grad factor is zero, value factor is omitted
	for (j=1; j<numVars; ++j)
	  accum_0[j] += t1_coeff;
      if (key_i0 == max0) {
	// accumulate sums over variables with max key value
	for (j=1, cit=++subset_indices.begin(); j<num_subset_v; ++j, ++cit) {
	  vj = *cit; bij = basis_index[vj]; key_ij  = key_i[vj];
	  accum_j = accumulator[j]; accum_jm1 = accumulator[j-1];
	  BasisPolynomial& poly_vj = polynomialBasis[basis_index[vj]][vj];
	  if (bij) {
	    eij = poly_vj.exact_index();
	    accum_j[j] += accum_jm1[j] *
	      poly_vj.barycentric_gradient_factor(key_ij);
	    if (eij == _NPOS) { // bc_vf_jj has general value
	      bc_vf_jj = poly_vj.barycentric_value_factor(key_ij);
	      for (k=0; k<numVars; ++k) {
		if (k != j) accum_j[k] += accum_jm1[k] * bc_vf_jj;
		accum_jm1[k] = 0.;
	      }
	    }
	    else if (eij == key_ij) { // bc_vf_jj is 1
	      for (k=0; k<numVars; ++k) {
		if (k != j) accum_j[k] += accum_jm1[k];
		accum_jm1[k] = 0.;
	      }
	    }
	    else // bc_vf_jj is 0
	      for (k=0; k<numVars; ++k)
		accum_jm1[k] = 0.;
	  }
	  else { // grad factor is zero, value factor is omitted
	    for (k=0; k<numVars; ++k) {
	      if (k != j) accum_j[k] += accum_jm1[k];
	      accum_jm1[k] = 0.;
	    }
	  }
	  if (key_ij + 1 != poly_vj.interpolation_size())
	    break;
	}
      }
    }
    Real bcg_scale = barycentric_gradient_scaling(basis_index, subset_indices);
    Real* accum = accumulator[num_subset_v-1];
    for (j=0; j<numVars; ++j)
      tpGradient[j] = accum[j] * bcg_scale;

    /*
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = subset_key[i];
      Real t1_coeff_i = (subset_colloc_index.empty()) ?
	subset_t1_coeffs[i] : subset_t1_coeffs[subset_colloc_index[i]];
      for (j=0; j<numVars; ++j)
	tpGradient[j] += t1_coeff_i *
	  barycentric_gradient_factor(j, key_i, basis_index, subset_indices);
    }
    tpGradient.scale(barycentric_gradient_scaling(basis_index, subset_indices));
    */
  }
  else if (subset_t2_coeffs.empty()) {
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = subset_key[i];
      Real t1_coeff_i = (subset_colloc_index.empty()) ?
	subset_t1_coeffs[i] : subset_t1_coeffs[subset_colloc_index[i]];
      for (j=0; j<numVars; ++j)
	tpGradient[j] += t1_coeff_i *
	  type1_interpolant_gradient(x, j, key_i, basis_index, subset_indices);
    }
  }
  else {
    size_t k;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = subset_key[i];
      Real t1_coeff_i = (subset_colloc_index.empty()) ?
	subset_t1_coeffs[i] : subset_t1_coeffs[subset_colloc_index[i]];
      const Real* t2_coeff_i = (subset_colloc_index.empty()) ?
	subset_t2_coeffs[i] : subset_t2_coeffs[subset_colloc_index[i]];
      for (j=0; j<numVars; ++j) { // ith contribution to jth grad component
	tpGradient[j] += t1_coeff_i *
	  type1_interpolant_gradient(x, j, key_i, basis_index, subset_indices);
	for (k=0; k<numVars; ++k) // type2 interpolant for kth grad comp
	  tpGradient[j] += t2_coeff_i[k] *
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
  if (!num_deriv_vars) return tpGradient;
  tpGradient = 0.;

  if (barycentricFlag) {

    // For barycentric interpolation: track x != newPoint within 1D basis
    set_new_point(x, basis_index, 3); // value+gradient factors needed

    // Note: for basis gradients, one cannot eliminate exact index dimensions
    // as is managed in the barycentric case of tensor_product_value()

    unsigned short key_i0, key_ij, bi0 = basis_index[0], bij;
    BasisPolynomial&   poly_0 = polynomialBasis[bi0][0];
    const RealVector& bc_vf_0 = poly_0.barycentric_value_factors();
    const RealVector& bc_gf_0 = poly_0.barycentric_gradient_factors();
    size_t k, ei0 = poly_0.exact_index(), eij, d0_index = dvv[0] - 1, start;
    unsigned short max0 = poly_0.interpolation_size() - 1;
    RealMatrix accumulator(num_deriv_vars, numVars); // init to 0.
    Real *accum_0 = accumulator[0], *accum_j, *accum_jm1,
      t1_coeff_bc_vf_0, bc_vf_j, bc_gf_j, t1_coeff;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i]; key_i0 = key_i[0];
      t1_coeff = (colloc_index.empty()) ? exp_t1_coeffs[i] :
	exp_t1_coeffs[colloc_index[i]];
      // mirror version w/o DVV, only logging a subset of accum components
      start = 0;
      if (bi0) {
	if (d0_index == 0)
	  { start = 1; accum_0[0] += t1_coeff * bc_gf_0[key_i0]; }
        if (ei0 == _NPOS) {
	  t1_coeff_bc_vf_0 = t1_coeff * bc_vf_0[key_i0];
	  for (j=start; j<num_deriv_vars; ++j)
	    accum_0[j] += t1_coeff_bc_vf_0;
	}
	else if (ei0 == key_i0) // value factor is 1, else value factor is 0
	  for (j=start; j<num_deriv_vars; ++j)
	    accum_0[j] += t1_coeff;
      }
      else { // grad factor is 0., value factor is omitted
	if (d0_index == 0) start = 1;
	for (j=start; j<num_deriv_vars; ++j)
	  accum_0[j] += t1_coeff;
      }
      if (key_i0 == max0) {
	// accumulate sums over variables with max key value
	for (j=1; j<numVars; ++j) {
	  bij = basis_index[j]; key_ij = key_i[j];
	  accum_j = accumulator[j]; accum_jm1 = accumulator[j-1];
	  BasisPolynomial& poly_j = polynomialBasis[bij][j];
	  if (bij) {
	    eij = poly_j.exact_index();
	    bc_gf_j = poly_j.barycentric_gradient_factor(key_ij);
	    if (eij == _NPOS) { // bc_vf_j has general value
	      bc_vf_j = poly_j.barycentric_value_factor(key_ij);
	      for (k=0; k<num_deriv_vars; ++k) {
		accum_j[k] += (j == dvv[k] - 1) ? accum_jm1[k] * bc_gf_j
		                                : accum_jm1[k] * bc_vf_j;
		accum_jm1[k] = 0.;
	      }
	    }
	    else if (eij == key_ij) // bc_vf_j is 1
	      for (k=0; k<num_deriv_vars; ++k) {
		accum_j[k] += (j == dvv[k] - 1) ? accum_jm1[k] * bc_gf_j
		                                : accum_jm1[k];
		accum_jm1[k] = 0.;
	      }
	    else // bc_vf_j is 0
	      for (k=0; k<num_deriv_vars; ++k) {
		if (j == dvv[k] - 1) accum_j[k] += accum_jm1[k] * bc_gf_j;
		accum_jm1[k] = 0.;
	      }
	  }
	  else { // grad factor is zero, value factor is omitted
	    for (k=0; k<num_deriv_vars; ++k) {
	      if (j != dvv[k] - 1) accum_j[k] += accum_jm1[k];
	      accum_jm1[k] = 0.;
	    }
	  }
	  if (key_ij + 1 != poly_j.interpolation_size())
	    break;
	}
      }
    }
    Real bcg_scale = barycentric_gradient_scaling(basis_index);
    Real* accum = accumulator[numVars-1];
    for (j=0; j<num_deriv_vars; ++j)
      tpGradient[j] = accum[j] * bcg_scale;

    /*
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray&   key_i = key[i];
      Real t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      for (j=0; j<num_deriv_vars; ++j) {
	deriv_index = dvv[j] - 1; // requires an "All" view
	tpGradient[j] += t1_coeff_i *
	  barycentric_gradient_factor(deriv_index, key_i, basis_index);
      }
    }
    tpGradient.scale(barycentric_gradient_scaling(basis_index));
    */
  }
  else if (exp_t2_coeffs.empty()) {
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray&   key_i = key[i];
      Real t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      for (j=0; j<num_deriv_vars; ++j) {
	deriv_index = dvv[j] - 1; // requires an "All" view
	tpGradient[j] += t1_coeff_i *
	  type1_interpolant_gradient(x, deriv_index, key_i, basis_index);
      }
    }
  }
  else {
    size_t k;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray&   key_i = key[i];
      Real t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      const Real* exp_t2_coeff_i = (colloc_index.empty()) ?
	exp_t2_coeffs[i] : exp_t2_coeffs[colloc_index[i]];
      for (j=0; j<num_deriv_vars; ++j) {
	deriv_index = dvv[j] - 1; // requires an "All" view
	tpGradient[j] += t1_coeff_i *
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
  size_t num_deriv_vars = exp_t1_coeff_grads.numRows();
  if (tpGradient.length() != num_deriv_vars)
    tpGradient.sizeUninitialized(num_deriv_vars);
  tpGradient = 0.;
  if (barycentricFlag) {

    // For barycentric interpolation: track x != newPoint within 1D basis
    set_new_point(x, basis_index, 1); // value factors needed for coeff grads

    size_t num_act_v = barycentric_active_variables(basis_index);
    if (num_act_v == 0) { // convert 1-D exact indices into n-D colloc index
      size_t pt_index = barycentric_exact_index(basis_index);
      if (colloc_index.empty())
	copy_data(exp_t1_coeff_grads[pt_index],(int)num_deriv_vars, tpGradient);
      else
	copy_data(exp_t1_coeff_grads[colloc_index[pt_index]],
		  (int)num_deriv_vars, tpGradient);
    }
    else if (num_act_v == numVars) { // interpolation over all variables
      RealMatrix accumulator(numVars, num_deriv_vars); // init to 0.
      BasisPolynomial&   poly_0 = polynomialBasis[basis_index[0]][0];
      const RealVector& bc_vf_0 = poly_0.barycentric_value_factors();
      size_t i, j, k, num_colloc_pts = key.size();
      unsigned short key_i0, key_ij, max0 = poly_0.interpolation_size() - 1;
      Real *accum_0, *accum_j, *accum_jm1, bc_vf_00, bc_vf_jj;
      for (i=0; i<num_colloc_pts; ++i) {
	const UShortArray& key_i = key[i]; key_i0 = key_i[0];
	const Real* grad = (colloc_index.empty()) ? exp_t1_coeff_grads[i] :
	  exp_t1_coeff_grads[colloc_index[i]];
	accum_0 = accumulator[0]; bc_vf_00 = bc_vf_0[key_i0];
	for (j=0; j<num_deriv_vars; ++j)
	  accum_0[j] += grad[j] * bc_vf_00;
	if (key_i0 == max0) {
	  // accumulate sums over variables with max key value
	  for (j=1; j<numVars; ++j) {
	    BasisPolynomial& poly_j = polynomialBasis[basis_index[j]][j];
	    key_ij = key_i[j];
	    accum_j = accumulator[j]; accum_jm1 = accumulator[j-1];
	    bc_vf_jj = poly_j.barycentric_value_factor(key_ij);
	    for (k=0; k<num_deriv_vars; ++k) {
	      accum_j[k]  += accum_jm1[k] * bc_vf_jj;
	      accum_jm1[k] = 0.;
	    }
	    if (key_ij + 1 != poly_j.interpolation_size())
	      break;
	  }
	}
      }
      Real scale = 1. / barycentric_value_scaling(basis_index);
      Real* accum = accumulator[numVars-1];
      for (j=0; j<num_deriv_vars; ++j)
	tpGradient[j] = accum[j] * scale;
    }
    else { // partial interpolation over active variables
      SizetArray pt_factors(num_act_v), act_v_set(num_act_v);
      unsigned short bi_j;
      size_t j, pts_vj, num_act_pts = 1, num_pts = 1, ej, av_cntr, pt_index = 0;
      for (j=0, av_cntr=0; j<numVars; ++j) {
	bi_j = basis_index[j];
	if (bi_j) { // else pts_vj = 1 and ej can be taken to be 0
	  BasisPolynomial& poly_j = polynomialBasis[bi_j][j];
	  ej = poly_j.exact_index(); pts_vj = poly_j.interpolation_size();
	  if (ej == _NPOS) { // active for interpolation
	    pt_factors[av_cntr] = num_pts; act_v_set[av_cntr] = j;
	    num_act_pts *= pts_vj; ++av_cntr;
	  }
	  else             // inactive for interpolation
	    pt_index += num_pts * ej;
	  num_pts *= pts_vj;
	}
      }
      // define initial pt_index offset
      RealMatrix accumulator(num_deriv_vars, num_act_v); // init to 0.
      size_t i, k, v0 = act_v_set[0], vj;
      BasisPolynomial&   poly_v0 = polynomialBasis[basis_index[v0]][v0];
      const RealVector& bc_vf_v0 = poly_v0.barycentric_value_factors();
      size_t pts_v0 = poly_v0.interpolation_size(), pf0 = pt_factors[0],
	 pts_v0_pf0 = pts_v0 * pf0, pfj, prev_pt_set;
      unsigned short key_i0, key_ij, max0 = pts_v0 - 1;
      Real *accum_0, *accum_j, *accum_jm1, bc_vf_00, bc_vf_jj;
      // loop over active pts, summing contributions from active variables
      for (i=0; i<num_act_pts; ++i) {
	const UShortArray& key_i = key[pt_index]; key_i0 = key_i[v0];
	const Real* grad = (colloc_index.empty()) ? exp_t1_coeff_grads[pt_index]
	  : exp_t1_coeff_grads[colloc_index[pt_index]];
	accum_0 = accumulator[0]; bc_vf_00 = bc_vf_v0[key_i0];
	for (j=0; j<num_deriv_vars; ++j)
	  accum_0[j] += grad[j] * bc_vf_00;
	pt_index += pf0;
	if (key_i0 == max0) {
	  // accumulate sums over variables with max key value
	  for (j=1, prev_pt_set=pts_v0_pf0; j<num_act_v; ++j) {
	    vj = act_v_set[j]; key_ij = key_i[vj]; pfj = pt_factors[j];
	    BasisPolynomial& poly_vj = polynomialBasis[basis_index[vj]][vj];
	    // update accumulators: push [j-1] entry up to [j] level
	    accum_j = accumulator[j]; accum_jm1 = accumulator[j-1];
	    bc_vf_jj = poly_vj.barycentric_value_factor(key_ij);
	    for (k=0; k<num_deriv_vars; ++k) {
	      accum_j[k]  += accum_jm1[k] * bc_vf_jj;
	      accum_jm1[k] = 0.;
	    }
	    // update pt_index: prev index rolls back to 0 and curr index +1.
	    // index increment is zero unless active vars are nonconsecutive.
	    if (pfj != prev_pt_set)
	      pt_index += pfj - prev_pt_set;
	    pts_vj = poly_vj.interpolation_size();
	    if (key_ij + 1 == pts_vj) prev_pt_set = pts_vj * pfj;
	    else                      break;
	  }
	}
      }
      Real scale = 1. / barycentric_value_scaling(basis_index);
      Real* accum = accumulator[num_act_v-1];
      for (j=0; j<num_deriv_vars; ++j)
	tpGradient[j] = accum[j] * scale;
    }

    /*
    size_t i, num_colloc_pts = key.size();
    for (i=0; i<num_colloc_pts; ++i) {
      const Real* t1_coeff_grad_i = (colloc_index.empty()) ?
	exp_t1_coeff_grads[i] : exp_t1_coeff_grads[colloc_index[i]];
      Real bc_fact = barycentric_value_factor(key[i], basis_index);
      for (j=0; j<num_deriv_vars; ++j)
	tpGradient[j] += t1_coeff_grad_i[j] * bc_fact;
    }
    // apply barycentric denominator
    tpGradient.scale(1. / barycentric_value_scaling(basis_index));
    */
  }
  else {
    size_t i, j, num_colloc_pts = key.size(); Real t1_val;
    for (i=0; i<num_colloc_pts; ++i) {
      const Real* t1_coeff_grad_i = (colloc_index.empty()) ?
	exp_t1_coeff_grads[i] : exp_t1_coeff_grads[colloc_index[i]];
      t1_val = type1_interpolant_value(x, key[i], basis_index);
      for (j=0; j<num_deriv_vars; ++j)
	tpGradient[j] += t1_coeff_grad_i[j] * t1_val;
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

/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HierarchInterpPolyApproximation
//- Description:  Implementation code for InterpPolyApproximation class
//-               
//- Owner:        Mike Eldred

#include "HierarchInterpPolyApproximation.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"
#include "pecos_stat_util.hpp"

//#define DEBUG
//#define VBD_DEBUG
//#define INTERPOLATION_TEST

namespace Pecos {


void HierarchInterpPolyApproximation::allocate_arrays()
{
  InterpPolyApproximation::allocate_arrays();

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  const ExpansionConfigOptions& ec_options = data_rep->expConfigOptions;
  const UShort4DArray& colloc_key = data_rep->hsg_driver()->collocation_key();
  size_t i, j, k, num_levels = colloc_key.size(), num_sets, num_tp_pts,
    num_deriv_v = surrData.num_derivative_variables();

  RealVector2DArray& exp_t1_coeffs = expT1CoeffsIter->second;
  RealMatrix2DArray& exp_t2_coeffs = expT2CoeffsIter->second;
  RealMatrix2DArray& exp_t1_coeff_grads = expT1CoeffGradsIter->second;
  
  if (exp_t1_coeffs.size() != num_levels)
    exp_t1_coeffs.resize(num_levels);
  if (exp_t2_coeffs.size() != num_levels)
    exp_t2_coeffs.resize(num_levels);
  if (exp_t1_coeff_grads.size() != num_levels)
    exp_t1_coeff_grads.resize(num_levels);
  for (i=0; i<num_levels; ++i) {
    const UShort3DArray& key_i = colloc_key[i];
    num_sets = key_i.size();
    if (exp_t1_coeffs[i].size() != num_sets)
      exp_t1_coeffs[i].resize(num_sets);
    if (exp_t2_coeffs[i].size() != num_sets)
      exp_t2_coeffs[i].resize(num_sets);
    if (exp_t1_coeff_grads[i].size() != num_sets)
      exp_t1_coeff_grads[i].resize(num_sets);
    for (j=0; j<num_sets; ++j) {
      num_tp_pts = key_i[j].size();
      for (k=0; k<num_tp_pts; ++k) {
	if (expansionCoeffFlag) {
	  exp_t1_coeffs[i][j].sizeUninitialized(num_tp_pts);
	  if (data_rep->basisConfigOptions.useDerivs)
	    exp_t2_coeffs[i][j].shapeUninitialized(num_deriv_v, num_tp_pts);
	}
	if (expansionCoeffGradFlag)
	  exp_t1_coeff_grads[i][j].shapeUninitialized(num_deriv_v, num_tp_pts);
      }
    }
  }

  // checking num_points is insufficient due to anisotropy --> changes in
  // anisotropic weights could move points around without changing the total
  //size_t num_points = surrData.points();
  //bool update_exp_form =
  //  ( (expansionCoeffFlag     && exp_t1_coeffs.length() != num_points) ||
  //    (expansionCoeffGradFlag && exp_t1_coeff_grads.numCols() != num_points));

  if (ec_options.refinementControl) {
    size_t num_moments = (data_rep->nonRandomIndices.empty()) ? 4 : 2;
    if (referenceMoments.empty())
      referenceMoments.sizeUninitialized(num_moments);
    if (deltaMoments.empty())
      deltaMoments.sizeUninitialized(num_moments);
  }
}


void HierarchInterpPolyApproximation::compute_coefficients()
{
  PolynomialApproximation::compute_coefficients();
  if (!expansionCoeffFlag && !expansionCoeffGradFlag)
    return;

  allocate_arrays();

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver   = data_rep->hsg_driver();
  const UShort3DArray&      sm_mi        = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      colloc_key   = hsg_driver->collocation_key();
  const Sizet3DArray&       colloc_index = hsg_driver->collocation_indices();
  size_t lev, set, pt, v, num_levels = colloc_key.size(), num_sets, num_tp_pts,
    cntr = 0, c_index, num_deriv_vars = surrData.num_derivative_variables();

  RealVector2DArray& exp_t1_coeffs = expT1CoeffsIter->second;
  RealMatrix2DArray& exp_t2_coeffs = expT2CoeffsIter->second;
  RealMatrix2DArray& exp_t1_coeff_grads = expT1CoeffGradsIter->second;
  const SDVArray& sdv_array = surrData.variables_data();
  const SDRArray& sdr_array = surrData.response_data();

  // level 0
  c_index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
  const SurrogateDataResp& sdr_0 = sdr_array[c_index];
  if (expansionCoeffFlag) {
    exp_t1_coeffs[0][0][0] = sdr_0.response_function();
    if (data_rep->basisConfigOptions.useDerivs)
      Teuchos::setCol(sdr_0.response_gradient(), 0, exp_t2_coeffs[0][0]);
  }
  if (expansionCoeffGradFlag)
    Teuchos::setCol(sdr_0.response_gradient(), 0, exp_t1_coeff_grads[0][0]);
  ++cntr;
  // levels 1 to num_levels
  for (lev=1; lev<num_levels; ++lev) {
    const UShort3DArray& key_l = colloc_key[lev];
    num_sets = key_l.size();
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = key_l[set].size();
      for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	c_index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	const RealVector& c_vars = sdv_array[c_index].continuous_variables();
	const SurrogateDataResp& sdr = sdr_array[c_index];
	// coefficients are hierarchical surpluses
	if (expansionCoeffFlag) {
	  exp_t1_coeffs[lev][set][pt] = sdr.response_function() -
	    value(c_vars, sm_mi, colloc_key, exp_t1_coeffs,exp_t2_coeffs,lev-1);
	  if (data_rep->basisConfigOptions.useDerivs) {
	    const RealVector& data_grad = sdr.response_gradient();
	    const RealVector& prev_grad = gradient_basis_variables(c_vars,
	      sm_mi, colloc_key, exp_t1_coeffs, exp_t2_coeffs, lev-1);
	    Real* hier_grad = exp_t2_coeffs[lev][set][pt];
	    for (v=0; v<num_deriv_vars; ++v)
	      hier_grad[v] = data_grad[v] - prev_grad[v];
	  }
	}
	if (expansionCoeffGradFlag) {
	  const RealVector& data_grad = sdr.response_gradient();
	  const RealVector& prev_grad = gradient_nonbasis_variables(c_vars,
	    sm_mi, colloc_key, exp_t1_coeff_grads, lev-1);
	  Real* hier_grad = exp_t1_coeff_grads[lev][set][pt];
	  for (v=0; v<num_deriv_vars; ++v)
	    hier_grad[v] = data_grad[v] - prev_grad[v];
	}
      }
    }
  }

#ifdef INTERPOLATION_TEST
  test_interpolation();
#endif

  clear_all_computed_bits();
}


void HierarchInterpPolyApproximation::increment_coefficients()
{
  // TO DO: partial sync for new TP data set, e.g. update_surrogate_data() ?
  synchronize_surrogate_data();

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();

  bool updated = update_active_iterators(data_rep->activeKey);
  if (updated) clear_all_computed_bits();
  else         increment_current_from_reference();

  switch (data_rep->expConfigOptions.refinementControl) {
  case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED: // generalized sparse grids
    increment_coefficients(hsg_driver->trial_set());
    break;
  default: {
    const UShort3DArray&   sm_mi = hsg_driver->smolyak_multi_index();
    const UShortArray& incr_sets = hsg_driver->increment_sets();
    size_t lev, num_lev = sm_mi.size(), set, start_set, num_sets;
    for (lev=0; lev<num_lev; ++lev) {
      start_set = incr_sets[lev]; num_sets = sm_mi[lev].size();
      for (set=start_set; set<num_sets; ++set)
	increment_coefficients(sm_mi[lev][set]);
    }
    break;
  }
  }

  // size sobolIndices based on shared sobolIndexMap
  allocate_component_sobol();
}


// ***************************************************************
// TO DO: verify that decrement/push is always valid for surpluses
// ***************************************************************


void HierarchInterpPolyApproximation::decrement_coefficients(bool save_data)
{
  // mirror changes to origSurrData for deep copied surrData
  if (deep_copied_surrogate_data())
    surrData.pop(save_data);

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  const UShortArray& key = data_rep->activeKey;

  bool updated = update_active_iterators(key);
  if (updated) clear_all_computed_bits();
  else         decrement_current_to_reference();

  const UShortArray& trial_set = data_rep->hsg_driver()->trial_set();
  size_t lev = l1_norm(trial_set);

  if (expansionCoeffFlag) {
    RealVector2DArray& exp_t1c = expT1CoeffsIter->second;
    poppedExpT1Coeffs[key][trial_set] = exp_t1c[lev].back();
    exp_t1c[lev].pop_back();
    if (data_rep->basisConfigOptions.useDerivs) {
      RealMatrix2DArray& exp_t2c = expT2CoeffsIter->second;
      poppedExpT2Coeffs[key][trial_set] = exp_t2c[lev].back();
      exp_t2c[lev].pop_back();
    }
  }
  if (expansionCoeffGradFlag) {
    RealMatrix2DArray& exp_t1cg = expT1CoeffGradsIter->second;
    poppedExpT1CoeffGrads[key][trial_set] = exp_t1cg[lev].back();
    exp_t1cg[lev].pop_back();
  }
}


void HierarchInterpPolyApproximation::push_coefficients()
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;

  // mirror changes to origSurrData for deep copied surrData
  if (deep_copied_surrogate_data())
    surrData.push(data_rep->retrieval_index());

  bool updated = update_active_iterators(data_rep->activeKey);
  if (updated) clear_all_computed_bits();
  else         increment_current_from_reference();

  push_coefficients(data_rep->hsg_driver()->trial_set());
}


void HierarchInterpPolyApproximation::finalize_coefficients()
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  const UShortArray& key = data_rep->activeKey;

  // mirror changes to origSurrData for deep copied surrData
  if (deep_copied_surrogate_data()) {
    size_t i, num_popped = surrData.popped_sets(); // # of popped trials
    for (i=0; i<num_popped; ++i)
      surrData.push(data_rep->finalization_index(i), false);
    surrData.clear_active_popped(); // only after process completed
  }

  // synchronize expansionCoeff{s,Grads} and approxData
  bool updated = update_active_iterators(key);
  if (updated) clear_all_computed_bits();
  else         clear_computed_bits(); //clear_delta_computed_bits();

  const UShort3DArray& sm_mi = data_rep->hsg_driver()->smolyak_multi_index();
  size_t lev, set, num_sets, num_levels = sm_mi.size(), num_smolyak_sets,
    num_coeff_sets;
  const RealVector2DArray& exp_t1_coeffs      = expT1CoeffsIter->second;
  const RealMatrix2DArray& exp_t1_coeff_grads = expT1CoeffGradsIter->second;
  for (lev=0; lev<num_levels; ++lev) {
    const UShort2DArray& sm_mi_l = sm_mi[lev];
    num_smolyak_sets = sm_mi_l.size();
    num_coeff_sets = (expansionCoeffFlag) ? exp_t1_coeffs[lev].size() :
      exp_t1_coeff_grads[lev].size();
    for (set=num_coeff_sets; set<num_smolyak_sets; ++set)
      push_coefficients(sm_mi_l[set]);
  }
  poppedExpT1Coeffs[key].clear(); poppedExpT2Coeffs[key].clear();
  poppedExpT1CoeffGrads[key].clear();
}


void HierarchInterpPolyApproximation::clear_inactive()
{
  // mirror changes to origSurrData for deep copied surrData
  if (deep_copied_surrogate_data())
    surrData.clear_inactive();

  std::map<UShortArray, RealVector2DArray>::iterator e1c_it
    = expansionType1Coeffs.begin();
  std::map<UShortArray, RealMatrix2DArray>::iterator e2c_it
    = expansionType2Coeffs.begin();
  std::map<UShortArray, RealMatrix2DArray>::iterator e1g_it
    = expansionType1CoeffGrads.begin();
  while (e1c_it != expansionType1Coeffs.end())
    if (e1c_it == expT1CoeffsIter) // preserve active
      { ++e1c_it; ++e2c_it; ++e1g_it; }
    else { // clear inactive: postfix increments manage iterator invalidations
      expansionType1Coeffs.erase(e1c_it++);
      expansionType2Coeffs.erase(e2c_it++);
      expansionType1CoeffGrads.erase(e1g_it++);
    }
}


void HierarchInterpPolyApproximation::combine_coefficients()
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  if (deep_copied_surrogate_data())
    surrData.active_key(data_rep->activeKey);

  // Coefficient combination is not dependent on active state
  //bool updated = update_active_iterators(data_rep->activeKey);
  //if (updated) clear_all_computed_bits();
  //else (combined stats managed separately)

  allocate_component_sobol(); // size sobolIndices from shared sobolIndexMap

  const UShort4DArray&   comb_key = data_rep->combinedCollocKey;
  const Sizet3DArray& comb_sm_map = data_rep->combinedSmolyakMultiIndexMap;
  size_t i, lev, set, pt, num_lev = comb_key.size(), num_sets, num_tp_pts,
    num_v = surrData.num_derivative_variables();
  bool use_derivs = data_rep->basisConfigOptions.useDerivs;

  // Resize combined expansion arrays and initialize to zero.
  // Note: STL resize() is no-op if already correct size, but Teuchos
  // resize()/reshape() is not; therefore protect the latter.
  combinedExpT1Coeffs.resize(num_lev); combinedExpT2Coeffs.resize(num_lev);
  combinedExpT1CoeffGrads.resize(num_lev);
  for (lev=0; lev<num_lev; ++lev) {
    const UShort3DArray& comb_key_l = comb_key[lev];
    num_sets = comb_key_l.size();
    RealVectorArray& comb_t1c_l = combinedExpT1Coeffs[lev];
    RealMatrixArray& comb_t2c_l = combinedExpT2Coeffs[lev];
    RealMatrixArray& comb_t1g_l = combinedExpT1CoeffGrads[lev];
    comb_t1c_l.resize(num_sets);  comb_t2c_l.resize(num_sets);
    comb_t1g_l.resize(num_sets);
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = comb_key_l[set].size();
      RealVector& comb_t1c_ls = comb_t1c_l[set];
      RealMatrix& comb_t2c_ls = comb_t2c_l[set];
      RealMatrix& comb_t1g_ls = comb_t1g_l[set];
      if (expansionCoeffFlag) {
	if (comb_t1c_ls.length() != num_tp_pts)
	  comb_t1c_ls.size(num_tp_pts);           // init to 0
	else comb_t1c_ls = 0.;
      }
      else if (!comb_t1c_ls.empty())
	comb_t1c_ls.size(0);
      if (use_derivs && expansionCoeffFlag) {
	if (comb_t2c_ls.numRows()!=num_v || comb_t2c_ls.numCols()!=num_tp_pts)
	  comb_t2c_ls.shape(num_v, num_tp_pts); // init to 0
	else comb_t2c_ls = 0.;
      }
      else if (!comb_t2c_ls.empty())
	comb_t2c_ls.shape(0, 0);
      if (expansionCoeffGradFlag) {
	if (comb_t1g_ls.numRows()!=num_v || comb_t1g_ls.numCols()!=num_tp_pts)
	  comb_t1g_ls.shape(num_v, num_tp_pts);   // init to 0
	else comb_t1g_ls = 0.;
      }
      else if (!comb_t1g_ls.empty())
	comb_t1g_ls.shape(0, 0);
    }
  }

  // Roll up hierarchical surplus increments for each level-set-point
  if (expansionCoeffFlag) {
    std::map<UShortArray, RealVector2DArray>::iterator ec1_it;
    std::map<UShortArray, RealMatrix2DArray>::iterator ec2_it;
    for (ec1_it  = expansionType1Coeffs.begin(),
	 ec2_it  = expansionType2Coeffs.begin(), i=0;
	 ec1_it != expansionType1Coeffs.end() &&
	 ec2_it != expansionType2Coeffs.end(); ++ec1_it, ++ec2_it, ++i) {
      const RealVector2DArray& t1c = ec1_it->second;
      const RealMatrix2DArray& t2c = ec2_it->second;
      num_lev = t1c.size();
      for (lev=0; lev<num_lev; ++lev) {
	const RealVectorArray& t1c_l = t1c[lev];
	const RealMatrixArray& t2c_l = t2c[lev];
	RealVectorArray&  comb_t1c_l = combinedExpT1Coeffs[lev];
	RealMatrixArray&  comb_t2c_l = combinedExpT2Coeffs[lev];
	const SizetArray& comb_sm_map_il = comb_sm_map[i][lev];
	num_sets = t1c_l.size();
	for (set=0; set<num_sets; ++set) {
	  // map from smolyakMultiIndex for this key to corresponding set
	  // in combinedSmolyakMultiIndex
	  size_t comb_set = comb_sm_map_il[set];
	  comb_t1c_l[comb_set]                 += t1c_l[set];
	  if (use_derivs) comb_t2c_l[comb_set] += t2c_l[set];
	}
      }
    }
  }
  if (expansionCoeffGradFlag) { // sum up expansionType1CoeffGrads
    std::map<UShortArray, RealMatrix2DArray>::iterator eg1_it;
    for (eg1_it  = expansionType1CoeffGrads.begin(), i=0;
	 eg1_it != expansionType1CoeffGrads.end(); ++eg1_it, ++i) {
      const RealMatrix2DArray& t1g = eg1_it->second;
      num_lev = t1g.size();
      for (lev=0; lev<num_lev; ++lev) {
	const RealMatrixArray& t1g_l = t1g[lev];
	RealMatrixArray&  comb_t1g_l = combinedExpT1CoeffGrads[lev];
	const SizetArray& comb_sm_map_il = comb_sm_map[i][lev];
	num_sets = t1g_l.size();
	for (set=0; set<num_sets; ++set) {
	  // map from smolyakMultiIndex for this key to corresponding set
	  // in combinedSmolyakMultiIndex
	  size_t comb_set = comb_sm_map_il[set];
	  comb_t1g_l[comb_set] += t1g_l[set];
	}
      }
    }
  }
}


void HierarchInterpPolyApproximation::combined_to_active()
{
  // replace active expansions with combined expansion arrays
  // Note: clear_inactive() takes care of the auxilliary inactive expansions
  //       that are now assimilated within the active expansion

  if (expansionCoeffFlag) {
    std::swap(expT1CoeffsIter->second, combinedExpT1Coeffs);
    combinedExpT1Coeffs.clear();
    SharedHierarchInterpPolyApproxData* data_rep
      = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
    if (data_rep->basisConfigOptions.useDerivs) {
      std::swap(expT2CoeffsIter->second, combinedExpT2Coeffs);
      combinedExpT2Coeffs.clear();
    }
  }
  if (expansionCoeffGradFlag) {
    std::swap(expT1CoeffGradsIter->second, combinedExpT1CoeffGrads);
    combinedExpT1CoeffGrads.clear();
  }

  clear_all_computed_bits();
}


/** Lower level helper function to process a single index set. */
void HierarchInterpPolyApproximation::
increment_coefficients(const UShortArray& index_set)
{
  RealVector2DArray& exp_t1_coeffs = expT1CoeffsIter->second;
  RealMatrix2DArray& exp_t2_coeffs = expT2CoeffsIter->second;
  RealMatrix2DArray& exp_t1_coeff_grads = expT1CoeffGradsIter->second;
  const SDVArray& sdv_array = surrData.variables_data();
  const SDRArray& sdr_array = surrData.response_data();

  size_t lev, old_levels = exp_t1_coeffs.size(), set, old_sets,
    pt, old_pts = 0;
  for (lev=0; lev<old_levels; ++lev) {
    old_sets = exp_t1_coeffs[lev].size();
    for (set=0; set<old_sets; ++set)
      old_pts += (expansionCoeffFlag) ?	exp_t1_coeffs[lev][set].length() :
	exp_t1_coeff_grads[lev][set].numCols();
  }
  lev = l1_norm(index_set);
  if (lev >= old_levels) {
    exp_t1_coeffs.resize(lev+1);
    exp_t2_coeffs.resize(lev+1);
    exp_t1_coeff_grads.resize(lev+1);
  }
  set = exp_t1_coeffs[lev].size();
  // append empty and update in place
  RealVector fns; RealMatrix grads;
  exp_t1_coeffs[lev].push_back(fns);
  exp_t2_coeffs[lev].push_back(grads);
  exp_t1_coeff_grads[lev].push_back(grads);
  RealVector& t1_coeffs      = exp_t1_coeffs[lev][set];
  RealMatrix& t2_coeffs      = exp_t2_coeffs[lev][set];
  RealMatrix& t1_coeff_grads = exp_t1_coeff_grads[lev][set];

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  const UShort4DArray& colloc_key = hsg_driver->collocation_key();
  size_t index, num_trial_pts = colloc_key[lev][set].size(),
    v, num_deriv_vars = 0;
  if (expansionCoeffFlag) {
    t1_coeffs.sizeUninitialized(num_trial_pts);
    if (data_rep->basisConfigOptions.useDerivs) {
      num_deriv_vars = exp_t2_coeffs[0][0].numRows();
      t2_coeffs.shapeUninitialized(num_deriv_vars, num_trial_pts);
    }
  }
  if (expansionCoeffGradFlag) {
    num_deriv_vars = exp_t1_coeff_grads[0][0].numRows();
    t1_coeff_grads.shapeUninitialized(num_deriv_vars, num_trial_pts);
  }
 
  for (pt=0, index=old_pts; pt<num_trial_pts; ++pt, ++index) {
    const RealVector& c_vars = sdv_array[index].continuous_variables();
    const SurrogateDataResp& sdr = sdr_array[index];
    if (expansionCoeffFlag) {
      t1_coeffs[pt] = sdr.response_function() -
	value(c_vars, sm_mi, colloc_key, exp_t1_coeffs, exp_t2_coeffs, lev-1);
      if (data_rep->basisConfigOptions.useDerivs) {
	const RealVector& data_grad = sdr.response_gradient();
	const RealVector& prev_grad = gradient_basis_variables(c_vars, sm_mi,
	  colloc_key, exp_t1_coeffs, exp_t2_coeffs, lev-1);
	Real* hier_grad = t2_coeffs[pt];
	for (v=0; v<num_deriv_vars; ++v)
	  hier_grad[v] = data_grad[v] - prev_grad[v];
      }
    }
    if (expansionCoeffGradFlag) {
      const RealVector& data_grad = sdr.response_gradient();
      const RealVector& prev_grad = gradient_nonbasis_variables(c_vars, sm_mi,
	colloc_key, exp_t1_coeff_grads, lev-1);
      Real* hier_grad = t1_coeff_grads[pt];
      for (v=0; v<num_deriv_vars; ++v)
	hier_grad[v] = data_grad[v] - prev_grad[v];
    }
  }
}


/** Lower level helper function to process a single index set. */
void HierarchInterpPolyApproximation::
push_coefficients(const UShortArray& push_set)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  const UShortArray& key = data_rep->activeKey;
  size_t lev = l1_norm(push_set);
  if (expansionCoeffFlag) {
    std::map<UShortArray, RealVector>& pop_t1c = poppedExpT1Coeffs[key];
    std::map<UShortArray, RealVector>::iterator itv = pop_t1c.find(push_set);
    expT1CoeffsIter->second[lev].push_back(itv->second);
    pop_t1c.erase(itv);
    SharedHierarchInterpPolyApproxData* data_rep
      = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
    if (data_rep->basisConfigOptions.useDerivs) {
      std::map<UShortArray, RealMatrix>& pop_t2c = poppedExpT2Coeffs[key];
      std::map<UShortArray, RealMatrix>::iterator itm = pop_t2c.find(push_set);
      expT2CoeffsIter->second[lev].push_back(itm->second);
      pop_t2c.erase(itm);
    }
  }
  if (expansionCoeffGradFlag) {
    std::map<UShortArray, RealMatrix>& pop_t1cg = poppedExpT1CoeffGrads[key];
    std::map<UShortArray, RealMatrix>::iterator itm = pop_t1cg.find(push_set);
    expT1CoeffGradsIter->second[lev].push_back(itm->second);
    pop_t1cg.erase(itm);
  }
}


Real HierarchInterpPolyApproximation::
value(const RealVector& x, const UShort3DArray& sm_mi,
      const UShort4DArray& colloc_key, const RealVector2DArray& t1_coeffs,
      const RealMatrix2DArray& t2_coeffs, unsigned short level,
      const UShort2DArray& set_partition)
{
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::value()" << std::endl;
    abort_handler(-1);
  }

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  Real approx_val = 0.;
  SizetArray colloc_index; // empty -> 2DArrays allow default indexing
  size_t lev, set, set_start = 0, set_end;
  bool partial = !set_partition.empty();
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = colloc_key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    if (partial)
      { set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
    else // sm_mi and key include all current index sets, whereas the t1/t2
         // coeffs may refect a partial state derived from a ref_key
      set_end = t1_coeffs_l.size();
    for (set=set_start; set<set_end; ++set)
      approx_val +=
	data_rep->tensor_product_value(x, t1_coeffs_l[set], t2_coeffs_l[set],
				       sm_mi_l[set], key_l[set], colloc_index);
  }
  return approx_val;
}


/** All variables version. */
Real HierarchInterpPolyApproximation::
value(const RealVector& x, const UShort3DArray& sm_mi,
      const UShort4DArray& colloc_key, const RealVector2DArray& t1_coeffs,
      const RealMatrix2DArray& t2_coeffs, unsigned short level,
      const SizetList& subset_indices, const UShort2DArray& set_partition)
{
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::value()" << std::endl;
    abort_handler(-1);
  }

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  Real approx_val = 0.;
  SizetArray colloc_index; // empty -> 2DArrays allow default indexing
  size_t lev, set, set_start = 0, set_end;
  bool partial = !set_partition.empty();
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = colloc_key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    if (partial)
      { set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
    else // sm_mi and key include all current index sets, whereas the t1/t2
         // coeffs may refect a partial state derived from a ref_key
      set_end = t1_coeffs_l.size();
    for (set=set_start; set<set_end; ++set)
      approx_val +=
	data_rep->tensor_product_value(x, t1_coeffs_l[set], t2_coeffs_l[set],
				       sm_mi_l[set], key_l[set], colloc_index,
				       subset_indices);
  }
  return approx_val;
}


const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const UShort3DArray& sm_mi,
			 const UShort4DArray& colloc_key,
			 const RealVector2DArray& t1_coeffs,
			 const RealMatrix2DArray& t2_coeffs,
			 unsigned short level,
			 const UShort2DArray& set_partition)
{
  // this could define a default_dvv and call gradient_basis_variables(x, dvv),
  // but we want this fn to be as fast as possible

  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in HierarchInterpPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  size_t num_v = sharedDataRep->numVars;
  if (approxGradient.length() != num_v)
    approxGradient.sizeUninitialized(num_v);
  approxGradient = 0.;

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  SizetArray colloc_index; // empty -> 2DArrays allow default indexing
  size_t lev, set, set_start = 0, set_end;
  bool partial = !set_partition.empty();
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = colloc_key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    if (partial)
      { set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
    else // sm_mi and key include all current index sets, whereas the t1/t2
         // coeffs may refect a partial state derived from a ref_key
      set_end = t1_coeffs_l.size();
    for (set=set_start; set<set_end; ++set)
      approxGradient +=
	data_rep->tensor_product_gradient_basis_variables(x, t1_coeffs_l[set],
	  t2_coeffs_l[set], sm_mi_l[set], key_l[set], colloc_index);
  }

  return approxGradient;
}


const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const UShort3DArray& sm_mi,
			 const UShort4DArray& colloc_key,
			 const RealVector2DArray& t1_coeffs,
			 const RealMatrix2DArray& t2_coeffs,
			 unsigned short level, const SizetList& subset_indices,
			 const UShort2DArray& set_partition)
{
  // this could define a default_dvv and call gradient_basis_variables(x, dvv),
  // but we want this fn to be as fast as possible

  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in HierarchInterpPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  size_t num_v = sharedDataRep->numVars;
  if (approxGradient.length() != num_v)
    approxGradient.sizeUninitialized(num_v);
  approxGradient = 0.;

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  SizetArray colloc_index; // empty -> 2DArrays allow default indexing
  size_t lev, set, set_start = 0, set_end;
  bool partial = !set_partition.empty();
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = colloc_key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    if (partial)
      { set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
    else // sm_mi and key include all current index sets, whereas the t1/t2
         // coeffs may refect a partial state derived from a ref_key
      set_end = t1_coeffs_l.size();
    for (set=set_start; set<set_end; ++set)
      approxGradient +=
	data_rep->tensor_product_gradient_basis_variables(x, t1_coeffs_l[set],
	  t2_coeffs_l[set], sm_mi_l[set], key_l[set], colloc_index,
	  subset_indices);
  }

  return approxGradient;
}


const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const UShort3DArray& sm_mi,
			 const UShort4DArray& colloc_key,
			 const RealVector2DArray& t1_coeffs,
			 const RealMatrix2DArray& t2_coeffs,
			 const SizetArray& dvv, unsigned short level,
			 const UShort2DArray& set_partition)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in HierarchInterpPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  size_t lev, set, set_start = 0, set_end, num_deriv_v = dvv.size();
  if (approxGradient.length() != num_deriv_v)
    approxGradient.sizeUninitialized(num_deriv_v);
  approxGradient = 0.;

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  SizetArray colloc_index; // empty -> 2DArrays allow default indexing
  bool partial = !set_partition.empty();
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = colloc_key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    if (partial)
      { set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
    else // sm_mi and key include all current index sets, whereas the t1/t2
         // coeffs may refect a partial state derived from a ref_key
      set_end = t1_coeffs_l.size();
    for (set=set_start; set<set_end; ++set)
      approxGradient +=
	data_rep->tensor_product_gradient_basis_variables(x, t1_coeffs_l[set],
	  t2_coeffs_l[set], sm_mi_l[set], key_l[set], colloc_index, dvv);
  }

  return approxGradient;
}


const RealVector& HierarchInterpPolyApproximation::
gradient_nonbasis_variables(const RealVector& x, const UShort3DArray& sm_mi,
			    const UShort4DArray& colloc_key,
			    const RealMatrix2DArray& t1_coeff_grads,
			    unsigned short level,
			    const UShort2DArray& set_partition)
{
  // Error check for required data
  size_t lev, set, set_start = 0, set_end, num_deriv_v;
  if (expansionCoeffGradFlag) {
    if (t1_coeff_grads.size() > level && t1_coeff_grads[level].size())
      num_deriv_v = t1_coeff_grads[level][0].numRows();
    else {
      PCerr << "Error: insufficient size in type1 expansion coefficient "
	    << "gradients in\n       HierarchInterpPolyApproximation::"
	    << "gradient_nonbasis_variables()" << std::endl;
      abort_handler(-1);
    }
  }
  else {
    PCerr << "Error: expansion coefficient gradients not defined in Hierarch"
	  << "InterpPolyApproximation::gradient_nonbasis_variables()"
	  << std::endl;
    abort_handler(-1);
  }

  if (approxGradient.length() != num_deriv_v)
    approxGradient.sizeUninitialized(num_deriv_v);
  approxGradient = 0.;

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  SizetArray colloc_index; // empty -> 2DArrays allow default indexing
  bool partial = !set_partition.empty();
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&            sm_mi_l = sm_mi[lev];
    const UShort3DArray&              key_l = colloc_key[lev];
    const RealMatrixArray& t1_coeff_grads_l = t1_coeff_grads[lev];
    if (partial)
      { set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
    else // sm_mi and key include all current index sets, whereas the t1/t2
         // coeffs may refect a partial state derived from a ref_key
      set_end = t1_coeff_grads_l.size();
    for (set=set_start; set<set_end; ++set)
      approxGradient +=
	data_rep->tensor_product_gradient_nonbasis_variables(x,
	  t1_coeff_grads_l[set], sm_mi_l[set], key_l[set], colloc_index);
  }

  return approxGradient;
}


const RealSymMatrix& HierarchInterpPolyApproximation::
hessian_basis_variables(const RealVector& x, const UShort3DArray& sm_mi,
			const UShort4DArray& colloc_key,
			const RealVector2DArray& t1_coeffs,
			unsigned short level,
			const UShort2DArray& set_partition)
{
  PCerr << "Error: HierarchInterpPolyApproximation::hessian_basis_variables() "
	<< "not yet implemented." << std::endl;
  abort_handler(-1);

  return approxHessian;
}


Real HierarchInterpPolyApproximation::mean()
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::mean()" << std::endl;
    abort_handler(-1);
  }

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedMean & 1))
    return numericalMoments[0];

  Real mean = expectation(expT1CoeffsIter->second, expT2CoeffsIter->second);
  if (std_mode)
    { numericalMoments[0] = mean; computedMean |= 1; }
  return mean;
}



Real HierarchInterpPolyApproximation::mean(const RealVector& x)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if (all_mode && (computedMean & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevMean))
    return numericalMoments[0];

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  Real mean = expectation(x, expT1CoeffsIter->second, expT2CoeffsIter->second);
  if (all_mode)
    { numericalMoments[0] = mean; computedMean |= 1; xPrevMean = x; }
  return mean;
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  In
    this case, the derivative of the expectation is the expectation of
    the derivative.  The mixed derivative case (some design variables
    are inserted and some are augmented) requires no special treatment. */
const RealVector& HierarchInterpPolyApproximation::mean_gradient()
{
  // Error check for required data
  if (!expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in Hierarch"
	  << "InterpPolyApproximation::mean_gradient()." << std::endl;
    abort_handler(-1);
  }

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedMean & 2))
    return meanGradient;

  meanGradient = expectation_gradient(expT1CoeffGradsIter->second);
  if (std_mode) computedMean |=  2; //   activate 2-bit
  else          computedMean &= ~2; // deactivate 2-bit: protect mixed usage
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
    are inserted (derivatives are obtained from expansionType1CoeffGrads). */
const RealVector& HierarchInterpPolyApproximation::
mean_gradient(const RealVector& x, const SizetArray& dvv)
{
  // if already computed, return previous result
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if ( all_mode && (computedMean & 2) &&
       data_rep->match_nonrandom_vars(x, xPrevMeanGrad) ) // && dvv == dvvPrev)
    return meanGradient;

  // ---------------------------------------------------------------------
  // For xi = ran vars, sa = augmented des vars, si = inserted design vars
  // Active variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(sa, si) L_i(xi)
  //   mu(sa, si)    = Sum_i r_i(sa, si) wt_prod_i
  //   dmu/ds        = Sum_i dr_i/ds wt_prod_i
  // All variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(si) L_i(xi, sa)
  //   mu(sa, si)    = Sum_i r_i(si) Lsa_i wt_prod_i
  //   dmu/dsa       = Sum_i r_i(si) dLsa_i/dsa wt_prod_i
  //   dmu/dsi       = Sum_i dr_i/dsi Lsa_i wt_prod_i
  // ---------------------------------------------------------------------
  size_t i, deriv_index, cntr = 0, num_deriv_vars = dvv.size();
  if (meanGradient.length() != num_deriv_vars)
    meanGradient.sizeUninitialized(num_deriv_vars);
  const RealVector2DArray& exp_t1_coeffs = expT1CoeffsIter->second;
  const RealMatrix2DArray& exp_t2_coeffs = expT2CoeffsIter->second;
  const RealMatrix2DArray& exp_t1_coeff_grads = expT1CoeffGradsIter->second;

  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    if (data_rep->randomVarsKey[deriv_index]) {
      // --------------------------------------------------------------------
      // derivative of All var expansion w.r.t. random var (design insertion)
      // --------------------------------------------------------------------
      if (!expansionCoeffGradFlag) { // required data check
	PCerr << "Error: expansion coefficient gradients not defined in "
	      << "HierarchInterpPolyApproximation::mean_gradient()."
	      << std::endl;
	abort_handler(-1);
      }
      if (data_rep->basisConfigOptions.useDerivs) {
	PCerr << "Error: combination of coefficient gradients and use_"
	      << "derivatives is not supported in HierarchInterpPoly"
	      << "Approximation::mean_gradient()." << std::endl;
	abort_handler(-1);
      }
      meanGradient[i] = expectation_gradient(x, exp_t1_coeff_grads, cntr);
      ++cntr;
    }
    else {
      // ---------------------------------------------------------------------
      // deriv of All var expansion w.r.t. nonrandom var (design augmentation)
      // ---------------------------------------------------------------------
      if (!expansionCoeffFlag) { // required data check
	PCerr << "Error: expansion coefficients not defined in HierarchInterp"
	      << "PolyApproximation::mean_gradient()." << std::endl;
	abort_handler(-1);
      }
      meanGradient[i]
	= expectation_gradient(x, exp_t1_coeffs, exp_t2_coeffs, deriv_index);
    }
  }
  if (all_mode) { computedMean |=  2; xPrevMeanGrad = x; }
  else            computedMean &= ~2; // deactivate 2-bit: protect mixed usage
  return meanGradient;
}


Real HierarchInterpPolyApproximation::
covariance(PolynomialApproximation* poly_approx_2)
{
  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool same = (this == hip_approx_2),
    std_mode = data_rep->nonRandomIndices.empty();

  // Error check for required data
  if ( !expansionCoeffFlag ||
       ( !same && !hip_approx_2->expansionCoeffFlag ) ) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  if (same && std_mode && (computedVariance & 1))
    return numericalMoments[1];

  RealVector2DArray cov_t1_coeffs; RealMatrix2DArray cov_t2_coeffs;
  Real mean_1 = mean(), mean_2 = (same) ? mean_1 : hip_approx_2->mean();
  central_product_interpolant(hip_approx_2, mean_1, mean_2,
			      cov_t1_coeffs, cov_t2_coeffs);

  // evaluate expectation of these t1/t2 coefficients
  Real covar = expectation(cov_t1_coeffs, cov_t2_coeffs);
  // Note: separation of reference and increment using cov_t{1,2}_coeffs
  // with {ref,incr}_key would provide an increment of a central moment
  // around an invariant center.  For hierarchical covariance, one must
  // also account for the change in mean as in delta_covariance().

  if (same && std_mode)
    { numericalMoments[1] = covar; computedVariance |= 1; }
  return covar;
}


Real HierarchInterpPolyApproximation::
covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool same = (this == hip_approx_2),
    all_mode = !data_rep->nonRandomIndices.empty();

  // Error check for required data
  if ( !expansionCoeffFlag ||
       ( !same && !hip_approx_2->expansionCoeffFlag ) ) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  if ( same && all_mode && (computedVariance & 1) &&
       data_rep->match_nonrandom_vars(x, xPrevVar) )
    return numericalMoments[1];

  RealVector2DArray cov_t1_coeffs; RealMatrix2DArray cov_t2_coeffs;
  Real mean_1 = mean(x), mean_2 = (same) ? mean_1 : hip_approx_2->mean(x);
  central_product_interpolant(hip_approx_2, mean_1, mean_2,
			      cov_t1_coeffs, cov_t2_coeffs);

  // evaluate expectation of these t1/t2 coefficients
  Real covar = expectation(x, cov_t1_coeffs, cov_t2_coeffs);

  if (same && all_mode)
    { numericalMoments[1] = covar; computedVariance |= 1; xPrevVar = x; }
  return covar;
}


Real HierarchInterpPolyApproximation::
combined_covariance(PolynomialApproximation* poly_approx_2)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;
  bool same = (this == hip_approx_2);

  const RealVector2DArray& comb_t1c_2 = hip_approx_2->combinedExpT1Coeffs;
  const RealMatrix2DArray& comb_t2c_2 = hip_approx_2->combinedExpT2Coeffs;
  const RealVector2DArray& comb_t1w   = data_rep->combinedT1WeightSets;
  const RealMatrix2DArray& comb_t2w   = data_rep->combinedT2WeightSets;
  Real mean_1 = expectation(combinedExpT1Coeffs, combinedExpT2Coeffs,
			    comb_t1w, comb_t2w),
       mean_2 = (same) ? mean_1 : expectation(comb_t1c_2, comb_t2c_2,
					      comb_t1w, comb_t2w);

  RealVector2DArray cov_t1_coeffs; RealMatrix2DArray cov_t2_coeffs;
  central_product_interpolant(data_rep->combinedVarSets,
			      data_rep->combinedSmolyakMultiIndex,
			      data_rep->combinedCollocKey, combinedExpT1Coeffs,
			      combinedExpT2Coeffs, comb_t1c_2, comb_t2c_2, same,
			      mean_1, mean_2, cov_t1_coeffs, cov_t2_coeffs);

  // evaluate expectation of these t1/t2 coefficients
  return expectation(cov_t1_coeffs, cov_t2_coeffs, comb_t1w, comb_t2w);
}


Real HierarchInterpPolyApproximation::
combined_covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;
  bool same = (this == hip_approx_2);

  const RealVector2DArray& comb_t1c_2 = hip_approx_2->combinedExpT1Coeffs;
  const RealMatrix2DArray& comb_t2c_2 = hip_approx_2->combinedExpT2Coeffs;
  const UShort3DArray&     sm_mi      = hsg_driver->smolyak_multi_index();
  const UShort4DArray&     colloc_key = hsg_driver->collocation_key();
  Real mean_1 =
    expectation(x, combinedExpT1Coeffs, combinedExpT2Coeffs, sm_mi, colloc_key),
       mean_2 = (same) ? mean_1 :
    expectation(x, comb_t1c_2, comb_t2c_2, sm_mi, colloc_key);

  RealVector2DArray cov_t1_coeffs; RealMatrix2DArray cov_t2_coeffs;
  central_product_interpolant(data_rep->combinedVarSets,
			      data_rep->combinedSmolyakMultiIndex,
			      data_rep->combinedCollocKey, combinedExpT1Coeffs,
			      combinedExpT2Coeffs, comb_t1c_2, comb_t2c_2, same,
			      mean_1, mean_2, cov_t1_coeffs, cov_t2_coeffs);

  // evaluate expectation of these t1/t2 coefficients
  return expectation(x, cov_t1_coeffs, cov_t2_coeffs, sm_mi, colloc_key);
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  The
    mixed derivative case (some design variables are inserted and some
    are augmented) requires no special treatment. */
const RealVector& HierarchInterpPolyApproximation::variance_gradient()
{
  // Error check for required data
  if (!expansionCoeffFlag ||
      !expansionCoeffGradFlag) {
    PCerr << "Error: insufficient expansion coefficient data in HierarchInterp"
	  << "PolyApproximation::variance_gradient()." << std::endl;
    abort_handler(-1);
  }

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedVariance & 2))
    return varianceGradient;

  Real mean_1 = mean(); const RealVector& mean1_grad = mean_gradient();
  RealMatrix2DArray cov_t1_coeff_grads;
  central_product_gradient_interpolant(this, mean_1, mean_1, mean1_grad,
				       mean1_grad, cov_t1_coeff_grads);
  varianceGradient = expectation_gradient(cov_t1_coeff_grads);
  if (std_mode) computedVariance |=  2;
  else          computedVariance &= ~2; // deactivate 2-bit: protect mixed usage
  return varianceGradient;
}


/** In this function, a subset of the expansion variables are random
    variables and any augmented design/state variables (i.e., not
    inserted as random variable distribution parameters) are included
    in the expansion.  This function must handle the mixed case, where
    some design/state variables are augmented (and are part of the
    expansion) and some are inserted (derivatives are obtained from
    expansionType1CoeffGrads). */
const RealVector& HierarchInterpPolyApproximation::
variance_gradient(const RealVector& x, const SizetArray& dvv)
{
  // if already computed, return previous result
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if ( all_mode && (computedVariance & 2) &&
       data_rep->match_nonrandom_vars(x, xPrevVarGrad) ) // && dvv == dvvPrev)
    return varianceGradient;

  // ---------------------------------------------------------------------
  // For xi = ran vars, sa = augmented des vars, si = inserted design vars
  // Active variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(sa, si) L_i(xi)
  //   mu(sa, si)    = Sum_i r_i(sa, si) wt_prod_i
  //   dmu/ds        = Sum_i dr_i/ds wt_prod_i
  // All variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(si) L_i(xi, sa)
  //   mu(sa, si)    = Sum_i r_i(si) Lsa_i wt_prod_i
  //   dmu/dsa       = Sum_i r_i(si) dLsa_i/dsa wt_prod_i
  //   dmu/dsi       = Sum_i dr_i/dsi Lsa_i wt_prod_i
  // ---------------------------------------------------------------------
  size_t i, deriv_index, cntr = 0, num_deriv_vars = dvv.size();
  bool insert = false, augment = false;
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    if (data_rep->randomVarsKey[deriv_index]) insert = true;
    else                                      augment = true;
  }

  RealVector2DArray cov_t1_coeffs;
  RealMatrix2DArray cov_t1_coeff_grads, cov_t2_coeffs;
  Real mean_1 = mean(x);
  if (insert) {
    const RealVector& mean1_grad = mean_gradient(x, dvv);
    central_product_gradient_interpolant(this, mean_1, mean_1, mean1_grad,
					 mean1_grad, cov_t1_coeff_grads);
  }
  if (augment)
    central_product_interpolant(this, mean_1, mean_1, cov_t1_coeffs,
				cov_t2_coeffs);

  if (varianceGradient.length() != num_deriv_vars)
    varianceGradient.sizeUninitialized(num_deriv_vars);
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    Real& grad_i = varianceGradient[i];
    if (data_rep->randomVarsKey[deriv_index]) {
      // --------------------------------------------------------------------
      // derivative of All var expansion w.r.t. random var (design insertion)
      // --------------------------------------------------------------------
      if (!expansionCoeffGradFlag) { // required data check
	PCerr << "Error: expansion coefficient gradients not defined in "
	      << "HierarchInterpPolyApproximation::variance_gradient()."
	      << std::endl;
	abort_handler(-1);
      }
      if (data_rep->basisConfigOptions.useDerivs) {
	PCerr << "Error: combination of coefficient gradients and use_"
	      << "derivatives is not supported in HierarchInterpPoly"
	      << "Approximation::variance_gradient()" << std::endl;
	abort_handler(-1);
      }
      varianceGradient[i] = expectation_gradient(x, cov_t1_coeff_grads, cntr);
      ++cntr;
    }
    else {
      // ---------------------------------------------------------------------
      // deriv of All var expansion w.r.t. nonrandom var (design augmentation)
      // ---------------------------------------------------------------------
      if (!expansionCoeffFlag) { // required data check
	PCerr << "Error: expansion coefficients not defined in Hierarch"
	      << "InterpPolyApproximation::variance_gradient()." << std::endl;
	abort_handler(-1);
      }
      varianceGradient[i] = expectation_gradient(x, cov_t1_coeffs,
						 cov_t2_coeffs, deriv_index);
    }
  }
  if (all_mode) { computedVariance |=  2; xPrevVarGrad = x; }
  else            computedVariance &= ~2;//deactivate 2-bit: protect mixed usage
  return varianceGradient;
}


Real HierarchInterpPolyApproximation::
reference_mean(const UShort2DArray& ref_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedRefMean & 1))
    return referenceMoments[0];

  Real ref_mean
    = expectation(expT1CoeffsIter->second, expT2CoeffsIter->second, ref_key);
  if (std_mode)
    { referenceMoments[0] = ref_mean; computedRefMean |= 1; }
  return ref_mean;
}


Real HierarchInterpPolyApproximation::
reference_mean(const RealVector& x, const UShort2DArray& ref_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if (all_mode && (computedRefMean & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevRefMean))
    return referenceMoments[0];

  Real ref_mean
    = expectation(x, expT1CoeffsIter->second, expT2CoeffsIter->second, ref_key);
  if (all_mode)
    { referenceMoments[0] = ref_mean; computedRefMean |= 1; xPrevRefMean = x; }
  return ref_mean;
}


Real HierarchInterpPolyApproximation::
reference_variance(const UShort2DArray& ref_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedRefVariance & 1))
    return referenceMoments[1];

  Real ref_mean = reference_mean(ref_key);
  RealVector2DArray cov_t1_coeffs; RealMatrix2DArray cov_t2_coeffs;
  central_product_interpolant(this, ref_mean, ref_mean, cov_t1_coeffs,
			      cov_t2_coeffs, ref_key);
  Real ref_var = expectation(cov_t1_coeffs, cov_t2_coeffs, ref_key);
  if (std_mode)
    { referenceMoments[1] = ref_var; computedRefVariance |= 1; }
  return ref_var;
}


Real HierarchInterpPolyApproximation::
reference_variance(const RealVector& x, const UShort2DArray& ref_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if (all_mode && (computedRefVariance & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevRefVar))
    return referenceMoments[1];

  Real ref_mean = reference_mean(x, ref_key);
  RealVector2DArray cov_t1_coeffs; RealMatrix2DArray cov_t2_coeffs;
  central_product_interpolant(this, ref_mean, ref_mean, cov_t1_coeffs,
			      cov_t2_coeffs, ref_key);
  Real ref_var = expectation(x, cov_t1_coeffs, cov_t2_coeffs, ref_key);
  if (all_mode) {
    referenceMoments[1] = ref_var;
    computedRefVariance |= 1; xPrevRefVar = x;
  }
  return ref_var;
}


Real HierarchInterpPolyApproximation::delta_mean(const UShort2DArray& incr_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedDeltaMean & 1))
    return deltaMoments[0];

  Real delta_mean
    = expectation(expT1CoeffsIter->second, expT2CoeffsIter->second, incr_key);
  if (std_mode)
    { deltaMoments[0] = delta_mean; computedDeltaMean |= 1; }
  return delta_mean;
}


Real HierarchInterpPolyApproximation::
delta_mean(const RealVector& x, const UShort2DArray& incr_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if (all_mode && (computedDeltaMean & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevDeltaMean))
    return deltaMoments[0];

  Real delta_mean = expectation(x, expT1CoeffsIter->second,
				expT2CoeffsIter->second, incr_key);
  if (all_mode) {
    deltaMoments[0] = delta_mean;
    computedDeltaMean |= 1; xPrevDeltaMean = x;
  }
  return delta_mean;
}


Real HierarchInterpPolyApproximation::
delta_variance(const UShort2DArray& ref_key, const UShort2DArray& incr_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedDeltaVariance & 1))
    return deltaMoments[1];

  // delta_variance() can leverage surrData for computing product interpolant:
  RealVector2DArray r1r2_t1_coeffs; RealMatrix2DArray r1r2_t2_coeffs;
  product_interpolant(this, r1r2_t1_coeffs, r1r2_t2_coeffs);//, incr_key);

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  Real delta_var =
    delta_covariance(expT1CoeffsIter->second, expT2CoeffsIter->second,
		     expT1CoeffsIter->second, expT2CoeffsIter->second, true,
		     r1r2_t1_coeffs, r1r2_t2_coeffs,
		     hsg_driver->type1_hierarchical_weight_sets(),
		     hsg_driver->type2_hierarchical_weight_sets(),
		     ref_key, incr_key);
  if (std_mode)
    { deltaMoments[1] = delta_var; computedDeltaVariance |= 1; }
  return delta_var;
}


Real HierarchInterpPolyApproximation::
delta_variance(const RealVector& x, const UShort2DArray& ref_key,
	       const UShort2DArray& incr_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if (all_mode && (computedDeltaVariance & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevDeltaVar))
    return deltaMoments[1];

  // delta_variance() can leverage surrData for computing product interpolant:
  RealVector2DArray r1r2_t1_coeffs; RealMatrix2DArray r1r2_t2_coeffs;
  product_interpolant(this, r1r2_t1_coeffs, r1r2_t2_coeffs);//, incr_key);

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  Real delta_var =
    delta_covariance(x, expT1CoeffsIter->second, expT2CoeffsIter->second,
		     expT1CoeffsIter->second, expT2CoeffsIter->second, true,
		     r1r2_t1_coeffs, r1r2_t2_coeffs,
		     hsg_driver->smolyak_multi_index(),
		     hsg_driver->collocation_key(), ref_key, incr_key);
  if (all_mode) {
    deltaMoments[1] = delta_var;
    computedDeltaVariance |= 1; xPrevDeltaVar = x;
  }
  return delta_var;
}


Real HierarchInterpPolyApproximation::
delta_std_deviation(const UShort2DArray& ref_key, const UShort2DArray& incr_key)
{
  // delta-sigma = sqrt( var0 + delta-var ) - sigma0
  //             = [ sqrt(1 + delta_var / var0) - 1 ] * sigma0
  //             = sqrt1pm1(delta_var / var0) * sigma0
  // where sqrt1pm1(x) = expm1[ log1p(x) / 2 ]

  Real delta_var = delta_variance(ref_key, incr_key),
       var0      = reference_variance(ref_key),
       sigma0    = std::sqrt(var0);

  return (delta_var < var0) ?
    bmth::sqrt1pm1(delta_var / var0) * sigma0 :            // preserve precision
    std::sqrt(var0 + delta_var) - sigma0; // precision OK; prevent div by var0=0
}


Real HierarchInterpPolyApproximation::
delta_std_deviation(const RealVector& x, const UShort2DArray& ref_key,
		    const UShort2DArray& incr_key)
{
  // delta-sigma = sqrt( var0 + delta-var ) - sigma0
  //             = [ sqrt(1 + delta_var / var0) - 1 ] * sigma0
  //             = sqrt1pm1(delta_var / var0) * sigma0
  // where sqrt1pm1(x) = expm1[ log1p(x) / 2 ]

  Real delta_var = delta_variance(x, ref_key, incr_key),
       var0      = reference_variance(x, ref_key),
       sigma0    = std::sqrt(var0);

  return (delta_var < var0) ?
    bmth::sqrt1pm1(delta_var / var0) * sigma0 :            // preserve precision
    std::sqrt(var0 + delta_var) - sigma0; // precision OK; prevent div by var0=0
}


Real HierarchInterpPolyApproximation::
delta_beta(bool cdf_flag, Real z_bar, const UShort2DArray& ref_key,
	   const UShort2DArray& incr_key)
{
  Real mu0      = reference_mean(ref_key),
    delta_mu    = delta_mean(incr_key),
    var0        = reference_variance(ref_key), 
    delta_sigma = delta_std_deviation(ref_key, incr_key);
  return delta_beta_map(mu0, delta_mu, var0, delta_sigma, cdf_flag, z_bar);
}


Real HierarchInterpPolyApproximation::
delta_beta(const RealVector& x, bool cdf_flag, Real z_bar,
	   const UShort2DArray& ref_key, const UShort2DArray& incr_key)
{
  Real mu0      = reference_mean(x, ref_key),
    delta_mu    = delta_mean(x, incr_key),
    var0        = reference_variance(x, ref_key), 
    delta_sigma = delta_std_deviation(x, ref_key, incr_key);
  return delta_beta_map(mu0, delta_mu, var0, delta_sigma, cdf_flag, z_bar);
}


Real HierarchInterpPolyApproximation::
delta_beta_map(Real mu0, Real delta_mu, Real var0, Real delta_sigma,
	       bool cdf_flag, Real z_bar)
{
  //  CDF delta-beta = (mu1 - z-bar)/sigma1 - (mu0 - z-bar)/sigma0
  //    = (mu1 sigma0 - z-bar sigma0 - mu0 sigma1 + z-bar sigma1)/sigma1/sigma0
  //    = (delta-mu sigma0 - mu0 delta-sigma + z-bar delta-sigma)/sigma1/sigma0
  //    = (delta-mu - delta-sigma beta0)/sigma1
  // CCDF delta-beta = (z-bar - mu1)/sigma1 - (z-bar - mu0)/sigma0
  //    = (z-bar sigma0 - mu1 sigma0 - z-bar sigma1 + mu0 sigma1)/sigma1/sigma0
  //    = (mu0 delta-sigma - delta-mu sigma0 - z-bar delta-sigma)/sigma1/sigma0
  //    = -delta-mu/sigma1 - delta_sigma (z-bar - mu0) / sigma0 / sigma1
  //    = (-delta-mu - delta-sigma beta0)/sigma1

  // Error traps are needed for zero variance: a single point ref grid
  // (level=0 sparse or m=1 tensor) has zero variance.  Unchanged response
  // values along an index set could then cause sigma1 also = 0.
  Real beta0, sigma0 = std::sqrt(var0), sigma1 = sigma0 + delta_sigma;
  if (cdf_flag) {
    if (sigma0 > SMALL_NUMBER && sigma1 > SMALL_NUMBER) {
      beta0 = (mu0 - z_bar) / sigma0;
      return ( delta_mu - delta_sigma * beta0) / sigma1;
    }
    else if (sigma1 > SMALL_NUMBER)// neglect beta0 term (zero init reliability)
      return delta_mu / sigma1; // or delta = beta1 = (mu1 - z_bar) / sigma1 ?
    else if (sigma0 > SMALL_NUMBER) // assume beta1 = 0 -> delta = -beta0
      return (z_bar - mu0) / sigma0;
    else                      // assume beta0 = beta1 = 0
      return 0;
  }
  else {
    if (sigma0 > SMALL_NUMBER && sigma1 > SMALL_NUMBER) {
      beta0 = (z_bar - mu0) / sigma0;
      return (-delta_mu - delta_sigma * beta0) / sigma1;
    }
    else if (sigma1 > SMALL_NUMBER)// neglect beta0 term (zero init reliability)
      return -delta_mu / sigma1;
    else if (sigma0 > SMALL_NUMBER) // assume beta1 = 0 -> delta = -beta0
      return (mu0 - z_bar) / sigma0;
    else                      // assume beta0 = beta1 = 0
      return 0;
  }
}


Real HierarchInterpPolyApproximation::
delta_z(bool cdf_flag, Real beta_bar, const UShort2DArray& ref_key,
	const UShort2DArray& incr_key)
{
  //  CDF delta-z = (mu1 - sigma1 beta-bar) - (mu0 - sigma0 beta-bar)
  //              = delta-mu - delta-sigma * beta-bar
  // CCDF delta-z = (mu1 + sigma1 beta-bar) - (mu0 + sigma0 beta-bar)
  //              = delta-mu + delta-sigma * beta-bar

  Real delta_mu = delta_mean(incr_key),
    delta_sigma = delta_std_deviation(ref_key, incr_key);
  return (cdf_flag) ? delta_mu - delta_sigma * beta_bar :
                      delta_mu + delta_sigma * beta_bar;
}


Real HierarchInterpPolyApproximation::
delta_z(const RealVector& x, bool cdf_flag, Real beta_bar,
	const UShort2DArray& ref_key, const UShort2DArray& incr_key)
{
  //  CDF delta-z = (mu1 - sigma1 beta-bar) - (mu0 - sigma0 beta-bar)
  //              = delta-mu - delta-sigma * beta-bar
  // CCDF delta-z = (mu1 + sigma1 beta-bar) - (mu0 + sigma0 beta-bar)
  //              = delta-mu + delta-sigma * beta-bar

  Real delta_mu = delta_mean(x, incr_key),
    delta_sigma = delta_std_deviation(x, ref_key, incr_key);
  return (cdf_flag) ? delta_mu - delta_sigma * beta_bar :
                      delta_mu + delta_sigma * beta_bar;
}


Real HierarchInterpPolyApproximation::
delta_covariance(PolynomialApproximation* poly_approx_2)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  UShort2DArray ref_key, incr_key;
  hsg_driver->partition_keys(ref_key, incr_key);

  // Supports multiple grid increments in discerning nominal from delta based
  // on isotropic/anisotropic/generalized index set increments.  In current
  // use, 2D keys with set ranges are sufficient: level -> {start,end} set.
  // In the future, may need 3D keys for level/set/point.
  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;
  bool same = (this == hip_approx_2);

  // Error check for required data
  if ( !expansionCoeffFlag ||
       ( !same && !hip_approx_2->expansionCoeffFlag )) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::delta_covariance()" << std::endl;
    abort_handler(-1);
  }

  // delta_covariance() can leverage surrData for computing product interpolant:
  RealVector2DArray r1r2_t1_coeffs; RealMatrix2DArray r1r2_t2_coeffs;
  product_interpolant(hip_approx_2,r1r2_t1_coeffs,r1r2_t2_coeffs);//incr_key);
  // only need expectation of R1R2 on incr_key, but product interpolant is
  // nonlinear and we form hierarchical evaluations based on deltaR1R2

  Real delta_covar =
    delta_covariance(expT1CoeffsIter->second, expT2CoeffsIter->second,
		     hip_approx_2->expT1CoeffsIter->second,
		     hip_approx_2->expT2CoeffsIter->second, same,
		     r1r2_t1_coeffs, r1r2_t2_coeffs,
		     hsg_driver->type1_hierarchical_weight_sets(),
		     hsg_driver->type2_hierarchical_weight_sets(),
		     ref_key, incr_key);
  if (same && data_rep->nonRandomIndices.empty()) // std mode
    { deltaMoments[1] = delta_covar; computedDeltaVariance |= 1; }
  return delta_covar;
}


Real HierarchInterpPolyApproximation::
delta_covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  UShort2DArray ref_key, incr_key;
  hsg_driver->partition_keys(ref_key, incr_key);

  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;
  bool same = (this == hip_approx_2);

  // Error check for required data
  if ( !expansionCoeffFlag ||
       ( !same && !hip_approx_2->expansionCoeffFlag )) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::delta_covariance()" << std::endl;
    abort_handler(-1);
  }

  // delta_covariance() can leverage surrData for computing product interpolant:
  RealVector2DArray r1r2_t1_coeffs; RealMatrix2DArray r1r2_t2_coeffs;
  product_interpolant(hip_approx_2, r1r2_t1_coeffs, r1r2_t2_coeffs);//incr_key);
  // only need expectation of R1R2 on incr_key, but product interpolant is
  // nonlinear and we form hierarchical evaluations based on deltaR1R2

  Real delta_covar =
    delta_covariance(x, expT1CoeffsIter->second, expT2CoeffsIter->second,
		     hip_approx_2->expT1CoeffsIter->second,
		     hip_approx_2->expT2CoeffsIter->second, same,
		     r1r2_t1_coeffs, r1r2_t2_coeffs,
		     hsg_driver->smolyak_multi_index(),
		     hsg_driver->collocation_key(), ref_key, incr_key);
  if (same && !data_rep->nonRandomIndices.empty()) { // all vars mode
    deltaMoments[1] = delta_covar;
    computedDeltaVariance |= 1; xPrevDeltaVar = x;
  }
  return delta_covar;
}


Real HierarchInterpPolyApproximation::
delta_combined_covariance(PolynomialApproximation* poly_approx_2)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  UShort2DArray ref_key, incr_key;
  data_rep->hsg_driver()->partition_keys(ref_key, incr_key); // *** TO DO ***

  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;
  bool same = (this == hip_approx_2);

  // delta_combined_covariance() cannot leverage surrData for r1r2 and must
  // form an interpolant for r1 * r2 using the combined coefficients:
  const RealVector2DArray& comb_t1c_2 = hip_approx_2->combinedExpT1Coeffs;
  const RealMatrix2DArray& comb_t2c_2 = hip_approx_2->combinedExpT2Coeffs;
  RealVector2DArray r1r2_t1_coeffs; RealMatrix2DArray r1r2_t2_coeffs;
  product_interpolant(data_rep->combinedVarSets,
		      data_rep->combinedSmolyakMultiIndex,
		      data_rep->combinedCollocKey, combinedExpT1Coeffs,
		      combinedExpT2Coeffs, comb_t1c_2, comb_t2c_2, same,
		      r1r2_t1_coeffs, r1r2_t2_coeffs);//, incr_key);
  // only need expectation of R1R2 on incr_key, but product interpolant is
  // nonlinear and we form hierarchical evaluations based on deltaR1R2

  return delta_covariance(combinedExpT1Coeffs, combinedExpT2Coeffs,
			  comb_t1c_2, comb_t2c_2, same, r1r2_t1_coeffs,
			  r1r2_t2_coeffs, data_rep->combinedT1WeightSets,
			  data_rep->combinedT2WeightSets, ref_key, incr_key);
}


Real HierarchInterpPolyApproximation::
delta_combined_covariance(const RealVector& x,
			  PolynomialApproximation* poly_approx_2)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  UShort2DArray ref_key, incr_key;
  data_rep->hsg_driver()->partition_keys(ref_key, incr_key); // *** TO DO ***

  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;
  bool same = (this == hip_approx_2);

  // delta_combined_covariance() cannot leverage surrData for r1r2 and must
  // form interpolants for r1 and r2 using the combined coefficients:
  const RealVector2DArray& comb_t1c_2 = hip_approx_2->combinedExpT1Coeffs;
  const RealMatrix2DArray& comb_t2c_2 = hip_approx_2->combinedExpT2Coeffs;
  const UShort3DArray& comb_sm_mi = data_rep->combinedSmolyakMultiIndex;
  const UShort4DArray& comb_key   = data_rep->combinedCollocKey;
  RealVector2DArray r1r2_t1_coeffs; RealMatrix2DArray r1r2_t2_coeffs;
  product_interpolant(data_rep->combinedVarSets, comb_sm_mi, comb_key,
		      combinedExpT1Coeffs, combinedExpT2Coeffs, comb_t1c_2,
		      comb_t2c_2, same, r1r2_t1_coeffs, r1r2_t2_coeffs);
                      //, incr_key);
  // only need expectation of R1R2 on incr_key, but product interpolant is
  // nonlinear and we form hierarchical evaluations based on deltaR1R2

  return delta_covariance(x, combinedExpT1Coeffs, combinedExpT2Coeffs,
			  comb_t1c_2, comb_t2c_2, same, r1r2_t1_coeffs,
			  r1r2_t2_coeffs, comb_sm_mi, comb_key,
			  ref_key, incr_key);
}


Real HierarchInterpPolyApproximation::
delta_covariance(const RealVector2DArray& r1_t1_coeffs,
		 const RealMatrix2DArray& r1_t2_coeffs,
		 const RealVector2DArray& r2_t1_coeffs,
		 const RealMatrix2DArray& r2_t2_coeffs, bool same,
		 const RealVector2DArray& r1r2_t1_coeffs,
		 const RealMatrix2DArray& r1r2_t2_coeffs,
		 const RealVector2DArray& t1_wts,
		 const RealMatrix2DArray& t2_wts, const UShort2DArray& ref_key,
		 const UShort2DArray& incr_key)
{
  // Compute surplus for r1, r2, and r1r2 and retrieve reference mean values
  Real  ref_mean_r1 =
    expectation(r1_t1_coeffs, r1_t2_coeffs, t1_wts, t2_wts, ref_key),
      delta_mean_r1 =
    expectation(r1_t1_coeffs, r1_t2_coeffs, t1_wts, t2_wts, incr_key),
        ref_mean_r2 = (same) ? ref_mean_r1 :
    expectation(r2_t1_coeffs, r2_t2_coeffs, t1_wts, t2_wts, ref_key),
      delta_mean_r2 = (same) ? delta_mean_r1 :
    expectation(r2_t1_coeffs, r2_t2_coeffs, t1_wts, t2_wts, incr_key),
    delta_mean_r1r2 =
    expectation(r1r2_t1_coeffs, r1r2_t2_coeffs, t1_wts, t2_wts, incr_key);

  // Hierarchical increment to covariance:
  // \Delta\Sigma_ij = \Sigma^1_ij - \Sigma^0_ij
  //   = ( E[Ri Rj]^1 - E[Ri]^1 E[Rj]^1 ) - ( E[Ri Rj]^0 - E[Ri]^0 E[Rj]^0 )
  //   = E[Ri Rj]^0 + \DeltaE[Ri Rj]
  //     - (E[Ri]^0 + \DeltaE[Ri]) (E[Rj]^0 + \DeltaE[Rj])
  //     - E[Ri Rj]^0 + E[Ri]^0 E[Rj]^0
  //   = \DeltaE[Ri Rj] - \DeltaE[Ri] E[Rj]^0 - E[Ri]^0 \DeltaE[Rj]
  //     - \DeltaE[Ri] \DeltaE[Rj]
  return delta_mean_r1r2 - ref_mean_r1 * delta_mean_r2 -
    ref_mean_r2 * delta_mean_r1 - delta_mean_r1 * delta_mean_r2;
}


Real HierarchInterpPolyApproximation::
delta_covariance(const std::map<UShortArray, RealVector2DArray>& r1_t1c_map,
		 const std::map<UShortArray, RealMatrix2DArray>& r1_t2c_map,
		 const std::map<UShortArray, RealVector2DArray>& r2_t1c_map,
		 const std::map<UShortArray, RealMatrix2DArray>& r2_t2c_map,
		 bool same,
		 const std::map<UShortArray, RealVector2DArray>& r1r2_t1c_map,
		 const std::map<UShortArray, RealMatrix2DArray>& r1r2_t2c_map,
		 const std::map<UShortArray, RealVector2DArray>& t1_wts_map,
		 const std::map<UShortArray, RealMatrix2DArray>& t2_wts_map,
		 const UShortArray& active_key,
		 const std::map<UShortArray, UShort2DArray>& ref_key_map,
		 const std::map<UShortArray, UShort2DArray>& incr_key_map)
{
  // Compute surplus for r1, r2, and r1r2 and retrieve reference mean values
  std::map<UShortArray, RealVector2DArray>::const_iterator r1t1_cit
    = r1_t1c_map.find(active_key), t1w_cit = t1_wts_map.find(active_key);
  std::map<UShortArray, RealMatrix2DArray>::const_iterator r1t2_cit
    = r1_t2c_map.find(active_key), t2w_cit = t2_wts_map.find(active_key);
  std::map<UShortArray, UShort2DArray>::const_iterator incr_cit
    = incr_key_map.find(active_key);
  if (r1t1_cit == r1_t1c_map.end() || t1w_cit == t1_wts_map.end() ||
      r1t2_cit == r1_t2c_map.end() || t2w_cit == t2_wts_map.end() ||
      incr_cit == incr_key_map.end()) {
    PCerr << "Error: failure in key lookup in HierarchInterpPolyApproximation"
	  << "::delta_covariance()" << std::endl;
    abort_handler(-1);
  }
  const RealVector2DArray& active_t1_wts =  t1w_cit->second;
  const RealMatrix2DArray& active_t2_wts =  t2w_cit->second;
  const UShort2DArray&   active_incr_key = incr_cit->second;

  Real  ref_mean_r1 =
    expectation(r1_t1c_map, r1_t2c_map, t1_wts_map, t2_wts_map, ref_key_map),
      delta_mean_r1 =
  //expectation(r1_t1c_map, r1_t2c_map, t1_wts_map, t2_wts_map, incr_key_map),
    expectation(r1t1_cit->second, r1t2_cit->second, active_t1_wts,
		active_t2_wts, active_incr_key);
  Real ref_mean_r2, delta_mean_r2;
  if (same)
    { ref_mean_r2 = ref_mean_r1; delta_mean_r2 = delta_mean_r1; }
  else {
    std::map<UShortArray, RealVector2DArray>::const_iterator r2t1_cit
      = r2_t1c_map.find(active_key);
    std::map<UShortArray, RealMatrix2DArray>::const_iterator r2t2_cit
      = r2_t2c_map.find(active_key);
    ref_mean_r2 =
      expectation(r2_t1c_map, r2_t2c_map, t1_wts_map, t2_wts_map, ref_key_map);
    delta_mean_r2 =
    //expectation(r2_t1c_map, r2_t2c_map, t1_wts_map, t2_wts_map, incr_key_map),
      expectation(r2t1_cit->second, r2t2_cit->second, active_t1_wts,
		  active_t2_wts, active_incr_key);
  }
  Real delta_mean_r1r2 =
    expectation(r1r2_t1c_map, r1r2_t2c_map, t1_wts_map, t2_wts_map,
		incr_key_map);

  // Hierarchical increment to covariance:
  // \Delta\Sigma_ij = \Sigma^1_ij - \Sigma^0_ij
  //   = ( E[Ri Rj]^1 - E[Ri]^1 E[Rj]^1 ) - ( E[Ri Rj]^0 - E[Ri]^0 E[Rj]^0 )
  //   = E[Ri Rj]^0 + \DeltaE[Ri Rj]
  //     - (E[Ri]^0 + \DeltaE[Ri]) (E[Rj]^0 + \DeltaE[Rj])
  //     - E[Ri Rj]^0 + E[Ri]^0 E[Rj]^0
  //   = \DeltaE[Ri Rj] - \DeltaE[Ri] E[Rj]^0 - E[Ri]^0 \DeltaE[Rj]
  //     - \DeltaE[Ri] \DeltaE[Rj]
  return delta_mean_r1r2 - ref_mean_r1 * delta_mean_r2 -
    ref_mean_r2 * delta_mean_r1 - delta_mean_r1 * delta_mean_r2;
}


Real HierarchInterpPolyApproximation::
delta_covariance(const RealVector& x, const RealVector2DArray& r1_t1_coeffs,
		 const RealMatrix2DArray& r1_t2_coeffs,
		 const RealVector2DArray& r2_t1_coeffs,
		 const RealMatrix2DArray& r2_t2_coeffs, bool same,
		 const RealVector2DArray& r1r2_t1_coeffs,
		 const RealMatrix2DArray& r1r2_t2_coeffs,
		 const UShort3DArray& sm_mi,   const UShort4DArray& colloc_key,
		 const UShort2DArray& ref_key, const UShort2DArray& incr_key)
{
  // Compute surplus for r1, r2, and r1r2 and retrieve reference mean values
  Real  ref_mean_r1 =
    expectation(x, r1_t1_coeffs, r1_t2_coeffs, sm_mi, colloc_key, ref_key),
      delta_mean_r1 =
    expectation(x, r1_t1_coeffs, r1_t2_coeffs, sm_mi, colloc_key, incr_key),
        ref_mean_r2 = (same) ? ref_mean_r1 :
    expectation(x, r2_t1_coeffs, r2_t2_coeffs, sm_mi, colloc_key, ref_key),
      delta_mean_r2 = (same) ? delta_mean_r1 :
    expectation(x, r2_t1_coeffs, r2_t2_coeffs, sm_mi, colloc_key, incr_key),
    delta_mean_r1r2 =
    expectation(x, r1r2_t1_coeffs, r1r2_t2_coeffs, sm_mi, colloc_key, incr_key);

  // same expression as standard expansion mode case above
  return delta_mean_r1r2 - ref_mean_r1 * delta_mean_r2 -
    ref_mean_r2 * delta_mean_r1 - delta_mean_r1 * delta_mean_r2;
}


Real HierarchInterpPolyApproximation::
expectation(const RealVector2DArray& t1_coeffs,
	    const RealMatrix2DArray& t2_coeffs, const RealVector2DArray& t1_wts,
	    const RealMatrix2DArray& t2_wts, const UShort2DArray& set_partition)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  Real integral = 0.;
  size_t lev, set, pt, num_levels = t1_coeffs.size(), set_start = 0, set_end,
    num_tp_pts;
  bool partial = !set_partition.empty();
  if (data_rep->basisConfigOptions.useDerivs) {
    size_t v, num_v = sharedDataRep->numVars;
    for (lev=0; lev<num_levels; ++lev) {
      const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
      if (partial)
	{ set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
      else
	set_end = t1_coeffs_l.size();
      for (set=set_start; set<set_end; ++set) {
	const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
	const RealMatrix& t2_coeffs_ls = t2_coeffs[lev][set];
	const RealVector&    t1_wts_ls = t1_wts[lev][set];
	const RealMatrix&    t2_wts_ls = t2_wts[lev][set];
	num_tp_pts = t1_coeffs_ls.length();
	for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	  integral += t1_coeffs_ls[pt] * t1_wts_ls[pt];
	  const Real* t2_coeffs_lsp = t2_coeffs_ls[pt];
	  const Real* t2_wts_lsp    = t2_wts_ls[pt];
	  for (v=0; v<num_v; ++v)
	    integral += t2_coeffs_lsp[v] * t2_wts_lsp[v];
	}
      }
    }
  }
  else {
    for (lev=0; lev<num_levels; ++lev) {
      const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
      if (partial)
	{ set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
      else
	set_end = t1_coeffs_l.size();
      for (set=set_start; set<set_end; ++set) {
	const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
	const RealVector&    t1_wts_ls = t1_wts[lev][set];
	num_tp_pts = t1_coeffs_ls.length();
	for (pt=0; pt<num_tp_pts; ++pt) // omitted if empty surplus vector
	  integral += t1_coeffs_ls[pt] * t1_wts_ls[pt];
      }
    }
  }
  return integral;
}


/*
Real HierarchInterpPolyApproximation::
expectation(const RealVector2DArray& t1_coeffs,
	    const RealMatrix2DArray& t2_coeffs,
	    const UShort3DArray& pt_partition)
{
  Real integral = 0.;
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const RealVector2DArray& t1_wts
    = hsg_driver->type1_hierarchical_weight_sets();
  size_t lev, set, pt, num_levels = t1_coeffs.size(), num_sets,
    tp_pt_start = 0, tp_pt_end;
  bool partial = !pt_partition.empty();
  switch (data_rep->basisConfigOptions.useDerivs) {
  case false:
    for (lev=0; lev<num_levels; ++lev) {
      const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
      num_sets = t1_coeffs_l.size();
      for (set=0; set<num_sets; ++set) {
	const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
	const RealVector&    t1_wts_ls = t1_wts[lev][set];
	if (partial) {
	  tp_pt_start = pt_partition[lev][set][0];
	  tp_pt_end   = pt_partition[lev][set][1];
	}
	else
	  tp_pt_end   = t1_coeffs_ls.length();
	for (pt=tp_pt_start; pt<tp_pt_end; ++pt)
	  integral += t1_coeffs_ls[pt] * t1_wts_ls[pt];
      }
    }
    break;
  case true: {
    const RealMatrix2DArray& t2_wts
      = hsg_driver->type2_hierarchical_weight_sets();
    size_t v, num_v = sharedDataRep->numVars;
    for (lev=0; lev<num_levels; ++lev) {
      const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
      num_sets = t1_coeffs_l.size();
      for (set=0; set<num_sets; ++set) {
	const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
	const RealMatrix& t2_coeffs_ls = t2_coeffs[lev][set];
	const RealVector&    t1_wts_ls = t1_wts[lev][set];
	const RealMatrix&    t2_wts_ls = t2_wts[lev][set];
	if (partial) {
	  tp_pt_start = pt_partition[lev][set][0];
	  tp_pt_end   = pt_partition[lev][set][1];
	}
	else
	  tp_pt_end   = t1_coeffs_ls.length();
	for (pt=tp_pt_start; pt<tp_pt_end; ++pt) {
	  integral += t1_coeffs_ls[pt] * t1_wts_ls[pt];
	  const Real* t2_coeffs_lsp = t2_coeffs_ls[pt];
	  const Real* t2_wts_lsp    = t2_wts_ls[pt];
	  for (v=0; v<num_v; ++v)
	    integral += t2_coeffs_lsp[v] * t2_wts_lsp[v];
	}
      }
    }
    break;
  }
  }

  return integral;
}
*/


Real HierarchInterpPolyApproximation::
expectation(const RealVector& x, const RealVector2DArray& t1_coeffs,
	    const RealMatrix2DArray& t2_coeffs, const UShort3DArray& sm_mi,
	    const UShort4DArray& colloc_key, const UShort2DArray& set_partition)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  Real integral = 0.;
  size_t lev, set, pt, num_levels = t1_coeffs.size(), set_start = 0, set_end,
    num_tp_pts;
  bool partial = !set_partition.empty();
  const SizetList&  r_ind = data_rep->randomIndices;
  const SizetList& nr_ind = data_rep->nonRandomIndices;
  if (data_rep->basisConfigOptions.useDerivs) {
    size_t v, num_v = sharedDataRep->numVars;
    for (lev=0; lev<num_levels; ++lev) {
      const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
      if (partial)
	{ set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
      else
	set_end = t1_coeffs_l.size();
      for (set=set_start; set<set_end; ++set) {
	const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
	const RealMatrix& t2_coeffs_ls = t2_coeffs[lev][set];
	const UShortArray&    sm_mi_ls = sm_mi[lev][set];
	num_tp_pts = t1_coeffs_ls.length();
	for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	  const UShortArray& key_lsp = colloc_key[lev][set][pt];
	  integral += t1_coeffs_ls[pt]
	    * data_rep->type1_interpolant_value(x, key_lsp, sm_mi_ls, nr_ind)
	    * data_rep->type1_weight(key_lsp, sm_mi_ls, r_ind);
	  const Real* t2_coeffs_lsp = t2_coeffs_ls[pt];
	  for (v=0; v<num_v; ++v)
	    integral += t2_coeffs_lsp[v]
	      * data_rep->type2_interpolant_value(x, v, key_lsp,sm_mi_ls,nr_ind)
	      * data_rep->type2_weight(v, key_lsp, sm_mi_ls, r_ind);
	}
      }
    }
  }
  else {
    for (lev=0; lev<num_levels; ++lev) {
      const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
      if (partial)
	{ set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
      else
	set_end = t1_coeffs_l.size();
      for (set=set_start; set<set_end; ++set) {
	const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
	const UShortArray&    sm_mi_ls = sm_mi[lev][set];
	num_tp_pts = t1_coeffs_ls.length();
	for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	  const UShortArray& key_lsp = colloc_key[lev][set][pt];
	  integral += t1_coeffs_ls[pt]
	    * data_rep->type1_interpolant_value(x, key_lsp, sm_mi_ls, nr_ind)
	    * data_rep->type1_weight(key_lsp, sm_mi_ls, r_ind);
	}
      }
    }
  }
  return integral;
}


Real HierarchInterpPolyApproximation::
expectation(const std::map<UShortArray, RealVector2DArray>& t1_coeffs_map,
	    const std::map<UShortArray, RealMatrix2DArray>& t2_coeffs_map,
	    const std::map<UShortArray, RealVector2DArray>& t1_wts_map,
	    const std::map<UShortArray, RealMatrix2DArray>& t2_wts_map,
	    const std::map<UShortArray, UShort2DArray>& set_partition_map)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  Real integral = 0.;
  size_t lev, set, pt, num_lev, start, end, num_tp_pts, v,
    num_v = sharedDataRep->numVars;
  std::map<UShortArray, RealVector2DArray>::const_iterator t1c_cit, t1w_cit;
  std::map<UShortArray, RealMatrix2DArray>::const_iterator t2c_cit, t2w_cit;
  std::map<UShortArray, UShort2DArray>::const_iterator p_cit;
  bool partial, use_derivs = data_rep->basisConfigOptions.useDerivs;
  for (t1c_cit = t1_coeffs_map.begin(),      t2c_cit  = t2_coeffs_map.begin(),
       t1w_cit = t1_wts_map.begin(),         t2w_cit  = t2_wts_map.begin(),
       p_cit   = set_partition_map.begin();  t1c_cit != t1_coeffs_map.end();
       ++t1c_cit, ++t2c_cit, ++t1w_cit, ++t2w_cit, ++p_cit) {
    const RealVector2DArray&  t1c = t1c_cit->second;
    const RealMatrix2DArray&  t2c = t2c_cit->second;
    const RealVector2DArray&  t1w = t1w_cit->second;
    const RealMatrix2DArray&  t2w = t2w_cit->second;
    const UShort2DArray& set_part = p_cit->second;
    num_lev = t1c.size();  partial = !set_part.empty();
    if (use_derivs) {
      for (lev=0; lev<num_lev; ++lev) {
	const RealVectorArray& t1c_l = t1c[lev];
	if (partial)  { start = set_part[lev][0]; end = set_part[lev][1]; }
	else          { start = 0;                end = t1c_l.size(); }
	for (set=start; set<end; ++set) {
	  const RealVector& t1c_ls = t1c_l[set];
	  const RealMatrix& t2c_ls = t2c[lev][set];
	  const RealVector& t1w_ls = t1w[lev][set];
	  const RealMatrix& t2w_ls = t2w[lev][set];
	  num_tp_pts = t1c_ls.length();
	  for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	    integral += t1c_ls[pt] * t1w_ls[pt];
	    const Real* t2c_lsp = t2c_ls[pt];
	    const Real* t2w_lsp = t2w_ls[pt];
	    for (v=0; v<num_v; ++v)
	      integral += t2c_lsp[v] * t2w_lsp[v];
	  }
	}
      }
    }
    else {
      for (lev=0; lev<num_lev; ++lev) {
	const RealVectorArray& t1c_l = t1c[lev];
	if (partial)  { start = set_part[lev][0]; end = set_part[lev][1]; }
	else          { start = 0;                end = t1c_l.size(); }
	for (set=start; set<end; ++set) {
	  const RealVector& t1c_ls = t1c_l[set];
	  const RealVector& t1w_ls = t1w[lev][set];
	  num_tp_pts = t1c_ls.length();
	  for (pt=0; pt<num_tp_pts; ++pt) // omitted if empty surplus vector
	    integral += t1c_ls[pt] * t1w_ls[pt];
	}
      }
    }
  }
  return integral;
}


Real HierarchInterpPolyApproximation::
expectation(const RealVector& x,
	    const std::map<UShortArray, RealVector2DArray>& t1_coeffs_map,
	    const std::map<UShortArray, RealMatrix2DArray>& t2_coeffs_map,
	    const std::map<UShortArray, UShort3DArray>& sm_mi_map,
	    const std::map<UShortArray, UShort4DArray>& colloc_key_map,
	    const std::map<UShortArray, UShort2DArray>& set_partition_map)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  Real integral = 0.;
  size_t lev, set, pt, num_lev, start, end, num_tp_pts, v,
    num_v = sharedDataRep->numVars;
  std::map<UShortArray, RealVector2DArray>::const_iterator t1c_cit;
  std::map<UShortArray, RealMatrix2DArray>::const_iterator t2c_cit;
  std::map<UShortArray, UShort3DArray>::const_iterator sm_cit;
  std::map<UShortArray, UShort4DArray>::const_iterator ck_cit;
  std::map<UShortArray, UShort2DArray>::const_iterator p_cit;
  bool partial, use_derivs = data_rep->basisConfigOptions.useDerivs;
  const SizetList&  r_ind = data_rep->randomIndices;
  const SizetList& nr_ind = data_rep->nonRandomIndices;
  for (t1c_cit = t1_coeffs_map.begin(),      t2c_cit  = t2_coeffs_map.begin(),
       sm_cit  = sm_mi_map.begin(),           ck_cit  = colloc_key_map.begin(),
       p_cit   = set_partition_map.begin();  t1c_cit != t1_coeffs_map.end();
       ++t1c_cit, ++t2c_cit, ++sm_cit, ++ck_cit, ++p_cit) {
    const RealVector2DArray&    t1c = t1c_cit->second;
    const RealMatrix2DArray&    t2c = t2c_cit->second;
    const UShort3DArray&      sm_mi =  sm_cit->second;
    const UShort4DArray& colloc_key =  ck_cit->second;
    const UShort2DArray&   set_part =   p_cit->second;
    num_lev = t1c.size();  partial = !set_part.empty();
    if (use_derivs) {
      for (lev=0; lev<num_lev; ++lev) {
	const RealVectorArray& t1c_l = t1c[lev];
	if (partial)  { start = set_part[lev][0]; end = set_part[lev][1]; }
	else          { start = 0;                end = t1c_l.size(); }
	for (set=start; set<end; ++set) {
	  const RealVector&    t1c_ls = t1c_l[set];
	  const RealMatrix&    t2c_ls = t2c[lev][set];
	  const UShortArray& sm_mi_ls = sm_mi[lev][set];
	  num_tp_pts = t1c_ls.length();
	  for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	    const UShortArray& key_lsp = colloc_key[lev][set][pt];
	    integral += t1c_ls[pt]
	      * data_rep->type1_interpolant_value(x, key_lsp, sm_mi_ls, nr_ind)
	      * data_rep->type1_weight(key_lsp, sm_mi_ls, r_ind);
	    const Real* t2c_lsp = t2c_ls[pt];
	    for (v=0; v<num_v; ++v)
	      integral += t2c_lsp[v]
		* data_rep->type2_interpolant_value(x,v,key_lsp,sm_mi_ls,nr_ind)
		* data_rep->type2_weight(v, key_lsp, sm_mi_ls, r_ind);
	  }
	}
      }
    }
    else {
      for (lev=0; lev<num_lev; ++lev) {
	const RealVectorArray& t1c_l = t1c[lev];
	if (partial)  { start = set_part[lev][0]; end = set_part[lev][1]; }
	else          { start = 0;                end = t1c_l.size(); }
	for (set=start; set<end; ++set) {
	  const RealVector&    t1c_ls = t1c_l[set];
	  const UShortArray& sm_mi_ls = sm_mi[lev][set];
	  num_tp_pts = t1c_ls.length();
	  for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	    const UShortArray& key_lsp = colloc_key[lev][set][pt];
	    integral += t1c_ls[pt]
	      * data_rep->type1_interpolant_value(x, key_lsp, sm_mi_ls, nr_ind)
	      * data_rep->type1_weight(key_lsp, sm_mi_ls, r_ind);
	  }
	}
      }
    }
  }
  return integral;
}


/** For inserted/augmented design/epistemic variables in standard mode. */
const RealVector& HierarchInterpPolyApproximation::
expectation_gradient(const RealMatrix2DArray& t1_coeff_grads,
		     const RealVector2DArray& t1_wts)
{
  size_t lev, num_levels = t1_coeff_grads.size(), set, num_sets, pt, num_tp_pts,
    v, num_deriv_vars = t1_coeff_grads[0][0].numRows();
  if (approxGradient.length() != num_deriv_vars)
    approxGradient.sizeUninitialized(num_deriv_vars);
  approxGradient = 0.;

  for (lev=0; lev<num_levels; ++lev) {
    const RealMatrixArray& t1_coeff_grads_l = t1_coeff_grads[lev];
    num_sets = t1_coeff_grads_l.size();
    for (set=0; set<num_sets; ++set) {
      const RealMatrix& t1_coeff_grads_ls = t1_coeff_grads_l[set];
      num_tp_pts = t1_coeff_grads_ls.numCols();
      for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	const Real* t1_coeff_grads_lsp = t1_coeff_grads_ls[pt];
	Real t1_wt_lsp = t1_wts[lev][set][pt];
	for (v=0; v<num_deriv_vars; ++v)
	  approxGradient[v] += t1_coeff_grads_lsp[v] * t1_wt_lsp;
      }
    }
  }

  return approxGradient;
}


/** For inserted design/epistemic variables in all_variables mode. */
Real HierarchInterpPolyApproximation::
expectation_gradient(const RealVector& x,
		     const RealMatrix2DArray& t1_coeff_grads,
		     const UShort3DArray& sm_mi,
		     const UShort4DArray& colloc_key, size_t t1cg_index)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  size_t lev, num_levels = t1_coeff_grads.size(), set, num_sets, pt, num_tp_pts;
  const SizetList&  r_ind = data_rep->randomIndices;
  const SizetList& nr_ind = data_rep->nonRandomIndices;

  Real grad = 0.;
  for (lev=0; lev<num_levels; ++lev) {
    const RealMatrixArray& t1_coeff_grads_l = t1_coeff_grads[lev];
    num_sets = t1_coeff_grads_l.size();
    for (set=0; set<num_sets; ++set) {
      const RealMatrix& t1_coeff_grads_ls = t1_coeff_grads_l[set];
      const UShortArray&         sm_mi_ls = sm_mi[lev][set];
      num_tp_pts = t1_coeff_grads_ls.numCols();
      for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	const UShortArray& key_lsp = colloc_key[lev][set][pt];
	grad += t1_coeff_grads_ls(t1cg_index, pt)
	  * data_rep->type1_interpolant_value(x, key_lsp, sm_mi_ls, nr_ind)
	  * data_rep->type1_weight(key_lsp, sm_mi_ls, r_ind);
      }
    }
  }
  return grad;
}


/** For augmented design/epistemic variables in all_variables mode. */
Real HierarchInterpPolyApproximation::
expectation_gradient(const RealVector& x, const RealVector2DArray& t1_coeffs,
		     const RealMatrix2DArray& t2_coeffs,
		     const UShort3DArray& sm_mi,
		     const UShort4DArray& colloc_key, size_t deriv_index)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  size_t lev, num_levels = t1_coeffs.size(), set, num_sets, pt, num_tp_pts, v,
    num_v = sharedDataRep->numVars;
  const SizetList&  r_ind = data_rep->randomIndices;
  const SizetList& nr_ind = data_rep->nonRandomIndices;

  Real grad = 0.;
  for (lev=0; lev<num_levels; ++lev) {
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    num_sets = t1_coeffs_l.size();
    for (set=0; set<num_sets; ++set) {
      const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
      const UShortArray&    sm_mi_ls = sm_mi[lev][set];
      num_tp_pts = t1_coeffs_ls.length();
      for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	const UShortArray& key_lsp = colloc_key[lev][set][pt];
	grad += t1_coeffs_ls[pt] * data_rep->
	  type1_interpolant_gradient(x, deriv_index, key_lsp, sm_mi_ls,
				     data_rep->nonRandomIndices) *
	  data_rep->type1_weight(key_lsp, sm_mi_ls, data_rep->randomIndices);
	if (data_rep->basisConfigOptions.useDerivs) {
	  const Real *t2_coeff_lsp = t2_coeffs[lev][set][pt];
	  for (v=0; v<num_v; ++v)
	    grad += t2_coeff_lsp[v] *
	      data_rep->type2_interpolant_gradient(x, deriv_index, v, key_lsp,
						   sm_mi_ls, nr_ind) *
	      data_rep->type2_weight(v, key_lsp, sm_mi_ls, r_ind);
	}
      }
    }
  }

  return grad;
}


/** This version accesses the surrogate data for efficiency in cases where
    the r1 and r2 interpolants correspond to the data and do not need to be
    evaluated.  Note: whereas expectation() supports either a reference or
    increment key (passed as generic set_partition), functions forming
    hierarchical product interpolants only support a reference key due to
    nonlinearity (only end point can be controlled).  Note: could consider 
    computing deltaR1R2 = R1 deltaR2 + R2 deltaR1 + deltaR1 deltaR2. */
void HierarchInterpPolyApproximation::
product_interpolant(HierarchInterpPolyApproximation* hip_approx_2,
		    RealVector2DArray& r1r2_t1_coeffs,
		    RealMatrix2DArray& r1r2_t2_coeffs,
		    const UShort2DArray& ref_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver   = data_rep->hsg_driver();
  const UShort3DArray&      sm_mi        = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      colloc_key   = hsg_driver->collocation_key();
  const Sizet3DArray&       colloc_index = hsg_driver->collocation_indices();
  size_t lev, set, pt, num_levels = expT1CoeffsIter->second.size(),
    num_sets, num_tp_pts, cntr = 0, c_index, v, num_v = sharedDataRep->numVars;
  bool partial = !ref_key.empty();
  const SDVArray& sdv_array = surrData.variables_data();
  const SDRArray& sdr_array_1 = surrData.response_data();
  const SDRArray& sdr_array_2 = hip_approx_2->surrData.response_data();

  // form hierarchical t1/t2 coeffs for raw moment R1 R2
  r1r2_t1_coeffs.resize(num_levels); r1r2_t1_coeffs[0].resize(1);
  r1r2_t2_coeffs.resize(num_levels); r1r2_t2_coeffs[0].resize(1);
  r1r2_t1_coeffs[0][0].sizeUninitialized(1);
  
  // level 0 (assume this is always contained in a partial ref_key)
  c_index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
  const SurrogateDataResp& sdr0_1 = sdr_array_1[c_index];
  const SurrogateDataResp& sdr0_2 = sdr_array_2[c_index];
  Real data_fn1 = sdr0_1.response_function(),
       data_fn2 = sdr0_2.response_function();
  r1r2_t1_coeffs[0][0][0] = data_fn1 * data_fn2;
  ++cntr;

  if (data_rep->basisConfigOptions.useDerivs) {
    // level 0 continued
    r1r2_t2_coeffs[0][0].shapeUninitialized(num_v, 1);
    Real *r1r2_t2_coeffs_000 = r1r2_t2_coeffs[0][0][0];
    const RealVector& data_grad1 = sdr0_1.response_gradient();
    const RealVector& data_grad2 = sdr0_2.response_gradient();
    for (v=0; v<num_v; ++v)
      r1r2_t2_coeffs_000[v]
	= data_fn1 * data_grad2[v] + data_fn2 * data_grad1[v];
  
    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = (partial) ? ref_key[lev][1] : colloc_key[lev].size();
      r1r2_t1_coeffs[lev].resize(num_sets);
      r1r2_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	RealVector& r1r2_t1_coeffs_ls = r1r2_t1_coeffs[lev][set];
	RealMatrix& r1r2_t2_coeffs_ls = r1r2_t2_coeffs[lev][set];
	num_tp_pts = colloc_key[lev][set].size();
	r1r2_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	r1r2_t2_coeffs_ls.shapeUninitialized(num_v, num_tp_pts);
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	  c_index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	  const RealVector& c_vars = sdv_array[c_index].continuous_variables();
	  const SurrogateDataResp& sdr_1 = sdr_array_1[c_index];
	  const SurrogateDataResp& sdr_2 = sdr_array_2[c_index];
	  // type1 hierarchical interpolation of R1 R2
	  data_fn1 = sdr_1.response_function();
	  data_fn2 = sdr_2.response_function();
	  r1r2_t1_coeffs_ls[pt] = data_fn1 * data_fn2
	    - value(c_vars, sm_mi, colloc_key, r1r2_t1_coeffs,
		    r1r2_t2_coeffs, lev-1, ref_key);
	  // type2 hierarchical interpolation of R1 R2
	  // --> interpolated grads are R1 * R2' + R2 * R1'
	  Real* r1r2_t2_coeffs_lsp = r1r2_t2_coeffs_ls[pt];
	  const RealVector& data_grad1 = sdr_1.response_gradient();
	  const RealVector& data_grad2 = sdr_2.response_gradient();
	  const RealVector& prev_grad  = gradient_basis_variables(c_vars,
	    sm_mi, colloc_key, r1r2_t1_coeffs, r1r2_t2_coeffs, lev-1, ref_key);
	  for (v=0; v<num_v; ++v)
	    r1r2_t2_coeffs_lsp[v] = data_fn1 * data_grad2[v]
	      + data_fn2 * data_grad1[v] - prev_grad[v];
	}
      }
    }
  }
  else {
    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = (partial) ? ref_key[lev][1] : colloc_key[lev].size();
      r1r2_t1_coeffs[lev].resize(num_sets);
      r1r2_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	num_tp_pts = colloc_key[lev][set].size();
	RealVector& r1r2_t1_coeffs_ls = r1r2_t1_coeffs[lev][set];
	r1r2_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	// type1 hierarchical interpolation of R1 R2
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	  c_index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	  r1r2_t1_coeffs_ls[pt] = sdr_array_1[c_index].response_function()
	    * sdr_array_2[c_index].response_function()
	    - value(sdv_array[c_index].continuous_variables(), sm_mi,
		    colloc_key, r1r2_t1_coeffs, r1r2_t2_coeffs, lev-1, ref_key);
	}
      }
    }
  }
}


/** This version evaluates the r1 and r2 interpolants for cases where a
    particular key of the surrogate data does not correspond (e.g., after
    combine_coefficients).  Whereas expectation() supports either a
    reference or increment key (passed as generic set_partition), functions
    forming hierarchical product interpolants only support a reference key
    due to nonlinearity (only end point can be controlled).  Note: consider 
    computing deltaR1R2 = R1 deltaR2 + R2 deltaR1 + deltaR1 deltaR2. */
// *** For use with a single set of combined coefficients; this would be a
//     model discrepancy for key > first ***
void HierarchInterpPolyApproximation::
product_interpolant(const RealMatrix2DArray& var_sets,
		    const UShort3DArray&     sm_mi,
		    const UShort4DArray&     colloc_key,
		    const RealVector2DArray& r1_t1_coeffs,
		    const RealMatrix2DArray& r1_t2_coeffs,
		    const RealVector2DArray& r2_t1_coeffs,
		    const RealMatrix2DArray& r2_t2_coeffs, bool same,
		    RealVector2DArray& r1r2_t1_coeffs,
		    RealMatrix2DArray& r1r2_t2_coeffs,
		    const UShort2DArray& ref_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  size_t lev, set, pt, num_levels = r1_t1_coeffs.size(), num_sets, num_tp_pts,
    v, num_v = sharedDataRep->numVars;
  bool partial = !ref_key.empty();

  // form hierarchical t1/t2 coeffs for raw moment R1 R2
  r1r2_t1_coeffs.resize(num_levels); r1r2_t1_coeffs[0].resize(1);
  r1r2_t2_coeffs.resize(num_levels); r1r2_t2_coeffs[0].resize(1);
  r1r2_t1_coeffs[0][0].sizeUninitialized(1);

  // level 0 (assume this is always contained in a partial ref_key)
  RealVector c_vars(Teuchos::View, const_cast<Real*>(var_sets[0][0][0]),
		    (int)num_v);
  // type1 hierarchical interpolation of R1 R2
  Real r1_val =
    value(c_vars, sm_mi, colloc_key, r1_t1_coeffs, r1_t2_coeffs, 0, ref_key),
       r2_val = (same) ? r1_val :
    value(c_vars, sm_mi, colloc_key, r2_t1_coeffs, r2_t2_coeffs, 0, ref_key);
  r1r2_t1_coeffs[0][0][0] = r1_val * r2_val;

  if (data_rep->basisConfigOptions.useDerivs) {
    // level 0 continued: type2 hierarchical interpolation of R1 R2
    const RealVector& r1_grad = gradient_basis_variables(c_vars, sm_mi,
      colloc_key, r1_t1_coeffs, r1_t2_coeffs, 0, ref_key);
    const RealVector& r2_grad = (same) ? r1_grad : gradient_basis_variables(
      c_vars, sm_mi, colloc_key, r2_t1_coeffs, r2_t2_coeffs, 0, ref_key);
    r1r2_t2_coeffs[0][0].shapeUninitialized(num_v, 1);
    Real *r1r2_t2_coeffs_000 = r1r2_t2_coeffs[0][0][0];
    for (v=0; v<num_v; ++v)
      r1r2_t2_coeffs_000[v] = r1_val * r2_grad[v] + r2_val * r1_grad[v];

    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = (partial) ? ref_key[lev][1] : colloc_key[lev].size();
      r1r2_t1_coeffs[lev].resize(num_sets);
      r1r2_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	RealVector& r1r2_t1_coeffs_ls = r1r2_t1_coeffs[lev][set];
	RealMatrix& r1r2_t2_coeffs_ls = r1r2_t2_coeffs[lev][set];
	num_tp_pts = colloc_key[lev][set].size();
	r1r2_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	r1r2_t2_coeffs_ls.shapeUninitialized(num_v, num_tp_pts);
	for (pt=0; pt<num_tp_pts; ++pt) {
	  RealVector c_vars(Teuchos::View,
	    const_cast<Real*>(var_sets[lev][set][pt]), (int)num_v);
	  // type1 hierarchical interpolation of R1 R2
	  r1_val =
	    value(c_vars, sm_mi, colloc_key, r1_t1_coeffs, r1_t2_coeffs,
		  lev, ref_key);
	  r2_val = (same) ? r1_val :
	    value(c_vars, sm_mi, colloc_key, r2_t1_coeffs, r2_t2_coeffs,
		  lev, ref_key);
	  r1r2_t1_coeffs_ls[pt] = r1_val * r2_val - value(c_vars, sm_mi,
	    colloc_key, r1r2_t1_coeffs, r1r2_t2_coeffs, lev-1, ref_key);
	  // type2 hierarchical interpolation of R1 R2
	  // --> interpolated grads are R1 * R2' + R2 * R1'
	  Real* r1r2_t2_coeffs_lsp = r1r2_t2_coeffs_ls[pt];
	  const RealVector& r1_grad = gradient_basis_variables(c_vars, sm_mi,
	    colloc_key, r1_t1_coeffs, r1_t2_coeffs, lev, ref_key);
	  const RealVector& r2_grad = (same) ? r1_grad :
	    gradient_basis_variables(c_vars, sm_mi, colloc_key, r2_t1_coeffs,
	    r2_t2_coeffs, lev, ref_key);
	  const RealVector& prev_grad = gradient_basis_variables(c_vars, sm_mi,
	    colloc_key, r1r2_t1_coeffs, r1r2_t2_coeffs, lev-1, ref_key);
	  for (v=0; v<num_v; ++v)
	    r1r2_t2_coeffs_lsp[v]
	      = r1_val * r2_grad[v] + r2_val * r1_grad[v] - prev_grad[v];
	}
      }
    }
  }
  else {
    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = (partial) ? ref_key[lev][1] : colloc_key[lev].size();
      r1r2_t1_coeffs[lev].resize(num_sets);
      r1r2_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	num_tp_pts = colloc_key[lev][set].size();
	RealVector& r1r2_t1_coeffs_ls = r1r2_t1_coeffs[lev][set];
	r1r2_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	// type1 hierarchical interpolation of R1 R2
	for (pt=0; pt<num_tp_pts; ++pt) {
	  RealVector c_vars(Teuchos::View,
	    const_cast<Real*>(var_sets[lev][set][pt]), (int)num_v);
	  // type1 hierarchical interpolation of R1 R2
	  r1_val = value(c_vars, sm_mi, colloc_key, r1_t1_coeffs,
	    r1_t2_coeffs, lev, ref_key);
	  r2_val = (same) ? r1_val : value(c_vars, sm_mi, colloc_key,
	    r2_t1_coeffs, r2_t2_coeffs, lev, ref_key);
	  r1r2_t1_coeffs_ls[pt] = r1_val * r2_val - value(c_vars, sm_mi,
	    colloc_key, r1r2_t1_coeffs, r1r2_t2_coeffs, lev-1, ref_key);
	}
      }
    }
  }
}


/*
void HierarchInterpPolyApproximation::
product_interpolant(HierarchInterpPolyApproximation* hip_approx_2,
  std::map<UShortArray, RealVector2DArray>& r1r2_t1c_map,
  std::map<UShortArray, RealMatrix2DArray>& r1r2_t2c_map,
  const std::map<UShortArray, UShort2DArray>& ref_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver   = data_rep->hsg_driver();
  const std::map<UShortArray, UShort3DArray>& sm_mi_map
    = hsg_driver->smolyak_multi_index_map();
  const std::map<UShortArray, UShort4DArray>& colloc_key_map
    = hsg_driver->collocation_key_map();
  const std::map<UShortArray, Sizet3DArray>&  colloc_index_map
    = hsg_driver->collocation_indices_map();
  const std::map<UShortArray, SDVArray>& sdv_array_map
    = surrData.variables_data_map();
  const std::map<UShortArray, Sizet3DArray>& c_index_map
    = hsg_driver->collocation_indices_map();
  size_t mod, lev, set, pt, num_levels, num_sets, num_tp_pts, v,
    num_v = sharedDataRep->numVars;
  bool same = (hip_approx_2 == this), partial = !ref_key.empty();

  // form hierarchical t1/t2 coeffs for raw moment R1 R2
  RealVector2DArray r1r2_t1c(num_levels); r1r2_t1c[0].resize(1);
  RealMatrix2DArray r1r2_t2c(num_levels); r1r2_t2c[0].resize(1);
  r1r2_t1c[0][0].sizeUninitialized(1);

  std::map<UShortArray, RealVector2DArray>::iterator r1_t1_it, r2_t1_it;
  std::map<UShortArray, RealMatrix2DArray>::iterator r1_t2_it, r2_t2_it;
  std::map<UShortArray, UShort4DArray>::iterator        ck_it;
  std::map<UShortArray, Sizet3DArray>::iterator         ci_it;
  std::map<UShortArray, SDVArray>::iterator           sdva_it;
  for (r1_t1_it  =               expansionType1Coeffs.begin(),
       r1_t2_it  =               expansionType2Coeffs.begin(),
       r2_t1_it  = hip_approx_2->expansionType1Coeffs.begin(),
       r2_t2_it  = hip_approx_2->expansionType2Coeffs.begin(),
       ck_it     = colloc_key_map.begin(), ci_it = c_index_map.begin(),
       sdva_it = sdv_array_map.begin(), mod=0;
       r1_t1_it !=               expansionType1Coeffs.end();
       ++r1_t1_it, ++r1_t2_it, ++r2_t1_it, ++r2_t2_it,
       ++ci_it, ++ck_it, ++sdva_it, ++mod) {
    const UShortArray&           key = r1_t1_it->first;
    const RealVector2DArray&  r1_t1c = r1_t1_it->second;
    const SDVArray&        sdv_array =  sdva_it->second;
    const Sizet3DArray& colloc_index =    ci_it->second;
    const UShort4DArray&  colloc_key =    ck_it->second;
    // level 0, model 0 could use surrogate data directly for R ...
    // *** TO DO ?

    // levels 1:L and models 1:M need to aggregate discrpenacies to predict R
    // --> use value() to reliably aggregate
    num_lev = r1_t1c.size();
    for (lev=0, cntr=0; lev<num_lev; ++lev) {
      const RealVectorArray& r1_t1c_l = r1_t1c[lev];
      const RealMatrixArray& r1_t2c_l = r1_t2_it->second[lev];
      const RealVectorArray& r2_t1c_l = r2_t1_it->second[lev];
      const RealMatrixArray& r2_t2c_l = r2_t2_it->second[lev];
      num_sets = r1_t1c_l.size();
      for (set=0; set<num_sets; ++set) {
	num_tp_pts = colloc_key[lev][set].size();
	RealVector& r1r2_t1_coeffs_ls = r1r2_t1_coeffs[lev][set];
	r1r2_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	// type1 hierarchical interpolation of R1 R2
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	  c_index = (c_index.empty()) ? cntr : colloc_index[lev][set][pt];
	  const RealVector& c_vars = sdv_array[c_index].continuous_variables();
	  // type1 hierarchical interpolation of R1 R2
	  r1_val = value(c_vars, sm_mi, colloc_key, r1_t1_coeffs,
	    r1_t2_coeffs, mod, lev, ref_key);
	  r2_val = (same) ? r1_val : value(c_vars, sm_mi, colloc_key,
	    r2_t1_coeffs, r2_t2_coeffs, mod, lev, ref_key);
	  r1r2_t1_coeffs_ls[pt] = r1_val * r2_val - value(c_vars, sm_mi,
	    colloc_key, r1r2_t1_coeffs, r1r2_t2_coeffs, mod, lev-1, ref_key);
	}
      }

      }
    }
    r1r2_t1c_map[key] = r1r2_t1c;
    r1r2_t2c_map[key] = r1r2_t2c;
  }

  // level 0 (assume this is always contained in a partial ref_key)
  RealVector c_vars(Teuchos::View, const_cast<Real*>(var_sets[0][0][0]),
		    (int)num_v);
  // type1 hierarchical interpolation of R1 R2
  Real r1_val =
    value(c_vars, sm_mi, colloc_key, r1_t1_coeffs, r1_t2_coeffs, 0, ref_key),
       r2_val = (same) ? r1_val :
    value(c_vars, sm_mi, colloc_key, r2_t1_coeffs, r2_t2_coeffs, 0, ref_key);
  r1r2_t1_coeffs[0][0][0] = r1_val * r2_val;

  if (data_rep->basisConfigOptions.useDerivs) {
    // level 0 continued: type2 hierarchical interpolation of R1 R2
    const RealVector& r1_grad = gradient_basis_variables(c_vars, sm_mi,
      colloc_key, r1_t1_coeffs, r1_t2_coeffs, 0, ref_key);
    const RealVector& r2_grad = (same) ? r1_grad : gradient_basis_variables(
      c_vars, sm_mi, colloc_key, r2_t1_coeffs, r2_t2_coeffs, 0, ref_key);
    r1r2_t2_coeffs[0][0].shapeUninitialized(num_v, 1);
    Real *r1r2_t2_coeffs_000 = r1r2_t2_coeffs[0][0][0];
    for (v=0; v<num_v; ++v)
      r1r2_t2_coeffs_000[v] = r1_val * r2_grad[v] + r2_val * r1_grad[v];

    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = (partial) ? ref_key[lev][1] : colloc_key[lev].size();
      r1r2_t1_coeffs[lev].resize(num_sets);
      r1r2_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	RealVector& r1r2_t1_coeffs_ls = r1r2_t1_coeffs[lev][set];
	RealMatrix& r1r2_t2_coeffs_ls = r1r2_t2_coeffs[lev][set];
	num_tp_pts = colloc_key[lev][set].size();
	r1r2_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	r1r2_t2_coeffs_ls.shapeUninitialized(num_v, num_tp_pts);
	for (pt=0; pt<num_tp_pts; ++pt) {
	  RealVector c_vars(Teuchos::View,
	    const_cast<Real*>(var_sets[lev][set][pt]), (int)num_v);
	  // type1 hierarchical interpolation of R1 R2
	  r1_val =
	    value(c_vars, sm_mi, colloc_key, r1_t1_coeffs, r1_t2_coeffs,
		  lev, ref_key);
	  r2_val = (same) ? r1_val :
	    value(c_vars, sm_mi, colloc_key, r2_t1_coeffs, r2_t2_coeffs,
		  lev, ref_key);
	  r1r2_t1_coeffs_ls[pt] = r1_val * r2_val - value(c_vars, sm_mi,
	    colloc_key, r1r2_t1_coeffs, r1r2_t2_coeffs, lev-1, ref_key);
	  // type2 hierarchical interpolation of R1 R2
	  // --> interpolated grads are R1 * R2' + R2 * R1'
	  Real* r1r2_t2_coeffs_lsp = r1r2_t2_coeffs_ls[pt];
	  const RealVector& r1_grad = gradient_basis_variables(c_vars, sm_mi,
	    colloc_key, r1_t1_coeffs, r1_t2_coeffs, lev, ref_key);
	  const RealVector& r2_grad = (same) ? r1_grad :
	    gradient_basis_variables(c_vars, sm_mi, colloc_key, r2_t1_coeffs,
	    r2_t2_coeffs, lev, ref_key);
	  const RealVector& prev_grad = gradient_basis_variables(c_vars, sm_mi,
	    colloc_key, r1r2_t1_coeffs, r1r2_t2_coeffs, lev-1, ref_key);
	  for (v=0; v<num_v; ++v)
	    r1r2_t2_coeffs_lsp[v]
	      = r1_val * r2_grad[v] + r2_val * r1_grad[v] - prev_grad[v];
	}
      }
    }
  }
  else {
    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = (partial) ? ref_key[lev][1] : colloc_key[lev].size();
      r1r2_t1_coeffs[lev].resize(num_sets);
      r1r2_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	num_tp_pts = colloc_key[lev][set].size();
	RealVector& r1r2_t1_coeffs_ls = r1r2_t1_coeffs[lev][set];
	r1r2_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	// type1 hierarchical interpolation of R1 R2
	for (pt=0; pt<num_tp_pts; ++pt) {
	  RealVector c_vars(Teuchos::View,
	    const_cast<Real*>(var_sets[lev][set][pt]), (int)num_v);
	  // type1 hierarchical interpolation of R1 R2
	  r1_val = value(c_vars, sm_mi, colloc_key, r1_t1_coeffs,
	    r1_t2_coeffs, lev, ref_key);
	  r2_val = (same) ? r1_val : value(c_vars, sm_mi, colloc_key,
	    r2_t1_coeffs, r2_t2_coeffs, lev, ref_key);
	  r1r2_t1_coeffs_ls[pt] = r1_val * r2_val - value(c_vars, sm_mi,
	    colloc_key, r1r2_t1_coeffs, r1r2_t2_coeffs, lev-1, ref_key);
	}
      }
    }
  }
}
*/
  

/** This version accesses the surrogate data for efficiency in cases where
    the r1 and r2 interpolants correspond to the data and do not need to be
    evaluated.  Note: whereas expectation() supports either a reference or
    increment key (passed as generic set_partition), functions forming
    hierarchical product interpolants only support a reference key due to
    nonlinearity (only end point can be controlled).  Note: could consider 
    computing deltaR1R2 = R1 deltaR2 + R2 deltaR1 + deltaR1 deltaR2. */
void HierarchInterpPolyApproximation::
central_product_interpolant(HierarchInterpPolyApproximation* hip_approx_2,
			    Real mean_1, Real mean_2,
			    RealVector2DArray& cov_t1_coeffs,
			    RealMatrix2DArray& cov_t2_coeffs,
			    const UShort2DArray& ref_key)
{
  // form hierarchical t1/t2 coeffs for (R_1 - \mu_1) (R_2 - \mu_2)
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver   = data_rep->hsg_driver();
  const UShort3DArray&      sm_mi        = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      colloc_key   = hsg_driver->collocation_key();
  const Sizet3DArray&       colloc_index = hsg_driver->collocation_indices();
  size_t lev, set, pt, num_levels = expT1CoeffsIter->second.size(),
    num_sets, num_tp_pts, cntr = 0, c_index, v, num_v = sharedDataRep->numVars;
  bool partial = !ref_key.empty();
  const SDVArray& sdv_array = surrData.variables_data();
  const SDRArray& sdr_array_1 = surrData.response_data();
  const SDRArray& sdr_array_2 = hip_approx_2->surrData.response_data();

  cov_t1_coeffs.resize(num_levels); cov_t1_coeffs[0].resize(1);
  cov_t2_coeffs.resize(num_levels); cov_t2_coeffs[0].resize(1);
  cov_t1_coeffs[0][0].sizeUninitialized(1);

  // level 0 (assume this is always contained in a partial ref_key)
  c_index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
  const SurrogateDataResp& sdr0_1 = sdr_array_1[c_index];
  const SurrogateDataResp& sdr0_2 = sdr_array_2[c_index];
  Real data_fn1_mm1 = sdr0_1.response_function() - mean_1,
       data_fn2_mm2 = sdr0_2.response_function() - mean_2;
  cov_t1_coeffs[0][0][0] = data_fn1_mm1 * data_fn2_mm2;
  ++cntr;

  if (data_rep->basisConfigOptions.useDerivs) {
    // level 0 continued
    cov_t2_coeffs[0][0].shapeUninitialized(num_v, 1);
    Real *cov_t2_coeffs_000 = cov_t2_coeffs[0][0][0];
    const RealVector& data_grad1 = sdr0_1.response_gradient();
    const RealVector& data_grad2 = sdr0_2.response_gradient();
    for (v=0; v<num_v; ++v)
      cov_t2_coeffs_000[v]
	= data_fn1_mm1 * data_grad2[v] + data_fn2_mm2 * data_grad1[v];

    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = (partial) ? ref_key[lev][1] : colloc_key[lev].size();
      cov_t1_coeffs[lev].resize(num_sets); cov_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	RealVector& cov_t1_coeffs_ls = cov_t1_coeffs[lev][set];
	RealMatrix& cov_t2_coeffs_ls = cov_t2_coeffs[lev][set];
	num_tp_pts = colloc_key[lev][set].size();
	cov_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	cov_t2_coeffs_ls.shapeUninitialized(num_v, num_tp_pts);
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	  c_index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	  const RealVector& c_vars = sdv_array[c_index].continuous_variables();
	  const SurrogateDataResp& sdr_1 = sdr_array_1[c_index];
	  const SurrogateDataResp& sdr_2 = sdr_array_2[c_index];
	  // type1 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	  data_fn1_mm1 = sdr_1.response_function() - mean_1;
	  data_fn2_mm2 = sdr_2.response_function() - mean_2;
	  cov_t1_coeffs_ls[pt] = data_fn1_mm1 * data_fn2_mm2 -
	    value(c_vars, sm_mi, colloc_key, cov_t1_coeffs, cov_t2_coeffs,
		  lev-1, ref_key);
	  // type2 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	  // --> interpolated grads are (R_1-\mu_1) * R_2' + (R_2-\mu_2) * R_1'
	  Real* cov_t2_coeffs_lsp = cov_t2_coeffs_ls[pt];
	  const RealVector& data_grad1 = sdr_1.response_gradient();
	  const RealVector& data_grad2 = sdr_2.response_gradient();
	  const RealVector& prev_grad  = gradient_basis_variables(c_vars,
	    sm_mi, colloc_key, cov_t1_coeffs, cov_t2_coeffs, lev-1, ref_key);
	  for (v=0; v<num_v; ++v)
	    cov_t2_coeffs_lsp[v] = data_fn1_mm1 * data_grad2[v]
	      + data_fn2_mm2 * data_grad1[v] - prev_grad[v];
	}
      }
    }
  }
  else {
    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = (partial) ? ref_key[lev][1] : colloc_key[lev].size();
      cov_t1_coeffs[lev].resize(num_sets); cov_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	num_tp_pts = colloc_key[lev][set].size();
	RealVector& cov_t1_coeffs_ls = cov_t1_coeffs[lev][set];
	cov_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	// type1 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	  c_index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	  cov_t1_coeffs_ls[pt]
	    = (sdr_array_1[c_index].response_function() - mean_1)
	    * (sdr_array_2[c_index].response_function() - mean_2)
	    - value(sdv_array[c_index].continuous_variables(), sm_mi,
		    colloc_key, cov_t1_coeffs, cov_t2_coeffs, lev-1, ref_key);
	}
      }
    }
  }
}


/** This version evaluates the r1 and r2 interpolants for cases where a
    particular key of the surrogate data does not correspond (e.g., after
    combine_coefficients).  Note: whereas expectation() supports either a
    reference or increment key (passed as generic set_partition), functions
    forming hierarchical product interpolants only support a reference key
    due to nonlinearity (only end point can be controlled).  Note: consider 
    computing deltaR1R2 = R1 deltaR2 + R2 deltaR1 + deltaR1 deltaR2. */
void HierarchInterpPolyApproximation::
central_product_interpolant(const RealMatrix2DArray& var_sets,
			    const UShort3DArray&     sm_mi,
			    const UShort4DArray&     colloc_key,
			    const RealVector2DArray& r1_t1_coeffs,
			    const RealMatrix2DArray& r1_t2_coeffs,
			    const RealVector2DArray& r2_t1_coeffs,
			    const RealMatrix2DArray& r2_t2_coeffs,
			    bool same, Real mean_r1, Real mean_r2,
			    RealVector2DArray& cov_t1_coeffs,
			    RealMatrix2DArray& cov_t2_coeffs,
			    const UShort2DArray& ref_key)
{
  // form hierarchical t1/t2 coeffs for (R_1 - \mu_1) (R_2 - \mu_2)
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  size_t lev, set, pt, num_levels = r1_t1_coeffs.size(), num_sets, num_tp_pts,
    v, num_v = sharedDataRep->numVars;
  bool partial = !ref_key.empty();

  cov_t1_coeffs.resize(num_levels); cov_t1_coeffs[0].resize(1);
  cov_t2_coeffs.resize(num_levels); cov_t2_coeffs[0].resize(1);
  cov_t1_coeffs[0][0].sizeUninitialized(1);

  // level 0 (assume this is always contained in a partial ref_key)
  RealVector c_vars(Teuchos::View, const_cast<Real*>(var_sets[0][0][0]),
		    (int)num_v);
  Real r1_val_mm1 = value(c_vars, sm_mi, colloc_key, r1_t1_coeffs,
	 r1_t2_coeffs, 0, ref_key) - mean_r1,
       r2_val_mm2 = (same) ? r1_val_mm1 : value(c_vars, sm_mi, colloc_key,
	 r2_t1_coeffs, r2_t2_coeffs, 0, ref_key) - mean_r2;
  cov_t1_coeffs[0][0][0] = r1_val_mm1 * r2_val_mm2;

  if (data_rep->basisConfigOptions.useDerivs) {
    // level 0 continued
    cov_t2_coeffs[0][0].shapeUninitialized(num_v, 1);
    Real *cov_t2_coeffs_000 = cov_t2_coeffs[0][0][0];
    const RealVector& r1_grad = gradient_basis_variables(c_vars, sm_mi,
      colloc_key, r1_t1_coeffs, r1_t2_coeffs, 0, ref_key);
    const RealVector& r2_grad = (same) ? r1_grad : gradient_basis_variables(
      c_vars, sm_mi, colloc_key, r2_t1_coeffs, r2_t2_coeffs, 0, ref_key);
    for (v=0; v<num_v; ++v)
      cov_t2_coeffs_000[v] = r1_val_mm1 * r2_grad[v] + r2_val_mm2 * r1_grad[v];

    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = (partial) ? ref_key[lev][1] : colloc_key[lev].size();
      cov_t1_coeffs[lev].resize(num_sets); cov_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	RealVector& cov_t1_coeffs_ls = cov_t1_coeffs[lev][set];
	RealMatrix& cov_t2_coeffs_ls = cov_t2_coeffs[lev][set];
	num_tp_pts = colloc_key[lev][set].size();
	cov_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	cov_t2_coeffs_ls.shapeUninitialized(num_v, num_tp_pts);
	for (pt=0; pt<num_tp_pts; ++pt) {
	  RealVector c_vars(Teuchos::View,
	    const_cast<Real*>(var_sets[lev][set][pt]), (int)num_v);
	  // type1 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	  r1_val_mm1 = value(c_vars, sm_mi, colloc_key, r1_t1_coeffs,
	    r1_t2_coeffs, lev, ref_key) - mean_r1;
	  r2_val_mm2 = (same) ? r1_val_mm1 : value(c_vars, sm_mi, colloc_key,
	    r2_t1_coeffs, r2_t2_coeffs, lev, ref_key) - mean_r2;
	  cov_t1_coeffs_ls[pt] = r1_val_mm1 * r2_val_mm2 - value(c_vars, sm_mi,
	    colloc_key, cov_t1_coeffs, cov_t2_coeffs, lev-1, ref_key);
	  // type2 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	  // --> interpolated grads are (R_1-\mu_1) * R_2' + (R_2-\mu_2) * R_1'
	  Real* cov_t2_coeffs_lsp = cov_t2_coeffs_ls[pt];
	  const RealVector& r1_grad = gradient_basis_variables(c_vars, sm_mi,
	    colloc_key, r1_t1_coeffs, r1_t2_coeffs, lev, ref_key);
	  const RealVector& r2_grad = (same) ? r1_grad :
	    gradient_basis_variables(c_vars, sm_mi, colloc_key, r2_t1_coeffs,
	    r2_t2_coeffs, lev, ref_key);
	  const RealVector& prev_grad = gradient_basis_variables(c_vars, sm_mi,
	    colloc_key, cov_t1_coeffs, cov_t2_coeffs, lev-1, ref_key);
	  for (v=0; v<num_v; ++v)
	    cov_t2_coeffs_lsp[v] = r1_val_mm1 * r2_grad[v]
	      + r2_val_mm2 * r1_grad[v] - prev_grad[v];
	}
      }
    }
  }
  else {
    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = (partial) ? ref_key[lev][1] : colloc_key[lev].size();
      cov_t1_coeffs[lev].resize(num_sets); cov_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	num_tp_pts = colloc_key[lev][set].size();
	RealVector& cov_t1_coeffs_ls = cov_t1_coeffs[lev][set];
	cov_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	// type1 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	for (pt=0; pt<num_tp_pts; ++pt) {
	  RealVector c_vars(Teuchos::View,
	    const_cast<Real*>(var_sets[lev][set][pt]), (int)num_v);
	  r1_val_mm1 = value(c_vars, sm_mi, colloc_key, r1_t1_coeffs,
	    r1_t2_coeffs, lev, ref_key) - mean_r1;
	  r2_val_mm2 = (same) ? r1_val_mm1 : value(c_vars, sm_mi, colloc_key,
	    r2_t1_coeffs, r2_t2_coeffs, lev, ref_key) - mean_r2;
	  cov_t1_coeffs_ls[pt] = r1_val_mm1 * r2_val_mm2 - value(c_vars, sm_mi,
	    colloc_key, cov_t1_coeffs, cov_t2_coeffs, lev-1, ref_key);
	}
      }
    }
  }
}


void HierarchInterpPolyApproximation::
central_product_gradient_interpolant(
  HierarchInterpPolyApproximation* hip_approx_2, Real mean_r1, Real mean_r2,
  const RealVector& mean1_grad, const RealVector& mean2_grad, 
  RealMatrix2DArray& cov_t1_coeff_grads, const UShort2DArray& ref_key)
{
  // form hierarchical t1 coeff grads for (R_1 - \mu_1) (R_2 - \mu_2)
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray&           sm_mi = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      colloc_key = hsg_driver->collocation_key();
  const Sizet3DArray&     colloc_index = hsg_driver->collocation_indices();
  size_t lev, set, pt, num_levels = colloc_key.size(), num_sets, num_tp_pts,
    cntr = 0, c_index, v,
    num_deriv_vars = expT1CoeffGradsIter->second[0][0].numRows();
  bool partial = !ref_key.empty();
  const SDVArray& sdv_array = surrData.variables_data();
  const SDRArray& sdr_array_1 = surrData.response_data();
  const SDRArray& sdr_array_2 = hip_approx_2->surrData.response_data();

  // level 0 (assume this is always contained in a partial ref_key)
  cov_t1_coeff_grads.resize(num_levels); cov_t1_coeff_grads[0].resize(1);
  cov_t1_coeff_grads[0][0].shapeUninitialized(num_deriv_vars, 1);
  Real* cov_t1_coeff_grads_000 = cov_t1_coeff_grads[0][0][0];
  c_index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
  const SurrogateDataResp& sdr0_1 = sdr_array_1[c_index];
  const SurrogateDataResp& sdr0_2 = sdr_array_2[c_index];
  Real r1_mm = sdr0_1.response_function() - mean_r1,
       r2_mm = sdr0_2.response_function() - mean_r2;
  const RealVector& r1_grad = sdr0_1.response_gradient();
  const RealVector& r2_grad = sdr0_2.response_gradient();
  for (v=0; v<num_deriv_vars; ++v)
    cov_t1_coeff_grads_000[v] = r1_mm * (r2_grad[v] - mean2_grad[v])
                              + r2_mm * (r1_grad[v] - mean1_grad[v]);
  ++cntr;
  // levels 1:w
  for (lev=1; lev<num_levels; ++lev) {
    num_sets = (partial) ? ref_key[lev][1] : colloc_key[lev].size();
    cov_t1_coeff_grads[lev].resize(num_sets);
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = colloc_key[lev][set].size();
      RealMatrix& cov_t1_coeff_grads_ls = cov_t1_coeff_grads[lev][set];
      cov_t1_coeff_grads_ls.shapeUninitialized(num_deriv_vars, num_tp_pts);
      // type1 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
      for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	c_index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	const SurrogateDataResp& sdr_1 = sdr_array_1[c_index];
	const SurrogateDataResp& sdr_2 = sdr_array_2[c_index];
	Real r1_mm = sdr_1.response_function() - mean_r1,
             r2_mm = sdr_2.response_function() - mean_r2;
	const RealVector& r1_grad = sdr_1.response_gradient();
	const RealVector& r2_grad = sdr_2.response_gradient();
	const RealVector& prev_grad = gradient_nonbasis_variables(
	  sdv_array[c_index].continuous_variables(), sm_mi, colloc_key,
	  cov_t1_coeff_grads, lev-1, ref_key);
	Real* cov_t1_coeff_grads_lsp = cov_t1_coeff_grads_ls[pt];
	for (v=0; v<num_deriv_vars; ++v)
	  cov_t1_coeff_grads_lsp[v] = r1_mm * (r2_grad[v] - mean2_grad[v])
	    + r2_mm * (r1_grad[v] - mean1_grad[v]) - prev_grad[v];
      }
    }
  }
}


void HierarchInterpPolyApproximation::
central_product_gradient_interpolant(const RealMatrix2DArray& var_sets,
				     const UShort3DArray&     sm_mi,
				     const UShort4DArray&     colloc_key,
				     const RealVector2DArray& r1_t1_coeffs,
				     const RealMatrix2DArray& r1_t2_coeffs,
				     const RealMatrix2DArray& r1_t1_coeff_grads,
				     const RealVector2DArray& r2_t1_coeffs,
				     const RealMatrix2DArray& r2_t2_coeffs,
				     const RealMatrix2DArray& r2_t1_coeff_grads,
				     bool same, Real mean_r1, Real mean_r2,
				     const RealVector& mean1_grad,
				     const RealVector& mean2_grad, 
				     RealMatrix2DArray& cov_t1_coeff_grads,
				     const UShort2DArray& ref_key)
{
  // form hierarchical t1 coeff grads for (R_1 - \mu_1) (R_2 - \mu_2)
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  size_t lev, set, pt, num_levels = colloc_key.size(), num_sets, num_tp_pts,
    v, num_deriv_vars = expT1CoeffGradsIter->second[0][0].numRows();
  bool partial = !ref_key.empty();

  // level 0 (assume this is always contained in a partial ref_key)
  cov_t1_coeff_grads.resize(num_levels); cov_t1_coeff_grads[0].resize(1);
  cov_t1_coeff_grads[0][0].shapeUninitialized(num_deriv_vars, 1);
  Real* cov_t1_coeff_grads_000 = cov_t1_coeff_grads[0][0][0];
  RealVector c_vars(Teuchos::View, const_cast<Real*>(var_sets[0][0][0]),
		    (int)num_deriv_vars);
  Real r1_mm = value(c_vars, sm_mi, colloc_key, r1_t1_coeffs, r1_t2_coeffs, 0,
	 ref_key) - mean_r1,
       r2_mm = (same) ? r1_mm : value(c_vars, sm_mi, colloc_key, r2_t1_coeffs,
	 r2_t2_coeffs, 0, ref_key) - mean_r2;
  const RealVector& r1_grad = gradient_nonbasis_variables(c_vars, sm_mi,
    colloc_key, r1_t1_coeff_grads, 0, ref_key);
  const RealVector& r2_grad = (same) ? r1_grad : gradient_nonbasis_variables(
    c_vars, sm_mi, colloc_key, r2_t1_coeff_grads, 0, ref_key);
  for (v=0; v<num_deriv_vars; ++v)
    cov_t1_coeff_grads_000[v] = r1_mm * (r2_grad[v] - mean2_grad[v])
                              + r2_mm * (r1_grad[v] - mean1_grad[v]);
  // levels 1:w
  for (lev=1; lev<num_levels; ++lev) {
    num_sets = (partial) ? ref_key[lev][1] : colloc_key[lev].size();
    cov_t1_coeff_grads[lev].resize(num_sets);
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = colloc_key[lev][set].size();
      RealMatrix& cov_t1_coeff_grads_ls = cov_t1_coeff_grads[lev][set];
      cov_t1_coeff_grads_ls.shapeUninitialized(num_deriv_vars, num_tp_pts);
      // type1 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
      for (pt=0; pt<num_tp_pts; ++pt) {
	RealVector c_vars(Teuchos::View,
	  const_cast<Real*>(var_sets[lev][set][pt]), (int)num_deriv_vars);
	Real r1_mm = value(c_vars, sm_mi, colloc_key, r1_t1_coeffs,
	       r1_t2_coeffs, lev, ref_key) - mean_r1,
             r2_mm = (same) ? r1_mm : value(c_vars, sm_mi, colloc_key,
	       r2_t1_coeffs, r2_t2_coeffs, lev, ref_key) - mean_r2;
	const RealVector& r1_grad = gradient_nonbasis_variables(c_vars, sm_mi,
	  colloc_key, r1_t1_coeff_grads, lev, ref_key);
	const RealVector& r2_grad = gradient_nonbasis_variables(c_vars, sm_mi,
	  colloc_key, r2_t1_coeff_grads, lev, ref_key);
	const RealVector& prev_grad = gradient_nonbasis_variables(c_vars, sm_mi,
	  colloc_key, cov_t1_coeff_grads, lev-1, ref_key);
	Real* cov_t1_coeff_grads_lsp = cov_t1_coeff_grads_ls[pt];
	for (v=0; v<num_deriv_vars; ++v)
	  cov_t1_coeff_grads_lsp[v]
	    = r1_mm * (r2_grad[v] - mean2_grad[v])
	    + r2_mm * (r1_grad[v] - mean1_grad[v]) - prev_grad[v];
      }
    }
  }
}


void HierarchInterpPolyApproximation::
integrate_response_moments(size_t num_moments, const UShort3DArray& sm_mi,
			   const UShort4DArray& colloc_key,
			   const Sizet3DArray&  colloc_index,
			   const SDVArray& sdv_array, const SDRArray& sdr_array)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in HierarchInterpPoly"
	  << "Approximation::integrate_response_moments()" << std::endl;
    abort_handler(-1);
  }

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  size_t lev, set, pt, num_levels = colloc_key.size(), num_sets, num_tp_pts,
    cntr, c_index, v, num_v = sharedDataRep->numVars;
  int m_index, moment;
  bool use_derivs = data_rep->basisConfigOptions.useDerivs;

  if (numericalMoments.length() != num_moments)
    numericalMoments.sizeUninitialized(num_moments);
  Real& mean = numericalMoments[0];
  mean = expectation(expT1CoeffsIter->second, expT2CoeffsIter->second);

  // size moment coefficient arrays
  RealVector2DArray mom_t1_coeffs(num_levels);
  RealMatrix2DArray mom_t2_coeffs(num_levels);
  for (lev=0; lev<num_levels; ++lev) {
    num_sets = colloc_key[lev].size();
    mom_t1_coeffs[lev].resize(num_sets);
    mom_t2_coeffs[lev].resize(num_sets);
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = colloc_key[lev][set].size();
      mom_t1_coeffs[lev][set].sizeUninitialized(num_tp_pts);
      if (use_derivs)
	mom_t2_coeffs[lev][set].shapeUninitialized(num_v, num_tp_pts);
    }
  }

  for (m_index=1; m_index<num_moments; ++m_index) {
    moment = m_index+1; cntr = 0;
    if (use_derivs) {
      // level 0
      c_index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
      const SurrogateDataResp& sdr0 = sdr_array[c_index];
      Real data_fn_mm         = sdr0.response_function() - mean;
      mom_t1_coeffs[0][0][0]  = std::pow(data_fn_mm, moment);
      Real* mom_t2_coeffs_000 = mom_t2_coeffs[0][0][0];
      Real deriv = moment * std::pow(data_fn_mm, m_index);
      const RealVector& data_grad = sdr0.response_gradient();
      for (v=0; v<num_v; ++v)
	mom_t2_coeffs_000[v] = deriv * data_grad[v];
      ++cntr;
      // levels 1:w
      for (lev=1; lev<num_levels; ++lev) {
	num_sets = colloc_key[lev].size();
	for (set=0; set<num_sets; ++set) {
	  num_tp_pts = colloc_key[lev][set].size();
	  RealVector& mom_t1_coeffs_ls = mom_t1_coeffs[lev][set];
	  RealMatrix& mom_t2_coeffs_ls = mom_t2_coeffs[lev][set];
	  for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	    c_index = (colloc_index.empty()) ? cntr :
	      colloc_index[lev][set][pt];
	    const RealVector& c_vars
	      = sdv_array[c_index].continuous_variables();
	    const SurrogateDataResp& sdr = sdr_array[c_index];
	    // type1 hierarchical interpolation of (R - \mu)^moment
	    data_fn_mm = sdr.response_function() - mean;
	    mom_t1_coeffs_ls[pt] = std::pow(data_fn_mm, moment)
	      - value(c_vars, sm_mi, colloc_key, mom_t1_coeffs,
		      mom_t2_coeffs, lev-1);
	    // type2 hierarchical interpolation of (R - \mu)^moment
	    // --> interpolated grads are moment(R-\mu)^{moment-1} R'
	    Real* mom_t2_coeffs_lsp = mom_t2_coeffs_ls[pt];
	    deriv = moment * std::pow(data_fn_mm, m_index);
	    const RealVector& data_grad = sdr.response_gradient();
	    const RealVector& prev_grad = gradient_basis_variables(c_vars,
	      sm_mi, colloc_key, mom_t1_coeffs, mom_t2_coeffs, lev-1);
	    for (v=0; v<num_v; ++v)
	      mom_t2_coeffs_lsp[v] = deriv * data_grad[v] - prev_grad[v];
	  }
	}
      }
    }
    else {
      // level 0
      c_index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
      mom_t1_coeffs[0][0][0]
	= std::pow(sdr_array[c_index].response_function() - mean, moment);
      ++cntr;
      // levels 1:w
      for (lev=1; lev<num_levels; ++lev) {
	num_sets = colloc_key[lev].size();
	for (set=0; set<num_sets; ++set) {
	  RealVector& mom_t1_coeffs_ls = mom_t1_coeffs[lev][set];
	  num_tp_pts = colloc_key[lev][set].size();
	  // type1 hierarchical interpolation of (R - \mu)^moment
	  for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	    c_index = (colloc_index.empty()) ? cntr :
	      colloc_index[lev][set][pt];
	    mom_t1_coeffs_ls[pt]
	      = std::pow(sdr_array[c_index].response_function() - mean, moment)
	      - value(sdv_array[c_index].continuous_variables(), sm_mi,
		      colloc_key, mom_t1_coeffs, mom_t2_coeffs, lev-1);
	  }
	}
      }
    }

    // pass these exp coefficients into a general expectation fn
    numericalMoments[m_index] = expectation(mom_t1_coeffs, mom_t2_coeffs);
  }

  // standardize third and higher central moments, if present
  //standardize_moments(numericalMoments);

  /*
  if (numericalMoments.size() != num_moments)
    numericalMoments.size(num_moments);
  if (use_derivs)
    integrate_moments(expT1CoeffsIter->second, expT2CoeffsIter->second,
		      hsg_driver->type1_hierarchical_weight_sets(),
		      hsg_driver->type2_hierarchical_weight_sets(),
		      numericalMoments);
  else
    integrate_moments(expT1CoeffsIter->second,
                      hsg_driver->type1_hierarchical_weight_sets(),
		      numericalMoments);
  */
}


void HierarchInterpPolyApproximation::
integrate_response_moments(size_t num_moments,
			   const RealMatrix2DArray& var_sets,
			   const UShort3DArray&     sm_mi,
			   const UShort4DArray&     colloc_key,
			   const RealVector2DArray& t1_coeffs,
			   const RealMatrix2DArray& t2_coeffs,
			   const RealVector2DArray& t1_wts,
			   const RealMatrix2DArray& t2_wts)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in HierarchInterpPoly"
	  << "Approximation::integrate_response_moments()" << std::endl;
    abort_handler(-1);
  }

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  size_t lev, set, pt, num_levels = colloc_key.size(), num_sets, num_tp_pts,
    v, num_v = sharedDataRep->numVars;
  int m_index, moment;
  bool use_derivs = data_rep->basisConfigOptions.useDerivs;

  if (numericalMoments.length() != num_moments)
    numericalMoments.sizeUninitialized(num_moments);
  Real& mean = numericalMoments[0];
  mean = expectation(t1_coeffs, t2_coeffs, t1_wts, t2_wts);

  // size moment coefficient arrays
  RealVector2DArray mom_t1_coeffs(num_levels),  r_vals(num_levels);
  RealMatrix2DArray mom_t2_coeffs(num_levels), r_grads(num_levels);
  for (lev=0; lev<num_levels; ++lev) {
    num_sets = colloc_key[lev].size();
    mom_t1_coeffs[lev].resize(num_sets);  r_vals[lev].resize(num_sets);
    mom_t2_coeffs[lev].resize(num_sets);  r_grads[lev].resize(num_sets);
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = colloc_key[lev][set].size();
      mom_t1_coeffs[lev][set].sizeUninitialized(num_tp_pts);
      r_vals[lev][set].sizeUninitialized(num_tp_pts);
      if (use_derivs) {
	mom_t2_coeffs[lev][set].shapeUninitialized(num_v, num_tp_pts);
	r_grads[lev][set].shapeUninitialized(num_v, num_tp_pts);
      }
    }
  }

  // precompute all value() and gradient_basis_variables() for {t1,t2}_coeffs
  RealVector c_vars(Teuchos::View, const_cast<Real*>(var_sets[0][0][0]),
		    (int)num_v);
  r_vals[0][0][0] =
    value(c_vars, sm_mi, colloc_key, t1_coeffs, t2_coeffs, 0) - mean;
  if (use_derivs)
    Teuchos::setCol(gradient_basis_variables(c_vars, sm_mi, colloc_key,
      t1_coeffs, t2_coeffs, 0), 0, r_grads[0][0]);
  for (lev=1; lev<num_levels; ++lev) {
    num_sets = colloc_key[lev].size();
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = colloc_key[lev][set].size();
      RealVector&         r_vals_ls =   r_vals[lev][set];
      RealMatrix&        r_grads_ls =  r_grads[lev][set];
      const RealMatrix& var_sets_ls = var_sets[lev][set];
      for (pt=0; pt<num_tp_pts; ++pt) {
	RealVector
	  c_vars(Teuchos::View, const_cast<Real*>(var_sets_ls[pt]), (int)num_v);
	r_vals_ls[pt] =
	  value(c_vars, sm_mi, colloc_key, t1_coeffs, t2_coeffs, lev) - mean;
	if (use_derivs)
	  Teuchos::setCol(gradient_basis_variables(c_vars, sm_mi, colloc_key,
	    t1_coeffs, t2_coeffs, lev), (int)pt, r_grads_ls);
      }
    }
  }

  for (m_index=1; m_index<num_moments; ++m_index) {
    moment = m_index+1;

    // level 0
    RealVector c_vars(Teuchos::View, const_cast<Real*>(var_sets[0][0][0]),
		      (int)num_v);
    Real r_val_mm = r_vals[0][0][0];
    mom_t1_coeffs[0][0][0] = std::pow(r_val_mm, moment);
    if (use_derivs) {
      // level 0 continued
      Real *mom_t2_coeffs_000 = mom_t2_coeffs[0][0][0],
	   *r_grad = r_grads[0][0][0],
	   deriv = moment * std::pow(r_val_mm, m_index);
      for (v=0; v<num_v; ++v)
	mom_t2_coeffs_000[v] = deriv * r_grad[v];
    
      // levels 1:w
      for (lev=1; lev<num_levels; ++lev) {
	num_sets = colloc_key[lev].size();
	for (set=0; set<num_sets; ++set) {
	  num_tp_pts = colloc_key[lev][set].size();
	  RealVector&  mom_t1_coeffs_ls = mom_t1_coeffs[lev][set];
	  RealMatrix&  mom_t2_coeffs_ls = mom_t2_coeffs[lev][set];
	  RealVector&         r_vals_ls =        r_vals[lev][set];
	  RealMatrix&        r_grads_ls =       r_grads[lev][set];
	  const RealMatrix& var_sets_ls =      var_sets[lev][set];
	  for (pt=0; pt<num_tp_pts; ++pt) {
	    RealVector c_vars(Teuchos::View,
	      const_cast<Real*>(var_sets_ls[pt]), (int)num_v);
	    // type1 hierarchical interpolation of (R - \mu)^moment
	    r_val_mm = r_vals_ls[pt];
	    mom_t1_coeffs_ls[pt] = std::pow(r_val_mm, moment) - value(c_vars,
	      sm_mi, colloc_key, mom_t1_coeffs, mom_t2_coeffs, lev-1);
	    // type2 hierarchical interpolation of (R - \mu)^moment
	    // --> interpolated grads are moment(R-\mu)^{moment-1} R'
	    Real *mom_t2_coeffs_lsp = mom_t2_coeffs_ls[pt],
	         *r_grad = r_grads_ls[pt],
	         deriv = moment * std::pow(r_val_mm, m_index);
	    const RealVector& prev_grad = gradient_basis_variables(c_vars,
	      sm_mi, colloc_key, mom_t1_coeffs, mom_t2_coeffs, lev-1);
	    for (v=0; v<num_v; ++v)
	      mom_t2_coeffs_lsp[v] = deriv * r_grad[v] - prev_grad[v];
	  }
	}
      }
    }
    else {
      // levels 1:w
      for (lev=1; lev<num_levels; ++lev) {
	num_sets = colloc_key[lev].size();
	for (set=0; set<num_sets; ++set) {
	  RealVector&  mom_t1_coeffs_ls = mom_t1_coeffs[lev][set];
	  RealVector&         r_vals_ls =        r_vals[lev][set];
	  const RealMatrix& var_sets_ls =      var_sets[lev][set];
	  num_tp_pts = colloc_key[lev][set].size();
	  // type1 hierarchical interpolation of (R - \mu)^moment
	  for (pt=0; pt<num_tp_pts; ++pt) {
	    RealVector c_vars(Teuchos::View,
	      const_cast<Real*>(var_sets_ls[pt]), (int)num_v);
	    mom_t1_coeffs_ls[pt] = std::pow(r_vals_ls[pt], moment) - value(
	      c_vars, sm_mi, colloc_key, mom_t1_coeffs, mom_t2_coeffs, lev-1);
	  }
	}
      }
    }

    // pass these exp coefficients into a general expectation fn
    numericalMoments[m_index]
      = expectation(mom_t1_coeffs, mom_t2_coeffs, t1_wts, t2_wts);
  }

  // standardize third and higher central moments, if present
  //standardize_moments(numericalMoments);
}


/** Computes the variance of component functions.  Assumes that all
    subsets of set_value have been computed in advance which will be
    true so long as the partial_variance is called following
    appropriate enumeration of set value. */
void HierarchInterpPolyApproximation::
compute_partial_variance(const BitArray& set_value)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  // Compute the partial integral corresponding to set_value
  Real& variance = partialVariance[data_rep->sobolIndexMap[set_value]];

  // Follows Tang, Iaccarino, Eldred (conference paper AIAA-2010-2922),
  // but adapted for hierarchical interpolation.

  // Perform inner integral over complementary set u' to form new weighted
  // coefficients h (stored as member_coeffs), stored hierarchically over
  // the member dimensions.
  RealVector2DArray member_t1_coeffs, member_t1_wts;
  RealMatrix2DArray member_t2_coeffs, member_t2_wts;
  UShort4DArray member_colloc_key; Sizet3DArray member_colloc_index;
  member_coefficients_weights(set_value, member_t1_coeffs, member_t1_wts,
			      member_t2_coeffs, member_t2_wts,
			      member_colloc_key, member_colloc_index);

  // re-interpolate h^2 (0 passed for mean to central product) using a
  // hierarchical formulation over the reduced member dimensions
  RealVector2DArray prod_member_t1_coeffs;
  RealMatrix2DArray prod_member_t2_coeffs;
  central_product_member_coefficients(set_value, member_t1_coeffs,
				      member_t2_coeffs, member_colloc_key,
				      member_colloc_index, 0.,
				      prod_member_t1_coeffs,
				      prod_member_t2_coeffs);

  // integrate h^2 over the reduced member dimensions
  variance = expectation(prod_member_t1_coeffs, prod_member_t2_coeffs,
			 member_t1_wts,         member_t2_wts);

#ifdef VBD_DEBUG
  PCout << "Partial variance = " << variance;
#endif // VBD_DEBUG
  // compute proper subsets and subtract their contributions
  InterpPolyApproximation::compute_partial_variance(set_value);
#ifdef VBD_DEBUG
  PCout << " (raw) " << variance << " (minus subsets)\n";
#endif // VBD_DEBUG
}


void HierarchInterpPolyApproximation::compute_total_sobol_indices()
{
  // Compute the total expansion mean and variance.  For standard mode, these
  // are likely already available, as managed by computed{Mean,Variance} data
  // reuse trackers.  For all vars mode, they are computed without partial
  // integration restricted to the random indices.  If negligible variance
  // (deterministic test fn) or negative variance (poor sparse grid resolution),
  // then attribution of this variance is suspect.  Note: zero is a good choice
  // since it drops out from anisotropic refinement based on a response-average
  // of total Sobol' indices.
  Real total_variance = variance();
  if (total_variance <= SMALL_NUMBER)
    { totalSobolIndices = 0.; return; }
  Real total_mean = mean();

  UShortArray quad_order; Real complement_variance;
  size_t v, num_v = sharedDataRep->numVars;
  BitArray complement_set(num_v);
  RealVector2DArray member_t1_coeffs, member_t1_wts, cprod_member_t1_coeffs;
  RealMatrix2DArray member_t2_coeffs, member_t2_wts, cprod_member_t2_coeffs;
  UShort4DArray member_colloc_key; Sizet3DArray member_colloc_index;
  // iterate each variable 
  for (v=0; v<num_v; ++v) {
    // define complement_set that includes all but index of interest
    complement_set.set(); complement_set.flip(v);

    // Perform inner integral over complementary set u' to form new weighted
    // coefficients h (stored as member_coeffs), stored hierarchically over
    // the member dimensions.
    member_coefficients_weights(complement_set, member_t1_coeffs, member_t1_wts,
				member_t2_coeffs, member_t2_wts,
				member_colloc_key, member_colloc_index);

    // re-interpolate (h-\mu)^2 using a hierarchical formulation over the
    // reduced member dimensions
    central_product_member_coefficients(complement_set, member_t1_coeffs,
					member_t2_coeffs, member_colloc_key,
					member_colloc_index, total_mean,
					cprod_member_t1_coeffs,
					cprod_member_t2_coeffs);

    // integrate (h-\mu)^2 over the reduced member dimensions using member wts
    complement_variance = expectation(cprod_member_t1_coeffs,
				      cprod_member_t2_coeffs, member_t1_wts,
				      member_t2_wts);

    // define total Sobol' index for this var from complement & total variance
    totalSobolIndices[v] = 1. - complement_variance / total_variance;
  }
}


/** Forms a lower dimensional interpolant over variables that are
    members of the given set. */
void HierarchInterpPolyApproximation::
member_coefficients_weights(const BitArray& member_bits,
  RealVector2DArray& member_t1_coeffs,  RealVector2DArray& member_t1_wts,
  RealMatrix2DArray& member_t2_coeffs,  RealMatrix2DArray& member_t2_wts,
  UShort4DArray&     member_colloc_key, Sizet3DArray&      member_colloc_index)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver   = data_rep->hsg_driver();
  const UShort3DArray&      sm_mi        = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      colloc_key   = hsg_driver->collocation_key();
  const Sizet3DArray&       colloc_index = hsg_driver->collocation_indices();

  const RealVector2DArray& exp_t1_coeffs = expT1CoeffsIter->second;
  const RealMatrix2DArray& exp_t2_coeffs = expT2CoeffsIter->second;

  // Perform inner integral over complementary set u' (non-member vars) to
  // form new weighted expansion coefficients h (stored as member_t1_coeffs)
  size_t v, num_v = sharedDataRep->numVars, lev, set, pt,
    num_levels = exp_t1_coeffs.size(), num_sets, num_tp_pts,
    num_member_coeffs, v_cntr, p_cntr = 0, m_index;
  Real member_wt, nonmember_wt;
  SizetArray indexing_factor; // for indexing member coeffs,wts
  UShortArray delta_sizes; UShort2DArray delta_keys;
  member_t1_coeffs.resize(num_levels);  member_t1_wts.resize(num_levels);
  member_t2_coeffs.resize(num_levels);  member_t2_wts.resize(num_levels);
  member_colloc_key.resize(num_levels); member_colloc_index.resize(num_levels);
  for (lev=0; lev<num_levels; ++lev) {
    const RealVectorArray& t1_coeffs_l = exp_t1_coeffs[lev];
    num_sets = t1_coeffs_l.size();
    member_t1_coeffs[lev].resize(num_sets); member_t1_wts[lev].resize(num_sets);
    member_t2_coeffs[lev].resize(num_sets); member_t2_wts[lev].resize(num_sets);
    member_colloc_key[lev].resize(num_sets);
    member_colloc_index[lev].resize(num_sets);
    for (set=0; set<num_sets; ++set) {
      const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
      const RealMatrix& t2_coeffs_ls = exp_t2_coeffs[lev][set];
      const UShortArray&    sm_mi_ls = sm_mi[lev][set];
      RealVector&     m_t1_coeffs_ls = member_t1_coeffs[lev][set];
      RealVector&        m_t1_wts_ls = member_t1_wts[lev][set];
      RealMatrix&     m_t2_coeffs_ls = member_t2_coeffs[lev][set];
      RealMatrix&        m_t2_wts_ls = member_t2_wts[lev][set];
      UShort2DArray&        m_key_ls = member_colloc_key[lev][set];
      SizetArray&         m_index_ls = member_colloc_index[lev][set];

      // precompute sizes and indexing factors and init member arrays
      num_tp_pts = t1_coeffs_ls.length();
      if (num_tp_pts) { // can be zero in case of growth restriction
	hsg_driver->levels_to_delta_keys(sm_mi_ls, delta_keys);
	hsg_driver->levels_to_delta_sizes(sm_mi_ls, delta_sizes);
	indexing_factor.clear();
	for (v=0, num_member_coeffs=1; v<num_v; ++v)
	  if (member_bits[v]) {
	    indexing_factor.push_back(num_member_coeffs); // for m_index below
	    num_member_coeffs *= delta_sizes[v];
	  }
	m_t1_coeffs_ls.size(num_member_coeffs);
	m_t1_wts_ls.size(num_member_coeffs);
	if (data_rep->basisConfigOptions.useDerivs) {
	  m_t2_coeffs_ls.shape(num_v, num_member_coeffs);
	  m_t2_wts_ls.shape(num_v, num_member_coeffs);
	}
	m_key_ls.resize(num_member_coeffs);
	m_index_ls.resize(num_member_coeffs);
      }

      for (pt=0; pt<num_tp_pts; ++pt, ++p_cntr) {
	// convert key_lsp to a corresponding index on member_t1_{coeffs,wts}.
	// We must also define a mapping to a member-variable collocation
	// key/index (from collapsing non-members) for use downstream.
	const UShortArray& key_lsp = colloc_key[lev][set][pt];
	for (v=0, m_index=0, v_cntr=0; v<num_v; ++v)
	  if (member_bits[v]) // key_lsp spans all pts, not deltas
	    m_index += find_index(delta_keys[v], key_lsp[v])
	            *  indexing_factor[v_cntr++];

	// integrate out nonmember dimensions and aggregate with type1 coeffs
	data_rep->type1_weight(key_lsp, sm_mi_ls, member_bits,
			       member_wt, nonmember_wt);
	m_t1_coeffs_ls[m_index] += nonmember_wt * t1_coeffs_ls[pt];
	// reduced dimension data updated more times than necessary, but
	// tracking this redundancy would be more expensive/complicated.
	// Note: key and index->c_vars components are the same only for the
	// member dimensions later used in value()/gradient_basis_variables().
	m_t1_wts_ls[m_index] = member_wt;
	m_index_ls[m_index]  = (colloc_index.empty()) ? p_cntr :
	  colloc_index[lev][set][pt];   // links back to surrData c_vars
	m_key_ls[m_index]    = key_lsp; // links back to interp polynomials

	// now do the same for the type2 coeffs and weights
	if (data_rep->basisConfigOptions.useDerivs) {
	  Real     *m_t2_coeffs_lsp = m_t2_coeffs_ls[m_index],
	              *m_t2_wts_lsp = m_t2_wts_ls[m_index];
	  const Real *t2_coeffs_lsp = t2_coeffs_ls[pt];
	  for (v=0; v<num_v; ++v) {
	    data_rep->type2_weight(v, key_lsp, sm_mi_ls, member_bits,
				   member_wt, nonmember_wt);
	    m_t2_coeffs_lsp[v] += nonmember_wt * t2_coeffs_lsp[v];
	    m_t2_wts_lsp[v]    =  member_wt;
	  }
	}
      }
#ifdef VBD_DEBUG
      PCout << "member_bits: " << member_bits // MSB->LSB: order reversed
	    << "\nmember_t1_coeffs[" << lev << "][" << set << "]:\n";
      write_data(PCout, member_t1_coeffs[lev][set]);
      PCout << "member_t1_wts[" << lev << "][" << set << "]:\n";
      write_data(PCout, member_t1_wts[lev][set]);
      if (data_rep->basisConfigOptions.useDerivs) {
	PCout << "member_t2_coeffs[" << lev << "][" << set << "]:\n";
	write_data(PCout, member_t2_coeffs[lev][set], false, true, true);
	PCout << "member_t2_wts[" << lev << "][" << set << "]:\n";
	write_data(PCout, member_t2_wts[lev][set],    false, true, true);
      }
      PCout << std::endl;
#endif // VBD_DEBUG
    }
  }
}


void HierarchInterpPolyApproximation::
central_product_member_coefficients(const BitArray& m_bits,
  const RealVector2DArray& m_t1_coeffs, const RealMatrix2DArray& m_t2_coeffs,
  const UShort4DArray& m_colloc_key,    const Sizet3DArray& m_colloc_index,
  Real mean, RealVector2DArray& cprod_m_t1_coeffs,
  RealMatrix2DArray& cprod_m_t2_coeffs)
{
  // We now have a lower dimensional hierarchical interpolant, for which
  // (h - mean)^2 must be formed hierarchically.  We use value() to both
  // reconstruct h from its surpluses and to evaluate new surpluses for
  // (h - mean)^2.  Note that we can simply pass mean = 0 for an h^2
  // product interpolant (used by compute_partial_variance()).

  // while colloc_{key,index} are redefined for member vars, sm_mi is not
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  const UShort3DArray& sm_mi = data_rep->hsg_driver()->smolyak_multi_index();
  const SDVArray& sdv_array = surrData.variables_data();

  size_t v, num_v = sharedDataRep->numVars, lev, set, pt,
    num_levels = m_t1_coeffs.size(), num_sets, num_tp_pts, c_index,
    h_level = num_levels - 1;
  SizetList member_indices;
  for (v=0; v<num_v; ++v)
    if (m_bits[v])
      member_indices.push_back(v);

  // form hierarchical t1/t2 coeffs for (h - mean)^2
  cprod_m_t1_coeffs.resize(num_levels); cprod_m_t1_coeffs[0].resize(1);
  cprod_m_t2_coeffs.resize(num_levels); cprod_m_t2_coeffs[0].resize(1);
  // level 0 type1
  cprod_m_t1_coeffs[0][0].sizeUninitialized(1);
  c_index = m_colloc_index[0][0][0];
  const RealVector& c_vars0 = sdv_array[c_index].continuous_variables();
  Real h_val_mm = value(c_vars0, sm_mi, m_colloc_key, m_t1_coeffs, m_t2_coeffs,
			h_level, member_indices) - mean;
  cprod_m_t1_coeffs[0][0][0] = h_val_mm * h_val_mm;
  if (data_rep->basisConfigOptions.useDerivs) {
    // level 0 type2
    const RealVector& h_grad_000 = gradient_basis_variables(c_vars0, sm_mi,
      m_colloc_key, m_t1_coeffs, m_t2_coeffs, h_level, member_indices);
    cprod_m_t2_coeffs[0][0].shapeUninitialized(num_v, 1);
    Real *cprod_m_t2_coeffs_000 = cprod_m_t2_coeffs[0][0][0];
    for (v=0; v<num_v; ++v)
      cprod_m_t2_coeffs_000[v] = 2. * h_val_mm * h_grad_000[v];
    // levels 1:w type1 and type2
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = m_colloc_key[lev].size();
      cprod_m_t1_coeffs[lev].resize(num_sets);
      cprod_m_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	RealVector& cprod_m_t1_coeffs_ls = cprod_m_t1_coeffs[lev][set];
	RealMatrix& cprod_m_t2_coeffs_ls = cprod_m_t2_coeffs[lev][set];
	num_tp_pts = m_colloc_key[lev][set].size();
	cprod_m_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	cprod_m_t2_coeffs_ls.shapeUninitialized(num_v,num_tp_pts);
	for (pt=0; pt<num_tp_pts; ++pt) {
	  c_index = m_colloc_index[lev][set][pt];
	  const RealVector& c_vars = sdv_array[c_index].continuous_variables();
	  // type1 hierarchical interpolation of h^2
	  h_val_mm = value(c_vars, sm_mi, m_colloc_key, m_t1_coeffs,
			   m_t2_coeffs, h_level, member_indices) - mean;
	  cprod_m_t1_coeffs_ls[pt] = h_val_mm * h_val_mm -
	    value(c_vars, sm_mi, m_colloc_key, cprod_m_t1_coeffs,
		  cprod_m_t2_coeffs, lev-1, member_indices);
	  // type2 hierarchical interpolation of h^2
	  // --> interpolated grads are 2 h h'
	  // --> make a copy of h_grad since prev_grad reuses approxGradient
	  Real* cprod_m_t2_coeffs_lsp = cprod_m_t2_coeffs_ls[pt];
	  RealVector h_grad = gradient_basis_variables(c_vars, sm_mi,
	    m_colloc_key, m_t1_coeffs, m_t2_coeffs, h_level, member_indices);
	  const RealVector& prev_grad =
	    gradient_basis_variables(c_vars, sm_mi, m_colloc_key,
	    cprod_m_t1_coeffs, cprod_m_t2_coeffs, lev-1, member_indices);
	  for (v=0; v<num_v; ++v)
	    cprod_m_t2_coeffs_lsp[v] = 2. * h_val_mm * h_grad[v] - prev_grad[v];
	}
#ifdef VBD_DEBUG
	PCout << "cprod_m_t1_coeffs[" << lev << "][" << set << "]:\n";
	write_data(PCout, cprod_m_t1_coeffs[lev][set]);
	PCout << "cprod_m_t2_coeffs[" << lev << "][" << set << "]:\n";
	write_data(PCout, cprod_m_t2_coeffs[lev][set], false, true, true);
	PCout << std::endl;
#endif // VBD_DEBUG
      }
    }
  }
  else {
    // levels 1:w type1
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = m_colloc_key[lev].size();
      cprod_m_t1_coeffs[lev].resize(num_sets);
      cprod_m_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	num_tp_pts = m_colloc_key[lev][set].size();
	RealVector& cprod_m_t1_coeffs_ls = cprod_m_t1_coeffs[lev][set];
	cprod_m_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	// type1 hierarchical interpolation of (h - mean)^2
	for (pt=0; pt<num_tp_pts; ++pt) {
	  c_index = m_colloc_index[lev][set][pt];
	  const RealVector& c_vars = sdv_array[c_index].continuous_variables();
	  h_val_mm = value(c_vars, sm_mi, m_colloc_key, m_t1_coeffs,
			   m_t2_coeffs, h_level, member_indices) - mean;
	  cprod_m_t1_coeffs_ls[pt] = h_val_mm * h_val_mm -
	    value(c_vars, sm_mi, m_colloc_key, cprod_m_t1_coeffs,
		  cprod_m_t2_coeffs, lev-1, member_indices);
	}
#ifdef VBD_DEBUG
	PCout << "cprod_m_t1_coeffs[" << lev << "][" << set << "]:\n";
	write_data(PCout, cprod_m_t1_coeffs[lev][set]);
	PCout << std::endl;
#endif // VBD_DEBUG
      }
    }
  }
}

}

/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        SharedHierarchInterpPolyApproxData
//- Description:  Implementation code for InterpPolyApproximation class
//-               
//- Owner:        Mike Eldred

#include "SharedHierarchInterpPolyApproxData.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"
#include "pecos_stat_util.hpp"

//#define DEBUG
//#define VBD_DEBUG

namespace Pecos {


void SharedHierarchInterpPolyApproxData::allocate_component_sobol()
{
  // Allocate memory specific to output control
  if (expConfigOptions.vbdFlag) {
    if (expConfigOptions.vbdOrderLimit == 1) // main effects only
      allocate_main_sobol();
    else { // main + interaction effects

      // One approach is to leverage PCE equivalence.  The exact order of
      // the interpolation polynomial (e.g., for nested rules, local or
      // gradient-enhanced interpolants) is not critical for defining
      // interactions; the issue is more the presence of constant dimensions.
      // While the collocation key has a very different meaning from the PCE
      // multi-index, the presence of non-zero's still indicates multi-point
      // interpolation and the presence of dimension effects.
      sobolIndexMap.clear();
      HierarchSparseGridDriver* hsg_driver
	= (HierarchSparseGridDriver*)driverRep;
      const UShort4DArray& key = hsg_driver->collocation_key();
      size_t lev, num_lev = key.size(), set, num_sets;
      for (lev=0; lev<num_lev; ++lev) {
	const UShort3DArray& key_l = key[lev];
	num_sets = key_l.size();
	for (set=0; set<num_sets; ++set)
	  multi_index_to_sobol_index_map(key_l[set]);
      }
      assign_sobol_index_map_values();

      // another approach: interrogate polynomialBasis[].interpolation_size()
      // or the quadrature/sparse level indices, again focusing on the presence
      // of constant dimensions (size = 1, level = 0).  But given the need to
      // regenerate the effect combinations from this reduced order data, the
      // collocation key idea seems preferable since it's already available.
      //polynomial_basis_to_sobol_indices();
    }
  }
}


void SharedHierarchInterpPolyApproxData::increment_component_sobol()
{
  // Allocate memory specific to output control
  if (expConfigOptions.vbdFlag && expConfigOptions.vbdOrderLimit != 1) {

    // return existing sobolIndexMap values to interaction counts
    reset_sobol_index_map_values();

    HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
    const UShort4DArray&      key        = hsg_driver->collocation_key();
    switch (expConfigOptions.refinementControl) {
    case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED: { // generalized sparse grids
      size_t lev = l1_norm(hsg_driver->trial_set());
      multi_index_to_sobol_index_map(key[lev].back());
      break;
    }
    default: { // isotropic/anisotropic refinement
      const UShortArray& incr_sets = hsg_driver->increment_sets();
      size_t lev, num_lev = key.size(), set, start_set, num_sets;
      for (lev=0; lev<num_lev; ++lev) {
	const UShort3DArray& key_l = key[lev];
	start_set = incr_sets[lev]; num_sets = key_l.size();
	for (set=start_set; set<num_sets; ++set)
	  multi_index_to_sobol_index_map(key_l[set]);
      }
      break;
    }
    }

    // update aggregated sobolIndexMap to indices into sobolIndices array
    assign_sobol_index_map_values();
  }
}


void SharedHierarchInterpPolyApproxData::pre_push_data()
{
  // pushIndex for Hierarch corresponds to candidate for current trial level

  // Note: pushIndex just caches result, avoiding need to call
  //       data_rep->retrieval_index() for each QoI
  const UShortArray& tr_set = ((SparseGridDriver*)driverRep)->trial_set();
  //unsigned short tr_lev = l1_norm(tr_set);
  pushIndex = find_index(poppedLevMultiIndex[activeKey]/*[tr_lev]*/, tr_set);
  // same as retrieval_index(), but eliminates key lookup for trial_set()
}


void SharedHierarchInterpPolyApproxData::post_push_data()
{
  // leave polynomialBasis as is (a previous increment is being restored)

  switch (expConfigOptions.refinementControl) {
  case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED: { // generalized sparse grids
    UShortArrayDeque& popped_lev_mi = poppedLevMultiIndex[activeKey]/*[tr_lev]*/;
    popped_lev_mi.erase(popped_lev_mi.begin() + pushIndex);
    break;
  }
  default:

    // Interpolation basis and component sobol already in incremented state

    break;
  }
}


void SharedHierarchInterpPolyApproxData::pre_combine_data()
{
  // Don't mix additive approach for hierarchical interpolation with
  // multiplicative/combined across multilevel-multifidelity
  if (expConfigOptions.combineType == MULT_COMBINE) {
    PCerr << "Error: only additive combinations supported in SharedHierarch"
	  << "InterpPolyApproxData::pre_combine_data()." << std::endl;
    abort_handler(-1);
  }

  // combine all multiIndex keys into combinedMultiIndex{,Map}
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const std::map<UShortArray, UShort3DArray>& sm_mi_map
    = hsg_driver->smolyak_multi_index_map();
  size_t i, num_combine = sm_mi_map.size(), lev, num_lev, combine_sm_map_ref;
  combinedSmolyakMultiIndex.clear();
  combinedSmolyakMultiIndexMap.resize(num_combine);
  std::map<UShortArray, UShort3DArray>::const_iterator sm_cit;
  for (sm_cit=sm_mi_map.begin(), i=0; sm_cit!=sm_mi_map.end(); ++sm_cit, ++i) {
    const UShort3DArray& sm_mi = sm_cit->second;
    num_lev = sm_mi.size();
    if (num_lev > combinedSmolyakMultiIndex.size()) // don't contract
      combinedSmolyakMultiIndex.resize(num_lev);
    Sizet2DArray& comb_sm_map_i = combinedSmolyakMultiIndexMap[i];
    comb_sm_map_i.resize(num_lev);
    for (lev=0; lev<num_lev; ++lev)
      append_multi_index(sm_mi[lev], combinedSmolyakMultiIndex[lev],
			 comb_sm_map_i[lev], combine_sm_map_ref);
  }

  // Flag indicates unstructured set ordering (can't assume a no-op return if
  // level/set counts are consistent).  Recompute combinedCollocKey from
  // combinedSmolyakMultiIndex, rather than overlay collocKey.
  hsg_driver->
    assign_collocation_key(combinedSmolyakMultiIndex, combinedCollocKey, false);
  // Can't use collocation indices for combined expansions without adding
  // additional tracking of model indices (no thanks), so create combined
  // points and weights to support expectation() calls.
  hsg_driver->compute_points_weights(combinedSmolyakMultiIndex,
    combinedCollocKey, combinedVarSets, combinedT1WeightSets,
    combinedT2WeightSets);
}


void SharedHierarchInterpPolyApproxData::combined_to_active(bool clear_combined)
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;

  hsg_driver->smolyak_multi_index(combinedSmolyakMultiIndex);
  hsg_driver->collocation_key(combinedCollocKey);
  // combinedSmolyakMultiIndexMap not used by hsg_driver

  //hsg_driver->combine_weight_sets(combinedSmolyakMultiIndexMap);
  //hsg_driver->combined_to_active();
  // Replace active weights with combined weights.
  // Update type2 wts even if inactive, so that 2D array sizes are correct
  // Note: inactive weight sets to be removed by clear_inactive()
  hsg_driver->type1_hierarchical_weight_sets(combinedT1WeightSets);
  hsg_driver->type2_hierarchical_weight_sets(combinedT2WeightSets);

  if (clear_combined) {
    combinedCollocKey.clear();
    combinedSmolyakMultiIndex.clear();  combinedSmolyakMultiIndexMap.clear();
    //combinedVarSets.clear();
    combinedT1WeightSets.clear();       combinedT2WeightSets.clear();
  }
}


size_t SharedHierarchInterpPolyApproxData::
barycentric_exact_index(const UShortArray& basis_index)
{
  size_t j, pt_index = 0, prod = 1, edi_j; unsigned short bi_j;
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  for (j=0; j<numVars; ++j) {
    bi_j = basis_index[j];
    // Note: if bi_j == 0, then constant interp with 1 point: we can replace
    // this constant interpolation with the value at the 1 colloc index (ei=0)
    if (bi_j) {
      edi_j = polynomialBasis[bi_j][j].exact_delta_index();
      if (edi_j == _NPOS) // manage exactIndex match w/o exactDeltaIndex match
	{ pt_index = _NPOS; break; }
      else {
	pt_index += edi_j * prod;
	prod     *= hsg_driver->level_to_delta_size(j, bi_j);
      }
    }
  }
  return pt_index;
}


size_t SharedHierarchInterpPolyApproxData::
barycentric_exact_index(const UShortArray& basis_index,
			const SizetList& subset_indices)
{
  size_t j, pt_index = 0, prod = 1, edi_j; unsigned short bi_j;
  SizetList::const_iterator cit;
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    j = *cit; bi_j = basis_index[j];
    // Note: if bi_j == 0, then constant interp with 1 point: we can replace
    // this constant interpolation with the value at the 1 colloc index (ei=0)
    if (bi_j) {
      edi_j = polynomialBasis[bi_j][j].exact_delta_index();
      if (edi_j == _NPOS) // manage exactIndex match w/o exactDeltaIndex match
	{ pt_index = _NPOS; break; }
      else {
	pt_index += edi_j * prod;
	prod     *= hsg_driver->level_to_delta_size(j, bi_j);
      }
    }
  }
  return pt_index;
}

}

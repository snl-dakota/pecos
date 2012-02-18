/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 HierarchSparseGridDriver
//- Description: Implementation code for HierarchSparseGridDriver class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "HierarchSparseGridDriver.hpp"
#include "PolynomialApproximation.hpp"
#include "sandia_sgmga.hpp"
#include "DistributionParams.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: HierarchSparseGridDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


int HierarchSparseGridDriver::grid_size()
{
  if (updateGridSize) {
    numCollocPts = 0;
    size_t i, j, num_lev = collocKey.size();
    if (num_lev != ssgLevel+1) {
      assign_smolyak_multi_index();
      assign_collocation_key();
      num_lev = collocKey.size();
    }
    for (i=0; i<num_lev; ++i) {
      const UShort3DArray& key_i = collocKey[i];
      size_t num_sets = key_i.size();
      for (j=0; j<num_sets; ++j)
	numCollocPts += key_i[j].size(); // hierarchical point increments
    }
    updateGridSize = false;
  }
  return numCollocPts;
}


/*
size_t HierarchSparseGridDriver::smolyak_size()
{
  size_t num_sm_mi = 0;
  unsigned short i, num_lev = smolyakMultiIndex.size();
  for (i=0; i<num_lev; ++i)
    num_sm_mi += smolyakMultiIndex[i].size();
  return num_sm_mi;
}
*/


void HierarchSparseGridDriver::assign_smolyak_multi_index()
{
  // Populate smolyakMultiIndex based on {iso,aniso}tropic index set constraint:
  // w = q - N = dimension-independent level.  For isotropic,
  //   w + 1 <= |i| <= w + N for i starts at 1 (used for index set defn.)
  //   w - N + 1 <= |j| <= w for j = i - 1 starts at 0 (used for generation)
  // For anisotropic, a weighted linear index set constraint is used.

  if (smolyakMultiIndex.size() == ssgLevel+1)
    return;

  smolyakMultiIndex.resize(ssgLevel+1);
  if (dimIsotropic) // initialize Smolyak multi_index
    for (unsigned short lev=0; lev<=ssgLevel; ++lev)
      PolynomialApproximation::total_order_multi_index(lev, numVars,
						       smolyakMultiIndex[lev]);
  else { // utilize webbur::sandia_sgmga_vcn_ordered
    for (unsigned short lev=0; lev<=ssgLevel; ++lev)
      smolyakMultiIndex[lev].clear();
    // w*alpha_min-|alpha| < |alpha . j| <= w*alpha_min for 0-based j index sets
    // With scaling alpha_min = 1: w-|alpha| < |alpha . j| <= w.
    // In the isotropic case, reduces to w-N < |j| <= w, which is the same as
    // w-N+1 <= |j| <= w.
    IntArray x(numVars), x_max(numVars); //x_max = ssgLevel;
    UShortArray index_set(numVars);
    Real wt_sum = 0., q_max = ssgLevel;
    for (size_t i=0; i<numVars; ++i) {
      const Real& wt_i = anisoLevelWts[i];
      wt_sum += wt_i;
      // minimum nonzero weight is scaled to 1, so just catch special case of 0
      x_max[i] = (wt_i > 1.e-10) ? (int)std::ceil(q_max/wt_i) : 0;
    }
    Real q_min = ssgLevel - wt_sum;
#ifdef DEBUG
    PCout << "q_min = " << q_min << " q_max = " << q_max;
#endif // DEBUG

    bool more = false;
    Real *aniso_wts = anisoLevelWts.values();
    int  *x0 = &x[0], *xm0 = &x_max[0];
    webbur::sandia_sgmga_vcn_ordered(numVars, aniso_wts, xm0, x0,
				     q_min, q_max, &more);
    while (more) {
      for (size_t i=0; i<numVars; ++i)
	index_set[i] = (unsigned short)x[i];
      unsigned int norm = index_norm(index_set);
      smolyakMultiIndex[norm].push_back(index_set);
      webbur::sandia_sgmga_vcn_ordered(numVars, aniso_wts, xm0, x0,
				       q_min, q_max, &more);
    }
  }

#ifdef DEBUG
  for (unsigned short i=0; i<=ssgLevel; ++i) {
    size_t j, num_sets = smolyakMultiIndex[i].size();
    for (j=0; j<num_sets; ++j)
      PCout << "Smolyak multi_index[" << i << "]:\n" << smolyakMultiIndex[i][j]
	    << "\n\n";
#endif // DEBUG
}


void HierarchSparseGridDriver::assign_collocation_key()
{
  if (collocKey.size() == ssgLevel+1)
    return;

  collocKey.resize(ssgLevel+1);
  if (nestedGrid) {
    size_t set, num_sets;
    UShort2DArray delta_quad(numVars);
    for (unsigned short lev=0; lev<=ssgLevel; ++lev) {
      const UShort2DArray& sm_mi_l = smolyakMultiIndex[lev];
      UShort3DArray&         key_l = collocKey[lev];
      num_sets = sm_mi_l.size();
      key_l.resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	level_to_delta_order(sm_mi_l[set], delta_quad);
	PolynomialApproximation::hierarchical_tensor_product_multi_index(
	  delta_quad, key_l[set]);
      }
    }
  }
  //else
  //  SparseGridDriver::assign_collocation_key();
}


void HierarchSparseGridDriver::update_collocation_key()
{
  UShort2DArray delta_quad(numVars);
  level_to_delta_order(trialSet, delta_quad);

  size_t lev = index_norm(trialSet);
  if (lev >= collocKey.size())
    collocKey.resize(lev+1);

  UShort3DArray& key_l = collocKey[lev];
  size_t set = key_l.size();
  UShort2DArray key_ls; key_l.push_back(key_ls); // update in place
  PolynomialApproximation::hierarchical_tensor_product_multi_index(
    delta_quad, key_l[set]);
}


void HierarchSparseGridDriver::assign_collocation_indices()
{
  size_t i, j, k, cntr = 0, num_levels = collocKey.size(), num_sets, num_tp_pts;
  if (collocIndices.size() < num_levels)
    collocIndices.resize(num_levels);
  for (i=0; i<num_levels; ++i) {
    const UShort3DArray& key_i = collocKey[i];
    num_sets = key_i.size();
    Sizet2DArray& indices_i = collocIndices[i];
    indices_i.resize(num_sets);
    for (j=0; j<num_sets; ++j) {
      const UShort2DArray& key_ij = key_i[j];
      num_tp_pts = key_ij.size();
      SizetArray& indices_ij = indices_i[j];
      indices_ij.resize(num_tp_pts);
      for (k=0; k<num_tp_pts; ++k, ++cntr)
	indices_ij[k] = cntr; // simple sequential ordering for unique pts
    }
  }
}


void HierarchSparseGridDriver::update_collocation_indices()
{
  size_t cntr = numCollocPts, num_levels = collocKey.size();
  if (collocIndices.size() < num_levels)
    collocIndices.resize(num_levels);

  size_t lev = index_norm(trialSet), num_tp_pts = collocKey[lev].back().size();
  Sizet2DArray& indices_l = collocIndices[lev];
  SizetArray indices; indices_l.push_back(indices); // update in place

  SizetArray& trial_indices = indices_l.back();
  trial_indices.resize(num_tp_pts);
  for (size_t k=0; k<num_tp_pts; ++k, ++cntr)
    trial_indices[k] = cntr;
  numCollocPts += num_tp_pts;
}


void HierarchSparseGridDriver::
level_to_delta_order(const UShortArray& levels, UShort2DArray& delta_quad)
{
  size_t i, j, num_lev = levels.size(), num_delta;
  if (delta_quad.size() != num_lev)
    delta_quad.resize(num_lev);
  for (i=0; i<num_lev; ++i) {
    unsigned short lev_i = levels[i], ord_i, ord_lm1 = 0;
    level_to_order(i, lev_i, ord_i);
    if (lev_i > 0) level_to_order(i, lev_i-1, ord_lm1);
    num_delta = ord_i - ord_lm1;
    UShortArray& delta_quad_i = delta_quad[i];
    delta_quad_i.resize(num_delta);
    switch(collocRules[i]) {
    case GAUSS_PATTERSON: // open nested
      for (j=0; j<num_delta; ++j)
	delta_quad_i[j] = 2*j; // new ends + new interior: 0,2,4,6,8,...
      break;
    case NEWTON_COTES: case CLENSHAW_CURTIS: // closed nested
      switch (lev_i) {
      case 0: delta_quad_i[0] = 0;                      break; // center of 1 pt
      case 1: delta_quad_i[0] = 0; delta_quad_i[1] = 2; break; // ends of 3 pt
      default:
	for (j=0; j<num_delta; ++j)
	  delta_quad_i[j] = 2*j+1; // new interior: 1,3,5,7,9,...
	break;
      }
      break;
    case GENZ_KEISTER: // open nested table lookup
      switch (lev_i) {
      case 0: delta_quad_i[0] = 0;                      break; // center of 1 pt
      case 1: delta_quad_i[0] = 0; delta_quad_i[1] = 2; break; // ends of 3 pt
      case 2:
	delta_quad_i[0] = 0; delta_quad_i[1] = 1; delta_quad_i[2] = 3;
	delta_quad_i[3] = 5; delta_quad_i[4] = 7; delta_quad_i[5] = 8;
	break; // 9 pt rule reusing 2, 4 (center), 6
      case 3:
	delta_quad_i[0] =  0; delta_quad_i[1] =  1; delta_quad_i[2] =  3;
	delta_quad_i[3] =  5; delta_quad_i[4] =  7; delta_quad_i[5] = 11;
	delta_quad_i[6] = 13; delta_quad_i[7] = 15; delta_quad_i[8] = 17;
	delta_quad_i[9] = 18;
	break; // 19 pt rule reusing 
      case 4:
	delta_quad_i[0]  =  0; delta_quad_i[1]  =  1; delta_quad_i[2]  =  2;
	delta_quad_i[3]  =  4; delta_quad_i[4]  =  6; delta_quad_i[5]  =  8;
	delta_quad_i[6]  = 12; delta_quad_i[7]  = 16; delta_quad_i[8]  = 18;
	delta_quad_i[9]  = 22; delta_quad_i[10] = 26; delta_quad_i[11] = 28;
	delta_quad_i[12] = 30; delta_quad_i[13] = 32; delta_quad_i[14] = 33;
	delta_quad_i[15] = 34;
	break; // 35 pt rule reusing
               // 3,5,7,9,10,11,13,14,15,17,19,20,21,23,24,25,27,29,31
      //case 5:  // 43 pt rule augments 19 pt rule, not 35 pt rule
      //  break; // disallow for hierarchical interpolation
      default:
	PCerr << "Error: out of range for hierarchical Genz-Keister rules in "
	      << "HierarchSparseGridDriver::level_to_delta_order()"<< std::endl;
	abort_handler(-1);
	break;
      }
      break;
    default:
      PCerr << "Error: bad rule type in level_to_delta_order()" << std::endl;
      abort_handler(-1);
      break;
    }
  }
}


void HierarchSparseGridDriver::compute_grid(RealMatrix& var_sets)
{
  assign_smolyak_multi_index();           // compute smolyakMultiIndex
  assign_collocation_key();               // compute collocKey
  assign_1d_collocation_points_weights(); // define 1-D point/weight sets

  // For efficiency reasons, incremental sparse grid definition uses different
  // point orderings than sgmg/sgmga.  Therefore, the reference grid
  // computations are kept completely separate.

  if (nestedGrid) {
    compute_points_weights(var_sets, type1WeightSets, type2WeightSets);
    //assign_collocation_indices();
  }
  /*
  else {
    // TO DO: hierarchical interpolation must difference interpolants among
    // full point sets rather than evaluating surpluses at point increments
    reference_unique(var_sets); // define reference grid

#ifdef DEBUG
    PCout << "HierarchSparseGridDriver::compute_grid() results:\n"
	  << "uniqueIndexMapping:\n" << uniqueIndexMapping << "\nvar_sets:\n";
    write_data(PCout, var_sets, false, true, true);
    if (trackUniqueProdWeights) {
      PCout << "\ntype1WeightSets:\n"; write_data(PCout, type1WeightSets);
      if (computeType2Weights) {
	PCout << "\ntype2WeightSets:\n";
	write_data(PCout, type2WeightSets, false, true, true);
      }
    }
#endif
  }
  */
}


void HierarchSparseGridDriver::
compute_points_weights(RealMatrix& pts, RealVector& t1_wts, RealMatrix& t2_wts,
		       const UShortArray& sm_index,
		       const UShort2DArray& colloc_key)
{
  size_t k, l, m, num_tp_pts = colloc_key.size();
  if (pts.numCols() != num_tp_pts)
    pts.shapeUninitialized(numVars, num_tp_pts);
  if (t1_wts.length() != num_tp_pts)
    t1_wts.sizeUninitialized(num_tp_pts);
  if (computeType2Weights && t2_wts.numCols() != num_tp_pts)
    t2_wts.shapeUninitialized(numVars, num_tp_pts);

  // update collocPts1D, type1CollocWts1D, and type2CollocWts1D
  UShortArray total_order;
  level_to_order(sm_index, total_order);
  update_1d_collocation_points_weights(total_order, sm_index);

  // define points and type 1/2 weights; weights are products of 1D weights
  for (k=0; k<num_tp_pts; ++k) {
    const UShortArray& key_k = colloc_key[k];
    Real* pt    =    pts[k]; // column vector
    Real& t1_wt = t1_wts[k]; t1_wt = 1.;
    for (l=0; l<numVars; ++l) {
      pt[l]  =      collocPts1D[sm_index[l]][l][key_k[l]];
      t1_wt *= type1CollocWts1D[sm_index[l]][l][key_k[l]];
    }
    if (computeType2Weights) {
      Real* t2_wt = t2_wts[k]; // column vector
      for (l=0; l<numVars; ++l) {
	Real& t2_wt_l = t2_wt[l]; t2_wt_l = 1.;
	for (m=0; m<numVars; ++m)
	  t2_wt_l *= (m==l) ? type2CollocWts1D[sm_index[m]][m][key_k[m]] :
	                      type1CollocWts1D[sm_index[m]][m][key_k[m]];
      }
    }
  }

#ifdef DEBUG
  PCout << "Tensor product points =\n"; write_data(PCout,pts,false,true,true);
  PCout << "Tensor product weights =\ntype1:\n"; write_data(PCout, t1_wts);
  PCout << "type2:\n"; write_data(PCout, t2_wts, false, true, true);
#endif // DEBUG
}


void HierarchSparseGridDriver::
compute_points_weights(RealMatrix& pts, RealVector& t1_wts, RealMatrix& t2_wts)
{
  size_t lev = index_norm(trialSet);
  compute_points_weights(pts, t1_wts, t2_wts, smolyakMultiIndex[lev].back(),
			 collocKey[lev].back());
}


/** Points are collapsed as required for compute_grid(var_sets), but t1/t2
    weights are hierarchical 2D arrays. */
void HierarchSparseGridDriver::
compute_points_weights(RealMatrix& pts, RealVector2DArray& t1_wts,
		       RealMatrix2DArray& t2_wts)
{
  size_t i, j, cntr = 0, num_colloc_pts = 0, num_tp_pts,
    num_levels = collocKey.size(), num_sets;
  if (t1_wts.size() != num_levels) t1_wts.resize(num_levels);
  if (t2_wts.size() != num_levels) t2_wts.resize(num_levels);
  // define num_colloc_pts
  for (i=0; i<num_levels; ++i) {
    const UShort3DArray& key_i = collocKey[i];
    num_sets = key_i.size();
    if (t1_wts[i].size() != num_sets) t1_wts[i].resize(num_sets);
    if (t2_wts[i].size() != num_sets) t2_wts[i].resize(num_sets);
    for (j=0; j<num_sets; ++j)
      num_colloc_pts += key_i[j].size();
  }
  if (pts.numCols() != num_colloc_pts)
    pts.shapeUninitialized(numVars, num_colloc_pts);

  // define points and type 1/2 weights; weights are products of 1D weights
  for (i=0; i<num_levels; ++i) {
    const UShort3DArray& key_i = collocKey[i];
    num_sets = key_i.size();
    for (j=0; j<num_sets; ++j) {
      const UShort2DArray& key_ij = key_i[j];
      num_tp_pts = key_ij.size();
      // take pts_ij sub-matrix view of full sample matrix pts
      RealMatrix pts_ij(Teuchos::View, pts, numVars, num_tp_pts, 0, cntr);
      compute_points_weights(pts_ij, t1_wts[i][j], t2_wts[i][j],
			     smolyakMultiIndex[i][j], key_ij);
      cntr += num_tp_pts;
    }
  }
}


void HierarchSparseGridDriver::initialize_sets()
{
  // define set O (old) from smolyakMultiIndex and smolyakCoeffs:
  //oldMultiIndex = smolyakMultiIndex;
  oldMultiIndex.clear();
  for (unsigned short lev=0; lev<=ssgLevel; ++lev)
    oldMultiIndex.insert(smolyakMultiIndex[lev].begin(),
			 smolyakMultiIndex[lev].end());

  // compute initial set A (active) by applying add_active_neighbors()
  // to the frontier of smolyakMultiIndex:
  if (dimIsotropic) {
    const UShort2DArray& sm_mi_l = smolyakMultiIndex[ssgLevel];
    size_t i, num_old_sets = sm_mi_l.size();
    for (i=0; i<num_old_sets; ++i)
      add_active_neighbors(sm_mi_l[i]);
  }
  else { // TO DO
    // For anisotropic, need to compute Pareto set.
  }

#ifdef DEBUG
  PCout << "SparseGridDriver::initialize_sets():\nold sets:\n" << oldMultiIndex
	<< "active sets:\n" << activeMultiIndex << std::endl;
#endif // DEBUG
}


void HierarchSparseGridDriver::push_trial_set(const UShortArray& set)
{
  trialSet = set;
  size_t lev = index_norm(set);
  if (smolyakMultiIndex.size() <= lev)
    smolyakMultiIndex.resize(lev+1);
  smolyakMultiIndex[lev].push_back(set);

  // collocKey, collocIndices, and uniqueIndexMapping updated within
  // either restore_set() or compute_trial_grid()
}


void HierarchSparseGridDriver::restore_set()
{
  // SparseGridDriver currently retains no memory, so updates are recomputed

  // update collocKey from trialSet
  update_collocation_key();
  if (nestedGrid) {
    //update_collocation_indices();
    // This approach stores less history than WeightSetsRef approach
    size_t lev = index_norm(trialSet);
    type1WeightSets[lev].push_back(savedT1WtSets[trialSet]);
    type2WeightSets[lev].push_back(savedT2WtSets[trialSet]);
  }
  /*
  else { // compute a2
    // update collocIndices and uniqueIndexMapping, but don't update pt/wt sets
    RealMatrix dummy_set;
    increment_unique(true, false, dummy_set);
    merge_unique(); // reset a1 --> INC3
  }
  */
}


void HierarchSparseGridDriver::compute_trial_grid(RealMatrix& var_sets)
{
  // track trial sets that have been evaluated (do here since
  // push_trial_set() used for both new trials and restorations)
  computedTrialSets.insert(trialSet);

  // update collocKey and compute trial variable/weight sets
  update_collocation_key();
  if (nestedGrid) {
    size_t lev = index_norm(trialSet);
    if (type1WeightSets.size() <= lev || type2WeightSets.size() <= lev)
      { type1WeightSets.resize(lev+1); type2WeightSets.resize(lev+1); }
    RealVectorArray& t1_wts_l = type1WeightSets[lev];
    RealMatrixArray& t2_wts_l = type2WeightSets[lev];
    size_t set = t1_wts_l.size();
    RealVector t1_wts_ls; t1_wts_l.push_back(t1_wts_ls); // update in place
    RealMatrix t2_wts_ls; t2_wts_l.push_back(t2_wts_ls); // update in place
    compute_points_weights(var_sets, t1_wts_l[set], t2_wts_l[set]);
    //update_collocation_indices();
  }
  /*
  else {
    compute_points_weights(a2Points, a2Type1Weights, a2Type2Weights);
    // update collocIndices, uniqueIndexMapping, and var_sets,
    // but don't recompute a2 data
    increment_unique(false, true, var_sets);
  }
  */

#ifdef DEBUG
  PCout << "compute_trial_grid() increment:\nunique variable sets:\n"
	<< var_sets;
#endif // DEBUG
}


void HierarchSparseGridDriver::pop_trial_set()
{
  size_t lev = index_norm(trialSet);
  if (nestedGrid)
    numCollocPts -= collocKey[lev].back().size(); // subtract # of trial pts
  /*
  else {
    numCollocPts -= numUnique2; // subtract number of trial points
    uniqueIndexMapping.resize(numCollocPts); // prune trial set from end
  }
  */

  smolyakMultiIndex[lev].pop_back();
  collocKey[lev].pop_back();
  //collocIndices[lev].pop_back();
  savedT1WtSets[trialSet] = type1WeightSets[lev].back();
  type1WeightSets[lev].pop_back();
  savedT2WtSets[trialSet] = type2WeightSets[lev].back();
  type2WeightSets[lev].pop_back();
}


void HierarchSparseGridDriver::merge_set()
{
  if (nestedGrid) {
    // no-op
  }
  //else
  //  merge_unique();
}


void HierarchSparseGridDriver::finalize_sets()
{
  // For final answer, push all evaluated sets into old and clear active.
  // Multiple trial insertion approach must be compatible with
  // Dakota::Approximation::savedSDPSet behavior (i.e., inc2/inc3 set 
  // insertions must occur one at a time without mixing).

  // don't insert activeMultiIndex, as this may include sets which have not
  // been evaluated (due to final update_sets() call); use computedTrialSets
  if (nestedGrid) {
    std::set<UShortArray>::iterator it;
    for (it=computedTrialSets.begin(); it!=computedTrialSets.end(); ++it) {
      trialSet = *it;
      size_t lev = index_norm(trialSet);
      smolyakMultiIndex[lev].push_back(trialSet);
      update_collocation_key();     // update collocKey
      //update_collocation_indices(); // update collocIndices and numCollocPts
      type1WeightSets[lev].push_back(savedT1WtSets[trialSet]);
      type2WeightSets[lev].push_back(savedT2WtSets[trialSet]);
    }
  }
  /*
  else {
    // ...
    // update a2 data, uniqueIndexMapping, collocIndices, numCollocPts
    finalize_unique(start_index);// assure no mixing of discrete a2's
    //merge_unique(); // a1 reference update not needed, no addtnl increments
    //update_reference();
  }
  */
  activeMultiIndex.clear(); computedTrialSets.clear();
  savedT1WtSets.clear();    savedT2WtSets.clear();
}


void HierarchSparseGridDriver::print_final_sets(bool converged_within_tol) const
{
  // this call should precede finalize_sets()
  size_t i, j, k, num_lev = smolyakMultiIndex.size();
  if (converged_within_tol) {
    size_t trial_lev = index_norm(trialSet);
    PCout << "Above tolerance index sets:\n";
    for (i=0; i<num_lev; ++i) {
      const UShort2DArray& sm_mi_i = smolyakMultiIndex[i];
      size_t num_sets = sm_mi_i.size();
      if (i==trial_lev) --num_sets; // omit trial set
      for (j=0; j<num_sets; ++j) {
	const UShortArray& sm_mi_ij = sm_mi_i[j];
	for (k=0; k<numVars; ++k)
	  PCout << std::setw(5) << sm_mi_ij[k];
	PCout << '\n';
      }
    }
    PCout << "Below tolerance index sets:\n";
    const UShort2DArray& sm_mi_l = smolyakMultiIndex[trial_lev];
    size_t last_set = sm_mi_l.size(); --last_set;
    const UShortArray& sm_mi_ls = sm_mi_l[last_set];
    for (k=0; k<numVars; ++k)
      PCout << std::setw(5) << sm_mi_ls[k];
    PCout << '\n';
  }
  else {
    PCout << "Final index sets:\n";
    for (i=0; i<num_lev; ++i) {
      const UShort2DArray& sm_mi_i = smolyakMultiIndex[i];
      size_t num_sets = sm_mi_i.size();
      for (j=0; j<num_sets; ++j) {
	const UShortArray& sm_mi_ij = sm_mi_i[j];
	for (k=0; k<numVars; ++k)
	  PCout << std::setw(5) << sm_mi_ij[k];
	PCout << '\n';
      }
    }
  }
  // now print the trial sets
  SparseGridDriver::print_final_sets(converged_within_tol);
}


void HierarchSparseGridDriver::print_smolyak_multi_index() const
{
  size_t i, j, k, num_lev = smolyakMultiIndex.size(), cntr = 1;
  for (i=0; i<num_lev; ++i) {
    const UShort2DArray& sm_mi_i = smolyakMultiIndex[i];
    size_t num_sets = sm_mi_i.size();
    for (j=0; j<num_sets; ++j, ++cntr) {
      PCout << "Smolyak index set " << cntr << ':';
      const UShortArray& sm_mi_ij = sm_mi_i[j];
      for (k=0; k<numVars; ++k)
	PCout << ' ' << sm_mi_ij[k];
      PCout << '\n';
    }
  }
}

} // namespace Pecos

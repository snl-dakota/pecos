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
#include "SharedPolyApproxData.hpp"
#include "sandia_sgmga.hpp"
#include "DistributionParams.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: HierarchSparseGridDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


void HierarchSparseGridDriver::
initialize_grid(unsigned short ssg_level, const RealVector& dim_pref,
		const ShortArray& u_types,
		const ExpansionConfigOptions& ec_options,
		BasisConfigOptions& bc_options, short growth_rate,
		bool track_colloc_indices)
{
  SparseGridDriver::initialize_grid(ssg_level, dim_pref, u_types, ec_options,
				    bc_options, growth_rate);
  trackCollocIndices = track_colloc_indices;
}


void HierarchSparseGridDriver::clear_inactive()
{
  // list erase could be done in two passes...
  //smolyakMultiIndex.erase(smolyakMultiIndex.begin(), smolMIIter);
  //smolyakMultiIndex.erase(++sm_it, smolyakMultiIndex.end());
  
  std::map<UShortArray, UShort3DArray>::iterator sm_it
    = smolyakMultiIndex.begin();
  std::map<UShortArray, UShort4DArray>::iterator ck_it = collocKey.begin();
  std::map<UShortArray, Sizet3DArray>::iterator  ci_it = collocIndices.begin();
  std::map<UShortArray, RealVector2DArray>::iterator t1_it
    = type1WeightSets.begin();
  std::map<UShortArray, RealMatrix2DArray>::iterator t2_it
    = type2WeightSets.begin();
  while (sm_it != smolyakMultiIndex.end())
    if (sm_it == smolMIIter) // preserve active
      { ++sm_it; ++ck_it; ++ci_it; ++t1_it; ++t2_it; }
    else { // clear inactive: postfix increments manage iterator invalidations
      smolyakMultiIndex.erase(sm_it++);
      collocKey.erase(ck_it++);       collocIndices.erase(ci_it++);
      type1WeightSets.erase(t1_it++); type2WeightSets.erase(t2_it++);
    }
}


/*
const UShortArray& HierarchSparseGridDriver::maximal_grid() const
{
  std::map<UShortArray, RealVector2DArray>::const_iterator
    w_cit = type1WeightSets.begin(), max_cit = w_cit;
  size_t l, s, num_lev, num_sets, num_wts, max_wts = 0;
  for (; w_cit!=type1WeightSets.end(); ++w_cit) {
    const RealVector2DArray& t1_wts = w_cit->second;
    num_lev = t1_wts.size(); num_wts = 0;
    for (l=0; l<num_lev; ++l) {
      const RealVectorArray& t1_wts_l = t1_wts[l];
      num_sets = t1_wts_l.size();
      for (s=0; s<num_sets; ++s)
	num_wts += t1_wts_l[s].length();
    }
    if (num_wts > max_wts)
      { max_wts = num_wts; max_cit = w_cit; }
  }
  return max_cit->first;
}
*/


int HierarchSparseGridDriver::grid_size()
{
  // This fn is currently used to ensure change in the number of grid points
  // under a coarse-grain refinement of incrementing ssgLevel.  As such, the
  // logic below does not currently support fine-grained index set increments.

  int& num_colloc_pts = numPtsIter->second;
  if (num_colloc_pts == 0) { // special value indicated update required
    unsigned short ssg_lev = ssgLevIter->second;
    if (collocKeyIter->second.size() == ssg_lev + 1) // coarse-grained updating
      update_collocation_points();
    else {
      update_smolyak_multi_index();
      // rather than full collocKey update, just sum grid sizes:
      UShortArray delta_sizes(numVars);
      unsigned short lev, set, num_sets;
      const UShort3DArray& sm_mi = smolMIIter->second;
      for (lev=0; lev<=ssg_lev; ++lev) {
	const UShort2DArray& sm_mi_l = sm_mi[lev];
	num_sets = sm_mi_l.size();
	for (set=0; set<num_sets; ++set) {
	  levels_to_delta_sizes(sm_mi_l[set], delta_sizes);
	  num_colloc_pts +=
	    SharedPolyApproxData::tensor_product_terms(delta_sizes, false);
	}
      }
    }
  }
  return num_colloc_pts;
}


void HierarchSparseGridDriver::update_smolyak_multi_index(bool clear_sm_mi)
{
  UShort3DArray& sm_mi     = smolMIIter->second;
  unsigned short ssg_lev   = ssgLevIter->second;
  RealVector&    aniso_wts = anisoWtsIter->second;

  if (clear_sm_mi) sm_mi.clear();

  size_t prev_sm_len = sm_mi.size();
  // anisotropic weight updates always accompanied by a level increment, so
  // we are already up to date for both iso and aniso cases if equal:
  if (prev_sm_len == ssg_lev + 1)
    return;

  // this function is for use with isotropic/anisotropic grids, including
  // the initial starting point for a generalized sparse grid adaptation
  bool from_scratch = (prev_sm_len == 0);
  if (!from_scratch && refineControl == DIMENSION_ADAPTIVE_CONTROL_GENERALIZED){
    PCerr << "Error: HierarchSparseGridDriver::update_smolyak_multi_index() "
	  << "intended for use with isotropic and anisotropic grid refinements."
	  << std::endl;
    abort_handler(-1);
  }

  // Populate smolyakMultiIndex based on {iso,aniso}tropic index set constraint:
  // w = q - N = dimension-independent level.  For isotropic,
  //   w + 1 <= |i| <= w + N for i starts at 1 (used for index set defn.)
  //   w - N + 1 <= |j| <= w for j = i - 1 starts at 0 (used for generation)
  // For anisotropic, a weighted linear index set constraint is used.

  size_t lev;
  sm_mi.resize(ssg_lev + 1);
  if (aniso_wts.empty())
    for (lev=prev_sm_len; lev<=ssg_lev; ++lev)
      SharedPolyApproxData::total_order_multi_index(lev, numVars, sm_mi[lev]);
  else { // utilize webbur::sandia_sgmga_vcn_ordered

    // With scaling alpha_min = 1: q_min < |alpha . j| <= q_max.
    IntArray x(numVars), x_max(numVars); //x_max = ssg_lev;
    Real wt_i, q_min = -1., q_max = ssg_lev; // no lower bound for hierarchical
    for (size_t i=0; i<numVars; ++i) {
      wt_i = aniso_wts[i];
      // minimum nonzero weight is scaled to 1, so just catch special case of 0
      x_max[i] = (wt_i > 1.e-10) ? (int)std::ceil(q_max/wt_i) : 0;
    }

    // We take the simple approach and restart vcn iteration from scratch,
    // checking each index_set for current inclusion in smolyakMultiIndex.
    // Warm starting the vcn iteration from the initial smolyakMultiIndex
    // state would be a more efficient option, although also more complex.
    bool more = false;
    UShortArray index_set(numVars);
    Real *aniso_wt_vals = aniso_wts.values();
    int  *x0 = &x[0], *xm0 = &x_max[0];
    webbur::sandia_sgmga_vcn_ordered(numVars, aniso_wt_vals, xm0, x0, q_min,
				     q_max, &more);
    while (more) {
      for (size_t i=0; i<numVars; ++i)
	index_set[i] = (unsigned short)x[i];
      lev = l1_norm(index_set);
      UShort2DArray& sm_mi_l = sm_mi[lev];
      if (from_scratch ||
	  std::find(sm_mi_l.begin(), sm_mi_l.end(), index_set) == sm_mi_l.end())
	sm_mi_l.push_back(index_set);
      webbur::sandia_sgmga_vcn_ordered(numVars, aniso_wt_vals, xm0, x0, q_min,
				       q_max, &more);
    }
  }

#ifdef DEBUG
  PCout << "HierarchSparseGridDriver::update_smolyak_multi_index():\n";
  size_t set, num_sets;
  for (lev=0; lev<=ssg_lev; ++lev) {
    num_sets = sm_mi[lev].size();
    for (set=0; set<num_sets; ++set)
      PCout << "Smolyak multi_index[" << lev << "][" << set << "]:\n"
	    << sm_mi[lev][set];
  }
#endif // DEBUG
}


void HierarchSparseGridDriver::
assign_collocation_key(const UShort3DArray& sm_mi, UShort4DArray& colloc_key,
		       bool ordered)
{
  size_t lev, num_lev = sm_mi.size();
  if (ordered && colloc_key.size() == num_lev) {
    bool rtn = true;
    for (lev=0; lev<num_lev; ++lev)
      if (sm_mi[lev].size() != colloc_key[lev].size())
	{ rtn = false; break; }
    if (rtn) return;
  }

  colloc_key.resize(num_lev);
  size_t set, num_sets;
  if (nestedGrid) {
    UShort2DArray delta_keys(numVars);
    for (lev=0; lev<num_lev; ++lev) {
      const UShort2DArray& sm_mi_l = sm_mi[lev];
      UShort3DArray&         key_l = colloc_key[lev];
      num_sets = sm_mi_l.size();
      key_l.resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	levels_to_delta_keys(sm_mi_l[set], delta_keys);
	SharedPolyApproxData::
	  hierarchical_tensor_product_multi_index(delta_keys, key_l[set]);
      }
    }
  }
  //else
  //  SparseGridDriver::assign_collocation_key();

#ifdef DEBUG
  PCout << "HierarchSparseGridDriver::assign_collocation_key():\n";
  size_t pt, num_tp_pts;
  for (lev=0; lev<num_lev; ++lev) {
    num_sets = colloc_key[lev].size();
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = colloc_key[lev][set].size();
      for (pt=0; pt<num_tp_pts; ++pt)
	PCout << "Collocation key[" << lev << "][" << set << "][" << pt
	      << "]:\n" << colloc_key[lev][set][pt];
    }
  }
#endif // DEBUG
}


void HierarchSparseGridDriver::
update_collocation_key_from_trial(const UShortArray& trial_set,
				  const UShort3DArray& sm_mi,
				  UShort4DArray& colloc_key)
{
  size_t num_lev = sm_mi.size();
  colloc_key.resize(num_lev);
  UShort2DArray delta_keys(numVars), key_ls;
  levels_to_delta_keys(trial_set, delta_keys);

  unsigned short trial_lev = l1_norm(trial_set);
  UShort3DArray& key_l = colloc_key[trial_lev];
  size_t set = key_l.size();
  key_l.push_back(key_ls); // update in place
  SharedPolyApproxData::
    hierarchical_tensor_product_multi_index(delta_keys, key_l[set]);

#ifdef DEBUG
  PCout << "HierarchSparseGridDriver::update_collocation_key_from_trial():\n";
  size_t pt, num_tp_pts = key_l[set].size();
  for (pt=0; pt<num_tp_pts; ++pt)
    PCout << "Collocation key[" << trial_lev << "][" << set << "][" << pt
	  << "]:\n" << key_l[set][pt];
#endif // DEBUG
}


void HierarchSparseGridDriver::
update_collocation_key_from_increment(UShortArray& incr_sets,
				      const UShort3DArray& sm_mi,
				      UShort4DArray& colloc_key)
{
  size_t num_lev = sm_mi.size();
  colloc_key.resize(num_lev);
  UShort2DArray delta_keys(numVars), key_ls;

  // isotropic and anisotropic grid refinements
  // define incr_sets to track iso/aniso grid refinement increment
  size_t lev, set, start_set, num_sets;
  incr_sets.resize(num_lev);
  for (lev=0; lev<num_lev; ++lev)
    incr_sets[lev] = colloc_key[lev].size();
  // update collocKey to correspond to smolyakMultiIndex
  for (lev=0; lev<num_lev; ++lev) {
    const UShort2DArray& sm_mi_l =      sm_mi[lev];
    UShort3DArray&         key_l = colloc_key[lev];
    start_set = incr_sets[lev]; num_sets = sm_mi_l.size();
    for (set=start_set; set<num_sets; ++set) {
      levels_to_delta_keys(sm_mi_l[set], delta_keys);
      key_l.push_back(key_ls); // update in place
      SharedPolyApproxData::
	hierarchical_tensor_product_multi_index(delta_keys, key_l[set]);
    }
  }

#ifdef DEBUG
  PCout<<"HierarchSparseGridDriver::update_collocation_key_from_increment():\n";
  size_t pt, num_tp_pts;
  for (lev=0; lev<num_lev; ++lev) {
    start_set = incr_sets[lev]; num_sets = colloc_key[lev].size();
    for (set=start_set; set<num_sets; ++set) {
      num_tp_pts = colloc_key[lev][set].size();
      for (pt=0; pt<num_tp_pts; ++pt)
	PCout << "Collocation key[" << lev << "][" << set << "][" << pt
	      << "]:\n" << colloc_key[lev][set][pt];
    }
  }
#endif // DEBUG
}


void HierarchSparseGridDriver::
assign_collocation_indices(const UShort4DArray& colloc_key,
			   Sizet3DArray& colloc_indices,
			   int& num_colloc_pts, bool ordered)
{
  size_t lev, num_lev = colloc_key.size();
  if (ordered && colloc_indices.size() == num_lev) {
    bool rtn = true;
    for (lev=0; lev<num_lev; ++lev)
      if (colloc_indices[lev].size() != colloc_key[lev].size())
	{ rtn = false; break; }
    if (rtn) return;
  }

  colloc_indices.resize(num_lev);
  size_t set, pt, cntr = 0, num_sets, num_tp_pts;
  for (lev=0; lev<num_lev; ++lev) {
    const UShort3DArray& key_l = colloc_key[lev];
    num_sets = key_l.size();
    Sizet2DArray& indices_l = colloc_indices[lev];
    indices_l.resize(num_sets);
    for (set=0; set<num_sets; ++set) {
      const UShort2DArray& key_ls = key_l[set];
      num_tp_pts = key_ls.size();
      SizetArray& indices_ls = indices_l[set];
      indices_ls.resize(num_tp_pts);
      for (pt=0; pt<num_tp_pts; ++pt, ++cntr)
	indices_ls[pt] = cntr; // simple sequential ordering for unique pts
    }
  }
  num_colloc_pts = cntr;

#ifdef DEBUG
  PCout << "HierarchSparseGridDriver::assign_collocation_indices():\n"
	<< "num collocation pts = " << num_colloc_pts << '\n';
  for (lev=0; lev<num_lev; ++lev) {
    num_sets = colloc_indices[lev].size();
    for (set=0; set<num_sets; ++set)
      PCout << "Collocation indices[" << lev << "][" << set << "]:\n"
	    << colloc_indices[lev][set];
  }
#endif // DEBUG
}


void HierarchSparseGridDriver::
update_collocation_indices_from_trial(const UShortArray& trial_set,
				      const UShort4DArray& colloc_key,
				      Sizet3DArray& colloc_indices,
				      int& num_colloc_pts)
{
  colloc_indices.resize(colloc_key.size());

  /*
  // TO DO: don't currently need to recompute reference point count since
  // grid_size() not accessible within generalized sparse grid loop, but
  // would be better to harden this as in _from_increment() case...

  unsigned short trial_lev = l1_norm(trial_set);
  size_t pt, num_tp_pts = colloc_key[trial_lev].back().size();
  update_collocation_points(colloc_key, num_colloc_pts);
  size_t cntr = num_colloc_pts - num_tp_pts;

  Sizet2DArray& indices_l = colloc_indices[trial_lev];
  SizetArray indices; indices_l.push_back(indices); // update in place
  SizetArray& trial_indices = indices_l.back();
  trial_indices.resize(num_tp_pts);
  for (pt=0; pt<num_tp_pts; ++pt, ++cntr)
    trial_indices[pt] = cntr;
  */

  unsigned short trial_lev = l1_norm(trial_set);
  size_t pt, num_tp_pts = colloc_key[trial_lev].back().size(),
    cntr = num_colloc_pts;
  Sizet2DArray& indices_l = colloc_indices[trial_lev];
  SizetArray indices; indices_l.push_back(indices); // update in place
  SizetArray& trial_indices = indices_l.back();
  trial_indices.resize(num_tp_pts);
  for (pt=0; pt<num_tp_pts; ++pt, ++cntr)
    trial_indices[pt] = cntr;
  num_colloc_pts += num_tp_pts;

#ifdef DEBUG
  PCout << "HierarchSparseGridDriver::update_collocation_indices_from_trial():"
	<< "\nnum collocation pts = " << num_colloc_pts << '\n';
  size_t set = colloc_indices[trial_lev].size() - 1;
  PCout << "Collocation indices[" << trial_lev << "][" << set << "]:\n"
	<< trial_indices;
#endif // DEBUG
}


void HierarchSparseGridDriver::
update_collocation_indices_from_increment(const UShortArray& incr_sets,
					  const UShort4DArray& colloc_key,
					  Sizet3DArray& colloc_indices,
					  int& num_colloc_pts)
{
  size_t lev, num_lev = colloc_key.size(), set, start_set, num_sets;
  colloc_indices.resize(num_lev);

  // recompute reference pt count due to possible use of grid_size()
  num_colloc_pts = 0;
  for (lev=0; lev<num_lev; ++lev) {
    const UShort3DArray& key_l = colloc_key[lev];
    start_set = incr_sets[lev];
    for (set=0; set<start_set; ++set)
      num_colloc_pts += key_l[set].size();
  }

  // update colloc_indices, num_colloc_pts
  size_t cntr = num_colloc_pts, pt, num_tp_pts;  SizetArray indices;
  for (lev=0; lev<num_lev; ++lev) {
    const UShort3DArray& key_l = colloc_key[lev];
    Sizet2DArray&    indices_l = colloc_indices[lev];
    start_set = incr_sets[lev]; num_sets = key_l.size();
    for (set=start_set; set<num_sets; ++set) {
      indices_l.push_back(indices); // update in place
      SizetArray& trial_indices = indices_l.back();
      num_tp_pts = key_l[set].size();
      trial_indices.resize(num_tp_pts);
      for (pt=0; pt<num_tp_pts; ++pt, ++cntr)
	trial_indices[pt] = cntr;
      num_colloc_pts += num_tp_pts;
    }
  }

#ifdef DEBUG
  PCout << "HierarchSparseGridDriver::update_collocation_indices_from_increment"
	<< "():\nnum collocation pts = " << num_colloc_pts << '\n';
  for (lev=0; lev<num_lev; ++lev) {
    start_set = incr_sets[lev]; num_sets = colloc_indices[lev].size();
    for (set=start_set; set<num_sets; ++set)
      PCout << "Collocation indices[" << lev << "][" << set << "]:\n"
	    << colloc_indices[lev][set];
  }
#endif // DEBUG
}


void HierarchSparseGridDriver::
update_collocation_points(const UShort4DArray& colloc_key, int& num_colloc_pts)
{
  size_t i, num_lev = colloc_key.size();
  num_colloc_pts = 0;
  for (i=0; i<num_lev; ++i) {
    const UShort3DArray& key_i = colloc_key[i];
    size_t j, num_sets = key_i.size();
    for (j=0; j<num_sets; ++j)
      num_colloc_pts += key_i[j].size(); // hierarchical point increments
  }
}


unsigned short HierarchSparseGridDriver::
level_to_delta_size(size_t i, unsigned short level)
{
  switch (level) { // growth restriction should not occur for lev_i = 0 or 1
  case 0: return 1;   break; // 1 new pt
  case 1: return 2;   break; // new ends of 3 pt (open or closed)
  default: { // difference point counts for level & level-1
    unsigned short num_delta, ord_lm1;
    level_to_order(i, level, num_delta); level_to_order(i, level-1, ord_lm1);
    num_delta -= ord_lm1; // Note: num_delta = 0 in case of growth restriction
    return num_delta; break;
  }
  }
}


void HierarchSparseGridDriver::
level_to_delta_key(size_t i, unsigned short lev_i, UShortArray& delta_key_i)
{
  unsigned short num_delta = level_to_delta_size(i, lev_i);
  delta_key_i.resize(num_delta);
  if (!num_delta) return; // possible due to growth restrictions

  switch(collocRules[i]) {
  case GAUSS_PATTERSON: // open nested
    for (size_t j=0; j<num_delta; ++j)
      delta_key_i[j] = 2*j; // new ends + new interior: 0,2,4,6,8,...
    break;
  case NEWTON_COTES: case CLENSHAW_CURTIS: // closed nested
    switch (lev_i) { // growth restriction should not occur for lev_i = 0 or 1
    case 0: delta_key_i[0] = 0;                     break; // center of 1 pt
    case 1: delta_key_i[0] = 0; delta_key_i[1] = 2; break; // ends of 3 pt
    default:
      for (size_t j=0; j<num_delta; ++j)
	delta_key_i[j] = 2*j+1; // new interior: 1,3,5,7,9,...
      break;
    }
    break;
  case GENZ_KEISTER: // open nested table lookup
    // switch on num_delta i/o lev_i due to possibility of growth restriction
    switch (num_delta) {
    case 1: delta_key_i[0] = 0;                     break; // center of 1 pt
    case 2: delta_key_i[0] = 0; delta_key_i[1] = 2; break; // ends of 3 pt
    case 6:
      delta_key_i[0] = 0; delta_key_i[1] = 1; delta_key_i[2] = 3;
      delta_key_i[3] = 5; delta_key_i[4] = 7; delta_key_i[5] = 8;
      break; // 9 pt rule reusing 2, 4 (center), 6
    case 10:
      delta_key_i[0] =  0; delta_key_i[1] =  1; delta_key_i[2] =  3;
      delta_key_i[3] =  5; delta_key_i[4] =  7; delta_key_i[5] = 11;
      delta_key_i[6] = 13; delta_key_i[7] = 15; delta_key_i[8] = 17;
      delta_key_i[9] = 18;
      break; // 19 pt rule reusing 
    case 16:
      delta_key_i[0]  =  0; delta_key_i[1]  =  1; delta_key_i[2]  =  2;
      delta_key_i[3]  =  4; delta_key_i[4]  =  6; delta_key_i[5]  =  8;
      delta_key_i[6]  = 12; delta_key_i[7]  = 16; delta_key_i[8]  = 18;
      delta_key_i[9]  = 22; delta_key_i[10] = 26; delta_key_i[11] = 28;
      delta_key_i[12] = 30; delta_key_i[13] = 32; delta_key_i[14] = 33;
      delta_key_i[15] = 34;
      break; // 35 pt rule reusing
      // 3,5,7,9,10,11,13,14,15,17,19,20,21,23,24,25,27,29,31
    //case 5:  // 43 pt rule augments 19 pt rule, not 35 pt rule
    //  break; // disallow for hierarchical interpolation
    default:
      PCerr << "Error: out of range for hierarchical Genz-Keister rules in "
	    << "HierarchSparseGridDriver::level_to_delta_key()"
	    << std::endl;
      abort_handler(-1);
      break;
    }
    break;
  default:
    PCerr << "Error: bad rule type in level_to_delta_key()" << std::endl;
    abort_handler(-1);
    break;
  }
}


UShortUShortPair HierarchSparseGridDriver::
level_to_delta_pair(size_t i, unsigned short level)
{
  switch (level) { // growth restriction should not occur for level = 0 or 1
  case 0: // +1 pt,  max index = 0 for 1 pt rule
    return UShortUShortPair(1,0); break;
  case 1: // +2 pts, max index = 2 for right end of 3 pt rule (open or closed)
    return UShortUShortPair(2,2); break;
  default: {
    unsigned short max_key, num_delta = level_to_delta_size(i, level);
    if (num_delta)
      switch (collocRules[i]) {
      case GAUSS_PATTERSON:                    // open nested
	max_key = 2 * num_delta - 2; break; // new exterior
      case NEWTON_COTES: case CLENSHAW_CURTIS: // closed nested
	max_key = 2 * num_delta - 1; break; // new interior
      case GENZ_KEISTER: // open nested table lookup
	// switch on num_delta due to possibility of growth restriction
	switch (num_delta) {
	case  6: max_key =  8; break; //  9 pt rule
	case 10: max_key = 18; break; // 19 pt rule
	case 16: max_key = 34; break; // 35 pt rule
	//case 5: // 43 pt rule augments 19 pt rule, not 35 pt rule --> disallow
	default:
	  PCerr << "Error: num_delta (" << num_delta << ") out of range for "
		<< "hierarchical Genz-Keister rules in\n       HierarchSparse"
		<< "GridDriver::level_to_delta_pair()" << std::endl;
	  abort_handler(-1);
	  break;
	}
	break;
      default:
	PCerr << "Error: bad collocation rule type in HierarchSparseGridDriver"
	      << "::level_to_delta_pair()" << std::endl;
	abort_handler(-1);
	break;
      }
    else // num_delta == 0 can occur due to growth restriction
      max_key = std::numeric_limits<unsigned short>::max(); // like _NPOS
    
    return UShortUShortPair(num_delta, max_key);
    break;
  }
  }
}


void HierarchSparseGridDriver::compute_grid(RealMatrix& var_sets)
{
  bool clear = (refineControl != NO_CONTROL); // restore prev state if refined
  update_smolyak_multi_index(clear);          // compute smolyakMultiIndex
  assign_collocation_key();                   // compute collocKey
  assign_1d_collocation_points_weights();     // define 1-D point/weight sets

  // For efficiency reasons, incremental sparse grid definition uses different
  // point orderings than sgmg/sgmga.  Therefore, the reference grid
  // computations are kept completely separate.

  if (nestedGrid) {
    compute_points_weights(var_sets, t1WtIter->second, t2WtIter->second);
    if (trackCollocIndices)
      assign_collocation_indices();
  }
  /*
  else {
    // TO DO: hierarchical interpolation must difference interpolants among
    // full point sets rather than evaluating surpluses at point increments
    reference_unique(var_sets); // define reference grid

#ifdef DEBUG
    PCout << "HierarchSparseGridDriver::compute_grid() results:\n"
	  << "uniqueIndexMapping:\n" << uniqueIndexMapping[activeKey]
	  << "\nvar_sets:\n";
    write_data(PCout, var_sets, false, true, true);
    if (trackUniqueProdWeights) {
      PCout << "\ntype1WeightSets:\n";
      write_data(PCout, t1WtIter->second);
      if (computeType2Weights) {
	PCout << "\ntype2WeightSets:\n";
	write_data(PCout, t2WtIter->second, false, true, true);
      }
    }
#endif
  }
  */

  //update_reference(); // not currently implemented for HSGD
}


void HierarchSparseGridDriver::compute_trial_grid(RealMatrix& var_sets)
{
  // track trial sets that have been evaluated (do here since
  // push_trial_set() used for both new trials and restorations)
  const UShortArray& tr_set = trial_set();
  computedTrialSets[activeKey].insert(tr_set);

  // update collocKey and compute trial variable/weight sets
  update_collocation_key_from_trial(tr_set);
  if (nestedGrid) {
    unsigned short tr_lev = trialLevIter->second;
    RealVector2DArray& t1_wts = t1WtIter->second;
    RealMatrix2DArray& t2_wts = t2WtIter->second;
    if (t1_wts.size() <= tr_lev || t2_wts.size() <= tr_lev)
      { t1_wts.resize(tr_lev+1); t2_wts.resize(tr_lev+1); }
    RealVectorArray& t1_wts_l = t1_wts[tr_lev];
    RealMatrixArray& t2_wts_l = t2_wts[tr_lev];
    size_t set = t1_wts_l.size();
    RealVector t1_wts_ls; t1_wts_l.push_back(t1_wts_ls); // update in place
    RealMatrix t2_wts_ls; t2_wts_l.push_back(t2_wts_ls); // update in place
    compute_points_weights(var_sets, t1_wts_l[set], t2_wts_l[set]);
    if (trackCollocIndices)
      update_collocation_indices_from_trial(tr_set);
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


void HierarchSparseGridDriver::compute_increment(RealMatrix& var_sets)
{
  // update collocKey and compute trial variable/weight sets

  update_smolyak_multi_index();
  UShortArray& incr_sets = incrSetsIter->second;
  update_collocation_key_from_increment(incr_sets);
  size_t lev, set, num_lev = incr_sets.size(), num_sets;
  if (nestedGrid) {
    RealVector2DArray& t1_wts = t1WtIter->second;
    RealMatrix2DArray& t2_wts = t2WtIter->second;
    if (t1_wts.size() < num_lev || t2_wts.size() < num_lev)
      { t1_wts.resize(num_lev); t2_wts.resize(num_lev); }
    // compute total increment evaluations and size var_sets
    size_t num_incr_pts = 0, start_set;
    const UShort4DArray& colloc_key = collocKeyIter->second;
    const UShort3DArray&      sm_mi =    smolMIIter->second;
    for (lev=0; lev<num_lev; ++lev) {
      const UShort3DArray& key_l = colloc_key[lev];
      start_set = incr_sets[lev]; num_sets = key_l.size();
      for (set=start_set; set<num_sets; ++set)
	num_incr_pts += key_l[set].size();
    }
    if (var_sets.numCols() != num_incr_pts)
      var_sets.shapeUninitialized(numVars, num_incr_pts);
    // update type1/2 weights and subset view of points
    size_t cntr = 0, num_tp_pts;
    for (lev=0; lev<num_lev; ++lev) {
      const UShort2DArray& sm_mi_l = sm_mi[lev];
      const UShort3DArray&   key_l = colloc_key[lev];
      RealVectorArray&    t1_wts_l = t1_wts[lev];
      RealMatrixArray&    t2_wts_l = t2_wts[lev];
      start_set = incr_sets[lev]; num_sets = sm_mi_l.size();
      for (set=start_set; set<num_sets; ++set) {
	RealVector t1_wts_ls; t1_wts_l.push_back(t1_wts_ls); // update in place
	RealMatrix t2_wts_ls; t2_wts_l.push_back(t2_wts_ls); // update in place
	const UShort2DArray& key_ls = key_l[set];
	num_tp_pts = key_ls.size();
	RealMatrix pts_ls(Teuchos::View, var_sets, numVars, num_tp_pts, 0,cntr);
	compute_points_weights(sm_mi_l[set], key_ls, pts_ls, t1_wts_l[set],
			       t2_wts_l[set]);
	cntr += num_tp_pts;
      }
    }
    if (trackCollocIndices)
      update_collocation_indices_from_increment(incr_sets);
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


void HierarchSparseGridDriver::
compute_points_weights(const UShortArray& sm_index,
		       const UShort2DArray& colloc_key, RealMatrix& pts,
		       RealVector& t1_wts, RealMatrix& t2_wts)
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
  unsigned short trial_lev = trialLevIter->second;
  compute_points_weights(smolMIIter->second[trial_lev].back(),
			 collocKeyIter->second[trial_lev].back(),
			 pts, t1_wts, t2_wts);
}


void HierarchSparseGridDriver::
compute_points_weights(const UShort3DArray& sm_mi,
		       const UShort4DArray& colloc_key, RealMatrix2DArray& pts,
		       RealVector2DArray& t1_wts, RealMatrix2DArray& t2_wts)
{
  // size consolidated weights according to greatest interpolation depth
  size_t lev, num_lev = sm_mi.size(), set, num_sets;
  pts.resize(num_lev);  t1_wts.resize(num_lev);  t2_wts.resize(num_lev);
  for (lev=0; lev<num_lev; ++lev) {
    const UShort3DArray&   key_l = colloc_key[lev];
    const UShort2DArray& sm_mi_l =  sm_mi[lev];   num_sets = sm_mi_l.size();
    RealMatrixArray&       pts_l =    pts[lev];      pts_l.resize(num_sets);
    RealVectorArray&    t1_wts_l = t1_wts[lev];   t1_wts_l.resize(num_sets);
    RealMatrixArray&    t2_wts_l = t2_wts[lev];   t2_wts_l.resize(num_sets);
    for (set=0; set<num_sets; ++set)
      compute_points_weights(sm_mi_l[set], key_l[set], pts_l[set],
			     t1_wts_l[set], t2_wts_l[set]);
  }
}


/** Points are collapsed as required for compute_grid(var_sets), but t1/t2
    weights are hierarchical 2D arrays. */
void HierarchSparseGridDriver::
compute_points_weights(RealMatrix& pts, RealVector2DArray& t1_wts,
		       RealMatrix2DArray& t2_wts)
{
  const UShort4DArray& colloc_key = collocKeyIter->second;
  const UShort3DArray&      sm_mi =    smolMIIter->second;
  //update_collocation_points(); // should be up to date
  int              num_colloc_pts =    numPtsIter->second;
  size_t i, j, cntr = 0, num_tp_pts, num_lev = colloc_key.size(), num_sets;

  // define points and type 1/2 weights; weights are products of 1D weights
  if (pts.numCols() != num_colloc_pts)
    pts.shapeUninitialized(numVars, num_colloc_pts);
  t1_wts.resize(num_lev);  t2_wts.resize(num_lev);
  for (i=0; i<num_lev; ++i) {
    const UShort3DArray& key_i = colloc_key[i];
    num_sets = key_i.size();
    t1_wts[i].resize(num_sets);  t2_wts[i].resize(num_sets);
    for (j=0; j<num_sets; ++j) {
      const UShort2DArray& key_ij = key_i[j];
      num_tp_pts = key_ij.size();
      // take pts_ij sub-matrix view of full sample matrix pts
      RealMatrix pts_ij(Teuchos::View, pts, numVars, num_tp_pts, 0, cntr);
      compute_points_weights(sm_mi[i][j], key_ij, pts_ij,
			     t1_wts[i][j], t2_wts[i][j]);
      cntr += num_tp_pts;
    }
  }
}


void HierarchSparseGridDriver::initialize_sets()
{
  // define set O (old) from smolyakMultiIndex and smolyakCoeffs:
  const UShort3DArray& sm_mi = smolMIIter->second;
  unsigned short     ssg_lev = ssgLevIter->second;
  UShortArraySet&     old_mi = oldMultiIndex[activeKey];
  old_mi.clear();
  for (unsigned short lev=0; lev<=ssg_lev; ++lev)
    old_mi.insert(sm_mi[lev].begin(), sm_mi[lev].end());

  // computedTrialSets no longer cleared in finalize_sets(), so do on init
  computedTrialSets[activeKey].clear();

  // compute initial set A (active) by applying add_active_neighbors()
  // to the frontier of smolyakMultiIndex:
  if (isotropic()) {
    const UShort2DArray& sm_mi_l = sm_mi[ssg_lev];
    size_t i, num_old_sets = sm_mi_l.size();
    for (i=0; i<num_old_sets; ++i)
      add_active_neighbors(sm_mi_l[i], true); // on frontier
  }
  else { // TO DO
    // For anisotropic, need to compute Pareto set.
  }

#ifdef DEBUG
  PCout << "HierarchSparseGridDriver::initialize_sets():\nold sets:\n" << old_mi
	<< "active sets:\n" << activeMultiIndex[activeKey] << std::endl;
#endif // DEBUG
}


void HierarchSparseGridDriver::push_trial_set(const UShortArray& set)
{
  unsigned short tr_lev = trialLevIter->second = l1_norm(set);
  UShort3DArray&  sm_mi =   smolMIIter->second;
  if (sm_mi.size() <= tr_lev)
    sm_mi.resize(tr_lev+1);
  sm_mi[tr_lev].push_back(set);

  // collocKey, collocIndices, and uniqueIndexMapping updated within
  // either restore_set() or compute_trial_grid()
}


void HierarchSparseGridDriver::restore_set()
{
  // recompute collocKey from trial set
  const UShortArray& tr_set = trial_set();
  update_collocation_key_from_trial(tr_set);
  if (nestedGrid) {
    if (trackCollocIndices)
      update_collocation_indices_from_trial(tr_set);
    // This approach stores less history than WeightSetsRef approach
    unsigned short tr_lev = trialLevIter->second;
    std::map<UShortArray, RealVector>& pop_t1_wts = poppedT1WtSets[activeKey];
    t1WtIter->second[tr_lev].push_back(pop_t1_wts[tr_set]);
    pop_t1_wts.erase(tr_set);
    if (computeType2Weights) {
      std::map<UShortArray, RealMatrix>& pop_t2_wts = poppedT2WtSets[activeKey];
      t2WtIter->second[tr_lev].push_back(pop_t2_wts[tr_set]);
      pop_t2_wts.erase(tr_set);
    }
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


void HierarchSparseGridDriver::pop_trial_set()
{
  UShort4DArray& colloc_key = collocKeyIter->second;
  Sizet3DArray&  colloc_ind = collocIndIter->second;
  UShort3DArray&      sm_mi =    smolMIIter->second;
  int&       num_colloc_pts =    numPtsIter->second;
  unsigned short     tr_lev =  trialLevIter->second;
  if (nestedGrid)
    num_colloc_pts -= colloc_key[tr_lev].back().size(); // subtract trial pts
  /*
  else {
    num_colloc_pts -= numUnique2; // subtract number of trial points
    uniqueIndexMapping.resize(num_colloc_pts); // prune trial set from end
  }
  */

  // migrate weights from popped to active status
  const UShortArray& tr_set = trial_set(); // valid prior to smolyakMI pop
  RealVectorArray& t1_wts_l = t1WtIter->second[tr_lev];
  poppedT1WtSets[activeKey][tr_set] = t1_wts_l.back();
  t1_wts_l.pop_back();
  if (computeType2Weights) {
    RealMatrixArray& t2_wts_l = t2WtIter->second[tr_lev];
    poppedT2WtSets[activeKey][tr_set] = t2_wts_l.back();
    t2_wts_l.pop_back();
  }
  // pop trailing set from smolyakMultiIndex, collocKey, collocIndices
  sm_mi[tr_lev].pop_back(); // tr_set no longer valid
  colloc_key[tr_lev].pop_back();
  if (trackCollocIndices)
    colloc_ind[tr_lev].pop_back();
}


/*
void HierarchSparseGridDriver::merge_increment()
{
  if (nestedGrid) {
    // no-op
  }
  //else
  //  merge_unique();
}
*/


void HierarchSparseGridDriver::
finalize_sets(bool output_sets, bool converged_within_tol, bool reverted)
{
  UShort3DArray& sm_mi = smolMIIter->second;
  unsigned short trial_lev = trialLevIter->second;
  std::map<UShortArray, RealVector>& pop_t1_wts = poppedT1WtSets[activeKey];
  std::map<UShortArray, std::map<UShortArray, RealMatrix> >::iterator pop2_it;
  if (computeType2Weights) pop2_it = poppedT2WtSets.find(activeKey);

  if (output_sets && converged_within_tol) {
    size_t l, s, num_lev = sm_mi.size();
    PCout << "Above tolerance index sets:\n";
    for (l=0; l<num_lev; ++l) {
      const UShort2DArray& sm_mi_l = sm_mi[l];
      size_t num_sets = sm_mi_l.size();
      // omit trial set if not reverted (was below tolerance at convergence)
      if (!reverted && l == trial_lev) --num_sets;
      for (s=0; s<num_sets; ++s)
	print_index_set(PCout, sm_mi_l[s]);
    }
    PCout << "Below tolerance index sets:\n";
    if (!reverted)
      print_index_set(PCout, sm_mi[trial_lev].back());
  }

  // For final answer, push all evaluated sets into old and clear active.
  // Multiple trial insertion approach must be compatible with bookkeeping
  // elsewhere (e.g., Dakota::Approximation), i.e., inc2/inc3 set insertions
  // occur one at a time without mixing.

  // don't insert activeMultiIndex, as this may include sets which have not
  // been evaluated (due to final update_sets() call); use computedTrialSets
  if (nestedGrid) {
    UShortArraySet& computed_trials = computedTrialSets[activeKey];
    UShortArraySet::iterator it;
    for (it=computed_trials.begin(); it!=computed_trials.end(); ++it) {
      const UShortArray& trial_set = *it;
      trial_lev = l1_norm(trial_set);
      sm_mi[trial_lev].push_back(trial_set);
      update_collocation_key_from_trial(trial_set); // update collocKey
      if (trackCollocIndices) // update collocIndices & numCollocPts
	update_collocation_indices_from_trial(trial_set);
      t1WtIter->second[trial_lev].push_back(pop_t1_wts[trial_set]);
      if (computeType2Weights)
	t2WtIter->second[trial_lev].push_back(pop2_it->second[trial_set]);
      if (output_sets && converged_within_tol) // print trials below tol
	print_index_set(PCout, trial_set);
    }
  }
  /*
  else {
    // ...
    // generate final grid, uniqueIndexMapping, collocIndices, numCollocPts
    increment_unique(start_index, false);
    merge_unique();
  }
  */

  if (output_sets && !converged_within_tol) { // print all together in order
    PCout << "Final index sets:\n";
    size_t l, s, num_lev = sm_mi.size();
    for (l=0; l<num_lev; ++l) {
      const UShort2DArray& sm_mi_l = sm_mi[l];
      size_t num_sets = sm_mi_l.size();
      for (s=0; s<num_sets; ++s)
	print_index_set(PCout, sm_mi_l[s]);
    }
  }

  activeMultiIndex[activeKey].clear();
  pop_t1_wts.clear();
  if (computeType2Weights) pop2_it->second.clear();
  // defer since needed for SharedPolyApproxData::finalization_index()
  //computed_trials.clear();
}


void HierarchSparseGridDriver::
partition_keys(UShort2DArray& reference_set_range,
	       UShort2DArray& increment_set_range) const
{
  const UShort3DArray&   sm_mi =   smolMIIter->second;
  unsigned short     trial_lev = trialLevIter->second;
  const UShortArray& incr_sets = incrSetsIter->second;
  size_t lev, num_lev = sm_mi.size(), num_sets;
  reference_set_range.resize(num_lev); increment_set_range.resize(num_lev);
  for (lev=0; lev<num_lev; ++lev) {
    UShortArray&  ref_l = reference_set_range[lev];  ref_l.resize(2);
    UShortArray& incr_l = increment_set_range[lev]; incr_l.resize(2);
    const UShort2DArray& sm_mi_l = sm_mi[lev];
    num_sets = sm_mi_l.size(); ref_l[0] = 0; incr_l[1] = num_sets;
    if (refineControl == DIMENSION_ADAPTIVE_CONTROL_GENERALIZED) {
      if (lev == trial_lev) ref_l[1] = incr_l[0] = num_sets - 1;
      else                  ref_l[1] = incr_l[0] = num_sets;
    }
    else
      ref_l[1] = incr_l[0] = incr_sets[lev];
  }
}


void HierarchSparseGridDriver::
partition_keys(std::map<UShortArray, UShort2DArray>& reference_range_map,
	       std::map<UShortArray, UShort2DArray>& increment_range_map) const
{
  unsigned short     active_trial_lev = trialLevIter->second;
  const UShortArray& active_incr_sets = incrSetsIter->second;

  reference_range_map.clear(); increment_range_map.clear();
  std::map<UShortArray, UShort3DArray>::const_iterator cit;
  size_t lev, num_lev, num_sets;
  UShort2DArray reference_set_range, increment_set_range;
  for (cit=smolyakMultiIndex.begin(); cit!=smolyakMultiIndex.end(); ++cit) {
    const UShort3DArray& sm_mi_i = cit->second;
    num_lev = sm_mi_i.size();
    reference_set_range.resize(num_lev); increment_set_range.resize(num_lev);
    for (lev=0; lev<num_lev; ++lev) {
      UShortArray&  ref_l = reference_set_range[lev];  ref_l.resize(2);
      UShortArray& incr_l = increment_set_range[lev]; incr_l.resize(2);
      const UShort2DArray& sm_mi_il = sm_mi_i[lev];
      num_sets = sm_mi_il.size(); ref_l[0] = 0; incr_l[1] = num_sets;
      if (cit != smolMIIter) // not the active key (no increment)
	ref_l[1] = incr_l[0] = num_sets;
      else if (refineControl == DIMENSION_ADAPTIVE_CONTROL_GENERALIZED) {
	if (lev == active_trial_lev) ref_l[1] = incr_l[0] = num_sets - 1;
	else                         ref_l[1] = incr_l[0] = num_sets;
      }
      else
	ref_l[1] = incr_l[0] = active_incr_sets[lev];
    }
    const UShortArray& key = cit->first;
    reference_range_map[key] = reference_set_range;
    increment_range_map[key] = increment_set_range;
  }
}


void HierarchSparseGridDriver::
partition_keys(UShort3DArray& reference_pt_range,
	       UShort3DArray& increment_pt_range) const
{
  if (refineControl != DIMENSION_ADAPTIVE_CONTROL_GENERALIZED) {
    PCerr << "Error: point set partitioning only supported in HierarchSparse"
	  << "GridDriver::partition_keys() for generalized sparse grids."
	  << std::endl;
    abort_handler(-1);
  }

  const UShort3DArray&      sm_mi =    smolMIIter->second;
  const UShort4DArray& colloc_key = collocKeyIter->second;
  size_t lev, num_lev = colloc_key.size(), set, num_sets, num_tp_pts;
  reference_pt_range.resize(num_lev); increment_pt_range.resize(num_lev);
  for (lev=0; lev<num_lev; ++lev) {
    num_sets = colloc_key[lev].size();
    reference_pt_range[lev].resize(num_sets);
    increment_pt_range[lev].resize(num_sets);
    for (set=0; set<num_sets; ++set) {
      const UShortArray& sm_mi_ls = sm_mi[lev][set];
      UShortArray&  ref_ls = reference_pt_range[lev][set];
      UShortArray& incr_ls = increment_pt_range[lev][set];
      ref_ls.resize(2); incr_ls.resize(2);
      num_tp_pts = colloc_key[lev][set].size();
      ref_ls[0] = 0; incr_ls[1] = num_tp_pts;
      /*
      if (set == trial_set())
	ref_ls[1] = incr_ls[0] = num_tp_pts-1;
      else
      */
	ref_ls[1] = incr_ls[0] = num_tp_pts;
    }
  }
}


/*
void HierarchSparseGridDriver::
combine_weight_sets(const Sizet3DArray& combined_sm_mi_map,
		    RealVector2DArray& comb_t1_wts,
		    RealMatrix2DArray& comb_t2_wts)
{
  // size consolidated weights according to greatest interpolation depth
  size_t i, lev, num_lev, set, num_sets, max_lev = 0, max_sets_il,
    num_map = combined_sm_mi_map.size();
  SizetArray max_sets;
  for (i=0; i<num_map; ++i) {
    const Sizet2DArray& comb_sm_map_i = combined_sm_mi_map[i];
    num_lev = comb_sm_map_i.size();
    if (num_lev > max_lev)                  max_lev = num_lev;
    for (lev=0; lev<num_lev; ++lev) {
      max_sets_il = find_max(comb_sm_map_i[lev]);
      if (max_sets.size() <= lev)           max_sets.push_back(max_sets_il);
      else if (max_sets[lev] < max_sets_il) max_sets[lev] = max_sets_il;
    }
  }
  comb_t1_wts.resize(max_lev);
  comb_t2_wts.resize(max_lev);
  for (lev=0; lev<max_lev; ++lev) {
    comb_t1_wts[lev].resize(max_sets[lev]);
    // be consistent with compute_points_weights(): size for num_{lev,sets}
    // but leave RealMatrix empty if type2 weights are inactive.
    comb_t2_wts.resize(max_sets[lev]);
  }

  std::map<UShortArray, RealVector2DArray>::iterator t1w_it;
  std::map<UShortArray, RealMatrix2DArray>::iterator t2w_it;
  for (t1w_it =type1WeightSets.begin(), t2w_it =type2WeightSets.begin(), i=0;
       t1w_it!=type1WeightSets.end() && t2w_it!=type2WeightSets.end();
       ++t1w_it, ++t2w_it, ++i) {
    const RealVector2DArray& t1w = t1w_it->second;  num_lev  = t1w.size();
    for (lev=0; lev<num_lev; ++lev) {
      const RealVectorArray& t1w_l = t1w[lev];      num_sets = t1w_l.size();
      const RealMatrixArray& t2w_l = t2w_it->second[lev];
      RealVectorArray&  comb_t1w_l = comb_t1_wts[lev];
      RealMatrixArray&  comb_t2w_l = comb_t2_wts[lev];
      const SizetArray& comb_sm_map_il = combined_sm_mi_map[i][lev];
      for (set=0; set<num_sets; ++set) {
	// map from weight sets for this key to corresponding sets in
	// combined t{1,2} weight sets
	size_t comb_set = comb_sm_map_il[set];
	// overlay but no need to accumulate or repeatedly overwrite
	if (comb_t1w_l[comb_set].empty()) comb_t1w_l[comb_set] = t1w_l[set];
	if (computeType2Weights && comb_t2w_l[comb_set].empty())
	  comb_t2w_l[comb_set] = t2w_l[set];
      }
    }
  }
}
*/


/*
const RealVector& HierarchSparseGridDriver::type1_weight_sets() // const
{
  int& num_colloc_pts = numPtsIter->second;
  if (concatT1WeightSets.length() != num_colloc_pts) {
    concatT1WeightSets.sizeUninitialized(num_colloc_pts);
    size_t lev, set, pt, cntr = 0, num_levels = type1WeightSets.size(),
      num_sets, num_tp_pts;
    for (lev=0; lev<num_levels; ++lev) {
      num_sets = type1WeightSets[lev].size();
      for (set=0; set<num_sets; ++set) {
	num_tp_pts = type1WeightSets[lev][set].length();
	const SizetArray& colloc_index = collocIndices[lev][set];
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr)
	  concatT1WeightSets[colloc_index[cntr]]
	    = type1WeightSets[lev][set][pt];
      }
    }
  }
  return concatT1WeightSets;
}


const RealMatrix& HierarchSparseGridDriver::type2_weight_sets() // const
{
  int& num_colloc_pts = numPtsIter->second;
  if (concatT2WeightSets.numCols() != num_colloc_pts) {
    concatT2WeightSets.shapeUninitialized(numVars, num_colloc_pts);
    size_t lev, set, pt, v, cntr = 0, num_levels = type2WeightSets.size(),
      num_sets, num_tp_pts;
    for (lev=0; lev<num_levels; ++lev) {
      num_sets = type2WeightSets[lev].size();
      for (set=0; set<num_sets; ++set) {
	num_tp_pts = type2WeightSets[lev][set].numCols();
	const SizetArray& colloc_index = collocIndices[lev][set];
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	  Real* concat_t2_wts = concatT2WeightSets[colloc_index[cntr]];
	  const Real*  t2_wts = type2WeightSets[lev][set][pt];
	  for (v=0; v<numVars; ++v)
	    concat_t2_wts[v] = t2_wts[v];
	}
      }
    }
  }
  return concatT2WeightSets;
}
*/

} // namespace Pecos

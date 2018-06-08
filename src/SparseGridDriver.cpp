/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 SparseGridDriver
//- Description: Implementation code for SparseGridDriver class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "SparseGridDriver.hpp"
#include "PolynomialApproximation.hpp"
#include "sandia_sgmga.hpp"
#include "sandia_sgmgg.hpp"
#include "DistributionParams.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: SparseGridDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


void SparseGridDriver::assign_1d_collocation_points_weights()
{
  // resize arrays
  unsigned short ssg_lev = ssgLevIter->second;
  size_t i, num_levels = ssg_lev + 1, curr_lev;
  curr_lev = collocPts1D.size();
  if (num_levels > curr_lev) {
    collocPts1D.resize(num_levels);
    for (i=curr_lev; i<num_levels; ++i)
      collocPts1D[i].resize(numVars);
  }
  curr_lev = type1CollocWts1D.size();
  if (num_levels > curr_lev) {
    type1CollocWts1D.resize(num_levels);
    for (i=curr_lev; i<num_levels; ++i)
      type1CollocWts1D[i].resize(numVars);
  }
  curr_lev = type2CollocWts1D.size();
  if (computeType2Weights && num_levels > curr_lev) {
    type2CollocWts1D.resize(num_levels);
    for (i=curr_lev; i<num_levels; ++i)
      type2CollocWts1D[i].resize(numVars);
  }
  // assign values
  // level_index (j indexing) range is 0:w, level (i indexing) range is 1:w+1
  unsigned short l_index, q_order;
  for (i=0; i<numVars; i++)
    for (l_index=0; l_index<num_levels; ++l_index) {
      level_to_order(i, l_index, q_order);
      IntegrationDriver::
	assign_1d_collocation_points_weights(i, q_order, l_index);
    }
}


void SparseGridDriver::dimension_preference(const RealVector& dim_pref)
{
  RealVector aniso_wts;
  if (!dim_pref.empty()) {
    size_t num_pref = dim_pref.length();
    aniso_wts.sizeUninitialized(num_pref);
#ifdef DEBUG
    PCout << "dimension preference:\n"; write_data(PCout, dim_pref);
#endif
    webbur::sandia_sgmga_importance_to_aniso(num_pref, dim_pref.values(),
					     aniso_wts.values());
#ifdef DEBUG
    PCout << "anisotropic weights after sandia_sgmga_importance_to_aniso():\n";
    write_data(PCout, aniso_wts);
#endif
  }
  anisotropic_weights(aniso_wts);
}


void SparseGridDriver::anisotropic_weights(const RealVector& aniso_wts)
{
  if (aniso_wts.empty()) {
    if (!dimIsotropic)          dimIsotropic = updateGridSize = true;
    if (!anisoLevelWts.empty()) anisoLevelWts.sizeUninitialized(0);
  }
  else {
    if (aniso_wts.length() != numVars) {
      PCerr << "Error: length of sparse grid anisotropic weights "
	    << "specification is inconsistent with\n       number of variables "
	    << "in SparseGridDriver::anisotropic_weights()." << std::endl;
      abort_handler(-1);
    }

    size_t i;
    // detect anisotropy
    bool prev_dim_iso = dimIsotropic; // history for updateGridSize
    dimIsotropic = true;
    const Real& wt0 = aniso_wts[0];
    for (i=1; i<numVars; ++i)
      if (std::abs(aniso_wts[i] - wt0) > DBL_EPSILON)
	{ dimIsotropic = false; break; }
    // define updateGridSize and anisoLevelWts
    if (dimIsotropic) {
      if (!prev_dim_iso)          updateGridSize = true;
      if (!anisoLevelWts.empty()) anisoLevelWts.sizeUninitialized(0);
    }
    else {
      RealVector prev_aniso_wts = anisoLevelWts; // history for updateGridSize
      // truncate any negative values
      anisoLevelWts.resize(numVars);
      for (i=0; i<numVars; ++i)
	anisoLevelWts[i] = std::max(aniso_wts[i], 0.);
      // normalize and enforce axis lower bounds/weight upper bounds
      int option = 1; // weights scaled so that minimum nonzero entry is 1
      webbur::sandia_sgmga_aniso_normalize(option, numVars,
					   anisoLevelWts.values());
#ifdef DEBUG
      PCout << "anisoLevelWts after sandia_sgmga_aniso_normalize():\n";
      write_data(PCout, anisoLevelWts);
#endif
      // enforce axis lower bounds, if present, for current ssgLevel.  An axis
      // lower bound defines a weight upper bound based on the current ssgLevel:
      // LB_i = level*wt_min/wt_i --> wt_i = level*wt_min/LB_i and wt_min=1.
      // Catch special case of dim_pref_i = 0 --> wt_i = LB_i = 0.
      if (!axisLowerBounds.empty()) {
	Real ssg_lev = (Real)(ssgLevIter->second);
	for (i=0; i<numVars; ++i)
	  if (axisLowerBounds[i] > 1.e-10) {                 // nonzero LB
	    Real wt_u_bnd = ssg_lev / axisLowerBounds[i];
	    anisoLevelWts[i] = (anisoLevelWts[i] > 1.e-10) ? // nonzero wt
	      std::min(wt_u_bnd, anisoLevelWts[i]) : wt_u_bnd;
	  }
#ifdef DEBUG
	PCout << "anisoLevelWts after axisLowerBounds enforcement:\n";
	write_data(PCout, anisoLevelWts);
#endif
      }
      // define updateGridSize
      if (anisoLevelWts != prev_aniso_wts)
	updateGridSize = true;
    }
  }
}


void SparseGridDriver::update_axis_lower_bounds()
{
  if (axisLowerBounds.empty())
    axisLowerBounds.sizeUninitialized(numVars);
  // An axisLowerBound is the maximum index coverage achieved on a coordinate
  // axis (when all other indices are zero); it defines a constraint for
  // minimum coordinate coverage in future refinements.  The linear index set
  // constraint is level*wt_min-|wt| < j.wt <= level*wt_min, which becomes
  // level-|wt| < j_i w_i <= level for wt_min=1 and all other indices=0.
  // The max feasible j_i is then level/w_i (except for special case w_i=0).
  Real ssg_lev = (Real)(ssgLevIter->second);
  if (dimIsotropic)
    axisLowerBounds = ssg_lev; // all weights = 1
  else // min nonzero weight scaled to 1 --> just catch special case w_i=0
    for (size_t i=0; i<numVars; ++i)
      axisLowerBounds[i] = (anisoLevelWts[i] > 1.e-10) ? // nonzero wt
	ssg_lev / anisoLevelWts[i] : 0.;
}


void SparseGridDriver::
initialize_grid(unsigned short ssg_level, const RealVector& dim_pref,
		const ShortArray& u_types,
		const ExpansionConfigOptions& ec_options,
		BasisConfigOptions& bc_options, short growth_rate)
{
  growthRate             = growth_rate;
  //refineType           = ec_options.refinementType;
  refineControl          = ec_options.refinementControl;

  // For unrestricted exponential growth, use of nested rules is restricted
  // to uniform/normal in order to enforce similar growth rates:
  if (bc_options.nestedRules && growthRate == UNRESTRICTED_GROWTH) {
    size_t i, num_u_types = u_types.size(); // numVars not yet defined
    for (i=0; i<num_u_types; ++i)
      if (u_types[i] != STD_UNIFORM && u_types[i] != STD_NORMAL)
	{ bc_options.nestedRules = false; break; }
  }
  // For MODERATE and SLOW restricted exponential growth, nested rules
  // can be used heterogeneously and synchronized with STANDARD and SLOW
  // linear growth, respectively.

  IntegrationDriver::initialize_grid(u_types, ec_options, bc_options);

  level(ssg_level);
  dimension_preference(dim_pref);
}


void SparseGridDriver::precompute_rules()
{
  unsigned short l, m, ssg_lev = ssgLevIter->second;
  if (dimIsotropic)
    for (size_t i=0; i<numVars; ++i) {
      level_to_order(i, ssg_lev, m); // max order is full level in this dim
      polynomialBasis[i].precompute_rules(m);
    }
  else
    for (size_t i=0; i<numVars; ++i) {
      Real& wt_i = anisoLevelWts[i];
      l = (wt_i > 0.) ? (unsigned short)((Real)ssg_lev / wt_i) : 0;
      level_to_order(i, l, m); // max order is full aniso level[dim]
      polynomialBasis[i].precompute_rules(m);
    }
}


void SparseGridDriver::initialize_sets()
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "initialize_sets()." << std::endl;
  abort_handler(-1);
}


void SparseGridDriver::push_trial_set(const UShortArray& set)
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "push_trial_set()." << std::endl;
  abort_handler(-1);
}


void SparseGridDriver::restore_set()
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "()." << std::endl;
  abort_handler(-1);
}


void SparseGridDriver::pop_trial_set()
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "pop_trial_set()." << std::endl;
  abort_handler(-1);
}


void SparseGridDriver::
finalize_sets(bool output_sets, bool converged_within_tol)
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "finalize_sets()." << std::endl;
  abort_handler(-1);
}


void SparseGridDriver::update_reference()
{
  // Not needed for HierarchSparseGridDriver, so use no-op as default

  /*
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "update_reference()." << std::endl;
  abort_handler(-1);
  */
}


void SparseGridDriver::compute_trial_grid(RealMatrix& var_sets)
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "compute_trial_grid()." << std::endl;
  abort_handler(-1);
}


void SparseGridDriver::compute_increment(RealMatrix& var_sets)
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "compute_increment()." << std::endl;
  abort_handler(-1);
}


void SparseGridDriver::push_increment()
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "push_increment()." << std::endl;
  abort_handler(-1);
}


void SparseGridDriver::pop_increment()
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "pop_increment()." << std::endl;
  abort_handler(-1);
}


void SparseGridDriver::merge_increment()
{
  // Not needed for HierarchSparseGridDriver, so use no-op as default

  /*
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "merge_increment()." << std::endl;
  abort_handler(-1);
  */
}


const UShortArray& SparseGridDriver::trial_set() const
{
  PCerr << "Error: no default implementation for SparseGridDriver::trial_set()."
	<< std::endl;
  abort_handler(-1);
  return activeKey; // dummy UShortArray
}


int SparseGridDriver::unique_trial_points() const
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "unique_trial_points()." << std::endl;
  abort_handler(-1);
  return 0;
}


void SparseGridDriver::update_smolyak_arrays()
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "update_smolyak_arrays()." << std::endl;
  abort_handler(-1);
}


void SparseGridDriver::update_sets(const UShortArray& set_star)
{
  // set_star is passed as *cit_star from the best entry in activeMultiIndex.
  // Therefore, we must use caution in updates to activeMultiIndex that can
  // invalidate cit_star.

  // update evaluation set smolyakMultiIndex (permanently, will not be popped)
  push_trial_set(set_star);
  restore_set();     // calls increment_unique() --> INC2
  merge_increment(); // calls merge_unique()     --> INC3

  // use trial set rather than incoming set_star due to iterator invalidation
  const UShortArray&       tr_set = trial_set();
  UShortArraySet& computed_trials = computedTrialSets[activeKey];
  UShortArraySet&       active_mi =  activeMultiIndex[activeKey];
  UShortArraySet&          old_mi =     oldMultiIndex[activeKey];

  // update set O by adding the trial set to oldMultiIndex:
  old_mi.insert(tr_set);
  // remove the trial set from set A by erasing from activeMultiIndex:
  active_mi.erase(tr_set); // invalidates cit_star -> set_star
  // update subset of A that have been evaluated as trial sets
  if (!computed_trials.empty()) // not tracked for LightweightSparseGridDriver
    computed_trials.erase(tr_set);

  // update set A (activeMultiIndex) based on neighbors of trial set
  add_active_neighbors(tr_set, false);//dimIsotropic);

  // Note: pruning irrelevant sets that have Coeff = 0 would be tricky,
  //       since a 0 close to the frontier can become nonzero

#ifdef DEBUG
  PCout << "Sets updated: (Smolyak,Old,Active,Trial) = (" << smolyak_size()
	<< ',' << old_mi.size() << ',' << active_mi.size() << ','
	<< computed_trials.size() << ')' << std::endl;
#endif // DEBUG
}


void SparseGridDriver::
add_active_neighbors(const UShortArray& set, bool frontier)
{
  UShortArray trial_set = set;
  UShortArraySet&    old_mi =    oldMultiIndex[activeKey];
  UShortArraySet& active_mi = activeMultiIndex[activeKey];
  UShortArraySet::const_iterator cit;
  size_t i, j, num_v = set.size();
  for (i=0; i<num_v; ++i) {
    // i^{th} candidate for set A (active) computed from forward neighbor:
    // increment by 1 in dimension i
    unsigned short& trial_set_i = trial_set[i];
    ++trial_set_i;
    // if !frontier, then candidates could exist in oldMultiIndex
    if (frontier || old_mi.find(trial_set) == old_mi.end()) {
      // test all backwards neighbors for membership in set O (old)
      bool backward_old = true;
      for (j=0; j<num_v; ++j) {
	unsigned short& trial_set_j = trial_set[j];
	if (trial_set_j) { // if 0, then admissible by default
	  --trial_set_j;
	  cit = old_mi.find(trial_set);
	  ++trial_set_j; // restore
	  if (cit == old_mi.end())
	    { backward_old = false; break; }
	}
      }
      if (backward_old) // std::set<> will discard any active duplicates
	active_mi.insert(trial_set);
    }
    --trial_set_i; // restore
  }
}


int SparseGridDriver::level_to_order_exp_hgk_interp(int level, int growth)
{
  if (level == 0) return 1;

  switch (growth) {
  case SLOW_RESTRICTED_GROWTH: {
    unsigned short level = 0, max_level = 5;
    int m = 1, m_goal = level + 1;
    while (level <= max_level && m < m_goal)
      { ++level; m = orderGenzKeister[level]; }
    return m; break;
  }
  case MODERATE_RESTRICTED_GROWTH: {
    unsigned short level = 0, max_level = 5;
    int m = 1, m_goal = 2 * level + 1;
    while (level <= max_level && m < m_goal)
      { ++level; m = orderGenzKeister[level]; }
    return m; break;
  }
  case UNRESTRICTED_GROWTH:
    return orderGenzKeister[std::min(level, 5)]; break;
  }
}


int SparseGridDriver::level_to_order_exp_closed_interp(int level, int growth)
{
  if (level == 0) return 1;

  switch (growth) {
  case SLOW_RESTRICTED_GROWTH: {
    int m = 1, lev_pow = 1, m_goal = level + 1;
    while (m < m_goal)
      { lev_pow *= 2; m = lev_pow + 1; } //std::pow(2.,level) + 1;
    return m; break;
  }
  case MODERATE_RESTRICTED_GROWTH: {
    int m = 1, lev_pow = 1, m_goal = 2 * level + 1;
    while (m < m_goal)
      { lev_pow *= 2; m = lev_pow + 1; } //std::pow(2.,level) + 1;
    return m; break;
  }
  case UNRESTRICTED_GROWTH:
    return (int)std::pow(2., level) + 1; break;
  }
}


int SparseGridDriver::level_to_order_exp_open_interp(int level, int growth)
{
  if (level == 0) return 1;

  switch (growth) {
  case SLOW_RESTRICTED_GROWTH: {
    int m = 1, lev_pow = 2, m_goal = level + 1;
    while (m < m_goal)
      { lev_pow *= 2; m = lev_pow - 1; } //std::pow(2.,level+1) - 1;
    return m; break;
  }
  case MODERATE_RESTRICTED_GROWTH: {
    int m = 1, lev_pow = 2, m_goal = 2 * level + 1;
    while (m < m_goal)
      { lev_pow *= 2; m = lev_pow - 1; } //std::pow(2.,level+1) - 1;
    return m; break;
  }
  case UNRESTRICTED_GROWTH:
    return (int)std::pow(2., level+1) - 1; break;
  }
}

} // namespace Pecos

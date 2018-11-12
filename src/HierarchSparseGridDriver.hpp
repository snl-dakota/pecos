/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 HierarchSparseGridDriver
//- Description: Wrapper class for C++ code from packages/quadrature/sparse_grid
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef HIERARCH_SPARSE_GRID_DRIVER_HPP
#define HIERARCH_SPARSE_GRID_DRIVER_HPP

#include "SparseGridDriver.hpp"

namespace Pecos {


/// Derived integration driver class that generates N-dimensional
/// Smolyak sparse grids for numerical evaluation of expectation
/// integrals over independent standard random variables.

/** This class is used by Dakota::NonDSparseGrid, but could also be
    used for general numerical integration of moments.  It employs 1-D
    Clenshaw-Curtis, Newton-Cotes, and Gaussian quadrature rules
    within Smolyak sparse grids. */

class HierarchSparseGridDriver: public SparseGridDriver
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  HierarchSparseGridDriver();
  /// constructor
  HierarchSparseGridDriver(unsigned short ssg_level,
			   const RealVector& dim_pref = RealVector(),
			   short growth_rate = MODERATE_RESTRICTED_GROWTH,
			   short refine_control = NO_CONTROL);
  /// destructor
  ~HierarchSparseGridDriver();

  //
  //- Heading: Virtual function redefinitions
  //

  void compute_grid(RealMatrix& var_sets);
  int grid_size();

  void clear_inactive();
  void clear_keys();

  //const UShortArray& maximal_grid() const;

  void initialize_sets();
  void increment_smolyak_multi_index(const UShortArray& set);
  bool push_trial_available(const UShortArray& key, const UShortArray& tr_set);
  bool push_trial_available(const UShortArray& key);
  bool push_trial_available();
  size_t push_trial_index(const UShortArray& key, const UShortArray& tr_set);
  size_t push_trial_index(const UShortArray& key);
  size_t push_trial_index();
  size_t push_index() const;
  void push_set();
  void compute_trial_grid(RealMatrix& var_sets);
  void pop_set();
  void finalize_sets(bool output_sets, bool converged_within_tol,
		     bool reverted);

  const UShortArray& trial_set(const UShortArray& key) const;
  const UShortArray& trial_set() const;
  unsigned short trial_level() const;
  int unique_trial_points() const;

  void compute_increment(RealMatrix& var_sets);
  void push_increment();
  void pop_increment();
  //void merge_increment();

  void print_smolyak_multi_index() const;

  // concatenate type1WeightSets for use in abstract integration functions
  //const RealVector& type1_weight_sets(); // const;
  // concatenate type2WeightSets for use in abstract integration functions
  //const RealMatrix& type2_weight_sets(); // const;

  UShortUShortPair level_to_delta_pair(size_t i, unsigned short lev_i);
  unsigned short level_to_delta_size(size_t i, unsigned short lev_i);
  void level_to_delta_key(size_t i, unsigned short lev_i,
			  UShortArray& delta_key_i);

  /// convert a Smolyak index set into hierarchical quadrature keys
  void levels_to_delta_keys(const UShortArray& levels,
			    UShort2DArray& delta_keys);
  /// convert a Smolyak index set into the sizes of hierarchical
  /// quadrature increments
  void levels_to_delta_sizes(const UShortArray& levels,
			     UShortArray& delta_sizes);

  //
  //- Heading: Member functions
  //

  /// initialize all sparse grid settings except for distribution params
  void initialize_grid(unsigned short ssg_level, const RealVector& dim_pref,
		       const ShortArray& u_types,
		       const ExpansionConfigOptions& ec_options,
		       BasisConfigOptions& bc_options,
		       short growth_rate = MODERATE_RESTRICTED_GROWTH,
		       bool track_colloc_indices = true);

  void assign_collocation_key();
  void assign_collocation_key(const UShort3DArray& sm_mi,
			      UShort4DArray& colloc_key, bool ordered = true);
  void update_collocation_key_from_trial(const UShortArray& trial_set);
  void update_collocation_key_from_trial(const UShortArray& trial_set,
					 const UShort3DArray& sm_mi,
					 UShort4DArray& colloc_key);
  void update_collocation_key_from_increment(UShortArray& incr_sets);
  void update_collocation_key_from_increment(UShortArray& incr_sets,
					     const UShort3DArray& sm_mi,
					     UShort4DArray& colloc_key);

  void assign_collocation_indices();
  void assign_collocation_indices(const UShort4DArray& colloc_key,
				  Sizet3DArray& colloc_indices,
				  int& num_colloc_pts, bool ordered = true);
  void update_collocation_indices_from_trial(const UShortArray& trial_set);
  void update_collocation_indices_from_trial(const UShortArray& trial_set,
					     const UShort4DArray& colloc_key,
					     Sizet3DArray& colloc_indices,
					     int& num_colloc_pts);
  void update_collocation_indices_from_increment(const UShortArray& incr_sets);
  void update_collocation_indices_from_increment(const UShortArray& incr_sets,
					     const UShort4DArray& colloc_key,
					     Sizet3DArray& colloc_indices,
					     int& num_colloc_pts);

  /// update active numCollocPts from active collocKey (contains unique points)
  void update_collocation_points();
  /// update num_colloc_pts from colloc_key (contains unique points)
  void update_collocation_points(const UShort4DArray& colloc_key,
				 int& num_colloc_pts);

  /// return active entry in incrementSets
  const UShortArray& increment_sets() const;

  // return number of sets within popped{T1,T2}WtSets identified by key
  //size_t popped_sets(const UShortArray& key) const;
  // return number of sets within popped{T1,T2}WtSets identified by activeKey
  //size_t popped_sets() const;

  /// return active entry in smolyakMultiIndex
  const UShort3DArray& smolyak_multi_index() const;
  /// set active entry in smolyakMultiIndex
  void smolyak_multi_index(const UShort3DArray& sm_mi);
  /// return smolyakMultiIndex[key]
  const UShort3DArray& smolyak_multi_index(const UShortArray& key) const;
  /// return smolyakMultiIndex
  const std::map<UShortArray, UShort3DArray>& smolyak_multi_index_map() const;

  /// set trackCollocIndices
  void track_collocation_indices(bool track_colloc_indices);
  /// get trackCollocIndices
  bool track_collocation_indices() const;

  /// return active entry in collocKey
  const UShort4DArray& collocation_key() const;
  /// set active entry in collocKey
  void collocation_key(const UShort4DArray& key);
  /// return collocKey[key]
  const UShort4DArray& collocation_key(const UShortArray& key) const;
  /// return collocKey
  const std::map<UShortArray, UShort4DArray>& collocation_key_map() const;

  /// return active entry in collocIndices
  const Sizet3DArray& collocation_indices() const;
  /// set active entry in collocIndices
  void collocation_indices(const Sizet3DArray& indices);
  /// return collocIndices[key]
  const Sizet3DArray& collocation_indices(const UShortArray& key) const;

  /// discriminate portions of the active level-set hierarchy that are
  /// reference sets from those in the current increment
  void partition_keys(UShort2DArray& reference_set_range,
		      UShort2DArray& increment_set_range) const;
  /// discriminate portions of the level-set hierarchy that are
  /// reference sets from those in the current increment
  void partition_keys(std::map<UShortArray, UShort2DArray>& reference_range_map,
    std::map<UShortArray, UShort2DArray>& increment_range_map) const;
  /// discriminate portions of the active level-set-point hierarchy that are
  /// in the reference grid from those in the current increment
  void partition_keys(UShort3DArray& reference_pt_range,
		      UShort3DArray& increment_pt_range) const;

  /// compute points and weights for all levels and sets of the hierarchical
  /// sparse grid indicated by Smolyak multi-index and collocation key
  void compute_points_weights(const UShort3DArray& sm_mi,
			      const UShort4DArray& colloc_key,
			      RealMatrix2DArray& pts, RealVector2DArray& t1_wts,
			      RealMatrix2DArray& t2_wts);
  /// compute points and weights for a trial set
  void compute_points_weights(RealMatrix& pts, RealVector& t1_wts,
			      RealMatrix& t2_wts);
  /// compute points and weights for all levels of the active sparse grid;
  /// in this case the points array is condensed to a single matrix but the
  /// weights remain organized by hierarchical levels and sets
  void compute_points_weights(RealMatrix& pts, RealVector2DArray& t1_wts,
			      RealMatrix2DArray& t2_wts);

  // overlay all type{1,2}WeightSets and store in active key
  //void combine_weight_sets(const Sizet3DArray& combined_sm_mi_map,
  //			     RealVector2DArray& comb_t1_wts,
  //			     RealMatrix2DArray& comb_t2_wts);

  /// return active type1WeightSets for use in hierarchical integration
  const RealVector2DArray& type1_hierarchical_weight_sets() const;
  /// set active type1WeightSets for use in hierarchical integration
  void type1_hierarchical_weight_sets(const RealVector2DArray& t1_wts);
  /// return active type2WeightSets for use in hierarchical integration
  const RealMatrix2DArray& type2_hierarchical_weight_sets() const;
  /// set active type2WeightSets for use in hierarchical integration
  void type2_hierarchical_weight_sets(const RealMatrix2DArray& t2_wts);

  /// return type1WeightSets for use in hierarchical integration
  const std::map<UShortArray, RealVector2DArray>& type1_weight_sets_map() const;
  /// return type2WeightSets for use in hierarchical integration
  const std::map<UShortArray, RealMatrix2DArray>& type2_weight_sets_map() const;

private:

  //
  //- Heading: Convenience functions
  //

  /// update {smolMI,collocKey,collocInd}Iter from activeKey
  void update_active_iterators();

  /// update active smolyakMultiIndex for change in level and/or aniso weights
  void update_smolyak_multi_index(bool clear_sm_mi = false);

  /// moves all data from popped weights to active arrays
  void push_popped_weights();

  /// kernel routine used for computing points and weights for a tensor grid
  /// corresponding to a single index set
  void compute_points_weights(const UShortArray& sm_index,
			      const UShort2DArray& colloc_key, RealMatrix& pts,
			      RealVector& t1_wts, RealMatrix& t2_wts);

  //
  //- Heading: Data
  //

  /// flag for use of fully nested 1D rules, allowing formulation using
  /// collocation point increments
  bool nestedGrid;

  /// due to the hierarchical structure, collocation indices only need
  /// to be defined in special cases (e.g., generalized sparse grids
  /// for which index sets can appear in different orders).
  bool trackCollocIndices;

  /// interpolation depth by index set by numVars array for identifying
  /// the index to use within the polynomialBasis for a particular variable
  /** The index sets correspond to j (0-based) for use as indices, which
      are offset from the i indices (1-based) normally used in the Smolyak
      expressions.  The indices correspond to levels, one within each
      anisotropic tensor-product integration of a Smolyak recursion. */
  std::map<UShortArray, UShort3DArray> smolyakMultiIndex;
  /// iterator for active entry within smolyakMultiIndex
  std::map<UShortArray, UShort3DArray>::iterator smolMIIter;

  /// level of trial evaluation set passed to increment_smolyak_multi_index()
  /// during an index set-based refinement
  std::map<UShortArray, unsigned short> trialLevel;
  /// iterator for active entry within trialLevel
  std::map<UShortArray, unsigned short>::iterator trialLevIter;
  /// identifies the trailing index set increments within smolyakMultiIndex
  /// due to an isotropic/anistropic grid refinement
  std::map<UShortArray, UShortArray> incrementSets;
  /// iterator for active entry within incrementSets
  std::map<UShortArray, UShortArray>::iterator incrSetsIter;

  /// levels-by-index sets-by-numDeltaPts-by-numVars array for identifying
  /// the 1-D point indices for sets of tensor-product collocation points
  std::map<UShortArray, UShort4DArray> collocKey;
  /// iterator for active entry within collocKey
  std::map<UShortArray, UShort4DArray>::iterator collocKeyIter;

  /// levels-by-index sets-by-numTensorProductPts array for linking the
  /// set of tensor products to the unique collocation points evaluated
  std::map<UShortArray, Sizet3DArray> collocIndices;
  /// iterator for active entry within collocIndices
  std::map<UShortArray, Sizet3DArray>::iterator collocIndIter;

  /// the set of type1 weights (for integration of value interpolants)
  /// associated with each point in the sparse grid
  std::map<UShortArray, RealVector2DArray> type1WeightSets;
  /// iterator for active entry within type1WeightSets
  std::map<UShortArray, RealVector2DArray>::iterator t1WtIter;
  /// the set of type2 weights (for integration of gradient interpolants)
  /// for each derivative component and for each point in the sparse grid
  std::map<UShortArray, RealMatrix2DArray> type2WeightSets;
  /// iterator for active entry within type2WeightSets
  std::map<UShortArray, RealMatrix2DArray>::iterator t2WtIter;

  // concatenation of type1WeightSets RealVector2DArray into a RealVector
  //RealVector concatT1WeightSets;
  // concatenation of type2WeightSets RealMatrix2DArray into a RealMatrix
  //RealMatrix concatT2WeightSets;

  /// popped trial sets that were computed but not selected
  std::map<UShortArray, UShortArrayDequeArray> poppedLevMultiIndex;
  /// index into poppedLevMultiIndex[trial_lev] for data to be restored
  size_t pushIndex;

  /// type 1 weight sets popped during decrement for later restoration to
  /// type1WeightSets
  std::map<UShortArray, RealVectorDequeArray> poppedT1WtSets;
  /// type 2 weight sets popped during decrement for later restoration to
  /// type2WeightSets
  std::map<UShortArray, RealMatrixDequeArray> poppedT2WtSets;
};


inline HierarchSparseGridDriver::HierarchSparseGridDriver():
  SparseGridDriver(), nestedGrid(true), trackCollocIndices(true),
  smolMIIter(smolyakMultiIndex.end()), pushIndex(_NPOS)
{ update_active_iterators(); }


inline HierarchSparseGridDriver::
HierarchSparseGridDriver(unsigned short ssg_level, const RealVector& dim_pref,
			 short growth_rate, short refine_control):
  SparseGridDriver(ssg_level, dim_pref, growth_rate, refine_control),
  nestedGrid(true), trackCollocIndices(true),
  smolMIIter(smolyakMultiIndex.end()), pushIndex(_NPOS)
{ update_active_iterators(); }


inline HierarchSparseGridDriver::~HierarchSparseGridDriver()
{ }


inline void HierarchSparseGridDriver::update_active_iterators()
{
  // Test for change
  if (smolMIIter != smolyakMultiIndex.end() && smolMIIter->first == activeKey)
    return;

  smolMIIter = smolyakMultiIndex.find(activeKey);
  if (smolMIIter == smolyakMultiIndex.end()) {
    std::pair<UShortArray, UShort3DArray> u3a_pair(activeKey, UShort3DArray());
    smolMIIter = smolyakMultiIndex.insert(u3a_pair).first;
  }
  trialLevIter = trialLevel.find(activeKey);
  if (trialLevIter == trialLevel.end()) {
    std::pair<UShortArray, unsigned short> us_pair(activeKey, 0);
    trialLevIter = trialLevel.insert(us_pair).first;
  }
  incrSetsIter = incrementSets.find(activeKey);
  if (incrSetsIter == incrementSets.end()) {
    std::pair<UShortArray, UShortArray> ua_pair(activeKey, UShortArray());
    incrSetsIter = incrementSets.insert(ua_pair).first;
  }
  collocKeyIter = collocKey.find(activeKey);
  if (collocKeyIter == collocKey.end()) {
    std::pair<UShortArray, UShort4DArray> u4a_pair(activeKey, UShort4DArray());
    collocKeyIter = collocKey.insert(u4a_pair).first;
  }
  collocIndIter = collocIndices.find(activeKey);
  if (collocIndIter == collocIndices.end()) {
    std::pair<UShortArray, Sizet3DArray> s3a_pair(activeKey, Sizet3DArray());
    collocIndIter = collocIndices.insert(s3a_pair).first;
  }
  t1WtIter = type1WeightSets.find(activeKey);
  if (t1WtIter == type1WeightSets.end()) {
    std::pair<UShortArray, RealVector2DArray>
      rv2_pair(activeKey, RealVector2DArray());
    t1WtIter = type1WeightSets.insert(rv2_pair).first;
  }
  t2WtIter = type2WeightSets.find(activeKey);
  if (t2WtIter == type2WeightSets.end()) {
    std::pair<UShortArray, RealMatrix2DArray>
      rm2_pair(activeKey, RealMatrix2DArray());
    t2WtIter = type2WeightSets.insert(rm2_pair).first;
  }

  SparseGridDriver::update_active_iterators();
}


inline void HierarchSparseGridDriver::assign_collocation_key()
{ assign_collocation_key(smolMIIter->second, collocKeyIter->second); }


inline void HierarchSparseGridDriver::
update_collocation_key_from_trial(const UShortArray& trial_set)
{
  update_collocation_key_from_trial(trial_set, smolMIIter->second,
				    collocKeyIter->second);
}


inline void HierarchSparseGridDriver::
update_collocation_key_from_increment(UShortArray& incr_sets)
{
  update_collocation_key_from_increment(incr_sets, smolMIIter->second,
					collocKeyIter->second);
}


inline void HierarchSparseGridDriver::assign_collocation_indices()
{
  assign_collocation_indices(collocKeyIter->second, collocIndIter->second,
			     numPtsIter->second);
}


inline void HierarchSparseGridDriver::
update_collocation_indices_from_trial(const UShortArray& trial_set)
{
  update_collocation_indices_from_trial(trial_set, collocKeyIter->second,
					collocIndIter->second,
					numPtsIter->second);
}


inline void HierarchSparseGridDriver::
update_collocation_indices_from_increment(const UShortArray& incr_sets)
{
  update_collocation_indices_from_increment(incr_sets, collocKeyIter->second,
					    collocIndIter->second,
					    numPtsIter->second);
}


inline void HierarchSparseGridDriver::update_collocation_points()
{ update_collocation_points(collocKeyIter->second, numPtsIter->second); }


inline void HierarchSparseGridDriver::clear_keys()
{
  SparseGridDriver::clear_keys();

  smolyakMultiIndex.clear();  smolMIIter    = smolyakMultiIndex.end();
  trialLevel.clear();         trialLevIter  = trialLevel.end();
  incrementSets.clear();      incrSetsIter  = incrementSets.end();
  collocKey.clear();          collocKeyIter = collocKey.end();
  collocIndices.clear();      collocIndIter = collocIndices.end();
  type1WeightSets.clear();    t1WtIter      = type1WeightSets.end();
  type2WeightSets.clear();    t2WtIter      = type2WeightSets.end();

  poppedLevMultiIndex.clear(); poppedT1WtSets.clear(); poppedT2WtSets.clear();
}


inline const UShortArray& HierarchSparseGridDriver::
trial_set(const UShortArray& key) const
{
  std::map<UShortArray, UShort3DArray>::const_iterator sm_cit
    = smolyakMultiIndex.find(key);
  std::map<UShortArray, unsigned short>::const_iterator tl_cit
    = trialLevel.find(key);
  if (sm_cit == smolyakMultiIndex.end() || tl_cit == trialLevel.end()) {
    PCerr << "Error: key not found in HierarchSparseGridDriver::trial_set()"
	  << std::endl;
    abort_handler(-1);
  }
  return sm_cit->second[tl_cit->second].back();
}


inline const UShortArray& HierarchSparseGridDriver::trial_set() const
{ return smolMIIter->second[trialLevIter->second].back(); }


inline unsigned short HierarchSparseGridDriver::trial_level() const
{ return trialLevIter->second; }


inline int HierarchSparseGridDriver::unique_trial_points() const
{ return collocKeyIter->second[trialLevIter->second].back().size(); }


/** identify if newly-pushed trial set exists within stored data sets */
inline bool HierarchSparseGridDriver::
push_trial_available(const UShortArray& key, const UShortArray& tr_set)
{
  size_t tr_lev = l1_norm(tr_set);
  const UShortArrayDeque& pop_mi_l = poppedLevMultiIndex[key][tr_lev];
  return
    (std::find(pop_mi_l.begin(), pop_mi_l.end(), tr_set) != pop_mi_l.end());
}


/** identify if newly-pushed trial set exists within stored data sets */
inline bool HierarchSparseGridDriver::
push_trial_available(const UShortArray& key)
{
  const UShortArray& tr_set = trial_set(key);
  size_t tr_lev = l1_norm(tr_set);
  const UShortArrayDeque& pop_mi_l = poppedLevMultiIndex[key][tr_lev];
  return
    (std::find(pop_mi_l.begin(), pop_mi_l.end(), tr_set) != pop_mi_l.end());
}


/** identify if newly-pushed trial set exists within stored data sets */
inline bool HierarchSparseGridDriver::push_trial_available()
{ return push_trial_available(activeKey, trial_set()); }


/** identify where newly-pushed trial set exists within stored data sets */
inline size_t HierarchSparseGridDriver::
push_trial_index(const UShortArray& key, const UShortArray& tr_set)
{
  size_t tr_lev = l1_norm(tr_set);
  return find_index(poppedLevMultiIndex[key][tr_lev], tr_set);
}


/** identify where newly-pushed trial set exists within stored data sets */
inline size_t HierarchSparseGridDriver::push_trial_index(const UShortArray& key)
{
  const UShortArray& tr_set = trial_set(key);
  size_t tr_lev = l1_norm(tr_set);
  return find_index(poppedLevMultiIndex[key][tr_lev], tr_set);
}


/** identify where newly-pushed trial set exists within stored data sets */
inline size_t HierarchSparseGridDriver::push_trial_index()
{ return push_trial_index(activeKey, trial_set()); }


inline size_t HierarchSparseGridDriver::push_index() const
{ return pushIndex; }


inline const UShortArray& HierarchSparseGridDriver::increment_sets() const
{ return incrSetsIter->second; }


/*
inline size_t HierarchSparseGridDriver::
popped_sets(const UShortArray& key) const
{
  // Avoid double lookup since T2 cannot currently exist w/o T1
  //return std::max(poppedT1WtSets[key].size(), poppedT2WtSets[key].size());

  std::map<UShortArray, RealVectorDequeArray>::const_iterator cit
    = poppedT1WtSets.find(key);
  if (cit == poppedT1WtSets.end())
    return 0;
  else {
    const RealVectorDequeArray& pop_t1w = cit->second;
    size_t lev, num_lev = pop_t1w.size(), num_sets = 0;
    for (lev=0; lev<num_lev; ++lev)
      num_sets += pop_t1w[lev].size();
    return num_sets;
  }
}


inline size_t HierarchSparseGridDriver::popped_sets() const
{ return popped_sets(activeKey); }
*/


inline void HierarchSparseGridDriver::print_smolyak_multi_index() const
{
  const UShort3DArray& sm_mi = smolMIIter->second;
  size_t i, j, k, num_lev = sm_mi.size(), cntr = 1;
  for (i=0; i<num_lev; ++i) {
    const UShort2DArray& sm_mi_i = sm_mi[i];
    size_t num_sets = sm_mi_i.size();
    for (j=0; j<num_sets; ++j, ++cntr) {
      PCout << "Smolyak index set " << cntr << ':';
      print_index_set(PCout, sm_mi_i[j]);
    }
  }
}


inline const UShort3DArray& HierarchSparseGridDriver::
smolyak_multi_index() const
{ return smolMIIter->second; }


inline void HierarchSparseGridDriver::
smolyak_multi_index(const UShort3DArray& sm_mi)
{ smolMIIter->second = sm_mi; }


inline const UShort3DArray& HierarchSparseGridDriver::
smolyak_multi_index(const UShortArray& key) const
{
  std::map<UShortArray, UShort3DArray>::const_iterator cit
    = smolyakMultiIndex.find(key);
  if (cit == smolyakMultiIndex.end()) {
    PCerr << "Error: key not found in HierarchSparseGridDriver::"
	  << "smolyak_multi_index()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const std::map<UShortArray, UShort3DArray>& HierarchSparseGridDriver::
smolyak_multi_index_map() const
{ return smolyakMultiIndex; }


inline void HierarchSparseGridDriver::
track_collocation_indices(bool track_colloc_indices)
{ trackCollocIndices = track_colloc_indices; }


inline bool HierarchSparseGridDriver::track_collocation_indices() const
{ return trackCollocIndices; }


inline const UShort4DArray& HierarchSparseGridDriver::collocation_key() const
{ return collocKeyIter->second; }


inline void HierarchSparseGridDriver::collocation_key(const UShort4DArray& key)
{ collocKeyIter->second = key; }


inline const UShort4DArray& HierarchSparseGridDriver::
collocation_key(const UShortArray& key) const
{
  std::map<UShortArray, UShort4DArray>::const_iterator cit
    = collocKey.find(key);
  if (cit == collocKey.end()) {
    PCerr << "Error: key not found in HierarchSparseGridDriver::"
	  << "collocation_key()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const std::map<UShortArray, UShort4DArray>& HierarchSparseGridDriver::
collocation_key_map() const
{ return collocKey; }


inline const Sizet3DArray& HierarchSparseGridDriver::collocation_indices() const
{ return collocIndIter->second; }


inline void HierarchSparseGridDriver::
collocation_indices(const Sizet3DArray& indices)
{ collocIndIter->second = indices; }


inline const Sizet3DArray& HierarchSparseGridDriver::
collocation_indices(const UShortArray& key) const
{
  std::map<UShortArray, Sizet3DArray>::const_iterator cit
    = collocIndices.find(key);
  if (cit == collocIndices.end()) {
    PCerr << "Error: key not found in HierarchSparseGridDriver::"
	  << "collocation_indices()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const RealVector2DArray& HierarchSparseGridDriver::
type1_hierarchical_weight_sets() const
{ return t1WtIter->second; }


inline void HierarchSparseGridDriver::
type1_hierarchical_weight_sets(const RealVector2DArray& rv2)
{ t1WtIter->second = rv2; }


inline const RealMatrix2DArray& HierarchSparseGridDriver::
type2_hierarchical_weight_sets() const
{ return t2WtIter->second; }


inline void HierarchSparseGridDriver::
type2_hierarchical_weight_sets(const RealMatrix2DArray& rm2)
{ t2WtIter->second = rm2; }


inline const std::map<UShortArray, RealVector2DArray>&
HierarchSparseGridDriver::type1_weight_sets_map() const
{ return type1WeightSets; }


inline const std::map<UShortArray, RealMatrix2DArray>&
HierarchSparseGridDriver::type2_weight_sets_map() const
{ return type2WeightSets; }


inline void HierarchSparseGridDriver::
levels_to_delta_sizes(const UShortArray& levels, UShortArray& delta_sizes)
{
  size_t i, num_lev = levels.size();
  if (delta_sizes.size() != num_lev)
    delta_sizes.resize(num_lev);
  for (i=0; i<num_lev; ++i)
    delta_sizes[i] = level_to_delta_size(i, levels[i]);
}


inline void HierarchSparseGridDriver::
levels_to_delta_keys(const UShortArray& levels, UShort2DArray& delta_keys)
{
  size_t i, num_lev = levels.size();
  if (delta_keys.size() != num_lev)
    delta_keys.resize(num_lev);
  for (i=0; i<num_lev; ++i)
    level_to_delta_key(i, levels[i], delta_keys[i]);
}

} // namespace Pecos

#endif

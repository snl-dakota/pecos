/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 SparseGridDriver
//- Description: Wrapper class for C++ code from packages/quadrature/sparse_grid
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef INCREMENTAL_SPARSE_GRID_DRIVER_HPP
#define INCREMENTAL_SPARSE_GRID_DRIVER_HPP

#include "CombinedSparseGridDriver.hpp"


namespace Pecos {

/// Derived integration driver class that generates N-dimensional
/// Smolyak sparse grids for numerical evaluation of expectation
/// integrals over independent standard random variables.

/** This class is used by Dakota::NonDSparseGrid, but could also be
    used for general numerical integration of moments.  It employs 1-D
    Clenshaw-Curtis, Newton-Cotes, and Gaussian quadrature rules
    within Smolyak sparse grids. */

class IncrementalSparseGridDriver: public CombinedSparseGridDriver
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  IncrementalSparseGridDriver();
  /// constructor
  IncrementalSparseGridDriver(unsigned short ssg_level,
			      const RealVector& dim_pref = RealVector(),
			      short growth_rate = MODERATE_RESTRICTED_GROWTH,
			      short refine_control = NO_CONTROL);
  /// destructor
  ~IncrementalSparseGridDriver();

  //
  //- Heading: Virtual function redefinitions
  //

  /// update {smolMI,smolCoeffs,collocKey,collocInd}Iter from activeKey
  void update_active_iterators();

  void compute_grid(RealMatrix& variable_sets);
  /// compute (if updateGridSize) and return number of collocation
  /// points with duplicates removed
  int grid_size();

  /// initialize all sparse grid settings (distribution params already
  /// set within poly_basis)
  void initialize_grid(const std::vector<BasisPolynomial>& poly_basis);

  void clear_inactive();
  void clear_keys();

  void initialize_sets();
  void push_trial_set(const UShortArray& set);
  void restore_set();
  void compute_trial_grid(RealMatrix& var_sets);
  void pop_trial_set();
  void finalize_sets(bool output_sets, bool converged_within_tol);

  void compute_increment(RealMatrix& var_sets);
  void push_increment();
  void merge_increment();
  
  /// return smolyakCoeffsRef
  const IntArray& smolyak_coefficients_reference() const;
  /// update smolyakCoeffsRef and type{1,2}WeightSetsRef for use within the
  /// adaptive grid refinement procedures
  void update_reference();

  /// return trialSet
  const UShortArray& trial_set() const;
  /// return num_unique2
  int unique_trial_points() const;

  //
  //- Heading: Member functions
  //

  /// initialize all sparse grid settings except for distribution params
  void initialize_grid(unsigned short ssg_level, const RealVector& dim_pref,
    const ShortArray& u_types, const ExpansionConfigOptions& ec_options,
    BasisConfigOptions& bc_options,
    short growth_rate = MODERATE_RESTRICTED_GROWTH,
    bool track_uniq_prod_wts = true);

  /// update smolyakMultiIndex and smolyakCoeffs
  void update_smolyak_arrays();
  /// overloaded form updates smolyakCoeffs from smolyakMultiIndex
  void update_smolyak_coefficients(size_t start_index);
  /// update collocKey for the trailing index sets within smolyakMultiIndex
  void update_collocation_key();

  /// define a1{Points,Type1Weights,Type2Weights} based on the reference grid
  void reference_unique(RealMatrix& var_sets);
  /// define a2Points and update collocIndices and uniqueIndexMapping
  /// for trailing index sets within smolyakMultiIndex
  void increment_unique(size_t start_index, bool update_1d_pts_wts = true);
  /// update a1Points by merging with unique a2Points
  void merge_unique();
  // apply all remaining trial sets
  //void finalize_unique(size_t start_index);

private:

  //
  //- Heading: Convenience functions
  //

  /// overloaded form updates sm_mi from sm_coeffs
  void update_smolyak_arrays(UShort2DArray& sm_mi, IntArray& sm_coeffs);
  /// increment sm_{mi,coeffs} to sync with new_sm_{mi,coeffs}
  void increment_smolyak_arrays(const UShort2DArray& new_sm_mi,
				const IntArray& new_sm_coeffs,
				UShort2DArray& sm_mi, IntArray& sm_coeffs);
  /// decrement sm_{mi,coeffs} to sync with new_sm_{mi,coeffs}
  void decrement_smolyak_arrays(const UShort2DArray& new_sm_mi,
				const IntArray& new_sm_coeffs,
				UShort2DArray& sm_mi, IntArray& sm_coeffs);
  /// update the coeffs array based on new trailing index sets within
  /// multi_index for incrementally generated generalized sparse grids
  void update_smolyak_coefficients(size_t start_index,
				   const UShort2DArray& sm_mi,
				   IntArray& sm_coeffs);

  /// aggregate point and weight sets across one or more tensor products
  void compute_tensor_points_weights(size_t start_index, size_t num_indices,
				     bool update_1d_pts_wts, RealMatrix& pts,
				     RealVector& t1_wts, RealMatrix& t2_wts);

  /// convenience function for updating sparse points from a set of
  /// aggregated tensor points
  void update_sparse_points(size_t start_index, const BitArray& is_unique,
			    int index_offset, const RealMatrix& tensor_pts,
			    RealMatrix& unique_pts);
  /// convenience function for updating sparse weights from a set of
  /// aggregated tensor weights
  void update_sparse_weights(size_t start_index,
			     const RealVector& tensor_t1_wts,
			     const RealMatrix& tensor_t2_wts,
			     RealVector& updated_t1_wts,
			     RealMatrix& updated_t2_wts);

  /// define the reference collocation indices
  void assign_collocation_indices();
  /// define an increment to the collocation indices
  void update_collocation_indices(size_t start_index);

  //
  //- Heading: Data
  //

  /// trial evaluation set from push_trial_set()
  UShortArray trialSet;

  /// reference values for the Smolyak combinatorial coefficients;
  /// used in incremental approaches that update smolyakCoeffs
  std::map<UShortArray, IntArray> smolyakCoeffsRef;
  /// reference values for the type1 weights corresponding to the current
  /// reference grid; used in incremental approaches that update type1WeightSets
  std::map<UShortArray, RealVector> type1WeightSetsRef;
  /// reference values for the type2 weights corresponding to the current
  /// reference grid; used in incremental approaches that update type2WeightSets
  std::map<UShortArray, RealMatrix> type2WeightSetsRef;

  /// number of unique points in set 1 (reference)
  std::map<UShortArray, int> numUnique1;
  /// active entry within numUnique1
  std::map<UShortArray, int>::iterator numUniq1Iter;
  /// number of unique points in set 2 (increment)
  std::map<UShortArray, int> numUnique2;
  /// active entry within numUnique2
  std::map<UShortArray, int>::iterator numUniq2Iter;

  /// random vector used within sgmgg for sorting
  std::map<UShortArray, RealVector> zVec;
  /// distance values for sorting in set 1 (reference)
  std::map<UShortArray, RealVector> r1Vec;
  /// distance values for sorting in set 2 (increment)
  std::map<UShortArray, RealVector> r2Vec;

  /// array of collocation points in set 1 (reference)
  std::map<UShortArray, RealMatrix> a1Points;
  /// active entry within a1Points
  std::map<UShortArray, RealMatrix>::iterator a1PIter;
  /// vector of type1 weights in set 1 (reference)
  std::map<UShortArray, RealVector> a1Type1Weights;
  /// active entry within a1Type1Weights
  std::map<UShortArray, RealVector>::iterator a1T1WIter;
  /// matrix of type2 weights in set 1 (reference)
  std::map<UShortArray, RealMatrix> a1Type2Weights;
  /// active entry within a1Type2Weights
  std::map<UShortArray, RealMatrix>::iterator a1T2WIter;

  /// array of collocation points in set 2 (increment)
  std::map<UShortArray, RealMatrix> a2Points;
  /// active entry within a2Points
  std::map<UShortArray, RealMatrix>::iterator a2PIter;
  /// vector of type1 weights in set 2 (increment)
  std::map<UShortArray, RealVector> a2Type1Weights;
  /// active entry within a2Type1Weights
  std::map<UShortArray, RealVector>::iterator a2T1WIter;
  /// matrix of type2 weights in set 2 (increment)
  std::map<UShortArray, RealMatrix> a2Type2Weights;
  /// active entry within a2Type2Weights
  std::map<UShortArray, RealMatrix>::iterator a2T2WIter;

  /// ascending sort index for set 1 (reference)
  std::map<UShortArray, IntArray> sortIndex1;
  /// ascending sort index for set 2 (increment)
  std::map<UShortArray, IntArray> sortIndex2;

  /// index within a1 (reference set) of unique points
  std::map<UShortArray, IntArray> uniqueSet1;
  /// active entry within uniqueSet1
  std::map<UShortArray, IntArray>::iterator uniqSet1Iter;
  /// index within a2 (increment set) of unique points
  std::map<UShortArray, IntArray> uniqueSet2;
  /// active entry within uniqueSet2
  std::map<UShortArray, IntArray>::iterator uniqSet2Iter;

  /// index within uniqueSet1 corresponding to all of a1
  std::map<UShortArray, IntArray> uniqueIndex1;
  /// active entry within uniqueIndex1
  std::map<UShortArray, IntArray>::iterator uniqInd1Iter;
  /// index within uniqueSet2 corresponding to all of a2
  std::map<UShortArray, IntArray> uniqueIndex2;
  /// active entry within uniqueIndex2
  std::map<UShortArray, IntArray>::iterator uniqInd2Iter;

  /// key to unique points in set 1 (reference)
  std::map<UShortArray, BitArray> isUnique1;
  /// active entry within isUnique1
  std::map<UShortArray, BitArray>::iterator isUniq1Iter;
  /// key to unique points in set 2 (increment)
  std::map<UShortArray, BitArray> isUnique2;
  /// active entry within isUnique2
  std::map<UShortArray, BitArray>::iterator isUniq2Iter;

  /// unrolled array of collocation indices
  std::map<UShortArray, IntArray> uniqueIndexMapping;
  /// active entry within uniqueIndexMapping
  std::map<UShortArray, IntArray>::iterator uniqIndMapIter;
};


inline IncrementalSparseGridDriver::IncrementalSparseGridDriver():
  CombinedSparseGridDriver()
{ update_active_iterators(); }


inline IncrementalSparseGridDriver::
IncrementalSparseGridDriver(unsigned short ssg_level,
			    const RealVector& dim_pref, short growth_rate,
			    short refine_control):
  CombinedSparseGridDriver(ssg_level, dim_pref, growth_rate, refine_control)
{ update_active_iterators(); }


inline IncrementalSparseGridDriver::~IncrementalSparseGridDriver()
{ }


inline void IncrementalSparseGridDriver::update_active_iterators()
{
  CombinedSparseGridDriver::update_active_iterators();

  numUniq1Iter = numUnique1.find(activeKey);
  if (numUniq1Iter == numUnique1.end()) {
    std::pair<UShortArray, int> ua_pair(activeKey, 0);
    numUniq1Iter = numUnique1.insert(ua_pair).first;
  }
  numUniq2Iter = numUnique2.find(activeKey);
  if (numUniq2Iter == numUnique2.end()) {
    std::pair<UShortArray, int> ua_pair(activeKey, 0);
    numUniq2Iter = numUnique2.insert(ua_pair).first;
  }
  a1PIter = a1Points.find(activeKey);
  if (a1PIter == a1Points.end()) {
    std::pair<UShortArray, RealMatrix> ua_pair(activeKey, RealMatrix());
    a1PIter = a1Points.insert(ua_pair).first;
  }
  a1T1WIter = a1Type1Weights.find(activeKey);
  if (a1T1WIter == a1Type1Weights.end()) {
    std::pair<UShortArray, RealVector> ua_pair(activeKey, RealVector());
    a1T1WIter = a1Type1Weights.insert(ua_pair).first;
  }
  a1T2WIter = a1Type2Weights.find(activeKey);
  if (a1T2WIter == a1Type2Weights.end()) {
    std::pair<UShortArray, RealMatrix> ua_pair(activeKey, RealMatrix());
    a1T2WIter = a1Type2Weights.insert(ua_pair).first;
  }
  a2PIter = a2Points.find(activeKey);
  if (a2PIter == a2Points.end()) {
    std::pair<UShortArray, RealMatrix> ua_pair(activeKey, RealMatrix());
    a2PIter = a2Points.insert(ua_pair).first;
  }
  a2T1WIter = a2Type1Weights.find(activeKey);
  if (a2T1WIter == a2Type1Weights.end()) {
    std::pair<UShortArray, RealVector> ua_pair(activeKey, RealVector());
    a2T1WIter = a2Type1Weights.insert(ua_pair).first;
  }
  a2T2WIter = a2Type2Weights.find(activeKey);
  if (a2T2WIter == a2Type2Weights.end()) {
    std::pair<UShortArray, RealMatrix> ua_pair(activeKey, RealMatrix());
    a2T2WIter = a2Type2Weights.insert(ua_pair).first;
  }
  uniqSet1Iter = uniqueSet1.find(activeKey);
  if (uniqSet1Iter == uniqueSet1.end()) {
    std::pair<UShortArray, IntArray> ua_pair(activeKey, IntArray());
    uniqSet1Iter = uniqueSet1.insert(ua_pair).first;
  }
  uniqSet2Iter = uniqueSet2.find(activeKey);
  if (uniqSet2Iter == uniqueSet2.end()) {
    std::pair<UShortArray, IntArray> ua_pair(activeKey, IntArray());
    uniqSet2Iter = uniqueSet2.insert(ua_pair).first;
  }
  uniqInd1Iter = uniqueIndex1.find(activeKey);
  if (uniqInd1Iter == uniqueIndex1.end()) {
    std::pair<UShortArray, IntArray> ua_pair(activeKey, IntArray());
    uniqInd1Iter = uniqueIndex1.insert(ua_pair).first;
  }
  uniqInd2Iter = uniqueIndex2.find(activeKey);
  if (uniqInd2Iter == uniqueIndex2.end()) {
    std::pair<UShortArray, IntArray> ua_pair(activeKey, IntArray());
    uniqInd2Iter = uniqueIndex2.insert(ua_pair).first;
  }
  isUniq1Iter = isUnique1.find(activeKey);
  if (isUniq1Iter == isUnique1.end()) {
    std::pair<UShortArray, BitArray> ua_pair(activeKey, BitArray());
    isUniq1Iter = isUnique1.insert(ua_pair).first;
  }
  isUniq2Iter = isUnique2.find(activeKey);
  if (isUniq2Iter == isUnique2.end()) {
    std::pair<UShortArray, BitArray> ua_pair(activeKey, BitArray());
    isUniq2Iter = isUnique2.insert(ua_pair).first;
  }
  uniqIndMapIter = uniqueIndexMapping.find(activeKey);
  if (uniqIndMapIter == uniqueIndexMapping.end()) {
    std::pair<UShortArray, IntArray> ua_pair(activeKey, IntArray());
    uniqIndMapIter = uniqueIndexMapping.insert(ua_pair).first;
  }
}


inline void IncrementalSparseGridDriver::clear_keys()
{
  CombinedSparseGridDriver::clear_keys();

  smolyakCoeffsRef.clear();
  type1WeightSetsRef.clear(); type2WeightSetsRef.clear();

  zVec.clear();            r1Vec.clear();          r2Vec.clear();
  sortIndex1.clear();      sortIndex2.clear();

  numUnique1.clear();      numUniq1Iter = numUnique1.end();
  numUnique2.clear();      numUniq2Iter = numUnique2.end();
  a1Points.clear();        a1PIter      = a1Points.end();
  a1Type1Weights.clear();  a1T1WIter    = a1Type1Weights.end();
  a1Type2Weights.clear();  a1T2WIter    = a1Type2Weights.end();
  a2Points.clear();        a2PIter      = a2Points.end();
  a2Type1Weights.clear();  a2T1WIter    = a2Type1Weights.end();
  a2Type2Weights.clear();  a2T2WIter    = a2Type2Weights.end();
  uniqueSet1.clear();      uniqSet1Iter = uniqueSet1.end();
  uniqueSet2.clear();      uniqSet2Iter = uniqueSet2.end();
  uniqueIndex1.clear();    uniqInd1Iter = uniqueIndex1.end();
  uniqueIndex2.clear();    uniqInd2Iter = uniqueIndex2.end();
  isUnique1.clear();       isUniq1Iter  = isUnique1.end();
  isUnique2.clear();       isUniq2Iter  = isUnique2.end();

  uniqueIndexMapping.clear();  uniqIndMapIter = uniqueIndexMapping.end();
}


inline const UShortArray& IncrementalSparseGridDriver::trial_set() const
{ return trialSet; }


inline void IncrementalSparseGridDriver::update_reference()
{
  smolyakCoeffsRef[activeKey] = smolCoeffsIter->second;
  if (trackUniqueProdWeights) {
    type1WeightSetsRef[activeKey] = type1WeightSets[activeKey];
    if (computeType2Weights)
      type2WeightSetsRef[activeKey] = type2WeightSets[activeKey];
  }
}


inline const IntArray& IncrementalSparseGridDriver::
smolyak_coefficients_reference() const
{
  std::map<UShortArray, IntArray>::const_iterator cit
    = smolyakCoeffsRef.find(activeKey);
  if (cit == smolyakCoeffsRef.end()) {
    PCerr << "Error: active key not found in CombinedSparseGridDriver::"
	  << "smolyak_coefficients_reference()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline int IncrementalSparseGridDriver::unique_trial_points() const
{ return numUniq2Iter->second; }


/** Start from scratch rather than incur incremental coefficient update. */
inline void IncrementalSparseGridDriver::update_smolyak_arrays()
{ update_smolyak_arrays(smolMIIter->second, smolCoeffsIter->second); }


inline void IncrementalSparseGridDriver::
update_smolyak_coefficients(size_t start_index)
{
  update_smolyak_coefficients(start_index, smolMIIter->second,
			      smolCoeffsIter->second);
}


inline void IncrementalSparseGridDriver::merge_increment()
{ merge_unique(); } // form a3 and promote to new a1

} // namespace Pecos

#endif

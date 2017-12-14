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

#ifndef COMBINED_SPARSE_GRID_DRIVER_HPP
#define COMBINED_SPARSE_GRID_DRIVER_HPP

#include "SparseGridDriver.hpp"

namespace Pecos {

/// pointer to a collocation point or weight evaluation function, matching
/// the GWPointer prototype required by Pecos/packages/VPISparseGrid
typedef void ( *CollocFnPtr ) ( int order, int index, double* data );
/// pointer to a level-growth-to-order mapping function, matching the
/// GWPointer2 prototype required by Pecos/packages/VPISparseGrid
typedef int ( *LevGrwOrdFnPtr ) ( int level, int growth );


/// Derived integration driver class that generates N-dimensional
/// Smolyak sparse grids for numerical evaluation of expectation
/// integrals over independent standard random variables.

/** This class is used by Dakota::NonDSparseGrid, but could also be
    used for general numerical integration of moments.  It employs 1-D
    Clenshaw-Curtis, Newton-Cotes, and Gaussian quadrature rules
    within Smolyak sparse grids. */

class CombinedSparseGridDriver: public SparseGridDriver
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  CombinedSparseGridDriver();
  /// constructor
  CombinedSparseGridDriver(unsigned short ssg_level,
			   const RealVector& dim_pref = RealVector(),
			   short growth_rate = MODERATE_RESTRICTED_GROWTH,
			   short refine_control = NO_CONTROL);
  /// destructor
  ~CombinedSparseGridDriver();

  //
  //- Heading: Virtual function redefinitions
  //

  void compute_grid(RealMatrix& variable_sets);
  /// compute (if updateGridSize) and return number of collocation
  /// points with duplicates removed
  int grid_size();
  void reinterpolated_tensor_grid(const UShortArray& lev_index,
				  const SizetList& reinterp_indices);

  /// initialize all sparse grid settings (distribution params already
  /// set within poly_basis)
  void initialize_grid(const std::vector<BasisPolynomial>& poly_basis);

  /*
  void store_grid(size_t index = _NPOS);
  void restore_grid(size_t index = _NPOS);
  void remove_stored_grid(size_t index = _NPOS);
  void clear_stored();

  void swap_grid(size_t index);
  */

  size_t maximal_grid() const;

  void initialize_sets();
  void push_trial_set(const UShortArray& set);
  void restore_set();
  void compute_trial_grid(RealMatrix& var_sets);
  void pop_trial_set();
  void merge_set();
  void finalize_sets(bool output_sets, bool converged_within_tol);

  /// update smolyakCoeffsRef and type{1,2}WeightSetsRef for use within the
  /// generalized sparse grid procedure
  void update_reference();

  /// return trialSet
  const UShortArray& trial_set() const;
  /// return num_unique2
  int unique_trial_points() const;

  void compute_grid_increment(RealMatrix& var_sets);

  void print_smolyak_multi_index() const;

  //
  //- Heading: Member functions
  //

  /// initialize all sparse grid settings except for distribution params
  void initialize_grid(unsigned short ssg_level, const RealVector& dim_pref,
    const ShortArray& u_types, const ExpansionConfigOptions& ec_options,
    BasisConfigOptions& bc_options,
    short growth_rate = MODERATE_RESTRICTED_GROWTH, bool track_colloc = false,
    bool track_uniq_prod_wts = true);

  /// overloaded form initializes smolyakMultiIndex and smolyakCoeffs
  void assign_smolyak_arrays();
  /// initialize Smolyak multi-index (index sets defining the set of tensor
  /// products) and Smolyak combinatorial coefficients using an isotropic or
  /// anisotropic index set constraint.  For anisotropic, webbur::sgmga_vcn_*
  /// functions are used to compute index sets satisfying the anisotropic
  /// index set constraint, along with their corresponding coefficients.
  void assign_smolyak_arrays(UShort2DArray& multi_index, IntArray& coeffs);
  /// overloaded form updates smolyakCoeffs from smolyakMultiIndex
  void update_smolyak_coefficients(size_t start_index);
  /// update the coeffs array based on new trailing index sets within
  /// multi_index for incrementally generated generalized sparse grids
  void update_smolyak_coefficients(size_t start_index,
				   const UShort2DArray& multi_index,
				   IntArray& coeffs);

  /// initialize collocKey from smolyakMultiIndex
  void assign_collocation_key();
  /// update collocKey for the trailing index sets within smolyakMultiIndex
  void update_collocation_key(size_t start_index);
  /// initialize collocIndices from collocKey and uniqueIndexMapping
  void assign_collocation_indices();

  /// define a1{Points,Type1Weights,Type2Weights} based on the reference grid
  void reference_unique(RealMatrix& var_sets);
  /// define a2Points and update collocIndices and uniqueIndexMapping
  /// for the trailing index set within smolyakMultiIndex
  void increment_unique(bool compute_a2, bool update_sets,
			RealMatrix& var_sets);
  /// update a1Points by merging with unique a2Points
  void merge_unique();
  /// apply all remaining trial sets
  void finalize_unique(size_t start_index);

  void compute_tensor_points_weights(size_t start_index,
    size_t num_indices, RealMatrix& pts, RealVector& t1_wts,
    RealMatrix& t2_wts);
  void assign_tensor_collocation_indices(size_t start_index,
					 const IntArray& unique_index);

  /// return smolyakMultiIndex[activeKey]
  const UShort2DArray& smolyak_multi_index() const;
  /// return smolyakMultiIndex[key]
  const UShort2DArray& smolyak_multi_index(const UShortArray& key) const;
  /// return smolyakCoeffs[activeKey]
  const IntArray& smolyak_coefficients() const;
  /// return smolyakCoeffs[key]
  const IntArray& smolyak_coefficients(const UShortArray& key) const;
  /// return smolyakCoeffsRef
  const IntArray& smolyak_coefficients_reference() const;

  /// set trackCollocDetails
  void track_collocation_details(bool track_colloc);
  /// get trackCollocDetails
  bool track_collocation_details() const;
  /// set trackUniqueProdWeights
  void track_unique_product_weights(bool track_uniq_prod_wts);
  /// get trackUniqueProdWeights
  bool track_unique_product_weights() const;

  /// return collocKey[activeKey]
  const UShort3DArray& collocation_key() const;
  /// return collocKey[key]
  const UShort3DArray& collocation_key(const UShortArray& key) const;
  /// return collocIndices[activeKey]
  const Sizet2DArray& collocation_indices() const;
  /// return collocIndices[key]
  const Sizet2DArray& collocation_indices(const UShortArray& key) const;
  /// return uniqueIndexMapping
  const IntArray& unique_index_mapping() const;
  // return duplicateTol
  //Real duplicate_tolerance() const;

  /// return type1WeightSets
  const RealVector& type1_weight_sets() const;
  /// return type2WeightSets
  const RealMatrix& type2_weight_sets() const;

private:

  //
  //- Heading: Convenience functions
  //

  /// convenience function for updating sparse points from a set of
  /// aggregated tensor points
  void update_sparse_points(size_t start_index, int new_index_offset,
			    const RealMatrix& tensor_pts,
			    const BitArray& is_unique,
			    const IntArray& unique_index,
			    RealMatrix& new_sparse_pts);
  /// convenience function for updating sparse weights from a set of
  /// aggregated tensor weights
  void update_sparse_weights(size_t start_index,
			     const RealVector& tensor_t1_wts,
			     const RealMatrix& tensor_t2_wts,
			     const IntArray& unique_index,
			     RealVector& updated_t1_wts,
			     RealMatrix& updated_t2_wts);

  /// function passed by pointer for computing collocation points for
  /// polynomialBasis[index]
  static void basis_collocation_points(int order, int index, double* data);
  /// function passed by pointer for computing type 1 collocation
  /// weights for polynomialBasis[index]
  static void basis_type1_collocation_weights(int order,int index,double* data);
  /// function passed by pointer for computing type 2 collocation
  /// weights for polynomialBasis[index]
  static void basis_type2_collocation_weights(int order,int index,double* data);

  /// set duplicateTol based on the content of collocRules: table lookups will
  /// generally be more precise/repeatable than numerically-generated rules
  void initialize_duplicate_tolerance();
  /// initialize compute1D{Points,Type1Weights,Type2Weights} function pointer
  /// arrays for use within webbur::sgmg() and webbur::sgmga() routines
  void initialize_rule_pointers();
  /// initialize levelGrowthToOrder function pointer arrays for use within
  /// webbur::sgmg() and webbur::sgmga() routines
  void initialize_growth_pointers();

  //
  //- Heading: Data
  //

  /// pointer to instance of this class for use in static member functions
  static CombinedSparseGridDriver* sgdInstance;

  /// numSmolyakIndices-by-numVars array for identifying the index to use
  /// within the polynomialBasis for a particular variable
  /** The index sets correspond to j (0-based) for use as indices, which
      are offset from the i indices (1-based) normally used in the Smolyak
      expressions.  The indices correspond to levels, one for each
      anisotropic tensor-product grid within a Smolyak recursion. */
  std::map<UShortArray, UShort2DArray> smolyakMultiIndex;
  /// iterator for active entry within smolyakMultiIndex
  std::map<UShortArray, UShort2DArray>::iterator smolMIIter;

  /// array of Smolyak combinatorial coefficients, one for each tensor
  /// product index set; order is synchronized with smolyakMultiIndex
  std::map<UShortArray, IntArray> smolyakCoeffs;
  /// iterator for active entry within smolyakCoeffs
  std::map<UShortArray, IntArray>::iterator smolCoeffsIter;

  /// reference values for the Smolyak combinatorial coefficients;
  /// used in incremental approaches that update smolyakCoeffs
  IntArray smolyakCoeffsRef;

  /// flag controls conditional population of collocKey, collocIndices,
  /// collocPts1D and type{1,2}CollocWts1D
  bool trackCollocDetails;
  /// flag indicating need to track {type1,type2}WeightSets (product weights for
  /// each unique grid point) as opposed to relying on collections of 1D weights
  bool trackUniqueProdWeights;

  /// numSmolyakIndices-by-numTensorProductPts-by-numVars array for identifying
  /// the 1-D point indices for sets of tensor-product collocation points
  std::map<UShortArray, UShort3DArray> collocKey;
  /// iterator for active entry within collocKey
  std::map<UShortArray, UShort3DArray>::iterator collocKeyIter;

  /// of tensor products to the unique collocation points evaluated
  std::map<UShortArray, Sizet2DArray> collocIndices;
  /// iterator for active entry within collocIndices
  std::map<UShortArray, Sizet2DArray>::iterator collocIndIter;

  // maps indices and bases from sgmga_index() to collocation point index
  //IntArraySizetMap ssgIndexMap;

  /// trial evaluation set from push_trial_set()
  UShortArray trialSet;

  /// the set of type1 weights (for integration of value interpolants)
  /// associated with each point in the sparse grid
  std::map<UShortArray, RealVector> type1WeightSets;
  /// the set of type2 weights (for integration of gradient interpolants)
  /// for each derivative component and for each point in the sparse grid
  std::map<UShortArray, RealMatrix> type2WeightSets;

  /// reference values for the type1 weights corresponding to the current
  /// reference grid; used in incremental approaches that update type1WeightSets
  RealVector type1WeightSetsRef;
  /// reference values for the type2 weights corresponding to the current
  /// reference grid; used in incremental approaches that update type2WeightSets
  RealMatrix type2WeightSetsRef;

  /// array of pointers to collocation point evaluation functions
  std::vector<CollocFnPtr> compute1DPoints;
  /// array of pointers to type1 collocation weight evaluation functions
  std::vector<CollocFnPtr> compute1DType1Weights;
  // 2D array of pointers to type2 collocation weight evaluation functions
  //std::vector<std::vector<CollocFnPtr> > compute1DType2Weights;

  /// array of pointers to webbur::level_to_growth functions
  std::vector<LevGrwOrdFnPtr> levelGrowthToOrder;

  /// output from sgmga_unique_index()
  std::map<UShortArray, IntArray> uniqueIndexMapping;
  /// duplication tolerance used in sgmga routines
  Real duplicateTol;

  int numUnique1;       ///< number of unique points in set 1 (reference)
  int numUnique2;       ///< number of unique points in set 2 (increment)
  RealVector zVec;      ///< random vector used within sgmgg for sorting
  RealVector r1Vec;     ///< distance values for sorting in set 1 (reference)
  RealVector r2Vec;     ///< distance values for sorting in set 2 (increment)
  RealMatrix a1Points;  ///< array of collocation points in set 1 (reference)
  RealMatrix a2Points;  ///< array of collocation points in set 2 (increment)
  RealVector a1Type1Weights; ///< vector of type1 weights in set 1 (reference)
  RealMatrix a1Type2Weights; ///< matrix of type2 weights in set 1 (reference)
  RealVector a2Type1Weights; ///< vector of type1 weights in set 2 (increment)
  RealMatrix a2Type2Weights; ///< matrix of type2 weights in set 2 (increment)
  IntArray sortIndex1;  ///< ascending sort index for set 1 (reference)
  IntArray sortIndex2;  ///< ascending sort index for set 2 (increment)
  IntArray uniqueSet1;  ///< index within a1 (reference set) of unique points
  IntArray uniqueSet2;  ///< index within a2 (increment set) of unique points
  IntArray uniqueIndex1;///< index within uniqueSet1 corresponding to all of a1
  IntArray uniqueIndex2;///< index within uniqueSet2 corresponding to all of a2
  BitArray isUnique1;   ///< key to unique points in set 1 (reference)
  BitArray isUnique2;   ///< key to unique points in set 2 (increment)
};


inline CombinedSparseGridDriver::CombinedSparseGridDriver():
  SparseGridDriver(), trackCollocDetails(false), trackUniqueProdWeights(false),
  duplicateTol(1.e-15)
{ }


inline CombinedSparseGridDriver::
CombinedSparseGridDriver(unsigned short ssg_level, const RealVector& dim_pref,
			 short growth_rate, short refine_control):
  SparseGridDriver(ssg_level, dim_pref, growth_rate, refine_control),
  trackCollocDetails(false), trackUniqueProdWeights(false), duplicateTol(1.e-15)
{ }


inline CombinedSparseGridDriver::~CombinedSparseGridDriver()
{ }


inline void CombinedSparseGridDriver::create_active_iterators()
{
  std::pair<UShortArray, UShort2DArray> u2a_pair(activeKey, UShort2DArray());
  std::pair<UShortArray, IntArray>       ia_pair(activeKey, IntArray());
  std::pair<UShortArray, UShort3DArray> u3a_pair(activeKey, UShort3DArray());
  std::pair<UShortArray, Sizet2DArray>  s2a_pair(activeKey, Sizet2DArray());

  // returned iterator points to existing instance or new insertion
  smolMIIter     = smolyakMultiIndex.insert(u2a_pair).first;
  smolCoeffsIter =     smolyakCoeffs.insert(ia_pair).first;
  collocKeyIter  =         collocKey.insert(u3a_pair).first;
  collocIndIter  =     collocIndices.insert(s2a_pair).first;
}


inline void CombinedSparseGridDriver::update_active_iterators()
{
  smolMIIter     = smolyakMultiIndex.find(activeKey);
  smolCoeffsIter =     smolyakCoeffs.find(activeKey);
  collocKeyIter  =         collocKey.find(activeKey);
  collocIndIter  =     collocIndices.find(activeKey);
}


inline const UShort2DArray& CombinedSparseGridDriver::
smolyak_multi_index() const
{ return smolMIIter->second; }


inline const UShort2DArray& CombinedSparseGridDriver::
smolyak_multi_index(const UShortArray& key) const
{
  std::map<UShortArray, UShort2DArray>::const_iterator cit
    = smolyakMultiIndex.find(key);
  if (cit == smolyakMultiIndex.end()) {
    PCerr << "Error: key not found in CombinedSparseGridDriver::"
	  << "smolyak_multi_index()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const IntArray& CombinedSparseGridDriver::smolyak_coefficients() const
{ return smolCoeffsIter->second; }


inline const IntArray& CombinedSparseGridDriver::
smolyak_coefficients(const UShortArray& key) const
{
  std::map<UShortArray, IntArray>::const_iterator cit = smolyakCoeffs.find(key);
  if (cit == smolyakCoeffs.end()) {
    PCerr << "Error: key not found in CombinedSparseGridDriver::"
	  << "smolyak_coefficients()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline void CombinedSparseGridDriver::
track_collocation_details(bool track_colloc)
{ trackCollocDetails = track_colloc; }


inline bool CombinedSparseGridDriver::track_collocation_details() const
{ return trackCollocDetails; }


inline void CombinedSparseGridDriver::
track_unique_product_weights(bool track_uniq_prod_wts)
{ trackUniqueProdWeights = track_uniq_prod_wts; }


inline bool CombinedSparseGridDriver::track_unique_product_weights() const
{ return trackUniqueProdWeights; }


inline const UShort3DArray& CombinedSparseGridDriver::collocation_key() const
{ return collocKeyIter->second; }


inline const UShort3DArray& CombinedSparseGridDriver::
collocation_key(const UShortArray& key) const
{
  std::map<UShortArray, UShort3DArray>::const_iterator cit
    = collocKey.find(key);
  if (cit == collocKey.end()) {
    PCerr << "Error: key not found in CombinedSparseGridDriver::"
	  << "collocation_key()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const Sizet2DArray& CombinedSparseGridDriver::collocation_indices() const
{ return collocIndIter->second; }


inline const Sizet2DArray& CombinedSparseGridDriver::
collocation_indices(const UShortArray& key) const
{
  std::map<UShortArray, Sizet2DArray>::const_iterator cit
    = collocIndices.find(key);
  if (cit == collocIndices.end()) {
    PCerr << "Error: key not found in CombinedSparseGridDriver::"
	  << "collocation_indices()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const IntArray& CombinedSparseGridDriver::unique_index_mapping() const
{
  std::map<UShortArray, >::const_iterator cit
    = uniqueIndexMapping.find(activeKey);
  if (cit == uniqueIndexMapping.end()) {
    PCerr << "Error: active key not found in CombinedSparseGridDriver::"
	  << "unique_index_mapping()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const UShortArray& CombinedSparseGridDriver::trial_set() const
{ return trialSet; }


inline int CombinedSparseGridDriver::unique_trial_points() const
{ return numUnique2; }


//inline Real CombinedSparseGridDriver::duplicate_tolerance() const
//{ return duplicateTol; }


inline void CombinedSparseGridDriver::print_smolyak_multi_index() const
{
  const UShort2DArray& sm_mi = smolyakMultiIndex[activeKey];
  const IntArray&  sm_coeffs = smolyakCoeffs[activeKey];
  size_t i, sm_mi_len = sm_mi.size(), cntr = 0;
  for (i=0; i<sm_mi_len; ++i) {
    if (sm_coeffs[i]) {
      PCout << "Smolyak index set " << ++cntr << ':';
      print_index_set(PCout, sm_mi[i]);
    }
  }
}


inline void CombinedSparseGridDriver::assign_smolyak_arrays()
{
  assign_smolyak_arrays(smolyakMultiIndex[activeKey],
			smolyakCoeffs[activeKey]);
}


inline void CombinedSparseGridDriver::
update_smolyak_coefficients(size_t start_index)
{
  update_smolyak_coefficients(start_index, smolyakMultiIndex[activeKey],
			      smolyakCoeffs[activeKey]);
}


inline void CombinedSparseGridDriver::merge_set()
{ merge_unique(); }


inline void CombinedSparseGridDriver::update_reference()
{
  smolyakCoeffsRef = smolyakCoeffs[activeKey];
  if (trackUniqueProdWeights) {
    type1WeightSetsRef = type1WeightSets[activeKey];
    if (computeType2Weights)
      type2WeightSetsRef = type2WeightSets[activeKey];
  }
}


inline const IntArray& CombinedSparseGridDriver::
smolyak_coefficients_reference() const
{ return smolyakCoeffsRef; }


inline const RealVector& CombinedSparseGridDriver::type1_weight_sets() const
{
  std::map<UShortArray, RealVector>::const_iterator cit
    = type1WeightSets.find(activeKey);
  if (cit == type1WeightSets.end()) {
    PCerr << "Error: active key not found in CombinedSparseGridDriver::"
	  << "type1_weight_sets()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const RealMatrix& CombinedSparseGridDriver::type2_weight_sets() const
{
  std::map<UShortArray, RealMatrix>::const_iterator cit
    = type2WeightSets.find(activeKey);
  if (cit == type2WeightSets.end()) {
    PCerr << "Error: active key not found in CombinedSparseGridDriver::"
	  << "type2_weight_sets()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline void CombinedSparseGridDriver::
basis_collocation_points(int order, int index, double* data)
{
  const RealArray& colloc_pts
    = sgdInstance->polynomialBasis[index].collocation_points(order);
  std::copy(colloc_pts.begin(), colloc_pts.begin()+order, data);
}


inline void CombinedSparseGridDriver::
basis_type1_collocation_weights(int order, int index, double* data)
{
  const RealArray& colloc_wts
    = sgdInstance->polynomialBasis[index].type1_collocation_weights(order);
  std::copy(colloc_wts.begin(), colloc_wts.begin()+order, data);
}


inline void CombinedSparseGridDriver::
basis_type2_collocation_weights(int order, int index, double* data)
{
  const RealArray& colloc_wts
    = sgdInstance->polynomialBasis[index].type2_collocation_weights(order);
  std::copy(colloc_wts.begin(), colloc_wts.begin()+order, data);
}

} // namespace Pecos

#endif

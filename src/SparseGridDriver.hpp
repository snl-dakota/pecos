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

#ifndef SPARSE_GRID_DRIVER_HPP
#define SPARSE_GRID_DRIVER_HPP

#include "IntegrationDriver.hpp"

namespace Pecos {

/// pointer to a collocation point or weight evaluation function, matching
/// the prototype required by Pecos/packages/VPISparseGrid
typedef void ( *FPType ) ( int order, int index, double* data );


/// Derived nondeterministic class that generates N-dimensional
/// Smolyak sparse grids for numerical evaluation of expectation
/// integrals over independent standard random variables.

/** This class is used by Dakota::NonDSparseGrid, but could also be
    used for general numerical integration of moments.  It employs 1-D
    Clenshaw-Curtis, Newton-Cotes, and Gaussian quadrature rules
    within Smolyak sparse grids. */

class SparseGridDriver: public IntegrationDriver
{
public:

  //
  //- Heading: Constructors and destructor
  //

  SparseGridDriver();  ///< default constructor
  ~SparseGridDriver(); ///< destructor

  //
  //- Heading: Virtual function redefinitions
  //

  /// compute scaled variable and weight sets for the sparse grid
  void compute_grid(RealMatrix& variable_sets);
  /// number of collocation points with duplicates removed
  int grid_size();

  //
  //- Heading: Member functions
  //

  /// overloaded form initializes smolyakMultiIndex and smolyakCoeffs
  void allocate_smolyak_arrays();
  /// initialize Smolyak multi-index (index sets defining the set of tensor
  /// products) and Smolyak combinatorial coefficients using an isotropic or
  /// anisotropic index set constraint.  For anisotropic, webbur::sgmga_vcn_*
  /// functions are used to compute index sets satisfying the anisotropic
  /// index set constraint, along with their corresponding coefficients.
  void allocate_smolyak_arrays(UShort2DArray& multi_index, IntArray& coeffs);
  /// overloaded form updates smolyakCoeffs from smolyakMultiIndex
  void allocate_smolyak_coefficients(size_t start_index);
  /// update the coeffs array based on new trailing index sets within
  /// multi_index for incrementally generated generalized sparse grids
  void allocate_smolyak_coefficients(size_t start_index,
				     const UShort2DArray& multi_index,
				     IntArray& coeffs);
  /// initialize collocKey from smolyakMultiIndex
  void allocate_collocation_key();
  /// update collocKey for the trailing index sets within smolyakMultiIndex
  void update_collocation_key(size_t start_index);
  /// initialize collocIndices from collocKey and uniqueIndexMapping
  void allocate_collocation_indices();
  /// initialize collocPts1D and collocWts1D
  void allocate_1d_collocation_points_weights();
  /// update collocPts1D and collocWts1D from pts_1d and wts_1d
  void update_1d_collocation_points_weights(const UShortArray& trial_set,
					    const Real2DArray& pts_1d,
					    const Real2DArray& wts_1d);

  /// define a1Points/a1Weights based on the reference grid
  void reference_unique();
  /// define a2Points and update collocIndices and uniqueIndexMapping
  /// for the trailing index set within smolyakMultiIndex
  void increment_unique(bool compute_a2 = true);
  /// update a1Points by merging with unique a2Points
  void merge_unique();
  /// apply all remaining trial sets
  void finalize_unique(size_t start_index);

  /// initialize all sparse grid settings except for distribution params
  void initialize_grid(const ShortArray& u_types, unsigned short ssg_level,
    const RealVector& dim_pref, //short refine_type = NO_REFINEMENT,
    short refine_control = NO_CONTROL, bool store_colloc = false,
    bool track_ensemble_wts = true, bool nested_rules = true,
    bool equidistant_rules = true,
    short growth_rate = MODERATE_RESTRICTED_GROWTH,
    short nested_uniform_rule = GAUSS_PATTERSON);
  /// initialize all sparse grid settings (distribution params already
  /// set within poly_basis)
  void initialize_grid(const std::vector<BasisPolynomial>& poly_basis,
    unsigned short ssg_level, const RealVector& dim_pref,
    //short refine_type = NO_REFINEMENT,
    short refine_control = NO_CONTROL, bool store_colloc = false,
    bool track_ensemble_wts = true,
    short growth_rate = MODERATE_RESTRICTED_GROWTH);

  /// update axisLowerBounds
  void update_axis_lower_bounds();

  /// initializes old/active/evaluation sets for use within the 
  /// generalized sparse grid procedure
  void initialize_sets();
  /// update activeMultiIndex from the passed trial set for use within
  /// the generalized sparse grid procedure
  void add_active_neighbors(const UShortArray& set);
  /// update smolyakCoeffsRef and weightSetsRef for use within the
  /// generalized sparse grid procedure
  void update_reference();
  /// update smolyakMultiIndex with a new trial set for use within the
  /// generalized sparse grid procedure
  void push_trial_set(const UShortArray& set);
  /// calls compute_tensor_grid() for the index set from push_trial_set()
  void compute_trial_grid(RealMatrix& unique_variable_sets);
  /// update collocKey, collocIndices, and uniqueIndexMapping based on
  /// restoration of previous trial to smolyakMultiIndex
  void restore_set();
  /// remove the previously pushed trial set from smolyakMultiIndex
  /// during the course of the generalized sparse grid procedure
  void pop_trial_set();
  /// accept the best of several trial sets and update old/active
  /// within the generalized sparse grid procedure
  void update_sets(const UShortArray& set_star);
  /// print the final accepted and trial sets prior to set finalization
  void print_final_sets(bool converged_within_tol);
  /// accept all remaining trial sets within the generalized sparse
  /// grid procedure
  void finalize_sets();

  /// converts an array of sparse grid levels to an array of
  /// quadrature orders based on integrationRules/growthRules
  void level_to_order(size_t i, unsigned short level,
		      unsigned short& order);
  /// converts an array of sparse grid levels to an array of
  /// quadrature orders based on integrationRules/growthRules
  void level_to_order(const UShortArray& levels, UShortArray& orders);

  /// set ssgLevel
  void level(unsigned short ssg_level);
  /// return ssgLevel
  unsigned short level() const;

  /// set anisoLevelWts
  void dimension_preference(const RealVector& dim_pref);
  /// set anisoLevelWts
  void anisotropic_weights(const RealVector& aniso_wts);
  /// return anisoLevelWts
  const RealVector& anisotropic_weights() const;
  /// return dimIsotropic
  bool isotropic() const;

  // return refineType
  //short refine_type()    const;
  /// return refineControl
  short refine_control() const;

  /// return smolyakMultiIndex
  const UShort2DArray& smolyak_multi_index() const;
  /// return smolyakCoeffs
  const IntArray& smolyak_coefficients() const;
  /// return collocKey
  const UShort3DArray& collocation_key() const;
  /// return collocIndices
  const Sizet2DArray& collocation_indices() const;
  /// return uniqueIndexMapping
  const IntArray& unique_index_mapping() const;
  /// return num_unique2
  int unique_trial_points() const;

  // return duplicateTol
  //const Real& duplicate_tolerance() const;

  /// return activeMultiIndex
  const std::set<UShortArray>& active_multi_index() const;
  /// return oldMultiIndex
  const std::set<UShortArray>& old_multi_index() const;
  /// return the trial index set from push_trial_set()
  const UShortArray& trial_index_set() const;
  /// return smolyakCoeffsRef
  const IntArray& smolyak_coefficients_reference() const;

  /// return collocPts1D
  const Real3DArray& collocation_points_array()  const;
  /// return collocWts1D
  const Real3DArray& collocation_weights_array() const;

private:

  //
  //- Heading: Convenience functions
  //

  /// function for computing collocation points for polynomialBasis[index]
  static void basis_collocation_points(int order, int index, double* data);
  /// function for computing collocation weights for polynomialBasis[index]
  static void basis_collocation_weights(int order, int index, double* data);

  /// function for computing collocation points for ChebyshevOrthogPolynomial
  static void chebyshev_points(int order, int index, double* data);
  /// function for computing collocation weights for ChebyshevOrthogPolynomial
  static void chebyshev_weights(int order, int index, double* data);

  /// convenience function for defining {a1,a2}{Points,Weights}
  void compute_tensor_points_weights(size_t start_index, size_t num_indices,
				     RealMatrix& pts, RealVector& wts);
  /// convenience function for updating sparse points/weights from a set of
  /// aggregated tensor points/weights
  void update_sparse_points(size_t start_index, int new_index_offset,
			    const RealMatrix& tensor_pts,
			    const BoolDeque& is_unique,
			    const IntArray& unique_index,
			    RealMatrix& new_sparse_pts);
  /// convenience function for updating sparse points/weights from a set of
  /// aggregated tensor points/weights
  void update_sparse_weights(size_t start_index, const RealVector& tensor_wts,
			     const IntArray& unique_index,
			     RealVector& updated_sparse_wts);
  ///convenience function for assigning collocIndices from uniqueIndex{1,2,3}
  void assign_tensor_collocation_indices(size_t start_index,
					 const IntArray& unique_index);

  /// compute 1-norm |i| (sum of indices) for the given index_set
  unsigned int index_norm(const UShortArray& index_set) const;

  //
  //- Heading: Data
  //

  /// pointer to instance of this class for use in statis member functions
  static SparseGridDriver* sgdInstance;

  /// pointer to a ChebyshevOrthogPolynomial instance for access to Fejer2 and
  /// Clenshaw-Curtis integration points/weights (avoid repeated instantiations)
  BasisPolynomial* chebyPolyPtr;

  /// the Smolyak sparse grid level
  unsigned short ssgLevel;

  /// flag indicating a dimension isotropic grid
  bool dimIsotropic;
  // vector of dimension preference levels for dimension anisotropic grids
  //RealVector dimPref;
  /// weighting vector for dimension anisotropic grids
  RealVector anisoLevelWts;

  // type of expansion refinement
  //short refineType;
  /// algorithm control governing expansion refinement
  short refineControl;
  /// controls conditional population of collocPts1D and collocWts1D
  bool storeCollocDetails;

  /// refinement constraints that ensure that level/anisotropic weight updates
  /// contain all previous multi-index sets
  RealVector axisLowerBounds;

  /// number of parameters for each polynomial; input to sgmga routines
  /// (corresponds to set of variables defined by integrationRules)
  IntArray numPolyParams;
  /// concatenated array of polynomial parameters for input to sgmga routines
  /// (corresponds to set of variables defined by integrationRules)
  RealArray polyParams;
  /// duplication tolerance used in sgmga routines
  Real duplicateTol;

  /// numSmolyakIndices-by-numVars array for identifying the index to use
  /// within the polynomialBasis for a particular variable
  /** The index sets correspond to j (0-based) for use as indices, which
      are offset from the i indices (1-based) normally used in the Smolyak
      expressions.  The indices correspond to levels, one within each
      anisotropic tensor-product integration of a Smolyak recursion. */
  UShort2DArray smolyakMultiIndex;
  /// array of Smolyak combinatorial coefficients, one for each tensor
  /// product index set; order is synchronized with smolyakMultiIndex
  IntArray smolyakCoeffs;
  /// numSmolyakIndices-by-numTensorProductPts-by-numVars array for identifying
  /// the 1-D point indices for sets of tensor-product collocation points
  UShort3DArray collocKey;
  /// numSmolyakIndices-by-numTensorProductPts array for linking the
  /// set of tensor products to the unique collocation points evaluated
  Sizet2DArray collocIndices;
  /// output from sgmga_unique_index()
  IntArray uniqueIndexMapping;
  /// the current number of unique points in the grid; used for incrementing
  /// and decrementing uniqueIndexMapping within generalized sparse grids
  int numCollocPts;
  // maps indices and bases from sgmga_index() to collocation point index
  //IntArraySizetMap ssgIndexMap;

  /// old reference index sets for generalized sparse grids
  std::set<UShortArray> oldMultiIndex; // or UShort2DArray
  /// active index sets under current consideration for inclusion in a
  /// generalized sparse grid
  std::set<UShortArray> activeMultiIndex; // or UShort2DArray
  /// subset of active set that have been evaluated as trial sets
  /// (incremented in compute_trial_grid() and decremented in update_sets())
  std::set<UShortArray> trialSets; // or UShort2DArray
  /// reference values for the Smolyak combinatorial coefficients;
  /// used in incremental approaches that update smolyakCoeffs
  IntArray smolyakCoeffsRef;
  /// reference values for the sparse grid weights corresponding to the current
  /// reference grid; used in incremental approaches that update weightSets
  RealVector weightSetsRef;
  /// flag indicating need to track weightSets for an ensemble sparse grid
  /// computed incrementally
  bool trackEnsembleWeights;

  int numUnique1;       ///< number of unique points in set 1 (reference)
  int numUnique2;       ///< number of unique points in set 2 (increment)
  RealVector zVec;      ///< random vector used within sgmgg for sorting
  RealVector r1Vec;     ///< distance values for sorting in set 1 (reference)
  RealVector r2Vec;     ///< distance values for sorting in set 2 (increment)
  RealMatrix a1Points;  ///< array of collocation points in set 1 (reference)
  RealMatrix a2Points;  ///< array of collocation points in set 2 (increment)
  RealVector a1Weights; ///< array of integration weights in set 1 (reference)
  RealVector a2Weights; ///< array of integration weights in set 2 (increment)
  IntArray sortIndex1;  ///< ascending sort index for set 1 (reference)
  IntArray sortIndex2;  ///< ascending sort index for set 2 (increment)
  IntArray uniqueSet1;  ///< index within a1 (reference set) of unique points
  IntArray uniqueSet2;  ///< index within a2 (increment set) of unique points
  IntArray uniqueIndex1;///< index within uniqueSet1 corresponding to all of a1
  IntArray uniqueIndex2;///< index within uniqueSet2 corresponding to all of a2
  BoolDeque isUnique1;  ///< key to unique points in set 1 (reference)
  BoolDeque isUnique2;  ///< key to unique points in set 2 (increment)

  /// num_levels_per_var x numContinuousVars sets of 1D collocation points
  Real3DArray collocPts1D;
  /// num_levels_per_var x numContinuousVars sets of 1D collocation weights
  Real3DArray collocWts1D;

  /// array of pointers to collocation point evaluation functions
  std::vector<FPType> compute1DPoints;
  /// array of pointers to collocation weight evaluation functions
  std::vector<FPType> compute1DWeights;
};


inline SparseGridDriver::SparseGridDriver():
  IntegrationDriver(BaseConstructor()), chebyPolyPtr(NULL), ssgLevel(0),
  storeCollocDetails(false), duplicateTol(1.e-15)
{ }


inline SparseGridDriver::~SparseGridDriver()
{ if (chebyPolyPtr) delete chebyPolyPtr; }


inline unsigned short SparseGridDriver::level() const
{ return ssgLevel; }


inline void SparseGridDriver::level(unsigned short ssg_level)
{ ssgLevel = ssg_level; }


inline const RealVector& SparseGridDriver::anisotropic_weights() const
{ return anisoLevelWts; }


inline bool SparseGridDriver::isotropic() const
{ return dimIsotropic; }


//inline short SparseGridDriver::refine_type() const
//{ return refineType; }


inline short SparseGridDriver::refine_control() const
{ return refineControl; }


inline const UShort2DArray& SparseGridDriver::smolyak_multi_index() const
{ return smolyakMultiIndex; }


inline const IntArray& SparseGridDriver::smolyak_coefficients() const
{ return smolyakCoeffs; }


inline const UShort3DArray& SparseGridDriver::collocation_key() const
{ return collocKey; }


inline const Sizet2DArray& SparseGridDriver::collocation_indices() const
{ return collocIndices; }


inline const IntArray& SparseGridDriver::unique_index_mapping() const
{ return uniqueIndexMapping; }


inline int SparseGridDriver::unique_trial_points() const
{ return numUnique2; }


//inline const Real& SparseGridDriver::duplicate_tolerance() const
//{ return duplicateTol; }


inline const Real3DArray& SparseGridDriver::collocation_points_array() const
{ return collocPts1D; }


inline const Real3DArray& SparseGridDriver::collocation_weights_array() const
{ return collocWts1D; }


inline void SparseGridDriver::allocate_smolyak_arrays()
{ allocate_smolyak_arrays(smolyakMultiIndex, smolyakCoeffs); }


inline void SparseGridDriver::allocate_smolyak_coefficients(size_t start_index)
{
  allocate_smolyak_coefficients(start_index, smolyakMultiIndex, smolyakCoeffs);
}


inline void SparseGridDriver::update_reference()
{
  smolyakCoeffsRef = smolyakCoeffs;
  if (trackEnsembleWeights)
    weightSetsRef = weightSets;
}


inline const std::set<UShortArray>& SparseGridDriver::active_multi_index() const
{ return activeMultiIndex; }


inline const std::set<UShortArray>& SparseGridDriver::old_multi_index() const
{ return oldMultiIndex; }


inline const UShortArray& SparseGridDriver::trial_index_set() const
{ return smolyakMultiIndex.back(); }


inline const IntArray& SparseGridDriver::smolyak_coefficients_reference() const
{ return smolyakCoeffsRef; }


inline void SparseGridDriver::
level_to_order(const UShortArray& levels, UShortArray& orders)
{
  size_t i, num_levels = levels.size();
  if (orders.size() != num_levels)
    orders.resize(num_levels);
  for (i=0; i<num_levels; i++)
    level_to_order(i, levels[i], orders[i]);
}


inline void SparseGridDriver::
basis_collocation_points(int order, int index, double* data)
{
  const RealArray& colloc_pts
    = sgdInstance->polynomialBasis[index].collocation_points(order);
  std::copy(colloc_pts.begin(), colloc_pts.begin()+order, data);
}


inline void SparseGridDriver::
basis_collocation_weights(int order, int index, double* data)
{
  const RealArray& colloc_wts
    = sgdInstance->polynomialBasis[index].collocation_weights(order);
  std::copy(colloc_wts.begin(), colloc_wts.begin()+order, data);
}


inline void SparseGridDriver::
chebyshev_points(int order, int index, double* data)
{
  sgdInstance->chebyPolyPtr->
    collocation_mode(sgdInstance->integrationRules[index]);
  const RealArray& colloc_pts
    = sgdInstance->chebyPolyPtr->collocation_points(order);
  std::copy(colloc_pts.begin(), colloc_pts.begin()+order, data);
}


inline void SparseGridDriver::
chebyshev_weights(int order, int index, double* data)
{
  sgdInstance->chebyPolyPtr->
    collocation_mode(sgdInstance->integrationRules[index]);
  const RealArray& colloc_wts
    = sgdInstance->chebyPolyPtr->collocation_weights(order);
  std::copy(colloc_wts.begin(), colloc_wts.begin()+order, data);
}


inline unsigned int SparseGridDriver::
index_norm(const UShortArray& index_set) const
{
  unsigned int i, norm = 0, len = index_set.size();
  for (i=0; i<len; ++i)
    norm += index_set[i];
  return norm;
}

} // namespace Pecos

#endif

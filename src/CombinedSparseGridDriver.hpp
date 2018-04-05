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

  /// update {smolMI,smolCoeffs,collocKey,collocInd}Iter from activeKey
  void update_active_iterators();

  void compute_grid(RealMatrix& variable_sets);
  /// compute (if updateGridSize) and return number of collocation
  /// points with duplicates removed
  int grid_size();
  void reinterpolated_tensor_grid(const UShortArray& lev_index,
				  const SizetList& reinterp_indices);

  /// initialize all sparse grid settings (distribution params already
  /// set within poly_basis)
  void initialize_grid(const std::vector<BasisPolynomial>& poly_basis);

  void clear_inactive();
  void clear_keys();

  const UShortArray& maximal_grid() const;

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
  /// initialize collocKey from smolyakMultiIndex
  void assign_collocation_key();
  /// initialize collocIndices from collocKey and unique_index_map
  void assign_collocation_indices(const IntArray& unique_index_map,
				  size_t start_index = 0);

  /// set duplicateTol based on the content of collocRules: table lookups will
  /// generally be more precise/repeatable than numerically-generated rules
  void initialize_duplicate_tolerance();

  /// return smolyakMultiIndex[activeKey]
  const UShort2DArray& smolyak_multi_index() const;
  /// return smolyakMultiIndex[key]
  const UShort2DArray& smolyak_multi_index(const UShortArray& key) const;
  /// return smolyakMultiIndex[key]
  const std::map<UShortArray, UShort2DArray>& smolyak_multi_index_map() const;

  /// return smolyakCoeffs[activeKey]
  const IntArray& smolyak_coefficients() const;
  /// return smolyakCoeffs[key]
  const IntArray& smolyak_coefficients(const UShortArray& key) const;

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

  // return duplicateTol
  //Real duplicate_tolerance() const;

  /// return type1WeightSets
  const RealVector& type1_weight_sets() const;
  /// return type2WeightSets
  const RealMatrix& type2_weight_sets() const;

protected:

  //
  //- Heading: Member functions
  //

  /// initialize Smolyak multi-index (index sets defining the set of tensor
  /// products) and Smolyak combinatorial coefficients using an isotropic or
  /// anisotropic index set constraint.  For anisotropic, webbur::sgmga_vcn_*
  /// functions are used to compute index sets satisfying the anisotropic
  /// index set constraint, along with their corresponding coefficients.
  void assign_smolyak_arrays(UShort2DArray& multi_index, IntArray& coeffs);

  //
  //- Heading: Data
  //

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

  /// the set of type1 weights (for integration of value interpolants)
  /// associated with each unique point in the sparse grid
  std::map<UShortArray, RealVector> type1WeightSets;
  /// the set of type2 weights (for integration of gradient interpolants) for
  /// each derivative component and for each unique point in the sparse grid
  std::map<UShortArray, RealMatrix> type2WeightSets;

  /// duplication tolerance used in sgmga routines
  Real duplicateTol;

private:

  //
  //- Heading: Convenience functions
  //

  /// function passed by pointer for computing collocation points for
  /// polynomialBasis[index]
  static void basis_collocation_points(int order, int index, double* data);
  /// function passed by pointer for computing type 1 collocation
  /// weights for polynomialBasis[index]
  static void basis_type1_collocation_weights(int order,int index,double* data);
  /// function passed by pointer for computing type 2 collocation
  /// weights for polynomialBasis[index]
  static void basis_type2_collocation_weights(int order,int index,double* data);

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

  /// flag controls conditional population of collocKey, collocIndices,
  /// collocPts1D and type{1,2}CollocWts1D
  bool trackCollocDetails;

  /// array of pointers to collocation point evaluation functions
  std::vector<CollocFnPtr> compute1DPoints;
  /// array of pointers to type1 collocation weight evaluation functions
  std::vector<CollocFnPtr> compute1DType1Weights;
  // 2D array of pointers to type2 collocation weight evaluation functions
  //std::vector<std::vector<CollocFnPtr> > compute1DType2Weights;

  /// array of pointers to webbur::level_to_growth functions
  std::vector<LevGrwOrdFnPtr> levelGrowthToOrder;
};


inline CombinedSparseGridDriver::CombinedSparseGridDriver():
  SparseGridDriver(), trackCollocDetails(false), trackUniqueProdWeights(false),
  duplicateTol(1.e-15)
{ update_active_iterators(); }


inline CombinedSparseGridDriver::
CombinedSparseGridDriver(unsigned short ssg_level, const RealVector& dim_pref,
			 short growth_rate, short refine_control):
  SparseGridDriver(ssg_level, dim_pref, growth_rate, refine_control),
  trackCollocDetails(false), trackUniqueProdWeights(false), duplicateTol(1.e-15)
{ update_active_iterators(); }


inline CombinedSparseGridDriver::~CombinedSparseGridDriver()
{ }


inline void CombinedSparseGridDriver::update_active_iterators()
{
  SparseGridDriver::update_active_iterators();

  smolMIIter = smolyakMultiIndex.find(activeKey);
  if (smolMIIter == smolyakMultiIndex.end()) {
    std::pair<UShortArray, UShort2DArray> u2a_pair(activeKey, UShort2DArray());
    smolMIIter = smolyakMultiIndex.insert(u2a_pair).first;
  }
  smolCoeffsIter = smolyakCoeffs.find(activeKey);
  if (smolCoeffsIter == smolyakCoeffs.end()) {
    std::pair<UShortArray, IntArray> ia_pair(activeKey, IntArray());
    smolCoeffsIter = smolyakCoeffs.insert(ia_pair).first;
  }
  collocKeyIter = collocKey.find(activeKey);
  if (collocKeyIter == collocKey.end()) {
    std::pair<UShortArray, UShort3DArray> u3a_pair(activeKey, UShort3DArray());
    collocKeyIter = collocKey.insert(u3a_pair).first;
  }
  collocIndIter = collocIndices.find(activeKey);
  if (collocIndIter == collocIndices.end()) {
    std::pair<UShortArray, Sizet2DArray> s2a_pair(activeKey, Sizet2DArray());
    collocIndIter = collocIndices.insert(s2a_pair).first;
  }
}


inline void CombinedSparseGridDriver::clear_keys()
{
  SparseGridDriver::clear_keys();

  smolyakMultiIndex.clear();  smolMIIter = smolyakMultiIndex.end();
  smolyakCoeffs.clear();  smolCoeffsIter = smolyakCoeffs.end();

  collocKey.clear();       collocKeyIter = collocKey.end();
  collocIndices.clear();   collocIndIter = collocIndices.end();

  type1WeightSets.clear();    type2WeightSets.clear();
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


inline const std::map<UShortArray, UShort2DArray>& CombinedSparseGridDriver::
smolyak_multi_index_map() const
{ return smolyakMultiIndex; }


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


//inline Real CombinedSparseGridDriver::duplicate_tolerance() const
//{ return duplicateTol; }


inline void CombinedSparseGridDriver::print_smolyak_multi_index() const
{
  const UShort2DArray& sm_mi = smolMIIter->second;
  const IntArray&  sm_coeffs = smolCoeffsIter->second;
  size_t i, sm_mi_len = sm_mi.size(), cntr = 0;
  for (i=0; i<sm_mi_len; ++i) {
    if (sm_coeffs[i]) {
      PCout << "Smolyak index set " << ++cntr << ':';
      print_index_set(PCout, sm_mi[i]);
    }
  }
}


inline void CombinedSparseGridDriver::assign_smolyak_arrays()
{ assign_smolyak_arrays(smolMIIter->second, smolCoeffsIter->second); }


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

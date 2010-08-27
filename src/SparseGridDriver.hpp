/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
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
#include "sandia_rules.H"

namespace Pecos {

/// pointer to a Gauss point or weight evaluation function, matching
/// the prototype required by Pecos/packages/VPISparseGrid
typedef void ( *FPType ) ( int order, int index, double* data );


/// Derived nondeterministic class that generates N-dimensional
/// Smolyak sparse grids for numerical evaluation of expectation
/// integrals over independent standard random variables.

/** This class is used by Dakota::NonDSparseGrid, but could also be
    used for general numerical integration of moments.  It employs 1-D
    Clenshaw-Curtis and Gaussian quadrature rules within Smolyak
    sparse grids. */

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
  void compute_grid();
  /// number of collocation points with duplicates removed
  int grid_size();

  /// set ssgAnisoLevelWts
  void dimension_preference(const RealVector& dim_pref);
  /// set ssgAnisoLevelWts
  void anisotropic_weights(const RealVector& aniso_wts);

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
  /// initialize the Smolyak combinatorial coefficients for the
  /// multi-indices defining a generalized sparse grid
  void allocate_generalized_coefficients(const UShort2DArray& multi_index,
					 IntArray& coeffs) const;
  /// initialize collocKey and expansionCoeffIndices
  void allocate_collocation_arrays();

  /// initialize all sparse grid settings except for distribution params
  void initialize_grid(const ShortArray& u_types, unsigned short ssg_level,
		       const RealVector& dim_pref, bool store_1d_gauss = false,
		       bool  nested_rules = true,
		       short growth_rate  = MODERATE_RESTRICTED_GROWTH,
		       short nested_uniform_rule = GAUSS_PATTERSON);
  /// initialize all sparse grid settings (distribution params already
  /// set within poly_basis)
  void initialize_grid(const std::vector<BasisPolynomial>& poly_basis,
		       unsigned short ssg_level, const RealVector& dim_pref,
		       bool  store_1d_gauss = false,
		       short growth_rate = MODERATE_RESTRICTED_GROWTH);

  /// total number of collocation points including duplicates
  int grid_size_total();

  /// update axisLowerBounds
  void update_axis_lower_bounds();

  /// generalized sparse grid function for 
  void initialize_sets();
  /// generalized sparse grid function for 
  void finalize_sets();
  /// generalized sparse grid function for 
  void push_trial_set(const UShortArray& set);
  /// generalized sparse grid function for 
  void pop_trial_set();
  /// generalized sparse grid function for 
  void update_sets(const UShortArray& set_star);
  /// generalized sparse grid function for 
  void add_active_neighbors(const UShortArray& set);

  /// calls compute_tensor_grid() for the index set from push_trial_set()
  void compute_trial_grid();

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

  /// return integrationRules
  const IntArray& integration_rules() const;
  /// return growthRules
  const IntArray& growth_rules() const;

  /// return smolyakMultiIndex
  const UShort2DArray& smolyak_multi_index() const;
  /// return smolyakCoeffs
  const IntArray& smolyak_coefficients() const;
  /// return collocKey
  const UShort3DArray& collocation_key() const;
  /// return expansionCoeffIndices
  const Sizet2DArray& expansion_coefficient_indices() const;
  /// return uniqueIndexMapping
  const IntArray& unique_index_mapping() const;

  // return duplicateTol
  //const Real& duplicate_tolerance() const;

  /// return activeMultiIndex
  const std::set<UShortArray>& active_multi_index() const;
  /// return oldMultiIndex
  const std::set<UShortArray>& old_multi_index() const;
  /// return the trial index set from push_trial_set()
  const UShortArray& trial_index_set() const;

  /// return gaussPts1D
  const Real3DArray& gauss_points_array()  const;
  /// return gaussWts1D
  const Real3DArray& gauss_weights_array() const;

private:

  //
  //- Heading: Convenience functions
  //

  /// function for computing Gauss points for polynomialBasis[index]
  static void basis_gauss_points(int order, int index, double* data);
  /// function for computing Gauss weights for polynomialBasis[index]
  static void basis_gauss_weights(int order, int index, double* data);

  /// function for computing collocation points for ChebyshevOrthogPolynomial
  static void chebyshev_points(int order, int index, double* data);
  /// function for computing collocation weights for ChebyshevOrthogPolynomial
  static void chebyshev_weights(int order, int index, double* data);

  // compute 1-norm |i| for index set i
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

  /// integer codes for integration rule options
  IntArray integrationRules;
  /// integer codes for growth rule options
  IntArray growthRules;

  /// controls conditional population of gaussPts1D and gaussWts1D
  bool store1DGauss;

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
  /// precomputed array of Smolyak combinatorial coefficients
  IntArray smolyakCoeffs;
  /// numSmolyakIndices-by-numTensorProductPts-by-numVars array for identifying
  /// the 1-D point indices for sets of tensor-product collocation points
  UShort3DArray collocKey;
  /// numSmolyakIndices-by-numTensorProductPts array for linking the
  /// set of tensor products to the expansionCoeffs array
  Sizet2DArray expansionCoeffIndices;
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

  /// num_levels_per_var x numContinuousVars sets of 1D Gauss points
  Real3DArray gaussPts1D;
  /// num_levels_per_var x numContinuousVars sets of 1D Gauss weights
  Real3DArray gaussWts1D;

  /// array of pointers to Gauss point evaluation functions
  std::vector<FPType> compute1DPoints;
  /// array of pointers to Gauss weight evaluation functions
  std::vector<FPType> compute1DWeights;
};


inline SparseGridDriver::SparseGridDriver():
  IntegrationDriver(BaseConstructor()), chebyPolyPtr(NULL), ssgLevel(0),
  store1DGauss(false), duplicateTol(1.e-15)
{ }


inline SparseGridDriver::~SparseGridDriver()
{ if (chebyPolyPtr) delete chebyPolyPtr; }


inline unsigned short SparseGridDriver::level() const
{ return ssgLevel; }


inline void SparseGridDriver::level(unsigned short ssg_level)
{ ssgLevel = ssg_level; }


inline const IntArray& SparseGridDriver::integration_rules() const
{ return integrationRules; }


inline const IntArray& SparseGridDriver::growth_rules() const
{ return growthRules; }


inline const UShort2DArray& SparseGridDriver::smolyak_multi_index() const
{ return smolyakMultiIndex; }


inline const IntArray& SparseGridDriver::smolyak_coefficients() const
{ return smolyakCoeffs; }


inline const UShort3DArray& SparseGridDriver::collocation_key() const
{ return collocKey; }


inline const Sizet2DArray& SparseGridDriver::
expansion_coefficient_indices() const
{ return expansionCoeffIndices; }


inline const IntArray& SparseGridDriver::unique_index_mapping() const
{ return uniqueIndexMapping; }


//inline const Real& SparseGridDriver::duplicate_tolerance() const
//{ return duplicateTol; }


inline const Real3DArray& SparseGridDriver::gauss_points_array() const
{ return gaussPts1D; }


inline const Real3DArray& SparseGridDriver::gauss_weights_array() const
{ return gaussWts1D; }


inline void SparseGridDriver::allocate_smolyak_arrays()
{ allocate_smolyak_arrays(smolyakMultiIndex, smolyakCoeffs); }


inline const std::set<UShortArray>& SparseGridDriver::active_multi_index() const
{ return activeMultiIndex; }


inline const std::set<UShortArray>& SparseGridDriver::old_multi_index() const
{ return oldMultiIndex; }


inline const UShortArray& SparseGridDriver::trial_index_set() const
{ return smolyakMultiIndex.back(); }


inline void SparseGridDriver::
level_to_order(size_t i, unsigned short level, unsigned short& order)
{
  int ilevel = level, iorder;
  webbur::level_growth_to_order(1, &ilevel, &integrationRules[i],
				&growthRules[i], &iorder);
  order = iorder;
}


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
basis_gauss_points(int order, int index, double* data)
{
  const RealArray& gauss_pts
    = sgdInstance->polynomialBasis[index].gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


inline void SparseGridDriver::
basis_gauss_weights(int order, int index, double* data)
{
  const RealArray& gauss_wts
    = sgdInstance->polynomialBasis[index].gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


inline void SparseGridDriver::
chebyshev_points(int order, int index, double* data)
{
  sgdInstance->chebyPolyPtr->gauss_mode(sgdInstance->integrationRules[index]);
  const RealArray& gauss_pts = sgdInstance->chebyPolyPtr->gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


inline void SparseGridDriver::
chebyshev_weights(int order, int index, double* data)
{
  sgdInstance->chebyPolyPtr->gauss_mode(sgdInstance->integrationRules[index]);
  const RealArray& gauss_wts = sgdInstance->chebyPolyPtr->gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
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

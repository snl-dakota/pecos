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
			   const RealVector& dim_pref);
  /// destructor
  ~HierarchSparseGridDriver();

  //
  //- Heading: Virtual function redefinitions
  //

  void compute_grid(RealMatrix& var_sets);
  int grid_size();

  void initialize_sets();
  void push_trial_set(const UShortArray& set);
  void restore_set();
  void compute_trial_grid(RealMatrix& var_sets);
  void pop_trial_set();
  void merge_set();
  void finalize_sets();

  int unique_trial_points() const;

  void print_final_sets(bool converged_within_tol) const;
  void print_smolyak_multi_index() const;

  //
  //- Heading: Member functions
  //

  /// return smolyakMultiIndex
  const UShort3DArray& smolyak_multi_index() const;
  /// return collocKey
  const UShort4DArray& collocation_key() const;
  /// return collocIndices
  const Sizet3DArray& collocation_indices() const;

  /// return type1WeightSets
  const RealVector2DArray& type1_weight_set_arrays() const;
  /// return type2WeightSets
  const RealMatrix2DArray& type2_weight_set_arrays() const;

private:

  //
  //- Heading: Convenience functions
  //

  void assign_smolyak_multi_index();
  void assign_collocation_key();
  void update_collocation_key();
  void assign_collocation_indices();
  void update_collocation_indices();

  /// kernel routine used for trialSet and full sparse grid computations
  void compute_points_weights(RealMatrix& pts, RealVector& t1_wts,
			      RealMatrix& t2_wts, const UShortArray& sm_index,
			      const UShort2DArray& colloc_key);
  /// compute points and weights for trialSet
  void compute_points_weights(RealMatrix& pts, RealVector& t1_wts,
			      RealMatrix& t2_wts);
  /// compute points and weights for all levels of the (initial) sparse grid
  void compute_points_weights(RealMatrix& pts, RealVector2DArray& t1_wts,
			      RealMatrix2DArray& t2_wts);

  void level_to_delta_order(const UShortArray& levels,
			    UShort2DArray& delta_quad);

  //
  //- Heading: Data
  //

  /// flag for use of fully nested 1D rules, allowing formulation using
  /// collocation point increments
  bool nestedGrid;

  /// interpolation depth by index set by numVars array for identifying
  /// the index to use within the polynomialBasis for a particular variable
  /** The index sets correspond to j (0-based) for use as indices, which
      are offset from the i indices (1-based) normally used in the Smolyak
      expressions.  The indices correspond to levels, one within each
      anisotropic tensor-product integration of a Smolyak recursion. */
  UShort3DArray smolyakMultiIndex;
  /// levels-by-index sets-by-numDeltaPts-by-numVars array for identifying
  /// the 1-D point indices for sets of tensor-product collocation points
  UShort4DArray collocKey;
  /// levels-by-index sets-by-numTensorProductPts array for linking the
  /// set of tensor products to the unique collocation points evaluated
  Sizet3DArray collocIndices;

  /// the set of type1 weights (for integration of value interpolants)
  /// associated with each point in the sparse grid
  RealVector2DArray type1WeightSets;
  /// the set of type2 weights (for integration of gradient interpolants)
  /// for each derivative component and for each point in the sparse grid
  RealMatrix2DArray type2WeightSets;

  /// stored type 1 weight sets for restoration to type1WeightSets
  std::map<UShortArray, RealVector> savedT1WtSets;
  /// stored type 2 weight sets for restoration to type2WeightSets
  std::map<UShortArray, RealMatrix> savedT2WtSets;
};


inline HierarchSparseGridDriver::HierarchSparseGridDriver():
  SparseGridDriver(), nestedGrid(true)
{ }


inline HierarchSparseGridDriver::
HierarchSparseGridDriver(unsigned short ssg_level, const RealVector& dim_pref):
  SparseGridDriver(ssg_level, dim_pref), nestedGrid(true)
{ }


inline HierarchSparseGridDriver::~HierarchSparseGridDriver()
{ }


inline int HierarchSparseGridDriver::unique_trial_points() const
{
  size_t lev = index_norm(trialSet);
  return collocKey[lev].back().size();
}


inline const UShort3DArray& HierarchSparseGridDriver::
smolyak_multi_index() const
{ return smolyakMultiIndex; }


inline const UShort4DArray& HierarchSparseGridDriver::collocation_key() const
{ return collocKey; }


inline const Sizet3DArray& HierarchSparseGridDriver::collocation_indices() const
{ return collocIndices; }


inline const RealVector2DArray& HierarchSparseGridDriver::
type1_weight_set_arrays() const
{ return type1WeightSets; }


inline const RealMatrix2DArray& HierarchSparseGridDriver::
type2_weight_set_arrays() const
{ return type2WeightSets; }


inline void HierarchSparseGridDriver::
level_to_delta_order(const UShortArray& levels, UShort2DArray& delta_quad)
{
  size_t i, j, num_lev = levels.size(), num_delta;
  if (delta_quad.size() != num_lev)
    delta_quad.resize(num_lev);
  for (i=0; i<num_lev; ++i) {
    unsigned short lev_i = levels[i], ord_i, ord_im1 = 0;
    level_to_order(i, lev_i, ord_i);
    if (lev_i > 0) level_to_order(i, --lev_i, ord_im1);
    num_delta = ord_i - ord_im1;
    UShortArray& delta_quad_i = delta_quad[i];
    delta_quad_i.resize(num_delta);
    switch(collocRules[i]) {
    case GAUSS_PATTERSON: // open nested
      for (j=0; j<num_delta; ++j)
	delta_quad_i[j] = 2*j;   // 0,2,4,6,8,...
      break;
    case NEWTON_COTES: case CLENSHAW_CURTIS: // closed nested
      for (j=0; j<num_delta; ++j)
	delta_quad_i[j] = 2*j+1; // 1,3,5,7,9,...
      break;
    case GENZ_KEISTER: // open nested table lookup
      // TO DO
      break;
    default:
      PCerr << "Error: bad rule type in level_to_delta_order()" << std::endl;
      abort_handler(-1);
      break;
    }
  }
}

} // namespace Pecos

#endif

/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        SharedProjectOrthogPolyApproxData
//- Description:  Class for Multivariate Orthogonal Polynomial Approximations
//-               
//- Owner:        Mike Eldred

#ifndef SHARED_PROJECT_ORTHOG_POLY_APPROX_DATA_HPP
#define SHARED_PROJECT_ORTHOG_POLY_APPROX_DATA_HPP

#include "SharedOrthogPolyApproxData.hpp"

namespace Pecos {

class TensorProductDriver;
class CombinedSparseGridDriver;
class CubatureDriver;


/// Derived approximation class for multivariate orthogonal polynomial
/// approximation with coefficient estimation via numerical integration.

/** The SharedProjectOrthogPolyApproxData class provides a global approximation
    based on multivariate orthogonal polynomials, where the coefficients are
    computed using numerical integration approaches such as quadrature,
    cubature, sparse grids, and random sampling.  It is used primarily for
    polynomial chaos expansion aproaches to UQ. */

class SharedProjectOrthogPolyApproxData: public SharedOrthogPolyApproxData
{
  //
  //- Heading: Friends
  //

  friend class ProjectOrthogPolyApproximation;

public:

  //
  //- Heading: Constructor and destructor
  //

  /// lightweight constructor
  SharedProjectOrthogPolyApproxData(short basis_type,
				    const UShortArray& approx_order,
				    size_t num_vars);
  /// full constructor
  SharedProjectOrthogPolyApproxData(short basis_type,
				    const UShortArray& approx_order,
				    size_t num_vars,
				    const ExpansionConfigOptions& ec_options,
				    const BasisConfigOptions& bc_options);
  /// destructor
  ~SharedProjectOrthogPolyApproxData();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void allocate_data();

  void increment_data();
  void decrement_data();
  void pre_restore_data();
  void post_restore_data();
  void pre_finalize_data();
  void post_finalize_data();

  void store_data();
  void pre_combine_data(short combine_type);
  void post_combine_data(short combine_type);

  void increment_component_sobol();

  //
  //- Heading: Member functions
  //

  /// return driverRep cast to requested derived type
  TensorProductDriver*      tpq_driver();
  /// return driverRep cast to requested derived type
  CombinedSparseGridDriver* csg_driver();
  /// return driverRep cast to requested derived type
  CubatureDriver*           cub_driver();

private:

  //
  //- Heading: Member functions
  //

  /// initialize multi_index using a sparse grid expansion
  void sparse_grid_multi_index(UShort2DArray& multi_index);
  // initialize tp_multi_index from tpMultiIndexMap
  //void map_tensor_product_multi_index(UShort2DArray& tp_multi_index,
  //				        size_t tp_index);

  /// convert quadrature orders to integrand orders using rigorous mappings
  void quadrature_order_to_integrand_order(const UShortArray& quad_order,
					   UShortArray& int_order);
  /// convert integrand orders to expansion orders using rigorous mappings
  void integrand_order_to_expansion_order(const UShortArray& int_order,
					  UShortArray& exp_order);
  /// convert a sparse grid index set and a growth setting to an integrand_order
  void sparse_grid_level_to_expansion_order(const UShortArray& levels,
					    UShortArray& exp_order);
                       //, short growth_rate = UNRESTRICTED_GROWTH);

  /// update the total Pareto set with new Pareto-optimal polynomial indices
  void update_pareto(const UShort2DArray& new_pareto,
		     UShort2DArray& total_pareto);
  /// assess whether new_pareto is dominated by total_pareto
  bool assess_dominance(const UShort2DArray& new_pareto,
			const UShort2DArray& total_pareto);
  /// assess bi-directional dominance for a new polynomial index set 
  /// against an incumbent polynomial index set
  void assess_dominance(const UShortArray& new_order,
			const UShortArray& existing_order,
			bool& new_dominated, bool& existing_dominated);

  /// Perform efficient calculation of tensor-product value via Horner's rule
  Real tensor_product_value(const RealVector& x, const RealVector& tp_coeffs,
			    const UShortArray& approx_order,
			    const UShort2DArray& tp_mi,
			    RealVector& accumulator);

  //
  //- Heading: Data
  //

  /// storage of level multi-index (levels for tensor or sparse grids)
  /// for subsequent restoration
  UShort2DArray storedLevMultiIndex;
  /// combination type for stored expansions; cached in class to bridge
  /// combine_coefficients() and compute_numerical_response_moments()
  short storedExpCombineType;

  /// numSmolyakIndices-by-numTensorProductPts-by-numVars array for
  /// identifying the orders of the one-dimensional orthogonal polynomials
  /// contributing to each of the multivariate orthogonal polynomials.
  /** For nested rules (GP, CC, or GK), the integration driver's collocKey
      is insufficient and we must track expansion orders separately. */
  UShort3DArray tpMultiIndex;
  /// sparse grid bookkeeping: mapping from num tensor-products by 
  /// tensor-product multi-indices into aggregated multiIndex
  Sizet2DArray tpMultiIndexMap;
  /// sparse grid bookkeeping: reference points for tpMultiIndexMap
  SizetArray tpMultiIndexMapRef;

  /// saved instances of tpMultiIndex that were computed but not selected
  std::deque<UShort2DArray> savedTPMultiIndex;
  /// saved instances of tpMultiIndexMap that were computed but not selected
  std::deque<SizetArray> savedTPMultiIndexMap;
  /// saved instances of tpMultiIndexMapRef that were computed but not selected
  std::deque<size_t> savedTPMultiIndexMapRef;

  /// index into saved sets of data to be restored (stored in this
  /// class for used by each ProjectOrthogPolyApproximation)
  size_t restoreIndex;
};


inline SharedProjectOrthogPolyApproxData::
SharedProjectOrthogPolyApproxData(short basis_type,
				  const UShortArray& approx_order,
				  size_t num_vars):
  SharedOrthogPolyApproxData(basis_type, approx_order, num_vars),
  storedExpCombineType(NO_COMBINE)
{ }


inline SharedProjectOrthogPolyApproxData::
SharedProjectOrthogPolyApproxData(short basis_type,
				  const UShortArray& approx_order,
				  size_t num_vars,
				  const ExpansionConfigOptions& ec_options,
				  const BasisConfigOptions& bc_options):
  SharedOrthogPolyApproxData(basis_type, approx_order, num_vars, ec_options,
			     bc_options), storedExpCombineType(NO_COMBINE)
{ }


inline SharedProjectOrthogPolyApproxData::~SharedProjectOrthogPolyApproxData()
{ }


inline TensorProductDriver* SharedProjectOrthogPolyApproxData::tpq_driver()
{ return (TensorProductDriver*)driverRep; }


inline CombinedSparseGridDriver* SharedProjectOrthogPolyApproxData::csg_driver()
{ return (CombinedSparseGridDriver*)driverRep; }


inline CubatureDriver* SharedProjectOrthogPolyApproxData::cub_driver()
{ return (CubatureDriver*)driverRep; }

} // namespace Pecos

#endif

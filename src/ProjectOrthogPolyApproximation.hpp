/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        ProjectOrthogPolyApproximation
//- Description:  Class for Multivariate Orthogonal Polynomial Approximations
//-               
//- Owner:        Mike Eldred

#ifndef PROJECT_ORTHOG_POLY_APPROXIMATION_HPP
#define PROJECT_ORTHOG_POLY_APPROXIMATION_HPP

#include "OrthogPolyApproximation.hpp"

namespace Pecos {


/// Derived approximation class for multivariate orthogonal polynomial
/// approximation with coefficient estimation via numerical integration.

/** The ProjectOrthogPolyApproximation class provides a global approximation
    based on multivariate orthogonal polynomials, where the coefficients are
    computed using numerical integration approaches such as quadrature,
    cubature, sparse grids, and random sampling.  It is used primarily for
    polynomial chaos expansion aproaches to UQ. */

class ProjectOrthogPolyApproximation: public OrthogPolyApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  ProjectOrthogPolyApproximation(const UShortArray& approx_order,
				 size_t num_vars, bool use_derivs);
  /// destructor
  ~ProjectOrthogPolyApproximation();

  //
  //- Heading: Member functions
  //

protected:

  //
  //- Heading: Virtual function redefinitions and member functions
  //

  int  min_coefficients() const;
  void compute_coefficients();

  void increment_coefficients();
  void decrement_coefficients();
  void restore_coefficients();
  void finalize_coefficients();

  void store_coefficients();
  void combine_coefficients(short combine_type);

  /// initialize polynomialBasis, multiIndex, et al.
  void allocate_arrays();

  Real value(const RealVector& x);
  Real stored_value(const RealVector& x);

  /// compute numerical moments to order 4 and expansion moments to order 2
  void compute_moments();
  /// compute expansion moments in all-variables mode to order 2
  void compute_moments(const RealVector& x);

private:

  //
  //- Heading: Member functions
  //

  /// Perform efficient calculation of tensor-product value via Horner's rule
  Real tensor_product_value(const RealVector& x, const RealVector& tp_coeffs,
			    const UShortArray& approx_order,
			    const UShort2DArray& tp_mi,
			    RealVector& accumulator);

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

  /// extract tp_data_points from surrData and tp_weights from
  /// driverRep->type1CollocWts1D
  void integration_data(size_t tp_index, SDVArray& tp_data_vars,
			SDRArray& tp_data_resp, RealVector& tp_weights);
  /// computes the chaosCoeffs via numerical integration (expCoeffsSolnApproach
  /// can be QUADRATURE, CUBATURE, or COMBINED_SPARSE_GRID)
  void integrate_expansion(const UShort2DArray& multi_index,
			   const SDVArray& data_vars, const SDRArray& data_resp,
			   const RealVector& wt_sets, RealVector& exp_coeffs,
			   RealMatrix& exp_coeff_grads);

  /// computes the chaosCoeffs via averaging of samples
  /// (expCoeffsSolnApproach is SAMPLING)
  void expectation();

  /// update expansion{Coeffs,CoeffGrads} by adding one or more tensor-product
  /// expansions and updating all Smolyak coefficients
  void append_tensor_expansions(size_t start_index);

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

  //
  //- Heading: Data
  //

  /// previous expansionCoeffs (aggregated total, not tensor-product
  /// contributions) prior to append_tensor_expansions()
  RealVector prevExpCoeffs;
  /// previous expansionCoeffGrads (aggregated total, not tensor-product
  /// contributions) prior to append_tensor_expansions()
  RealMatrix prevExpCoeffGrads;

  /// storage of level multi-index (levels for tensor or sparse grids)
  /// for subsequent restoration
  UShort2DArray storedLevMultiIndex;

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
  /// the set of tensor-product contributions to expansionCoeffs
  RealVectorArray tpExpansionCoeffs;
  /// the set of tensor-product contributions to expansionCoeffGrads
  RealMatrixArray tpExpansionCoeffGrads;

  /// saved instances of tpMultiIndex that were computed but not selected
  std::deque<UShort2DArray> savedTPMultiIndex;
  /// saved instances of tpMultiIndexMap that were computed but not selected
  std::deque<SizetArray> savedTPMultiIndexMap;
  /// saved instances of tpMultiIndexMapRef that were computed but not selected
  std::deque<size_t> savedTPMultiIndexMapRef;
  /// saved instances of tpExpansionCoeffs that were computed but not selected
  std::deque<RealVector> savedTPExpCoeffs;
  /// saved tpExpansionCoeffGrads instances that were computed but not selected
  std::deque<RealMatrix> savedTPExpCoeffGrads;
};


inline ProjectOrthogPolyApproximation::
ProjectOrthogPolyApproximation(const UShortArray& approx_order, size_t num_vars,
			       bool use_derivs):
  OrthogPolyApproximation(approx_order, num_vars, use_derivs)
{ }


inline ProjectOrthogPolyApproximation::~ProjectOrthogPolyApproximation()
{ }


inline void ProjectOrthogPolyApproximation::compute_moments()
{
  // standard variables mode supports two expansion and four numerical moments
  mean(); variance(); // updates expansionMoments[0] and [1]
  //standardize_moments(expansionMoments);
  if (expConfigOptions.expCoeffsSolnApproach != SAMPLING)
    compute_numerical_response_moments(4);
}


inline void ProjectOrthogPolyApproximation::compute_moments(const RealVector& x)
{
  // all variables mode only supports first two moments
  mean(x); variance(x); // updates expansionMoments[0] and [1]
  //standardize_moments(expansionMoments);
  //compute_numerical_response_moments(2, x); // TO DO
}

} // namespace Pecos

#endif

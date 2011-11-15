/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        OrthogPolyApproximation
//- Description:  Class for Multivariate Orthogonal Polynomial Approximations
//-               
//- Owner:        Mike Eldred

#ifndef ORTHOG_POLY_APPROXIMATION_HPP
#define ORTHOG_POLY_APPROXIMATION_HPP

#include "PolynomialApproximation.hpp"
#include "NumericGenOrthogPolynomial.hpp"


namespace Pecos {

// special values for quadratureExpansion and sparseGridExpansion
enum { TENSOR_INT_TOTAL_ORD_EXP,      TENSOR_INT_TENSOR_EXP,
       TENSOR_INT_TENSOR_SUM_EXP,     SPARSE_INT_TOTAL_ORD_EXP,
       SPARSE_INT_HEUR_TOTAL_ORD_EXP, SPARSE_INT_TENSOR_SUM_EXP,
       SPARSE_INT_RESTR_TENSOR_SUM_EXP };


/// Derived approximation class for orthogonal polynomials (global
/// approximation).

/** The OrthogPolyApproximation class provides a global approximation
    based on orthogonal polynomials.  It is used primarily for polynomial
    chaos expansions (for stochastic finite element approaches to
    uncertainty quantification). */

class OrthogPolyApproximation: public PolynomialApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  OrthogPolyApproximation(const UShortArray& approx_order, size_t num_vars,
			  bool use_derivs);
  /// destructor
  ~OrthogPolyApproximation();

  //
  //- Heading: Member functions
  //

  /// set numExpansionTerms
  void expansion_terms(int exp_terms);
  /// get numExpansionTerms
  int expansion_terms() const;

  /// invoke initialize_basis_types(), initialize_polynomial_basis(),
  /// and, if needed, update_basis_distribution_parameters()
  void construct_basis(const ShortArray& u_types, const DistributionParams& dp,
		       const BasisConfigOptions& bc_options);

  /// set basisTypes
  void basis_types(const ShortArray& basis_types);
  /// get basisTypes
  const ShortArray& basis_types() const;

  /// get polynomialBasis
  const std::vector<BasisPolynomial>& polynomial_basis() const;
  /// set polynomialBasis
  void polynomial_basis(const std::vector<BasisPolynomial>& poly_basis);

  /// set NumericGenOrthogPolynomial::coeffsNormsFlag
  void coefficients_norms_flag(bool flag);

protected:

  //
  //- Heading: Virtual function redefinitions and member functions
  //

  int min_coefficients() const;

  void compute_coefficients();
  void increment_coefficients();
  void decrement_coefficients();
  void restore_coefficients();
  void finalize_coefficients();
  void store_coefficients();
  void combine_coefficients(short combine_type);

  void print_coefficients(std::ostream& s) const;

  const RealVector& approximation_coefficients() const;
  void approximation_coefficients(const RealVector& approx_coeffs);

  /// initialize basis_types from u_types
  bool initialize_basis_types(const ShortArray& u_types,
			      ShortArray& basis_types);

  /// initialize polynomialBasis, multiIndex, et al.
  void allocate_arrays();

  /// Performs global sensitivity analysis using Sobol' Indices;
  /// computes component (main and interaction) effects
  void compute_component_effects();
  /// Performs global sensitivity analysis using Sobol' Indices;
  /// computes total effects
  void compute_total_effects();

  /// estimate chaos expansion coefficient decay rates for each random
  /// variable dimension using linear least squares in semilog space
  const RealVector& dimension_decay_rates();

  /// uniformly increment approxOrder
  void increment_order();

  /// retrieve the response PCE value for a given parameter vector
  Real value(const RealVector& x);
  /// retrieve the gradient of the response PCE with respect to all
  /// variables included in the polynomial basis (e.g., probabilistic
  /// variables) for a given parameter vector
  const RealVector& gradient_basis_variables(const RealVector& x);
  /// retrieve the gradient of the response PCE with respect to variables
  /// included in the polynomial basis (e.g., probabilistic or "all"
  /// variables) for a given parameter vector and a given DVV subset
  const RealVector& gradient_basis_variables(const RealVector& x,
					     const SizetArray& dvv);
  /// retrieve the gradient of the response PCE with respect to variables not
  /// included in the polynomial basis (nonprobabilistic variables such as
  /// design or epistemic when not in "all" mode) for a given parameter vector
  const RealVector& gradient_nonbasis_variables(const RealVector& x);

  Real stored_value(const RealVector& x);
  const RealVector& stored_gradient_basis_variables(const RealVector& x);
  const RealVector& stored_gradient_nonbasis_variables(const RealVector& x);

  Real mean();
  Real mean(const RealVector& x);
  const RealVector& mean_gradient();
  const RealVector& mean_gradient(const RealVector& x, const SizetArray& dvv);

  Real variance();
  Real variance(const RealVector& x);
  const RealVector& variance_gradient();
  const RealVector& variance_gradient(const RealVector& x,
				      const SizetArray& dvv);

  Real covariance(PolynomialApproximation* poly_approx_2);
  Real covariance(const RealVector& x, PolynomialApproximation* poly_approx_2);

  /// compute numerical moments to order 4 and expansion moments to order 2
  void compute_moments();
  /// compute expansion moments in all-variables mode to order 2
  void compute_moments(const RealVector& x);
  /// return expansionMoments
  const RealVector& moments() const;

  /// returns the norm-squared of a particular multivariate polynomial,
  /// treating all variables as probabilistic
  Real norm_squared(const UShortArray& indices);
  /// returns the norm-squared of a particular multivariate polynomial,
  /// treating a subset of the variables as probabilistic
  Real norm_squared_random(const UShortArray& indices);

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
  /// convert sparse grid level to expansion orders using available heuristics
  void heuristic_sparse_grid_level_to_expansion_order(unsigned short ssg_level,
						      UShortArray& exp_order);
  /// convert a sparse grid index set and a growth setting to an integrand_order
  void sparse_grid_levels_to_expansion_order(const UShortArray& levels,
    UShortArray& exp_order, short growth_rate = UNRESTRICTED_GROWTH);

  /// append multi-indices from app_multi_index that do not already
  /// appear in multi_index
  void append_multi_index(const UShort2DArray& app_multi_index,
			  UShort2DArray& multi_index);
  /// append multi-indices from app_multi_index that do not already
  /// appear in multi_index; define multi_index_map and multi_index_map_ref
  void append_multi_index(const UShort2DArray& app_multi_index,
			  UShort2DArray& multi_index,
			  Sizet2DArray& multi_index_map,
			  SizetArray& multi_index_map_ref);
  /// append multi-indices from app_multi_index that do not already
  /// appear in multi_index, using previously defined multi_index_map
  /// and multi_index_map_ref for mapping
  void append_multi_index(const UShort2DArray& app_multi_index,
			  SizetArray& multi_index_map,
			  size_t& multi_index_map_ref,
			  UShort2DArray& multi_index);

  /// overlay the passed tensor-product expansion with the aggregate
  /// expansion{Coeffs,CoeffGrads}
  void overlay_expansion(const SizetArray& multi_index_map,
			 const RealVector& exp_coeffs,
			 const RealMatrix& exp_grads, int coeff);
  /// multiply current expansion ("a") with incoming expansion ("b")
  /// and store in product expansion ("c")
  void multiply_expansion(const UShort2DArray& multi_index_b,
			  const RealVector& exp_coeffs_b,
			  const RealMatrix& exp_grads_b,
			  const UShort2DArray& multi_index_c);
  /// update expansion{Coeffs,CoeffGrads} by adding one or more tensor-product
  /// expansions and updating all Smolyak coefficients
  void append_tensor_expansions(size_t start_index);
  /// synchronize expansion{Coeffs,CoeffGrads} with an updated multiIndex
  void resize_expansion();

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

  /// calculate a particular multivariate orthogonal polynomial value
  /// evaluated at a particular parameter set
  Real multivariate_polynomial(const RealVector& xi,const UShortArray& indices);
  /// calculate a particular multivariate orthogonal polynomial gradient
  /// evaluated at a particular parameter set
  const RealVector& multivariate_polynomial_gradient(const RealVector& xi,
    const UShortArray& indices);
  /// calculate a particular multivariate orthogonal polynomial gradient with
  /// respect to specified dvv and evaluated at a particular parameter set
  const RealVector& multivariate_polynomial_gradient(const RealVector& xi,
    const UShortArray& indices, const SizetArray& dvv);

  /// computes the chaosCoeffs via linear regression
  /// (expCoeffsSolnApproach is REGRESSION)
  void regression();
  /// computes the chaosCoeffs via averaging of samples
  /// (expCoeffsSolnApproach is SAMPLING)
  void expectation();

  /// screen data sets for samples with Inf/Nan that should be excluded;
  /// defines failedIndices
  void sample_checks();
  /// perform sanity checks prior to numerical integration
  void integration_checks();
  /// extract tp_data_points from surrData and tp_weights from
  /// driverRep->collocWts1D
  void integration_data(size_t tp_index, SDVArray& tp_data_vars,
			SDRArray& tp_data_resp, RealVector& tp_weights);
  /// computes the chaosCoeffs via numerical integration
  /// (expCoeffsSolnApproach is QUADRATURE, CUBATURE, or SPARSE_GRID)
  void integrate_expansion(const UShort2DArray& multi_index,
			   const SDVArray& data_vars, const SDRArray& data_resp,
			   const RealVector& wt_sets, RealVector& exp_coeffs,
			   RealMatrix& exp_coeff_grads);

  /// update numericalMoments using numerical integration applied
  /// directly to surrData
  void compute_numerical_response_moments(size_t num_moments);

  /// cross-validates alternate gradient expressions
  void gradient_check();

  //
  //- Heading: Data
  //

  /// number of terms in orthogonal polynomial expansion (length of chaosCoeffs)
  int numExpansionTerms;
  /// order of orthogonal polynomial expansion
  UShortArray approxOrder;
  /// previous value of approxOrder; used for detecting when a multiIndex
  /// update is needed
  UShortArray approxOrderPrev;
  /// flag identifying use of a partial total-order expansion due to
  /// user specification of a numExpansionTerms value that does not
  /// define a full total-order expansion
  bool partialOrder;

  /// array of basis types for each one-dimensional orthogonal polynomial:
  /// HERMITE_ORTHOG, LEGENDRE_ORTHOG, LAGUERRE_ORTHOG, JACOBI_ORTHOG,
  /// GEN_LAGUERRE_ORTHOG, CHEBYSHEV_ORTHOG, or NUM_GEN_ORTHOG
  ShortArray basisTypes;

  /// array of one-dimensional basis polynomial objects which are used in
  /// constructing the multivariate orthogonal/interpolation polynomials
  std::vector<BasisPolynomial> polynomialBasis;

  /// numExpansionTerms-by-numVars array for identifying the orders of
  /// the one-dimensional orthogonal polynomials contributing to each
  /// of the multivariate orthogonal polynomials
  UShort2DArray multiIndex;

  /// the coefficients of the expansion
  RealVector expansionCoeffs;
  /// the gradients of the expansion coefficients
  /** may be interpreted as either the gradients of the expansion
      coefficients or the coefficients of expansions for the response
      gradients.  This array is used when sensitivities of moments are
      needed with respect to variables that do not appear in the
      expansion (e.g., with respect to design or epistemic variables
      for an expansion only over probabilistic variables). */
  RealMatrix expansionCoeffGrads;

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

  /// copy of multiIndex (aggregated total, not tensor-product contributions)
  /// stored in store_coefficients() for use in combine_coefficients()
  UShort2DArray storedMultiIndex;
  /// copy of expansionCoeffs (aggregated total, not tensor-product contribs)
  /// stored in store_coefficients() for use in combine_coefficients()
  RealVector storedExpCoeffs;
  /// copy of expansionCoeffGrads (aggregated total, not indiv tensor-products)
  /// stored in store_coefficients() for use in combine_coefficients()
  RealMatrix storedExpCoeffGrads;
  /// copy of approxOrder stored in store_coefficients() for use in
  /// combine_coefficients()
  UShortArray storedApproxOrder;
  /// combination type for stored expansions; cached in class to bridge
  /// combine_coefficients() and compute_numerical_response_moments()
  short storedExpCombineType;

  /// previous expansionCoeffs (aggregated total, not tensor-product
  /// contributions) prior to append_tensor_expansions()
  RealVector prevExpCoeffs;
  /// previous expansionCoeffGrads (aggregated total, not tensor-product
  /// contributions) prior to append_tensor_expansions()
  RealMatrix prevExpCoeffGrads;

  /// Data vector for storing the gradients of individual expansion term
  /// polynomials.  Called in multivariate_polynomial_gradient().
  RealVector mvpGradient;

  /// switch for formulation of orthogonal polynomial expansion
  /// integrated with tensor-product quadrature:
  /// TENSOR_INT_TOTAL_ORD_EXP or TENSOR_INT_TENSOR_EXP expansion.
  short quadratureExpansion;
  /// switch for formulation of orthogonal polynomial expansion for
  /// sparse grids: TENSOR_INT_TENSOR_SUM_EXP, SPARSE_INT_TENSOR_SUM_EXP,
  /// SPARSE_INT_RESTR_TENSOR_SUM_EXP, SPARSE_INT_TOTAL_ORD_EXP, or
  /// SPARSE_INT_HEUR_TOTAL_ORD_EXP expansion.
  short sparseGridExpansion;

  /// spectral coefficient decay rates estimated by LLS on log of
  /// univariate expansion coefficients
  RealVector decayRates;

  /// list of failed evaluation indices defined in sample_checks() and
  /// used for fault tolerance in regression() and expectation()
  SizetList failedIndices;
};


inline OrthogPolyApproximation::
OrthogPolyApproximation(const UShortArray& approx_order, size_t num_vars,
			bool use_derivs):
  PolynomialApproximation(num_vars, use_derivs), numExpansionTerms(0),
  approxOrder(approx_order), partialOrder(false),
  storedExpCombineType(NO_COMBINE), quadratureExpansion(TENSOR_INT_TENSOR_EXP),
//quadratureExpansion(TENSOR_INT_TOTAL_ORDER_EXP),
  sparseGridExpansion(TENSOR_INT_TENSOR_SUM_EXP)
//sparseGridExpansion(SPARSE_INT_TOTAL_ORDER_EXP)
//sparseGridExpansion(SPARSE_INT_HEUR_TOTAL_ORDER_EXP)
//sparseGridExpansion(SPARSE_INT_TENSOR_SUM_EXP)
//sparseGridExpansion(SPARSE_INT_RESTR_TENSOR_SUM_EXP)
// Note: for sparseGridExpansion == SPARSE_INT_*, all_variables mode requires
//       track_wts = true in Dakota::NonDExpansion::construct_sparse_grid().
{ }


inline OrthogPolyApproximation::~OrthogPolyApproximation()
{ }


inline void OrthogPolyApproximation::increment_order()
{
  // increment approxOrder
  for (size_t i=0; i<numVars; ++i)
    approxOrder[i] += 1;
  // need numExpansionTerms updated for use in NonDPolynomialChaos::
  // increment_expansion(), but multiIndex update can wait until
  // compute_coefficients() => allocate_arrays().
  numExpansionTerms = total_order_terms(approxOrder);
}


inline void OrthogPolyApproximation::compute_moments()
{
  // standard variables mode supports two expansion and four numerical moments
  expansionMoments.sizeUninitialized(2); mean(); variance();
  if (expConfigOptions.expCoeffsSolnApproach == QUADRATURE ||
      expConfigOptions.expCoeffsSolnApproach == CUBATURE   ||
      expConfigOptions.expCoeffsSolnApproach == SPARSE_GRID)
    compute_numerical_response_moments(4);
}


inline void OrthogPolyApproximation::compute_moments(const RealVector& x)
{
  // all variables mode only supports first two moments
  expansionMoments.sizeUninitialized(2); mean(x); variance(x);
  //compute_numerical_response_moments(2, x); // TO DO
}


inline const RealVector& OrthogPolyApproximation::moments() const
{ return expansionMoments; }


inline void OrthogPolyApproximation::expansion_terms(int exp_terms)
{ numExpansionTerms = exp_terms; }


inline int OrthogPolyApproximation::expansion_terms() const
{ return numExpansionTerms; }


/** This function is invoked to create basisTypes and polynomialBasis
    for cases where they have not already been created by an
    IntegrationDriver (i.e., expansion_samples or regression). */
inline void OrthogPolyApproximation::
construct_basis(const ShortArray& u_types, const DistributionParams& dp,
		const BasisConfigOptions& bc_options)
{
  bool dist_params = initialize_basis_types(u_types, basisTypes);
  ShortArray colloc_rules;
  initialize_collocation_rules(u_types, bc_options, colloc_rules);
  initialize_polynomial_basis(basisTypes, colloc_rules, polynomialBasis);
  if (dist_params)
    update_basis_distribution_parameters(u_types, dp, polynomialBasis);
}


inline void OrthogPolyApproximation::basis_types(const ShortArray& basis_types)
{ basisTypes = basis_types; }


inline const ShortArray& OrthogPolyApproximation::basis_types() const
{ return basisTypes; }


inline const std::vector<BasisPolynomial>& OrthogPolyApproximation::
polynomial_basis() const
{ return polynomialBasis; }


inline void OrthogPolyApproximation::
polynomial_basis(const std::vector<BasisPolynomial>& poly_basis)
{ polynomialBasis = poly_basis; }


inline void OrthogPolyApproximation::coefficients_norms_flag(bool flag)
{
  size_t i, num_basis = basisTypes.size();
  for (i=0; i<num_basis; ++i)
    if (basisTypes[i] == NUM_GEN_ORTHOG)
      ((NumericGenOrthogPolynomial*)polynomialBasis[i].polynomial_rep())
	->coefficients_norms_flag(flag);
}


inline const RealVector& OrthogPolyApproximation::
approximation_coefficients() const
{ return expansionCoeffs; }


inline void OrthogPolyApproximation::
approximation_coefficients(const RealVector& approx_coeffs)
{ expansionCoeffs = approx_coeffs; }


inline void OrthogPolyApproximation::resize_expansion()
{
  numExpansionTerms = multiIndex.size();
  if (expConfigOptions.expansionCoeffFlag)
    expansionCoeffs.resize(numExpansionTerms); // new terms initialized to 0
  if (expConfigOptions.expansionCoeffGradFlag) {
    size_t num_deriv_vars = expansionCoeffGrads.numRows();
    expansionCoeffGrads.reshape(num_deriv_vars, numExpansionTerms);
    // new terms initialized to 0
  }
}


inline Real OrthogPolyApproximation::
multivariate_polynomial(const RealVector& xi, const UShortArray& indices)
{
  unsigned short order_1d;
  Real mvp = 1.;
  for (size_t i=0; i<numVars; ++i) {
    order_1d = indices[i];
    if (order_1d)
      mvp *= polynomialBasis[i].type1_value(xi[i], order_1d);
  }
  return mvp;
}


inline const RealVector& OrthogPolyApproximation::
multivariate_polynomial_gradient(const RealVector& xi,
				 const UShortArray& indices)
{
  if (mvpGradient.length() != numVars)
    mvpGradient.sizeUninitialized(numVars);
  size_t i, j;
  for (i=0; i<numVars; ++i) {
    Real& mvp_grad_i = mvpGradient[i];
    mvp_grad_i = 1.0;
    // differentiation of product of 1D polynomials
    for (j=0; j<numVars; ++j)
      mvp_grad_i *= (j == i) ?
	polynomialBasis[j].type1_gradient(xi[j], indices[j]) :
	polynomialBasis[j].type1_value(xi[j],    indices[j]);
  }
  return mvpGradient;
}


inline const RealVector& OrthogPolyApproximation::
multivariate_polynomial_gradient(const RealVector& xi,
				 const UShortArray& indices,
				 const SizetArray& dvv)
{
  size_t i, j, deriv_index, num_deriv_vars = dvv.size();
  if (mvpGradient.length() != num_deriv_vars)
    mvpGradient.sizeUninitialized(num_deriv_vars);
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // *** requires an "All" view
    Real& mvp_grad_i = mvpGradient[i];
    mvp_grad_i = 1.0;
    // differentiation of product of 1D polynomials
    for (j=0; j<numVars; ++j)
      mvp_grad_i *= (j == deriv_index) ?
	polynomialBasis[j].type1_gradient(xi[j], indices[j]) :
	polynomialBasis[j].type1_value(xi[j],    indices[j]);
  }
  return mvpGradient;
}

} // namespace Pecos

#endif

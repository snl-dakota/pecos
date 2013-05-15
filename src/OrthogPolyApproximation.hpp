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

  /// get numExpansionTerms
  int expansion_terms() const;

  /// invoke initialize_basis_types(), initialize_polynomial_basis(),
  /// and, if needed, update_basis_distribution_parameters()
  void construct_basis(const ShortArray& u_types, const AleatoryDistParams& adp,
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

  void store_coefficients();
  void combine_coefficients(short combine_type);

  void print_coefficients(std::ostream& s) const;

  const RealVector& approximation_coefficients() const;
  void approximation_coefficients(const RealVector& approx_coeffs);

  void coefficient_labels(std::vector<std::string>& all_coeff_tags) const;

  /// initialize polynomialBasis, multiIndex, et al.
  void allocate_arrays();

  /// Performs global sensitivity analysis using Sobol' Indices;
  /// computes component (main and interaction) effects
  void compute_component_effects();
  /// Performs global sensitivity analysis using Sobol' Indices;
  /// computes total effects
  void compute_total_effects();

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

  /// compute expansion moments to order 2
  void compute_moments();
  /// compute expansion moments in all-variables mode to order 2
  void compute_moments(const RealVector& x);

  /// return expansionMoments
  const RealVector& moments() const;

  //
  //- Heading: Member functions
  //

  /// initialize basis_types from u_types
  bool initialize_basis_types(const ShortArray& u_types,
			      ShortArray& basis_types);

  /// size expansion{Coeffs,CoeffGrads} based on multiIndex
  void size_expansion();
  /// synchronize expansion{Coeffs,CoeffGrads} with an updated multiIndex
  void resize_expansion();

  /// returns the norm-squared of a particular multivariate polynomial,
  /// treating all variables as probabilistic
  Real norm_squared(const UShortArray& indices);
  /// returns the norm-squared of a particular multivariate polynomial,
  /// treating a subset of the variables as probabilistic
  Real norm_squared(const UShortArray& indices, const SizetList& rand_indices);

  /// calculate a particular multivariate orthogonal polynomial value
  /// evaluated at a particular parameter set
  Real multivariate_polynomial(const RealVector& x,const UShortArray& indices);
  /// calculate a particular multivariate orthogonal polynomial value over
  /// the nonrandom variable subset evaluated at a particular parameter set
  Real multivariate_polynomial(const RealVector& x,const UShortArray& indices,
			       const SizetList& non_rand_indices);

  /// compute multivariate orthogonal polynomial gradient for term
  /// corresponding to deriv_index, evaluated at x
  Real multivariate_polynomial_gradient(const RealVector& x, size_t deriv_index,
    const UShortArray& indices);
  /// compute multivariate orthogonal polynomial gradient for term corresponding
  /// to deriv_index, for nonrandom variable subset, and evaluated at x
  Real multivariate_polynomial_gradient(const RealVector& x, size_t deriv_index,
    const UShortArray& indices, const SizetList& non_rand_indices);

  /// calculate multivariate orthogonal polynomial gradient vector
  /// evaluated at a particular parameter set
  const RealVector& multivariate_polynomial_gradient_vector(const RealVector& x,
    const UShortArray& indices);
  /// calculate multivariate orthogonal polynomial gradient vector with
  /// respect to specified dvv and evaluated at a particular parameter set
  const RealVector& multivariate_polynomial_gradient_vector(const RealVector& x,
    const UShortArray& indices, const SizetArray& dvv);

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

  /// overlay the passed expansion with the aggregate
  /// expansion{Coeffs,CoeffGrads} as managed by the multi_index_map
  void overlay_expansion(const SizetArray& multi_index_map,
			 const RealVector& exp_coeffs,
			 const RealMatrix& exp_grads, int coeff);
  /// multiply current expansion ("a") with incoming expansion ("b")
  /// and store in product expansion ("c")
  void multiply_expansion(const UShort2DArray& multi_index_b,
			  const RealVector& exp_coeffs_b,
			  const RealMatrix& exp_grads_b,
			  const UShort2DArray& multi_index_c);

  /// estimate chaos expansion coefficient decay rates for each random
  /// variable dimension using linear least squares in semilog space
  const RealVector& dimension_decay_rates();

  /// perform sanity checks prior to numerical integration
  void integration_checks();
  /// tests 1D gradient computations (active in DEBUG compile mode)
  void gradient_check();

  /// update add_val and add_gradient based on surrData's failure map
  void fail_booleans(SizetShortMap::const_iterator& fit, size_t j,
		     bool& add_val, bool& add_grad);

  //
  //- Heading: Data
  //

  /// array of basis types for each one-dimensional orthogonal polynomial:
  /// HERMITE_ORTHOG, LEGENDRE_ORTHOG, LAGUERRE_ORTHOG, JACOBI_ORTHOG,
  /// GEN_LAGUERRE_ORTHOG, CHEBYSHEV_ORTHOG, or NUM_GEN_ORTHOG
  ShortArray basisTypes;

  /// array of one-dimensional basis polynomial objects which are used in
  /// constructing the multivariate orthogonal/interpolation polynomials
  std::vector<BasisPolynomial> polynomialBasis;

  /// order of orthogonal polynomial expansion
  UShortArray approxOrder;

  /// number of terms in orthogonal polynomial expansion (length of chaosCoeffs)
  int numExpansionTerms;

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

  /// copy of multiIndex (aggregated expansion) stored in
  /// store_coefficients() for use in combine_coefficients()
  UShort2DArray storedMultiIndex;
  /// copy of expansionCoeffs (aggregated expansion) stored in
  /// store_coefficients() for use in combine_coefficients()
  RealVector storedExpCoeffs;
  /// copy of expansionCoeffGrads (aggregated expansion) stored in
  /// store_coefficients() for use in combine_coefficients()
  RealMatrix storedExpCoeffGrads;
  /// copy of approxOrder stored in store_coefficients() for use in
  /// combine_coefficients()
  UShortArray storedApproxOrder;

private:

  //
  //- Heading: Member functions
  //

  /// test for nonzero indices in random variable subset
  bool zero_random(const UShortArray& indices) const;

  /// Generate the coefficient tag for expansion term i, variable j
  void get_tag(char* tag, size_t i, size_t j) const;

  //
  //- Heading: Data
  //

  /// previous value of approxOrder; used for detecting when a multiIndex
  /// update is needed
  UShortArray approxOrderPrev;

  /// Data vector for storing the gradients of individual expansion term
  /// polynomials (see multivariate_polynomial_gradient_vector())
  RealVector mvpGradient;

  /// spectral coefficient decay rates estimated by LLS on log of
  /// univariate expansion coefficients
  RealVector decayRates;
};


inline OrthogPolyApproximation::
OrthogPolyApproximation(const UShortArray& approx_order, size_t num_vars,
			bool use_derivs):
  PolynomialApproximation(num_vars, use_derivs), numExpansionTerms(0),
  approxOrder(approx_order)
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
  mean(); variance(); // updates expansionMoments[0] and [1]
  //standardize_moments(expansionMoments);
}


inline void OrthogPolyApproximation::compute_moments(const RealVector& x)
{
  // all variables mode only supports first two moments
  mean(x); variance(x); // updates expansionMoments[0] and [1]
  //standardize_moments(expansionMoments);
}


inline const RealVector& OrthogPolyApproximation::moments() const
{ return expansionMoments; }


inline int OrthogPolyApproximation::expansion_terms() const
{ return numExpansionTerms; }


/** This function is invoked to create basisTypes and polynomialBasis
    for cases where they have not already been created by an
    IntegrationDriver (i.e., expansion_samples or regression). */
inline void OrthogPolyApproximation::
construct_basis(const ShortArray& u_types, const AleatoryDistParams& adp,
		const BasisConfigOptions& bc_options)
{
  bool dist_params = initialize_basis_types(u_types, basisTypes);
  ShortArray colloc_rules;
  initialize_collocation_rules(u_types, bc_options, colloc_rules);
  initialize_polynomial_basis(basisTypes, colloc_rules, polynomialBasis);
  if (dist_params)
    update_basis_distribution_parameters(u_types, adp, polynomialBasis);
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


inline void OrthogPolyApproximation::size_expansion()
{
  numExpansionTerms = multiIndex.size();
  if (expConfigOptions.expansionCoeffFlag &&
      expansionCoeffs.length() != numExpansionTerms)
    expansionCoeffs.sizeUninitialized(numExpansionTerms);
  if (expConfigOptions.expansionCoeffGradFlag) {
    size_t num_deriv_vars = surrData.num_derivative_variables();
    if (expansionCoeffGrads.numRows() != num_deriv_vars ||
	expansionCoeffGrads.numCols() != numExpansionTerms)
      expansionCoeffGrads.shapeUninitialized(num_deriv_vars, numExpansionTerms);
  }
}


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
multivariate_polynomial(const RealVector& x, const UShortArray& indices)
{
  Real mvp = 1.; unsigned short order_1d;
  for (size_t i=0; i<numVars; ++i) {
    order_1d = indices[i];
    if (order_1d)
      mvp *= polynomialBasis[i].type1_value(x[i], order_1d);
  }
  return mvp;
}


/** All variables version. */
inline Real OrthogPolyApproximation::
multivariate_polynomial(const RealVector& x, const UShortArray& indices,
			const SizetList& non_rand_indices)
{
  Real mvp = 1.; SizetList::const_iterator cit;
  unsigned short order_1d; size_t i;
  for (cit=non_rand_indices.begin(); cit!=non_rand_indices.end(); ++cit) {
    i = *cit; order_1d = indices[i];
    if (order_1d)
      mvp *= polynomialBasis[i].type1_value(x[i], order_1d);
  }
  return mvp;
}


inline Real OrthogPolyApproximation::
multivariate_polynomial_gradient(const RealVector& x, size_t deriv_index,
				 const UShortArray& indices)
{
  Real mvp_grad = 1.;
  for (size_t k=0; k<numVars; ++k)
    mvp_grad *= (k == deriv_index) ?
      polynomialBasis[k].type1_gradient(x[k], indices[k]) :
      polynomialBasis[k].type1_value(x[k],    indices[k]);
  return mvp_grad;
}


/** All variables version. */
inline Real OrthogPolyApproximation::
multivariate_polynomial_gradient(const RealVector& x, size_t deriv_index,
				 const UShortArray& indices,
				 const SizetList& non_rand_indices)
{
  Real mvp_grad = 1.; SizetList::const_iterator cit; size_t k;
  for (cit=non_rand_indices.begin(); cit!=non_rand_indices.end(); ++cit) {
    k = *cit;
    mvp_grad *= (k == deriv_index) ?
      polynomialBasis[k].type1_gradient(x[k], indices[k]) :
      polynomialBasis[k].type1_value(x[k],    indices[k]);
  }
  return mvp_grad;
}


inline const RealVector& OrthogPolyApproximation::
multivariate_polynomial_gradient_vector(const RealVector& x,
					const UShortArray& indices)
{
  if (mvpGradient.length() != numVars)
    mvpGradient.sizeUninitialized(numVars);
  for (size_t i=0; i<numVars; ++i)
    mvpGradient[i] = multivariate_polynomial_gradient(x, i, indices);
  return mvpGradient;
}


inline const RealVector& OrthogPolyApproximation::
multivariate_polynomial_gradient_vector(const RealVector& x,
					const UShortArray& indices,
					const SizetArray& dvv)
{
  size_t i, j, deriv_index, num_deriv_vars = dvv.size();
  if (mvpGradient.length() != num_deriv_vars)
    mvpGradient.sizeUninitialized(num_deriv_vars);
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // *** requires an "All" view
    mvpGradient[i] = multivariate_polynomial_gradient(x, deriv_index, indices);
  }
  return mvpGradient;
}


inline bool OrthogPolyApproximation::
zero_random(const UShortArray& indices) const
{
  SizetList::const_iterator cit;
  for (cit=randomIndices.begin(); cit!=randomIndices.end(); ++cit)
    if (indices[*cit])
      return false;
  return true;
}


inline Real OrthogPolyApproximation::norm_squared(const UShortArray& indices)
{
  // the norm squared of a particular multivariate polynomial is the product of
  // the norms squared of the numVars univariate polynomials that comprise it.
  Real norm_sq = 1.; unsigned short order_1d;
  for (size_t i=0; i<numVars; ++i) {
    order_1d = indices[i];
    if (order_1d)
      norm_sq *= polynomialBasis[i].norm_squared(order_1d);
  }
  return norm_sq;
}


/** All variables version. */
inline Real OrthogPolyApproximation::
norm_squared(const UShortArray& indices, const SizetList& rand_indices)
{
  // the norm squared of a particular multivariate polynomial is the product of
  // the norms squared of the numVars univariate polynomials that comprise it.
  Real norm_sq = 1.; SizetList::const_iterator cit;
  unsigned short order_1d; size_t i;
  for (cit=rand_indices.begin(); cit!=rand_indices.end(); ++cit) {
    i = *cit; order_1d = indices[i];
    if (order_1d)
      norm_sq *= polynomialBasis[i].norm_squared(order_1d);
  }
  return norm_sq;
}


inline void OrthogPolyApproximation::
fail_booleans(SizetShortMap::const_iterator& fit, size_t j,
	      bool& add_val, bool& add_grad)
{
  if (fit != surrData.failed_response_data().end() && fit->first == j) {
    short fail_bits = fit->second;
    if (fail_bits & 1) add_val  = false;
    if (fail_bits & 2) add_grad = false;
    ++fit;
  }
}

} // namespace Pecos

#endif

/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
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
enum { TOTAL_ORDER, TENSOR_PRODUCT, HEURISTIC_TOTAL_ORDER, TENSOR_PRODUCT_SUM };

class DistributionParams;


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
  OrthogPolyApproximation(const UShortArray& approx_order, size_t num_vars);
  /// destructor
  ~OrthogPolyApproximation();

  //
  //- Heading: Member functions
  //

  /// set numExpansionTerms
  void expansion_terms(int exp_terms);
  /// get numExpansionTerms
  const int& expansion_terms() const;

  /// invoke distribution_types() and, if needed, distribution_parameters()
  void distributions(const ShortArray& u_types,
		     const IntArray& int_rules, const DistributionParams& dp);
  /// invoke distribution_types() and, if needed, distribution_parameters()
  static void distributions(const ShortArray& u_types,
			    const IntArray& int_rules,
			    const DistributionParams& dp,
			    std::vector<BasisPolynomial>& poly_basis,
			    ShortArray& basis_types, ShortArray& gauss_modes);

  /// allocate polynomialBasis and basisTypes based on u_types
  bool distribution_types(const ShortArray& u_types,
			  const IntArray& int_rules);
  /// allocate poly_basis and basis_types based on u_types
  static bool distribution_types(const ShortArray& u_types,
				 const IntArray& int_rules,
				 ShortArray& basis_types,
				 ShortArray& gauss_modes);

  /// allocate polynomialBasis based on basisTypes and gaussModes
  void distribution_basis();
  /// allocate poly_basis based on basis_types and gauss_modes
  static void distribution_basis(const ShortArray& basis_types,
				 const ShortArray& gauss_modes,
				 std::vector<BasisPolynomial>& poly_basis);

  /// pass distribution parameters from dp to polynomialBasis
  void distribution_parameters(const ShortArray& u_types,
			       const DistributionParams& dp);
  /// pass distribution parameters from dp to poly_basis
  static void distribution_parameters(const ShortArray& u_types,
				      const DistributionParams& dp,
				      std::vector<BasisPolynomial>& poly_basis);

  /// get polynomialBasis
  const std::vector<BasisPolynomial>& polynomial_basis() const;
  /// set polynomialBasis
  void polynomial_basis(const std::vector<BasisPolynomial>& poly_basis);

  /// set NumericGenOrthogPolynomial::coeffsNormsFlag
  void coefficients_norms_flag(bool flag);

  /// initialize polynomialBasis, multiIndex, et al.
  void allocate_arrays();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  int min_coefficients() const;

  /// find the coefficients for the expansion of multivariate
  /// orthogonal polynomials
  void find_coefficients();

  /// print the coefficients for the expansion
  void print_coefficients(std::ostream& s) const;

  /// Performs global sensitivity analysis using Sobol' Indices
  void compute_global_sensitivity();

  /// retrieve the response PCE value for a given parameter vector
  const Real& get_value(const RealVector& x);
  /// retrieve the response PCE gradient for a given parameter vector
  /// and default DVV
  const RealVector& get_gradient(const RealVector& x);
  /// retrieve the response PCE gradient for a given parameter vector
  /// and given DVV
  const RealVector& get_gradient(const RealVector& x,
				      const UIntArray& dvv);

  /// return the mean of the PCE, treating all variables as random
  const Real& get_mean();
  /// return the mean of the PCE for a given parameter vector,
  /// treating a subset of the variables as random
  const Real& get_mean(const RealVector& x);
  /// return the gradient of the PCE mean for a given parameter
  /// vector, treating all variables as random
  const RealVector& get_mean_gradient();
  /// return the gradient of the PCE mean for a given parameter vector
  /// and given DVV, treating a subset of the variables as random
  const RealVector& get_mean_gradient(const RealVector& x,
				      const UIntArray& dvv);

  /// return the variance of the PCE, treating all variables as random
  const Real& get_variance();
  /// return the variance of the PCE for a given parameter vector,
  /// treating a subset of the variables as random
  const Real& get_variance(const RealVector& x);
  /// return the gradient of the PCE variance for a given parameter
  /// vector, treating all variables as random
  const RealVector& get_variance_gradient();
  /// return the gradient of the PCE variance for a given parameter
  /// vector and given DVV, treating a subset of the variables as random
  const RealVector& get_variance_gradient(const RealVector& x,
					  const UIntArray& dvv);

  /// return the covariance of the PCE, treating all variables as random
  const Real& get_covariance(const RealVector& exp_coeffs_2);
  // return the covariance of the PCE for a given parameter vector,
  // treating a subset of the variables as random
  //const Real& get_covariance(const RealVector& x,
  //                           const RealVector& exp_coeffs_2);

  /// returns the norm-squared of a particular multivariate polynomial,
  /// treating all variables as random
  const Real& norm_squared(size_t expansion_index);
  /// returns the norm-squared of a particular multivariate polynomial,
  /// treating a subset of the variables as random
  const Real& norm_squared_random(size_t expansion_index);

private:

  //
  //- Heading: Member functions
  //

  /// initialize multiIndex using a sparse grid expansion
  void sparse_grid_multi_index(UShort2DArray& multi_index);

  /// convert quadrature orders to integrand orders using rigorous mappings
  void quadrature_order_to_integrand_order(const UShortArray& quad_order,
					   UShortArray& int_order);
  /// convert integrand orders to expansion orders using rigorous mappings
  void integrand_order_to_expansion_order(const UShortArray& int_order,
					  UShortArray& exp_order);
  /// convert sparse grid level to expansion orders using available heuristics
  void sparse_grid_level_to_expansion_order(unsigned short ssg_level,
					    UShortArray& exp_order);
  /// convert a level index set and a growth setting to an integrand_order
  void level_growth_to_integrand_order(const UShortArray& levels,
				       short exp_growth,
				       UShortArray& int_order);

  /// append multi-indices from tp_multi_index that do not already
  /// appear in multi_index
  void append_unique(const UShort2DArray& tp_multi_index,
		     UShort2DArray& multi_index);
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

  /// calculate a particular multivariate orthogonal polynomial
  /// evaluated at a particular parameter set
  Real multivariate_polynomial(const RealVector& x, size_t term);

  /// computes the chaosCoeffs via numerical integration
  /// (expCoeffsSolnApproach is QUADRATURE or SPARSE_GRID)
  void integration();
  /// computes the chaosCoeffs via linear regression
  /// (expCoeffsSolnApproach is REGRESSION)
  void regression();
  /// computes the chaosCoeffs via averaging of samples
  /// (expCoeffsSolnApproach is SAMPLING)
  void expectation();

  /// cross-validates alternate gradient expressions
  void gradient_check();

  //
  //- Heading: Data
  //

  /// number of terms in orthogonal polynomial expansion (length of chaosCoeffs)
  int numExpansionTerms;
  /// order of orthogonal polynomial expansion
  UShortArray approxOrder;

  /// array of basis types for each one-dimensional orthogonal polynomial:
  /// HERMITE, LEGENDRE, LAGUERRE, JACOBI, GENERALIZED_LAGUERRE, CHEBYSHEV,
  /// or NUMERICALLY_GENERATED
  ShortArray basisTypes;

  /// array of Gauss mode options for some derived orthogonal polynomial 
  /// types: Legendre supports GAUSS_LEGENDRE or GAUSS_PATTERSON,
  /// Chebyshev supports CLENSHAW_CURTIS or FEJER2.
  ShortArray gaussModes;

  /// array of one-dimensional basis polynomial objects which are used in
  /// constructing the multivariate orthogonal/interpolation polynomials
  std::vector<BasisPolynomial> polynomialBasis;

  /// numExpansionTerms-by-numVars array for identifying the orders of
  /// the one-dimensional orthogonal polynomials contributing to each
  /// of the multivariate orthogonal polynomials
  UShort2DArray multiIndex;

  /// norm-squared of one of the multivariate polynomial basis functions
  Real multiPolyNormSq;

  /// switch for formulation of orthogonal polynomial expansion
  /// integrated with tensor-product quadrature: TOTAL_ORDER or
  /// TENSOR_PRODUCT expansion.
  short quadratureExpansion;
  /// switch for formulation of orthogonal polynomial expansion
  /// integrated with sparse grids: TOTAL_ORDER,
  /// HEURISTIC_TOTAL_ORDER, or TENSOR_PRODUCT_SUM expansion.
  short sparseGridExpansion;
};


inline OrthogPolyApproximation::
OrthogPolyApproximation(const UShortArray& approx_order, size_t num_vars):
  PolynomialApproximation(num_vars), numExpansionTerms(0),
  approxOrder(approx_order), quadratureExpansion(TENSOR_PRODUCT) //TOTAL_ORDER
  //sparseGridExpansion()//TOTAL_ORDER,HEURISTIC_TOTAL_ORDER,TENSOR_PRODUCT_SUM
{ }


inline OrthogPolyApproximation::~OrthogPolyApproximation()
{ }


inline void OrthogPolyApproximation::expansion_terms(int exp_terms)
{ numExpansionTerms = exp_terms; }


inline const int& OrthogPolyApproximation::expansion_terms() const
{ return numExpansionTerms; }


inline bool OrthogPolyApproximation::
distribution_types(const ShortArray& u_types,
		   const IntArray& int_rules)
{ return distribution_types(u_types, int_rules, basisTypes, gaussModes); }


inline void OrthogPolyApproximation::distribution_basis()
{ distribution_basis(basisTypes, gaussModes, polynomialBasis); }


inline void OrthogPolyApproximation::
distribution_parameters(const ShortArray& u_types, const DistributionParams& dp)
{ distribution_parameters(u_types, dp, polynomialBasis); }


inline void OrthogPolyApproximation::
distributions(const ShortArray& u_types, const IntArray& int_rules,
	      const DistributionParams& dp,
	      std::vector<BasisPolynomial>& poly_basis,
	      ShortArray& basis_types, ShortArray& gauss_modes)
{
  bool dist_params
    = distribution_types(u_types, int_rules, basis_types, gauss_modes);
  distribution_basis(basis_types, gauss_modes, poly_basis);
  if (dist_params)
    distribution_parameters(u_types, dp, poly_basis);
}


inline void OrthogPolyApproximation::
distributions(const ShortArray& u_types, const IntArray& int_rules,
	      const DistributionParams& dp)
{
  if (u_types.size() != numVars) {
    PCerr << "Error: incoming u_types array length (" << u_types.size()
	  << ") does not match number of variables (" << numVars
	  << ") in OrthogPolyApproximation::distribution_types()." << std::endl;
    abort_handler(-1);
  }
  distributions(u_types, int_rules, dp, polynomialBasis, basisTypes,
		gaussModes);
}


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
    if (basisTypes[i] == NUMERICALLY_GENERATED)
      ((NumericGenOrthogPolynomial*)polynomialBasis[i].polynomial_rep())
	->coefficients_norms_flag(flag);
}


inline Real OrthogPolyApproximation::
multivariate_polynomial(const RealVector& xi, size_t term)
{
  unsigned short order_1d;
  Real mvp = 1.0;
  for (size_t i=0; i<numVars; ++i) {
    order_1d = multiIndex[term][i];
    if (order_1d)
      mvp *= polynomialBasis[i].get_value(xi[i], order_1d);
  }
  return mvp;
}

} // namespace Pecos

#endif

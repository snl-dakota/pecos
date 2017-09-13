#ifndef ORTHOG_POLY_BASIS_HPP
#define ORTHOG_POLY_BASIS_HPP
#include <PolynomialApproximation.hpp>

namespace Surrogates {
/**
\class OrthogonalPolynomialBasis()
\brief A multivariate orthogonal polynomial basis

This class was requires the PECOS library.
*/
class OrthogonalPolynomialBasis(){
 public:
  boost::shared_ptr<Surrogates::AleatoryVariableTransform> varTransform_;
  std::vector<BasisPolynomial> polynomialBasis_;
  bool nestedRules_;

  /// Constructor
  OrthogonalPolynomialBasis();

  /// Destructor
  virtual ~OrthogonalPolynomialBasis();

  /**\brief Get the orthogonal basis types from an aleatory variable type.
   */
  static short get_basis_type(short var_type);

  /**\brief Get the orthogonal basis types from the aleatory transformation.
   */
  static void get_basis_types(const boost::shared_ptr<Surrogates::AleatoryVariables> vars,
			      ShortArray& basis_types);

  /**\brief Initialize the polynomials to be orthogonal to the aleatory
     variables
   */
  static void initialize_polynomial_basis(const boost::shared_ptr<Surrogates::AleatoryVariables> vars, std::vector<BasisPolynomial> &polynomial_basis);

}; // class OrthogonalPolynomialBasis

}; // namespace Surrogates
#endif // ORTHOG_POLY_BASIS_HPP

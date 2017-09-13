#ifndef POLYNOMIAL_CHAOS_EXPANSION_HPP
#define POLYNOMIAL_CHAOS_EXPANSION_HPP
#include <PolynomialApproximation.hpp>

namespace Surrogates {
/**
\class PolynomialChaosExpansion
\brief A multivariate polynomial chaos expansion approximation.

This class was requires the PECOS library.
*/
class PolynomialChaosExpansion : public PolynomialApproximation {
public:
  OrthogonalPolynomialBasis orthogPolyBasis_;
  boost::shared_ptr<Surrogates::AleatoryVariableTransformation> aleatoryVarTransform_;

  PolynomialChaosExpansion();

  ~PolynomialChaosExpansion();

  void set_options(const OptionsList &opts);

  /** \copydoc PolynomialApproximation::generate_basis_matrix() */
  void generate_basis_matrix(const RealMatrix &samples, RealMatrix &values_out);

  /** \copydoc PolynomialApproximation::set_variable_transformation() */
  void set_variable_transformation(const boost::shared_ptr<Surrogates::VariableTransform> var_transform);

}; // class PolynomialChaosExpansion

}; // namespace Surrogates
#endif // POLYNOMIAL_CHAOS_EXPANSION_HPP

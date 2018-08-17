#ifndef POLYNOMIAL_CHAOS_EXPANSION_HPP
#define POLYNOMIAL_CHAOS_EXPANSION_HPP
#include "PolyApproximation.hpp"
#include "OrthogonalPolynomialBasis.hpp"
#include "AleatoryVariableTransformation.hpp"

namespace Surrogates {
/**
\class PolynomialChaosExpansion
\brief A multivariate polynomial chaos expansion approximation.

This class was requires the PECOS library.
*/
class PolynomialChaosExpansion : public PolyApproximation {
private:
  OrthogonalPolynomialBasis orthogPolyBasis_;
public:

  //std::shared_ptr<Surrogates::AleatoryVariableTransformation> aleatoryVarTransform_;

  PolynomialChaosExpansion();

  ~PolynomialChaosExpansion();

  void set_options(const OptionsList &opts);

  /** \copydoc PolyApproximation::generate_canonical_basis_matrix() */
  void generate_canonical_basis_matrix(const RealMatrix &samples, RealMatrix &result_0);

  /** \copydoc PolyApproximation::set_variable_transformation() */
  void set_variable_transformation(const std::shared_ptr<Surrogates::VariableTransformation> var_transform);

  void initialize_polynomial_basis_from_basis_types(const Pecos::ShortArray &basis_types);

}; // class PolynomialChaosExpansion

} // namespace Surrogates
#endif // POLYNOMIAL_CHAOS_EXPANSION_HPP

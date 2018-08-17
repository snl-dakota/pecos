#include "PolynomialChaosExpansion.hpp"
namespace Surrogates {

PolynomialChaosExpansion::PolynomialChaosExpansion(){}

PolynomialChaosExpansion::~PolynomialChaosExpansion(){}

void PolynomialChaosExpansion::set_options(const OptionsList &opts){
  PolyApproximation::set_options(opts);
}

void PolynomialChaosExpansion::
generate_canonical_basis_matrix(const RealMatrix &samples, RealMatrix &basis_matrix){
  orthogPolyBasis_.value(samples, basisIndices_, basis_matrix);
}

void PolynomialChaosExpansion::
initialize_polynomial_basis_from_basis_types(const Pecos::ShortArray &basis_types){
  orthogPolyBasis_.initialize_polynomial_basis_from_basis_types(basis_types);
}

void PolynomialChaosExpansion::
set_variable_transformation(const std::shared_ptr<Surrogates::VariableTransformation>var_transform){
  PolyApproximation::set_variable_transformation(var_transform);
  // Need the above line but also below eventually
  // aleatoryVarTransform_ =
  //   Teuchos::rcp_dynamic_cast<Surrogates::AleatoryVariableTransform>(varTransform,true);
  // if (aleatoryVarTransform_.is_null())
  //   throw(std::runtime_error("var_transform is not an object of type AleatoryVariableTransform"));

  // initialize_polynomial_basis();
}

// void PolynomialChaosExpansion::
// initialize_polynomial_basis(){
//   if (aleatoryVarTransform_.is_null())
//     throw(std::runtime_error("Aleatory Variable transform has not been set"));

//   std::shared_ptr<Surrogates::AleatoryVariables> transformed_vars;
//   aleatoryVariableTransform_.get_transformed_varables(transformed_vars)
//   orthogPolyBasis_.initialize_polynomial_basis(transformed_vars);
// }

} // namespace Surrogates

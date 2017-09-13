#include <PolynomialChaosExpansion.hpp>
#include <OrthogPolynomialApproximation.hpp>
namespace Surrogates {

PolynomialChaosExpansion::PolynomialChaosExpansion(){}

PolynomialChaosExpansion::~PolynomialChaosExpansion(){}

void PolynomialChaosExpansion::set_options(const OptionsList &opts){
  PolynomialApproximation::set_options(opts);
}

void PolynomialChaosExpansion::
generate_basis_matrix(const RealMatrix &samples, RealMatrix &basis_matrix){
  //todo: Replace this function with one that takes advantage of recursion
  //formulae of orthogonal polynomials
  size_t i, j, num_exp_terms = multi_index.size(), num_samples = x.numCols(),
    num_vars = x.numRows();
  basis_matrix.shapeUninitialized(num_samples,num_exp_terms);
  for (j=0; j<num_exp_terms; ++j){
    for (i=0; i<num_samples; ++i){
      RealVector sample( Teuchos::View, const_cast<double*>(x[i]), num_vars);
      basis_matrix(i,j) =
	SharedOrthogPolyApproxData::multivariate_polynomial(sample,
							    multi_index[j],
							    polynomial_basis);
    }
  }
}

void PolynomialChaosExpansion::
set_variable_transformation(const boost::shared_ptr<Surrogates::VariableTransformation>var_transform){
  PolynomialApproximation::set_variable_transformation(var_transform);
  aleatoryVarTransform_ =
    Teuchos::rcp_dynamic_cast<Surrogates::AleatoryVariableTransform>(varTransform,true);
  if (aleatoryVarTransform_.is_null())
    throw(std::runtime_error("var_transform is not an object of type AleatoryVariableTransform"));

  initialize_polynomial_basis();
}

void PolynomialChaosExpansion::initialize_polynomial_basis(){
  if (aleatoryVarTransform_.is_null())
    throw(std::runtime_error("Aleatory Variable transform has not been set"));

  boost::shared_ptr<Surrogates::AleatoryVariables> transformed_vars;
  aleatoryVariableTransform_.get_transformed_varables(transformed_vars)
  orthogPolyBasis_.initialize_polynomial_basis(transformed_vars);
}

} // namespace Surrogates

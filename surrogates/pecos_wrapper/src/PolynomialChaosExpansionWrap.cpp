#include "PolynomialChaosExpansionWrap.hpp"

namespace Surrogates {

PolynomialChaosExpansionWrap::PolynomialChaosExpansionWrap(){
  size_t num_vars=2, degree=2;
  short basis_type=Pecos::LEGENDRE_ORTHOG;
  Pecos::UShortArray approx_order(num_vars,degree);
  sharedData_=Pecos::SharedBasisApproxData(basis_type,approx_order,num_vars);
  poly_ = Pecos::BasisApproximation(sharedData_);
};

PolynomialChaosExpansionWrap::~PolynomialChaosExpansionWrap(){}

void PolynomialChaosExpansionWrap::
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
        data_rep->multivariate_polynomial(sample, basisIndices_(d,j);
    }
  }
}

}// namespace Surrogates

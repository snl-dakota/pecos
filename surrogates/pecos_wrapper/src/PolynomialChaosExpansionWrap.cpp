#include "PolynomialChaosExpansionWrap.hpp"
#include "RegressOrthogPolyApproximation.hpp"

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
  Pecos::RegressOrthogPolyApproximation* poly_rep =
    (Pecos::RegressOrthogPolyApproximation*) (&poly_);
  std::cout << "HERE\n";
  poly_rep->build_linear_system(basis_matrix);
}

}// namespace Surrogates

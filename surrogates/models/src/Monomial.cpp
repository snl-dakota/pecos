#include "Monomial.hpp"

namespace Surrogates {

Monomial::Monomial(){}

Monomial::~Monomial(){}

void Monomial::set_options(const OptionsList &opts){
  PolynomialApproximation::set_options(opts);
}

void Monomial::
generate_basis_matrix(const RealMatrix &samples, RealMatrix &basis_matrix){
  int num_vars = varTransform_->num_vars();
  resize_if_needed(basis_matrix, samples.numCols(), basisIndices_.numCols());
  for(int j=0; j<basisIndices_.numCols(); ++j){
    for(int i=0; i<samples.numCols(); ++i){
      Real val = 1.;
      for (int d=0; d<num_vars; d++)
	val *= std::pow(samples(d,i), basisIndices_(d,j));
      basis_matrix(i,j) = val;
    }
  }
}

} // namespace Surrogates

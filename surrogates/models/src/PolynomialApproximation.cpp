#include <PolynomialApproximation.hpp>
#include <teuchos_data_types.hpp>

namespace Surrogates {

PolynomialApproximation::PolynomialApproximation(){};

PolynomialApproximation::~PolynomialApproximation(){};

void PolynomialApproximation::value(const RealMatrix &samples, RealMatrix &approx_vals) {
  if (basisCoeffs_.numRows()!=basisIndices_.numCols())
    throw(std::runtime_error("basis indices and coefficients are inconsistent"));
  resize_if_needed(approx_vals,samples.numCols(),basisCoeffs_.numCols());
  RealMatrix basis_matrix;
  generate_basis_matrix(samples, basis_matrix);
  multiply(basis_matrix,basisCoeffs_,approx_vals,1.0,0.0);
}


void PolynomialApproximation::set_coefficients(const RealMatrix &coeffs) {
    basisCoeffs_.shapeUninitialized(coeffs.numRows(), coeffs.numCols());
    basisCoeffs_.assign(coeffs);  // TODO: compare dimensions to basisIndices?
}

void PolynomialApproximation::get_coefficients(RealMatrix &coeffs) const{
  coeffs = basisCoeffs_;
}

void PolynomialApproximation::set_basis_indices(const IntMatrix &basis_indices){
  if (num_vars()!=basis_indices.numRows())
    throw(std::runtime_error("basis indices must be num_vars x num_indices"));
  resize_if_needed(basisIndices_, basis_indices.numRows(),
		   basis_indices.numCols());
  basisIndices_.assign(basis_indices);
};

void PolynomialApproximation:: get_basis_indices(IntMatrix &basis_indices) const{
  basis_indices = basisIndices_;
}

void PolynomialApproximation::gradient(const RealMatrix &samples, int qoi, RealMatrix &gradients) {
  throw(std::string("This Function type does not provide gradients"));
}

void PolynomialApproximation::jacobian(const RealVector &sample, RealMatrix &jacobian) {
  throw(std::string("This Function type does not provide a jacobian"));
}

void PolynomialApproximation::hessian(const RealMatrix &samples, int qoi, RealMatrixList &hessians) {
  throw(std::string("This Function type does not provide gradients"));
}

} // namespace Surrogates

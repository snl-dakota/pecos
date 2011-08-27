/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        NodalInterpPolyApproximation
//- Description:  Class for Nodal Interpolation Polynomial Approximation
//-               
//- Owner:        Mike Eldred

#ifndef NODAL_INTERP_POLY_APPROXIMATION_HPP
#define NODAL_INTERP_POLY_APPROXIMATION_HPP

#include "InterpPolyApproximation.hpp"


namespace Pecos {

/// Derived approximation class for nodal interpolation polynomials
/// (global approximation interpolating function values and
/// potentially gradients at collocation points).

/** The NodalInterpPolyApproximation class provides a global polynomial
    approximation based on either Lagrange or Hermite interpolation
    polynomials using a nodal basis approach.  It is used primarily
    for stochastic collocation approaches to uncertainty quantification. */

class NodalInterpPolyApproximation: public InterpPolyApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  NodalInterpPolyApproximation(short basis_type, size_t num_vars,
			       bool use_derivs);
  /// destructor
  ~NodalInterpPolyApproximation();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// retrieve the response value for the current expansion using the
  /// given parameter vector
  Real value(const RealVector& x);
  /// retrieve the response gradient for the current expansion using
  /// the given parameter vector and default DVV
  const RealVector& gradient(const RealVector& x);
  /// retrieve the response gradient for the current expansion using
  /// the given parameter vector and given DVV
  const RealVector& gradient(const RealVector& x, const SizetArray& dvv);

  Real stored_value(const RealVector& x);
  const RealVector& stored_gradient(const RealVector& x);

  Real mean();
  Real mean(const RealVector& x);
  const RealVector& mean_gradient();
  const RealVector& mean_gradient(const RealVector& x, const SizetArray& dvv);

  Real variance();
  Real variance(const RealVector& x);
  const RealVector& variance_gradient();
  const RealVector& variance_gradient(const RealVector& x,
				      const SizetArray& dvv);

  Real covariance(PolynomialApproximation* poly_approx_2);
  Real covariance(const RealVector& x, PolynomialApproximation* poly_approx_2);

private:

  //
  //- Heading: Convenience functions
  //

  /// return value of type 1 interpolation polynomial (TPQ case)
  Real type1_interpolant_value(const RealVector& x, const UShortArray& key);
  /// return value of type 1 interpolation polynomial (SSG case)
  Real type1_interpolant_value(const RealVector& x, const UShortArray& key,
			       const UShortArray& sm_index);

  /// return gradient of type 1 interpolation polynomial (TPQ case)
  Real type1_interpolant_gradient(const RealVector& x, size_t deriv_index,
				  const UShortArray& key);
  /// return gradient of type 1 interpolation polynomial (SSG case)
  Real type1_interpolant_gradient(const RealVector& x, size_t deriv_index,
				  const UShortArray& key,
				  const UShortArray& sm_index);

  /// return value of type 2 interpolation polynomial (TPQ case)
  Real type2_interpolant_value(const RealVector& x, size_t interp_index,
			       const UShortArray& key);
  /// return value of type 2 interpolation polynomial (SSG case)
  Real type2_interpolant_value(const RealVector& x, size_t interp_index,
			       const UShortArray& key,
			       const UShortArray& sm_index);

  /// return gradient of type 2 interpolation polynomial (TPQ case)
  Real type2_interpolant_gradient(const RealVector& x, size_t deriv_index,
				  size_t interp_index, const UShortArray& key);
  /// return gradient of type 2 interpolation polynomial (SSG case)
  Real type2_interpolant_gradient(const RealVector& x, size_t deriv_index,
				  size_t interp_index, const UShortArray& key,
				  const UShortArray& sm_index);

  /// compute the value of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to value(x)
  Real tensor_product_value(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs,
    const UShort2DArray& key);
  /// compute the value of a sparse interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to value(x)
  Real tensor_product_value(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs,
    const UShortArray& sm_index, const UShort2DArray& key,
    const SizetArray& colloc_index);

  /// compute the gradient of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to gradient(x)
  const RealVector& tensor_product_gradient(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs,
    const UShort2DArray& key);
  /// compute the gradient of a sparse interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to gradient(x)
  const RealVector& tensor_product_gradient(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs,
    const UShortArray& sm_index, const UShort2DArray& key,
    const SizetArray& colloc_index);
  /// compute the gradient of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid for given DVV; contributes to gradient(x, dvv)
  const RealVector& tensor_product_gradient(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs,
    const UShort2DArray& key, const SizetArray& dvv);
  /// compute the gradient of a sparse interpolant on an isotropic/anisotropic
  /// tensor-product grid for given DVV; contributes to gradient(x, dvv)
  const RealVector& tensor_product_gradient(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs,
    const UShortArray& sm_index, const UShort2DArray& key,
    const SizetArray& colloc_index, const SizetArray& dvv);

  /// compute the mean of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to mean(x)
  Real tensor_product_mean(const RealVector& x, const UShort2DArray& key);
  /// compute the mean of a sparse interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to mean(x)
  Real tensor_product_mean(const RealVector& x, const UShortArray& sm_index,
    const UShort2DArray& key, const SizetArray& colloc_index);

  /// compute the mean of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to mean(x)
  const RealVector& tensor_product_mean_gradient(const RealVector& x,
    const UShort2DArray& key, const SizetArray& dvv);
  /// compute the mean of a sparse interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to mean(x)
  const RealVector& tensor_product_mean_gradient(const RealVector& x,
    const UShortArray& sm_index,    const UShort2DArray& key,
    const SizetArray& colloc_index, const SizetArray& dvv);

  /// compute the covariance of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to covariance(x, poly_approx_2)
  Real tensor_product_covariance(const RealVector& x, const UShort2DArray& key,
				 const RealVector& exp_coeffs_2);
  /// compute the covariance of a sparse interpolant on an isotropic/anisotropic
  /// sparse grid; contributes to covariance(x, poly_approx_2)
  Real tensor_product_covariance(const RealVector& x,
    const UShortArray& sm_index,    const UShort2DArray& key,
    const SizetArray& colloc_index, const RealVector& exp_coeffs_2);

  /// compute the variance of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to variance(x)
  const RealVector& tensor_product_variance_gradient(const RealVector& x,
    const UShort2DArray& key, const SizetArray& dvv);
  /// compute the variance of a sparse interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to variance(x)
  const RealVector& tensor_product_variance_gradient(const RealVector& x,
    const UShortArray& sm_index, const UShort2DArray& key,
    const SizetArray& colloc_index, const SizetArray& dvv);

  //
  //- Heading: Data
  //

};


inline NodalInterpPolyApproximation::
NodalInterpPolyApproximation(short basis_type, size_t num_vars,
			     bool use_derivs):
  InterpPolyApproximation(basis_type, num_vars, use_derivs)
{ }


inline NodalInterpPolyApproximation::~NodalInterpPolyApproximation()
{ }


inline Real NodalInterpPolyApproximation::
type1_interpolant_value(const RealVector& x, const UShortArray& key)
{
  Real L1 = 1.;
  std::vector<BasisPolynomial>& basis_0 = polynomialBasis[0];
  for (size_t j=0; j<numVars; ++j)
    L1 *= basis_0[j].type1_value(x[j], key[j]);
  return L1;
}


inline Real NodalInterpPolyApproximation::
type1_interpolant_value(const RealVector& x, const UShortArray& key,
			const UShortArray& sm_index)
{
  Real L1 = 1.;
  for (size_t j=0; j<numVars; ++j)
    L1 *= polynomialBasis[sm_index[j]][j].type1_value(x[j], key[j]);
  return L1;
}


inline Real NodalInterpPolyApproximation::
type1_interpolant_gradient(const RealVector& x, size_t deriv_index,
			   const UShortArray& key)
{
  Real L1_grad = 1.;
  std::vector<BasisPolynomial>& basis_0 = polynomialBasis[0];
  for (size_t k=0; k<numVars; ++k)
    L1_grad *= (k == deriv_index) ? basis_0[k].type1_gradient(x[k], key[k]) :
                                    basis_0[k].type1_value(x[k],    key[k]);
  return L1_grad;
}


inline Real NodalInterpPolyApproximation::
type1_interpolant_gradient(const RealVector& x,    size_t deriv_index,
			   const UShortArray& key, const UShortArray& sm_index)
{
  Real L1_grad = 1.;
  for (size_t k=0; k<numVars; ++k)
    L1_grad *= (k == deriv_index) ?
      polynomialBasis[sm_index[k]][k].type1_gradient(x[k], key[k]) :
      polynomialBasis[sm_index[k]][k].type1_value(x[k],    key[k]);
  return L1_grad;
}


inline Real NodalInterpPolyApproximation::
type2_interpolant_value(const RealVector& x, size_t interp_index,
			const UShortArray& key)
{
  Real L2 = 1.;
  std::vector<BasisPolynomial>& basis_0 = polynomialBasis[0];
  for (size_t k=0; k<numVars; ++k)
    L2 *= (interp_index == k) ? basis_0[k].type2_value(x[k], key[k]) :
                                basis_0[k].type1_value(x[k], key[k]);
  return L2;
}


inline Real NodalInterpPolyApproximation::
type2_interpolant_value(const RealVector& x,    size_t interp_index,
			const UShortArray& key, const UShortArray& sm_index)
{
  Real L2 = 1.;
  for (size_t k=0; k<numVars; ++k)
    L2 *= (interp_index == k) ?
      polynomialBasis[sm_index[k]][k].type2_value(x[k], key[k]) :
      polynomialBasis[sm_index[k]][k].type1_value(x[k], key[k]);
  return L2;
}


inline Real NodalInterpPolyApproximation::
type2_interpolant_gradient(const RealVector& x, size_t deriv_index,
			   size_t interp_index, const UShortArray& key)
{
  // deriv_index  = desired gradient component
  // interp_index = index of gradient component used in type2 interpolation
  Real L2_grad = 1.;
  std::vector<BasisPolynomial>& basis_0 = polynomialBasis[0];
  if (interp_index == deriv_index) // match in grad and type2 interp components
    for (size_t l=0; l<numVars; ++l)
      L2_grad *= (l == interp_index) ? basis_0[l].type2_gradient(x[l], key[l]) :
                                       basis_0[l].type1_value(x[l],    key[l]);
  else                          // mismatch in grad and type2 interp components
    for (size_t l=0; l<numVars; ++l)
      L2_grad *= (l == interp_index) ? basis_0[l].type2_value(x[l],    key[l]) :
	                               basis_0[l].type1_gradient(x[l], key[l]);
  return L2_grad;
}


inline Real NodalInterpPolyApproximation::
type2_interpolant_gradient(const RealVector& x, size_t deriv_index,
			   size_t interp_index, const UShortArray& key,
			   const UShortArray& sm_index)
{
  // deriv_index  = desired gradient component
  // interp_index = index of gradient component used in type2 interpolation
  Real L2_grad = 1.;
  if (interp_index == deriv_index) // match in grad and type2 interp components
    for (size_t l=0; l<numVars; ++l)
      L2_grad *= (l == interp_index) ?
	polynomialBasis[sm_index[l]][l].type2_gradient(x[l], key[l]) :
	polynomialBasis[sm_index[l]][l].type1_value(x[l],    key[l]);
  else                          // mismatch in grad and type2 interp components
    for (size_t l=0; l<numVars; ++l)
      L2_grad *= (l == interp_index) ?
	polynomialBasis[sm_index[l]][l].type2_value(x[l],    key[l]) :
	polynomialBasis[sm_index[l]][l].type1_gradient(x[l], key[l]);
  return L2_grad;
}

} // namespace Pecos

#endif

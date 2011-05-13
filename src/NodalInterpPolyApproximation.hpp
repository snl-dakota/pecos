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

  /// retrieve the response expansion value for a given parameter vector
  const Real& value(const RealVector& x);
  /// retrieve the response expansion gradient for a given parameter vector
  /// and default DVV
  const RealVector& gradient(const RealVector& x);
  /// retrieve the response expansion gradient for a given parameter vector
  /// and given DVV
  const RealVector& gradient(const RealVector& x, const SizetArray& dvv);

  /// return the mean of the expansion, treating all variables as random
  const Real& mean();
  /// return the mean of the expansion for a given parameter vector,
  /// treating a subset of the variables as random
  const Real& mean(const RealVector& x);
  /// return the gradient of the expansion mean for a given parameter vector,
  /// treating all variables as random
  const RealVector& mean_gradient();
  /// return the gradient of the expansion mean for a given parameter vector
  /// and given DVV, treating a subset of the variables as random
  const RealVector& mean_gradient(const RealVector& x, const SizetArray& dvv);

  /// return the variance of the expansion, treating all variables as random
  const Real& variance();
  /// return the variance of the expansion for a given parameter vector,
  /// treating a subset of the variables as random
  const Real& variance(const RealVector& x);
  /// return the gradient of the expansion variance for a given parameter
  /// vector, treating all variables as random
  const RealVector& variance_gradient();
  /// return the gradient of the expansion variance for a given parameter
  /// vector and given DVV, treating a subset of the variables as random
  const RealVector& variance_gradient(const RealVector& x,
				      const SizetArray& dvv);

  /// return the covariance of the expansion, treating all variables as random
  Real covariance(PolynomialApproximation* poly_approx_2);
  /// return the covariance of the expansion for a given parameter vector,
  /// treating a subset of the variables as random
  Real covariance(const RealVector& x, PolynomialApproximation* poly_approx_2);

private:

  //
  //- Heading: Convenience functions
  //

  /// compute the value of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to value(x)
  const Real& tensor_product_value(const RealVector& x);
  /// compute the value of a sparse interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to value(x)
  const Real& tensor_product_value(const RealVector& x, size_t tp_index);

  /// compute the gradient of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to gradient(x)
  const RealVector& tensor_product_gradient(const RealVector& x);
  /// compute the gradient of a sparse interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to gradient(x)
  const RealVector& tensor_product_gradient(const RealVector& x,
					    size_t tp_index);
  /// compute the gradient of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid for given DVV; contributes to gradient(x, dvv)
  const RealVector& tensor_product_gradient(const RealVector& x,
					    const SizetArray& dvv);
  /// compute the gradient of a sparse interpolant on an isotropic/anisotropic
  /// tensor-product grid for given DVV; contributes to gradient(x, dvv)
  const RealVector& tensor_product_gradient(const RealVector& x,
    size_t tp_index, const SizetArray& dvv);

  /// compute the mean of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to mean(x)
  const Real& tensor_product_mean(const RealVector& x);
  /// compute the mean of a sparse interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to mean(x)
  const Real& tensor_product_mean(const RealVector& x, size_t tp_index);

  /// compute the mean of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to mean(x)
  const RealVector& tensor_product_mean_gradient(const RealVector& x,
						 const SizetArray& dvv);
  /// compute the mean of a sparse interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to mean(x)
  const RealVector& tensor_product_mean_gradient(const RealVector& x,
    size_t tp_index, const SizetArray& dvv);

  /// compute the covariance of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to covariance(x, poly_approx_2)
  const Real& tensor_product_covariance(const RealVector& x,
					const RealVector& exp_coeffs_2);
  /// compute the covariance of a sparse interpolant on an isotropic/anisotropic
  /// sparse grid; contributes to covariance(x, poly_approx_2)
  const Real& tensor_product_covariance(const RealVector& x,
					const RealVector& exp_coeffs_2,
					size_t tp_index);

  /// compute the variance of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to variance(x)
  const RealVector& tensor_product_variance_gradient(const RealVector& x,
						     const SizetArray& dvv);
  /// compute the variance of a sparse interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to variance(x)
  const RealVector& tensor_product_variance_gradient(const RealVector& x,
    size_t tp_index, const SizetArray& dvv);

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

} // namespace Pecos

#endif

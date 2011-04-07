/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        LagrangeInterpPolyApproximation
//- Description:  Class for Lagrange Interpolation Polynomial Approximation
//-               
//- Owner:        Mike Eldred

#ifndef LAGRANGE_INTERP_POLY_APPROXIMATION_HPP
#define LAGRANGE_INTERP_POLY_APPROXIMATION_HPP

#include "InterpPolyApproximation.hpp"


namespace Pecos {

/// Derived approximation class for Lagrange interpolation polynomials (global
/// approximation interpolating function values at collocation points).

/** The LagrangeInterpPolyApproximation class provides a global approximation
    based on Lagrange interpolation polynomials.  It is used primarily for
    stochastic collocation approaches to uncertainty quantification. */

class LagrangeInterpPolyApproximation: public InterpPolyApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  LagrangeInterpPolyApproximation(short basis_type, size_t num_vars);
  /// destructor
  ~LagrangeInterpPolyApproximation();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// retrieve the response expansion value for a given parameter vector
  const Real& get_value(const RealVector& x);
  /// retrieve the response expansion gradient for a given parameter vector
  /// and default DVV
  const RealVector& get_gradient(const RealVector& x);
  /// retrieve the response expansion gradient for a given parameter vector
  /// and given DVV
  const RealVector& get_gradient(const RealVector& x, const SizetArray& dvv);

  /// return the mean of the expansion, treating all variables as random
  const Real& get_mean();
  /// return the mean of the expansion for a given parameter vector,
  /// treating a subset of the variables as random
  const Real& get_mean(const RealVector& x);
  /// return the gradient of the expansion mean for a given parameter vector,
  /// treating all variables as random
  const RealVector& get_mean_gradient();
  /// return the gradient of the expansion mean for a given parameter vector
  /// and given DVV, treating a subset of the variables as random
  const RealVector& get_mean_gradient(const RealVector& x,
				      const SizetArray& dvv);

  /// return the variance of the expansion, treating all variables as random
  const Real& get_variance();
  /// return the variance of the expansion for a given parameter vector,
  /// treating a subset of the variables as random
  const Real& get_variance(const RealVector& x);
  /// return the gradient of the expansion variance for a given parameter
  /// vector, treating all variables as random
  const RealVector& get_variance_gradient();
  /// return the gradient of the expansion variance for a given parameter
  /// vector and given DVV, treating a subset of the variables as random
  const RealVector& get_variance_gradient(const RealVector& x,
					  const SizetArray& dvv);

  /// return the covariance of the expansion, treating all variables as random
  Real get_covariance(const RealVector& exp_coeffs_2);
  /// return the covariance of the expansion for a given parameter vector,
  /// treating a subset of the variables as random
  Real get_covariance(const RealVector& x, const RealVector& exp_coeffs_2);

private:

  //
  //- Heading: Convenience functions
  //

  /// compute the value of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to get_value(x)
  const Real& tensor_product_value(const RealVector& x);
  /// compute the value of a sparse interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to get_value(x)
  const Real& tensor_product_value(const RealVector& x, size_t tp_index);

  /// compute the gradient of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to get_gradient(x)
  const RealVector& tensor_product_gradient(const RealVector& x);
  /// compute the gradient of a sparse interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to get_gradient(x)
  const RealVector& tensor_product_gradient(const RealVector& x,
					    size_t tp_index);
  /// compute the gradient of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid for given DVV; contributes to get_gradient(x, dvv)
  const RealVector& tensor_product_gradient(const RealVector& x,
					    const SizetArray& dvv);
  /// compute the gradient of a sparse interpolant on an isotropic/anisotropic
  /// tensor-product grid for given DVV; contributes to get_gradient(x, dvv)
  const RealVector& tensor_product_gradient(const RealVector& x,
    size_t tp_index, const SizetArray& dvv);

  /// compute the mean of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to get_mean(x)
  const Real& tensor_product_mean(const RealVector& x);
  /// compute the mean of a sparse interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to get_mean(x)
  const Real& tensor_product_mean(const RealVector& x, size_t tp_index);

  /// compute the mean of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to get_mean(x)
  const RealVector& tensor_product_mean_gradient(const RealVector& x,
						 const SizetArray& dvv);
  /// compute the mean of a sparse interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to get_mean(x)
  const RealVector& tensor_product_mean_gradient(const RealVector& x,
    size_t tp_index, const SizetArray& dvv);

  /// compute the covariance of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to get_covariance(x, exp_coeffs_2)
  const Real& tensor_product_covariance(const RealVector& x,
					const RealVector& exp_coeffs_2);
  /// compute the covariance of a sparse interpolant on an isotropic/anisotropic
  /// sparse grid; contributes to get_covariance(x, exp_coeffs_2)
  const Real& tensor_product_covariance(const RealVector& x,
					const RealVector& exp_coeffs_2,
					size_t tp_index);

  /// compute the variance of a tensor interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to get_variance(x)
  const RealVector& tensor_product_variance_gradient(const RealVector& x,
						     const SizetArray& dvv);
  /// compute the variance of a sparse interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to get_variance(x)
  const RealVector& tensor_product_variance_gradient(const RealVector& x,
    size_t tp_index, const SizetArray& dvv);

  //
  //- Heading: Data
  //
};


inline LagrangeInterpPolyApproximation::
LagrangeInterpPolyApproximation(short basis_type, size_t num_vars):
  InterpPolyApproximation(basis_type, num_vars)
{ useDerivs = false; }


inline LagrangeInterpPolyApproximation::~LagrangeInterpPolyApproximation()
{ }

} // namespace Pecos

#endif

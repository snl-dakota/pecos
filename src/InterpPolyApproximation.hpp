/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        InterpPolyApproximation
//- Description:  Class for Lagrange Interpolation Polynomial Approximation
//-               
//- Owner:        Mike Eldred

#ifndef INTERP_POLY_APPROXIMATION_HPP
#define INTERP_POLY_APPROXIMATION_HPP

#include "PolynomialApproximation.hpp"
#include "BasisPolynomial.hpp"


namespace Pecos {

/// Derived approximation class for interpolation polynomials (global
/// approximation).

/** The InterpPolyApproximation class provides a global approximation
    based on interpolation polynomials.  It is used primarily for
    stochastic collocation approaches to uncertainty quantification. */

class InterpPolyApproximation: public PolynomialApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  InterpPolyApproximation(size_t num_vars, short output_level);
  /// destructor
  ~InterpPolyApproximation();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  int min_coefficients() const;

  /// find the coefficients for the expansion of multivariate
  /// interpolation polynomials
  void find_coefficients();

  /// size expansionCoeffs and expansionCoeffGrads
  void allocate_arrays();

  /// Performs global sensitivity analysis using Sobol' Indices
  void compute_global_sensitivity();

  /// retrieve the response expansion value for a given parameter vector
  const Real& get_value(const RealVector& x);
  /// retrieve the response expansion gradient for a given parameter vector
  /// and default DVV
  const RealVector& get_gradient(const RealVector& x);
  /// retrieve the response expansion gradient for a given parameter vector
  /// and given DVV
  const RealVector& get_gradient(const RealVector& x, const UIntArray& dvv);

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
				      const UIntArray& dvv);

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
					  const UIntArray& dvv);

  /// return the covariance of the expansion, treating all variables as random
  const Real& get_covariance(const RealVector& exp_coeffs_2);
  // return the covariance of the expansion for a given parameter vector,
  // treating a subset of the variables as random
  //const Real& get_covariance(const RealVector& x,
  //                           const RealVector& exp_coeffs_2);

private:

  //
  //- Heading: Convenience functions
  //

  /// compute the value of an interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to get_value(x)
  const Real& tensor_product_value(const RealVector& x, size_t tp_index);
  /// compute the gradient of an interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to get_gradient(x)
  const RealVector& tensor_product_gradient(const RealVector& x,
					    size_t tp_index);
  /// compute the gradient of an interpolant on an isotropic/anisotropic
  /// tensor-product grid for given DVV; contributes to get_gradient(x, dvv)
  const RealVector& tensor_product_gradient(const RealVector& x,
    size_t tp_index, const UIntArray& dvv);
  /// compute the mean of an interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to get_mean(x)
  const Real& tensor_product_mean(const RealVector& x, size_t tp_index);
  /// compute the mean of an interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to get_mean(x)
  const RealVector& tensor_product_mean_gradient(const RealVector& x,
    size_t tp_index, const UIntArray& dvv);

  /// compute the variance of an interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to get_variance(x)
  const Real& tensor_product_variance(const RealVector& x, size_t tp_index);
  /// compute the variance of an interpolant on an isotropic/anisotropic
  /// tensor-product grid; contributes to get_variance(x)
  const RealVector& tensor_product_variance_gradient(const RealVector& x,
    size_t tp_index, const UIntArray& dvv);

  /// performs sorting to store constituent subsets (constituentSets)
  void get_subsets();
  /// recursively identifies constituent subsets
  void lower_sets(int plus_one_set, IntSet &top_level_set);
  /// finds variance of interpolant with respect to variables in the set
  Real partial_variance_integral(const int &set_value, size_t tp_index,
				 UShortArray &quad_order);
  /// computes partialVariance
  void partial_variance(const int &set_value);

  //
  //- Heading: Data
  //

  /// the constituent subsets for each superset
  std::vector<IntSet> constituentSets;
  /// the partialVariances of subset functions f_alpha
  RealVector partialVariance; 

  /// 2D array of one-dimensional basis polynomial objects which are used in
  /// constructing the multivariate orthogonal/interpolation polynomials.
  /** Each variable (outer array size = numVars) may have multiple integration
      orders associated with it (inner array size = num_levels_per_var = 1 for
      quadrature, w + numVars for sparse grid). */
  std::vector< std::vector< BasisPolynomial > > polynomialBasis;

  /// total number of collocation points = number of terms in interpolation
  /// expansion (length of expansionCoeffs)
  int numCollocPts;

  /// the value of a tensor-product interpolant; a contributor to approxValue
  Real tpValue;
  /// the gradient of a tensor-product interpolant; a contributor to
  /// approxGradient
  RealVector tpGradient;
  /// the mean of a tensor-product interpolant; a contributor to expansionMean
  Real tpMean;
  /// the gradient of the mean of a tensor-product interpolant; a
  /// contributor to expansionMeanGrad
  RealVector tpMeanGrad;
  /// the variance of a tensor-product interpolant; a contributor to
  /// expansionVariance
  Real tpVariance;
  /// the gradient of the variance of a tensor-product interpolant; a
  /// contributor to expansionVarianceGrad 
  RealVector tpVarianceGrad;
};


inline InterpPolyApproximation::
InterpPolyApproximation(size_t num_vars, short output_level):
  PolynomialApproximation(num_vars, output_level), numCollocPts(0)
{ }


inline InterpPolyApproximation::~InterpPolyApproximation()
{ }

} // namespace Pecos

#endif

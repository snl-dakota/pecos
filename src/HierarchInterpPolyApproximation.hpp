/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HierarchInterpPolyApproximation
//- Description:  Class for Nodal Interpolation Polynomial Approximation
//-               
//- Owner:        Chris Miller

#ifndef HIERARCH_INTERP_POLY_APPROXIMATION_HPP
#define HIERARCH_INTERP_POLY_APPROXIMATION_HPP

#include "InterpPolyApproximation.hpp"
#include "LocalRefinableDriver.hpp"

namespace Pecos {

/// Derived approximation class for piecewise linear and cubic hierarchical 
/// interpolation polynomials (local approximation interpolating function values and
/// potentially gradients at collocation points).

/** The HierarchInterpPolyApproximation class provides a local piecewise polynomial
    approximation based on the hierarchical approach described in 
    X. Ma and N. Zabaras "An adaptive hierarchical sparse grid collocation algorithm
    for the solution of stochastic differential equations", Journal of Computational
    Physics, 228 (2009), 3084-3113.  Both piecewise linear basis functions using 
    function values at the collocation points and cubic Hermite basis functions
    using both values and derivatives are available.  It is used primarily
    for stochastic collocation approaches to uncertainty quantification. */

class HierarchInterpPolyApproximation: public InterpPolyApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  HierarchInterpPolyApproximation(short basis_type, 
				  size_t num_vars,
				  bool use_derivs);
  /// destructor
  ~HierarchInterpPolyApproximation();

  /// Compute the basis coefficients from the data.
  virtual void compute_coefficients();

  /// Update the coefficients using new data.
  virtual void increment_coefficients();

  //
  //- Heading: Virtual function redefinitions
  //


  /// retrieve the response expansion value for a given parameter vector
  Real value(const RealVector& x);
  /// retrieve the response expansion gradient for a given parameter vector
  /// and default DVV
  const RealVector& gradient(const RealVector& x);
  /// retrieve the response expansion gradient for a given parameter vector
  /// and given DVV
  const RealVector& gradient(const RealVector& x, const SizetArray& dvv);

  /// return the mean of the expansion, treating all variables as random
  Real mean();
  /// return the mean of the expansion for a given parameter vector,
  /// treating a subset of the variables as random
  Real mean(const RealVector& x);
  /// return the gradient of the expansion mean for a given parameter vector,
  /// treating all variables as random
  const RealVector& mean_gradient();
  /// return the gradient of the expansion mean for a given parameter vector
  /// and given DVV, treating a subset of the variables as random
  const RealVector& mean_gradient(const RealVector& x,
				  const SizetArray& dvv);

  /// return the variance of the expansion, treating all variables as random
  Real variance();
  /// return the variance of the expansion for a given parameter vector,
  /// treating a subset of the variables as random
  Real variance(const RealVector& x);
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
  Real covariance(const RealVector& x,
		      PolynomialApproximation* poly_approx_2);

protected:
  /// returns an int array containing the indices of the points whose support includes x
  const IntArray& in_support_of(const RealVector& x);

  /// compute the value at a point using a lower level than the full approximation.
  Real value(const RealVector& x, unsigned int max_level);

  /// compute the approximate gradient at a point using a lower level than the full approximation.
  const RealVector& gradient(const RealVector& x, unsigned int max_level);

  /// The largest computed coefficient.
  unsigned int maxComputedCoeff;
private:
  
  ///Pecos:PIECEWISE_INTERP_POLYNOMIAL or Pecos:PIECEWISE_CUBIC_INTERP
  short polyType;

  //Array of ints indicating which basis functions have support containing a given point.
  IntArray supportIndicator;

};


} // namespace Pecos

#endif

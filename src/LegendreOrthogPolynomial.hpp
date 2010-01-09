/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        LegendreOrthogPolynomial
//- Description:  Class for Legendre Orthogonal Polynomial
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef LEGENDRE_ORTHOG_POLYNOMIAL_HPP
#define LEGENDRE_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"


namespace Pecos {

/// Derived orthogonal polynomial class for Legendre polynomials

/** The LegendreOrthogPolynomial class evaluates a univariate Legendre
    polynomial of a particular order.  These polynomials are
    orthogonal with respect to the weight function 1 when integrated
    over the support range of [-1,+1].  This corresponds to the
    probability density function f(x) = 1/(U-L) = 1/2 for the uniform
    distribution for [L,U]=[-1,1].  It enables (mixed)
    multidimensional orthogonal polynomial basis functions within
    OrthogPolyApproximation.  Legendre polynomials are a special case
    (alpha = beta = 0) of the more general Jacobi polynomials
    (implemented separately) which correspond to the beta distribution. */

class LegendreOrthogPolynomial: public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  LegendreOrthogPolynomial();  ///< default constructor
  ~LegendreOrthogPolynomial(); ///< destructor

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// retrieve the Legendre polynomial value for a given parameter x 
  const Real& get_value(const Real& x, unsigned short order);
  /// retrieve the Legendre polynomial gradient for a given parameter x 
  const Real& get_gradient(const Real& x, unsigned short order);

  /// return the inner product <P_n,P_n> = ||P_n||^2
  const Real& norm_squared(unsigned short order);

  /// return the Gauss-Legendre quadrature points corresponding to
  /// polynomial order n
  const RealArray& gauss_points(unsigned short order);
  /// return the Gauss-Legendre quadrature weights corresponding to
  /// polynomial order n
  const RealArray& gauss_weights(unsigned short order);

private:

  //
  //- Heading: Data
  //

};


inline LegendreOrthogPolynomial::LegendreOrthogPolynomial()
{ wtFactor = 0.5; }


inline LegendreOrthogPolynomial::~LegendreOrthogPolynomial()
{ }

} // namespace Pecos

#endif

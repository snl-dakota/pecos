/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        LaguerreOrthogPolynomial
//- Description:  Class for Laguerre Orthogonal Polynomial
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef LAGUERRE_ORTHOG_POLYNOMIAL_HPP
#define LAGUERRE_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"


namespace Pecos {

/// Derived orthogonal polynomial class for Laguerre polynomials

/** The LaguerreOrthogPolynomial class evaluates a univariate Laguerre
    polynomial of a particular order.  These polynomials are
    orthogonal with respect to the weight function exp(-x) when
    integrated over the support range of [0,+infinity].  This
    corresponds to the probability density function for the standard
    exponential distribution.  It enables (mixed) multidimensional
    orthogonal polynomial basis functions within
    OrthogPolyApproximation.  Laguerre polynomials are a special case
    (alpha = 0) of the generalized Laguerre polynomials (implemented
    separately) which correspond to the standard gamma distribution. */

class LaguerreOrthogPolynomial: public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  LaguerreOrthogPolynomial();  ///< default constructor
  ~LaguerreOrthogPolynomial(); ///< destructor

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// retrieve the Laguerre polynomial value for a given parameter x 
  const Real& get_value(const Real& x, unsigned short order);
  /// retrieve the Laguerre polynomial gradient for a given parameter x 
  const Real& get_gradient(const Real& x, unsigned short order);

  /// return the inner product <L_n,L_n> = ||L_n||^2
  const Real& norm_squared(unsigned short order);

  /// return the Gauss-Laguerre quadrature points corresponding to
  /// polynomial order n
  const RealArray& collocation_points(unsigned short order);
  /// return the Gauss-Laguerre quadrature weights corresponding to
  /// polynomial order n
  const RealArray& collocation_weights(unsigned short order);

private:

  //
  //- Heading: Data
  //

};


inline LaguerreOrthogPolynomial::LaguerreOrthogPolynomial()
{ collocMode = GAUSS_LAGUERRE; }


inline LaguerreOrthogPolynomial::~LaguerreOrthogPolynomial()
{ }

} // namespace Pecos

#endif

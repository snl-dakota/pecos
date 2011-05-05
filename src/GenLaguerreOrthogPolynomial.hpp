/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        GenLaguerreOrthogPolynomial
//- Description:  Class for GenLaguerre Orthogonal Polynomial
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef GEN_LAGUERRE_ORTHOG_POLYNOMIAL_HPP
#define GEN_LAGUERRE_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"


namespace Pecos {

/// Derived orthogonal polynomial class for generalized Laguerre polynomials

/** The GenLaguerreOrthogPolynomial class evaluates a univariate
    generalized/associated Laguerre polynomial L^(alpha)_n of a
    particular order.  These polynomials are orthogonal with respect
    to the weight function x^alpha exp(-x) when integrated over the
    support range of [0,+infinity].  This corresponds to the
    probability density function f(x) = x^alpha exp(-x) / Gamma(alpha+1)
    for the standard gamma distribution, although common statistical
    PDF parameter conventions (see, e.g., the uncertain variables
    section in the DAKOTA Reference Manual) and the Abramowitz and
    Stegun orthogonal polynomial parameter conventions require an
    offset conversion in this case (alpha_poly = alpha_stat - 1 with
    the poly definition used in both cases above).  It enables (mixed)
    multidimensional orthogonal polynomial basis functions within
    OrthogPolyApproximation.  A special case is the
    LaguerreOrthogPolynomial (implemented separately), for which
    alpha_poly = 0 and weight function = exp(-x) (the standard
    exponential distribution). */

class GenLaguerreOrthogPolynomial: public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  GenLaguerreOrthogPolynomial();                       ///< default constructor
  GenLaguerreOrthogPolynomial(const Real& alpha_stat); ///< standard constructor
  ~GenLaguerreOrthogPolynomial();                      ///< destructor

  //
  //- Heading: Virtual function redefinitions
  //

  /// calculate and return wtFactor based on alphaPoly
  const Real& weight_factor();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// retrieve the generalized Laguerre polynomial value for a given
  /// parameter x
  const Real& get_value(const Real& x, unsigned short order);
  /// retrieve the generalized Laguerre polynomial gradient for a
  /// given parameter x
  const Real& get_gradient(const Real& x, unsigned short order);

  /// return the inner product <L^(alpha)_n,L^(alpha)_n> = ||L^(alpha)_n||^2
  const Real& norm_squared(unsigned short order);

  /// return the generalized Gauss-Laguerre quadrature points
  /// corresponding to polynomial order n
  const RealArray& collocation_points(unsigned short order);
  /// return the generalized Gauss-Laguerre quadrature weights
  /// corresponding to polynomial order n
  const RealArray& collocation_weights(unsigned short order);

  /// return alphaPoly
  const Real& alpha_polynomial() const;
  /// set alphaPoly using the conversion alphaPoly = alpha_stat-1.
  void alpha_stat(const Real& alpha);

private:

  //
  //- Heading: Data
  //

  /// the alpha parameter for the generalized Laguerre polynomial as defined
  /// by Abramowitz and Stegun (differs from statistical PDF notation)
  Real alphaPoly;
};


inline GenLaguerreOrthogPolynomial::GenLaguerreOrthogPolynomial(): alphaPoly(0.)
{ collocRule = GEN_GAUSS_LAGUERRE; }


// TO DO
inline GenLaguerreOrthogPolynomial::
GenLaguerreOrthogPolynomial(const Real& alpha_stat): alphaPoly(alpha_stat-1.)
{ collocRule = GEN_GAUSS_LAGUERRE; }


inline GenLaguerreOrthogPolynomial::~GenLaguerreOrthogPolynomial()
{ }


inline const Real& GenLaguerreOrthogPolynomial::alpha_polynomial() const
{ return alphaPoly; }


inline void GenLaguerreOrthogPolynomial::alpha_stat(const Real& alpha)
{ alphaPoly = alpha - 1.; reset_gauss(); }

} // namespace Pecos

#endif

/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        ChebyshevOrthogPolynomial
//- Description:  Class for Chebyshev Orthogonal Polynomial
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef CHEBYSHEV_ORTHOG_POLYNOMIAL_HPP
#define CHEBYSHEV_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"


namespace Pecos {

/// Derived orthogonal polynomial class for Chebyshev polynomials

/** The ChebyshevOrthogPolynomial class evaluates a univariate
    Chebyshev polynomial of the first kind (T_n(x)) of a particular
    order.  These polynomials are orthogonal with respect to the
    weight function 1/sqrt(1-x^2) when integrated over the support
    range of [-1,+1].  It enables (mixed) multidimensional orthogonal
    polynomial basis functions within OrthogPolyApproximation.
    Chebyshev polynomials are a special case of the more general
    Jacobi polynomials (implemented separately). */

class ChebyshevOrthogPolynomial: public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  ChebyshevOrthogPolynomial(short colloc_rule); ///< extended constructor
  ChebyshevOrthogPolynomial();                  ///< default constructor
  ~ChebyshevOrthogPolynomial();                 ///< destructor

  //
  //- Heading: Virtual function redefinitions
  //

  /// return the Chebyshev quadrature points corresponding to
  /// polynomial order n
  const RealArray& collocation_points(unsigned short order);
  /// return the Chebyshev quadrature weights corresponding to
  /// polynomial order n
  const RealArray& collocation_weights(unsigned short order);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// retrieve the Chebyshev polynomial value for a given parameter x 
  const Real& get_value(const Real& x, unsigned short order);
  /// retrieve the Chebyshev polynomial gradient for a given parameter x 
  const Real& get_gradient(const Real& x, unsigned short order);

  /// return the inner product <T_n,T_n> = ||T_n||^2
  const Real& norm_squared(unsigned short order);

private:

  //
  //- Heading: Data
  //
};


// collocRule may be CLENSHAW_CURTIS (default) or FEJER2
inline ChebyshevOrthogPolynomial::ChebyshevOrthogPolynomial(short colloc_rule)
{ collocRule = colloc_rule;     wtFactor = 0.5; }


inline ChebyshevOrthogPolynomial::ChebyshevOrthogPolynomial()
{ collocRule = CLENSHAW_CURTIS; wtFactor = 0.5; }


inline ChebyshevOrthogPolynomial::~ChebyshevOrthogPolynomial()
{ }

} // namespace Pecos

#endif

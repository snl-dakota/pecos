/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        KrawtchoukOrthogPolynomial
//- Description:  Class for Krawtchouk Orthogonal Polynomial
//-               
//- Owner:        Russell Hooper, Sandia National Laboratories

#ifndef KRAWTCHOUK_ORTHOG_POLYNOMIAL_HPP
#define KRAWTCHOUK_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"


namespace Pecos {

/// Derived orthogonal polynomial class for Krawtchouk polynomials

/** The KrawtchoukOrthogPolynomial class evaluates a univariate Krawtchouk
    polynomial K^(p,N)_n of a particular order.  These polynomials
    are orthogonal with respect to the weight function (N choose
    k)*p^k*(1-p)^(n-k) hen summed over the discrete points, N.
    This corresponds to the binomial probability mass function.
    See appendix in Xiu & Karniadakis, Siam J. Sci. Comp., v24, n2,
    pp. 619-644, 2002 for more details.  */

class KrawtchoukOrthogPolynomial: public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  KrawtchoukOrthogPolynomial();
  /// destructor
  ~KrawtchoukOrthogPolynomial();

  //
  //- Heading: Virtual function redefinitions
  //

  //
  //- Heading: Noninherited memeber functions
  //
  void set_N(short Nval);
  void set_p(Real pval);

  short get_N();
  Real get_p();

protected:

  //
  //- Heading: Virtual function redefinitions
  //
  Real type1_value(Real x, unsigned short order);


private:

  //
  //- Heading: Data
  //

  /// the number of discrete points on which to base the polynomial 
  short N;
  /// the probability of a "success" for each experiment
  Real p;
};


inline KrawtchoukOrthogPolynomial::KrawtchoukOrthogPolynomial() :
  N(0), p(-1.0)
{ }

inline KrawtchoukOrthogPolynomial::~KrawtchoukOrthogPolynomial()
{ }

inline void KrawtchoukOrthogPolynomial::set_N(short Nval)
{ N = Nval; }

inline void KrawtchoukOrthogPolynomial::set_p(Real pval)
{ p = pval; }

inline short KrawtchoukOrthogPolynomial::get_N()
{ return N; }

inline Real KrawtchoukOrthogPolynomial::get_p()
{ return p; }

} // namespace Pecos

#endif

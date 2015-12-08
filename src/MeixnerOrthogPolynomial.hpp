/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        MeixnerOrthogPolynomial
//- Description:  Class for Meixner Orthogonal Polynomial
//-               
//- Owner:        Russell Hooper, Sandia National Laboratories

#ifndef MEIXNER_ORTHOG_POLYNOMIAL_HPP
#define MEIXNER_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"


namespace Pecos {

/// Derived orthogonal polynomial class for Meixner polynomials

/** The MeixnerOrthogPolynomial class evaluates a univariate Meixner
    polynomial M^(c,Beta)_n of a particular order.  These polynomials
    are orthogonal with respect to the weight function 
    
    EDIT HERE:

    (N choose
    k)*p^k*(1-p)^(n-k) hen summed over the discrete points, N.
    This corresponds to the binomial probability mass function.
    See appendix in Xiu & Karniadakis, Siam J. Sci. Comp., v24, n2,
    pp. 619-644, 2002 for more details.  */

class MeixnerOrthogPolynomial: public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  MeixnerOrthogPolynomial();
  /// destructor
  ~MeixnerOrthogPolynomial();

  //
  //- Heading: Virtual function redefinitions
  //

  //
  //- Heading: Noninherited memeber functions
  //
  void set_c(Real cval);
  void set_beta(Real betaval);

  Real get_c();
  Real get_beta();

protected:

  //
  //- Heading: Virtual function redefinitions
  //
  Real type1_value(Real x, unsigned short order);


private:

  //
  //- Heading: Data
  //

  /// the c-parameter associated with a negative binomial distribution
  Real c;
  /// the beta-parameter associated with a negative binomial distribution
  Real beta;
};


inline MeixnerOrthogPolynomial::MeixnerOrthogPolynomial() :
  c(-1.0), beta(-1.0)
{ }

inline MeixnerOrthogPolynomial::~MeixnerOrthogPolynomial()
{ }

inline void MeixnerOrthogPolynomial::set_c(Real cval)
{ c = cval; }

inline void MeixnerOrthogPolynomial::set_beta(Real betaval)
{ beta = betaval; }

inline Real MeixnerOrthogPolynomial::get_c()
{ return c; }

inline Real MeixnerOrthogPolynomial::get_beta()
{ return beta; }

} // namespace Pecos

#endif

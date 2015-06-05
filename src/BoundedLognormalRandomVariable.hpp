/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 BoundedLognormalRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef BOUNDED_LOGNORMAL_RANDOM_VARIABLE_HPP
#define BOUNDED_LOGNORMAL_RANDOM_VARIABLE_HPP

#include "LognormalRandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for bounded lognormal random variables.

/** Manages lower and upper bounds. */

class BoundedLognormalRandomVariable: public LognormalRandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  BoundedLognormalRandomVariable();
  /// alternate constructor
  BoundedLognormalRandomVariable(Real lambda, Real zeta, Real lwr, Real upr);
  /// destructor
  ~BoundedLognormalRandomVariable();

  //
  //- Heading: Member functions
  //

  Real cdf(Real x) const;
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  //Real pdf_hessian(Real x) const;

  //Real log_pdf(Real x) const;

  void update(Real lambda, Real zeta, Real lwr, Real upr);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real lambda, Real zeta, Real lwr, Real upr);
  static Real cdf(Real x, Real lambda, Real zeta, Real lwr, Real upr);

protected:

  //
  //- Heading: Data
  //

  /// lower bound of bounded_lognormal random variable
  Real lowerBnd;
  /// upper bound of bounded_lognormal random variable
  Real upperBnd;
};


inline BoundedLognormalRandomVariable::BoundedLognormalRandomVariable():
  LognormalRandomVariable(), lowerBnd(0.),
  upperBnd(std::numeric_limits<Real>::infinity())
{ }


inline BoundedLognormalRandomVariable::
BoundedLognormalRandomVariable(Real lambda, Real zeta, Real lwr, Real upr):
  LognormalRandomVariable(lambda, zeta), lowerBnd(lwr), upperBnd(upr)
{ }


inline BoundedLognormalRandomVariable::~BoundedLognormalRandomVariable()
{ }


inline Real BoundedLognormalRandomVariable::
pdf(Real x, Real lambda, Real zeta, Real lwr, Real upr)
{
  if (x < lwr || x > upr) return 0.;
  else {
    Real Phi_lms = (lwr > 0.) ?
      NormalRandomVariable::std_cdf((std::log(lwr)-lambda)/zeta) : 0.;
    Real Phi_ums = (upr < std::numeric_limits<Real>::infinity()) ?
      NormalRandomVariable::std_cdf((std::log(upr)-lambda)/zeta) : 1.;
    return NormalRandomVariable::std_pdf((std::log(x)-lambda)/zeta) /
      (Phi_ums-Phi_lms)/x/zeta;
  }
}


inline Real BoundedLognormalRandomVariable::
cdf(Real x, Real lambda, Real zeta, Real lwr, Real upr)
{
  if      (x < lwr) return 0.;
  else if (x > upr) return 1.;
  else {
    Real Phi_lms = (lwr > 0.) ?
      NormalRandomVariable::std_cdf((std::log(lwr)-lambda)/zeta) : 0.;
    Real Phi_ums = (upr < std::numeric_limits<Real>::infinity()) ?
      NormalRandomVariable::std_cdf((std::log(upr)-lambda)/zeta) : 1.;
    return (NormalRandomVariable::std_cdf((std::log(x)-lambda)/zeta) - Phi_lms)
      / (Phi_ums - Phi_lms);
  }
}


inline Real BoundedLognormalRandomVariable::cdf(Real x) const
{ return cdf(x, lnLambda, lnZeta, lowerBnd, upperBnd); }


inline Real BoundedLognormalRandomVariable::ccdf(Real x) const
{
  if      (x < lowerBnd) return 1.;
  else if (x > upperBnd) return 0.;
  else {
    Real Phi_lms = (lowerBnd > 0.) ?
      NormalRandomVariable::std_cdf((std::log(lowerBnd)-lnLambda)/lnZeta) : 0.;
    Real Phi_ums = (upperBnd < std::numeric_limits<Real>::infinity()) ?
      NormalRandomVariable::std_cdf((std::log(upperBnd)-lnLambda)/lnZeta) : 1.;
    return (Phi_ums - NormalRandomVariable::
	    std_cdf((std::log(x)-lnLambda)/lnZeta)) / (Phi_ums - Phi_lms);
  }
}


inline Real BoundedLognormalRandomVariable::inverse_cdf(Real p_cdf) const
{
  if      (p_cdf <= 0.) return lowerBnd;
  else if (p_cdf >= 1.) return upperBnd;
  else {
    // p = (Phi((log(x)-lambda)/zeta) - Phi_lms)/(Phi_ums - Phi_lms)
    // log(x) = Phi_inverse[p * (Phi_ums - Phi_lms) + Phi_lms] * zeta + lambda
    Real Phi_lms = (lowerBnd > 0.) ?
      NormalRandomVariable::std_cdf((log(lowerBnd)-lnLambda)/lnZeta) : 0.;
    Real Phi_ums = (upperBnd < std::numeric_limits<Real>::infinity()) ?
      NormalRandomVariable::std_cdf((log(upperBnd)-lnLambda)/lnZeta) : 1.;
    return std::exp(lnLambda + lnZeta * NormalRandomVariable::
		    inverse_std_cdf(p_cdf * (Phi_ums - Phi_lms) + Phi_lms));
  }
}


inline Real BoundedLognormalRandomVariable::inverse_ccdf(Real p_ccdf) const
{
  if      (p_ccdf >= 1.) return lowerBnd;
  else if (p_ccdf <= 0.) return upperBnd;
  else {
    // p = (Phi_ums - Phi((log(x)-lambda)/zeta))/(Phi_ums - Phi_lms)
    // log(x) = lambda + zeta * Phi_inverse[Phi_ums - p * (Phi_ums - Phi_lms)]
    Real Phi_lms = (lowerBnd > 0.) ?
      NormalRandomVariable::std_cdf((log(lowerBnd)-lnLambda)/lnZeta) : 0.;
    Real Phi_ums = (upperBnd < std::numeric_limits<Real>::infinity()) ?
      NormalRandomVariable::std_cdf((log(upperBnd)-lnLambda)/lnZeta) : 1.;
    return std::exp(lnLambda + lnZeta * NormalRandomVariable::
		    inverse_std_cdf(Phi_ums - p_ccdf * (Phi_ums - Phi_lms)));
  }
}


inline Real BoundedLognormalRandomVariable::pdf(Real x) const
{ return pdf(x, lnLambda, lnZeta, lowerBnd, upperBnd); }


inline Real BoundedLognormalRandomVariable::pdf_gradient(Real x) const
{ return -pdf(x) * (1. + (log(x)-lnLambda)/(lnZeta*lnZeta)) / x; }


//inline Real BoundedLognormalRandomVariable::pdf_hessian(Real x) const
//{
//  return pdf(x) * ...; // TO DO
//}


inline void BoundedLognormalRandomVariable::
update(Real lambda, Real zeta, Real lwr, Real upr)
{ lnLambda = lambda; lnZeta = zeta; lowerBnd = lwr; upperBnd = upr; }

} // namespace Pecos

#endif

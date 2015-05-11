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
  BoundedLognormalRandomVariable(Real mean, Real stdev, Real lwr, Real upr);
  /// destructor
  ~BoundedLognormalRandomVariable();

  //
  //- Heading: Member functions
  //

  Real cdf(Real x) const;
  Real cdf_inverse(Real p) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  //Real pdf_hessian(Real x) const;

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
  LognormalRandomVariable(), lowerBnd(-1.), upperBnd(1.) // TO DO: +/- inf?
{ }


inline BoundedLognormalRandomVariable::
BoundedLognormalRandomVariable(Real mean, Real stdev, Real lwr, Real upr):
  LognormalRandomVariable(mean, stdev), lowerBnd(lwr), upperBnd(upr)
{ }


inline BoundedLognormalRandomVariable::~BoundedLognormalRandomVariable()
{ }


inline Real BoundedLognormalRandomVariable::cdf(Real x) const
{ return bounded_lognormal_cdf(x, gaussMean, gaussStdDev, lowerBnd, upperBnd); }


inline Real BoundedLognormalRandomVariable::cdf_inverse(Real p) const
{
  // p = (Phi((log(x)-lambda)/zeta) - Phi_lms)/(Phi_ums - Phi_lms)
  // log(x) = Phi_inverse[p * (Phi_ums - Phi_lms) + Phi_lms] * zeta + lambda
  Real Phi_lms = (lowerBnd > 0.) ? Phi((log(lowerBnd)-lambda)/zeta) : 0.;
  Real Phi_ums = (upperBnd < std::numeric_limits<Real>::infinity()) ?
    Phi((log(upperBnd)-lambda)/zeta) : 1.;
  return
    std::exp(Phi_inverse(p * (Phi_ums - Phi_lms) + Phi_lms) * zeta + lambda);
}


inline Real BoundedLognormalRandomVariable::pdf(Real x) const
{ return bounded_lognormal_pdf(x, gaussMean, gaussStdDev, lowerBnd, upperBnd); }


inline Real BoundedLognormalRandomVariable::pdf_gradient(Real x) const
{
  Real xms = (log(x)-lambda)/zeta,
    pdf = bounded_lognormal_pdf(x, gaussMean, gaussStdDev, lowerBnd, upperBnd);
  return -pdf * (1. + xms / zeta) / x;
}


//inline Real BoundedLognormalRandomVariable::pdf_hessian(Real x) const
//{
//  return bounded_lognormal_pdf(x, gaussMean, gaussStdDev, lowerBnd, upperBnd);
// * ...; // TO DO
//}

} // namespace Pecos

#endif

/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 BoundedNormalRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef BOUNDED_NORMAL_RANDOM_VARIABLE_HPP
#define BOUNDED_NORMAL_RANDOM_VARIABLE_HPP

#include "NormalRandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for bounded normal random variables.

/** Manages lower and upper bounds. */

class BoundedNormalRandomVariable: public NormalRandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  BoundedNormalRandomVariable();
  /// alternate constructor
  BoundedNormalRandomVariable(Real mean, Real stdev, Real lwr, Real upr);
  /// destructor
  ~BoundedNormalRandomVariable();

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

  /// lower bound of bounded_normal random variable
  Real lowerBnd;
  /// upper bound of bounded_normal random variable
  Real upperBnd;
};


inline BoundedNormalRandomVariable::BoundedNormalRandomVariable():
  NormalRandomVariable(), lowerBnd(-1.), upperBnd(1.) // TO DO: +/- inf?
{ }


inline BoundedNormalRandomVariable::
BoundedNormalRandomVariable(Real mean, Real stdev, Real lwr, Real upr):
  NormalRandomVariable(mean, stdev), lowerBnd(lwr), upperBnd(upr)
{ }


inline BoundedNormalRandomVariable::~BoundedNormalRandomVariable()
{ }


inline Real BoundedNormalRandomVariable::cdf(Real x) const
{ return bounded_normal_cdf(x, gaussMean, gaussStdDev, lowerBnd, upperBnd); }


inline Real BoundedNormalRandomVariable::cdf_inverse(Real p) const
{
  // p = (Phi((x-mean)/std_dev) - Phi_lms)/(Phi_ums - Phi_lms)
  Real dbl_inf = std::numeric_limits<Real>::infinity();
  Real Phi_lms
    = (lowerBnd > -dbl_inf) ? Phi((lowerBnd-gaussMean)/gaussStdDev) : 0.;
  Real Phi_ums
    = (upperBnd <  dbl_inf) ? Phi((upperBnd-gaussMean)/gaussStdDev) : 1.;
  return
    Phi_inverse(p * (Phi_ums - Phi_lms) + Phi_lms) * gaussStdDev + gaussMean;
}


inline Real BoundedNormalRandomVariable::pdf(Real x) const
{ return bounded_normal_pdf(x, gaussMean, gaussStdDev, lowerBnd, upperBnd); }


inline Real BoundedNormalRandomVariable::pdf_gradient(Real x) const
{
  return pdf(x) * (gaussMean - x) / (gaussStdDev * gaussStdDev);
}


//inline Real BoundedNormalRandomVariable::pdf_hessian(Real x) const
//{
//  return pdf(x); // * ...; // TO DO
//}

} // namespace Pecos

#endif

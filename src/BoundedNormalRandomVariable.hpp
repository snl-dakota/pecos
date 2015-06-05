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
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  //Real pdf_hessian(Real x) const;

  //Real log_pdf(Real x) const;

  void update(Real mean, Real stdev, Real lwr, Real upr);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real mean, Real std_dev, Real lwr, Real upr);
  static Real cdf(Real x, Real mean, Real std_dev, Real lwr, Real upr);

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
  NormalRandomVariable(), lowerBnd(-std::numeric_limits<Real>::infinity()),
  upperBnd(std::numeric_limits<Real>::infinity())
{ }


inline BoundedNormalRandomVariable::
BoundedNormalRandomVariable(Real mean, Real stdev, Real lwr, Real upr):
  NormalRandomVariable(mean, stdev), lowerBnd(lwr), upperBnd(upr)
{ }


inline BoundedNormalRandomVariable::~BoundedNormalRandomVariable()
{ }


inline Real BoundedNormalRandomVariable::
pdf(Real x, Real mean, Real std_dev, Real lwr, Real upr)
{
  if (x < lwr || x > upr) return 0.;
  else {
    Real dbl_inf = std::numeric_limits<Real>::infinity();
    Real Phi_lms = (lwr > -dbl_inf) ? std_cdf((lwr-mean)/std_dev) : 0.;
    Real Phi_ums = (upr <  dbl_inf) ? std_cdf((upr-mean)/std_dev) : 1.;
    return std_pdf((x-mean)/std_dev)/(Phi_ums - Phi_lms)/std_dev;
  }
}


inline Real BoundedNormalRandomVariable::
cdf(Real x, Real mean, Real std_dev, Real lwr, Real upr)
{
  if      (x < lwr) return 0.;
  else if (x > upr) return 1.;
  else {
    Real dbl_inf = std::numeric_limits<Real>::infinity();
    Real Phi_lms = (lwr > -dbl_inf) ? std_cdf((lwr-mean)/std_dev) : 0.;
    Real Phi_ums = (upr <  dbl_inf) ? std_cdf((upr-mean)/std_dev) : 1.;
    return (std_cdf((x-mean)/std_dev) - Phi_lms) / (Phi_ums - Phi_lms);
  }
}


inline Real BoundedNormalRandomVariable::cdf(Real x) const
{ return cdf(x, gaussMean, gaussStdDev, lowerBnd, upperBnd); }


inline Real BoundedNormalRandomVariable::ccdf(Real x) const
{
  if      (x < lowerBnd) return 1.;
  else if (x > upperBnd) return 0.;
  else {
    Real dbl_inf = std::numeric_limits<Real>::infinity();
    Real Phi_lms = (lowerBnd > -dbl_inf) ?
      std_cdf((lowerBnd-gaussMean)/gaussStdDev) : 0.;
    Real Phi_ums = (upperBnd <  dbl_inf) ?
      std_cdf((upperBnd-gaussMean)/gaussStdDev) : 1.;
    return (Phi_ums - std_cdf((x-gaussMean)/gaussStdDev)) / (Phi_ums - Phi_lms);
  }
}


inline Real BoundedNormalRandomVariable::inverse_cdf(Real p_cdf) const
{
  if      (p_cdf <= 0.) return lowerBnd;
  else if (p_cdf >= 1.) return upperBnd;
  else {
    // p = (Phi((x-mean)/std_dev) - Phi_lms)/(Phi_ums - Phi_lms)
    Real dbl_inf = std::numeric_limits<Real>::infinity();
    Real Phi_lms = (lowerBnd > -dbl_inf) ?
      std_cdf((lowerBnd-gaussMean)/gaussStdDev) : 0.;
    Real Phi_ums = (upperBnd <  dbl_inf) ?
      std_cdf((upperBnd-gaussMean)/gaussStdDev) : 1.;
    return gaussMean + gaussStdDev *
      inverse_std_cdf(p_cdf * (Phi_ums - Phi_lms) + Phi_lms);
  }
}


inline Real BoundedNormalRandomVariable::inverse_ccdf(Real p_ccdf) const
{
  if      (p_ccdf >= 1.) return lowerBnd;
  else if (p_ccdf <= 0.) return upperBnd;
  else {
    // p = (Phi_ums - Phi((x-mean)/std_dev))/(Phi_ums - Phi_lms)
    Real dbl_inf = std::numeric_limits<Real>::infinity();
    Real Phi_lms = (lowerBnd > -dbl_inf) ?
      std_cdf((lowerBnd-gaussMean)/gaussStdDev) : 0.;
    Real Phi_ums = (upperBnd <  dbl_inf) ?
      std_cdf((upperBnd-gaussMean)/gaussStdDev) : 1.;
    return gaussMean + gaussStdDev *
      inverse_std_cdf(Phi_ums - p_ccdf * (Phi_ums - Phi_lms));
  }
}


inline Real BoundedNormalRandomVariable::pdf(Real x) const
{ return pdf(x, gaussMean, gaussStdDev, lowerBnd, upperBnd); }


inline Real BoundedNormalRandomVariable::pdf_gradient(Real x) const
{ return pdf(x) * (gaussMean - x) / (gaussStdDev * gaussStdDev); }


//inline Real BoundedNormalRandomVariable::pdf_hessian(Real x) const
//{
//  return pdf(x) * ...; // TO DO
//}


inline void BoundedNormalRandomVariable::
update(Real mean, Real stdev, Real lwr, Real upr)
{ gaussMean = mean; gaussStdDev = stdev; lowerBnd = lwr; upperBnd = upr; }

} // namespace Pecos

#endif

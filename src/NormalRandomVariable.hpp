/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 NormalRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef NORMAL_RANDOM_VARIABLE_HPP
#define NORMAL_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for Gaussian random variables.

/** Manages mean and standard deviation. */

class NormalRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  NormalRandomVariable();                      ///< constructor
  NormalRandomVariable(Real mean, Real stdev); ///< alternate constructor
  ~NormalRandomVariable();                     ///< destructor

  //
  //- Heading: Member functions
  //

  Real cdf(Real x) const;
  Real cdf_inverse(Real p) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  Real pdf_hessian(Real x) const;

protected:

  //
  //- Heading: Data
  //

  /// mean of Gaussian random variable
  Real gaussMean;
  /// standard deviation of Gaussian random variable
  Real gaussStdDev;

  // normal distribution instance from Boost math
  //normal_dist* normDist;
};


inline NormalRandomVariable::NormalRandomVariable():
  RandomVariable(BaseConstructor()), gaussMean(0), gaussStdDev(1.)
{ }


inline NormalRandomVariable::NormalRandomVariable(Real mean, Real stdev):
  RandomVariable(BaseConstructor()), gaussMean(mean), gaussStdDev(stdev)
{ }


inline NormalRandomVariable::~NormalRandomVariable()
{ }


inline Real NormalRandomVariable::cdf(Real x) const
{
  normal_dist norm(gaussMean, gaussStdDev);
  return bmth::cdf(norm, x);
}


inline Real NormalRandomVariable::cdf_inverse(Real p) const
{
  normal_dist norm(gaussMean, gaussStdDev);
  return bmth::quantile(norm, p); 
}


inline Real NormalRandomVariable::pdf(Real x) const
{
  normal_dist norm(gaussMean, gaussStdDev);
  return bmth::pdf(norm, x);
}


inline Real NormalRandomVariable::pdf_gradient(Real x) const
{
  return pdf(x) * (gaussMean - x) / (gaussStdDev * gaussStdDev);
}


inline Real NormalRandomVariable::pdf_hessian(Real x) const
{
  Real var = gaussStdDev * gaussStdDev, mu_minus_x = gaussMean - x;
  return pdf(x) * ( mu_minus_x * mu_minus_x / var - 1. ) / var;
}

} // namespace Pecos

#endif

/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 LoguniformRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef LOGUNIFORM_RANDOM_VARIABLE_HPP
#define LOGUNIFORM_RANDOM_VARIABLE_HPP

#include "UniformRandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for loguniform random variables.

/** Manages lower and upper bounds.  See SAND98-0210 LHS manual, pp. 43-44. */

class LoguniformRandomVariable: public UniformRandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  LoguniformRandomVariable();                   ///< default constructor
  LoguniformRandomVariable(Real lwr, Real upr); ///< alternate constructor
  ~LoguniformRandomVariable();                  ///< destructor

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

  //void update(Real lwr, Real upr); // inherits from UniformRV

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real lwr, Real upr);
  static Real cdf(Real x, Real lwr, Real upr);

  static void moments_from_params(Real lwr, Real upr, Real& mean,
				  Real& std_dev);

protected:

  //
  //- Heading: Data
  //
};


inline LoguniformRandomVariable::LoguniformRandomVariable():
  UniformRandomVariable()
{ }


inline LoguniformRandomVariable::LoguniformRandomVariable(Real lwr, Real upr):
  UniformRandomVariable(lwr, upr)
{ }


inline LoguniformRandomVariable::~LoguniformRandomVariable()
{ }


inline Real LoguniformRandomVariable::cdf(Real x) const
{
  return (std::log(x)        - std::log(lowerBnd)) /
         (std::log(upperBnd) - std::log(lowerBnd));
}


inline Real LoguniformRandomVariable::ccdf(Real x) const
{
  return (std::log(upperBnd) - std::log(x)) /
         (std::log(upperBnd) - std::log(lowerBnd));
}


inline Real LoguniformRandomVariable::inverse_cdf(Real p_cdf) const
{
  // p = (ln x - ln L)/(ln U - ln L)
  return lowerBnd * std::exp(p_cdf * (std::log(upperBnd) - std::log(lowerBnd)));
}


inline Real LoguniformRandomVariable::inverse_ccdf(Real p_ccdf) const
{
  // p = (ln U - ln x)/(ln U - ln L)
  return upperBnd /
    std::exp(p_ccdf * (std::log(upperBnd) - std::log(lowerBnd)));
}


inline Real LoguniformRandomVariable::pdf(Real x) const
{ return 1./(std::log(upperBnd) - std::log(lowerBnd))/x; }


inline Real LoguniformRandomVariable::pdf_gradient(Real x) const
{ return -pdf(x, lowerBnd, upperBnd) / x; }


//inline Real LoguniformRandomVariable::pdf_hessian(Real x) const
//{ return pdf(x) * ; }

// static functions:

inline Real LoguniformRandomVariable::pdf(Real x, Real lwr, Real upr)
{ return 1./(std::log(upr) - std::log(lwr))/x; }


inline Real LoguniformRandomVariable::cdf(Real x, Real lwr, Real upr)
{ return (std::log(x) - std::log(lwr))/(std::log(upr) - std::log(lwr)); }


inline void LoguniformRandomVariable::
moments_from_params(Real lwr, Real upr, Real& mean, Real& std_dev)
{
  Real range = upr - lwr, log_range = std::log(upr) - std::log(lwr);
  mean       = range/log_range;
  std_dev    = std::sqrt(range*(log_range*(upr+lwr)-2.*range)/2.)/log_range;
}

} // namespace Pecos

#endif

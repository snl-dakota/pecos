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

/** Manages lower and upper bounds. */

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
  Real cdf_inverse(Real p) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  //Real pdf_hessian(Real x) const;

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
{ return loguniform_cdf(x, lowerBnd, upperBnd); }


inline Real LoguniformRandomVariable::cdf_inverse(Real p) const
{
  // p = (ln x - ln L)/(ln U - ln L)
  Real log_lwr = std::log(lowerBnd);
  return std::exp(p * (std::log(upperBnd) - log_lwr) + log_lwr);
}


inline Real LoguniformRandomVariable::pdf(Real x) const
{ return loguniform_pdf(x, lowerBnd, upperBnd); }


inline Real LoguniformRandomVariable::pdf_gradient(Real x) const
{ return -loguniform_pdf(x, lowerBnd, upperBnd) / x; }


//inline Real LoguniformRandomVariable::pdf_hessian(Real x) const
//{ return pdf(x); }// * ; }

} // namespace Pecos

#endif

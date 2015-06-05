/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 ExponentialRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef EXPONENTIAL_RANDOM_VARIABLE_HPP
#define EXPONENTIAL_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for exponential random variables.

/** Manages beta parameter.  Pecos employs the 1/beta exp(-x/beta)
    definition, which differs from the lambda exp(-lambda x) LHS
    definition (lambda_LHS = 1/beta_PECOS). */

class ExponentialRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  ExponentialRandomVariable();
  /// alternate constructor
  ExponentialRandomVariable(Real beta);
  /// destructor
  ~ExponentialRandomVariable();

  //
  //- Heading: Member functions
  //

  Real cdf(Real x) const;
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  Real pdf_hessian(Real x) const;

  Real inverse_log_ccdf(Real log_p_ccdf) const;
  Real log_pdf(Real x) const;

  Real to_std(Real x) const;
  Real from_std(Real z) const;

  void update(Real beta);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real std_pdf(Real x);
  static Real std_cdf(Real x);

  static Real pdf(Real x, Real beta);
  static Real cdf(Real x, Real beta);

  static void moments_from_params(Real beta, Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Data
  //

  /// beta parameter of exponential random variable
  Real betaStat;
};


inline ExponentialRandomVariable::ExponentialRandomVariable():
  RandomVariable(BaseConstructor()), betaStat(1.)
{ }


inline ExponentialRandomVariable::ExponentialRandomVariable(Real beta):
  RandomVariable(BaseConstructor()), betaStat(beta)
{ }


inline ExponentialRandomVariable::~ExponentialRandomVariable()
{ }


inline Real ExponentialRandomVariable::cdf(Real x) const
{ return -bmth::expm1(-x/betaStat); }


inline Real ExponentialRandomVariable::ccdf(Real x) const
{ return std::exp(-x/betaStat); }


inline Real ExponentialRandomVariable::inverse_cdf(Real p_cdf) const
{
  // p_cdf = 1 - exp(-x/beta)  -->  -x/beta = log(1-p_cdf)
  return -betaStat * bmth::log1p(-p_cdf);
}


inline Real ExponentialRandomVariable::inverse_ccdf(Real p_ccdf) const
{ return -betaStat * std::log(p_ccdf); }


inline Real ExponentialRandomVariable::inverse_log_ccdf(Real log_p_ccdf) const
{ return -betaStat * log_p_ccdf; }


inline Real ExponentialRandomVariable::pdf(Real x) const
{ return std::exp(-x/betaStat)/betaStat; }


inline Real ExponentialRandomVariable::pdf_gradient(Real x) const
{ return -pdf(x, betaStat) / betaStat; }


inline Real ExponentialRandomVariable::pdf_hessian(Real x) const
{ return pdf(x, betaStat) / (betaStat * betaStat); }


inline Real ExponentialRandomVariable::log_pdf(Real x) const
{ return -x / betaStat - std::log(betaStat); }


inline Real ExponentialRandomVariable::to_std(Real x) const
{ return x / betaStat; }


inline Real ExponentialRandomVariable::from_std(Real z) const
{ return z * betaStat; }


inline void ExponentialRandomVariable::update(Real beta)
{ betaStat = beta; }

// static functions:

inline Real ExponentialRandomVariable::std_pdf(Real x)
{ return std::exp(-x); }


inline Real ExponentialRandomVariable::std_cdf(Real x)
{
  // as with log1p(), avoid numerical probs when exp(~0) is ~ 1
  return -bmth::expm1(-x); //1. - std::exp(-x);
}


inline Real ExponentialRandomVariable::pdf(Real x, Real beta)
{ return std::exp(-x/beta)/beta; }


inline Real ExponentialRandomVariable::cdf(Real x, Real beta)
{
  // as with log1p(), avoid numerical probs when exp(~0) is ~ 1
  return -bmth::expm1(-x/beta);
}


inline void ExponentialRandomVariable::
moments_from_params(Real beta, Real& mean, Real& std_dev)
{ mean = std_dev = beta; }

} // namespace Pecos

#endif

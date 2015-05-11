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

/** Manages beta parameter. */

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
  Real cdf_inverse(Real p) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  Real pdf_hessian(Real x) const;

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
{ return exponential_cdf(x, betaStat); }


inline Real ExponentialRandomVariable::cdf_inverse(Real p) const
{
  // p = 1. - std::exp(-x/beta)
  // -x/beta = log(1.-p)
  return -betaStat * bmth::log1p(-p);
}


inline Real ExponentialRandomVariable::pdf(Real x) const
{ return exponential_pdf(x, betaStat); }


inline Real ExponentialRandomVariable::pdf_gradient(Real x) const
{ return -pdf(x) / betaStat; }


inline Real ExponentialRandomVariable::pdf_hessian(Real x) const
{ return pdf(x) / (betaStat * betaStat); }

} // namespace Pecos

#endif

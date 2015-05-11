/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 WeibullRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef WEIBULL_RANDOM_VARIABLE_HPP
#define WEIBULL_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for weibull random variables.

/** Manages alpha and beta parameters. */

class WeibullRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  WeibullRandomVariable();
  /// alternate constructor
  WeibullRandomVariable(Real alpha, Real beta);
  /// destructor
  ~WeibullRandomVariable();

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

  /// alpha parameter of weibull random variable
  Real alphaStat;
  /// beta parameter of weibull random variable
  Real betaStat;
};


inline WeibullRandomVariable::WeibullRandomVariable():
  RandomVariable(BaseConstructor()), alphaStat(0), betaStat(0.)
{ }


inline WeibullRandomVariable::WeibullRandomVariable(Real alpha, Real beta):
  RandomVariable(BaseConstructor()), alphaStat(alpha), betaStat(beta)
{ }


inline WeibullRandomVariable::~WeibullRandomVariable()
{ }


inline Real WeibullRandomVariable::cdf(Real x) const
{ return weibull_cdf(x, alphaStat, betaStat); }


inline Real WeibullRandomVariable::cdf_inverse(Real p) const
{
  // p = 1 - e^(-(x/beta)^alpha)
  return betaStat * std::pow(-bmth::log1p(-p), 1./alphaStat);
}


inline Real WeibullRandomVariable::pdf(Real x) const
{ return weibull_pdf(x, alphaStat, betaStat); }


inline Real WeibullRandomVariable::pdf_gradient(Real x) const
{
  Real num = x / betaStat, num2 = std::exp(-std::pow(num, alphaStat)),
    ab_ratio = alphaStat/betaStat,
    pdf = ab_ratio * num2 * std::pow(num, alphaStat - 1.);
  return ab_ratio * (num2 * (alphaStat-1.) / betaStat * 
		     std::pow(num, alphaStat - 2.) -
		     std::pow(num, alphaStat - 1.) * pdf);
}


//inline Real WeibullRandomVariable::pdf_hessian(Real x) const
//{
//  return weibull_pdf(x, alphaStat, betaStat);
// * ...; // TO DO
//}

} // namespace Pecos

#endif

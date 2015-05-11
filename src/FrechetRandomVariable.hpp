/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 FrechetRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef FRECHET_RANDOM_VARIABLE_HPP
#define FRECHET_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for frechet random variables.

/** Manages alpha and beta parameters. */

class FrechetRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  FrechetRandomVariable();
  /// alternate constructor
  FrechetRandomVariable(Real alpha, Real beta);
  /// destructor
  ~FrechetRandomVariable();

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

  /// alpha parameter of frechet random variable
  Real alphaStat;
  /// beta parameter of frechet random variable
  Real betaStat;
};


inline FrechetRandomVariable::FrechetRandomVariable():
  RandomVariable(BaseConstructor()), alphaStat(0), betaStat(0.)
{ }


inline FrechetRandomVariable::FrechetRandomVariable(Real alpha, Real beta):
  RandomVariable(BaseConstructor()), alphaStat(alpha), betaStat(beta)
{ }


inline FrechetRandomVariable::~FrechetRandomVariable()
{ }


inline Real FrechetRandomVariable::cdf(Real x) const
{ return frechet_cdf(x, alphaStat, betaStat); }


inline Real FrechetRandomVariable::cdf_inverse(Real p) const
{
  // p = std::exp(-std::pow(beta/x, alpha))
  // beta/x = std::pow( -log p, 1/alpha)
  return betaStat * std::pow(-std::log(p), -1./alphaStat);
}


inline Real FrechetRandomVariable::pdf(Real x) const
{ return frechet_pdf(x, alphaStat, betaStat); }


inline Real FrechetRandomVariable::pdf_gradient(Real x) const
{
  Real num = betaStat/x, cdf = std::exp(-std::pow(num, alphaStat)),
    ab_ratio = alphaStat/betaStat,
    pdf = ab_ratio * std::pow(num, alphaStat+1.) * cdf;
  return ab_ratio * (std::pow(num,alphaStat+1.) * pdf - 
		     cdf*(alphaStat+1.)/betaStat * std::pow(num,alphaStat+2.));
}


//inline Real FrechetRandomVariable::pdf_hessian(Real x) const
//{
//  return frechet_pdf(x, alphaStat, betaStat); // * ...; // TO DO
//}

} // namespace Pecos

#endif

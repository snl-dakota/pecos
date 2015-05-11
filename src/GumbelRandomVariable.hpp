/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 GumbelRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef GUMBEL_RANDOM_VARIABLE_HPP
#define GUMBEL_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for gumbel random variables.

/** Manages alpha and beta parameters. */

class GumbelRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  GumbelRandomVariable();
  /// alternate constructor
  GumbelRandomVariable(Real alpha, Real beta);
  /// destructor
  ~GumbelRandomVariable();

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

  /// alpha parameter of gumbel random variable
  Real alphaStat;
  /// beta parameter of gumbel random variable
  Real betaStat;
};


inline GumbelRandomVariable::GumbelRandomVariable():
  RandomVariable(BaseConstructor()), alphaStat(0), betaStat(0.)
{ }


inline GumbelRandomVariable::GumbelRandomVariable(Real alpha, Real beta):
  RandomVariable(BaseConstructor()), alphaStat(alpha), betaStat(beta)
{ }


inline GumbelRandomVariable::~GumbelRandomVariable()
{ }


inline Real GumbelRandomVariable::cdf(Real x) const
{ return gumbel_cdf(x, alphaStat, betaStat); }


inline Real GumbelRandomVariable::cdf_inverse(Real p) const
{
  // p = std::exp(-std::exp(-alpha * (x - beta)))
  return betaStat - std::log(-std::log(p)) / alphaStat;
}


inline Real GumbelRandomVariable::pdf(Real x) const
{ return gumbel_pdf(x, alphaStat, betaStat); }


inline Real GumbelRandomVariable::pdf_gradient(Real x) const
{
  Real num = std::exp(-alphaStat*(x-betaStat)),
    cdf = std::exp(-num), pdf = alphaStat * num * cdf;
  return alphaStat * num * (pdf - alphaStat * cdf);
}


//inline Real GumbelRandomVariable::pdf_hessian(Real x) const
//{
//  return gumbel_pdf(x, alphaStat, betaStat);
// * ...; // TO DO
//}

} // namespace Pecos

#endif

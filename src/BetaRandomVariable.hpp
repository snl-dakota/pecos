/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 BetaRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef BETA_RANDOM_VARIABLE_HPP
#define BETA_RANDOM_VARIABLE_HPP

#include "UniformRandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for beta random variables.

/** Manages alpha and beta and inherits bounds. */

class BetaRandomVariable: public UniformRandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  BetaRandomVariable();
  /// alternate constructor
  BetaRandomVariable(Real alpha, Real beta, Real lwr, Real upr);
  /// destructor
  ~BetaRandomVariable();

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

  /// alpha parameter of beta random variable (statistical PDF
  /// convention; differs from Jacobi polynomial convention)
  Real alphaStat;
  /// beta parameter of beta random variable (statistical PDF
  /// convention; differs from Jacobi polynomial convention)
  Real betaStat;
};


inline BetaRandomVariable::BetaRandomVariable():
  UniformRandomVariable(), alphaStat(0), betaStat(0.)
{ }


inline BetaRandomVariable::
BetaRandomVariable(Real alpha, Real beta, Real lwr, Real upr):
  UniformRandomVariable(lwr, upr), alphaStat(alpha), betaStat(beta)
{ }


inline BetaRandomVariable::~BetaRandomVariable()
{ }


inline Real BetaRandomVariable::cdf(Real x) const
{ return beta_cdf(x, alphaStat, betaStat, lowerBnd, upperBnd); }


inline Real BetaRandomVariable::cdf_inverse(Real p) const
{
  Real std_x = std_beta_cdf_inverse(p, alphaStat, betaStat); // [-1,1]
  // scale from [-1,1] to [L,U]:
  return lowerBnd + (upperBnd - lowerBnd) * (std_x + 1.) / 2.;
}


inline Real BetaRandomVariable::pdf(Real x) const
{ return beta_pdf(x, alphaStat, betaStat, lowerBnd, upperBnd); }


inline Real BetaRandomVariable::pdf_gradient(Real x) const
{
  return beta_pdf(x, alphaStat, betaStat, lowerBnd, upperBnd)
    * ( (alphaStat-1.) / (x-lowerBnd) - (betaStat-1.) / (upperBnd-x) );
}


//inline Real BetaRandomVariable::pdf_hessian(Real x) const
//{
//  return beta_pdf(x, alphaStat, betaStat, lowerBnd, upperBnd);// * ...; TO DO
//}

} // namespace Pecos

#endif

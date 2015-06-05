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

/** Manages alpha and beta and inherits bounds.  See Haldar and
    Mahadevan, p. 72. */

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
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  //Real pdf_hessian(Real x) const;

  // inherited from UniformRandomVariable
  //Real to_std(Real x) const;
  //Real from_std(Real z) const;

  void update(Real alpha, Real beta, Real lwr, Real upr);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real std_pdf(Real x, Real alpha, Real beta);
  static Real std_cdf(Real x, Real alpha, Real beta);
  static Real inverse_std_cdf(Real cdf, Real alpha, Real beta);

  static Real pdf(Real x, Real alpha, Real beta, Real lwr, Real upr);
  static Real cdf(Real x, Real alpha, Real beta, Real lwr, Real upr);

  static void moments_from_params(Real lwr, Real upr, Real alpha, Real beta,
				  Real& mean, Real& std_dev);

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

  /// pointer to the Boost beta_distribution instance
  beta_dist* betaDist;
};


inline BetaRandomVariable::BetaRandomVariable():
  UniformRandomVariable(), alphaStat(1.), betaStat(1.), betaDist(NULL)
{ }


inline BetaRandomVariable::
BetaRandomVariable(Real alpha, Real beta, Real lwr, Real upr):
  UniformRandomVariable(lwr, upr), alphaStat(alpha), betaStat(beta),
  betaDist(new beta_dist(alphaStat, betaStat))
{ }


inline BetaRandomVariable::~BetaRandomVariable()
{ if (betaDist) delete betaDist; }


inline Real BetaRandomVariable::cdf(Real x) const
{
  //return beta_cdf(x, alphaStat, betaStat, lowerBnd, upperBnd);

  Real std_x = (x - lowerBnd)/(upperBnd - lowerBnd);// scale from [L,U] to [0,1]
  return bmth::cdf(*betaDist, std_x);
}


inline Real BetaRandomVariable::ccdf(Real x) const
{
  //return beta_cdf(x, alphaStat, betaStat, lowerBnd, upperBnd);

  Real std_x = (x - lowerBnd)/(upperBnd - lowerBnd);// scale from [L,U] to [0,1]
  return bmth::cdf(complement(*betaDist, std_x));
}


inline Real BetaRandomVariable::inverse_cdf(Real p_cdf) const
{
  //Real std_x = inverse_std_cdf(p, alphaStat, betaStat); // [0,1]

  Real std_x = bmth::quantile(*betaDist, p_cdf);
  return lowerBnd + (upperBnd - lowerBnd) * std_x;  // scale from [0,1] to [L,U]
}


inline Real BetaRandomVariable::inverse_ccdf(Real p_ccdf) const
{
  //Real std_x = inverse_std_cdf(p, alphaStat, betaStat); // [0,1]

  Real std_x = bmth::quantile(complement(*betaDist, p_ccdf));
  return lowerBnd + (upperBnd - lowerBnd) * std_x;  // scale from [0,1] to [L,U]
}


inline Real BetaRandomVariable::pdf(Real x) const
{
  //return beta_pdf(x, alphaStat, betaStat, lowerBnd, upperBnd);

  Real scale = upperBnd - lowerBnd, scaled_x = (x - lowerBnd)/scale;
  return bmth::pdf(*betaDist, scaled_x) / scale;
}


inline Real BetaRandomVariable::pdf_gradient(Real x) const
{
  return pdf(x) *
    ( (alphaStat-1.) / (x-lowerBnd) - (betaStat-1.) / (upperBnd-x) );
}


//inline Real BetaRandomVariable::pdf_hessian(Real x) const
//{
//  return pdf(x) * ...; // TO DO
//}


inline void BetaRandomVariable::
update(Real alpha, Real beta, Real lwr, Real upr)
{
  lowerBnd = lwr; upperBnd = upr; // don't affect betaDist
  if (!betaDist || alphaStat != alpha || betaStat != beta) {
    alphaStat = alpha; betaStat = beta;
    if (betaDist) delete betaDist;
    betaDist = new beta_dist(alphaStat, betaStat);
  }
}


inline Real BetaRandomVariable::std_pdf(Real x, Real alpha, Real beta)
{
  beta_dist beta1(alpha, beta);
  return bmth::pdf(beta1, x);
}


inline Real BetaRandomVariable::std_cdf(Real x, Real alpha, Real beta)
{
  beta_dist beta1(alpha, beta);
  return bmth::cdf(beta1, x);
}


inline Real BetaRandomVariable::inverse_std_cdf(Real cdf, Real alpha, Real beta)
{
  beta_dist beta1(alpha, beta);
  return bmth::quantile(beta1, cdf);
}


inline Real BetaRandomVariable::
pdf(Real x, Real alpha, Real beta, Real lwr, Real upr)
{
  Real scale = upr - lwr, scaled_x = (x - lwr)/scale;
  return std_pdf(scaled_x, alpha, beta) / scale;
}


inline Real BetaRandomVariable::
cdf(Real x, Real alpha, Real beta, Real lwr, Real upr)
{
  Real scaled_x = (x - lwr)/(upr - lwr);
  return std_cdf(scaled_x, alpha, beta);
}


inline void BetaRandomVariable::
moments_from_params(Real lwr, Real upr, Real alpha, Real beta,
		    Real& mean, Real& std_dev)
{
  Real range = upr - lwr;
  mean       = lwr + alpha/(alpha+beta)*range;
  std_dev    = std::sqrt(alpha*beta/(alpha+beta+1.))/(alpha+beta)*range;
}

} // namespace Pecos

#endif

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
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  //Real pdf_hessian(Real x) const;

  Real inverse_log_cdf(Real log_p) const;
  Real log_pdf(Real x) const;

  void update(Real alpha, Real beta);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real alpha, Real beta);
  static Real cdf(Real x, Real alpha, Real beta);

  static void moments_from_params(Real alpha, Real beta,
				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Data
  //

  /// alpha parameter of frechet random variable (shape)
  Real alphaStat;
  /// beta parameter of frechet random variable (scale)
  Real betaStat;
};


inline FrechetRandomVariable::FrechetRandomVariable():
  RandomVariable(BaseConstructor()), alphaStat(1.), betaStat(1.)
{ }


inline FrechetRandomVariable::FrechetRandomVariable(Real alpha, Real beta):
  RandomVariable(BaseConstructor()), alphaStat(alpha), betaStat(beta)
{ }


inline FrechetRandomVariable::~FrechetRandomVariable()
{ }


inline Real FrechetRandomVariable::cdf(Real x) const
{ return std::exp(-std::pow(betaStat/x, alphaStat)); }


inline Real FrechetRandomVariable::ccdf(Real x) const
{ return -bmth::expm1(-std::pow(betaStat/x, alphaStat)); }


inline Real FrechetRandomVariable::inverse_cdf(Real p_cdf) const
{
  // p = std::exp(-std::pow(beta/x, alpha))
  // beta/x = std::pow( -log p, 1/alpha)
  return betaStat * std::pow(-std::log(p_cdf), -1./alphaStat);
}


inline Real FrechetRandomVariable::inverse_ccdf(Real p_ccdf) const
{ return betaStat * std::pow(-bmth::log1p(-p_ccdf), -1./alphaStat); }


inline Real FrechetRandomVariable::inverse_log_cdf(Real log_p) const
{ return betaStat * std::pow(-log_p, -1./alphaStat); }


inline Real FrechetRandomVariable::pdf(Real x) const
{
  Real num = std::pow(betaStat/x, alphaStat);
  return alphaStat/x*num*std::exp(-num);
}


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
//  return pdf(x, alphaStat, betaStat) * ...; // TO DO
//}


inline Real FrechetRandomVariable::log_pdf(Real x) const
{
  Real num = std::pow(betaStat/x, alphaStat);
  return std::log(alphaStat/x*num) - num; // fewer operations

  // more decomposed -> less likelihood of overflow?
  //Real num = betaStat/x;
  //return std::log(alphaStat/x) + alphaStat * std::log(num)
  //  - std::pow(num, alphaStat);
}


inline void FrechetRandomVariable::update(Real alpha, Real beta)
{ alphaStat = alpha; betaStat = beta; }

// static functions:

inline Real FrechetRandomVariable::pdf(Real x, Real alpha, Real beta)
{
  Real num = std::pow(beta/x, alpha);
  return alpha/x*num*std::exp(-num);
}


inline Real FrechetRandomVariable::cdf(Real x, Real alpha, Real beta)
{ return std::exp(-std::pow(beta/x, alpha)); }


inline void FrechetRandomVariable::
moments_from_params(Real alpha, Real beta, Real& mean, Real& std_dev)
{
  // See Haldar and Mahadevan, p. 91-92
  Real gam = bmth::tgamma(1.-1./alpha);
  mean     = beta*gam;
  std_dev  = beta*std::sqrt(bmth::tgamma(1.-2./alpha)-gam*gam);
}

} // namespace Pecos

#endif

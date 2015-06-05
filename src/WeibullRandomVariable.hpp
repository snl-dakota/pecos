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
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  //Real pdf_hessian(Real x) const;

  Real inverse_log_ccdf(Real log_p_ccdf) const;
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

  /// alpha parameter of weibull random variable (shape)
  Real alphaStat;
  /// beta parameter of weibull random variable (scale)
  Real betaStat;

  /// pointer to the Boost weibull_distribution instance
  weibull_dist* weibullDist;
};


inline WeibullRandomVariable::WeibullRandomVariable():
  RandomVariable(BaseConstructor()), alphaStat(1.), betaStat(1.),
  weibullDist(NULL)
{ }


inline WeibullRandomVariable::WeibullRandomVariable(Real alpha, Real beta):
  RandomVariable(BaseConstructor()), alphaStat(alpha), betaStat(beta),
  weibullDist(new weibull_dist(alpha, beta))
{ }


inline WeibullRandomVariable::~WeibullRandomVariable()
{ }


inline Real WeibullRandomVariable::cdf(Real x) const
{ return bmth::cdf(*weibullDist, x); }


inline Real WeibullRandomVariable::ccdf(Real x) const
{ return bmth::cdf(complement(*weibullDist, x)); }


inline Real WeibullRandomVariable::inverse_cdf(Real p_cdf) const
{
  return bmth::quantile(*weibullDist, p_cdf);
  //return betaStat * std::pow(-bmth::log1p(-p), 1./alphaStat);
}


inline Real WeibullRandomVariable::inverse_ccdf(Real p_ccdf) const
{
  return bmth::quantile(complement(*weibullDist, p_ccdf));
  //return betaStat * std::pow(-std::log(p_ccdf), 1./alphaStat);
}


inline Real WeibullRandomVariable::inverse_log_ccdf(Real log_p_ccdf) const
{ return betaStat * std::pow(-log_p_ccdf, 1./alphaStat); }


inline Real WeibullRandomVariable::pdf(Real x) const
{
  return bmth::pdf(*weibullDist, x);
  //return alpha/beta * std::pow(x/beta,alpha-1.) *
  //  std::exp(-std::pow(x/beta,alpha));
}


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
//  return pdf(x) * ...; // TO DO
//}


inline Real WeibullRandomVariable::log_pdf(Real x) const
{
  Real num = x/betaStat;
  return std::log(alphaStat/betaStat) + (alphaStat-1.) * std::log(num)
    - std::pow(num,alphaStat);
}


inline void WeibullRandomVariable::update(Real alpha, Real beta)
{
  if (!weibullDist || alphaStat != alpha || betaStat != beta) {
    alphaStat = alpha; betaStat = beta;
    if (weibullDist) delete weibullDist;
    weibullDist = new weibull_dist(alphaStat, betaStat);
  }
}

// Static member functions:

inline Real WeibullRandomVariable::pdf(Real x, Real alpha, Real beta)
{
  weibull_dist weibull1(alpha, beta);
  return bmth::pdf(weibull1, x);
  //return alpha/beta * std::pow(x/beta,alpha-1.) *
  //  std::exp(-std::pow(x/beta,alpha));
}


inline Real WeibullRandomVariable::cdf(Real x, Real alpha, Real beta)
{
  weibull_dist weibull1(alpha, beta);
  return bmth::cdf(weibull1, x);
  // avoid numerical probs when exp()~1
  //return -std::expm1(-std::pow(x/beta, alpha));
}


inline void WeibullRandomVariable::
moments_from_params(Real alpha, Real beta, Real& mean, Real& std_dev)
{
  // See Haldar and Mahadevan, p. 97
  Real gam = bmth::tgamma(1.+1./alpha),
    cf_var = std::sqrt(bmth::tgamma(1.+2./alpha)/gam/gam - 1.);
  mean     = beta*gam;
  std_dev  = cf_var*mean;
}

} // namespace Pecos

#endif

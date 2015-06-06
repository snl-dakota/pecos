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
  //- Heading: Virtual function redefinitions
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

  Real coefficient_of_variation() const;
  Real correlation_warping_factor(const RandomVariable& rv, Real corr) const;

  //
  //- Heading: Member functions
  //

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


inline Real WeibullRandomVariable::coefficient_of_variation() const
{
  Real mean, std_dev;
  moments_from_params(alphaStat, betaStat, mean, std_dev);
  return std_dev / mean;
}


inline Real WeibullRandomVariable::
correlation_warping_factor(const RandomVariable& rv, Real corr) const
{
  // correlation warping factor for transformations to STD_NORMAL space
  // Der Kiureghian and Liu, ASCE JEM 112:1, 1986
  Real COV = coefficient_of_variation(), COV_rv;
  switch (rv.type()) { // x-space types mapped to STD_NORMAL u-space

  // Der Kiureghian & Liu: Table 3 (quadratic approximations in COV)
  case NORMAL:      // Max Error 0.1%
    return 1.031 + (-0.195 + 0.328*COV)*COV; break;

  // Der Kiureghian & Liu: Table 5 (quadratic approximations in COV,corr)
  case UNIFORM:     // Max Error 0.5%
    return 1.061 + (-0.237 + 0.379*COV)*COV - 0.005*corr*corr; break;
  case EXPONENTIAL: // Max Error 0.4%
    return 1.147 + (0.145 + 0.010*corr)*corr
      + (-0.271 + 0.459*COV - 0.467*corr)*COV; break;
  case GUMBEL:      // Max Error 0.2%
    return 1.064 + (0.065 +  0.003*corr)*corr
      + (-0.21 + 0.356*COV - 0.211*corr)*COV; break;

  // Der Kiureghian & Liu: Table 6 (quadratic approximations in COV)
  case LOGNORMAL: // Max Error 2.4%
    COV_rv = rv.coefficient_of_variation();
    return 1.031 + (0.052 + 0.002*corr)*corr
      + (0.011 + 0.22*COV_rv + 0.005*corr)*COV_rv
      + (-0.21 + 0.35*COV + 0.009*COV_rv - 0.174*corr)*COV; break;
  case GAMMA:     // Max Error 4.0%
    COV = rv.coefficient_of_variation();
    return 1.032 + 0.034*corr + (-0.202 + 0.339*COV - 0.111*corr)*COV
      + (-0.007 + 0.121*COV_rv - 0.006*corr + 0.003*COV)*COV_rv; break;
  case FRECHET:   // Max Error 3.8%
    COV_rv = rv.coefficient_of_variation();
    return 1.065 + (0.146 + 0.013*corr)*corr
      + (-0.259 + 0.435*COV + 0.034*COV_rv - 0.481*corr)*COV
      + ( 0.241 + 0.372*COV_rv + 0.005*corr)*COV_rv; break;
  case WEIBULL:   // Max Error 2.6%
    COV_rv = rv.coefficient_of_variation();
    return 1.063 + (-0.004 - 0.001*corr)*corr - 0.007*COV*COV_rv
      + (COV + COV_rv) * (0.007*corr - 0.2) + 0.337*(COV*COV + COV_rv*COV_rv);
    break;

  default: // Unsupported warping (should be prevented upsteam)
    PCerr << "Error: unsupported correlation warping for WeibullRV."<<std::endl;
    abort_handler(-1); return 1.; break;
  }
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
       COV = std::sqrt(bmth::tgamma(1.+2./alpha)/gam/gam - 1.);
  mean     = beta*gam;
  std_dev  = COV*mean;
}

} // namespace Pecos

#endif

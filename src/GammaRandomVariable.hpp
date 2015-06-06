/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 GammaRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef GAMMA_RANDOM_VARIABLE_HPP
#define GAMMA_RANDOM_VARIABLE_HPP

#include "ExponentialRandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for gamma random variables.

/** Manages alpha and inherits beta.  This follows the
    Gamma(alpha,beta) definition in Law & Kelton, which differs from
    the LHS definition (beta_LK = 1/beta_LHS). */

class GammaRandomVariable: public ExponentialRandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  GammaRandomVariable();
  /// alternate constructor
  GammaRandomVariable(Real alpha, Real beta);
  /// destructor
  ~GammaRandomVariable();

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

  // inherited from ExponentialRandomVariable
  //Real to_std(Real x) const;
  //Real from_std(Real z) const;

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
  static Real inverse_cdf(Real cdf, Real alpha, Real beta);

  static void moments_from_params(Real alpha, Real beta,
				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Data
  //

  /// alpha parameter of gamma random variable (statistical PDF
  /// convention; differs from GenLaguerre polynomial convention)
  Real alphaStat;

  /// pointer to the Boost gamma_distribution instance
  gamma_dist* gammaDist;
};


inline GammaRandomVariable::GammaRandomVariable():
  ExponentialRandomVariable(), alphaStat(1.), gammaDist(NULL)
{ }


inline GammaRandomVariable::GammaRandomVariable(Real alpha, Real beta):
  ExponentialRandomVariable(beta), alphaStat(alpha),
  gammaDist(new gamma_dist(alphaStat, betaStat))
{ }


inline GammaRandomVariable::~GammaRandomVariable()
{ if (gammaDist) delete gammaDist; }


inline Real GammaRandomVariable::cdf(Real x) const
{ return bmth::cdf(*gammaDist, x); }


inline Real GammaRandomVariable::ccdf(Real x) const
{ return bmth::cdf(complement(*gammaDist, x)); }


inline Real GammaRandomVariable::inverse_cdf(Real p_cdf) const
{ return bmth::quantile(*gammaDist, p_cdf); }


inline Real GammaRandomVariable::inverse_ccdf(Real p_ccdf) const
{ return bmth::quantile(complement(*gammaDist, p_ccdf)); }


inline Real GammaRandomVariable::pdf(Real x) const
{ return bmth::pdf(*gammaDist, x); }


inline Real GammaRandomVariable::pdf_gradient(Real x) const
{
  return std::pow(betaStat,-alphaStat) / bmth::tgamma(alphaStat) *
    (std::exp(-x/betaStat)*(alphaStat-1.) * std::pow(x,alphaStat-2.) -
     std::pow(x,alphaStat-1.) * std::exp(-x/betaStat)/betaStat);
}


//inline Real GammaRandomVariable::pdf_hessian(Real x) const
//{
//  return gamma_pdf(x, alphaStat, betaStat); // * ...; // TO DO
//}


inline Real GammaRandomVariable::coefficient_of_variation() const
{
  Real mean, std_dev;
  moments_from_params(alphaStat, betaStat, mean, std_dev);
  return std_dev / mean;
}


inline Real GammaRandomVariable::
correlation_warping_factor(const RandomVariable& rv, Real corr) const
{
  // correlation warping factor for transformations to STD_NORMAL space
  // Der Kiureghian and Liu, ASCE JEM 112:1, 1986
  Real COV = coefficient_of_variation(), COV_rv;
  switch (rv.type()) { // x-space types mapped to STD_NORMAL u-space

  // Der Kiureghian & Liu: Table 3 (quadratic approximations in COV)
  case NORMAL: // Max Error 0.0%
    return 1.001 + (-0.007 + 0.118*COV)*COV; break;
  case UNIFORM: // Max Error 0.1%
    return 1.023 + (-0.007 + 0.127*COV)*COV + 0.002*corr*corr; break;
  case EXPONENTIAL: // Max Error 0.9%
    return 1.104 + (0.003 + 0.014*corr)*corr
      + (-0.008 + 0.173*COV - 0.296*corr)*COV; break;
  case GUMBEL: // Max Error 0.3%
    return 1.031 + (0.001 +  0.003*corr)*corr
      + (-0.007 + 0.131*COV - 0.132*corr)*COV; break;

  // Der Kiureghian & Liu: Table 6
  case LOGNORMAL: // Max Error 4.0%
    COV_rv = rv.coefficient_of_variation();
    return 1.001 + (0.033 + 0.002*corr)*corr
      + (0.004 + 0.223*COV_rv - 0.104*corr)*COV_rv
      + (0.016 + 0.13*COV + 0.029*COV_rv - 0.119*corr)*COV; break;
  case GAMMA: // Max Error 4.0%
    COV_rv = rv.coefficient_of_variation();
    return 1.002 + 0.022*corr - 0.012*(COV + COV_rv) + 0.001*corr*corr
      + 0.125*(COV*COV + COV_rv*COV_rv) - 0.077*corr*(COV + COV_rv)
      + 0.014*COV*COV_rv; break;
  case FRECHET: // Max Error 4.2%
    COV_rv = rv.coefficient_of_variation();
    return 1.029 + (0.056 +  0.012*corr)*corr
      + (-0.03 + 0.174*COV - 0.313*corr)*COV
      + (0.225 + 0.379*COV_rv +  0.075*COV - 0.182*corr)*COV_rv; break;
  case WEIBULL: // Max Error 4.0%
    COV_rv = rv.coefficient_of_variation();
    return 1.032 + 0.034*corr
      + (-0.007 + 0.121*COV - 0.006*corr + 0.003*COV_rv)*COV
      + (-0.202 + 0.339*COV_rv - 0.111*corr)*COV_rv; break;

  default: // Unsupported warping (should be prevented upsteam)
    PCerr << "Error: unsupported correlation warping for GammaRV."<< std::endl;
    abort_handler(-1); return 1.; break;
  }
}


inline void GammaRandomVariable::update(Real alpha, Real beta)
{
  if (!gammaDist || alphaStat != alpha || betaStat != beta) {
    alphaStat = alpha; betaStat = beta;
    if (gammaDist) delete gammaDist;
    gammaDist = new gamma_dist(alphaStat, betaStat);
  }
}


inline Real GammaRandomVariable::pdf(Real x, Real alpha, Real beta)
{
  gamma_dist gamma1(alpha, beta);
  return bmth::pdf(gamma1, x);
}


inline Real GammaRandomVariable::cdf(Real x, Real alpha, Real beta)
{
  gamma_dist gamma1(alpha, beta);
  return bmth::cdf(gamma1, x);
}


inline Real GammaRandomVariable::inverse_cdf(Real cdf, Real alpha, Real beta)
{
  gamma_dist gamma1(alpha, beta);
  return bmth::quantile(gamma1, cdf);
}


inline void GammaRandomVariable::
moments_from_params(Real alpha, Real beta, Real& mean, Real& std_dev)
{ mean = alpha*beta; std_dev = std::sqrt(alpha)*beta; }

} // namespace Pecos

#endif

/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 LognormalRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef LOGNORMAL_RANDOM_VARIABLE_HPP
#define LOGNORMAL_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for lognormal random variables.

/** Manages lambda and zeta (mean and std deviation of underlying
    normal distribution). */

class LognormalRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  LognormalRandomVariable();                      ///< constructor
  LognormalRandomVariable(Real lambda, Real zeta); ///< alternate constructor
  ~LognormalRandomVariable();                     ///< destructor

  //
  //- Heading: Virtual function redefinitions
  //

  Real cdf(Real x) const;
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;
  //Real pdf_gradient(Real x) const;
  //Real pdf_hessian(Real x) const;

  Real log_pdf(Real x) const;

  Real coefficient_of_variation() const;
  Real correlation_warping_factor(const RandomVariable& rv, Real corr) const;

  //
  //- Heading: Member functions
  //

  void update(Real lambda, Real zeta);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real lambda, Real zeta);
  static Real cdf(Real x, Real lambda, Real zeta);

  static void std_deviation_from_err_factor(Real mean, Real err_fact,
					    Real& std_dev);
  static void err_factor_from_std_deviation(Real mean, Real std_dev,
					    Real& err_fact);
  static void moments_from_params(Real lambda, Real zeta, Real& mean,
				  Real& std_dev);
  static void zeta_sq_from_moments(Real mean, Real std_dev, Real& zeta_sq);
  static void params_from_moments(Real mean, Real std_dev, Real& lambda,
				  Real& zeta);

protected:

  //
  //- Heading: Data
  //

  /// lambda parameter for lognormal random variable
  Real lnLambda;
  /// zeta parameter for lognormal random variable
  Real lnZeta;
};


inline LognormalRandomVariable::LognormalRandomVariable():
  RandomVariable(BaseConstructor())
{ params_from_moments(0., 1., lnLambda, lnZeta); }


inline LognormalRandomVariable::LognormalRandomVariable(Real lambda, Real zeta):
  RandomVariable(BaseConstructor()), lnLambda(lambda), lnZeta(zeta)
{ }


inline LognormalRandomVariable::~LognormalRandomVariable()
{ }


inline Real LognormalRandomVariable::cdf(Real x) const
{
  lognormal_dist logn1(lnLambda, lnZeta);
  return bmth::cdf(logn1, x);

  //return NormalRandomVariable::std_cdf((std::log(x) - lnLambda)/lnZeta);
}


inline Real LognormalRandomVariable::ccdf(Real x) const
{
  lognormal_dist logn1(lnLambda, lnZeta);
  return bmth::cdf(complement(logn1, x));

  //return NormalRandomVariable::std_ccdf((std::log(x) - lnLambda)/lnZeta);
}


inline Real LognormalRandomVariable::inverse_cdf(Real p_cdf) const
{
  lognormal_dist logn1(lnLambda, lnZeta);
  return bmth::quantile(logn1, p_cdf);

  // p = Phi((std::log(x) - lambda)/zeta)
  //return std::exp(NormalRandomVariable::inverse_std_cdf(p_cdf) *
  //		    lnZeta + lnLambda);
}


inline Real LognormalRandomVariable::inverse_ccdf(Real p_ccdf) const
{
  lognormal_dist logn1(lnLambda, lnZeta);
  return bmth::quantile(complement(logn1, p_ccdf));

  //return std::exp(NormalRandomVariable::inverse_std_ccdf(p_ccdf) *
  //		    lnZeta + lnLambda);
}


inline Real LognormalRandomVariable::pdf(Real x) const
{
  lognormal_dist logn1(lnLambda, lnZeta);
  return bmth::pdf(logn1, x);
}


//inline Real LognormalRandomVariable::pdf_gradient(Real x) const
//{
//  return pdf(x) * ...; // TO DO
//}


//inline Real LognormalRandomVariable::pdf_hessian(Real x) const
//{
//  return pdf(x) * ...; // TO DO
//}


inline Real LognormalRandomVariable::log_pdf(Real x) const
{
  Real num = (std::log(x) - lnLambda) / lnZeta;
  return -std::log(std::sqrt(2.*PI) * lnZeta * x) - num*num/2.;
}


inline Real LognormalRandomVariable::coefficient_of_variation() const
{
  Real mean, std_dev;
  moments_from_params(lnLambda, lnZeta, mean, std_dev);
  return std_dev / mean;
}


inline Real LognormalRandomVariable::
correlation_warping_factor(const RandomVariable& rv, Real corr) const
{
  // correlation warping factor for transformations to STD_NORMAL space
  // Der Kiureghian and Liu, ASCE JEM 112:1, 1986
  Real COV = coefficient_of_variation(), COV_rv;
  switch (rv.type()) { // x-space types mapped to STD_NORMAL u-space

  // Der Kiureghian & Liu: Table 3
  case NORMAL:
    return COV/std::sqrt(bmth::log1p(COV*COV)); break; // Exact

  // Der Kiureghian & Liu: Table 5 (quadratic approx in corr,COV)
  case UNIFORM: // Max Error 0.7%
    return 1.019 + (0.014 + 0.249*COV)*COV + 0.01*corr*corr; break;
  case EXPONENTIAL: // Max Error 1.6%
    return 1.098 + (0.003 + 0.025*corr)*corr
      + ( 0.019 + 0.303*COV - 0.437*corr)*COV; break;
  case GUMBEL: // Max Error 0.3%
    return 1.029 + (0.001 +  0.004*corr)*corr
      + ( 0.014 + 0.233*COV - 0.197*corr)*COV; break;

  // Der Kiureghian & Liu: Table 6 (quadratic approx in corr,COV,COV_rv)
  case LOGNORMAL: // Exact
    COV_rv = rv.coefficient_of_variation();
    return bmth::log1p(COV*COV_rv*corr) / corr /
      std::sqrt(bmth::log1p(COV*COV)*bmth::log1p(COV_rv*COV_rv)); break;
  case GAMMA: // Max Error 4.0%
    COV_rv = rv.coefficient_of_variation();
    return 1.001 + (0.033 + 0.002*corr)*corr
      + (0.004 + 0.223*COV - 0.104*corr)*COV
      + (0.016 + 0.13*COV_rv + 0.029*COV - 0.119*corr)*COV_rv; break;
  case FRECHET: // Max Error 4.3%
    COV_rv = rv.coefficient_of_variation();
    return  1.026 + (0.082 + 0.018*corr)*corr
      + (-0.019 + 0.288*COV - 0.441*corr)*COV
      + ( 0.222 + 0.379*COV_rv + 0.126*COV - 0.277*corr)*COV_rv; break;
  case WEIBULL: // Max Error 2.4%
    COV_rv = rv.coefficient_of_variation();
    return 1.031 + (0.052 + 0.002*corr)*corr
      + (0.011 + 0.22*COV + 0.005*corr)*COV
      + (-0.21 + 0.35*COV_rv +  0.009*COV - 0.174*corr)*COV_rv; break;

  default: // Unsupported warping (should be prevented upsteam)
    PCerr << "Error: unsupported correlation warping for LognormalRV."
	  << std::endl;
    abort_handler(-1); return 1.; break;
  }
}


inline void LognormalRandomVariable::update(Real lambda, Real zeta)
{ lnLambda = lambda; lnZeta = zeta; }

// static functions

inline Real LognormalRandomVariable::pdf(Real x, Real lambda, Real zeta)
{
  lognormal_dist logn1(lambda, zeta);
  return bmth::pdf(logn1, x);

  //return NormalRandomVariable::std_pdf((std::log(x) - lambda)/zeta);
}


inline Real LognormalRandomVariable::cdf(Real x, Real lambda, Real zeta)
{
  lognormal_dist logn1(lambda, zeta);
  return bmth::cdf(logn1, x);

  //return NormalRandomVariable::std_cdf((std::log(x) - lambda)/zeta);
}


inline void LognormalRandomVariable::
std_deviation_from_err_factor(Real mean, Real err_fact, Real& std_dev)
{
  // Phi^{-1}(0.95) ~= 1.645
  Real zeta = std::log(err_fact)/NormalRandomVariable::Phi_inverse(0.95);
  std_dev   = mean * std::sqrt(bmth::expm1(zeta*zeta));
}


inline void LognormalRandomVariable::
err_factor_from_std_deviation(Real mean, Real std_dev, Real& err_fact)
{
  Real COV = std_dev/mean, zeta = std::sqrt(bmth::log1p(COV*COV));
  err_fact = std::exp(NormalRandomVariable::Phi_inverse(0.95)*zeta);
}


inline void LognormalRandomVariable::
moments_from_params(Real lambda, Real zeta, Real& mean, Real& std_dev)
{
  Real zeta_sq = zeta*zeta;
  mean    = std::exp(lambda + zeta_sq/2.);
  std_dev = mean * std::sqrt(bmth::expm1(zeta_sq));
}


inline void LognormalRandomVariable::
zeta_sq_from_moments(Real mean, Real std_dev, Real& zeta_sq)
{
  Real COV = std_dev/mean;
  zeta_sq = bmth::log1p(COV*COV);
}


inline void LognormalRandomVariable::
params_from_moments(Real mean, Real std_dev, Real& lambda, Real& zeta)
{
  Real zeta_sq;
  zeta_sq_from_moments(mean, std_dev, zeta_sq);
  lambda = std::log(mean) - zeta_sq/2.;
  zeta   = std::sqrt(zeta_sq);
}

} // namespace Pecos

#endif

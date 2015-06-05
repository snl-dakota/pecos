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
  //- Heading: Member functions
  //

  Real cdf(Real x) const;
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;
  //Real pdf_gradient(Real x) const;
  //Real pdf_hessian(Real x) const;

  Real log_pdf(Real x) const;

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
  Real cf_var = std_dev/mean, zeta = std::sqrt(std::log(1. + cf_var*cf_var));
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
  Real cf_var = std_dev/mean;
  zeta_sq = std::log(1. + cf_var*cf_var);
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

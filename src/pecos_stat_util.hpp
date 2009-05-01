/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_STAT_UTIL_HPP
#define PECOS_STAT_UTIL_HPP

#include "pecos_data_types.hpp"
#ifdef HAVE_GSL
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_sf_gamma.h"
#endif
#ifdef HAVE_BOOST 
#define BOOST_MATH_PROMOTE_DOUBLE_POLICY false
// WJB - ToDo: investigate error in boost/math/tools/traits.hpp with SunProCC
//#include <boost/math/distributions.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/weibull.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/policies/policy.hpp>
namespace bmth = boost::math;
namespace bmp  = bmth::policies;
#endif
#include <math.h>


namespace Pecos {


// -----------------------------------
// Non-default boost math/policy types
// -----------------------------------
#ifdef HAVE_BOOST
typedef bmth::
  normal_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  normal_dist;
typedef bmth::
  lognormal_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  lognormal_dist;
typedef bmth::
  exponential_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  exponential_dist;
typedef bmth::
  beta_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  beta_dist;
typedef bmth::
  gamma_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  gamma_dist;
typedef bmth::
  weibull_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  weibull_dist;
typedef bmth::
  chi_squared_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  chi_squared_dist;
typedef bmth::
  students_t_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  students_t_dist;
typedef bmth::
  fisher_f_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  fisher_f_dist;
#endif


#ifndef HAVE_BOOST
#ifdef HAVE_GSL
/// Inverse of standard beta CDF by backtracking Newton (not supported by GSL)
Real cdf_beta_Pinv(const Real& cdf, const Real& alpha, const Real& beta);
#else
/// Inverse of error function used in Phi_inverse()
Real erf_inverse(const Real& p);
#endif // HAVE_GSL
#endif // HAVE_BOOST


inline Real phi(const Real& beta)
{
#ifdef HAVE_BOOST
  normal_dist norm(0., 1.);
  return bmth::pdf(norm, beta);
#elif HAVE_GSL
  return gsl_ran_ugaussian_pdf(beta);
#else
  return exp(-beta*beta/2.)/sqrt(2.*Pi);
#endif // HAVE_GSL or HAVE_BOOST
}


/** returns a probability < 0.5 for negative beta and a probability > 0.5
    for positive beta. */
inline Real Phi(const Real& beta)
{
#ifdef HAVE_BOOST
  normal_dist norm(0., 1.);
  return bmth::cdf(norm, beta);
#elif HAVE_GSL
  return gsl_cdf_ugaussian_P(beta);
#else
  return .5 + .5*erf(beta/sqrt(2.));
#endif // HAVE_GSL or HAVE_BOOST
}


/** returns a negative beta for probability < 0.5 and a positive beta for
    probability > 0.5. */
inline Real Phi_inverse(const Real& p)
{
#ifdef HAVE_BOOST
  normal_dist norm(0., 1.);
  return bmth::quantile(norm, p); 
#elif HAVE_GSL
  return gsl_cdf_ugaussian_Pinv(p);
#else
  return sqrt(2.)*erf_inverse(2.*p - 1.);
#endif // HAVE_GSL or HAVE_BOOST
}


inline Real lognormal_pdf(const Real& x, const Real& mean, const Real& std_dev)
{
  // convert mean/std_dev of lognormal to lambda/zeta of underlying normal
  Real cf_var = std_dev/mean, zeta_sq = std::log(1 + cf_var*cf_var),
    lambda = std::log(mean) - zeta_sq/2., zeta = sqrt(zeta_sq);
#ifdef HAVE_BOOST
  lognormal_dist logn1(lambda, zeta);
  return bmth::pdf(logn1, x);
#elif HAVE_GSL
  // GSL gamma passes alpha followed by beta
  return gsl_ran_lognormal_pdf(x, lambda, zeta);
#else
  return 1./sqrt(2.*Pi)/zeta/x *
    std::exp(-std::pow(std::log(x)-lambda, 2)/2./zeta_sq);
#endif // HAVE_GSL or HAVE_BOOST
}


inline Real lognormal_cdf(const Real& x, const Real& mean, const Real& std_dev)
{
  // convert mean/std_dev of lognormal to lambda/zeta of underlying normal
  Real cf_var = std_dev/mean, zeta_sq = std::log(1 + cf_var*cf_var),
    lambda = std::log(mean) - zeta_sq/2., zeta = sqrt(zeta_sq);
  return Phi((std::log(x) - lambda)/zeta);
}


inline Real loguniform_pdf(const Real& x, const Real& lwr, const Real& upr)
{ return 1./(std::log(upr) - std::log(lwr))/x; }


inline Real loguniform_cdf(const Real& x, const Real& lwr, const Real& upr)
{ return (std::log(x) - std::log(lwr))/(std::log(upr) - std::log(lwr)); }


inline Real std_uniform_pdf(const Real& x)
{ return 0.5; } // [-1,1]


inline Real std_uniform_cdf(const Real& x)
{ return (x + 1.)/2.; } // [-1,1]


inline Real uniform_pdf(const Real& x, const Real& lwr, const Real& upr)
{ return 1./(upr - lwr); } // [lwr,upr]


inline Real uniform_cdf(const Real& x, const Real& lwr, const Real& upr)
{ return (x - lwr)/(upr - lwr); } // [lwr,upr]


inline Real std_exponential_pdf(const Real& x)
{ return std::exp(-x); }


inline Real std_exponential_cdf(const Real& x)
{ return 1. - std::exp(-x); }


inline Real exponential_pdf(const Real& x, const Real& beta)
{ return std::exp(-x/beta)/beta; }


inline Real exponential_cdf(const Real& x, const Real& beta)
{ return 1. - std::exp(-x/beta); }


inline Real std_beta_pdf(const Real& x, const Real& alpha, const Real& beta)
{
#ifdef HAVE_BOOST
  beta_dist beta1(alpha, beta);
  return bmth::pdf(beta1, x);
#elif HAVE_GSL
  // GSL beta passes alpha followed by beta
  return gsl_ran_beta_pdf(x, alpha, beta);
#endif // HAVE_GSL or BOOST
}


inline Real std_beta_cdf(const Real& x, const Real& alpha, const Real& beta)
{
#ifdef HAVE_BOOST
  beta_dist beta1(alpha, beta);
  return bmth::cdf(beta1, x);
#elif HAVE_GSL
  // GSL beta passes alpha followed by beta
  return gsl_cdf_beta_P(x, alpha, beta);
#endif //HAVE_GSL or BOOST
}


inline Real std_beta_cdf_inverse(const Real& cdf, const Real& alpha,
				 const Real& beta)
{
#ifdef HAVE_BOOST
  beta_dist beta1(alpha, beta);
  return bmth::quantile(beta1, cdf);
#elif HAVE_GSL
  // GSL does not support beta CDF inverse, therefore a special Newton solve
  // has been implemented for the inversion (which uses the GSL CDF/PDF fns).
  return cdf_beta_Pinv(cdf, alpha, beta);
#endif // HAVE_BOOST
}


inline Real gamma_pdf(const Real& x, const Real& alpha, const Real& beta)
{
#ifdef HAVE_BOOST
  gamma_dist gamma1(alpha, beta);
  return bmth::pdf(gamma1, x);
#elif HAVE_GSL
  // GSL gamma passes alpha followed by beta
  return gsl_ran_gamma_pdf(x, alpha, beta);
#endif // HAVE_GSL
}


inline Real gamma_pdf_deriv(const Real& x, const Real& alpha, const Real& beta)
{
#ifdef HAVE_BOOST
  Real gamma_alpha = bmth::tgamma(alpha);
#elif HAVE_GSL
  Real gamma_alpha = gsl_sf_gamma(alpha);
#endif // HAVE_GSL or BOOST
  return pow(beta,-alpha) / gamma_alpha * (exp(-x/beta)*(alpha-1.) *
    pow(x,alpha-2.) - pow(x,alpha-1.) * exp(-x/beta)/beta);
}


inline Real gamma_cdf(const Real& x, const Real& alpha, const Real& beta)
{
#ifdef HAVE_BOOST
  gamma_dist gamma1(alpha, beta);
  return bmth::cdf(gamma1, x);
#elif HAVE_GSL
  // GSL gamma passes alpha followed by beta
  return gsl_cdf_gamma_P(x, alpha, beta);
#endif // HAVE_GSL
}


inline Real gamma_cdf_inverse(const Real& cdf, const Real& alpha,
			      const Real& beta)
{
#ifdef HAVE_BOOST
  gamma_dist gamma1(alpha, beta);
  return bmth::quantile(gamma1, cdf);
#elif HAVE_GSL
  // GSL gamma passes alpha followed by beta
  return gsl_cdf_gamma_Pinv(cdf, alpha, beta);
#endif // HAVE_GSL or HAVE_BOOST
}


inline Real gumbel_pdf(const Real& x, const Real& alpha, const Real& beta)
{
  Real num = std::exp(-alpha*(x-beta));
  return alpha*num*std::exp(-num);
}


inline Real gumbel_cdf(const Real& x, const Real& alpha, const Real& beta)
{ return std::exp(-std::exp(-alpha * (x - beta))); }


inline Real frechet_pdf(const Real& x, const Real& alpha, const Real& beta)
{
  Real num = beta/x;
  return alpha/beta*std::pow(num,alpha+1.)*std::exp(-std::pow(num,alpha));
}


inline Real frechet_cdf(const Real& x, const Real& alpha, const Real& beta)
{ return std::exp(-std::pow(beta/x, alpha)); }


inline Real weibull_pdf(const Real& x, const Real& alpha, const Real& beta)
{
#ifdef HAVE_BOOST
  weibull_dist weibull1(alpha, beta);
  return bmth::pdf(weibull1, x);
#elif HAVE_GSL
  // GSL weibull passes beta followed by alpha
  return gsl_ran_weibull_pdf(x, beta, alpha);
#else
  return alpha/beta * std::pow(x/beta,alpha-1.) *
    std::exp(-std::pow(x/beta,alpha));
#endif // HAVE_GSL or BOOST
}


inline Real weibull_cdf(const Real& x, const Real& alpha, const Real& beta)
{
#ifdef HAVE_BOOST
  weibull_dist weibull1(alpha, beta);
  return bmth::cdf(weibull1, x);
#elif HAVE_GSL
  // GSL weibull passes beta followed by alpha
  return gsl_cdf_weibull_P(x, beta, alpha);
#else
  // avoid numerical probs when exp()~1
  return -std::expm1(-std::pow(x/beta, alpha));
#endif // HAVE_GSL or BOOST
}


inline void moments_from_lognormal_params(const Real& mean,
					  const Real& err_fact, Real& std_dev)
{
  Real zeta = log(err_fact)/1.645;
  std_dev   = mean*sqrt(exp(zeta*zeta)-1.);
}


inline void moments_from_uniform_params(const Real& lwr, const Real& upr,
					Real& mean, Real& std_dev)
{ mean = (lwr + upr)/2.; std_dev = (upr - lwr)/sqrt(12.); }


inline void moments_from_loguniform_params(const Real& lwr, const Real& upr,
					   Real& mean, Real& std_dev)
{
  Real range = upr - lwr, log_range = log(upr) - log(lwr);
  mean       = range/log_range;
  std_dev    = sqrt(range*(log_range*(upr+lwr)-2.*range)/2.)/log_range;
}


inline void moments_from_triangular_params(const Real& lwr, const Real& upr,
					   const Real& mode, Real& mean,
					   Real& std_dev)
{
  mean    = (lwr + mode + upr)/3.;
  std_dev = sqrt((lwr*(lwr - mode) + mode*(mode - upr) + upr*(upr - lwr))/18.);
}


inline void moments_from_exponential_params(const Real& beta, Real& mean,
					    Real& std_dev)
{ mean = beta; std_dev = beta; }


inline void moments_from_beta_params(const Real& lwr, const Real& upr,
				     const Real& alpha, const Real& beta,
				     Real& mean, Real& std_dev)
{
  Real range = upr - lwr;
  mean       = lwr + alpha/(alpha+beta)*range;
  std_dev    = sqrt(alpha*beta/(alpha+beta+1.))/(alpha+beta)*range;
}


inline void moments_from_gamma_params(const Real& alpha, const Real& beta,
				      Real& mean, Real& std_dev)
{ mean = alpha*beta; std_dev = sqrt(alpha)*beta; }


inline void moments_from_gumbel_params(const Real& alpha, const Real& beta,
				       Real& mean, Real& std_dev)
{ mean = beta + 0.5772/alpha; std_dev = Pi/sqrt(6.)/alpha; }


inline void moments_from_frechet_params(const Real& alpha, const Real& beta,
					Real& mean, Real& std_dev)
{
  // See Haldar and Mahadevan, p. 91-92
#ifdef HAVE_BOOST
  Real gam = bmth::tgamma(1.-1./alpha);
  std_dev = beta*sqrt(bmth::tgamma(1.-2./alpha)-gam*gam);
#elif HAVE_GSL
  Real gam = gsl_sf_gamma(1.-1./alpha);
  std_dev = beta*sqrt(gsl_sf_gamma(1.-2./alpha)-gam*gam);
#else
  PCerr << "Error: frechet distributions only supported in executables "
	<< "configured with the GSL or Boost library." << std::endl;
  abort_handler(-1);
#endif // HAVE_GSL or HAVE_BOOST
  mean = beta*gam;
}


inline void moments_from_weibull_params(const Real& alpha, const Real& beta,
					Real& mean, Real& std_dev)
{
  // See Haldar and Mahadevan, p. 97
#ifdef HAVE_BOOST
  Real gam = bmth::tgamma(1.+1./alpha),
    cf_var = sqrt(bmth::tgamma(1.+2./alpha)/gam/gam - 1.);
#elif HAVE_GSL
  Real gam = gsl_sf_gamma(1.+1./alpha),
    cf_var = sqrt(gsl_sf_gamma(1.+2./alpha)/gam/gam - 1.);
#else
  PCerr << "Error: weibull distributions only supported in executables "
	<< "configured with the GSL or Boost library." << std::endl;
  abort_handler(-1);
#endif // HAVE_GSL or HAVE_BOOST
  mean    = beta*gam;
  std_dev = cf_var*mean;
}

} // namespace Pecos

#endif

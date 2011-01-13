/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_STAT_UTIL_HPP
#define PECOS_STAT_UTIL_HPP

#include "pecos_data_types.hpp"
/*
#ifdef HAVE_GSL
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_sf_gamma.h"
#endif 
*/

#define BOOST_MATH_PROMOTE_DOUBLE_POLICY false
#include <boost/math/distributions.hpp>
#include <boost/math/special_functions/expm1.hpp>
using boost::math::expm1;

namespace bmth = boost::math;
namespace bmp  = bmth::policies;


namespace Pecos {


// -----------------------------------
// Non-default boost math/policy types
// -----------------------------------

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


inline Real gamma_function(const Real& x)
{
//#ifdef HAVE_BOOST
  return bmth::tgamma(x);
/*
#elif HAVE_GSL
  return gsl_sf_gamma(x);
#else
  PCerr << "Error: gamma function only supported in executables configured "
	<< "with the GSL or Boost library." << std::endl;
  abort_handler(-1);
#endif // HAVE_GSL or HAVE_BOOST
*/
}


inline Real phi(const Real& beta)
{
//#ifdef HAVE_BOOST
  normal_dist norm(0., 1.);
  return bmth::pdf(norm, beta);
/*
#elif HAVE_GSL
  return gsl_ran_ugaussian_pdf(beta);
#else
  return std::exp(-beta*beta/2.)/std::sqrt(2.*PI);
#endif
*/
}


/** returns a probability < 0.5 for negative beta and a probability > 0.5
    for positive beta. */
inline Real Phi(const Real& beta)
{
//#ifdef HAVE_BOOST
  normal_dist norm(0., 1.);
  return bmth::cdf(norm, beta);
/*
#elif HAVE_GSL
  return gsl_cdf_ugaussian_P(beta);
#else
  return .5 + .5*erf(beta/std::sqrt(2.));
#endif
*/
}


/** returns a negative beta for probability < 0.5 and a positive beta for
    probability > 0.5. */
inline Real Phi_inverse(const Real& p)
{
//#ifdef HAVE_BOOST
  normal_dist norm(0., 1.);
  return bmth::quantile(norm, p); 
/*
#elif HAVE_GSL
  return gsl_cdf_ugaussian_Pinv(p);
#else
  return std::sqrt(2.)*erf_inverse(2.*p - 1.);
#endif // HAVE_GSL or HAVE_BOOST
*/
}


/*
#ifndef HAVE_BOOST
#ifdef HAVE_GSL
/// Inverse of standard beta CDF by backtracking Newton (not supported by GSL)
Real cdf_beta_Pinv(const Real& cdf, const Real& alpha, const Real& beta);
#else
/// Inverse of error function used in Phi_inverse()
Real erf_inverse(const Real& p);
#endif // HAVE_GSL
#endif // HAVE_BOOST
*/


inline void lognormal_std_deviation_from_err_factor(const Real& mean,
						    const Real& err_fact,
						    Real& std_dev)
{
  Real zeta = std::log(err_fact)/Phi_inverse(0.95); // Phi^{-1}(0.95) ~= 1.645
  std_dev   = mean * std::sqrt(expm1(zeta*zeta));
}


inline void lognormal_err_factor_from_std_deviation(const Real& mean,
						    const Real& std_dev,
						    Real& err_fact)
{
  Real cf_var = std_dev/mean, zeta = std::sqrt(std::log(1. + cf_var*cf_var));
  err_fact = std::exp(Phi_inverse(0.95)*zeta);
}


inline void moments_from_lognormal_params(const Real& lambda, const Real& zeta,
					  Real& mean, Real& std_dev)
{
  Real zeta_sq = zeta*zeta;
  mean    = std::exp(lambda + zeta_sq/2.);
  std_dev = mean * std::sqrt(expm1(zeta_sq));
}


inline void lognormal_zeta_sq_from_moments(const Real& mean,
					   const Real& std_dev,
					   Real& zeta_sq)
{
  Real cf_var = std_dev/mean;
  zeta_sq = std::log(1. + cf_var*cf_var);
}


inline void lognormal_params_from_moments(const Real& mean, const Real& std_dev,
					  Real& lambda, Real& zeta)
{
  Real zeta_sq;
  lognormal_zeta_sq_from_moments(mean, std_dev, zeta_sq);
  lambda = std::log(mean) - zeta_sq/2.;
  zeta   = std::sqrt(zeta_sq);
}


inline void moments_from_lognormal_spec(const RealVector& ln_means,
					const RealVector& ln_std_devs,
					const RealVector& ln_lambdas,
					const RealVector& ln_zetas,
					const RealVector& ln_err_facts, 
					size_t i, Real& mean, Real& std_dev)
{
  if (!ln_lambdas.empty())
    moments_from_lognormal_params(ln_lambdas[i], ln_zetas[i], mean, std_dev);
  else {
    mean = ln_means[i];
    if (!ln_std_devs.empty())
      std_dev = ln_std_devs[i];
    else
      lognormal_std_deviation_from_err_factor(mean, ln_err_facts[i], std_dev);
  }
}


inline void params_from_lognormal_spec(const RealVector& ln_means,
				       const RealVector& ln_std_devs,
				       const RealVector& ln_lambdas,
				       const RealVector& ln_zetas,
				       const RealVector& ln_err_facts, 
				       size_t i, Real& lambda, Real& zeta)
{
  if (!ln_lambdas.empty()) {
    lambda = ln_lambdas[i];
    zeta   = ln_zetas[i];
  }
  else {
    if (!ln_std_devs.empty())
      lognormal_params_from_moments(ln_means[i], ln_std_devs[i], lambda, zeta);
    else {
      Real mean = ln_means[i], stdev;
      lognormal_std_deviation_from_err_factor(mean, ln_err_facts[i], stdev);
      lognormal_params_from_moments(mean, stdev, lambda, zeta);
    }
  }
}


inline void all_from_lognormal_spec(const RealVector& ln_means,
				    const RealVector& ln_std_devs,
				    const RealVector& ln_lambdas,
				    const RealVector& ln_zetas,
				    const RealVector& ln_err_facts, 
				    size_t i, Real& mean, Real& std_dev,
				    Real& lambda, Real& zeta, Real& err_fact)
{
  if (!ln_lambdas.empty()) { // lambda/zeta -> mean/std_dev
    lambda = ln_lambdas[i]; zeta = ln_zetas[i];
    moments_from_lognormal_params(lambda, zeta, mean, std_dev);
    lognormal_err_factor_from_std_deviation(mean, std_dev, err_fact);
  }
  else if (!ln_std_devs.empty()) {
    mean = ln_means[i]; std_dev = ln_std_devs[i];
    lognormal_params_from_moments(mean, std_dev, lambda, zeta);
    lognormal_err_factor_from_std_deviation(mean, std_dev, err_fact);
  }
  else { // mean/err_fact -> mean/std_dev
    mean = ln_means[i]; err_fact = ln_err_facts[i];
    lognormal_std_deviation_from_err_factor(mean, err_fact, std_dev);
    lognormal_params_from_moments(mean, std_dev, lambda, zeta);
  }
}


inline void moments_from_uniform_params(const Real& lwr, const Real& upr,
					Real& mean, Real& std_dev)
{ mean = (lwr + upr)/2.; std_dev = (upr - lwr)/std::sqrt(12.); }


inline void moments_from_loguniform_params(const Real& lwr, const Real& upr,
					   Real& mean, Real& std_dev)
{
  Real range = upr - lwr, log_range = std::log(upr) - std::log(lwr);
  mean       = range/log_range;
  std_dev    = std::sqrt(range*(log_range*(upr+lwr)-2.*range)/2.)/log_range;
}


inline void moments_from_triangular_params(const Real& lwr, const Real& upr,
					   const Real& mode, Real& mean,
					   Real& std_dev)
{
  mean = (lwr + mode + upr)/3.;
  std_dev
    = std::sqrt((lwr*(lwr - mode) + mode*(mode - upr) + upr*(upr - lwr))/18.);
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
  std_dev    = std::sqrt(alpha*beta/(alpha+beta+1.))/(alpha+beta)*range;
}


inline void moments_from_gamma_params(const Real& alpha, const Real& beta,
				      Real& mean, Real& std_dev)
{ mean = alpha*beta; std_dev = std::sqrt(alpha)*beta; }


inline void moments_from_gumbel_params(const Real& alpha, const Real& beta,
				       Real& mean, Real& std_dev)
{ mean = beta + 0.57721566490153286/alpha; std_dev = PI/std::sqrt(6.)/alpha; }
/* Euler-Mascheroni constant is 0.5772... */


inline void moments_from_frechet_params(const Real& alpha, const Real& beta,
					Real& mean, Real& std_dev)
{
  // See Haldar and Mahadevan, p. 91-92
  Real gam = gamma_function(1.-1./alpha);
  mean     = beta*gam;
  std_dev  = beta*std::sqrt(gamma_function(1.-2./alpha)-gam*gam);
}


inline void moments_from_weibull_params(const Real& alpha, const Real& beta,
					Real& mean, Real& std_dev)
{
  // See Haldar and Mahadevan, p. 97
  Real gam = gamma_function(1.+1./alpha),
    cf_var = std::sqrt(gamma_function(1.+2./alpha)/gam/gam - 1.);
  mean     = beta*gam;
  std_dev  = cf_var*mean;
}


inline void moments_from_histogram_bin_params(const RealVector& hist_bin_prs,
					      Real& mean, Real& std_dev)
{
  // in bin case, (x,y) and (x,c) are not equivalent since bins have non-zero
  // width -> assume (x,c) count/frequency-based (already converted from (x,y)
  // skyline/density-based) with normalization (counts sum to 1.)
  mean = std_dev = 0.;
  size_t num_bins = hist_bin_prs.length() / 2 - 1;
  for (int i=0; i<num_bins; ++i) {
    const Real& lwr   = hist_bin_prs[2*i];
    const Real& count = hist_bin_prs[2*i+1];
    const Real& upr   = hist_bin_prs[2*i+2];
    //Real clu = count * (lwr + upr);
    mean    += count * (lwr + upr); // clu
    std_dev += count * (upr*upr + upr*lwr + lwr*lwr); // upr*clu + count*lwr*lwr
  }
  mean   /= 2.;
  std_dev = std::sqrt(std_dev/3. - mean*mean);
}


inline void moments_from_poisson_params(const Real& lambda,
					Real& mean, Real& std_dev)
{
  mean    = lambda;
  std_dev = std::sqrt(lambda);
}


inline void moments_from_binomial_params(const Real& prob_pertrial,
                                         const int& num_trials, 
                                         Real& mean, Real& std_dev)
{
  mean    = prob_pertrial * num_trials;
  std_dev = std::sqrt(prob_pertrial * num_trials *(1.-prob_pertrial));
}


inline void moments_from_negative_binomial_params(const Real& prob_pertrial,
			                          const int& num_trials,
                                                  Real& mean, Real& std_dev)
{
  mean    = (Real)num_trials * (1.-prob_pertrial)/prob_pertrial;
  std_dev = std::sqrt((Real)num_trials * (1.-prob_pertrial) /
		      std::pow(prob_pertrial,2));
}


inline void moments_from_geometric_params(const Real& prob_pertrial,
					  Real& mean, Real& std_dev)
{
  mean    = (1.-prob_pertrial)/prob_pertrial;
  std_dev = std::sqrt((1.-prob_pertrial)/std::pow(prob_pertrial,2));
}


inline void moments_from_hypergeometric_params(const int& num_total_pop,
					       const int& num_sel_pop, 
					       const int& num_fail, 
					       Real& mean, Real& std_dev)
{
  mean    = (Real)(num_fail*num_sel_pop)/(Real)num_total_pop;
  std_dev = std::sqrt((Real)(num_fail*num_sel_pop*(num_total_pop-num_fail)*
			     (num_total_pop-num_sel_pop))/
		      (Real)(num_total_pop*num_total_pop*(num_total_pop-1)));
}


inline void moments_from_histogram_pt_params(const RealVector& hist_pt_prs,
					     Real& mean, Real& std_dev)
{
  // in point case, (x,y) and (x,c) are equivalent since bins have zero-width.
  // assume normalization (counts sum to 1.).
  mean = std_dev = 0.;
  size_t num_pts = hist_pt_prs.length() / 2;
  for (int i=0; i<num_pts; ++i) {
    const Real& val = hist_pt_prs[2*i];
    Real cv = hist_pt_prs[2*i+1]*val; // count * val
    mean    += cv;
    std_dev += cv*val;
  }
  std_dev = std::sqrt(std_dev - mean*mean);
}


inline Real bounded_normal_pdf(const Real& x, const Real& mean, 
			       const Real& std_dev, const Real& lwr,
			       const Real& upr)
{
  Real Phi_lms = (lwr > -DBL_MAX) ? Phi((lwr-mean)/std_dev) : 0.;
  Real Phi_ums = (upr <  DBL_MAX) ? Phi((upr-mean)/std_dev) : 1.;
  return phi((x-mean)/std_dev)/(Phi_ums - Phi_lms)/std_dev;
}


inline Real bounded_normal_cdf(const Real& x, const Real& mean, 
			       const Real& std_dev, const Real& lwr,
			       const Real& upr)
{
  Real Phi_lms = (lwr > -DBL_MAX) ? Phi((lwr-mean)/std_dev) : 0.;
  Real Phi_ums = (upr <  DBL_MAX) ? Phi((upr-mean)/std_dev) : 1.;
  return (Phi((x-mean)/std_dev) - Phi_lms)/(Phi_ums - Phi_lms);
}


inline Real lognormal_pdf(const Real& x, const Real& mean, const Real& std_dev)
{
  Real lambda, zeta;
//#ifdef HAVE_BOOST
  lognormal_params_from_moments(mean, std_dev, lambda, zeta);
  lognormal_dist logn1(lambda, zeta);
  return bmth::pdf(logn1, x);
/*
#elif HAVE_GSL
  lognormal_params_from_moments(mean, std_dev, lambda, zeta);
  return gsl_ran_lognormal_pdf(x, lambda, zeta);
#else
  Real zeta_sq;
  lognormal_zeta_sq_from_moments(mean, std_dev, zeta_sq);
  lambda = std::log(mean) - zeta_sq/2.;
  zeta = std::sqrt(zeta_sq);
  return 1./std::sqrt(2.*PI)/zeta/x *
    std::exp(-std::pow(std::log(x)-lambda, 2)/2./zeta_sq);
#endif // HAVE_GSL or HAVE_BOOST
*/
}


inline Real lognormal_cdf(const Real& x, const Real& mean, const Real& std_dev)
{
  Real lambda, zeta;
  lognormal_params_from_moments(mean, std_dev, lambda, zeta);
  return Phi((std::log(x) - lambda)/zeta);
}


inline Real bounded_lognormal_pdf(const Real& x, const Real& mean, 
				  const Real& std_dev, const Real& lwr,
				  const Real& upr)
{
  Real lambda, zeta;
  lognormal_params_from_moments(mean, std_dev, lambda, zeta);
  Real Phi_lms = (lwr > 0.)      ? Phi((std::log(lwr)-lambda)/zeta) : 0.;
  Real Phi_ums = (upr < DBL_MAX) ? Phi((std::log(upr)-lambda)/zeta) : 1.;
  return phi((std::log(x)-lambda)/zeta)/(Phi_ums-Phi_lms)/x/zeta;
}


inline Real bounded_lognormal_cdf(const Real& x, const Real& mean, 
				  const Real& std_dev, const Real& lwr,
				  const Real& upr)
{
  Real lambda, zeta;
  lognormal_params_from_moments(mean, std_dev, lambda, zeta);
  Real Phi_lms = (lwr > 0.)      ? Phi((std::log(lwr)-lambda)/zeta) : 0.;
  Real Phi_ums = (upr < DBL_MAX) ? Phi((std::log(upr)-lambda)/zeta) : 1.;
  return (Phi((std::log(x)-lambda)/zeta) - Phi_lms)/(Phi_ums - Phi_lms);
}


inline Real std_uniform_pdf()
{ return 0.5; } // equal probability on [-1,1]


inline Real std_uniform_cdf(const Real& x)
{ return (x + 1.)/2.; } // linear [-1,1] -> [0,1]


inline Real std_uniform_cdf_inverse(const Real& cdf)
{ return 2.*cdf - 1.; } // linear [0,1] -> [-1,1]


inline Real uniform_pdf(const Real& lwr, const Real& upr)
{ return 1./(upr - lwr); } // equal probability on [lwr,upr]


inline Real uniform_cdf(const Real& x, const Real& lwr, const Real& upr)
{ return (x - lwr)/(upr - lwr); } // linear [lwr,upr] -> [0,1]


inline Real loguniform_pdf(const Real& x, const Real& lwr, const Real& upr)
{ return 1./(std::log(upr) - std::log(lwr))/x; }


inline Real loguniform_cdf(const Real& x, const Real& lwr, const Real& upr)
{ return (std::log(x) - std::log(lwr))/(std::log(upr) - std::log(lwr)); }


inline Real triangular_pdf(const Real& x, const Real& mode, const Real& lwr,
			   const Real& upr)
{
  return (x < mode) ? 2.*(x-lwr)/(upr-lwr)/(mode-lwr) :
                      2.*(upr-x)/(upr-lwr)/(upr-mode);
}


inline Real triangular_cdf(const Real& x, const Real& mode, const Real& lwr,
			   const Real& upr)
{
  return (x < mode) ? std::pow(x-lwr,2.)/(upr-lwr)/(mode-lwr) :
    ((mode-lwr) - (x+mode-2*upr)*(x-mode)/(upr-mode))/(upr-lwr);
}


inline Real std_exponential_pdf(const Real& x)
{ return std::exp(-x); }


inline Real std_exponential_cdf(const Real& x)
{
  // as with log1p(), avoid numerical probs when exp(~0) is ~ 1
  return -expm1(-x); //1. - std::exp(-x);
}


inline Real exponential_pdf(const Real& x, const Real& beta)
{ return std::exp(-x/beta)/beta; }


inline Real exponential_cdf(const Real& x, const Real& beta)
{
  // as with log1p(), avoid numerical probs when exp(~0) is ~ 1
  return -expm1(-x/beta); //1. - std::exp(-x/beta);
}


inline Real std_beta_pdf(const Real& x, const Real& alpha, const Real& beta)
{
//#ifdef HAVE_BOOST
  beta_dist beta1(alpha, beta);
  return bmth::pdf(beta1, x);
/*
#elif HAVE_GSL
  // GSL beta passes alpha followed by beta
  return gsl_ran_beta_pdf(x, alpha, beta);
#endif // HAVE_GSL or HAVE_BOOST
*/
}


inline Real std_beta_cdf(const Real& x, const Real& alpha, const Real& beta)
{
//#ifdef HAVE_BOOST
  beta_dist beta1(alpha, beta);
  return bmth::cdf(beta1, x);
/*
#elif HAVE_GSL
  // GSL beta passes alpha followed by beta
  return gsl_cdf_beta_P(x, alpha, beta);
#endif //HAVE_GSL or HAVE_BOOST
*/
}


inline Real std_beta_cdf_inverse(const Real& cdf, const Real& alpha,
				 const Real& beta)
{
//#ifdef HAVE_BOOST
  beta_dist beta1(alpha, beta);
  return bmth::quantile(beta1, cdf);
/*
#elif HAVE_GSL
  // GSL does not support beta CDF inverse, therefore a special Newton solve
  // has been implemented for the inversion (which uses the GSL CDF/PDF fns).
  return cdf_beta_Pinv(cdf, alpha, beta);
#endif // HAVE_BOOST
*/
}


inline Real gamma_pdf(const Real& x, const Real& alpha, const Real& beta)
{
//#ifdef HAVE_BOOST
  gamma_dist gamma1(alpha, beta);
  return bmth::pdf(gamma1, x);
/*
#elif HAVE_GSL
  // GSL gamma passes alpha followed by beta
  return gsl_ran_gamma_pdf(x, alpha, beta);
#endif // HAVE_GSL or HAVE_BOOST
*/
}


//inline Real gamma_pdf_deriv(const Real& x, const Real& alpha,const Real& beta)
//{
//  return std::pow(beta,-alpha) / gamma_function(alpha) *
//    (std::exp(-x/beta)*(alpha-1.) * std::pow(x,alpha-2.) -
//     std::pow(x,alpha-1.) * std::exp(-x/beta)/beta);
//}


inline Real gamma_cdf(const Real& x, const Real& alpha, const Real& beta)
{
//#ifdef HAVE_BOOST
  gamma_dist gamma1(alpha, beta);
  return bmth::cdf(gamma1, x);
/*
#elif HAVE_GSL
  // GSL gamma passes alpha followed by beta
  return gsl_cdf_gamma_P(x, alpha, beta);
#endif // HAVE_GSL or HAVE_BOOST
*/
}


inline Real gamma_cdf_inverse(const Real& cdf, const Real& alpha,
			      const Real& beta)
{
//#ifdef HAVE_BOOST
  gamma_dist gamma1(alpha, beta);
  return bmth::quantile(gamma1, cdf);
/*
#elif HAVE_GSL
  // GSL gamma passes alpha followed by beta
  return gsl_cdf_gamma_Pinv(cdf, alpha, beta);
#endif // HAVE_GSL or HAVE_BOOST
*/
}


inline Real gumbel_pdf(const Real& x, const Real& alpha, const Real& beta)
{
  Real num = std::exp(-alpha*(x-beta));
  // if num overflows, apply l'Hopital's rule
  return (num > DBL_MAX) ? 0. : alpha*num*std::exp(-num);
}


inline Real gumbel_cdf(const Real& x, const Real& alpha, const Real& beta)
{ return std::exp(-std::exp(-alpha * (x - beta))); }


inline Real frechet_pdf(const Real& x, const Real& alpha, const Real& beta)
{
  Real num = std::pow(beta/x, alpha);
  return alpha/x*num*std::exp(-num);
}


inline Real frechet_cdf(const Real& x, const Real& alpha, const Real& beta)
{ return std::exp(-std::pow(beta/x, alpha)); }


inline Real weibull_pdf(const Real& x, const Real& alpha, const Real& beta)
{
//#ifdef HAVE_BOOST
  weibull_dist weibull1(alpha, beta);
  return bmth::pdf(weibull1, x);
/*
#elif HAVE_GSL
  // GSL weibull passes beta followed by alpha
  return gsl_ran_weibull_pdf(x, beta, alpha);
#else
  return alpha/beta * std::pow(x/beta,alpha-1.) *
    std::exp(-std::pow(x/beta,alpha));
#endif // HAVE_GSL or HAVE_BOOST
*/
}


inline Real weibull_cdf(const Real& x, const Real& alpha, const Real& beta)
{
//#ifdef HAVE_BOOST
  weibull_dist weibull1(alpha, beta);
  return bmth::cdf(weibull1, x);
/*
#elif HAVE_GSL
  // GSL weibull passes beta followed by alpha
  return gsl_cdf_weibull_P(x, beta, alpha);
#else
  // avoid numerical probs when exp()~1
  return -expm1(-std::pow(x/beta, alpha));
#endif // HAVE_GSL or HAVE_BOOST
*/
}


inline Real histogram_bin_pdf(const Real& x, const RealVector& hist_bin_prs)
{
  // The PDF is discontinuous at the bin bounds.  To resolve this, this
  // function employs an (arbitrary) convention: closed/inclusive lower bound
  // and open/exclusive upper bound.
  size_t num_bins = hist_bin_prs.length() / 2 - 1;
  if (x < hist_bin_prs[0] || x >= hist_bin_prs[2*num_bins])
    return 0.;
  else {
    for (int i=0; i<num_bins; ++i) {
      const Real& upr = hist_bin_prs[2*i+2];
      if (x < upr) // return density = count / (upr - lwr);
	return hist_bin_prs[2*i+1] / (upr - hist_bin_prs[2*i]);
    }
    return 0.;
  }
}


inline Real histogram_bin_cdf(const Real& x, const RealVector& hist_bin_prs)
{
  size_t num_bins = hist_bin_prs.length() / 2 - 1;
  if (x <= hist_bin_prs[0])
    return 0.;
  else if (x >= hist_bin_prs[2*num_bins])
    return 1.;
  else {
    Real cdf = 0.;
    for (int i=0; i<num_bins; ++i) {
      const Real& count = hist_bin_prs[2*i+1];
      const Real& upr   = hist_bin_prs[2*i+2];
      if (x < upr) {
	const Real& lwr = hist_bin_prs[2*i];
	cdf += count * (x - lwr) / (upr - lwr);
	break;
      }
      else
	cdf += count;
    }
    return cdf;
  }
}


inline Real histogram_bin_cdf_inverse(const Real& cdf,
				      const RealVector& hist_bin_prs)
{
  size_t num_bins = hist_bin_prs.length() / 2 - 1;
  if (cdf <= 0.)
    return hist_bin_prs[0];          // lower bound abscissa
  else if (cdf >= 1.)
    return hist_bin_prs[2*num_bins]; // upper bound abscissa
  else {
    Real upr_cdf = 0.;
    for (int i=0; i<num_bins; ++i) {
      const Real& count = hist_bin_prs[2*i+1];
      upr_cdf += count;
      if (cdf < upr_cdf) {
	const Real& upr = hist_bin_prs[2*i+2];
	const Real& lwr = hist_bin_prs[2*i];
	return lwr + (upr_cdf - cdf) / count * (upr - lwr);
      }
    }
    // If still no return due to numerical roundoff, return upper bound
    return hist_bin_prs[2*num_bins]; // upper bound abscissa
  }
}

} // namespace Pecos

#endif

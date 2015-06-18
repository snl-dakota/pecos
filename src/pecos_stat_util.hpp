/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_STAT_UTIL_HPP
#define PECOS_STAT_UTIL_HPP

#include "pecos_data_types.hpp"
#include "BoundedNormalRandomVariable.hpp"
#include "BoundedLognormalRandomVariable.hpp"
#include "LoguniformRandomVariable.hpp"
#include "TriangularRandomVariable.hpp"
#include "BetaRandomVariable.hpp"
#include "GammaRandomVariable.hpp"
#include "GumbelRandomVariable.hpp"
#include "FrechetRandomVariable.hpp"
#include "WeibullRandomVariable.hpp"
#include "HistogramBinRandomVariable.hpp"


namespace Pecos {

inline Real gamma_function(Real x)
{ return bmth::tgamma(x); }


inline void moments_from_lognormal_spec(const RealVector& ln_means,
					const RealVector& ln_std_devs,
					const RealVector& ln_lambdas,
					const RealVector& ln_zetas,
					const RealVector& ln_err_facts, 
					size_t i, Real& mean, Real& std_dev)
{
  // spec options include the mean/std dev or mean/error factor of the actual
  // lognormal distribution or the mean/std dev (lambda/zeta) of the
  // corresponding normal distribution

  if (!ln_lambdas.empty())
    LognormalRandomVariable::
      moments_from_params(ln_lambdas[i], ln_zetas[i], mean, std_dev);
  else {
    mean = ln_means[i];
    if (!ln_std_devs.empty())
      std_dev = ln_std_devs[i];
    else
      LognormalRandomVariable::
	std_deviation_from_error_factor(mean, ln_err_facts[i], std_dev);
  }
}


inline void params_from_lognormal_spec(const RealVector& ln_means,
				       const RealVector& ln_std_devs,
				       const RealVector& ln_lambdas,
				       const RealVector& ln_zetas,
				       const RealVector& ln_err_facts, 
				       size_t i, Real& lambda, Real& zeta)
{
  // spec options include the mean/std dev or mean/error factor of the actual
  // lognormal distribution or the mean/std dev (lambda/zeta) of the
  // corresponding normal distribution

  if (!ln_lambdas.empty()) {
    lambda = ln_lambdas[i];
    zeta   = ln_zetas[i];
  }
  else {
    if (!ln_std_devs.empty())
      LognormalRandomVariable::
	params_from_moments(ln_means[i], ln_std_devs[i], lambda, zeta);
    else {
      Real mean = ln_means[i], stdev;
      LognormalRandomVariable::
	std_deviation_from_error_factor(mean, ln_err_facts[i], stdev);
      LognormalRandomVariable::params_from_moments(mean, stdev, lambda, zeta);
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
  // spec options include the mean/std dev or mean/error factor of the actual
  // lognormal distribution or the mean/std dev (lambda/zeta) of the
  // corresponding normal distribution

  if (!ln_lambdas.empty()) { // lambda/zeta -> mean/std_dev
    lambda = ln_lambdas[i]; zeta = ln_zetas[i];
    LognormalRandomVariable::moments_from_params(lambda, zeta, mean, std_dev);
    LognormalRandomVariable::
      error_factor_from_std_deviation(mean, std_dev, err_fact);
  }
  else if (!ln_std_devs.empty()) {
    mean = ln_means[i]; std_dev = ln_std_devs[i];
    LognormalRandomVariable::params_from_moments(mean, std_dev, lambda, zeta);
    LognormalRandomVariable::
      error_factor_from_std_deviation(mean, std_dev, err_fact);
  }
  else { // mean/err_fact -> mean/std_dev
    mean = ln_means[i]; err_fact = ln_err_facts[i];
    LognormalRandomVariable::
      std_deviation_from_error_factor(mean, err_fact, std_dev);
    LognormalRandomVariable::params_from_moments(mean, std_dev, lambda, zeta);
  }
}


inline void moments_from_poisson_params(Real lambda, Real& mean, Real& std_dev)
{
  mean    = lambda;
  std_dev = std::sqrt(lambda);
}


inline void moments_from_binomial_params(Real prob_pertrial, int num_trials, 
                                         Real& mean, Real& std_dev)
{
  mean    = prob_pertrial * num_trials;
  std_dev = std::sqrt(prob_pertrial * num_trials *(1.-prob_pertrial));
}


inline void moments_from_negative_binomial_params(Real prob_pertrial,
						  int num_trials, Real& mean,
						  Real& std_dev)
{
  mean    = (Real)num_trials * (1.-prob_pertrial)/prob_pertrial;
  std_dev = std::sqrt((Real)num_trials * (1.-prob_pertrial) /
		      std::pow(prob_pertrial,2));
}


inline void moments_from_geometric_params(Real prob_pertrial,
					  Real& mean, Real& std_dev)
{
  mean    = (1.-prob_pertrial)/prob_pertrial;
  std_dev = std::sqrt((1.-prob_pertrial)/std::pow(prob_pertrial,2));
}


inline void moments_from_hypergeometric_params(int num_total_pop,
					       int num_sel_pop, 
					       int num_fail, 
					       Real& mean, Real& std_dev)
{
  mean    = (Real)(num_fail*num_sel_pop)/(Real)num_total_pop;
  std_dev = std::sqrt((Real)(num_fail*num_sel_pop*(num_total_pop-num_fail)*
			     (num_total_pop-num_sel_pop))/
		      (Real)(num_total_pop*num_total_pop*(num_total_pop-1)));
}


/// for integer-valued histogram, return a real-valued mean and std dev
inline void moments_from_histogram_pt_params(const IntRealMap& hist_pt_prs,
					     Real& mean, Real& std_dev)
{
  // in point case, (x,y) and (x,c) are equivalent since bins have zero-width.
  // assume normalization (counts sum to 1.).
  mean = std_dev = 0.;
  Real val, count, prod;
  IRMCIter cit = hist_pt_prs.begin();
  IRMCIter cit_end = hist_pt_prs.end();
  for ( ; cit != cit_end; ++cit) {
    val = cit->first; count = cit->second; prod = count * val;
    mean    += prod;
    std_dev += prod * val;
  }
  std_dev = std::sqrt(std_dev - mean * mean);
}


/// for string variables, define the mean as the count-weighted mean
/// of a zero-based index
inline void moments_from_histogram_pt_params(const StringRealMap& hist_pt_prs,
					     Real& mean, Real& std_dev)
{
  // in point case, (x,y) and (x,c) are equivalent since bins have zero-width.
  // assume normalization (counts sum to 1.).
  mean = std_dev = 0.;
  Real val, count, prod;
  size_t index = 0;
  SRMCIter cit = hist_pt_prs.begin();
  SRMCIter cit_end = hist_pt_prs.end();
  for ( ; cit != cit_end; ++cit, ++index) {
    val = index; count = cit->second; prod = count * val;
    mean    += prod;
    std_dev += prod * val;
  }
  std_dev = std::sqrt(std_dev - mean * mean);
}


/// return the mean and standard deviation of a real-valued point histogram
inline void moments_from_histogram_pt_params(const RealRealMap& hist_pt_prs,
					     Real& mean, Real& std_dev)
{
  // in point case, (x,y) and (x,c) are equivalent since bins have zero-width.
  // assume normalization (counts sum to 1.).
  mean = std_dev = 0.;
  Real val, count, prod;
  RRMCIter cit = hist_pt_prs.begin();
  RRMCIter cit_end = hist_pt_prs.end();
  for ( ; cit != cit_end; ++cit) {
    val = cit->first; count = cit->second; prod = count * val;
    mean    += prod;
    std_dev += prod * val;
  }
  std_dev = std::sqrt(std_dev - mean * mean);
}

} // namespace Pecos

#endif

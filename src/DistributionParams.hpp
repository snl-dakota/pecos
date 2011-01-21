/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef DISTRIBUTION_PARAMS_HPP
#define DISTRIBUTION_PARAMS_HPP

#include "pecos_data_types.hpp"


namespace Pecos {

/// The representation of a surrogate data point.  This representation,
/// or body, may be shared by multiple DistributionParams handle instances.

/** The DistributionParams/DistributionParamsRep pairs utilize a
    handle-body idiom (Coplien, Advanced C++). */

class PECOS_EXPORT DistributionParamsRep{
  //
  //- Heading: Friends
  //

  /// the handle class can access attributes of the body class directly
  friend class DistributionParams;

private:

  //
  //- Heading: Private member functions
  //

  /// default constructor
  DistributionParamsRep();
  /// constructor
  DistributionParamsRep(const RealVector& nuv_means,
    const RealVector& nuv_std_devs,  const RealVector& nuv_l_bnds,
    const RealVector& nuv_u_bnds,    const RealVector& lnuv_means,
    const RealVector& lnuv_std_devs, const RealVector& lnuv_lambdas,
    const RealVector& lnuv_zetas,    const RealVector& lnuv_err_facts,
    const RealVector& lnuv_l_bnds,   const RealVector& lnuv_u_bnds,
    const RealVector& uuv_l_bnds,    const RealVector& uuv_u_bnds,
    const RealVector& luuv_l_bnds,   const RealVector& luuv_u_bnds,
    const RealVector& tuv_modes,     const RealVector& tuv_l_bnds,
    const RealVector& tuv_u_bnds,    const RealVector& euv_betas,
    const RealVector& beuv_alphas,   const RealVector& beuv_betas,
    const RealVector& beuv_l_bnds,   const RealVector& beuv_u_bnds,
    const RealVector& gauv_alphas,   const RealVector& gauv_betas,
    const RealVector& guuv_alphas,   const RealVector& guuv_betas,
    const RealVector& fuv_alphas,    const RealVector& fuv_betas,
    const RealVector& wuv_alphas,    const RealVector& wuv_betas,
    const RealVectorArray& hbuv_prs, const RealVector& puv_lambdas,
    const RealVector& biuv_p_per_tr, const IntVector&  biuv_num_trials, 
    const RealVector& nbuv_p_per_tr, const IntVector&  nbuv_num_trials, 
    const RealVector& geuv_p_per_tr, const IntVector&  hguv_tot_pop,
    const IntVector& hguv_sel_pop,   const IntVector& hguv_num_drawn,
    const RealVectorArray& hpuv_prs, const RealVectorArray& iuv_probs,
    const RealVectorArray& iuv_bnds, const RealSymMatrix& uv_corr);
  /// destructor
  ~DistributionParamsRep();

  //
  //- Heading: Private data members
  //

  /// normal uncertain variable means
  RealVector normalMeans;
  /// normal uncertain variable standard deviations
  RealVector normalStdDevs;
  /// normal uncertain variable lower bounds
  RealVector normalLowerBnds;
  /// normal uncertain variable upper bounds
  RealVector normalUpperBnds;
  /// lognormal uncertain variable means
  RealVector lognormalMeans;
  /// lognormal uncertain variable standard deviations
  RealVector lognormalStdDevs;
  /// lognormal uncertain variable lambdas
  RealVector lognormalLambdas;
  /// lognormal uncertain variable zetas
  RealVector lognormalZetas;
  /// lognormal uncertain variable error factors
  RealVector lognormalErrFacts;
  /// lognormal uncertain variable lower bounds
  RealVector lognormalLowerBnds;
  /// lognormal uncertain variable upper bounds
  RealVector lognormalUpperBnds;
  /// uniform uncertain variable lower bounds
  RealVector uniformLowerBnds;
  /// uniform uncertain variable upper bounds
  RealVector uniformUpperBnds;
  /// loguniform uncertain variable lower bounds
  RealVector loguniformLowerBnds;
  /// loguniform uncertain variable upper bounds
  RealVector loguniformUpperBnds;
  /// triangular uncertain variable modes
  RealVector triangularModes;
  /// triangular uncertain variable lower bounds
  RealVector triangularLowerBnds;
  /// triangular uncertain variable upper bounds
  RealVector triangularUpperBnds;
  /// exponential uncertain variable betas
  RealVector exponentialBetas;
  /// beta uncertain variable alphas
  RealVector betaAlphas;
  /// beta uncertain variable betas
  RealVector betaBetas;
  /// beta uncertain variable lower bounds
  RealVector betaLowerBnds;
  /// beta uncertain variable upper bounds
  RealVector betaUpperBnds;
  /// gamma uncertain variable alphas
  RealVector gammaAlphas;
  /// gamma uncertain variable betas
  RealVector gammaBetas;
  /// gumbel uncertain variable alphas
  RealVector gumbelAlphas;
  /// gumbel uncertain variable betas
  RealVector gumbelBetas;
  /// frechet uncertain variable alphas
  RealVector frechetAlphas;
  /// frechet uncertain variable betas
  RealVector frechetBetas;
  /// weibull uncertain variable alphas
  RealVector weibullAlphas;
  /// weibull uncertain variable betas
  RealVector weibullBetas;
  /// histogram uncertain (x,y) bin pairs (continuous linear histogram)
  RealVectorArray histogramBinPairs;
  /// poisson uncertain variable lambdas
  RealVector poissonLambdas;
  /// binomial uncertain variable probabilities per trial
  RealVector binomialProbPerTrial;
  /// binomial uncertain variable numbers of trials
  IntVector binomialNumTrials;
  /// negative binomial uncertain variable probabilities per trial
  RealVector negBinomialProbPerTrial;
  /// negative binomial uncertain variable numbers of trials
  IntVector negBinomialNumTrials;
  /// geometric uncertain variable probabilities per trial
  RealVector geometricProbPerTrial;
  /// hypergeometric uncertain variable numbers in total population
  IntVector hyperGeomTotalPopulation;
  /// hypergeometric uncertain variable numbers in selected population
  IntVector hyperGeomSelectedPopulation;
  /// hypergeometric uncertain variable numbers failed in population
  IntVector hyperGeomNumDrawn;
  /// histogram uncertain (x,y) point pairs (discrete histogram)
  RealVectorArray histogramPointPairs;
  /// basic probability values for interval uncertain variables
  RealVectorArray intervalBasicProbs;
  /// interval lower/upper bounds for interval uncertain variables
  RealVectorArray intervalBounds;

  /// uncertain variable correlation matrix (rank correlations for sampling
  /// and correlation coefficients for reliability)
  RealSymMatrix uncertainCorrelations;

  /// number of handle objects sharing dpRep
  int referenceCount;
};


inline DistributionParamsRep::DistributionParamsRep() : referenceCount(1)
{ }


inline DistributionParamsRep::
DistributionParamsRep(const RealVector& nuv_means,
  const RealVector& nuv_std_devs,  const RealVector& nuv_l_bnds,
  const RealVector& nuv_u_bnds,    const RealVector& lnuv_means,
  const RealVector& lnuv_std_devs, const RealVector& lnuv_lambdas,
  const RealVector& lnuv_zetas,    const RealVector& lnuv_err_facts,
  const RealVector& lnuv_l_bnds,   const RealVector& lnuv_u_bnds,
  const RealVector& uuv_l_bnds,    const RealVector& uuv_u_bnds,
  const RealVector& luuv_l_bnds,   const RealVector& luuv_u_bnds,
  const RealVector& tuv_modes,     const RealVector& tuv_l_bnds,
  const RealVector& tuv_u_bnds,    const RealVector& euv_betas,
  const RealVector& beuv_alphas,   const RealVector& beuv_betas,
  const RealVector& beuv_l_bnds,   const RealVector& beuv_u_bnds,
  const RealVector& gauv_alphas,   const RealVector& gauv_betas,
  const RealVector& guuv_alphas,   const RealVector& guuv_betas,
  const RealVector& fuv_alphas,    const RealVector& fuv_betas,
  const RealVector& wuv_alphas,    const RealVector& wuv_betas,
  const RealVectorArray& hbuv_prs, const RealVector& puv_lambdas,
  const RealVector& biuv_p_per_tr, const IntVector&  biuv_num_trials, 
  const RealVector& nbuv_p_per_tr, const IntVector&  nbuv_num_trials, 
  const RealVector& geuv_p_per_tr, const IntVector&  hguv_tot_pop,
  const IntVector& hguv_sel_pop,   const IntVector& hguv_num_drawn,
  const RealVectorArray& hpuv_prs, const RealVectorArray& iuv_probs,
  const RealVectorArray& iuv_bnds, const RealSymMatrix& uv_corr):
  normalMeans(nuv_means), normalStdDevs(nuv_std_devs),
  normalLowerBnds(nuv_l_bnds), normalUpperBnds(nuv_u_bnds),
  lognormalMeans(lnuv_means), lognormalStdDevs(lnuv_std_devs),
  lognormalLambdas(lnuv_lambdas), lognormalZetas(lnuv_zetas),
  lognormalErrFacts(lnuv_err_facts), lognormalLowerBnds(lnuv_l_bnds),
  lognormalUpperBnds(lnuv_u_bnds), uniformLowerBnds(uuv_l_bnds),
  uniformUpperBnds(uuv_u_bnds), loguniformLowerBnds(luuv_l_bnds),
  loguniformUpperBnds(luuv_u_bnds), triangularModes(tuv_modes),
  triangularLowerBnds(tuv_l_bnds), triangularUpperBnds(tuv_u_bnds),
  exponentialBetas(euv_betas), betaAlphas(beuv_alphas), betaBetas(beuv_betas),
  betaLowerBnds(beuv_l_bnds), betaUpperBnds(beuv_u_bnds),
  gammaAlphas(gauv_alphas), gammaBetas(gauv_betas), gumbelAlphas(guuv_alphas),
  gumbelBetas(guuv_betas), frechetAlphas(fuv_alphas), frechetBetas(fuv_betas),
  weibullAlphas(wuv_alphas), weibullBetas(wuv_betas),
  histogramBinPairs(hbuv_prs), poissonLambdas(puv_lambdas),
  binomialProbPerTrial(biuv_p_per_tr), binomialNumTrials(biuv_num_trials), 
  negBinomialProbPerTrial(nbuv_p_per_tr), negBinomialNumTrials(nbuv_num_trials),
  geometricProbPerTrial(geuv_p_per_tr), hyperGeomTotalPopulation(hguv_tot_pop),
  hyperGeomSelectedPopulation(hguv_sel_pop),  hyperGeomNumDrawn(hguv_num_drawn),
  histogramPointPairs(hpuv_prs), intervalBasicProbs(iuv_probs),
  intervalBounds(iuv_bnds), uncertainCorrelations(uv_corr), referenceCount(1)
{ }


inline DistributionParamsRep::~DistributionParamsRep()
{ }


/// Container class encapsulating basic parameter and response data
/// for defining a "truth" data point.

/** A list of these data points is contained in Approximation instances
    (e.g., Dakota::Approximation::currentPoints) and provides the data
    to build the approximation.  A handle-body idiom is used to avoid
    excessive data copying overhead. */

class PECOS_EXPORT DistributionParams{
public:

  //
  //- Heading: Constructors, destructor, and operators
  //

  /// default constructor
  DistributionParams();
  /// standard constructor
  DistributionParams(const RealVector& nuv_means,
    const RealVector& nuv_std_devs,  const RealVector& nuv_l_bnds,
    const RealVector& nuv_u_bnds,    const RealVector& lnuv_means,
    const RealVector& lnuv_std_devs, const RealVector& lnuv_lambdas,
    const RealVector& lnuv_zetas,    const RealVector& lnuv_err_facts,
    const RealVector& lnuv_l_bnds,   const RealVector& lnuv_u_bnds,
    const RealVector& uuv_l_bnds,    const RealVector& uuv_u_bnds,
    const RealVector& luuv_l_bnds,   const RealVector& luuv_u_bnds,
    const RealVector& tuv_modes,     const RealVector& tuv_l_bnds,
    const RealVector& tuv_u_bnds,    const RealVector& euv_betas,
    const RealVector& beuv_alphas,   const RealVector& beuv_betas,
    const RealVector& beuv_l_bnds,   const RealVector& beuv_u_bnds,
    const RealVector& gauv_alphas,   const RealVector& gauv_betas,
    const RealVector& guuv_alphas,   const RealVector& guuv_betas,
    const RealVector& fuv_alphas,    const RealVector& fuv_betas,
    const RealVector& wuv_alphas,    const RealVector& wuv_betas,
    const RealVectorArray& hbuv_prs, const RealVector& puv_lambdas,
    const RealVector& biuv_p_per_tr, const IntVector&  biuv_num_trials, 
    const RealVector& nbuv_p_per_tr, const IntVector&  nbuv_num_trials, 
    const RealVector& geuv_p_per_tr, const IntVector&  hguv_tot_pop,
    const IntVector& hguv_sel_pop,   const IntVector& hguv_num_drawn,
    const RealVectorArray& hpuv_prs, const RealVectorArray& iuv_probs,
    const RealVectorArray& iuv_bnds, const RealSymMatrix& uv_corr);
  /// copy constructor
  DistributionParams(const DistributionParams& dp);
  /// destructor
  ~DistributionParams();

  /// assignment operator
  DistributionParams& operator=(const DistributionParams& dp);

  //
  //- Heading: member functions
  //

  /// deep copy (as opposed to operator= shallow copy)
  void copy(const DistributionParams& dp);
  /// data update (no changes to representation (unless null))
  void update(const DistributionParams& dp);

  /// return the normal uncertain variable means
  const RealVector& normal_means() const;
  /// return the ith normal uncertain variable mean
  const Real& normal_mean(size_t i) const;
  /// set the normal uncertain variable means
  void normal_means(const RealVector& n_means);
  /// set the ith normal uncertain variable mean
  void normal_mean(const Real& n_mean, size_t i);
  /// return the normal uncertain variable standard deviations
  const RealVector& normal_std_deviations() const;
  /// return the ith normal uncertain variable standard deviation
  const Real& normal_std_deviation(size_t i) const;
  /// set the normal uncertain variable standard deviations
  void normal_std_deviations(const RealVector& n_std_devs);
  /// set the ith normal uncertain variable standard deviation
  void normal_std_deviation(const Real& n_std_dev, size_t i);
  /// return the normal uncertain variable lower bounds
  const RealVector& normal_lower_bounds() const;
  /// return the ith normal uncertain variable lower bound
  const Real& normal_lower_bound(size_t i) const;
  /// set the normal uncertain variable lower bounds
  void normal_lower_bounds(const RealVector& n_lower_bnds);
  /// set the ith normal uncertain variable lower bound
  void normal_lower_bound(const Real& n_lower_bnd, size_t i);
  /// return the normal uncertain variable upper bounds
  const RealVector& normal_upper_bounds() const;
  /// return the ith normal uncertain variable upper bound
  const Real& normal_upper_bound(size_t i) const;
  /// set the normal uncertain variable upper bounds
  void normal_upper_bounds(const RealVector& n_upper_bnds);
  /// set the ith normal uncertain variable upper bound
  void normal_upper_bound(const Real& n_upper_bnd, size_t i);

  /// return the lognormal uncertain variable means
  const RealVector& lognormal_means() const;
  /// return the ith lognormal uncertain variable mean
  const Real& lognormal_mean(size_t i) const;
  /// set the lognormal uncertain variable means
  void lognormal_means(const RealVector& ln_means);
  /// set the ith lognormal uncertain variable mean
  void lognormal_mean(const Real& ln_mean, size_t i);
  /// return the lognormal uncertain variable standard deviations
  const RealVector& lognormal_std_deviations() const;
  /// return the ith lognormal uncertain variable standard deviation
  const Real& lognormal_std_deviation(size_t i) const;
  /// set the lognormal uncertain variable standard deviations
  void lognormal_std_deviations(const RealVector& ln_std_devs);
  /// set the ith lognormal uncertain variable standard deviation
  void lognormal_std_deviation(const Real& ln_std_dev, size_t i);
  /// return the lognormal uncertain variable lambdas
  const RealVector& lognormal_lambdas() const;
  /// return the ith lognormal uncertain variable lambda
  const Real& lognormal_lambda(size_t i) const;
  /// set the lognormal uncertain variable lambdas
  void lognormal_lambdas(const RealVector& ln_lambdas);
  /// set the ith lognormal uncertain variable lambda
  void lognormal_lambda(const Real& ln_lambda, size_t i);
  /// return the lognormal uncertain variable zetas
  const RealVector& lognormal_zetas() const;
  /// return the ith lognormal uncertain variable zeta
  const Real& lognormal_zeta(size_t i) const;
  /// set the lognormal uncertain variable zetas
  void lognormal_zetas(const RealVector& ln_std_devs);
  /// set the ith lognormal uncertain variable zeta
  void lognormal_zeta(const Real& ln_std_dev, size_t i);
  /// return the lognormal uncertain variable error factors
  const RealVector& lognormal_error_factors() const;
  /// return the ith lognormal uncertain variable error factor
  const Real& lognormal_error_factor(size_t i) const;
  /// set the lognormal uncertain variable error factors
  void lognormal_error_factors(const RealVector& ln_err_facts);
  /// set the ith lognormal uncertain variable error factor
  void lognormal_error_factor(const Real& ln_err_fact, size_t i);
  /// return the lognormal uncertain variable lower bounds
  const RealVector& lognormal_lower_bounds() const;
  /// return the ith lognormal uncertain variable lower bound
  const Real& lognormal_lower_bound(size_t i) const;
  /// set the lognormal uncertain variable lower bounds
  void lognormal_lower_bounds(const RealVector& ln_lower_bnds);
  /// set the ith lognormal uncertain variable lower bound
  void lognormal_lower_bound(const Real& ln_lower_bnd, size_t i);
  /// return the lognormal uncertain variable upper bounds
  const RealVector& lognormal_upper_bounds() const;
  /// return the ith lognormal uncertain variable upper bound
  const Real& lognormal_upper_bound(size_t i) const;
  /// set the lognormal uncertain variable upper bounds
  void lognormal_upper_bounds(const RealVector& ln_upper_bnds);
  /// set the ith lognormal uncertain variable upper bound
  void lognormal_upper_bound(const Real& ln_upper_bnd, size_t i);

  /// return the uniform uncertain variable lower bounds
  const RealVector& uniform_lower_bounds() const;
  /// return the ith uniform uncertain variable lower bound
  const Real& uniform_lower_bound(size_t i) const;
  /// set the uniform uncertain variable lower bounds
  void uniform_lower_bounds(const RealVector& u_lower_bnds);
  /// set the ith uniform uncertain variable lower bound
  void uniform_lower_bound(const Real& u_lower_bnd, size_t i);
  /// return the uniform uncertain variable upper bounds
  const RealVector& uniform_upper_bounds() const;
  /// return the ith uniform uncertain variable upper bound
  const Real& uniform_upper_bound(size_t i) const;
  /// set the uniform uncertain variable upper bounds
  void uniform_upper_bounds(const RealVector& u_upper_bnds);
  /// set the ith uniform uncertain variable upper bound
  void uniform_upper_bound(const Real& u_upper_bnd, size_t i);

  /// return the loguniform uncertain variable lower bounds
  const RealVector& loguniform_lower_bounds() const;
  /// return the ith loguniform uncertain variable lower bound
  const Real& loguniform_lower_bound(size_t i) const;
  /// set the loguniform uncertain variable lower bounds
  void loguniform_lower_bounds(const RealVector& lu_lower_bnds);
  /// set the ith loguniform uncertain variable lower bound
  void loguniform_lower_bound(const Real& lu_lower_bnd, size_t i);
  /// return the loguniform uncertain variable upper bounds
  const RealVector& loguniform_upper_bounds() const;
  /// return the ith loguniform uncertain variable upper bound
  const Real& loguniform_upper_bound(size_t i) const;
  /// set the loguniform uncertain variable upper bounds
  void loguniform_upper_bounds(const RealVector& lu_upper_bnds);
  /// set the ith loguniform uncertain variable upper bound
  void loguniform_upper_bound(const Real& lu_upper_bnd, size_t i);

  /// return the triangular uncertain variable modes
  const RealVector& triangular_modes() const;
  /// return the ith triangular uncertain variable mode
  const Real& triangular_mode(size_t i) const;
  /// set the triangular uncertain variable modes
  void triangular_modes(const RealVector& t_modes);
  /// set the ith triangular uncertain variable mode
  void triangular_mode(const Real& t_mode, size_t i);
  /// return the triangular uncertain variable lower bounds
  const RealVector& triangular_lower_bounds() const;
  /// return the ith triangular uncertain variable lower bound
  const Real& triangular_lower_bound(size_t i) const;
  /// set the triangular uncertain variable lower bounds
  void triangular_lower_bounds(const RealVector& t_lower_bnds);
  /// set the ith triangular uncertain variable lower bound
  void triangular_lower_bound(const Real& t_lower_bnd, size_t i);
  /// return the triangular uncertain variable upper bounds
  const RealVector& triangular_upper_bounds() const;
  /// return the ith triangular uncertain variable upper bound
  const Real& triangular_upper_bound(size_t i) const;
  /// set the triangular uncertain variable upper bounds
  void triangular_upper_bounds(const RealVector& t_upper_bnds);
  /// set the ith triangular uncertain variable upper bound
  void triangular_upper_bound(const Real& t_upper_bnd, size_t i);

  /// return the exponential uncertain variable beta parameters
  const RealVector& exponential_betas() const;
  /// return the ith exponential uncertain variable beta parameter
  const Real& exponential_beta(size_t i) const;
  /// set the exponential uncertain variable beta parameters
  void exponential_betas(const RealVector& e_betas);
  /// set the ith exponential uncertain variable beta parameter
  void exponential_beta(const Real& e_beta, size_t i);

  /// return the beta uncertain variable alphas
  const RealVector& beta_alphas() const;
  /// return the ith beta uncertain variable alpha
  const Real& beta_alpha(size_t i) const;
  /// set the beta uncertain variable alphas
  void beta_alphas(const RealVector& b_alphas);
  /// set the ith beta uncertain variable alpha
  void beta_alpha(const Real& b_alpha, size_t i);
  /// return the beta uncertain variable betas
  const RealVector& beta_betas() const;
  /// return the ith beta uncertain variable beta
  const Real& beta_beta(size_t i) const;
  /// set the beta uncertain variable betas
  void beta_betas(const RealVector& b_betas);
  /// set the ith beta uncertain variable beta
  void beta_beta(const Real& b_beta, size_t i);
  /// return the beta uncertain variable lower bounds
  const RealVector& beta_lower_bounds() const;
  /// return the ith beta uncertain variable lower bound
  const Real& beta_lower_bound(size_t i) const;
  /// set the beta uncertain variable lower bounds
  void beta_lower_bounds(const RealVector& b_lower_bnds);
  /// set the ith beta uncertain variable lower bound
  void beta_lower_bound(const Real& b_lower_bnd, size_t i);
  /// return the beta uncertain variable upper bounds
  const RealVector& beta_upper_bounds() const;
  /// return the ith beta uncertain variable upper bound
  const Real& beta_upper_bound(size_t i) const;
  /// set the beta uncertain variable upper bounds
  void beta_upper_bounds(const RealVector& b_upper_bnds);
  /// set the ith beta uncertain variable upper bound
  void beta_upper_bound(const Real& b_upper_bnd, size_t i);

  /// return the gamma uncertain variable alpha parameters
  const RealVector& gamma_alphas() const;
  /// return the ith gamma uncertain variable alpha parameter
  const Real& gamma_alpha(size_t i) const;
  /// set the gamma uncertain variable alpha parameters
  void gamma_alphas(const RealVector& ga_alphas);
  /// set the ith gamma uncertain variable alpha parameter
  void gamma_alpha(const Real& ga_alpha, size_t i);
  /// return the gamma uncertain variable beta parameters
  const RealVector& gamma_betas() const;
  /// return the ith gamma uncertain variable beta parameter
  const Real& gamma_beta(size_t i) const;
  /// set the gamma uncertain variable beta parameters
  void gamma_betas(const RealVector& ga_betas);
  /// set the ith gamma uncertain variable beta parameter
  void gamma_beta(const Real& ga_beta, size_t i);

  /// return the gumbel uncertain variable alphas
  const RealVector& gumbel_alphas() const;
  /// return the ith gumbel uncertain variable alpha
  const Real& gumbel_alpha(size_t i) const;
  /// set the gumbel uncertain variable alphas
  void gumbel_alphas(const RealVector& gu_alphas);
  /// set the ith gumbel uncertain variable alpha
  void gumbel_alpha(const Real& gu_alpha, size_t i);
  /// return the gumbel uncertain variable betas
  const RealVector& gumbel_betas() const;
  /// return the ith gumbel uncertain variable beta
  const Real& gumbel_beta(size_t i) const;
  /// set the gumbel uncertain variable betas
  void gumbel_betas(const RealVector& gu_betas);
  /// set the ith gumbel uncertain variable beta
  void gumbel_beta(const Real& gu_beta, size_t i);

  /// return the frechet uncertain variable alpha parameters
  const RealVector& frechet_alphas() const;
  /// return the ith frechet uncertain variable alpha parameter
  const Real& frechet_alpha(size_t i) const;
  /// set the frechet uncertain variable alpha parameters
  void frechet_alphas(const RealVector& f_alphas);
  /// set the ith frechet uncertain variable alpha parameter
  void frechet_alpha(const Real& f_alpha, size_t i);
  /// return the frechet uncertain variable beta parameters
  const RealVector& frechet_betas() const;
  /// return the ith frechet uncertain variable beta parameter
  const Real& frechet_beta(size_t i) const;
  /// set the frechet uncertain variable beta parameters
  void frechet_betas(const RealVector& f_betas);
  /// set the ith frechet uncertain variable beta parameter
  void frechet_beta(const Real& f_beta, size_t i);

  /// return the weibull uncertain variable alpha parameters
  const RealVector& weibull_alphas() const;
  /// return the ith weibull uncertain variable alpha parameter
  const Real& weibull_alpha(size_t i) const;
  /// set the weibull uncertain variable alpha parameters
  void weibull_alphas(const RealVector& w_alphas);
  /// set the ith weibull uncertain variable alpha parameter
  void weibull_alpha(const Real& w_alpha, size_t i);
  /// return the weibull uncertain variable beta parameters
  const RealVector& weibull_betas() const;
  /// return the ith weibull uncertain variable beta parameter
  const Real& weibull_beta(size_t i) const;
  /// set the weibull uncertain variable beta parameters
  void weibull_betas(const RealVector& w_betas);
  /// set the ith weibull uncertain variable beta parameter
  void weibull_beta(const Real& w_beta, size_t i);

  /// return the histogram uncertain bin pairs
  const RealVectorArray& histogram_bin_pairs() const;
  /// return the ith histogram uncertain bin pair
  const RealVector& histogram_bin_pairs(size_t i) const;
  /// set the histogram uncertain bin pairs
  void histogram_bin_pairs(const RealVectorArray& h_bin_pairs);
  /// set the ith histogram uncertain bin pair
  void histogram_bin_pairs(const RealVector& h_bin_pairs_i, size_t i);

  /// return the poisson uncertain variable lambda parameters
  const RealVector& poisson_lambdas() const;
  /// return the ith poisson uncertain variable lambda parameter
  const Real& poisson_lambda(size_t i) const;
  /// set the poisson uncertain variable lambda parameters
  void poisson_lambdas(const RealVector& p_lambdas);
  /// set the ith poisson uncertain variable lambda parameter
  void poisson_lambda(const Real& p_lambda, size_t i);

  /// return the binomial probabilities per each trial (p) 
  const RealVector& binomial_probabilities_per_trial() const;
  /// return the ith binomial probabilities per each trial (p) 
  const Real& binomial_probabilities_per_trial(size_t i) const;
  /// set the binomial probabilities per each trial (p) 
  void binomial_probabilities_per_trial(const RealVector& probs_per_trial);
  /// set the ith binomial probabilities per each trial (p) 
  void binomial_probability_per_trial(const Real& prob_per_trial, size_t i);
  /// return the binomial number of trials (N)
  const IntVector& binomial_num_trials() const;
  /// return the ith binomial number of trials (N)
  int binomial_num_trials(size_t i) const;
  /// set the binomial number of trials (N)
  void binomial_num_trials(const IntVector& num_trials);
  /// set the ith binomial number of trials (N)
  void binomial_num_trials(int num_trials, size_t i);

  /// return the negative binomial probabilities per each trial (p) 
  const RealVector& negative_binomial_probabilities_per_trial() const;
  /// return the ith negative binomial probabilities per each trial (p) 
  const Real& negative_binomial_probabilities_per_trial(size_t i) const;
  /// set the negative binomial probabilities per each trial (p) 
  void negative_binomial_probabilities_per_trial(
    const RealVector& probs_per_trial);
  /// set the ith negative binomial probabilities per each trial (p) 
  void negative_binomial_probability_per_trial(
    const Real& prob_per_trial, size_t i);
  /// return the negative binomial number of trials (N)
  const IntVector& negative_binomial_num_trials() const;
  /// return the ith negative binomial number of trials (N)
  int negative_binomial_num_trials(size_t i) const;
  /// set the negative binomial number of trials (N)
  void negative_binomial_num_trials(const IntVector& num_trials);
  /// set the ith negative binomial number of trials (N)
  void negative_binomial_num_trials(int num_trials, size_t i);

  /// return the geometric probabilities per each trial (p) 
  const RealVector& geometric_probabilities_per_trial() const;
  /// return the ith geometric probabilities per each trial (p) 
  const Real& geometric_probabilities_per_trial(size_t i) const;
  /// set the geometric probabilities per each trial (p) 
  void geometric_probabilities_per_trial(const RealVector& probs_per_trial);
  /// set the ith geometric probabilities per each trial (p) 
  void geometric_probability_per_trial(const Real& prob_per_trial, size_t i);

  /// return the hypergeometric number in total population 
  const IntVector& hypergeometric_total_population() const;
  /// return the ith hypergeometric number in total population 
  int hypergeometric_total_population(size_t i) const;
  /// set the hypergeometric number in total population
  void hypergeometric_total_population(const IntVector& total_pop);
  /// set the ith hypergeometric number in total population
  void hypergeometric_total_population(int total_pop, size_t i);
  /// return the hypergeometric number in selected population
  const IntVector& hypergeometric_selected_population() const;
  /// return the ith hypergeometric number in selected population
  int hypergeometric_selected_population(size_t i) const;
  /// set the hypergeometric number in selected population
  void hypergeometric_selected_population(const IntVector& sel_pop);
  /// set the ith hypergeometric number in selected population
  void hypergeometric_selected_population(int sel_pop, size_t i);
  /// return the hypergeometric number failed
  const IntVector& hypergeometric_num_drawn() const;
  /// return the ith hypergeometric number failed
  int hypergeometric_num_drawn(size_t i) const;
  /// set the hypergeometric number in total population
  void hypergeometric_num_drawn(const IntVector& num_drawn);
  /// set the ith hypergeometric number in total population
  void hypergeometric_num_drawn(int num_drawn, size_t i);

  /// return the histogram uncertain point pairs
  const RealVectorArray& histogram_point_pairs() const;
  /// return the ith histogram uncertain point pair
  const RealVector& histogram_point_pairs(size_t i) const;
  /// set the histogram uncertain point pairs
  void histogram_point_pairs(const RealVectorArray& h_pt_pairs);
  /// set the ith histogram uncertain point pair
  void histogram_point_pairs(const RealVector& h_pt_pairs_i, size_t i);

  /// return the interval basic probability values
  const RealVectorArray& interval_probabilities() const;
  /// return the ith interval basic probability value
  const RealVector& interval_probabilities(size_t i) const;
  /// set the interval basic probability values
  void interval_probabilities(const RealVectorArray& int_probs);
  /// set the ith interval basic probability value
  void interval_probabilities(const RealVector& int_probs_i, size_t i);
  /// return the interval bounds
  const RealVectorArray& interval_bounds() const;
  /// return the ith interval bound
  const RealVector& interval_bounds(size_t i) const;
  /// set the interval bounds
  void interval_bounds(const RealVectorArray& int_bounds);
  /// set the ith interval bound
  void interval_bounds(const RealVector& int_bounds_i, size_t i);

  /// return the uncertain variable correlations
  const RealSymMatrix& uncertain_correlations() const;
  /// set the uncertain variable correlations
  void uncertain_correlations(const RealSymMatrix& uncertain_corr);

  /// function to check dpRep (does this handle contain a body)
  bool is_null() const;

private:

  //
  //- Heading: Private data members
  //

  /// pointer to the body (handle-body idiom)
  DistributionParamsRep* dpRep;
};


inline DistributionParams::DistributionParams():
  dpRep(new DistributionParamsRep())
{ }


inline DistributionParams::
DistributionParams(const RealVector& nuv_means,
  const RealVector& nuv_std_devs,  const RealVector& nuv_l_bnds,
  const RealVector& nuv_u_bnds,    const RealVector& lnuv_means,
  const RealVector& lnuv_std_devs, const RealVector& lnuv_lambdas,
  const RealVector& lnuv_zetas,    const RealVector& lnuv_err_facts,
  const RealVector& lnuv_l_bnds,   const RealVector& lnuv_u_bnds,
  const RealVector& uuv_l_bnds,    const RealVector& uuv_u_bnds,
  const RealVector& luuv_l_bnds,   const RealVector& luuv_u_bnds,
  const RealVector& tuv_modes,     const RealVector& tuv_l_bnds,
  const RealVector& tuv_u_bnds,    const RealVector& euv_betas,
  const RealVector& beuv_alphas,   const RealVector& beuv_betas,
  const RealVector& beuv_l_bnds,   const RealVector& beuv_u_bnds,
  const RealVector& gauv_alphas,   const RealVector& gauv_betas,
  const RealVector& guuv_alphas,   const RealVector& guuv_betas,
  const RealVector& fuv_alphas,    const RealVector& fuv_betas,
  const RealVector& wuv_alphas,    const RealVector& wuv_betas,
  const RealVectorArray& hbuv_prs, const RealVector& puv_lambdas,
  const RealVector& biuv_p_per_tr, const IntVector&  biuv_num_trials, 
  const RealVector& nbuv_p_per_tr, const IntVector&  nbuv_num_trials, 
  const RealVector& geuv_p_per_tr, const IntVector&  hguv_tot_pop,
  const IntVector& hguv_sel_pop,   const IntVector& hguv_num_drawn,
  const RealVectorArray& hpuv_prs, const RealVectorArray& iuv_probs,
  const RealVectorArray& iuv_bnds, const RealSymMatrix& uv_corr):
  dpRep(new DistributionParamsRep(nuv_means, nuv_std_devs, nuv_l_bnds,
	nuv_u_bnds, lnuv_means, lnuv_std_devs, lnuv_lambdas, lnuv_zetas,
	lnuv_err_facts, lnuv_l_bnds, lnuv_u_bnds, uuv_l_bnds, uuv_u_bnds,
	luuv_l_bnds, luuv_u_bnds, tuv_modes, tuv_l_bnds, tuv_u_bnds, euv_betas,
	beuv_alphas, beuv_betas, beuv_l_bnds, beuv_u_bnds, gauv_alphas,
	gauv_betas, guuv_alphas, guuv_betas, fuv_alphas, fuv_betas, wuv_alphas,
	wuv_betas, hbuv_prs, puv_lambdas, biuv_p_per_tr, biuv_num_trials,
	nbuv_p_per_tr, nbuv_num_trials, geuv_p_per_tr, hguv_tot_pop,
	hguv_sel_pop, hguv_num_drawn, hpuv_prs, iuv_probs, iuv_bnds, uv_corr))
{ }


inline DistributionParams::DistributionParams(const DistributionParams& dp)
{
  // Increment new (no old to decrement)
  dpRep = dp.dpRep;
  if (dpRep) // Check for an assignment of NULL
    dpRep->referenceCount++;
}


inline DistributionParams::~DistributionParams()
{
  if (dpRep) { // Check for NULL
    --dpRep->referenceCount; // decrement
    if (dpRep->referenceCount == 0)
      delete dpRep;
  }
}


inline DistributionParams& DistributionParams::
operator=(const DistributionParams& dp)
{
  // Decrement old
  if (dpRep) // Check for NULL
    if ( --dpRep->referenceCount == 0 ) 
      delete dpRep;
  // Increment new
  dpRep = dp.dpRep;
  if (dpRep) // Check for an assignment of NULL
    dpRep->referenceCount++;
  return *this;
}


inline void DistributionParams::copy(const DistributionParams& dp)
{ 
  // Decrement old
  if (dpRep) // Check for NULL
    if ( --dpRep->referenceCount == 0 ) 
      delete dpRep;
  // Create new
  dpRep = new DistributionParamsRep(dp.normal_means(),
    dp.normal_std_deviations(), dp.normal_lower_bounds(),
    dp.normal_upper_bounds(), dp.lognormal_means(),
    dp.lognormal_std_deviations(), dp.lognormal_lambdas(), dp.lognormal_zetas(),
    dp.lognormal_error_factors(), dp.lognormal_lower_bounds(),
    dp.lognormal_upper_bounds(), dp.uniform_lower_bounds(),
    dp.uniform_upper_bounds(), dp.loguniform_lower_bounds(),
    dp.loguniform_upper_bounds(), dp.triangular_modes(),
    dp.triangular_lower_bounds(), dp.triangular_upper_bounds(),
    dp.exponential_betas(), dp.beta_alphas(), dp.beta_betas(),
    dp.beta_lower_bounds(), dp.beta_upper_bounds(), dp.gamma_alphas(),
    dp.gamma_betas(), dp.gumbel_alphas(), dp.gumbel_betas(),
    dp.frechet_alphas(), dp.frechet_betas(), dp.weibull_alphas(),
    dp.weibull_betas(), dp.histogram_bin_pairs(), dp.poisson_lambdas(),
    dp.binomial_probabilities_per_trial(), dp.binomial_num_trials(),
    dp.negative_binomial_probabilities_per_trial(),
    dp.negative_binomial_num_trials(), dp.geometric_probabilities_per_trial(),
    dp.hypergeometric_total_population(),
    dp.hypergeometric_selected_population(), dp.hypergeometric_num_drawn(),
    dp.histogram_point_pairs(), dp.interval_probabilities(),
    dp.interval_bounds(), dp.uncertain_correlations());
}


inline const RealVector& DistributionParams::normal_means() const
{ return dpRep->normalMeans; }


inline const Real& DistributionParams::normal_mean(size_t i) const
{ return dpRep->normalMeans[i]; }


inline void DistributionParams::normal_means(const RealVector& n_means)
{ dpRep->normalMeans = n_means; }


inline void DistributionParams::normal_mean(const Real& n_mean, size_t i)
{ dpRep->normalMeans[i] = n_mean; }


inline const RealVector& DistributionParams::normal_std_deviations() const
{ return dpRep->normalStdDevs; }


inline const Real& DistributionParams::normal_std_deviation(size_t i) const
{ return dpRep->normalStdDevs[i]; }


inline void DistributionParams::
normal_std_deviations(const RealVector& n_std_devs)
{ dpRep->normalStdDevs = n_std_devs; }


inline void DistributionParams::
normal_std_deviation(const Real& n_std_dev, size_t i)
{ dpRep->normalStdDevs[i] = n_std_dev; }


inline const RealVector& DistributionParams::normal_lower_bounds() const
{ return dpRep->normalLowerBnds; }


inline const Real& DistributionParams::normal_lower_bound(size_t i) const
{ return dpRep->normalLowerBnds[i]; }


inline void DistributionParams::
normal_lower_bounds(const RealVector& n_lower_bnds)
{ dpRep->normalLowerBnds = n_lower_bnds; }


inline void DistributionParams::
normal_lower_bound(const Real& n_lower_bnd, size_t i)
{ dpRep->normalLowerBnds[i] = n_lower_bnd; }


inline const RealVector& DistributionParams::normal_upper_bounds() const
{ return dpRep->normalUpperBnds; }


inline const Real& DistributionParams::normal_upper_bound(size_t i) const
{ return dpRep->normalUpperBnds[i]; }


inline void DistributionParams::
normal_upper_bounds(const RealVector& n_upper_bnds)
{ dpRep->normalUpperBnds = n_upper_bnds; }


inline void DistributionParams::
normal_upper_bound(const Real& n_upper_bnd, size_t i)
{ dpRep->normalUpperBnds[i] = n_upper_bnd; }


inline const RealVector& DistributionParams::lognormal_means() const
{ return dpRep->lognormalMeans; }


inline const Real& DistributionParams::lognormal_mean(size_t i) const
{ return dpRep->lognormalMeans[i]; }


inline void DistributionParams::lognormal_means(const RealVector& ln_means)
{ dpRep->lognormalMeans = ln_means; }


inline void DistributionParams::lognormal_mean(const Real& ln_mean, size_t i)
{ dpRep->lognormalMeans[i] = ln_mean; }


inline const RealVector& DistributionParams::lognormal_std_deviations() const
{ return dpRep->lognormalStdDevs; }


inline const Real& DistributionParams::lognormal_std_deviation(size_t i) const
{ return dpRep->lognormalStdDevs[i]; }


inline void DistributionParams::
lognormal_std_deviations(const RealVector& ln_std_devs)
{ dpRep->lognormalStdDevs = ln_std_devs; }


inline void DistributionParams::
lognormal_std_deviation(const Real& ln_std_dev, size_t i)
{ dpRep->lognormalStdDevs[i] = ln_std_dev; }


inline const RealVector& DistributionParams::lognormal_lambdas() const
{ return dpRep->lognormalLambdas; }


inline const Real& DistributionParams::lognormal_lambda(size_t i) const
{ return dpRep->lognormalLambdas[i]; }


inline void DistributionParams::lognormal_lambdas(const RealVector& ln_lambdas)
{ dpRep->lognormalLambdas = ln_lambdas; }


inline void DistributionParams::
lognormal_lambda(const Real& ln_lambda, size_t i)
{ dpRep->lognormalLambdas[i] = ln_lambda; }


inline const RealVector& DistributionParams::lognormal_zetas() const
{ return dpRep->lognormalZetas; }


inline const Real& DistributionParams::lognormal_zeta(size_t i) const
{ return dpRep->lognormalZetas[i]; }


inline void DistributionParams::lognormal_zetas(const RealVector& ln_zetas)
{ dpRep->lognormalZetas = ln_zetas; }


inline void DistributionParams::lognormal_zeta(const Real& ln_zeta, size_t i)
{ dpRep->lognormalZetas[i] = ln_zeta; }


inline const RealVector& DistributionParams::lognormal_error_factors() const
{ return dpRep->lognormalErrFacts; }


inline const Real& DistributionParams::lognormal_error_factor(size_t i) const
{ return dpRep->lognormalErrFacts[i]; }


inline void DistributionParams::
lognormal_error_factors(const RealVector& ln_err_facts)
{ dpRep->lognormalErrFacts = ln_err_facts; }


inline void DistributionParams::
lognormal_error_factor(const Real& ln_err_fact, size_t i)
{ dpRep->lognormalErrFacts[i] = ln_err_fact; }


inline const RealVector& DistributionParams::lognormal_lower_bounds() const
{ return dpRep->lognormalLowerBnds; }


inline const Real& DistributionParams::lognormal_lower_bound(size_t i) const
{ return dpRep->lognormalLowerBnds[i]; }


inline void DistributionParams::
lognormal_lower_bounds(const RealVector& ln_lower_bnds)
{ dpRep->lognormalLowerBnds = ln_lower_bnds; }


inline void DistributionParams::
lognormal_lower_bound(const Real& ln_lower_bnd, size_t i)
{ dpRep->lognormalLowerBnds[i] = ln_lower_bnd; }


inline const RealVector& DistributionParams::lognormal_upper_bounds() const
{ return dpRep->lognormalUpperBnds; }


inline const Real& DistributionParams::lognormal_upper_bound(size_t i) const
{ return dpRep->lognormalUpperBnds[i]; }


inline void DistributionParams::
lognormal_upper_bounds(const RealVector& ln_upper_bnds)
{ dpRep->lognormalUpperBnds = ln_upper_bnds; }


inline void DistributionParams::
lognormal_upper_bound(const Real& ln_upper_bnd, size_t i)
{ dpRep->lognormalUpperBnds[i] = ln_upper_bnd; }


inline const RealVector& DistributionParams::uniform_lower_bounds() const
{ return dpRep->uniformLowerBnds; }


inline const Real& DistributionParams::uniform_lower_bound(size_t i) const
{ return dpRep->uniformLowerBnds[i]; }


inline void DistributionParams::
uniform_lower_bounds(const RealVector& u_lower_bnds)
{ dpRep->uniformLowerBnds = u_lower_bnds; }


inline void DistributionParams::
uniform_lower_bound(const Real& u_lower_bnd, size_t i)
{ dpRep->uniformLowerBnds[i] = u_lower_bnd; }


inline const RealVector& DistributionParams::uniform_upper_bounds() const
{ return dpRep->uniformUpperBnds; }


inline const Real& DistributionParams::uniform_upper_bound(size_t i) const
{ return dpRep->uniformUpperBnds[i]; }


inline void DistributionParams::
uniform_upper_bounds(const RealVector& u_upper_bnds)
{ dpRep->uniformUpperBnds = u_upper_bnds; }


inline void DistributionParams::
uniform_upper_bound(const Real& u_upper_bnd, size_t i)
{ dpRep->uniformUpperBnds[i] = u_upper_bnd; }


inline const RealVector& DistributionParams::loguniform_lower_bounds() const
{ return dpRep->loguniformLowerBnds; }


inline const Real& DistributionParams::loguniform_lower_bound(size_t i) const
{ return dpRep->loguniformLowerBnds[i]; }


inline void DistributionParams::
loguniform_lower_bounds(const RealVector& lu_lower_bnds)
{ dpRep->loguniformLowerBnds = lu_lower_bnds; }


inline void DistributionParams::
loguniform_lower_bound(const Real& lu_lower_bnd, size_t i)
{ dpRep->loguniformLowerBnds[i] = lu_lower_bnd; }


inline const RealVector& DistributionParams::loguniform_upper_bounds() const
{ return dpRep->loguniformUpperBnds; }


inline const Real& DistributionParams::loguniform_upper_bound(size_t i) const
{ return dpRep->loguniformUpperBnds[i]; }


inline void DistributionParams::
loguniform_upper_bounds(const RealVector& lu_upper_bnds)
{ dpRep->loguniformUpperBnds = lu_upper_bnds; }


inline void DistributionParams::
loguniform_upper_bound(const Real& lu_upper_bnd, size_t i)
{ dpRep->loguniformUpperBnds[i] = lu_upper_bnd; }


inline const RealVector& DistributionParams::triangular_modes() const
{ return dpRep->triangularModes; }


inline const Real& DistributionParams::triangular_mode(size_t i) const
{ return dpRep->triangularModes[i]; }


inline void DistributionParams::triangular_modes(const RealVector& t_modes)
{ dpRep->triangularModes = t_modes; }


inline void DistributionParams::triangular_mode(const Real& t_mode, size_t i)
{ dpRep->triangularModes[i] = t_mode; }


inline const RealVector& DistributionParams::triangular_lower_bounds() const
{ return dpRep->triangularLowerBnds; }


inline const Real& DistributionParams::triangular_lower_bound(size_t i) const
{ return dpRep->triangularLowerBnds[i]; }


inline void DistributionParams::
triangular_lower_bounds(const RealVector& t_lower_bnds)
{ dpRep->triangularLowerBnds = t_lower_bnds; }


inline void DistributionParams::
triangular_lower_bound(const Real& t_lower_bnd, size_t i)
{ dpRep->triangularLowerBnds[i] = t_lower_bnd; }


inline const RealVector& DistributionParams::triangular_upper_bounds() const
{ return dpRep->triangularUpperBnds; }


inline const Real& DistributionParams::triangular_upper_bound(size_t i) const
{ return dpRep->triangularUpperBnds[i]; }


inline void DistributionParams::
triangular_upper_bounds(const RealVector& t_upper_bnds)
{ dpRep->triangularUpperBnds = t_upper_bnds; }


inline void DistributionParams::
triangular_upper_bound(const Real& t_upper_bnd, size_t i)
{ dpRep->triangularUpperBnds[i] = t_upper_bnd; }


inline const RealVector& DistributionParams::exponential_betas() const
{ return dpRep->exponentialBetas; }


inline const Real& DistributionParams::exponential_beta(size_t i) const
{ return dpRep->exponentialBetas[i]; }


inline void DistributionParams::exponential_betas(const RealVector& e_betas)
{ dpRep->exponentialBetas = e_betas; }


inline void DistributionParams::exponential_beta(const Real& e_beta, size_t i)
{ dpRep->exponentialBetas[i] = e_beta; }


inline const RealVector& DistributionParams::beta_alphas() const
{ return dpRep->betaAlphas; }


inline const Real& DistributionParams::beta_alpha(size_t i) const
{ return dpRep->betaAlphas[i]; }


inline void DistributionParams::beta_alphas(const RealVector& b_alphas)
{ dpRep->betaAlphas = b_alphas; }


inline void DistributionParams::beta_alpha(const Real& b_alpha, size_t i)
{ dpRep->betaAlphas[i] = b_alpha; }


inline const RealVector& DistributionParams::beta_betas() const
{ return dpRep->betaBetas; }


inline const Real& DistributionParams::beta_beta(size_t i) const
{ return dpRep->betaBetas[i]; }


inline void DistributionParams::beta_betas(const RealVector& b_betas)
{ dpRep->betaBetas = b_betas; }


inline void DistributionParams::beta_beta(const Real& b_beta, size_t i)
{ dpRep->betaBetas[i] = b_beta; }


inline const RealVector& DistributionParams::beta_lower_bounds() const
{ return dpRep->betaLowerBnds; }


inline const Real& DistributionParams::beta_lower_bound(size_t i) const
{ return dpRep->betaLowerBnds[i]; }


inline void DistributionParams::
beta_lower_bounds(const RealVector& b_lower_bnds)
{ dpRep->betaLowerBnds = b_lower_bnds; }


inline void DistributionParams::
beta_lower_bound(const Real& b_lower_bnd, size_t i)
{ dpRep->betaLowerBnds[i] = b_lower_bnd; }


inline const RealVector& DistributionParams::beta_upper_bounds() const
{ return dpRep->betaUpperBnds; }


inline const Real& DistributionParams::beta_upper_bound(size_t i) const
{ return dpRep->betaUpperBnds[i]; }


inline void DistributionParams::
beta_upper_bounds(const RealVector& b_upper_bnds)
{ dpRep->betaUpperBnds = b_upper_bnds; }


inline void DistributionParams::
beta_upper_bound(const Real& b_upper_bnd, size_t i)
{ dpRep->betaUpperBnds[i] = b_upper_bnd; }


inline const RealVector& DistributionParams::gamma_alphas() const
{ return dpRep->gammaAlphas; }


inline const Real& DistributionParams::gamma_alpha(size_t i) const
{ return dpRep->gammaAlphas[i]; }


inline void DistributionParams::gamma_alphas(const RealVector& ga_alphas)
{ dpRep->gammaAlphas = ga_alphas; }


inline void DistributionParams::gamma_alpha(const Real& ga_alpha, size_t i)
{ dpRep->gammaAlphas[i] = ga_alpha; }


inline const RealVector& DistributionParams::gamma_betas() const
{ return dpRep->gammaBetas; }


inline const Real& DistributionParams::gamma_beta(size_t i) const
{ return dpRep->gammaBetas[i]; }


inline void DistributionParams::gamma_betas(const RealVector& ga_betas)
{ dpRep->gammaBetas = ga_betas; }


inline void DistributionParams::gamma_beta(const Real& ga_beta, size_t i)
{ dpRep->gammaBetas[i] = ga_beta; }


inline const RealVector& DistributionParams::gumbel_alphas() const
{ return dpRep->gumbelAlphas; }


inline const Real& DistributionParams::gumbel_alpha(size_t i) const
{ return dpRep->gumbelAlphas[i]; }


inline void DistributionParams::gumbel_alphas(const RealVector& gu_alphas)
{ dpRep->gumbelAlphas = gu_alphas; }


inline void DistributionParams::gumbel_alpha(const Real& gu_alpha, size_t i)
{ dpRep->gumbelAlphas[i] = gu_alpha; }


inline const RealVector& DistributionParams::gumbel_betas() const
{ return dpRep->gumbelBetas; }


inline const Real& DistributionParams::gumbel_beta(size_t i) const
{ return dpRep->gumbelBetas[i]; }


inline void DistributionParams::gumbel_betas(const RealVector& gu_betas)
{ dpRep->gumbelBetas = gu_betas; }


inline void DistributionParams::gumbel_beta(const Real& gu_beta, size_t i)
{ dpRep->gumbelBetas[i] = gu_beta; }


inline const RealVector& DistributionParams::frechet_alphas() const
{ return dpRep->frechetAlphas; }


inline const Real& DistributionParams::frechet_alpha(size_t i) const
{ return dpRep->frechetAlphas[i]; }


inline void DistributionParams::frechet_alphas(const RealVector& f_alphas)
{ dpRep->frechetAlphas = f_alphas; }


inline void DistributionParams::frechet_alpha(const Real& f_alpha, size_t i)
{ dpRep->frechetAlphas[i] = f_alpha; }


inline const RealVector& DistributionParams::frechet_betas() const
{ return dpRep->frechetBetas; }


inline const Real& DistributionParams::frechet_beta(size_t i) const
{ return dpRep->frechetBetas[i]; }


inline void DistributionParams::frechet_betas(const RealVector& f_betas)
{ dpRep->frechetBetas = f_betas; }


inline void DistributionParams::frechet_beta(const Real& f_beta, size_t i)
{ dpRep->frechetBetas[i] = f_beta; }


inline const RealVector& DistributionParams::weibull_alphas() const
{ return dpRep->weibullAlphas; }


inline const Real& DistributionParams::weibull_alpha(size_t i) const
{ return dpRep->weibullAlphas[i]; }


inline void DistributionParams::weibull_alphas(const RealVector& w_alphas)
{ dpRep->weibullAlphas = w_alphas; }


inline void DistributionParams::weibull_alpha(const Real& alpha, size_t i)
{ dpRep->weibullAlphas[i] = alpha; }


inline const RealVector& DistributionParams::weibull_betas() const
{ return dpRep->weibullBetas; }


inline const Real& DistributionParams::weibull_beta(size_t i) const
{ return dpRep->weibullBetas[i]; }


inline void DistributionParams::weibull_betas(const RealVector& w_betas)
{ dpRep->weibullBetas = w_betas; }


inline void DistributionParams::weibull_beta(const Real& beta, size_t i)
{ dpRep->weibullBetas[i] = beta; }


inline const RealVectorArray& DistributionParams::histogram_bin_pairs() const
{ return dpRep->histogramBinPairs; }


inline const RealVector& DistributionParams::histogram_bin_pairs(size_t i) const
{ return dpRep->histogramBinPairs[i]; }


inline void DistributionParams::
histogram_bin_pairs(const RealVectorArray& h_bin_pairs)
{ dpRep->histogramBinPairs = h_bin_pairs; }


inline void DistributionParams::
histogram_bin_pairs(const RealVector& h_bin_pr, size_t i)
{ dpRep->histogramBinPairs[i] = h_bin_pr; }


inline const RealVector& DistributionParams::poisson_lambdas() const
{ return dpRep->poissonLambdas; }


inline const Real& DistributionParams::poisson_lambda(size_t i) const
{ return dpRep->poissonLambdas[i]; }


inline void DistributionParams::poisson_lambdas(const RealVector& p_lambdas)
{ dpRep->poissonLambdas = p_lambdas; }


inline void DistributionParams::poisson_lambda(const Real& p_lambda, size_t i)
{ dpRep->poissonLambdas[i] = p_lambda; }


inline const RealVector& DistributionParams::
binomial_probabilities_per_trial() const
{ return dpRep->binomialProbPerTrial; }


inline const Real& DistributionParams::
binomial_probabilities_per_trial(size_t i) const
{ return dpRep->binomialProbPerTrial[i]; }


inline void DistributionParams::
binomial_probabilities_per_trial(const RealVector& probs_per_tr)
{ dpRep->binomialProbPerTrial = probs_per_tr; }


inline void DistributionParams::
binomial_probability_per_trial(const Real& prob_per_tr, size_t i)
{ dpRep->binomialProbPerTrial[i] = prob_per_tr; }


inline const IntVector& DistributionParams::binomial_num_trials() const
{ return dpRep->binomialNumTrials; }


inline int DistributionParams::binomial_num_trials(size_t i) const
{ return dpRep->binomialNumTrials[i]; }


inline void DistributionParams::binomial_num_trials(const IntVector& num_tr)
{ dpRep->binomialNumTrials = num_tr; }


inline void DistributionParams::binomial_num_trials(int num_tr, size_t i)
{ dpRep->binomialNumTrials[i] = num_tr; }


inline const RealVector& DistributionParams::
negative_binomial_probabilities_per_trial() const
{ return dpRep->negBinomialProbPerTrial; }


inline const Real& DistributionParams::
negative_binomial_probabilities_per_trial(size_t i) const
{ return dpRep->negBinomialProbPerTrial[i]; }


inline void DistributionParams::
negative_binomial_probabilities_per_trial(const RealVector& probs_per_tr)
{
  if (dpRep->negBinomialProbPerTrial.empty())//Teuchos operator= does not resize
    dpRep->negBinomialProbPerTrial.sizeUninitialized(probs_per_tr.length());
  dpRep->negBinomialProbPerTrial = probs_per_tr; }


inline void DistributionParams::
negative_binomial_probability_per_trial(const Real& prob_per_tr, size_t i)
{ dpRep->negBinomialProbPerTrial[i] = prob_per_tr; }


inline const IntVector& DistributionParams::negative_binomial_num_trials() const
{ return dpRep->negBinomialNumTrials; }


inline int DistributionParams::negative_binomial_num_trials(size_t i) const
{ return dpRep->negBinomialNumTrials[i]; }


inline void DistributionParams::
negative_binomial_num_trials(const IntVector& num_tr)
{ dpRep->negBinomialNumTrials = num_tr; }


inline void DistributionParams::
negative_binomial_num_trials(int num_tr, size_t i)
{ dpRep->negBinomialNumTrials[i] = num_tr; }


inline const RealVector& DistributionParams::
geometric_probabilities_per_trial() const
{ return dpRep->geometricProbPerTrial; }


inline const Real& DistributionParams::
geometric_probabilities_per_trial(size_t i) const
{ return dpRep->geometricProbPerTrial[i]; }


inline void DistributionParams::
geometric_probabilities_per_trial(const RealVector& probs_per_tr)
{ dpRep->geometricProbPerTrial = probs_per_tr; }


inline void DistributionParams::
geometric_probability_per_trial(const Real& prob_per_tr, size_t i)
{ dpRep->geometricProbPerTrial[i] = prob_per_tr; }


inline const IntVector& DistributionParams::
hypergeometric_total_population() const
{ return dpRep->hyperGeomTotalPopulation; }


inline int DistributionParams::hypergeometric_total_population(size_t i) const
{ return dpRep->hyperGeomTotalPopulation[i]; }


inline void DistributionParams::
hypergeometric_total_population(const IntVector& total_pop)
{ dpRep->hyperGeomTotalPopulation = total_pop; }


inline void DistributionParams::
hypergeometric_total_population(int total_pop, size_t i)
{ dpRep->hyperGeomTotalPopulation[i] = total_pop; }


inline const IntVector& DistributionParams::
hypergeometric_selected_population() const
{ return dpRep->hyperGeomSelectedPopulation; }


inline int DistributionParams::
hypergeometric_selected_population(size_t i) const
{ return dpRep->hyperGeomSelectedPopulation[i]; }


inline void DistributionParams::
hypergeometric_selected_population(const IntVector& sel_pop)
{ dpRep->hyperGeomSelectedPopulation = sel_pop; }


inline void DistributionParams::
hypergeometric_selected_population(int sel_pop, size_t i)
{ dpRep->hyperGeomSelectedPopulation[i] = sel_pop; }


inline const IntVector& DistributionParams::hypergeometric_num_drawn() const
{ return dpRep->hyperGeomNumDrawn; }


inline int DistributionParams::hypergeometric_num_drawn(size_t i) const
{ return dpRep->hyperGeomNumDrawn[i]; }


inline void DistributionParams::
hypergeometric_num_drawn(const IntVector& num_drawn)
{ dpRep->hyperGeomNumDrawn = num_drawn; }


inline void DistributionParams::
hypergeometric_num_drawn(int num_drawn, size_t i)
{ dpRep->hyperGeomNumDrawn[i] = num_drawn; }


inline const RealVectorArray& DistributionParams::histogram_point_pairs() const
{ return dpRep->histogramPointPairs; }


inline const RealVector& DistributionParams::
histogram_point_pairs(size_t i) const
{ return dpRep->histogramPointPairs[i]; }


inline void DistributionParams::
histogram_point_pairs(const RealVectorArray& h_pt_pairs)
{ dpRep->histogramPointPairs = h_pt_pairs; }


inline void DistributionParams::
histogram_point_pairs(const RealVector& h_pt_pairs_i, size_t i)
{ dpRep->histogramPointPairs[i] = h_pt_pairs_i; }


inline const RealVectorArray& DistributionParams::interval_probabilities() const
{ return dpRep->intervalBasicProbs; }


inline const RealVector& DistributionParams::
interval_probabilities(size_t i) const
{ return dpRep->intervalBasicProbs[i]; }


inline void DistributionParams::
interval_probabilities(const RealVectorArray& int_probs)
{ dpRep->intervalBasicProbs = int_probs; }


inline void DistributionParams::
interval_probabilities(const RealVector& int_probs_i, size_t i)
{ dpRep->intervalBasicProbs[i] = int_probs_i; }


inline const RealVectorArray& DistributionParams::interval_bounds() const
{ return dpRep->intervalBounds; }


inline const RealVector& DistributionParams::interval_bounds(size_t i) const
{ return dpRep->intervalBounds[i]; }


inline void DistributionParams::
interval_bounds(const RealVectorArray& int_bounds)
{ dpRep->intervalBounds = int_bounds; }


inline void DistributionParams::
interval_bounds(const RealVector& int_bounds_i, size_t i)
{ dpRep->intervalBounds[i] = int_bounds_i; }


inline const RealSymMatrix& DistributionParams::uncertain_correlations() const
{ return dpRep->uncertainCorrelations; }


inline void DistributionParams::
uncertain_correlations(const RealSymMatrix& uncertain_corr)
{ dpRep->uncertainCorrelations = uncertain_corr; }


inline void DistributionParams::update(const DistributionParams& dp)
{
  if (!dpRep) // if no rep, create a new instance
    copy(dp);
  else {      // update data of existing instance
    normal_means(dp.normal_means());
    normal_std_deviations(dp.normal_std_deviations());
    normal_lower_bounds(dp.normal_lower_bounds());
    normal_upper_bounds(dp.normal_upper_bounds());
    lognormal_means(dp.lognormal_means());
    lognormal_std_deviations(dp.lognormal_std_deviations());
    lognormal_lambdas(dp.lognormal_lambdas());
    lognormal_zetas(dp.lognormal_zetas());
    lognormal_error_factors(dp.lognormal_error_factors());
    lognormal_lower_bounds(dp.lognormal_lower_bounds());
    lognormal_upper_bounds(dp.lognormal_upper_bounds());
    uniform_lower_bounds(dp.uniform_lower_bounds());
    uniform_upper_bounds(dp.uniform_upper_bounds());
    loguniform_lower_bounds(dp.loguniform_lower_bounds());
    loguniform_upper_bounds(dp.loguniform_upper_bounds());
    triangular_modes(dp.triangular_modes());
    triangular_lower_bounds(dp.triangular_lower_bounds());
    triangular_upper_bounds(dp.triangular_upper_bounds());
    exponential_betas(dp.exponential_betas());
    beta_alphas(dp.beta_alphas());
    beta_betas(dp.beta_betas());
    beta_lower_bounds(dp.beta_lower_bounds());
    beta_upper_bounds(dp.beta_upper_bounds());
    gamma_alphas(dp.gamma_alphas());
    gamma_betas(dp.gamma_betas());
    gumbel_alphas(dp.gumbel_alphas());
    gumbel_betas(dp.gumbel_betas());
    frechet_alphas(dp.frechet_alphas());
    frechet_betas(dp.frechet_betas());
    weibull_alphas(dp.weibull_alphas());
    weibull_betas(dp.weibull_betas());
    histogram_bin_pairs(dp.histogram_bin_pairs());
    poisson_lambdas(dp.poisson_lambdas());
    binomial_probabilities_per_trial(dp.binomial_probabilities_per_trial());
    binomial_num_trials(dp.binomial_num_trials());
    negative_binomial_probabilities_per_trial(
      dp.negative_binomial_probabilities_per_trial());
    negative_binomial_num_trials(dp.negative_binomial_num_trials());
    geometric_probabilities_per_trial(dp.geometric_probabilities_per_trial());
    hypergeometric_total_population(dp.hypergeometric_total_population());
    hypergeometric_selected_population(dp.hypergeometric_selected_population());
    hypergeometric_num_drawn(dp.hypergeometric_num_drawn());
    histogram_point_pairs(dp.histogram_point_pairs());
    interval_probabilities(dp.interval_probabilities());
    interval_bounds(dp.interval_bounds());
    uncertain_correlations(dp.uncertain_correlations());
  }
}


inline bool DistributionParams::is_null() const
{ return (dpRep) ? false : true; }

} // namespace Pecos

#endif

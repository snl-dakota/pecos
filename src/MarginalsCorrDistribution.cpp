/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "MarginalsCorrDistribution.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: MarginalsCorrDistribution.cpp 4768 2007-12-17 17:49:32Z mseldre $";

//#define DEBUG


namespace Pecos {

void MarginalsCorrDistribution::initialize(const ShortArray& rv_types)
{
  ranVarTypes = rv_types;

  // construction of x-space random variables occurs once (updates to
  // distribution parameters can occur repeatedly)
  size_t i, num_v = rv_types.size();
  randomVars.resize(num_v);
  for (i=0; i<num_v; ++i)
    randomVars[i] = RandomVariable(rv_types[i]);
}


void MarginalsCorrDistribution::
parameters(size_t v, const ShortArray& dist_params, const RealVector& values)
{
  RandomVariable& random_var = randomVars[v];
  size_t i, num_params = std::min(dist_params.size(), values.length());
  for (i=0; i<num_params; ++i)
    random_var.parameter(dist_params[i], values[i]);
}


void MarginalsCorrDistribution::
parameters(size_t start_v, size_t num_v, short dist_param,
	   const RealVector& values)
{
  size_t i, v, num_updates = std::min(values.length(), num_v);
  for (i=0, v=start_v; i<num_updates; ++i, ++v)
    randomVars[v].parameter(dist_param, values[i]);
}


void MarginalsCorrDistribution::
parameters(short rv_type, short dist_param, const RealVector& values)
{
  size_t rv, num_rv = ranVarTypes.size(), cntr = 0, num_vals = values.length();
  for (rv=0; i<num_rv; ++rv)
    if (ranVarTypes[rv] == rv_type && cntr < num_vals)
      randomVars[rv].parameter(dist_param, values[cntr++]);
}


void MarginalsCorrDistribution::correlations(const RealSymMatrix& corr)
{
  corrMatrix = corr;

  size_t i, j, num_rv = corr.numRows();
  correlationFlag = false;
  for (i=1; i<num_rv; i++)
    for (j=0; j<i; j++)
      if (std::abs(corr(i,j)) > SMALL_NUMBER)
	{ correlationFlag = true; break; }
}


void MarginalsCorrDistribution::
expand_correlation_matrix(size_t num_lead_v, size_t num_prob_v,
			  size_t num_trail_v)
{
  if (!correlationFlag)
    return;
  size_t num_prev_v  = corrMatrix.numRows(),
         num_total_v = num_lead_v + num_prob_v + num_trail_v;
  if (num_prev_v == num_total_v)
    return;
  else if (num_prev_v != num_prob_v) { // old->new subset assumed below
    PCerr << "\nError: unexpected matrix size (" << num_prev_v
	  << ") in MarginalsCorrDistribution::expand_correlation_matrix()."
	  << std::endl;
    abort_handler(-1);
  }

  // expand from num_prob_v to num_total_v

  size_t i, j, offset;
  RealSymMatrix old_corr_matrix(corrMatrix); // copy
  corrMatrix.shape(num_total_v); // init to zero
  for (i=0; i<num_lead_v; i++)
    corrMatrix(i,i) = 1.;
  offset = num_lead_v;
  for (i=0; i<num_prob_v; i++)
    for (j=0; j<=i; j++)
      corrMatrix(i+offset, j+offset) = old_corr_matrix(i,j);
  offset += num_prob_v;
  for (i=0; i<num_trail_v; i++)
    corrMatrix(i+offset, i+offset) = 1.;
}


void MarginalsCorrDistribution::
contract_correlation_matrix(size_t num_lead_v, size_t num_prob_v,
			    size_t num_trail_v)
{
  if (!correlationFlag)
    return;
  size_t num_prev_v  = corrMatrix.numRows(),
         num_total_v = num_lead_v + num_prob_v + num_trail_v;
  if (num_prev_v == num_prob_v)
    return;
  else if (num_prev_v != num_total_v) {
    PCerr << "\nError: unexpected matrix size (" << num_prev_v
	  << ") in MarginalsCorrDistribution::contract_correlation_matrix()."
	  << std::endl;
    abort_handler(-1);
  }

  // contract from num_total_v to num_prob_v

  size_t i, j;
  RealSymMatrix old_corr_matrix(corrMatrix); // copy
  corrMatrix.shape(num_prob_v); // init to zero
  for (i=0; i<num_prob_v; i++)
    for (j=0; j<=i; j++)
      corrMatrix(i, j) = old_corr_matrix(i+num_lead_v,j+num_lead_v);
}


/*
void MarginalsCorrDistribution::
initialize_random_variable_parameters(const RealVector& cd_l_bnds,
				      const RealVector& cd_u_bnds,
				      const AleatoryDistParams& adp,
				      const EpistemicDistParams& edp,
				      const RealVector& cs_l_bnds,
				      const RealVector& cs_u_bnds)
{
  if (mvDistRep) // envelope fwd to letter
    mvDistRep->
      initialize_random_variable_parameters(cd_l_bnds, cd_u_bnds, adp, edp,
					    cs_l_bnds, cs_u_bnds);
  else {
    size_t i, num_v = randomVarsX.size(), cd_cntr = 0, n_cntr = 0, ln_cntr = 0,
      u_cntr = 0, lu_cntr = 0, t_cntr = 0, e_cntr = 0, b_cntr = 0, ga_cntr = 0,
      gu_cntr = 0, f_cntr = 0, w_cntr = 0, h_cntr = 0, p_cntr = 0, bi_cntr = 0,
      nbi_cntr = 0, ge_cntr = 0, hge_cntr = 0, hpi_cntr = 0, hps_cntr = 0,
      hpr_cntr = 0, ci_cntr = 0, cs_cntr = 0;
    Real dbl_inf = std::numeric_limits<Real>::infinity();
    RandomVariable* rv_rep_i;
    const RealVector& ln_means     = adp.lognormal_means();
    const RealVector& ln_std_devs  = adp.lognormal_std_deviations();
    const RealVector& ln_lambdas   = adp.lognormal_lambdas();
    const RealVector& ln_zetas     = adp.lognormal_zetas();
    const RealVector& ln_err_facts = adp.lognormal_error_factors();
    for (i=0; i<num_v; ++i) {
      rv_rep_i = randomVarsX[i].random_variable_rep();
      switch (rv_rep_i->type()) {
      case CONTINUOUS_DESIGN: {
	Real lwr = cd_l_bnds[cd_cntr], upr = cd_u_bnds[cd_cntr];
	if (lwr == -dbl_inf || upr == dbl_inf) {
	  PCerr << "Error: bounds specification required for design variable "
	       << "transformation to standard uniform." << std::endl;
	  abort_handler(-1);
	}
	((UniformRandomVariable*)rv_rep_i)->update(lwr, upr);
	++cd_cntr; break;
      }
      case NORMAL:
	((NormalRandomVariable*)rv_rep_i)->
	  update(adp.normal_mean(n_cntr), adp.normal_std_deviation(n_cntr));
	++n_cntr; break;
      case BOUNDED_NORMAL:
	((BoundedNormalRandomVariable*)rv_rep_i)->
	  update(adp.normal_mean(n_cntr), adp.normal_std_deviation(n_cntr),
		 adp.normal_lower_bound(n_cntr),adp.normal_upper_bound(n_cntr));
	++n_cntr; break;
      case LOGNORMAL: {
	Real lambda, zeta;
	params_from_lognormal_spec(ln_means, ln_std_devs, ln_lambdas, ln_zetas,
				   ln_err_facts, ln_cntr, lambda, zeta);
        ((LognormalRandomVariable*)rv_rep_i)->update(lambda, zeta);
	++ln_cntr; break;
      }
      case BOUNDED_LOGNORMAL: {
	Real lambda, zeta;
	params_from_lognormal_spec(ln_means, ln_std_devs, ln_lambdas, ln_zetas,
				   ln_err_facts, ln_cntr, lambda, zeta);
	((BoundedLognormalRandomVariable*)rv_rep_i)->
	  update(lambda, zeta, adp.lognormal_lower_bound(ln_cntr),
		 adp.lognormal_upper_bound(ln_cntr));
	++ln_cntr; break;
      }
      case UNIFORM:
	((UniformRandomVariable*)rv_rep_i)->
	  update(adp.uniform_lower_bound(u_cntr),
		 adp.uniform_upper_bound(u_cntr));
	++u_cntr; break;
      case LOGUNIFORM:
	((LoguniformRandomVariable*)rv_rep_i)->
	  update(adp.loguniform_lower_bound(lu_cntr),
		 adp.loguniform_upper_bound(lu_cntr));
	++lu_cntr; break;
      case TRIANGULAR:
	((TriangularRandomVariable*)rv_rep_i)->
	  update(adp.triangular_lower_bound(t_cntr),
		 adp.triangular_mode(t_cntr),
		 adp.triangular_upper_bound(t_cntr));
	++t_cntr; break;
      case EXPONENTIAL:
	((ExponentialRandomVariable*)rv_rep_i)->
	  update(adp.exponential_beta(e_cntr));
	++e_cntr; break;
      case BETA:
	((BetaRandomVariable*)rv_rep_i)->
	  update(adp.beta_alpha(b_cntr), adp.beta_beta(b_cntr),
		 adp.beta_lower_bound(b_cntr), adp.beta_upper_bound(b_cntr));
	++b_cntr; break;
      case GAMMA:
	((GammaRandomVariable*)rv_rep_i)->
	  update(adp.gamma_alpha(ga_cntr), adp.gamma_beta(ga_cntr));
	+ga_cntr; break;
      case GUMBEL:
	((GumbelRandomVariable*)rv_rep_i)->
	  update(adp.gumbel_alpha(gu_cntr), adp.gumbel_beta(gu_cntr));
	++gu_cntr; break;
      case FRECHET:
	((FrechetRandomVariable*)rv_rep_i)->
	  update(adp.frechet_alpha(f_cntr), adp.frechet_beta(f_cntr));
	++f_cntr; break;
      case WEIBULL:
	((WeibullRandomVariable*)rv_rep_i)->
	  update(adp.weibull_alpha(w_cntr), adp.weibull_beta(w_cntr));
	++w_cntr; break;
      case HISTOGRAM_BIN:
	((HistogramBinRandomVariable*)rv_rep_i)->
	  update(adp.histogram_bin_pairs(h_cntr));
	++h_cntr; break;

      // discrete int aleatory uncertain
      case POISSON:
	((PoissonRandomVariable*)rv_rep_i)->update(adp.poisson_lambda(p_cntr));
	++p_cntr; break;
      case BINOMIAL:
 	((BinomialRandomVariable*)rv_rep_i)->
	  update(adp.binomial_num_trials(bi_cntr),
		 adp.binomial_probability_per_trial(bi_cntr));
	++bi_cntr; break;
      case NEGATIVE_BINOMIAL:
 	((NegBinomialRandomVariable*)rv_rep_i)->
	  update(adp.negative_binomial_num_trials(nbi_cntr),
		 adp.negative_binomial_probability_per_trial(nbi_cntr));
	++nbi_cntr; break;
      case GEOMETRIC:
 	((GeometricRandomVariable*)rv_rep_i)->
	  update(adp.geometric_probability_per_trial(ge_cntr));
	++ge_cntr; break;
      case HYPERGEOMETRIC:
 	((HypergeometricRandomVariable*)rv_rep_i)->
	  update(adp.hypergeometric_total_population(hge_cntr),
		 adp.hypergeometric_selected_population(hge_cntr),
		 adp.hypergeometric_num_drawn(hge_cntr));
	++hge_cntr; break;
      case HISTOGRAM_PT_INT:
 	((HistogramPtRandomVariable*)rv_rep_i)->
	  update(adp.histogram_point_int_pairs(hpi_cntr));
	++hpi_cntr; break;

      // discrete string aleatory uncertain
      case HISTOGRAM_PT_STRING:
 	((HistogramPtRandomVariable*)rv_rep_i)->
	  update(adp.histogram_point_string_pairs(hps_cntr));
	++hps_cntr; break;

      // discrete real aleatory uncertain
      case HISTOGRAM_PT_REAL:
 	((HistogramPtRandomVariable*)rv_rep_i)->
	  update(adp.histogram_point_real_pairs(hpr_cntr));
	++hpr_cntr; break;

      case CONTINUOUS_INTERVAL: {
	const RealRealPairRealMap& ci_bp
	  = edp.continuous_interval_basic_probabilities(ci_cntr);
	// process interval sets for overall bounds; since intervals are sorted,
	// lb should be in 1st interval but go ahead and process lb same as ub
	Real lb = dbl_inf, ub = -dbl_inf;
	RealRealPairRealMap::const_iterator cit;
	for (cit=ci_bp.begin(); cit!=ci_bp.end(); ++cit) {
	  const RealRealPair& bnds = cit->first;
	  if (bnds.first  < lb) lb = bnds.first;
	  if (bnds.second > ub) ub = bnds.second;
	}
	((UniformRandomVariable*)rv_rep_i)->update(lb, ub);
	++ci_cntr; break;
      }

      // discrete int epistemic uncertain

      // discrete string epistemic uncertain

      // discrete real epistemic uncertain

      case CONTINUOUS_STATE: {
	Real lwr = cs_l_bnds[cs_cntr], upr = cs_u_bnds[cs_cntr];
	if (lwr == -dbl_inf || upr == dbl_inf) {
	  PCerr << "Error: bounds specification required for state variable "
	       << "transformation to standard uniform." << std::endl;
	  abort_handler(-1);
	}
	((UniformRandomVariable*)rv_rep_i)->update(lwr, upr);
	++cs_cntr; break;
      }
      case STD_NORMAL:      ++n_cntr; break;
      case STD_UNIFORM:     ++u_cntr; break;
      case STD_EXPONENTIAL: ++e_cntr; break;
      case STD_BETA:        ++b_cntr; break;
      case STD_GAMMA:      ++ga_cntr; break;
      //default:
      }
    }
  }
}


RealRealPairArray MarginalsCorrDistribution::u_moments() const
{
  if (mvDistRep) return mvDistRep->u_moments();
  else {
    size_t i, num_v = randomVars.size();
    RealRealPairArray u_mom(num_v);
    Real unif_stdev = 1./std::sqrt(3.);
    for (i=0; i<num_v; ++i)
      switch (ranVarTypesU[i]) {
      case STD_NORMAL:      u_mom[i] = RealRealPair(0., 1.);         break;
      case STD_UNIFORM:     u_mom[i] = RealRealPair(0., unif_stdev); break;
      case STD_EXPONENTIAL: u_mom[i] = RealRealPair(1., 1.);         break;
      case STD_BETA: {
	check_type(i, BETA);
	Real mean, stdev;
	BetaRandomVariable::
	  moments_from_params(randomVars[i].parameter(BE_ALPHA),
			      randomVars[i].parameter(BE_BETA), -1., 1.,
			      mean, stdev);
        u_mom[i] = RealRealPair(mean, stdev); break;
      }
      case STD_GAMMA: {
	check_type(i, GAMMA);
	Real mean, stdev;
	GammaRandomVariable::
	  moments_from_params(randomVars[i].parameter(GA_ALPHA), 1.,
			      mean, stdev);
        u_mom[i] = RealRealPair(mean, stdev); break;
      }
      default: // no transformation (e.g., PCE w/ numerically-generated basis)
	check_type(i, ranVarTypesU[i]);
	u_mom[i] = randomVars[i].moments(); break;
      }
    return u_mom;
  }
}


RealRealPairArray MarginalsCorrDistribution::u_bounds() const
{
  if (mvDistRep) return mvDistRep->u_bounds();
  else {
    size_t i, num_v = randomVars.size();
    RealRealPairArray u_bnds(num_v);
    Real dbl_inf = std::numeric_limits<Real>::infinity();
    for (i=0; i<num_v; ++i)
      switch (ranVarTypesU[i]) {
      case STD_NORMAL:
	u_bnds[i] = RealRealPair(-dbl_inf, dbl_inf); break;
      case STD_UNIFORM:     case STD_BETA:
	u_bnds[i] = RealRealPair(-1., 1.);           break;
      case STD_EXPONENTIAL: case STD_GAMMA:
	u_bnds[i] = RealRealPair(0., dbl_inf);       break;
      default: // no transformation (e.g., PCE w/ numerically-generated basis)
	check_type(i, ranVarTypesU[i]);
	u_bnds[i] = randomVars[i].bounds();         break;
      }
    return u_bnds;
  }
}


Real MarginalsCorrDistribution::u_pdf(Real u_val, size_t i) const
{
  if (mvDistRep) return mvDistRep->u_pdf(u_val, i);
  else {
    // can only use randomVars[i].standard_pdf() for cases where u_type is a
    // standardized form of the x_type.  For STD_NORMAL and STD_UNIFORM, many
    // x_types can be mapped to these u_types, so use global utility fns
    // whenever there are no auxilliary parameters to manage.
    switch (ranVarTypesU[i]) {
    case STD_NORMAL:  return  NormalRandomVariable::std_pdf(u_val); break;
    case STD_UNIFORM: return UniformRandomVariable::std_pdf();      break;
    case STD_EXPONENTIAL:
      return ExponentialRandomVariable::std_pdf(u_val);             break;
    case STD_BETA:
      check_type(i, BETA);  return randomVars[i].standard_pdf(u_val); break;
    case STD_GAMMA:
      check_type(i, GAMMA); return randomVars[i].standard_pdf(u_val); break;
    default: // no transformation (e.g., PCE with numerically-generated bases)
      check_type(i, ranVarTypesU[i]);
      return randomVars[i].pdf(u_val); break;
    }
  }
}


Real MarginalsCorrDistribution::u_log_pdf(Real u_val, size_t i) const
{
  if (mvDistRep) return mvDistRep->u_log_pdf(u_val, i);
  else {
    switch (ranVarTypesU[i]) {
    case STD_NORMAL:  return  NormalRandomVariable::log_std_pdf(u_val); break;
    case STD_UNIFORM: return UniformRandomVariable::log_std_pdf();      break;
    case STD_EXPONENTIAL:
      return ExponentialRandomVariable::log_std_pdf(u_val);             break;
    case STD_BETA:
      check_type(i, BETA);  // need alphaStat,betaStat for variable i
      return randomVars[i].log_standard_pdf(u_val); break;
    case STD_GAMMA:
      check_type(i, GAMMA); // need alphaStat for variable i
      return randomVars[i].log_standard_pdf(u_val); break;
    default: // no transformation (e.g., PCE with numerically-generated bases)
      check_type(i, ranVarTypesU[i]);
      return randomVars[i].log_pdf(u_val);          break;
    }
  }
}


Real MarginalsCorrDistribution::u_log_pdf_gradient(Real u_val, size_t i) const
{
  if (mvDistRep) return mvDistRep->u_log_pdf_gradient(u_val, i);
  else {
    switch (ranVarTypesU[i]) {
    case STD_NORMAL:
      return      NormalRandomVariable::log_std_pdf_gradient(u_val); break;
    case STD_UNIFORM:
      return     UniformRandomVariable::log_std_pdf_gradient(); break;
    case STD_EXPONENTIAL:
      return ExponentialRandomVariable::log_std_pdf_gradient(); break;
    case STD_BETA:
      check_type(i, BETA);  // need alphaStat,betaStat for variable i
      return randomVars[i].log_standard_pdf_gradient(u_val); break;
    case STD_GAMMA:
      check_type(i, GAMMA); // need alphaStat for variable i
      return randomVars[i].log_standard_pdf_gradient(u_val); break;
    default: // no transformation (e.g., PCE with numerically-generated bases)
      check_type(i, ranVarTypesU[i]);
      return randomVars[i].log_pdf_gradient(u_val);          break;
    }
  }
}


Real MarginalsCorrDistribution::u_log_pdf_hessian(Real u_val, size_t i) const
{
  if (mvDistRep) return mvDistRep->u_log_pdf_hessian(u_val, i);
  else {
    switch (ranVarTypesU[i]) {
    case STD_NORMAL:
      return      NormalRandomVariable::log_std_pdf_hessian(); break;
    case STD_UNIFORM:
      return     UniformRandomVariable::log_std_pdf_hessian(); break;
    case STD_EXPONENTIAL:
      return ExponentialRandomVariable::log_std_pdf_hessian(); break;
    case STD_BETA:
      check_type(i, BETA);  // need alphaStat,betaStat for variable i
      return randomVars[i].log_standard_pdf_hessian(u_val); break;
    case STD_GAMMA:
      check_type(i, GAMMA); // need alphaStat for variable i
      return randomVars[i].log_standard_pdf_hessian(u_val); break;
    default: // no transformation (e.g., PCE with numerically-generated bases)
      check_type(i, ranVarTypesU[i]);
      return randomVars[i].log_pdf_hessian(u_val);          break;
    }
  }
}


inline Real MarginalsCorrDistribution::u_pdf(const RealVector& u_pt) const
{
  if (probTransRep) return probTransRep->u_pdf(u_pt);
  else {
    // u-space is independent -> use product of marginals
    size_t i, num_v = randomVarsX.size();
    Real density = 1.;
    for (i=0; i<num_v; ++i)
      density *= u_pdf(u_pt[i], i);
    return density;
  }
}


inline Real MarginalsCorrDistribution::u_log_pdf(const RealVector& u_pt) const
{
  if (probTransRep) return probTransRep->u_log_pdf(u_pt);
  else {
    // u-space is independent -> use sum of log marginals
    size_t i, num_v = randomVarsX.size();
    Real log_density = 0.;
    for (i=0; i<num_v; ++i)
      log_density += u_log_pdf(u_pt[i], i);
    return log_density;
  }
}
*/

} // namespace Pecos

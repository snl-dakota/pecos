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

void MarginalsCorrDistribution::
initialize_types(const ShortArray& rv_types, const BitArray& active_vars)
{
  ranVarTypes = rv_types;
  activeVars  = active_vars;

  // construction of x-space random variables occurs once (updates to
  // distribution parameters can occur repeatedly)
  size_t i, num_v = rv_types.size();
  randomVars.resize(num_v);
  for (i=0; i<num_v; ++i)
    randomVars[i] = RandomVariable(rv_types[i]);
}


void MarginalsCorrDistribution::
initialize_correlations(const RealSymMatrix& corr, const BitArray& active_corr)
{
  corrMatrix = corr;
  activeCorr = active_corr; // active RV subset for correlation matrix

  correlationFlag = false;
  size_t num_corr = corr.numRows();
  if (num_corr == 0) return;

  size_t num_rv = randomVars.size();
  bool no_mask = activeCorr.empty();
  if (no_mask) {
    if (num_corr != num_rv) {
      PCerr << "Error: correlation matrix size (" << num_corr
	    << ") inconsistent with number of random variables (" << num_rv
	    << ")." << std::endl;
      abort_handler(-1);
    }
  }
  else {
    if (num_corr != activeCorr.count()) {
      PCerr << "Error: correlation matrix size (" << num_corr
	    << ") inconsistent with active correlation subset ("
	    << activeCorr.count() << ")." << std::endl;
      abort_handler(-1);
    }
  }

  size_t i, j, cntr_i, cntr_j;
  for (i=0, cntr_i=0; i<num_rv; ++i) {
    if (no_mask || activeCorr[i]) {
      for (j=0, cntr_j=0; j<i; ++j) {
	if (no_mask || activeCorr[j]) {
	  if (std::abs(corr(cntr_i, cntr_j)) > SMALL_NUMBER)
	    correlationFlag = true;
	  ++cntr_j;
	}
      }
      ++cntr_i;
    }
    if (correlationFlag) break;
  }
}


void MarginalsCorrDistribution::
pull_distribution_parameters(MultivariateDistribution* mv_dist_rep)
{
  const std::vector<RandomVariable>& rv_in = mv_dist_rep->random_variables();
  const ShortArray& rv_types_in = mv_dist_rep->random_variable_types();
  size_t i, num_rv = ranVarTypes.size();

  for (i=0; i<num_rv; ++i) {
    RandomVariable&       push_rv = randomVars[i];
    const RandomVariable& pull_rv = rv_in[i];
    short push_type = ranVarTypes[i], pull_type = rv_types_in[i];
    switch (push_type) {

    // push RV are fully standard: no updates to perform
    case STD_NORMAL:  case STD_UNIFORM:  case STD_EXPONENTIAL:           break;

    // push RV have standardized scale params; copy shape params
    case STD_BETA:
      push_rv.push_parameter(BE_ALPHA, pull_rv.pull_parameter<Real>(BE_ALPHA));
      push_rv.push_parameter(BE_BETA,  pull_rv.pull_parameter<Real>(BE_BETA));
      break;
    case STD_GAMMA:
      push_rv.push_parameter(GA_ALPHA, pull_rv.pull_parameter<Real>(GA_ALPHA));
      break;

    // push RV are non-standard, pull all non-standard data from rv_in
    default:
      switch (pull_type) {
      // pull RV are fully standard: no data to pull
      case STD_NORMAL:  case STD_UNIFORM:  case STD_EXPONENTIAL:         break;

      // pull RV have standardized scale params; pull shape params
      case STD_BETA:
	push_rv.push_parameter(BE_ALPHA,pull_rv.pull_parameter<Real>(BE_ALPHA));
	push_rv.push_parameter(BE_BETA, pull_rv.pull_parameter<Real>(BE_BETA));
	break;
      case STD_GAMMA:
	push_rv.push_parameter(GA_ALPHA,pull_rv.pull_parameter<Real>(GA_ALPHA));
	break;

      // pull and push RV are non-standard; copy all params
      default:
	push_rv.copy_parameters(rv_in[i]);                               break;
      }
      break;
    }
  }
}


void MarginalsCorrDistribution::copy_rep(MultivariateDistribution* source_rep)
{
  // copy base class data
  MultivariateDistribution::copy_rep(source_rep);
  // specialization for marginals + corr
  MarginalsCorrDistribution* mcd_rep = (MarginalsCorrDistribution*)source_rep;
  initialize_types(mcd_rep->ranVarTypes, mcd_rep->activeVars);
  initialize_correlations(mcd_rep->corrMatrix, mcd_rep->activeCorr);
  pull_distribution_parameters(source_rep);
}


/* Reshaping the correlation matrix should no longer be required
   (subsets now supported with activeCorr BitArray)
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
  // Note: a reshape is not helpful due to num_lead_v

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
  // Note: a reshape is not helpful due to num_lead_v

  size_t i, j;
  RealSymMatrix old_corr_matrix(corrMatrix); // copy
  corrMatrix.shape(num_prob_v); // init to zero
  for (i=0; i<num_prob_v; i++)
    for (j=0; j<=i; j++)
      corrMatrix(i, j) = old_corr_matrix(i+num_lead_v,j+num_lead_v);
}
*/

} // namespace Pecos

/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef LHS_DRIVER_H
#define LHS_DRIVER_H

#include "pecos_data_types.hpp"
#include "pecos_global_defs.hpp"
#include "DistributionParams.hpp"


namespace Pecos {

// LHS rank array processing modes:
enum { IGNORE_RANKS, SET_RANKS, GET_RANKS, SET_GET_RANKS };


/// Driver class for Latin Hypercube Sampling (LHS)

/** This class provides common code for sampling methods which
    employ the Latin Hypercube Sampling (LHS) package from Sandia
    Albuquerque's Risk and Reliability organization. */

class LHSDriver
{
public:

  //
  //- Heading: Constructors and destructor
  //

  LHSDriver();  ///< default constructor
  LHSDriver(const String& sample_type, short sample_ranks_mode = IGNORE_RANKS,
            bool reports = true);  ///< constructor
  ~LHSDriver(); ///< destructor

  //
  //- Heading: Member functions
  //

  /// populate data when not passed through ctor
  void initialize(const String& sample_type,
		  short sample_ranks_mode = IGNORE_RANKS, bool reports = true);

  /// set randomSeed
  void seed(int seed);
  /// return randomSeed
  int seed() const;

  /// set random number generator, passing the name of the uniform
  /// generator: rnum2 or mt19937 (default).  Passed value is
  /// superceded by environment variable DAKOTA_LHS_UNIFGEN, if
  /// present
  void rng(String unif_gen);
  // return name of uniform generator
  //String rng();

  /// reseed using a deterministic sequence
  void advance_seed_sequence();

  /// generates the desired set of parameter samples from within general
  /// user-specified probabilistic distributions
  void generate_samples(const std::vector<RandomVariables>& random_vars,
			const RealSymMatrix& correlations, int num_samples,
			RealMatrix& samples, RealMatrix& sample_ranks,
			const BitArray& active_vars,
			const BitArray& active_corr);

  /// generates the desired set of parameter samples from within
  /// AleatoryDistParams and EpistemicDistParams specifications
  void generate_uncertain_samples(
    const std::vector<RandomVariables>& random_vars,
    const RealSymMatrix& correlations, int num_samples,
    RealMatrix& samples_array, bool backfill_flag = false);
  /// generates the desired set of parameter samples from within an
  /// AleatoryDistParams specification
  void generate_aleatory_samples(
    const std::vector<RandomVariables>& random_vars,
    const RealSymMatrix& correlations, int num_samples,
    RealMatrix& samples_array, bool backfill_flag = false);
  /// generates the desired set of parameter samples from within a
  /// EpistemicDistParams specification
  void generate_epistemic_samples(
    const std::vector<RandomVariables>& random_vars,
    const RealSymMatrix& correlations, int num_samples,
    RealMatrix& samples_array, bool backfill_flag = false);

  /// generates the desired set of parameter samples from within
  /// uncorrelated normal distributions
  void generate_normal_samples(const RealVector& n_means,
    const RealVector& n_std_devs, const RealVector& n_l_bnds,
    const RealVector& n_u_bnds,	RealSymMatrix& corr, int num_samples,
    RealMatrix& samples_array);

  /// generates the desired set of parameter samples from within
  /// uncorrelated uniform distributions
  void generate_uniform_samples(const RealVector& u_l_bnds,
    const RealVector& u_u_bnds, RealSymMatrix& corr, int num_samples,
    RealMatrix& samples_array);

  /// generates integer index samples from within uncorrelated uniform
  /// distributions
  void generate_uniform_index_samples(const IntVector& index_l_bnds,
    const IntVector& index_u_bnds, int num_samples, IntMatrix& index_samples,
    bool backfill_flag = false);

  /// generates unique integer index samples from within uncorrelated
  /// uniform distributions (more expensive than non-unique case)
  //void generate_unique_index_samples(const IntVector& index_l_bnds,
  //  const IntVector& index_u_bnds, int num_samples,
  //  std::set<IntArray>& sorted_samples);
  void generate_unique_index_samples(const IntVector& index_l_bnds,
    const IntVector& index_u_bnds, int num_samples, IntMatrix& sorted_samples);

  /// Similar to generate_samples but this function ensures that all discrete 
  /// samples are unique
  void generate_unique_samples(const std::vector<RandomVariables>& random_vars,
			       const RealSymMatrix& correlations,
			       int num_samples, RealMatrix& samples,
			       RealMatrix& sample_ranks);

private:

  //
  //- Heading: Convenience functions
  //

  /// checks whether LHS is enabled in the build and aborts if not
  static void abort_if_no_lhs();

  /// checks the return codes from LHS routines and aborts if an
  /// error is returned
  void check_error(int err_code, const char* err_source) const;

  //
  //- Heading: Data
  //

  /// type of sampling: random, lhs, incremental_lhs, or incremental_random
  String sampleType;

  /// variable labels passed to LHS; required for specification of correlations
  StringArray lhsNames;

  /// the current random number seed
  int randomSeed;

  /// mode of sample ranks I/O: IGNORE_RANKS, SET_RANKS, GET_RANKS, or
  /// SET_GET_RANKS
  short sampleRanksMode;

  /// flag for generating LHS report output
  bool reportFlag;

  /// for honoring advance_seed_sequence() calls
  short allowSeedAdvance; // bit 1 = first-time flag
		          // bit 2 = allow repeated seed update
};


inline LHSDriver::LHSDriver() : sampleType("lhs"), randomSeed(0),
  sampleRanksMode(IGNORE_RANKS), reportFlag(true), allowSeedAdvance(1)
{ abort_if_no_lhs(); }


inline LHSDriver::~LHSDriver()
{ }


inline void LHSDriver::
initialize(const String& sample_type, short sample_ranks_mode, bool reports)
{
  sampleType      = sample_type;
  sampleRanksMode = sample_ranks_mode;
  reportFlag      = reports;
}


inline LHSDriver::LHSDriver(const String& sample_type,
			    short sample_ranks_mode, bool reports) :
  allowSeedAdvance(1)
{
  abort_if_no_lhs();
  initialize(sample_type, sample_ranks_mode, reports);
}


inline int LHSDriver::seed() const
{ return randomSeed; }


/** It would be preferable to call srand() only once and then call rand()
    for each LHS execution (the intended usage model), but possible
    interaction with other uses of rand() in other contexts is a concern.
    E.g., an srand(clock()) executed elsewhere would ruin the
    repeatability of the LHSDriver seed sequence.  The only way to insure
    isolation is to reseed each time.  Any induced correlation should be
    inconsequential for the intended use. */
inline void LHSDriver::advance_seed_sequence()
{
  if (allowSeedAdvance & 2) { // repeated seed updates allowed
    std::srand(randomSeed);
    seed(1 + std::rand()); // from 1 to RANDMAX+1
  }
}


inline void LHSDriver::
uncertain_subset(const std::vector<RandomVariables>& random_vars,
		 BitArray& subset)
{
  size_t i, num_rv = random_vars.size();
  subset.resize(num_rv, true); // init bits to true
  for (i=0; i<num_rv; ++i)
    // deactivate complement of uncertain vars
    switch (random_vars[i].type()) {
    case CONTINUOUS_RANGE:    case DISCRETE_RANGE: case DISCRETE_SET_INT:
    case DISCRETE_SET_STRING: case DISCRETE_SET_REAL:
      subset.reset(i); break;
    }
}


inline void LHSDriver::
aleatory_uncertain_subset(const std::vector<RandomVariables>& random_vars,
			  BitArray& subset)
{
  size_t i, num_rv = random_vars.size();
  subset.resize(num_rv, true); // init bits to true
  for (i=0; i<num_rv; ++i)
    // deactivate complement of aleatory uncertain vars
    switch (random_vars[i].type()) {
    case CONTINUOUS_RANGE:    case DISCRETE_RANGE: case DISCRETE_SET_INT:
    case DISCRETE_SET_STRING: case DISCRETE_SET_REAL:
    case CONTINUOUS_INTERVAL_UNCERTAIN: case DISCRETE_INTERVAL_UNCERTAIN:
    case DISCRETE_UNCERTAIN_SET_INT:    case DISCRETE_UNCERTAIN_SET_STRING:
    case DISCRETE_UNCERTAIN_SET_REAL:
      subset.reset(i); break;
    }
}


inline void LHSDriver::
epistemic_uncertain_subset(const std::vector<RandomVariables>& random_vars,
			   BitArray& subset)
{
  size_t i, num_rv = random_vars.size();
  subset.resize(num_rv, false); // init bits to false
  for (i=0; i<num_rv; ++i)
    // activate epistemic uncertain vars
    switch (random_vars[i].type()) {
    case CONTINUOUS_INTERVAL_UNCERTAIN: case DISCRETE_INTERVAL_UNCERTAIN:
    case DISCRETE_UNCERTAIN_SET_INT:    case DISCRETE_UNCERTAIN_SET_STRING:
    case DISCRETE_UNCERTAIN_SET_REAL:
      subset.set(i); break;
    }
}


inline void LHSDriver::
generate_uncertain_samples(const std::vector<RandomVariables>& random_vars,
			   const RealSymMatrix& correlations, int num_samples,
			   RealMatrix& samples_array, bool backfill_flag)
{
  if (sampleRanksMode) {
    PCerr << "Error: LHSDriver::generate_uncertain_samples() does not support "
	  << "sample rank input/output." << std::endl;
    abort_handler(-1);
  }
  RealMatrix ranks;
  BitArray active_rv;            uncertain_subset(random_vars, active_rv);
  BitArray active_corr; aleatory_uncertain_subset(random_vars, active_corr);
  if (backfill_flag)
    generate_unique_samples(random_vars, correlations, num_samples,
			    samples_array, ranks, active_rv, active_corr);
  else
    generate_samples(random_vars, correlations, num_samples,
		     samples_array, ranks, active_rv, active_corr);
}


inline void LHSDriver::
generate_aleatory_samples(const std::vector<RandomVariables>& random_vars,
			  const RealSymMatrix& correlations, int num_samples,
			  RealMatrix& samples_array, bool backfill_flag)
{
  if (sampleRanksMode) {
    PCerr << "Error: generate_aleatory_samples() does not support sample rank "
	  << "input/output." << std::endl;
    abort_handler(-1);
  }
  RealMatrix ranks;
  BitArray active_rv; aleatory_uncertain_subset(random_vars, active_rv);
  if (backfill_flag)
    generate_unique_samples(random_vars, correlations, num_samples,
			    samples_array, ranks, active_rv, active_rv);
  else
    generate_samples(random_vars, correlations, num_samples,
		     samples_array, ranks, active_rv, active_rv);
}


inline void LHSDriver::
generate_epistemic_samples(const std::vector<RandomVariables>& random_vars,
			   const RealSymMatrix& correlations, int num_samples,
			   RealMatrix& samples_array, bool backfill_flag)
{
  if (sampleRanksMode) {
    PCerr << "Error: generate_epistemic_samples() does not support sample rank "
	  << "input/output." << std::endl;
    abort_handler(-1);
  }
  RealMatrix ranks;
  BitArray active_rv;  epistemic_uncertain_subset(random_vars, active_rv);
  BitArray active_corr; aleatory_uncertain_subset(random_vars, active_corr);
  if (backfill_flag)
    generate_unique_samples(random_vars, correlations, num_samples,
			    samples_array, ranks, active_rv, active_corr);
  else
    generate_samples(random_vars, correlations, num_samples,
		     samples_array, ranks, active_rv, active_corr);
}


inline void LHSDriver::
generate_normal_samples(const RealVector& n_means, const RealVector& n_std_devs,
			const RealVector& n_l_bnds, const RealVector& n_u_bnds,
			RealSymMatrix& corr, int num_samples,
			RealMatrix& samples_array)
{
  if (sampleRanksMode) {
    PCerr << "Error: generate_normal_samples() does not support sample rank "
  	  << "input/output." << std::endl;
    abort_handler(-1);
  }
  size_t i, num_rv = n_means.length();
  std::vector<RandomVariables> random_vars(num_rv);
  bool l_bnd = !n_l_bnds.empty(), u_bnd = !n_u_bnds.empty();
  for (i=0; i<num_rv; ++i) {
    RandomVariable& rv_i = randomVars[i];
    rv_i = RandomVariable(NORMAL);
    rv_i.push_parameter(N_MEAN,    n_means[i]);
    rv_i.push_parameter(N_STD_DEV, n_std_devs[i]);
    if (l_bnd) rv_i.push_parameter(N_LWR_BND, n_l_bnds[i]);
    if (u_bnd) rv_i.push_parameter(N_UPR_BND, n_u_bnds[i]);
  }
  RealMatrix ranks;
  generate_samples(random_vars, corr, num_samples, samples_array, ranks);
}


inline void LHSDriver::
generate_uniform_samples(const RealVector& u_l_bnds, const RealVector& u_u_bnds,
			 RealSymMatrix& corr, int num_samples,
			 RealMatrix& samples_array)
{
  if (sampleRanksMode) {
    PCerr << "Error: generate_uniform_samples() does not support sample rank "
	  << "input/output." << std::endl;
    abort_handler(-1);
  }
  size_t i, num_rv = u_l_bnds.length();
  std::vector<RandomVariables> random_vars(num_rv);
  for (i=0; i<num_rv; ++i) {
    RandomVariable& rv_i = randomVars[i];
    rv_i = RandomVariable(UNIFORM);
    rv_i.push_parameter(U_LWR_BND, u_l_bnds[i]);
    rv_i.push_parameter(U_UPR_BND, u_u_bnds[i]);
  }
  RealMatrix ranks;
  generate_samples(random_vars, corr, num_samples, samples_array, ranks);
}


inline void LHSDriver::
generate_uniform_index_samples(const IntVector& index_l_bnds,
			       const IntVector& index_u_bnds, int num_samples,
			       IntMatrix& index_samples, bool backfill_flag)
{
  if (sampleRanksMode) {
    PCerr << "Error: generate_uniform_index_samples() does not support sample "
	  << "rank input/output." << std::endl;
    abort_handler(-1);
  }
  // For    uniform probability, model as discrete range (this fn).
  // For nonuniform probability, model as discrete uncertain set integer.
  size_t i, num_rv = index_l_bnds.length();
  std::vector<RandomVariables> random_vars(num_rv);
  for (i=0; i<num_rv; ++i) {
    RandomVariable& rv_i = randomVars[i];
    rv_i = RandomVariable(DISCRETE_RANGE);
    rv_i.push_parameter(DR_LWR_BND, index_l_bnds[i]);
    rv_i.push_parameter(DR_UPR_BND, index_u_bnds[i]);
  }
  RealMatrix ranks, samples_rm;  RealSymMatrix corr; // assume uncorrelated
  if (backfill_flag)
    generate_unique_samples(random_vars, corr, num_samples, samples_rm, ranks);
  else
    generate_samples(random_vars, corr, num_samples, samples_rm, ranks);
  copy_data(samples_rm, index_samples);
}


/** Helper function to create labels with appended tags that Fortran
    will see as character*16 values which are NOT null-terminated. */
inline void LHSDriver::f77name16(const char* name, size_t index, String& label)
{
  label = name + boost::lexical_cast<String>(index+1);
  label.resize(16, ' '); // NOTE: no NULL terminator
}


/** Helper function to create labels that Fortran will see as character*32
    values which are NOT null-terminated. */
inline void LHSDriver::f77name32(const char* name, String& label)
{
  label = name;
  label.resize(32, ' '); // NOTE: no NULL terminator
}


template <typename T>
void LHSDriver::check_range(T l_bnd, T u_bnd, bool allow_equal) const
{
  if (l_bnd > u_bnd) {
    PCerr << "\nError: lower bound exceeds upper bound in Pecos::LHSDriver."
	  << std::endl;
    abort_handler(-1);
  }
  else if (!allow_equal && l_bnd == u_bnd) {
    PCerr << "\nError: Pecos::LHSDriver requires non-zero range between lower "
	  << "and upper bounds." << std::endl;
    abort_handler(-1);
  }
}


template <typename T>
void LHSDriver::check_finite<T>(T l_bnd, T u_bnd) const
{
  if (std::numeric_limits<T>::has_infinity) { // floating point types
    typename T T_inf = std::numeric_limits<T>::infinity();
    if (l_bnd <= -T_inf || u_bnd >= T_inf) {
      PCerr << "\nError: Pecos::LHSDriver requires finite bounds to sample a "
	    << "continuous range." << std::endl;
      abort_handler(-1);
    }
  }
  else if (l_bnd <= std::numeric_limits<T>::min() ||
	   u_bnd >= std::numeric_limits<T>::max()) { // INT_{MIN,MAX}
    PCerr << "\nError: Pecos::LHSDriver requires finite bounds to sample a "
	  << "discrete range." << std::endl;
    abort_handler(-1);
  }
}


inline void LHSDriver::
check_error(int err_code, const char* err_source, const char* err_case) const
{
  if (err_code) {
    PCerr << "Error: code " << err_code << " in LHSDriver";
    if (err_source != NULL) PCerr << " returned from " << err_source;
    if (err_case   != NULL) PCerr << " for case "      << err_case;
    PCerr << "." << std::endl;
    abort_handler(-1);
  }
}


inline void LHSDriver::
bins_to_udist_params(const RealRealMap& h_bin_prs,
		     RealArray& x_val, RealArray& y_val)
{
  // histogram bins: pairs are defined from an abscissa in the first field
  // and a count (not a density) in the second field.  This distinction is
  // important for unequal bin widths.

  size_t i, num_params = h_bin_prs.size(), end = num_params - 1;
  RRMCIter cit = h_bin_prs.begin();

  // Assume already normalized with sum = 1
  //Real sum = 0.;
  //RRMCIter end = --h_bin_prs.end(); // last y from DAKOTA must be zero
  //for (; cit!=end; ++cit)
  //  sum += cit->second;

  // LHS requires accumulation of CDF with first y at 0 and last y at 1
  x_val.resize(num_params);  y_val.resize(num_params);
  y_val[0] = 0.;
  for (i=0; i<end; ++i, ++cit) {
    x_val[i]   = cit->first;
    y_val[i+1] = y_val[i] + cit->second/* /sum */;
  }
  x_val[end] = cit->first; // last prob value (cit->second) must be zero
}


inline void LHSDriver::
int_range_to_udist_params(int l_bnd,        int u_bnd,
			  RealArray& x_val, RealArray& y_val)
{
  // supports either discrete integer range or range of set indices

  int i, num_params = ub_i - lb_i + 1;
  x_val.resize(num_params);  y_val.assign(num_params, 1.);
  for (i=0; i<num_params; ++i)
    x_val[i] = (Real)(l_bnd + i);
}


inlinetemplate <typename T>
void LHSDriver::
set_to_udist_params(const std::set<T>& values,
		    RealArray& x_val, RealArray& y_val)
{
  int i, num_params = values.size();
  x_val.resize(num_params);  y_val.assign(num_params, 1.);
  typename std::set<T>::const_iterator cit;
  for (cit=values.begin(), i=0; cit!=values.end(); ++cit, ++i)
    x_val[i] = (Real)(*cit); // value
}


template <typename T>
void LHSDriver::
map_to_udist_params(const std::map<T, Real>& vals_probs,
		    RealArray& x_val, RealArray& y_val)
{
  int i, num_params = vals_probs.size();
  x_val.resize(num_params);  y_val.resize(num_params);
  typename std::set<T>::const_iterator cit;
  for (cit=vals_probs.begin(), i=0; cit!=vals_probs.end(); ++cit, ++i) {
    x_val[i] = (Real)cit->first;  // value
    y_val[i] =       cit->second; // probability
  }
}


template <typename T>
void LHSDriver::
map_indices_to_udist_params(const std::map<T, Real>& vals_probs,
			    RealArray& x_val, RealArray& y_val)
{
  int i, num_params = vals_probs.size();
  x_val.resize(num_params);  y_val.resize(num_params);
  typename std::set<T>::const_iterator cit;
  for (cit=vals_probs.begin(), i=0; cit!=vals_probs.end(); ++cit, ++i) {
    x_val[i] = (Real)i;     // index rather than value
    y_val[i] = cit->second; // probability
  }
}


void LHSDriver::
intervals_to_udist_params(const RealRealPairRealMap& ci_bpa,
			  RealArray& x_val, RealArray& y_val)
{
  // x_sort_unique is a set with ALL of the interval bounds for this variable
  // in increasing order and unique.  For example, if there are 2 intervals
  // for a variable, and the bounds are (1,4) and (3,6), x_sorted will be
  // (1, 3, 4, 6).  If the intervals are contiguous, e.g. one interval is
  // (1,3) and the next is (3,5), x_sort_unique is (1,3,5).
  RRPRMCIter cit;  RealSet x_sort_unique;
  for (cit=ci_bpa.begin(); cit!=ci_bpa.end(); ++cit) {
    const RealRealPair& bounds = cit->first;
    x_sort_unique.insert(bounds.first);
    x_sort_unique.insert(bounds.second);
  }
  // convert sorted RealSet to x_val
  num_params = x_sort_unique.size();
  x_val.reshape(num_params);  y_val.reshape(num_params);
  RSIter it = x_sort_unique.begin();
  for (j=0; j<num_params; ++j, ++it)
    x_val[j] = *it;

  // Calculate the probability densities, and account for the cases where
  // there are intervals that are overlapping.  This section of code goes
  // through the original intervals and see where they fall relative to the
  // new, sorted intervals for the density calculation.
  RealVector prob_dens(num_params); // initialize to 0.
  for (cit=ci_bpa.begin(); cit!=ci_bpa.end(); ++cit) {
    const RealRealPair& bounds = cit->first;
    Real l_bnd = bounds.first, u_bnd = bounds.second;
    Real ci_density = cit->second / (u_bnd - l_bnd);
    int cum_int_index = 0;
    while (l_bnd > x_val[cum_int_index])
      ++cum_int_index;
    ++cum_int_index;
    while (cum_int_index < num_params && x_val[cum_int_index] <= u_bnd)
      { prob_dens[cum_int_index] += ci_density; ++cum_int_index; }
  }

  // put the densities in a cumulative format necessary for LHS histograms.
  // Note that x_val and y_val are defined as Real* for input to f77.
  y_val[0] = 0.;
  for (j=1; j<num_params; ++j)
    y_val[j] = (prob_dens[j] > 0.0) ?
      y_val[j-1] + prob_dens[j] * (x_val[j] - x_val[j-1]) :
      y_val[j-1] + 0.0001; // handle case where there is a gap
  // normalize if necessary
  if (y_val[num_params-1] != 1.) {
    Real y_total = y_val[num_params-1];
    for (j=1; j<num_params; ++j)
      y_val[j] /= y_total;
  }
#ifdef DEBUG
  for (j=0; j<num_params; ++j)
    PCout << "ciuv: x_val[" << j << "] is " << x_val[j]
	  << " y_val[" << j << "] is " << y_val[j] << '\n';
#endif // DEBUG
}


void LHSDriver::
intervals_to_udist_params(const IntIntPairRealMap& di_bpa,
			  RealArray& x_val, RealArray& y_val)
{
  // x_sort_unique contains ALL of the unique integer values for this
  // x_sort_unique contains ALL of the unique integer values for this
  // discrete interval variable in increasing order.  For example, if
  // there are 3 intervals for a variable and the bounds are (1,4),
  // (3,6), and (9,10), x_sorted will be (1,2,3,4,5,6,9,10).
  IIPRMCIter cit; IntSet x_sort_unique;
  for (cit=di_bpa.begin(); cit!=di_bpa.end(); ++cit) {
    const RealRealPair& bounds = cit->first;
    int val, u_bnd = bounds.second;
    for (val=bounds.first; val<=u_bnd; ++val)
      x_sort_unique.insert(val);
  }
  // copy sorted IntSet to x_val
  num_params = x_sort_unique.size();
  x_val.reshape(num_params);  y_val.reshape(num_params);
  ISIter it = x_sort_unique.begin();
  for (j=0; j<num_params; ++j, ++it)
    x_val[j] = *it;

  // Calculate probability densities and account for overlapping intervals.
  // Loop over the original intervals and see where they fall relative to
  // the new, sorted intervals for the density calculation.
  for (j=0; j<num_params; ++j) y_val[j] = 0.;
  int l_bnd, u_bnd; size_t index;
  for (cit=di_bpa.begin(); cit!=di_bpa.end(); ++cit) {
    const RealRealPair& bounds = cit->first;
    int val, l_bnd = bounds.first, u_bnd = bounds.second;
    Real di_density = cit->second / (u_bnd - l_bnd + 1); // prob/#integers
    it = x_sort_unique.find(l_bnd);
    if (it == x_sort_unique.end()) {
      PCerr << "Error: lower bound not found in sorted set within LHSDriver "
	    << "mapping of discrete interval uncertain variable."<< std::endl;
      abort_handler(-1);
    }
    index = std::distance(x_sort_unique.begin(), it);
    for (val=l_bnd; val<=u_bnd; ++val, ++index)
      y_val[index] += di_density;
  }
#ifdef DEBUG
  for (j=0; j<num_params; ++j)
    PCout << "diuv: x_val[" << j << "] is " << x_val[j]
	  << " y_val[" << j << "] is " << y_val[j] << '\n';
#endif // DEBUG
}


inline void LHSDriver::
lhs_dist_register(const char* var_name, const char* dist_name, size_t rv,
		  const RealArray& dist_params) const
{
  String dist_string;                 f77name32(dist_name,   dist_string);
  String& var_string = lhsNames[rv];  f77name16(var_name, rv, var_string);
  int num_params = dist_params.size(), err_code = 0, ptval_flag = 0, // inputs
      dist_num, pv_num; // outputs (not used)
  Real ptval = 0.;

  LHS_DIST2_FC(var_string.data(), ptval_flag, ptval, dist_string.data(),
	       &dist_params[0], num_params, err_code, dist_num, pv_num);
  check_error(err_code, "lhs_dist()", var_string.data());
}


inline void LHSDriver::
lhs_udist_register(const char* var_name, const char* dist_name, size_t rv,
		   const RealArray& x_val, const RealArray& y_val) const
{
  String dist_string;                 f77name32(dist_name,   dist_string);
  String& var_string = lhsNames[rv];  f77name16(var_name, rv, var_string);
  int num_params = std::min(x_val.size(), y_val.size()), err_code = 0,
    ptval_flag = 0, dist_num, pv_num;
  Real ptval = 0.;

  LHS_UDIST2_FC(var_string.data(), ptval_flag, ptval, dist_string.data(),
		num_params, &x_val[0], &y_val[0], err_code, dist_num, pv_num);
  check_error(err_code, "lhs_udist()", var_string.data());
}


template <typename T>
void LHSDriver::lhs_const_register(const char* var_name, size_t rv, T val)
{
  String& var_string = lhsNames[rv];  f77name16(var_name, rv, var_string);
  int err_code = 0, pv_num;           Real pt_val = (Real)val;
  LHS_CONST2_FC(var_string.data(), pt_val, err_code, pv_num);
  check_error(err_code, "lhs_const()", var_string.data());
}

} // namespace Pecos

#endif

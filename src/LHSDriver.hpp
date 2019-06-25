/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef LHS_DRIVER_H
#define LHS_DRIVER_H

#include "pecos_data_types.hpp"
#include "RandomVariable.hpp"


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

  /// generates a set of parameter samples from RandomVariable distributions
  void generate_samples(const std::vector<RandomVariable>& random_vars,
			const RealSymMatrix& corr, int num_samples,
			RealMatrix& samples, RealMatrix& sample_ranks,
			const BitArray& active_vars = BitArray(),
			const BitArray& active_corr = BitArray());
  /// Similar to generate_samples, but this function seeks uniqueness in
  /// discrete samples through iterative replacement of duplicates
  void generate_unique_samples(const std::vector<RandomVariable>& random_vars,
			       const RealSymMatrix& corr, int num_samples,
			       RealMatrix& samples, RealMatrix& sample_ranks,
			       const BitArray& active_vars = BitArray(),
			       const BitArray& active_corr = BitArray());

  /// generates a set of parameter samples from RandomVariable distributions
  /// for the uncertain variable subset
  void generate_uncertain_samples(
    const std::vector<RandomVariable>& random_vars, const RealSymMatrix& corr,
    int num_samples, RealMatrix& samples_array, bool backfill_flag = false);
  /// generates a set of parameter samples from RandomVariable distributions
  /// for the aleatory uncertain variable subset
  void generate_aleatory_samples(
    const std::vector<RandomVariable>& random_vars,
    const RealSymMatrix& corr, int num_samples, RealMatrix& samples_array,
    bool backfill_flag = false);
  /// generates a set of parameter samples from RandomVariable distributions
  /// for the epistemic uncertain variable subset
  void generate_epistemic_samples(
    const std::vector<RandomVariable>& random_vars, const RealSymMatrix& corr,
    int num_samples, RealMatrix& samples_array, bool backfill_flag = false);

  /// generates a set of parameter samples for a set of (correlated)
  /// normal distributions
  void generate_normal_samples(const RealVector& n_means,
    const RealVector& n_std_devs, const RealVector& n_l_bnds,
    const RealVector& n_u_bnds,	RealSymMatrix& corr, int num_samples,
    RealMatrix& samples_array);
  /// generates a set of parameter samples for a set of (correlated)
  /// uniform distributions
  void generate_uniform_samples(const RealVector& u_l_bnds,
    const RealVector& u_u_bnds, RealSymMatrix& corr, int num_samples,
    RealMatrix& samples_array);

  /// generates integer index samples from within uncorrelated uniform
  /// discrete range distributions
  void generate_uniform_index_samples(const IntVector& index_l_bnds,
    const IntVector& index_u_bnds, int num_samples, IntMatrix& index_samples,
    bool backfill_flag = false);

  /// define the subset of design or state variables from RV type within
  /// the range [start_set,start_set+num_set)
  void design_state_subset(const std::vector<RandomVariable>& random_vars,
			   BitArray& subset, size_t start_set,
			   size_t num_set) const;
  /// define the subset of uncertain variables from RV type
  void uncertain_subset(const std::vector<RandomVariable>& random_vars,
			BitArray& subset) const;
  /// define the subset of aleatory uncertain variables from RV type
  void aleatory_uncertain_subset(
    const std::vector<RandomVariable>& random_vars, BitArray& subset) const;
  /// define the subset of epistemic uncertain variables from RV type
  void epistemic_uncertain_subset(
    const std::vector<RandomVariable>& random_vars, BitArray& subset) const;

private:

  //
  //- Heading: Convenience functions
  //

  /// checks whether LHS is enabled in the build and aborts if not
  static void abort_if_no_lhs();

  /// check for valid range for type T
  template <typename T>
  void check_range(T l_bnd, T u_bnd, bool allow_equal) const;
  /// check for finite bounds for type T
  template <typename T>
  void check_finite(T l_bnd, T u_bnd) const;
  /// checks the return codes from LHS routines and aborts with an
  /// error message if an error was returned
  void check_error(int err_code, const char* err_source = NULL,
		   const char* err_case = NULL) const;

  /// convert histogram bin pairs to LHS udist x,y inputs
  void bins_to_udist_params(const RealRealMap& h_bin_prs,
			    RealArray& x_val, RealArray& y_val);
  /// convert int range to LHS udist x,y inputs
  void int_range_to_udist_params(int l_bnd,        int u_bnd,
				 RealArray& x_val, RealArray& y_val);
  /// convert set values to LHS udist x,y inputs
  template <typename T>
  void set_to_udist_params(const std::set<T>& values,
			   RealArray& x_val, RealArray& y_val);
  /// convert map pairs (value T to Real prob) to LHS udist x,y inputs
  template <typename T>
  void map_to_udist_params(const std::map<T, Real>& vals_probs,
			   RealArray& x_val, RealArray& y_val);
  /// convert map pairs (string index to Real prob) to LHS udist x,y inputs
  template <typename T>
  void map_indices_to_udist_params(const std::map<T, Real>& vals_probs,
				   RealArray& x_val, RealArray& y_val);
  /// convert epistemic BPAs to LHS udist x,y inputs
  void intervals_to_udist_params(const RealRealPairRealMap& ci_bpa,
				 RealArray& x_val, RealArray& y_val);
  /// convert epistemic BPAs to LHS udist x,y inputs
  void intervals_to_udist_params(const IntIntPairRealMap& di_bpa,
				 RealArray& x_val, RealArray& y_val);

  /// register a random variable with LHS using the DIST format
  void lhs_dist_register(const char* var_name, const char* dist_name,
			 size_t rv, const RealArray& dist_params);
  /// register a random variable with LHS using the UDIST format
  void lhs_udist_register(const char* var_name, const char* dist_name,
			  size_t rv, const RealArray& x_val,
			  const RealArray& y_val);
  /// register a random variable with LHS using the CONST format
  void lhs_const_register(const char* var_name, size_t rv, Real val);

  /// create F77 string label from base name + tag + padding (field = 16 chars)
  void f77name16(const char* name, size_t index, String& label);
  /// create F77 string label from base name + padding (field = 32 chars)
  void f77name32(const char* name, String& label);

  //
  //- Heading: Data
  //

  /// type of sampling: random, lhs, incremental_lhs, or incremental_random
  String sampleType;

  /// variable labels passed to LHS; required for specification of correlations
  StringArray lhsNames;

  /// mode of sample ranks I/O: IGNORE_RANKS, SET_RANKS, GET_RANKS, or
  /// SET_GET_RANKS
  short sampleRanksMode;

  /// flag for generating LHS report output
  bool reportFlag;

  /// the current random number seed
  int randomSeed;
  /// for honoring advance_seed_sequence() calls
  short allowSeedAdvance; // bit 1 = first-time flag
		          // bit 2 = allow repeated seed update
};


inline LHSDriver::LHSDriver():
  sampleType("lhs"), sampleRanksMode(IGNORE_RANKS), reportFlag(true),
  randomSeed(0), allowSeedAdvance(1)
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
design_state_subset(const std::vector<RandomVariable>& random_vars,
		    BitArray& subset, size_t start_set, size_t num_set) const
{
  size_t i, num_rv = random_vars.size(), end_set = start_set + num_set;
  subset.resize(num_rv, false); // init bits to false
  for (i=start_set; i<end_set; ++i)
    // activate design + state vars
    switch (random_vars[i].type()) {
    case CONTINUOUS_RANGE:    case DISCRETE_RANGE: case DISCRETE_SET_INT:
    case DISCRETE_SET_STRING: case DISCRETE_SET_REAL:
      subset.set(i); break;
    }
}


inline void LHSDriver::
uncertain_subset(const std::vector<RandomVariable>& random_vars,
		 BitArray& subset) const
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
aleatory_uncertain_subset(const std::vector<RandomVariable>& random_vars,
			  BitArray& subset) const
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
epistemic_uncertain_subset(const std::vector<RandomVariable>& random_vars,
			   BitArray& subset) const
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
generate_uncertain_samples(const std::vector<RandomVariable>& random_vars,
			   const RealSymMatrix& corr, int num_samples,
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
    generate_unique_samples(random_vars, corr, num_samples, samples_array,
			    ranks, active_rv, active_corr);
  else
    generate_samples(random_vars, corr, num_samples, samples_array, ranks,
		     active_rv, active_corr);
}


inline void LHSDriver::
generate_aleatory_samples(const std::vector<RandomVariable>& random_vars,
			  const RealSymMatrix& corr, int num_samples,
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
    generate_unique_samples(random_vars, corr, num_samples, samples_array,
			    ranks, active_rv, active_rv);
  else
    generate_samples(random_vars, corr, num_samples, samples_array, ranks,
		     active_rv, active_rv);
}


inline void LHSDriver::
generate_epistemic_samples(const std::vector<RandomVariable>& random_vars,
			   const RealSymMatrix& corr, int num_samples,
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
    generate_unique_samples(random_vars, corr, num_samples, samples_array,
			    ranks, active_rv, active_corr);
  else
    generate_samples(random_vars, corr, num_samples,
		     samples_array, ranks, active_rv, active_corr);
}


inline void LHSDriver::
generate_normal_samples(const RealVector& n_means, const RealVector& n_std_devs,
			const RealVector& n_l_bnds, const RealVector& n_u_bnds,
			RealSymMatrix& corr, int num_samples,
			RealMatrix& samples_array)//, bool backfill_flag)
{
  if (sampleRanksMode) {
    PCerr << "Error: generate_normal_samples() does not support sample rank "
  	  << "input/output." << std::endl;
    abort_handler(-1);
  }
  size_t i, num_rv = n_means.length();
  std::vector<RandomVariable> random_vars(num_rv);
  bool l_bnd = !n_l_bnds.empty(), u_bnd = !n_u_bnds.empty();
  for (i=0; i<num_rv; ++i) {
    RandomVariable& rv_i = random_vars[i];
    rv_i = RandomVariable(NORMAL);
    rv_i.push_parameter(N_MEAN,    n_means[i]);
    rv_i.push_parameter(N_STD_DEV, n_std_devs[i]);
    if (l_bnd) rv_i.push_parameter(N_LWR_BND, n_l_bnds[i]);
    if (u_bnd) rv_i.push_parameter(N_UPR_BND, n_u_bnds[i]);
  }
  RealMatrix ranks;
  // All variables are continuous normal --> no need for backfill
  //if (backfill_flag)
  //  generate_unique_samples(random_vars,corr,num_samples,samples_array,ranks);
  //else
  generate_samples(random_vars, corr, num_samples, samples_array, ranks);
}


inline void LHSDriver::
generate_uniform_samples(const RealVector& u_l_bnds, const RealVector& u_u_bnds,
			 RealSymMatrix& corr, int num_samples,
			 RealMatrix& samples_array)//, bool backfill_flag)
{
  if (sampleRanksMode) {
    PCerr << "Error: generate_uniform_samples() does not support sample rank "
	  << "input/output." << std::endl;
    abort_handler(-1);
  }
  size_t i, num_rv = u_l_bnds.length();
  std::vector<RandomVariable> random_vars(num_rv);
  for (i=0; i<num_rv; ++i) {
    RandomVariable& rv_i = random_vars[i];
    rv_i = RandomVariable(UNIFORM);
    rv_i.push_parameter(U_LWR_BND, u_l_bnds[i]);
    rv_i.push_parameter(U_UPR_BND, u_u_bnds[i]);
  }
  RealMatrix ranks;
  // All variables are continuous uniform --> no need for backfill
  //if (backfill_flag)
  //  generate_unique_samples(random_vars,corr,num_samples,samples_array,ranks);
  //else
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
  std::vector<RandomVariable> random_vars(num_rv);
  for (i=0; i<num_rv; ++i) {
    RandomVariable& rv_i = random_vars[i];
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
void LHSDriver::check_finite(T l_bnd, T u_bnd) const
{
  if (std::numeric_limits<T>::has_infinity) { // floating point types
    T T_inf = std::numeric_limits<T>::infinity();
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

  int i, num_params = u_bnd - l_bnd + 1;
  x_val.resize(num_params);  y_val.assign(num_params, 1.);
  for (i=0; i<num_params; ++i)
    x_val[i] = (Real)(l_bnd + i);
}


template <typename T>
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
  typename std::map<T, Real>::const_iterator cit;
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
  typename std::map<T, Real>::const_iterator cit;
  for (cit=vals_probs.begin(), i=0; cit!=vals_probs.end(); ++cit, ++i) {
    x_val[i] = (Real)i;     // index rather than value
    y_val[i] = cit->second; // probability
  }
}

} // namespace Pecos

#endif

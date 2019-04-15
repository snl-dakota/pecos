/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "LHSDriver.hpp"
#include "BoostRNG_Monostate.hpp"
#include "LognormalRandomVariable.hpp"
#include <boost/lexical_cast.hpp>

static const char rcsId[]="@(#) $Id: LHSDriver.cpp 5248 2008-09-05 18:51:52Z wjbohnh $";

//#define DEBUG

#ifdef HAVE_LHS
#ifdef HAVE_CONFIG_H
// Use the classic, autotools Fortran name mangling macros in pecos_config.h
#define LHS_INIT_MEM_FC FC_FUNC_(lhs_init_mem,LHS_INIT_MEM)
#define LHS_PREP_FC     FC_FUNC_(lhs_prep,LHS_PREP)
#define LHS_RUN_FC      FC_FUNC_(lhs_run,LHS_RUN)
#define LHS_CLOSE_FC    FC_FUNC_(lhs_close,LHS_CLOSE)
#define LHS_OPTIONS2_FC FC_FUNC_(lhs_options2,LHS_OPTIONS2)
#define LHS_DIST2_FC    FC_FUNC_(lhs_dist2,LHS_DIST2)
#define LHS_UDIST2_FC   FC_FUNC_(lhs_udist2,LHS_UDIST2)
#define LHS_CONST2_FC   FC_FUNC_(lhs_const2,LHS_CONST2)
#define LHS_CORR2_FC    FC_FUNC_(lhs_corr2,LHS_CORR2)
#define LHS_FILES2_FC   FC_FUNC_(lhs_files2,LHS_FILES2)

#define defaultrnum1       FC_FUNC(defaultrnum1,DEFAULTRNUM1)
#define defaultrnum2       FC_FUNC(defaultrnum2,DEFAULTRNUM2)

#else
// Use the CMake-generated PREFIXED, fortran name mangling macros (no warnings)
#define LHS_INIT_MEM_FC LHS_GLOBAL_(lhs_init_mem,LHS_INIT_MEM)
#define LHS_PREP_FC     LHS_GLOBAL_(lhs_prep,LHS_PREP)
#define LHS_RUN_FC      LHS_GLOBAL_(lhs_run,LHS_RUN)
#define LHS_CLOSE_FC    LHS_GLOBAL_(lhs_close,LHS_CLOSE)
//#define LHS_OPTIONS2_FC LHS_GLOBAL_(lhs_options2,LHS_OPTIONS2)
//#define LHS_DIST2_FC    LHS_GLOBAL_(lhs_dist2,LHS_DIST2)
//#define LHS_UDIST2_FC   LHS_GLOBAL_(lhs_udist2,LHS_UDIST2)
//#define LHS_CONST2_FC   LHS_GLOBAL_(lhs_const2,LHS_CONST2)
//#define LHS_CORR2_FC    LHS_GLOBAL_(lhs_corr2,LHS_CORR2)
//#define LHS_FILES2_FC   LHS_GLOBAL_(lhs_files2,LHS_FILES2)
#define defaultrnum1    LHS_GLOBAL(defaultrnum1,DEFAULTRNUM1)
#define defaultrnum2    LHS_GLOBAL(defaultrnum2,DEFAULTRNUM2)

// BMA (20160315): Changed to use Fortran 2003 ISO C bindings.
// The Fortran symbol will be lowercase with same name as if in C
#define  LHS_OPTIONS2_FC lhs_options2
#define  LHS_DIST2_FC lhs_dist2
#define  LHS_UDIST2_FC lhs_udist2
#define  LHS_CONST2_FC lhs_const2
#define  LHS_CORR2_FC lhs_corr2
#define  LHS_FILES2_FC lhs_files2


#endif // HAVE_CONFIG_H

extern "C" {

// for these functions, a straight interface to F90 can be used
void LHS_INIT_MEM_FC( int& obs, int& seed, int& max_obs, int& max_samp_size,
		      int& max_var, int& max_interval, int& max_corr,
		      int& max_table, int& print_level, int& output_width,
		      int& ierror );

void LHS_PREP_FC( int& ierror, int& num_names, int& num_vars );

void LHS_RUN_FC( int& max_var, int& max_obs, int& max_names, int& ierror,
		 char* dist_names, int* name_order, Pecos::Real* ptvals,
		 int& num_names, Pecos::Real* sample_matrix, int& num_vars,
                 Pecos::Real* rank_matrix, int& rank_flag );

void LHS_CLOSE_FC( int& ierror );

// for these functions, bridge-function interfaces to F90 are required
void LHS_OPTIONS2_FC( int& num_replications, int& ptval_option,
		      const char* sampling_options, int& ierror );

void LHS_DIST2_FC( const char* label, int& ptval_flag, Pecos::Real& ptval,
		   const char* dist_type, Pecos::Real* dist_params,
		   int& num_params, int& ierror, int& dist_id, int& ptval_id );

void LHS_UDIST2_FC( const char* label, int& ptval_flag, Pecos::Real& ptval,
		    const char* dist_type, int& num_pts, Pecos::Real* x,
		    Pecos::Real* y, int& ierror, int& dist_id, int& ptval_id );

void LHS_CONST2_FC( const char* label, Pecos::Real& ptval, int& ierror,
		    int& ptval_id );

void LHS_CORR2_FC( char* label1, char* label2, Pecos::Real& corr, int& ierror );

void LHS_FILES2_FC( const char* lhsout, const char* lhsmsg, const char* lhstitl,
		    const char* lhsopts, int& ierror );

//void lhs_run2( int* max_var, int* max_obs, int* max_names, int* ierror,
//               char* dist_names, int* name_order, double* ptvals,
//               int* num_names, double* sample_matrix, int* num_vars );

Pecos::Real defaultrnum1( void );
Pecos::Real defaultrnum2( void );

//Pecos::Real rnumlhs10( void );
//Pecos::Real rnumlhs20( void );

}
#endif // HAVE_LHS


namespace Pecos {

void LHSDriver::seed(int seed)
{
  randomSeed = seed;
  // The Boost RNG is not set by LHS_INIT_MEM, so must be done here.
  if (BoostRNG_Monostate::randomNum == BoostRNG_Monostate::mt19937)
    BoostRNG_Monostate::seed(seed);
  // This would be redundant since the f77 ISeed is set in LHS_INIT_MEM:
  //else
  // lhs_setseed(&seed);
}


void LHSDriver::rng(String unif_gen)
{

  // check the environment (once only) for an RNG preference and cache it
  static bool first_entry = true;
  static const char *env_unifgen;
  if (first_entry) {
    env_unifgen = std::getenv("DAKOTA_LHS_UNIFGEN");
    first_entry = false;
  }

  // the environment overrides the passed rng specification
  if (env_unifgen) {
    unif_gen = env_unifgen;
    if (unif_gen != "rnum2" && unif_gen != "mt19937") {
      PCerr << "Error: LHSDriver::rng() expected $DAKOTA_LHS_UNIFGEN to be "
	    << "\"rnum2\" or \"mt19937\", not \"" << env_unifgen << "\".\n"
	    << std::endl;
      abort_handler(-1);
    }
  }

  // now point the monostate RNG to the desired generator function
  if (unif_gen == "mt19937" || unif_gen.empty()) {
    BoostRNG_Monostate::randomNum  = BoostRNG_Monostate::mt19937;
    BoostRNG_Monostate::randomNum2 = BoostRNG_Monostate::mt19937;
    allowSeedAdvance &= ~2; // drop 2 bit: disallow repeated seed update
  }
  else if (unif_gen == "rnum2") {
#ifdef HAVE_LHS
    BoostRNG_Monostate::randomNum  = (Rfunc)defaultrnum1;
    BoostRNG_Monostate::randomNum2 = (Rfunc)defaultrnum2;
    allowSeedAdvance |= 2;  // add 2 bit: allow repeated seed update
#else
    PCerr << "Error: LHSDriver::rng() Use of rnum2 for RNG selection is NOT "
	  << "supported in current (without-lhs) configuration" << std::endl;
    abort_handler(-1);
#endif
  }
  else {
    PCerr << "Error: LHSDriver::rng() expected string to be \"rnum2\" or "
	  << "\"mt19937\", not \"" << unif_gen << "\".\n" << std::endl;
    abort_handler(-1);
  }
}


void LHSDriver::abort_if_no_lhs()
{
#ifndef HAVE_LHS
  PCerr << "Error: LHSDriver not available as PECOS was configured without LHS."
        << std::endl;
  abort_handler(-1);
#endif
}


//String LHSDriver::rng()
//{
//  if (BoostRNG_Monostate::randomNum == BoostRNG_Monostate::random_num1)
//    return "mt19937";
//  if (BoostRNG_Monostate::randomNum == (Rfunc)rnumlhs10)
//    return "rnum2";
//  else
//    return "";
//}


/** While it would be desirable in some cases to carve this function
    into smaller parts and allow multiple invocations of LHS_RUN
    following a single initialization of types and arrays, the LHS
    code does not currently allow this: it will return an error if
    LHS_INIT/LHS_INIT_MEM, at least one distribution call (i.e.,
    LHS_DIST, LHS_UDIST or LHS_SDIST), and LHS_PREP are not called
    prior to each invocation of LHS_RUN.  Since LHS_INIT/LHS_INIT_MEM
    require input of a seed, the approach to computing multiple
    distinct sample sets must employ advance_seed_sequence() to
    re-seed multiple generate_samples() calls, rather than continuing
    an existing random number sequence. */
void LHSDriver::
generate_samples(const std::vector<RandomVariables>& random_vars,
		 const RealSymMatrix& correlations, int num_samples,
		 RealMatrix& samples, RealMatrix& sample_ranks,
		 const BitArray& active_vars, const BitArray& active_corr)
{
#ifdef HAVE_LHS
  // generate samples within user-specified parameter distributions

  // error check on program parameters
  if (!num_samples) {
    PCerr << "\nError: number of samples in LHSDriver::generate_samples() must "
	  << "be nonzero." << std::endl;
    abort_handler(-1);
  }

  // active_vars identifies the active subset of random_vars for which we will
  // generate samples; it will vary based on sampling context.
  size_t i, num_rv = random_vars.size(), num_active_rv;
  bool subset_rv = false;
  if (active_vars.empty())
    num_active_rv = num_rv;
  else {
    num_active_rv = active_vars.count();
    if (num_active_rv < num_rv) subset_rv = true;
  }
  lhsNames.resize(num_active_rv);

  // active_corr identifies the subset of random_vars for which we specify
  // correlations; it need not vary based on sampling context, potentially
  // indicating a static subset that admits correlations in general (e.g.,
  // cont/disc aleatory RV in the current Dakota XML spec).
  // > Both active subsets are relative to the full RV vector (i.e., not
  //   nested subsets) such that they overlay when specifying correlations
  //   to LHS, e.g., uncertain vars are (currently) active for sampling and
  //   correlations are (always) bound to aleatory uncertain vars, so both
  //   bits must be on to call LHS_CORR2.
  // > avoid resizing the correlation matrix (inflating with 1's on diagonal)
  //   for active sampling context: the size of the corr matrix is defined by
  //   the (static) RV subset that admits correlations (correlations.numRows()
  //   == active_corr.count()).  Since this is independent of the active
  //   sampling context, a subset of corr matrix could be passed to LHS_CORR2.
  bool correlation_flag = !correlations.empty(), subset_corr = false;
  int  num_corr = 0, num_active_corr = 0;
  if (correlation_flag) {
    if (active_corr.empty())
      num_corr = correlations.numRows();
    else {
      for (i=0; i<num_rv; ++i)
	if (active_corr[i])
	  { ++num_active_corr; if (active_vars[i]) ++num_corr; }
      if (num_active_corr < num_rv) subset_corr = true;
    }
  }

  int max_corr_size = (num_corr > 1) ? num_corr * (num_corr - 1) / 2 : -1,
      max_var = num_active_rv, max_samp_size = max_var * num_samples,
      max_interval = -1, max_table = -1, err_code = 0, print_level = 0,
      output_width = 1;
  // randomSeed passed below propagates to ISeed in the f77 rnum2, but does
  // not propagate to Boost RNGs (LHSDriver::seed() must be used for that).
  LHS_INIT_MEM_FC(num_samples, randomSeed, num_samples, max_samp_size, max_var,
		  max_interval, max_corr_size, max_table, print_level,
		  output_width, err_code);
  check_error(err_code, "lhs_init_mem");

  // set sample type to either LHS (default) or random Monte Carlo (optional)
  bool call_lhs_option = false;
  String option_string("              ");
  if (sampleType == "random" || sampleType == "incremental_random")
    { option_string = "RANDOM SAMPLE "; call_lhs_option = true; }
  // set mixing option to either restricted pairing (default) or random pairing
  // (optional).  For enforcing user-specified correlation, restricted pairing
  // is required.  And for uncorrelated variables, restricted pairing results
  // in lower correlation values than random pairing.  For these reasons, the
  // random pairing option is not currently active, although a specification
  // option for it could be added in the future if a use arises.
  bool random_pairing_flag = false; // this option hard-wired off for now
  if (!correlation_flag && random_pairing_flag)
    { option_string += "RANDOM PAIRING"; call_lhs_option = true; }
  // else // use default of restricted pairing
  option_string.resize(32, ' ');
  if (call_lhs_option) {
    // Don't null-terminate the string since the '\0' is not used in Fortran
    int num_replic = 1, ptval_option = 1;
    LHS_OPTIONS2_FC(num_replic, ptval_option, option_string.data(), err_code);
    check_error(err_code, "lhs_options");
  }

  // Create files showing distributions and associated statistics.  Avoid
  // issues with null-terminated strings from C++ (which mess up the Fortran
  // output) by using std::string::data().
  String output_string("LHS_samples.out");
  output_string.resize(32, ' ');
  String message_string("LHS_distributions.out");
  message_string.resize(32, ' ');
  String title_string("Pecos::LHSDriver");
  title_string.resize(32, ' ');
  // From the LHS manual (p. 100): LHSRPTS is used to specify which reports LHS
  // will print in the message output file. If LHSRPTS is omitted, the message
  // file will contain only the title, run specifications, and descriptions of
  // the distributions sampled. If LHSRPTS is included, it must be followed by
  // one or more of the following three additional keywords:
  //   > CORR: Print both the achieved raw correlation matrix and the achieved
  //           rank correlation matrix.
  //   > HIST: Print a text-based histogram for each random variable.
  //   > DATA: Print the complete set of all data samples and their ranks.
  // Pecos::LHSDriver::reportFlag is set from Dakota::Iterator::subIteratorFlag,
  // which accomplishes two things: (1) it reduces some output when generating
  // multiple sample sets (the report files get overwritten anyway), and (2) it
  // avoids numerical problems with generating input variable histogram plots
  // as trust regions become small in SBO (mainly an issue before conversion of
  // f90 LHS to double precision).
  String options_string = (reportFlag) ? "LHSRPTS CORR HIST DATA" : " ";
  options_string.resize(32, ' ');
  LHS_FILES2_FC(output_string.data(), message_string.data(),
                title_string.data(), options_string.data(), err_code);
  check_error(err_code, "lhs_files");

  //////////////////////////////////////////////////////////
  // Register RandomVariables with lhs_{dist,udist,const} //
  //////////////////////////////////////////////////////////
  RealArray dist_params;
  for (i=0; i<num_rv; ++i) {
    if (subset_rv && !active_vars[i]) continue; // skip this RV if not active

    const RandomVariable& rv_i = random_vars[i];
    switch (rv_i.type()) {
    case CONTINUOUS_RANGE: {
      Real l_bnd;  rv_i.pull_parameter(CR_LWR_BND, l_bnd);
      Real u_bnd;  rv_i.pull_parameter(CR_UPR_BND, u_bnd);
      check_finite(l_bnd, u_bnd);
      if (u_bnd > l_bnd) {
	dist_params.resize(2);
	dist_params[0] = l_bnd;  dist_params[1] = u_bnd;
	lhs_dist_register("ContRange", "uniform", i, dist_params);
      }
      else {
	check_range(l_bnd, u_bnd, true); // allow equal
	lhs_const_register("ContRange", i, l_bnd);
      }
      break;
    }
    case DISCRETE_RANGE: {
      int l_bnd;  rv_i.pull_parameter(DR_LWR_BND, l_bnd);
      int u_bnd;  rv_i.pull_parameter(DR_UPR_BND, u_bnd);
      check_finite(l_bnd, u_bnd);
      if (u_bnd > l_bnd) {
	RealArray x_val, y_val;
	int_range_to_udist_params(l_bnd, u_bnd, x_val, y_val);
	lhs_udist_register("DiscRange", "discrete histogram", i, x_val, y_val);
      }
      else {
	check_range(l_bnd, u_bnd, true); // allow equal
	lhs_const_register("DiscRange", i, l_bnd);
      }
      break;
    }
    case DISCRETE_SET_INT: {
      IntSet int_set;  rv_i.pull_parameter(DSI_VALUES, int_set);
      size_t set_size = int_set.size();
      if (set_size > 1) {
	RealArray x_val, y_val;
	set_to_udist_params(int_set, x_val, y_val); // set values
	lhs_udist_register("DiscSetI", "discrete histogram", i, x_val, y_val);
      }
      else if (set_size)
	lhs_const_register("DiscSetI", i, *int_set.begin());
      break;
    }
    case DISCRETE_SET_STRING: {
      StringSet str_set;  rv_i.pull_parameter(DSS_VALUES, str_set);
      int set_size = str_set.size();
      if (set_size > 1) {
	RealArray x_val, y_val;
	int_range_to_udist_params(0, set_size - 1, x_val, y_val); // indices
	lhs_udist_register("DiscSetS", "discrete histogram", i, x_val, y_val);
      }
      else if (set_size)
	lhs_const_register("DiscSetS", i, 0);
      break;
    }
    case DISCRETE_SET_REAL: {
      RealSet real_set;  rv_i.pull_parameter(DSR_VALUES, real_set);
      size_t set_size = real_set.size();
      if (set_size > 1) {
	RealArray x_val, y_val;
	set_to_udist_params(real_set, x_val, y_val); // set values
	lhs_udist_register("DiscSetR", "discrete histogram", i, x_val, y_val);
      }
      else if (set_size)
	lhs_const_register("DiscSetR", i, *real_set.begin());
      break;
    }
    case NORMAL: case STD_NORMAL:
      dist_params.resize(2);
      rv_i.pull_parameter(N_MEAN,    dist_params[0]);
      rv_i.pull_parameter(N_STD_DEV, dist_params[1]);
      lhs_dist_register("Normal", "normal", i, dist_params);
      break;
    case BOUNDED_NORMAL:
      dist_params.resize(4);
      rv_i.pull_parameter(N_MEAN,    dist_params[0]);
      rv_i.pull_parameter(N_STD_DEV, dist_params[1]);
      rv_i.pull_parameter(N_LWR_BND, dist_params[2]);
      rv_i.pull_parameter(N_UPR_BND, dist_params[3]);
      check_range(dist_params[2], dist_params[3], false);
      lhs_dist_register("Normal", "bounded normal", i, dist_params);
      break;
    case LOGNORMAL:     // LognormalRandomVariable standardizes on Lambda/Zeta
      dist_params.resize(2);
      rv_i.pull_parameter(LN_LAMBDA,  dist_params[0]);
      rv_i.pull_parameter(LN_ZETA,    dist_params[1]);
      lhs_dist_register("Lognormal", "lognormal-n", i, dist_params);
      break;
    case BOUNDED_LOGNORMAL: // BoundedLognormalRandomVariable uses Lambda/Zeta
      dist_params.resize(4);
      rv_i.pull_parameter(LN_LAMBDA,  dist_params[0]);
      rv_i.pull_parameter(LN_ZETA,    dist_params[1]);
      rv_i.pull_parameter(LN_LWR_BND, dist_params[2]);
      rv_i.pull_parameter(LN_UPR_BND, dist_params[3]);
      check_range(dist_params[2], dist_params[3], false);
      lhs_dist_register("Lognormal", "bounded lognormal-n", i, dist_params);
      break;
    case UNIFORM: case STD_UNIFORM:
      dist_params.resize(2);
      rv_i.pull_parameter(U_LWR_BND, dist_params[0]);
      rv_i.pull_parameter(U_UPR_BND, dist_params[1]);
      check_range(dist_params[0],  dist_params[1], false);
      check_finite(dist_params[0], dist_params[1]);
      lhs_dist_register("Uniform", "uniform", i, dist_params);
      break;
    case LOGUNIFORM:
      dist_params.resize(2);
      rv_i.pull_parameter(LU_LWR_BND, dist_params[0]);
      rv_i.pull_parameter(LU_UPR_BND, dist_params[1]);
      check_range(dist_params[0],  dist_params[1], false);
      check_finite(dist_params[0], dist_params[1]);
      lhs_dist_register("Loguniform", "loguniform", i, dist_params);
      break;
    case TRIANGULAR:
      dist_params.resize(2);
      rv_i.pull_parameter(T_LWR_BND, dist_params[0]);
      rv_i.pull_parameter(T_MODE,    dist_params[1]);
      rv_i.pull_parameter(T_UPR_BND, dist_params[2]);
      check_range(dist_params[0],  dist_params[2], false);
      check_finite(dist_params[0], dist_params[2]);
      lhs_dist_register("Triangular", "triangular", i, dist_params);
      break;
    case EXPONENTIAL: case STD_EXPONENTIAL: {
      dist_params.resize(1);
      Real e_beta; rv_i.pull_parameter(E_BETA, e_beta);
      dist_params[0] = 1./e_beta; // convert to LHS convention
      lhs_dist_register("Exponential", "exponential", i, dist_params);
      break;
    }
    case BETA: case STD_BETA:
      dist_params.resize(4);
      rv_i.pull_parameter(BE_LWR_BND, dist_params[0]);
      rv_i.pull_parameter(BE_UPR_BND, dist_params[1]);
      rv_i.pull_parameter(BE_ALPHA,   dist_params[2]);
      rv_i.pull_parameter(BE_BETA,    dist_params[3]);
      check_range(dist_params[0],  dist_params[1], false);
      check_finite(dist_params[0], dist_params[1]);
      lhs_dist_register("Beta", "beta", i, dist_params);
      break;
    case GAMMA: case STD_GAMMA: {
      dist_params.resize(2);
      rv_i.pull_parameter(GA_ALPHA, dist_params[0]);
      Real ga_beta; rv_i.pull_parameter(GA_BETA, ga_beta);
      dist_params[1] = 1./ga_beta; // convert to LHS convention
      lhs_dist_register("Gamma", "gamma", i, dist_params);
      break;
    }
    case GUMBEL:
      dist_params.resize(2);
      rv_i.pull_parameter(GU_ALPHA, dist_params[0]);
      rv_i.pull_parameter(GU_BETA,  dist_params[1]);
      lhs_dist_register("Gumbel", "gumbel", i, dist_params);
      break;
    case FRECHET:
      dist_params.resize(2);
      rv_i.pull_parameter(F_ALPHA, dist_params[0]);
      rv_i.pull_parameter(F_BETA,  dist_params[1]);
      lhs_dist_register("Frechet", "frechet", i, dist_params);
      break;
    case WEIBULL:
      dist_params.resize(2);
      rv_i.pull_parameter(W_ALPHA, dist_params[0]);
      rv_i.pull_parameter(W_BETA,  dist_params[1]);
      lhs_dist_register("Weibull", "weibull", i, dist_params);
      break;
    case HISTOGRAM_BIN: {
      RealRealMap h_bin_pairs;  rv_i.pull_parameter(H_BIN_PAIRS, h_bin_pairs);
      bins_to_udist_params(h_bin_pairs, x_val, y_val);
      lhs_udist_register("HistBin", "continuous linear", i, x_val, y_val);
      break;
    }
    case POISSON:
      dist_params.resize(1);
      rv_i.pull_parameter(P_LAMBDA, dist_params[0]);
      lhs_dist_register("Poisson", "poisson", i, dist_params);
      break;
    case BINOMIAL: {
      dist_params.resize(2);  int num_tr;
      rv_i.pull_parameter(BI_P_PER_TRIAL, dist_params[0]);
      rv_i.pull_parameter(BI_TRIALS, num_tr); dist_params[1] = (Real)num_tr;
      lhs_dist_register("Binomial", "binomial", i, dist_params);
      break;
    }
    case NEGATIVE_BINOMIAL: {
      dist_params.resize(2);  int num_tr;
      rv_i.pull_parameter(NBI_P_PER_TRIAL, dist_params[0]);
      rv_i.pull_parameter(NBI_TRIALS, num_tr); dist_params[1] = (Real)num_tr;
      lhs_dist_register("NegBinomial", "negative binomial", i, dist_params);
      break;
    }
    case GEOMETRIC:
      dist_params.resize(1);
      rv_i.pull_parameter(GE_P_PER_TRIAL, dist_params[0]);
      lhs_dist_register("Geometric", "geometric", i, dist_params);
      break;
    case HYPERGEOMETRIC: {
      dist_params.resize(3);  int tot_pop, num_drw, sel_pop;
      rv_i.pull_parameter(HGE_TOT_POP, tot_pop); dist_params[0] = (Real)tot_pop;
      rv_i.pull_parameter(HGE_DRAWN,   num_drw); dist_params[1] = (Real)num_drw;
      rv_i.pull_parameter(HGE_SEL_POP, sel_pop); dist_params[2] = (Real)sel_pop;
      lhs_dist_register("Hypergeom", "hypergeometric", i, dist_params);
      break;
    }
    case HISTOGRAM_PT_INT: {
      IntRealMap ir_map;  rv_i.pull_parameter(H_PT_INT_PAIRS, ir_map);
      size_t map_size = ir_map.size();
      if (map_size > 1) {
	RealArray x_val, y_val;
	map_to_udist_params(ir_map, x_val, y_val);
	lhs_udist_register("HistPtInt", "discrete histogram", i, x_val, y_val);
      }
      else if (map_size)
	lhs_const_register("HistPtInt", i, ir_map.begin()->first);
      break;
    }
    case HISTOGRAM_PT_STRING: {
      StringRealMap sr_map;  rv_i.pull_parameter(H_PT_STR_PAIRS, sr_map);
      int map_size = sr_map.size();
      if (map_size > 1) {
	RealArray x_val, y_val;
	map_indices_to_udist_params(sr_map, x_val, y_val);
	lhs_udist_register("HistPtString","discrete histogram",i, x_val, y_val);
      }
      else if (map_size)
	lhs_const_register("HistPtString", i, 0);
      break;
    }
    case HISTOGRAM_PT_REAL: {
      RealRealMap rr_map;  rv_i.pull_parameter(H_PT_REAL_PAIRS, rr_map);
      size_t map_size = rr_map.size();
      if (map_size > 1) {
	RealArray x_val, y_val;
	map_to_udist_params(rr_map, x_val, y_val);
	lhs_udist_register("HistPtReal", "discrete histogram", i, x_val, y_val);
      }
      else if (map_size)
	lhs_const_register("HistPtReal", i, rr_map.begin()->first);
      break;
    }
    case CONTINUOUS_INTERVAL_UNCERTAIN: {
      RealRealPairRealMap ci_bpa;  RealArray x_val, y_val;
      intervals_to_udist_params(ci_bpa, x_val, y_val);
      lhs_udist_register("ContInterval", "continuous linear", i, x_val, y_val);
      break;
    }
    case DISCRETE_INTERVAL_UNCERTAIN: {
      IntIntPairRealMap di_bpa;  RealArray x_val, y_val;
      intervals_to_udist_params(di_bpa, x_val, y_val);
      lhs_udist_register("DiscInterval", "discrete histogram", i, x_val, y_val);
      break;
    }
    case DISCRETE_UNCERTAIN_SET_INT: {
      IntRealMap ir_map;  rv_i.pull_parameter(DUSI_VALUES_PROBS, ir_map);
      size_t map_size = ir_map.size();
      if (map_size > 1) {
	RealArray x_val, y_val;
	map_to_udist_params(ir_map, x_val, y_val);
	lhs_udist_register("DiscUncSetI","discrete histogram", i, x_val, y_val);
      }
      else if (map_size)
	lhs_const_register("DiscUncSetI", i, ir_map.begin()->first);
      break;
    }
    case DISCRETE_UNCERTAIN_SET_STRING: {
      StringRealMap sr_map;  rv_i.pull_parameter(DUSS_VALUES_PROBS, sr_map);
      int map_size = sr_map.size();
      if (map_size > 1) {
	RealArray x_val, y_val;
	map_indices_to_udist_params(sr_map, x_val, y_val);
	lhs_udist_register("DiscUncSetS","discrete histogram", i, x_val, y_val);
      }
      else if (map_size)
	lhs_const_register("DiscUncSetS", i, 0);
      break;
    }
    case DISCRETE_UNCERTAIN_SET_REAL: {
      RealRealMap rr_map;  rv_i.pull_parameter(DUSR_VALUES_PROBS, rr_map);
      size_t map_size = rr_map.size();
      if (map_size > 1) {
	RealArray x_val, y_val;
	map_to_udist_params(rr_map, x_val, y_val);
	lhs_udist_register("DiscUncSetR","discrete histogram", i, x_val, y_val);
      }
      else if (map_size)
	lhs_const_register("DiscUncSetR", i, rr_map.begin()->first);
      break;
    }
    }
  }

  /////////////////////////////////////////
  // Register correlations with lhs_corr //
  /////////////////////////////////////////
  // specify the rank correlations among the RV.  Only non-zero values in the
  // lower triangular portion of the rank correlation matrix are specified.
  if (correlation_flag) {
    // Spec order: {cdv, ddv}, {cauv, dauv, corr}, {ceuv, deuv}, {csv, dsv}
    // > pass in a bit array for RV's to sample + another for keying corr's
    // > Default empty arrays --> all RVs active; corr matrix applies to all RVs
    size_t j, cntr_i, cntr_j;
    for (i=1, cntr_i=1; i<num_rv; ++i) {
      if (!subset_corr || active_corr[i]) {
	for (j=0, cntr_j=0; j<i; ++j) {
	  if (!subset_corr || active_corr[j]) {
	    if (!subset_rv || (active_vars[i] && active_vars[j])) {
	      Real corr_val = correlations(cntr_i,cntr_j);
	      if (std::abs(corr_val) > SMALL_NUMBER) {
		LHS_CORR2_FC(const_cast<char*>(lhsNames[i].data()),
			     const_cast<char*>(lhsNames[j].data()),
			     corr_val, err_code);
		check_error(err_code, "lhs_corr");
	      }
	    }
	    ++cntr_j;
	  }
	}
	++cntr_i;
      }
    }
  }

  /////////////////////
  // RUN THE SAMPLER //
  /////////////////////
  // perform internal checks on input to LHS
  int num_nam = num_rv, num_var = num_rv;
  LHS_PREP_FC(err_code, num_nam, num_var);
  check_error(err_code, "lhs_prep");

  // allocate the memory to hold samples, pt values, variable names, etc.
  int*      index_list = new  int    [num_nam];  // output
  Real*     ptval_list = new Real    [num_nam];  // output
  char* dist_name_list = new char [16*num_nam];  // output
  // dist_name_list is a bit tricky since the f90 array is declared as
  // CHARACTER(LEN=16) :: LSTNAME(MAXNAM), which would seem to be a
  // noncontiguous memory model.  However, a char** does not map correctly to
  // f90.  Rather, f90 takes the contiguous memory block from the C++ char*
  // allocation and indexes into it as if it were a vector of 16 char arrays
  // arranged head to tail.

  // The matrix of parameter samples from Fortran 90 is arranged in column-major
  // order with all variables for sample 1, followed by all variables for
  // sample 2, etc.  Teuchos::SerialDenseMatrix using column-major memory layout
  // as well, so use samples(var#,sample#) or samples[sample#][var#] for access.
  if (samples.numRows() != num_var || samples.numCols() != num_samples)
    samples.shapeUninitialized(num_var, num_samples);
  if (sampleRanksMode && sample_ranks.empty()) {
    if (sampleRanksMode == SET_RANKS || sampleRanksMode == SET_GET_RANKS) {
      PCerr << "Error: empty sample ranks array cannot be set in Pecos::"
	    << "LHSDriver::get_parameter_sets()" << std::endl;
      abort_handler(-1);
    }
    else if (sampleRanksMode == GET_RANKS)
      sample_ranks.shapeUninitialized(num_var, num_samples);
  }

  // generate the samples
  int rflag = sampleRanksMode; // short -> int
  LHS_RUN_FC(max_var, num_samples, num_nam, err_code, dist_name_list,
	     index_list, ptval_list, num_nam, samples.values(), num_var,
	     sample_ranks.values(), rflag);
  check_error(err_code, "lhs_run");

  // deallocate LHS memory
  LHS_CLOSE_FC(err_code);
  check_error(err_code, "lhs_close");

  // clean up memory
  delete [] index_list;
  delete [] ptval_list;
  delete [] dist_name_list;
#endif // HAVE_LHS
}


void LHSDriver::
generate_unique_samples(const std::vector<RandomVariables>& random_vars,
			const RealSymMatrix& correlations, int num_samples,
			RealMatrix& samples, RealMatrix& sample_ranks )
{
  // NonDSampling ordering of variables
  // Design    - continuous, discrete range, discrete set integer, 
  //             discrete set string, discrete set real
  // Aleatory  - continuous, discrete range, discrete set integer, 
  //             discrete set string, discrete set real
  // Epistemic - continuous, discrete range, discrete set integer, 
  //             discrete set string, discrete set real
  // State     - continuous, discrete range, discrete set integer, 
  //             discrete set string, discrete set real
 
  // LHSDriver ordering of variables returned in allSamples matrix:
  // Continuous - design, aleatory, epistemic, state
  // Discrete range - design,    Discrete set int - design
  // Discrete range - aleatory,  Discrete set int - aleatory
  // Discrete range - epistemic, Discrete set int - epistemic
  // Discrete range - state,     Discrete set int - state
  // Discrete set string - design, aleatory, epistemic, state
  // Discrete set real   - design, aleatory, epistemic, state

  // LHSDriver ordering of discrete uncertain variables returned in 
  // allSamples matrix
  // poisson uncertain, binomial uncertain, negative binomial uncertain, 
  // geometric uncertain, hypergeometric uncertain, histogram point int, 
  // discrete interval uncertain, discrete uncertain set integer

  // allSamples is a num_dims x num_samples RealMatrix
  // the values are stored in the following manner:
  // continuous: value
  // discrete set string: integer (mapped to Real) 
  //   I think an alphabetical ordering is applied internally by dakota so
  //   string values input with order 'b' 'a' are assigned int values of 1 0
  // discrete set int: value (mapped to Real)
  // discrete integer range: value (mapped to Real)
  // discrete set real: value

  // get total number of variables and number of discrete variables
  size_t num_cd_vars = cd_l_bnds.length(), num_cs_vars = cs_l_bnds.length(),
    num_cau_vars = adp.cauv(), num_ceu_vars = edp.ceuv(), 
    num_dau_vars = adp.dauv(), num_deu_vars = edp.deuv(),
    //num_daui_vars=adp.dauiv(),num_dausv_vars=adp.dausv(),num_daur_vars=daurv(),
    num_ddri_vars = ddri_l_bnds.length(), num_ddsi_vars = ddsi_values.size(),
    num_ddss_vars = ddss_values.size(), num_ddsr_vars = ddsr_values.size(),
    num_dsri_vars = dsri_l_bnds.length(), num_dssi_vars = dssi_values.size(),
    num_dsss_vars = dsss_values.size(), num_dssr_vars = dssr_values.size();
  size_t num_continuous_vars = num_cd_vars + num_cs_vars + num_cau_vars + 
    num_ceu_vars;
  size_t num_discrete_vars = num_dau_vars + num_deu_vars + num_ddri_vars + 
    num_ddsi_vars + num_ddss_vars + num_ddsr_vars + num_dsri_vars + 
    num_dssi_vars + num_dsss_vars + num_dssr_vars;
  size_t num_vars = num_continuous_vars + num_discrete_vars;

  // compute the number of total possible combinations of discrete variables
  // TODO Must look over all data structures passed into function
  // for range variables use ub-lb+1
  // for set variables use sum_i (ddsi_values[i].size())
  // if num_values**d < num_samples then call generate_samples
  int k=0;
  IntVector num_discrete_strata_1d( num_discrete_vars, false );
  // Discrete design variables
  for (int i=0; i<num_ddri_vars; i++)
    {  num_discrete_strata_1d[k] = ddri_u_bnds[i]-ddri_l_bnds[i]+1; k++; }
  for (int i=0; i<num_ddsi_vars; i++)
    {  num_discrete_strata_1d[k] = ddsi_values[i].size(); k++; }
  for (int i=0; i<num_ddss_vars; i++)
    {  num_discrete_strata_1d[k] = ddss_values[i].size(); k++; }
  for (int i=0; i<num_ddsr_vars; i++)
    {  num_discrete_strata_1d[k] = ddsr_values[i].size(); k++; }
  for (int i=0; i<num_dsri_vars; i++)

  // Discrete state variables
    {  num_discrete_strata_1d[k] = dsri_u_bnds[i]-dsri_l_bnds[i]+1; k++; }
  for (int i=0; i<num_dssi_vars; i++)
    {  num_discrete_strata_1d[k] = dssi_values.size(); k++; }
  for (int i=0; i<num_dsss_vars; i++)
    {  num_discrete_strata_1d[k] = dsss_values.size(); k++; }
  for (int i=0; i<num_dssr_vars; i++)
    {  num_discrete_strata_1d[k] = dssr_values.size(); k++; }

  // Epistemic discrete variables
  IntIntPairRealMapArray di_bpa = edp.discrete_interval_basic_probabilities();
  int num_diuv = di_bpa.size();
  for (int i=0; i<num_diuv; ++i){
    const IntIntPairRealMap& di_bpa_i = di_bpa[i];
    IIPRMCIter cit;
    // x_sort_unique contains ALL of the unique integer values for this
    // discrete interval variable in increasing order.  For example, if
    // there are 3 intervals for a variable and the bounds are (1,4),
    // (3,6), and [9,10], x_sorted will be (1,2,3,4,5,6,9,10).
    IntSet x_sort_unique;
    for (cit=di_bpa_i.begin(); cit!=di_bpa_i.end(); ++cit) {
      const RealRealPair& bounds = cit->first;
      int val, u_bnd = bounds.second;
      for (val=bounds.first; val<=u_bnd; ++val)
	x_sort_unique.insert(val);
    }
    int num_vals = x_sort_unique.size();
    num_discrete_strata_1d[k] = num_vals;
    k++;
  }
  IntRealMapArray dusi_vals_probs = edp.discrete_set_int_values_probabilities();
  int num_dusiv = dusi_vals_probs.size();
  for (int i=0; i<num_dusiv; ++i){
    const IntRealMap& dusi_v_p_i = dusi_vals_probs[i];
    int num_vals = dusi_v_p_i.size();
    num_discrete_strata_1d[k] = num_vals;
    k++;
  }
  StringRealMapArray duss_vals_probs = 
    edp.discrete_set_string_values_probabilities();
  int num_dussv = duss_vals_probs.size();
  for (int i=0; i<num_dussv; ++i){
    const StringRealMap& duss_v_p_i = duss_vals_probs[i];
    int num_vals = duss_v_p_i.size();
    num_discrete_strata_1d[k] = num_vals;
    k++;
  }
  RealRealMapArray dusr_vals_probs = 
    edp.discrete_set_real_values_probabilities();
  int num_dusrv = dusr_vals_probs.size();
  for (int i=0; i<num_dusrv; ++i){
    const RealRealMap& dusr_v_p_i = dusr_vals_probs[i];
    int num_vals = dusr_v_p_i.size();
    num_discrete_strata_1d[k] = num_vals;
    k++;
  }

  // Aleatory discrete variables with finite support
  IntVector bnt = adp.binomial_num_trials();
  for ( int i=0; i < bnt.length(); i++ )
    { num_discrete_strata_1d[k] = bnt[i]+1; k++;}
  IntVector htp = adp.hypergeometric_total_population();
  IntVector hsp = adp.hypergeometric_selected_population();
  IntVector hnd = adp.hypergeometric_num_drawn();
  for ( int i=0; i < htp.length(); i++ )
    { int num_fail=hnd[i], num_total_pop=htp[i], num_sel_pop=hsp[i]; 
      // Todo confirm this
      num_discrete_strata_1d[k] = std::min(num_fail,num_sel_pop)-
	std::max(0,num_sel_pop+num_fail-num_total_pop)+1;
      k++;
    }

  // Aleatory discrete variables with infinite support.  If any of these
  // variables are present then backfill can always be used.
  int max_num_unique_discrete_samples = 1;
  RealVector pl = adp.poisson_lambdas();
  RealVector nbppt = adp.negative_binomial_probability_per_trial();
  RealVector gppt = adp.geometric_probability_per_trial();
  if ( pl.length() > 0 ) max_num_unique_discrete_samples = INT_MAX;
  else if ( nbppt.length() > 0 ) max_num_unique_discrete_samples = INT_MAX;
  else if ( gppt.length() > 0 ) max_num_unique_discrete_samples = INT_MAX;
  else
    {
      num_discrete_strata_1d.resize( num_discrete_vars-pl.length()-
				     nbppt.length()- gppt.length());
      for (int k=0; k < num_discrete_strata_1d.length(); k++)
	  max_num_unique_discrete_samples *= num_discrete_strata_1d[k];
    }
  
  if ( max_num_unique_discrete_samples >= num_samples )
    // If the number of samples requested is greater than the maximum possible
    // number of discrete samples then we must allow replicates of the 
    // disscrete variables to obtain the desired number of variables. 
    // If not then we can proceed with generating a unique set of discrete 
    // samples.
    {
      if (samples.numRows() != num_vars || samples.numCols() != num_samples)
	samples.shapeUninitialized(num_vars, num_samples);
      // Currently sample_ranks will always be returned empty. It should only be
      // filled when NonDSampling.sampleRanksMode>0. But I cannot see anywhere
      // in the code where this is true.
      //sample_ranks.shapeUninitialized( num_vars, num_samples );

      RealMatrix sample_ranks_rm, samples_rm;

      // unique index of all discrete variables if any
      std::set<RealArray>::iterator it;
      std::set<RealArray> sorted_discrete_samples; 
      RealArray discrete_sample( num_discrete_vars );

      // Determine the columns in samples_rm that contain discrete variables
      IntVector discrete_samples_map( num_discrete_vars, false );
      for (int i=0; i<num_discrete_vars; i++)
	discrete_samples_map[i]=num_continuous_vars+i;
  
      // Eliminate redundant samples by resampling if necessary.  Could pad
      // num_samples in anticipation of duplicates, but this would alter LHS
      // stratification that could be intended, so use num_samples for now.
      bool complete = false, initial = true;
      int num_unique_samples = 0;
      while (!complete) {
	generate_samples(random_vars, correlations, num_samples,
			 samples_rm, sample_ranks_rm);

	if (initial) { // pack initial sample set
	  for (int i=0; i<num_samples; ++i) { // or matrix->set<vector> ?
	    //PCout << "[";
	    for (int j=0; j<num_discrete_vars; j++) {
	      int index = discrete_samples_map[j];
	      discrete_sample[j] = samples_rm(index,i);
	      //PCout << discrete_sample[j] << ",";
	    }
	    //PCout << "]\n";
	    sorted_discrete_samples.insert(discrete_sample);
	    if ( sorted_discrete_samples.size() > num_unique_samples ){
	      // copy sample into samples matrix
	      for (int j=0; j<num_vars; j++)
		samples(j,num_unique_samples) = samples_rm(j,i);
	      num_unique_samples++;
	    }
	  }
	  if (num_unique_samples == num_samples) complete = true;
	  else initial = false;
	}
	else { // backfill duplicates with new samples
	  //PCout << num_unique_samples << "," << sorted_discrete_samples.size()
	  //      << "," << num_discrete_vars << "," << num_vars << ","
	  //      << num_continuous_vars << std::endl;

	  for (int i=0; i<num_samples; ++i) {
	    if (num_unique_samples < num_samples) {
	      //PCout << "[";
	      for (int j=0; j<num_discrete_vars; j++) {
		int index = discrete_samples_map[j];
		discrete_sample[j] = samples_rm(index,i);
		//PCout << discrete_sample[j] << ",";
	      }
	      //PCout << "]\n";
	      sorted_discrete_samples.insert(discrete_sample);
	      if ( sorted_discrete_samples.size() > num_unique_samples ){
		// copy sample into samples matrix
		for (int j=0; j<num_vars; j++)
		  samples(j,num_unique_samples) = samples_rm(j,i);
		num_unique_samples++;
	      }
	    }
	    else
	      { complete = true; break; }
	  }
	}
      }
    }
  else {
    PCout << "LHS backfill was requested, but the discrete variables "
	  << "specified do not have enough unique values ("
	  << max_num_unique_discrete_samples
	  << ") to obtain the number of samples requested so replicated "
	  << "discrete samples have been allowed.\n";
    generate_samples(random_vars, correlations, num_samples, samples,
		     sample_ranks);
  }
}


void LHSDriver::
generate_unique_index_samples(const IntVector& index_l_bnds,
			      const IntVector& index_u_bnds, int num_samples,
			      IntMatrix& index_samples )
{
  // For    uniform probability, model as discrete design range (this fn).
  // For nonuniform probability, model as discrete uncertain set integer.

  RealVector  empty_rv;  RealMatrix empty_rm, samples_rm;
  IntVector   empty_iv;  IntArray sample;
  IntSetArray empty_isa; StringSetArray empty_ssa; RealSetArray empty_rsa;
  AleatoryDistParams adp; EpistemicDistParams edp;
  generate_unique_samples( empty_rv, empty_rv, index_l_bnds, index_u_bnds,
			   empty_isa, empty_ssa, empty_rsa, empty_rv, empty_rv,
			   empty_iv, empty_iv, empty_isa, empty_ssa, empty_rsa,
			   adp, edp, num_samples, samples_rm, empty_rm);
  
  index_samples.shapeUninitialized( samples_rm.numRows(), samples_rm.numCols() );
  for (int j=0; j<samples_rm.numCols(); j++){
    for (int i=0; i<samples_rm.numRows(); i++)
      index_samples(i,j) = (int)samples_rm(i,j);
  }
}

} // namespace Pecos

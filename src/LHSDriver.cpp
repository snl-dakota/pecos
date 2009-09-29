/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "LHSDriver.hpp"
#include "pecos_stat_util.hpp"
#include <algorithm>

static const char rcsId[]="@(#) $Id: LHSDriver.cpp 5248 2008-09-05 18:51:52Z wjbohnh $";


#ifdef HAVE_LHS
#define LHS_INIT_MEM_FC FC_FUNC_(lhs_init_mem,LHS_INIT_MEM)
#define LHS_PREP_FC     FC_FUNC_(lhs_prep,LHS_PREP)
#define LHS_RUN_FC      FC_FUNC_(lhs_run,LHS_RUN)
#define LHS_CLOSE_FC    FC_FUNC_(lhs_close,LHS_CLOSE)
#define LHS_OPTIONS2_FC FC_FUNC_(lhs_options2,LHS_OPTIONS2)
#define LHS_DIST2_FC    FC_FUNC_(lhs_dist2,LHS_DIST2)
#define LHS_UDIST2_FC   FC_FUNC_(lhs_udist2,LHS_UDIST2)
#define LHS_CORR2_FC    FC_FUNC_(lhs_corr2,LHS_CORR2)
#define LHS_FILES2_FC   FC_FUNC_(lhs_files2,LHS_FILES2)

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
		      char* sampling_options, int& ierror );

void LHS_DIST2_FC( char* label, int& ptval_flag, Pecos::Real& ptval,
		   char* dist_type, Pecos::Real* dist_params, int& num_params,
		   int& ierror, int& dist_id, int& ptval_id );

void LHS_UDIST2_FC( char* label, int& ptval_flag, Pecos::Real& ptval,
		    char* dist_type, int& num_pts, Pecos::Real* x,
		    Pecos::Real* y, int& ierror, int& dist_id, int& ptval_id );

void LHS_CORR2_FC( char* label1, char* label2, Pecos::Real& corr, int& ierror );

void LHS_FILES2_FC( char* lhsout, char* lhsmsg, char* lhstitl, char* lhsopts,
		    int& ierror );

//void lhs_run2( int* max_var, int* max_obs, int* max_names, int* ierror,
//               char* dist_names, int* name_order, double* ptvals,
//               int* num_names, double* sample_matrix, int* num_vars );

}
#endif // HAVE_LHS


namespace Pecos {


/** Helper function to copy names that Fortran will see as
    character*16 values, and are null-terminated for assignment
    to String values.  Fortran won't see the null at the end. */
 static void
f77name16(char buf[17], const char *name, StringArray &lhs_names, int cntr)
{
	char *b, *be;
	for(b = buf, be = buf + 16; b < be; ++b)
		if (!(*b = *name++)) {
			b += snprintf(b, be-b, "%d", cntr+1);
			while(b < be)
				*b++ = ' ';
			break;
			}
	*be = 0;
	lhs_names[cntr] = buf;
	}


 static void
f77dist32(char buf[32], const char *name) // blank-filled but not null-terminated
{
	char *b, *be;
	for(b = buf, be = buf + 32; b < be; ++b)
		if (!(*b = *name++)) {
			do *b++ = ' ';
			   while(b < be);
			break;
			}
	}


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
generate_samples(const RealVector& d_l_bnds,     const RealVector& d_u_bnds,
		 const RealVector& s_l_bnds,     const RealVector& s_u_bnds,
		 const RealVector& n_means,      const RealVector& n_std_devs,
		 const RealVector& n_l_bnds,     const RealVector& n_u_bnds,
		 const RealVector& ln_means,     const RealVector& ln_std_devs,
		 const RealVector& ln_lambdas,   const RealVector& ln_zetas,
		 const RealVector& ln_err_facts, const RealVector& ln_l_bnds,
		 const RealVector& ln_u_bnds,	 const RealVector& u_l_bnds,
		 const RealVector& u_u_bnds,	 const RealVector& lu_l_bnds,
		 const RealVector& lu_u_bnds,    const RealVector& t_modes,
		 const RealVector& t_l_bnds,     const RealVector& t_u_bnds,
		 const RealVector& e_betas,      const RealVector& b_alphas,
		 const RealVector& b_betas,      const RealVector& b_l_bnds,
		 const RealVector& b_u_bnds,     const RealVector& ga_alphas,
		 const RealVector& ga_betas,     const RealVector& w_alphas,
		 const RealVector& w_betas,      const RealVector& gu_alphas,
		 const RealVector& gu_betas,     const RealVector& f_alphas,
		 const RealVector& f_betas,
		 const RealVectorArray& h_bin_prs,
		 const RealVectorArray& h_pt_prs,
		 const RealVectorArray& i_probs, const RealVectorArray& i_bnds,
		 const RealSymMatrix& correlations, int num_samples,
		 RealMatrix& samples, RealMatrix& sample_ranks)
{
#ifdef HAVE_LHS
  // generate samples within user-specified parameter distributions

  // error check on program parameters
  if (!num_samples) {
    PCerr << "\nError: number of samples in LHSDriver::generate_samples() must "
	  << "be nonzero." << std::endl;
    abort_handler(-1);
  }

  bool correlation_flag = !correlations.empty();
  size_t i, j, num_dv = d_l_bnds.length(), num_sv = s_l_bnds.length(),
    num_nuv = n_means.length(),
    num_lnuv = std::max(ln_means.length(), ln_lambdas.length()),
    num_uuv = u_l_bnds.length(), num_luuv = lu_l_bnds.length(),
    num_tuv = t_modes.length(),  num_euv  = e_betas.length(),
    num_buv = b_alphas.length(), num_gauv = ga_alphas.length(),
    num_wuv = w_alphas.length(), num_guuv = gu_alphas.length(),
    num_fuv = f_alphas.length(), num_huv  = h_bin_prs.size() + h_pt_prs.size(),
    num_iuv = i_probs.size(),
    num_uv  = num_nuv + num_lnuv + num_uuv + num_luuv + num_tuv + num_euv
            + num_buv + num_gauv + num_wuv + num_guuv + num_fuv
            + num_huv + num_iuv,
    num_av  = num_dv + num_uv + num_sv;

  int err_code = 0, max_var = num_av, max_obs = num_samples,
      max_samp_size = num_av*num_samples, max_interval = -1,
      max_unc_corr = (num_uv*num_uv - num_uv)/2,
      max_table = -1, print_level = 0, output_width = 1;
  int max_corr = (num_uv > 1) ? max_unc_corr : -1;
  //lhs_init(num_samples, randomSeed, err_code);
  LHSDriver::seed(randomSeed);
  LHS_INIT_MEM_FC(num_samples, randomSeed, max_obs, max_samp_size, max_var,
		  max_interval, max_corr, max_table, print_level, output_width,
		  err_code);
  check_error(err_code, "lhs_init_mem");

  // set sample type to either LHS (default) or random Monte Carlo (optional)
  bool call_lhs_option = false;
  char option_string[32];
  std::ostringstream option_stream;
  if (sampleType == "random" || sampleType == "incremental_random") {
    option_stream << "RANDOM SAMPLE ";
    call_lhs_option = true;
  }
  else
    option_stream << "              ";
  // set mixing option to either restricted pairing (default) or random pairing
  // (optional).  For enforcing user-specified correlation, restricted pairing
  // is required.  And for uncorrelated variables, restricted pairing results
  // in lower correlation values than random pairing.  For these reasons, the
  // random pairing option is not currently active, although a specification
  // option for it could be added in the future if a use arises.
  bool random_pairing_flag = false; // this option hard-wired off for now
  if (!correlation_flag && random_pairing_flag) {
    option_stream << "RANDOM PAIRING    ";
    call_lhs_option = true;
  }
  else // use default of restricted pairing
    option_stream << "                  ";
  if (call_lhs_option) {
    // Don't null-terminate the string since the '\0' is not used in Fortran
    //option_stream << ends;
    int num_replic = 1, ptval_option = 1;
    String option_s(option_stream.str());
    std::copy(option_s.begin(), option_s.end(), option_string);
    LHS_OPTIONS2_FC(num_replic, ptval_option, option_string, err_code);
    check_error(err_code, "lhs_options");
  }

  int num_params, cntr = 0, ptval_flag = 0;
  int dist_num, pv_num; // outputs (ignored)
  Real ptval = 0., dist_params[4];
  StringArray lhs_names(num_av);
  char dist_string[32], name_string[17];
  const char *distname;

  // design (treated as uniform)
  for (i=0; i<num_dv; i++, cntr++) {
    f77name16(name_string, "Design", lhs_names, cntr);
    f77dist32(dist_string, "uniform");
    num_params = 2;
    if (d_l_bnds[i] > -DBL_MAX && d_u_bnds[i] < DBL_MAX) {
      if (d_l_bnds[i] >= d_u_bnds[i]) {
	PCerr << "\nError: Pecos::LHSDriver requires lower bounds strictly "
	      << "less than upper bounds to\n       sample design variables "
	      << "using uniform distributions." << std::endl;
	abort_handler(-1);
      }
      dist_params[0] = d_l_bnds[i];
      dist_params[1] = d_u_bnds[i];
    }
    else {
      PCerr << "\nError: Pecos::LHSDriver requires bounds to sample design "
	    << "variables using uniform\n       distributions." << std::endl;
      abort_handler(-1);
    }
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string, dist_params,
		 num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(design)");
  }

  // normal uncertain
  bool n_bnd_spec = (!n_l_bnds.empty() && !n_u_bnds.empty());
  for (i=0; i<num_nuv; i++, cntr++) {
    f77name16(name_string, "Normal", lhs_names, cntr);
    dist_params[0] = n_means[i];
    dist_params[1] = n_std_devs[i];
    // check for bounded normal
    if ( n_bnd_spec && (n_l_bnds[i] > -DBL_MAX || n_u_bnds[i] < DBL_MAX) ) {
      if (n_l_bnds[i] >= n_u_bnds[i]) {
	PCerr << "\nError: Pecos::LHSDriver requires lower bounds strictly "
	      << "less than upper bounds to\n       sample using bounded "
	      << "normal distributions." << std::endl;
	abort_handler(-1);
      }
      num_params = 4;
      dist_params[2] = n_l_bnds[i];
      dist_params[3] = n_u_bnds[i];
      distname = "bounded normal";
    }
    else { // normal
      num_params = 2;
      distname = "normal";
    }
    f77dist32(dist_string, distname);
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string, dist_params,
		 num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(normal)");
  }

  // lognormal uncertain
  bool ln_bnd_spec = (!ln_l_bnds.empty()  && !ln_u_bnds.empty());
  bool n_dist      = (!ln_lambdas.empty() || !ln_std_devs.empty());
  for (i=0; i<num_lnuv; i++, cntr++) {
    f77name16(name_string, "Lognormal", lhs_names, cntr);
    if (n_dist) {
      // In the mean/std dev specification case, LHS expects the mean/std dev
      // of the underlying normal distribution (LHS manual, SAND#98-0210, p.39).
      // Therefore, map from the DAKOTA spec (mean/std_dev or lambda/zeta) to
      // the LHS input requirements (lambda/zeta), if required.
      if (!ln_lambdas.empty()) {
	dist_params[0] = ln_lambdas[i];
	dist_params[1] = ln_zetas[i];
      }
      else
	lognormal_params_from_moments(ln_means[i], ln_std_devs[i],
				      dist_params[0], dist_params[1]);
    }
    else {
      // In the mean/err factor specification case, DAKOTA and LHS are
      // consistent (LHS manual, SAND#98-0210, p.39) and no mapping is required.
      dist_params[0] = ln_means[i];
      dist_params[1] = ln_err_facts[i];
    }
    // check for bounded lognormal
    if (ln_bnd_spec && (ln_l_bnds[i] > 0.0 || ln_u_bnds[i] < DBL_MAX) ) {
      if (ln_l_bnds[i] >= ln_u_bnds[i]) {
	PCerr << "\nError: Pecos::LHSDriver requires lower bounds strictly "
	      << "less than upper bounds to\n       sample using bounded "
	      << "lognormal distributions." << std::endl;
	abort_handler(-1);
      }
      num_params = 4;
      dist_params[2] = ln_l_bnds[i];
      dist_params[3] = ln_u_bnds[i];
      if (n_dist)
	distname = "bounded lognormal-n";
      else
	distname = "bounded lognormal";
    }
    else {
      num_params = 2;
      if (n_dist)
	distname = "lognormal-n";
      else
	distname = "lognormal";
    }
    f77dist32(dist_string, distname);
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string, dist_params,
		 num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(lognormal)");
  }

  // uniform uncertain
  for (i=0; i<num_uuv; i++, cntr++) {
    f77name16(name_string, "Uniform", lhs_names, cntr);
    f77dist32(dist_string, "uniform");
    num_params = 2;
    if (u_l_bnds[i] >= u_u_bnds[i]) {
	PCerr << "\nError: Pecos::LHSDriver requires lower bounds strictly "
	      << "less than upper bounds to\n       sample using uniform "
	      << "distributions." << std::endl;
	abort_handler(-1);
    }
    dist_params[0] = u_l_bnds[i];
    dist_params[1] = u_u_bnds[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string, dist_params,
		 num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(uniform)");
  }

  // loguniform uncertain
  for (i=0; i<num_luuv; i++, cntr++) {
    f77name16(name_string, "Loguniform", lhs_names, cntr);
    f77dist32(dist_string, "loguniform");
    num_params = 2;
    if (lu_l_bnds[i] >= lu_u_bnds[i]) {
	PCerr << "\nError: Pecos::LHSDriver requires lower bounds strictly "
	      << "less than upper bounds to\n       sample using loguniform "
	      << "distributions." << std::endl;
	abort_handler(-1);
    }
    dist_params[0] = lu_l_bnds[i];
    dist_params[1] = lu_u_bnds[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string, dist_params,
		 num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(loguniform)");
  }

  // triangular uncertain
  for (i=0; i<num_tuv; i++, cntr++) {
    f77name16(name_string, "Triangular", lhs_names, cntr);
    f77dist32(dist_string, "triangular");
    num_params = 3;
    if (t_l_bnds[i] >= t_u_bnds[i]) {
      PCerr << "\nError: Pecos::LHSDriver requires lower bounds strictly less "
	    << "than upper bounds to\n       sample using triangular "
	    << "distributions." << std::endl;
      abort_handler(-1);
    }
    dist_params[0] = t_l_bnds[i];
    dist_params[1] = t_modes[i];
    dist_params[2] = t_u_bnds[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string, dist_params,
		 num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(triangular)");
  }

  // exponential uncertain
  for (i=0; i<num_euv; i++, cntr++) {
    f77name16(name_string, "Exponential", lhs_names, cntr);
    f77dist32(dist_string, "exponential");
    num_params = 1;
    dist_params[0] = 1./e_betas[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string, dist_params,
		 num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(exponential)");
  }

  // beta uncertain
  for (i=0; i<num_buv; i++, cntr++) {
    f77name16(name_string, "Beta", lhs_names, cntr);
    f77dist32(dist_string, "beta");
    num_params = 4;
    if (b_l_bnds[i] >= b_u_bnds[i]) {
      PCerr << "\nError: Pecos::LHSDriver requires lower bounds strictly less "
	    << "than upper bounds to\n       sample using beta "
	    << "distributions." << std::endl;
      abort_handler(-1);
    }
    dist_params[0] = b_l_bnds[i];
    dist_params[1] = b_u_bnds[i];
    dist_params[2] = b_alphas[i];
    dist_params[3] = b_betas[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string, dist_params,
		 num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(beta)");
  }

  // gamma uncertain
  for (i=0; i<num_gauv; i++, cntr++) {
    f77name16(name_string, "Gamma", lhs_names, cntr);
    f77dist32(dist_string, "gamma");
    num_params = 2;
    dist_params[0] = ga_alphas[i];
    dist_params[1] = 1./ga_betas[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string, dist_params,
		 num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(gamma)");
  }

  // gumbel uncertain
  for (i=0; i<num_guuv; i++, cntr++) {
    f77name16(name_string, "Gumbel", lhs_names, cntr);
    f77dist32(dist_string, "gumbel");
    num_params = 2;
    dist_params[0] = gu_alphas[i];
    dist_params[1] = gu_betas[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string, dist_params,
		 num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(gumbel)");
  }

  // frechet uncertain
  for (i=0; i<num_fuv; i++, cntr++) {
    f77name16(name_string, "Frechet", lhs_names, cntr);
    f77dist32(dist_string, "frechet");
    num_params = 2;
    dist_params[0] = f_alphas[i];
    dist_params[1] = f_betas[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string, dist_params,
		 num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(frechet)");
  }

  // weibull uncertain
  for (i=0; i<num_wuv; i++, cntr++) {
    f77name16(name_string, "Weibull", lhs_names, cntr);
    f77dist32(dist_string, "weibull");
    num_params = 2;
    dist_params[0] = w_alphas[i];
    dist_params[1] = w_betas[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string, dist_params,
		 num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(weibull)");
  }

  // histogram uncertain: in both formats, pairs are defined from an abscissa
  // in the first field and a count (not a density) in the second field.  The
  // distinction in the second field is only important for unequal bin widths.
  size_t num_c_huv = h_bin_prs.size();
  for (i=0; i<num_huv; i++, cntr++) {
    f77name16(name_string, "Histogram", lhs_names, cntr);
     num_params = (i < num_c_huv) ? h_bin_prs[i].length()/2
                                 : h_pt_prs[i-num_c_huv].length()/2;
    Real* x_val = new Real [num_params];
    Real* y_val = new Real [num_params];
    if (i < num_c_huv) {
      // LHS requires accumulation of CDF with first y at 0 and last y at 1
      distname = "continuous linear";
      for (j=0; j<num_params; j++)
        x_val[j] = h_bin_prs[i][2*j];
      // Assume already normalized with sum = 1
      //Real sum = 0.;
      //for (j=1; j<num_params; j++)
      //  sum += h_bin_prs[i][2*j-1]; // last y from DAKOTA must be zero
      y_val[0] = 0.;
      for (j=1; j<num_params; j++)
        y_val[j] = y_val[j-1] + h_bin_prs[i][2*j-1];// / sum;
    }
    else {
      // LHS can use discrete frequency information directly
      distname = "discrete histogram";
      for (j=0; j<num_params; j++) {
        x_val[j] = h_pt_prs[i-num_c_huv][2*j];
        y_val[j] = h_pt_prs[i-num_c_huv][2*j+1];
      }
    }
#ifdef DEBUG
    for (j=0; j<num_params; j++)
      PCout << "Histogram " << i+1 << ": " << x_val[j] << ' ' << y_val[j]
	    << std::endl;
#endif //DEBUG
    f77dist32(dist_string, distname);
    LHS_UDIST2_FC(name_string, ptval_flag, ptval, dist_string, num_params,
		  x_val, y_val, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_udist(histogram)");
    delete [] x_val;
    delete [] y_val;
  }

  // interval uncertain: convert to histogram for sampling
  for (i=0; i<num_iuv; i++, cntr++) {
    f77name16(name_string, "Interval", lhs_names, cntr);
    f77dist32(dist_string, "continuous linear");

    // x_sort_unique is a set with ALL of the interval bounds for this variable
    // in increasing order and unique.  For example, if there are 2 intervals
    // for a variable, and the bounds are (1,4) and (3,6), x_sorted will be
    // (1, 3, 4, 6).  If the intervals are contiguous, e.g. one interval is
    // (1,3) and the next is (3,5), x_sort_unique is (1,3,5).
    const RealVector& interval_probs_i = i_probs[i];
    const RealVector& interval_bnds_i  = i_bnds[i];
    RealSet x_sort_unique;
    int num_intervals_i = interval_probs_i.length(),
        num_bounds_i    = 2*num_intervals_i;
    for (j=0; j<num_bounds_i; j++)
      x_sort_unique.insert(interval_bnds_i[j]);
    // convert RealSet to Real*
    num_params = x_sort_unique.size();
    Real* x_val = new Real [num_params];
    RealSet::iterator it = x_sort_unique.begin();
    for (j=0; j<num_params; j++, it++)
      x_val[j] = *it;

    // Calculate the probability densities, and account for the cases where
    // there are intervals that are overlapping.  This section of code goes
    // through the original intervals and see where they fall relative to the
    // new, sorted intervals for the density calculation.
    RealVector prob_dens(num_params); // initialize to 0.
    for (j=0; j<num_intervals_i; j++) {
      const Real& lower_value = interval_bnds_i[2*j];
      const Real& upper_value = interval_bnds_i[2*j+1];
      Real interval_density = interval_probs_i[j] / (upper_value - lower_value);
      int cum_int_index = 0;
      while (lower_value > x_val[cum_int_index])
	cum_int_index++;
      cum_int_index++;
      while (cum_int_index < num_params && x_val[cum_int_index] <= upper_value){
	prob_dens[cum_int_index] += interval_density;
	cum_int_index++;
      }
    }

    // put the densities in a cumulative format necessary for
    // the LHS histogram variable.  Note that x_val and y_val
    // are defined as Real* for input to the Fortran call.
    Real* y_val = new Real [num_params];
    y_val[0] = 0.;
    for (j=1; j<num_params; j++) {
      if (prob_dens[j] > 0.0)
	y_val[j] = y_val[j-1] + prob_dens[j] * (x_val[j] - x_val[j-1]);
      else // handle case where there is a gap
	y_val[j] = y_val[j-1] + 0.0001;
    }
    // normalize if necessary
    if (y_val[num_params-1] != 1.0) {
      Real y_total = y_val[num_params-1];
      for (j=1; j<num_params; j++)
	y_val[j] /= y_total;
    }

#ifdef DEBUG
    for (j=0;j<num_params;j++)
      PCout << "x_val " << j << "is " << x_val[j] << "\ny_val " << j << "is "
	    << y_val[j]<< '\n';
#endif //DEBUG
    LHS_UDIST2_FC(name_string, ptval_flag, ptval, dist_string, num_params,
		  x_val, y_val, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_udist(interval)");
    delete [] x_val;
    delete [] y_val;
  }

  // state (treated as uniform)
  for (i=0; i<num_sv; i++, cntr++) {
    f77name16(name_string, "State", lhs_names, cntr);
    f77dist32(dist_string, "uniform");
    num_params = 2;
    if (s_l_bnds[i] > -DBL_MAX && s_u_bnds[i] < DBL_MAX) {
      if (s_l_bnds[i] >= s_u_bnds[i]) {
	PCerr << "\nError: Pecos::LHSDriver requires lower bounds strictly "
	      << "less than upper bounds to\n       sample state variables "
	      << "using uniform distributions." << std::endl;
	abort_handler(-1);
      }
      dist_params[0] = s_l_bnds[i];
      dist_params[1] = s_u_bnds[i];
    }
    else {
      PCerr << "\nError: Pecos::LHSDriver requires bounds to sample state "
	    << "variables using uniform\n       distributions." << std::endl;
      abort_handler(-1);
    }
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string, dist_params,
		 num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(state)");
  }

  // specify the rank correlations among the uncertain vars (no correlation
  // currently supported for design and state vars in allVars mode).  Only
  // non-zero values in the lower triangular portion of the rank correlation
  // matrix are specified.
  if (correlation_flag) {
    for (i=1; i<num_uv; i++) {
      for (j=0; j<i; j++) {
	Real corr_val = correlations(i,j);
	if (fabs(corr_val) > 1.e-25) {
	  LHS_CORR2_FC(const_cast<char*>(lhs_names[i+num_dv].c_str()),
		       const_cast<char*>(lhs_names[j+num_dv].c_str()),
		       corr_val, err_code);
	  check_error(err_code, "lhs_corr");
	}
      }
    }
  }

  // Create files showing distributions and associated statistics.  Avoid
  // issues with null-terminated strings from C++ (which mess up the Fortran
  // output) by using ostrstreams without ends insertions.
  char output_string[32], message_string[32], title_string[32],
    options_string[32];
  std::ostringstream output_stream, message_stream, title_stream,
    options_stream;
  output_stream  << "LHS_samples.out                 ";
  message_stream << "LHS_distributions.out           ";
  title_stream   << "DAKOTA LHS                      ";
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
  if (reportFlag)
    options_stream << "LHSRPTS CORR HIST DATA          ";
  else
    options_stream << "                                ";
  String output_s(output_stream.str());
  std::copy(output_s.begin(), output_s.end(), output_string);
  String message_s(message_stream.str());
  std::copy(message_s.begin(), message_s.end(), message_string);
  String title_s(title_stream.str());
  std::copy(title_s.begin(), title_s.end(), title_string);
  String options_s(options_stream.str());
  std::copy(options_s.begin(), options_s.end(), options_string);
  LHS_FILES2_FC(output_string, message_string, title_string, options_string,
		err_code);
  check_error(err_code, "lhs_files");

  // perform internal checks on input to LHS
  int num_nam = num_av, num_var = num_av;
  LHS_PREP_FC(err_code, num_nam, num_var);
  check_error(err_code, "lhs_prep");

  // allocate the memory to hold samples, pt values, variable names, etc.
  int   max_nam        = num_av;
  int*  index_list     = new int    [max_nam];       // output
  Real* ptval_list     = new Real   [max_nam];       // output
  char* dist_name_list = new char   [16*max_nam];    // output
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
  if (samples.empty())
    samples.shapeUninitialized(num_av, num_samples);
  if (sampleRanksMode && sample_ranks.empty()) {
    if (sampleRanksMode == SET_RANKS || sampleRanksMode == SET_GET_RANKS) {
      PCerr << "Error: empty sample ranks array cannot be set in Pecos::"
	    << "LHSDriver::get_parameter_sets()" << std::endl;
      abort_handler(-1);
    }
    else if (sampleRanksMode == GET_RANKS)
      sample_ranks.shapeUninitialized(num_av, num_samples);
  }

  // generate the samples
  int rflag = sampleRanksMode; // short -> int
  LHS_RUN_FC(max_var, max_obs, max_nam, err_code, dist_name_list,
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

} // namespace Pecos

/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef LHS_DRIVER_H
#define LHS_DRIVER_H

#include "pecos_data_types.hpp"


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
	    bool reportFlag = true);  ///< constructor
  ~LHSDriver(); ///< destructor

  //
  //- Heading: Member functions
  //

  /// populate data when not passed through ctor
  void initialize(const String& sample_type,
		  short sample_ranks_mode = IGNORE_RANKS,
		  bool reportFlag = true);

  /// generates the desired set of parameter samples from within general
  /// user-specified probabilistic distributions
  void generate_samples(const RealArray& d_l_bnds,  const RealArray& d_u_bnds,
			const RealArray& s_l_bnds,  const RealArray& s_u_bnds,
			const RealArray& n_means,   const RealArray& n_std_devs,
			const RealArray& n_l_bnds,  const RealArray& n_u_bnds,
			const RealArray& ln_means,
			const RealArray& ln_std_devs,
			const RealArray& ln_err_facts,
			const RealArray& ln_l_bnds,
			const RealArray& ln_u_bnds, const RealArray& u_l_bnds,
			const RealArray& u_u_bnds,  const RealArray& lu_l_bnds,
			const RealArray& lu_u_bnds, const RealArray& t_modes,
			const RealArray& t_l_bnds,  const RealArray& t_u_bnds,
			const RealArray& e_betas,   const RealArray& b_alphas,
			const RealArray& b_betas,   const RealArray& b_l_bnds,
			const RealArray& b_u_bnds,  const RealArray& ga_alphas,
			const RealArray& ga_betas,  const RealArray& w_alphas,
			const RealArray& w_betas,   const RealArray& gu_alphas,
			const RealArray& gu_betas,  const RealArray& f_alphas,
			const RealArray& f_betas, const Real2DArray& h_bin_prs,
			const Real2DArray& h_pt_prs, const Real2DArray& i_probs,
			const Real2DArray& i_bnds,
			const Real2DArray& correlations,
			int num_samples,            int seed, 
			Real2DArray& samples_array, Real2DArray& rank_array);

  /// generates the desired set of parameter samples from within
  /// uncorrelated normal distributions
  void generate_normal_samples(const RealArray& n_means,
			       const RealArray& n_std_devs,
			       const RealArray& n_l_bnds,
			       const RealArray& n_u_bnds, int num_samples,
			       int seed, Real2DArray& samples_array);

  /// generates the desired set of parameter samples from within
  /// uncorrelated uniform distributions
  void generate_uniform_samples(const RealArray& u_l_bnds,
				const RealArray& u_u_bnds, int num_samples,
				int seed, Real2DArray& samples_array);

private:

  //
  //- Heading: Convenience functions
  //

  /// checks the return codes from LHS routines and aborts if an
  /// error is returned
  void check_error(const int& err_code, const char* err_source) const;

  //
  //- Heading: Data
  //

  /// type of sampling: random, lhs, incremental_lhs, or incremental_random
  String sampleType;

  /// mode of sample ranks I/O: IGNORE_RANKS, SET_RANKS, GET_RANKS, or
  /// SET_GET_RANKS
  short sampleRanksMode;

  /// flag for generating LHS report output
  bool reportFlag;
};


inline LHSDriver::LHSDriver() :
  sampleType("lhs"), sampleRanksMode(IGNORE_RANKS), reportFlag(true)
{}


inline LHSDriver::~LHSDriver()
{}


inline void LHSDriver::
initialize(const String& sample_type, short sample_ranks_mode, bool reports)
{
  sampleType      = sample_type;
  sampleRanksMode = sample_ranks_mode;
  reportFlag      = reports;
}


inline LHSDriver::
LHSDriver(const String& sample_type, short sample_ranks_mode, bool reports)
{ initialize(sample_type, sample_ranks_mode, reports); }


inline void LHSDriver::
check_error(const int&  err_code, const char* err_source) const
{
  if (err_code) {
    PCerr << "Error: code " << err_code << " returned from " << err_source 
	  << " in LHSDriver." << std::endl;
    abort_handler(-1);
  }
}

} // namespace Pecos

#endif

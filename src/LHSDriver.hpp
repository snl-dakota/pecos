/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
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
	    bool reportFlag = true);  ///< constructor
  ~LHSDriver(); ///< destructor

  //
  //- Heading: Member functions
  //

  /// populate data when not passed through ctor
  void initialize(const String& sample_type,
		  short sample_ranks_mode = IGNORE_RANKS,
		  bool reportFlag = true);

  /// set randomSeed
  void seed(int seed);
  /// return randomSeed
  int seed() const;

  /// set random number generator
  void rng(const String& unif_gen);
  // return name of uniform generator
  //String rng();

  /// reseed using a deterministic sequence
  void advance_seed_sequence();

  /// generates the desired set of parameter samples from within general
  /// user-specified probabilistic distributions (DistributionParams
  /// specification augmented with design and state variables)
  void generate_samples(const RealVector&   cd_l_bnds,
			const RealVector&   cd_u_bnds,
			const IntVector&    ddr_l_bnds,
			const IntVector&    ddr_u_bnds,
			const IntSetArray&  ddsi_values,
			const RealSetArray& ddsr_values,
			const RealVector&   cs_l_bnds,
			const RealVector&   cs_u_bnds,
			const IntVector&    dsr_l_bnds,
			const IntVector&    dsr_u_bnds,
			const IntSetArray&  dssi_values,
			const RealSetArray& dssr_values,
			const DistributionParams& dp, int num_samples,
			RealMatrix& samples_array, RealMatrix& rank_array);

  /// generates the desired set of parameter samples from within a
  /// DistributionParams specification
  void generate_samples(const DistributionParams& dp, int num_samples,
			RealMatrix& samples_array);

  /// generates the desired set of parameter samples from within
  /// uncorrelated normal distributions
  void generate_normal_samples(const RealVector& n_means,
			       const RealVector& n_std_devs,
			       const RealVector& n_l_bnds,
			       const RealVector& n_u_bnds, int num_samples,
			       RealMatrix& samples_array);

  /// generates the desired set of parameter samples from within
  /// uncorrelated uniform distributions
  void generate_uniform_samples(const RealVector& u_l_bnds,
				const RealVector& u_u_bnds, int num_samples,
				RealMatrix& samples_array);

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

  /// the current random number seed
  int randomSeed;

  /// mode of sample ranks I/O: IGNORE_RANKS, SET_RANKS, GET_RANKS, or
  /// SET_GET_RANKS
  short sampleRanksMode;

  /// flag for generating LHS report output
  bool reportFlag;

  /// for honoring advance_seed_sequence() calls
  int allowSeedAdvance; // bit 1 = first-time flag
		        // bit 2 = allow repeated seed update
};


inline LHSDriver::LHSDriver() : sampleType("lhs"), randomSeed(0),
  sampleRanksMode(IGNORE_RANKS), reportFlag(true), allowSeedAdvance(1)
{
#ifndef HAVE_LHS
  PCerr << "Error: LHSDriver not available as PECOS was configured without LHS."
	<< std::endl;
  abort_handler(-1);
#endif
}


inline LHSDriver::~LHSDriver()
{}


inline void LHSDriver::
initialize(const String& sample_type, short sample_ranks_mode, bool reports)
{
  sampleType      = sample_type;
  sampleRanksMode = sample_ranks_mode;
  reportFlag      = reports;
}


inline LHSDriver::LHSDriver(const String& sample_type, short sample_ranks_mode,
			    bool reports): allowSeedAdvance(1)
{
#ifndef HAVE_LHS
  PCerr << "Error: LHSDriver not available as PECOS was configured without LHS."
	<< std::endl;
  abort_handler(-1);
#endif
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
generate_samples(const DistributionParams& dp, int num_samples,
		 RealMatrix& samples_array)
{
  if (sampleRanksMode) {
    PCerr << "Error: generate_samples(DistributionParams&) does not support "
	  << "sample rank input/output." << std::endl;
    abort_handler(-1);
  }
  RealVector  empty_rv;  IntVector    empty_iv; RealMatrix empty_rm;
  IntSetArray empty_isa; RealSetArray empty_rsa;
  generate_samples(empty_rv, empty_rv, empty_iv, empty_iv, empty_isa, empty_rsa,
		   empty_rv, empty_rv, empty_iv, empty_iv, empty_isa, empty_rsa,
		   dp, num_samples, samples_array, empty_rm);
}


inline void LHSDriver::
generate_normal_samples(const RealVector& n_means, const RealVector& n_std_devs,
			const RealVector& n_l_bnds, const RealVector& n_u_bnds,
			int num_samples, RealMatrix& samples_array)
{
  if (sampleRanksMode) {
    PCerr << "Error: generate_normal_samples() does not support sample rank "
	  << "input/output." << std::endl;
    abort_handler(-1);
  }
  RealVector  empty_rv;  RealVectorArray empty_rva; IntVector empty_iv;
  RealMatrix  empty_rm;  RealSymMatrix   empty_rsm;
  IntSetArray empty_isa; RealSetArray    empty_rsa;
  DistributionParams dp(n_means, n_std_devs, n_l_bnds, n_u_bnds, empty_rv,
			empty_rv, empty_rv, empty_rv, empty_rv, empty_rv,
			empty_rv, empty_rv, empty_rv, empty_rv, empty_rv,
			empty_rv, empty_rv, empty_rv, empty_rv, empty_rv,
			empty_rv, empty_rv, empty_rv, empty_rv, empty_rv,
			empty_rv, empty_rv, empty_rv, empty_rv, empty_rv,
			empty_rv, empty_rva, empty_rv, empty_rv, empty_iv,
			empty_rv, empty_iv, empty_rv, empty_iv, empty_iv,
			empty_iv, empty_rva, empty_rva, empty_rva, empty_rsm);
  generate_samples(empty_rv, empty_rv, empty_iv, empty_iv, empty_isa, empty_rsa,
		   empty_rv, empty_rv, empty_iv, empty_iv, empty_isa, empty_rsa,
		   dp, num_samples, samples_array, empty_rm);
}


inline void LHSDriver::
generate_uniform_samples(const RealVector& u_l_bnds, const RealVector& u_u_bnds,
			 int num_samples, RealMatrix& samples_array)
{
  if (sampleRanksMode) {
    PCerr << "Error: generate_uniform_samples() does not support sample rank "
	  << "input/output." << std::endl;
    abort_handler(-1);
  }
  RealVector empty_rv;   RealVectorArray empty_rva; IntVector empty_iv;
  RealMatrix empty_rm;   RealSymMatrix   empty_rsm; 
  IntSetArray empty_isa; RealSetArray    empty_rsa;
  DistributionParams dp(empty_rv, empty_rv, empty_rv, empty_rv, empty_rv,
			empty_rv, empty_rv, empty_rv, empty_rv, empty_rv,
			empty_rv, u_l_bnds, u_u_bnds, empty_rv, empty_rv,
			empty_rv, empty_rv, empty_rv, empty_rv, empty_rv,
			empty_rv, empty_rv, empty_rv, empty_rv, empty_rv,
			empty_rv, empty_rv, empty_rv, empty_rv, empty_rv,
			empty_rv, empty_rva, empty_rv, empty_rv, empty_iv,
			empty_rv, empty_iv, empty_rv, empty_iv, empty_iv,
			empty_iv, empty_rva, empty_rva, empty_rva, empty_rsm);
  generate_samples(empty_rv, empty_rv, empty_iv, empty_iv, empty_isa, empty_rsa,
		   empty_rv, empty_rv, empty_iv, empty_iv, empty_isa, empty_rsa,
		   dp, num_samples, samples_array, empty_rm);
}


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

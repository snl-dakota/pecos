/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        SharedRegressPolyApproxData
//- Description:  Class for Multivariate Orthogonal Polynomial Approximations
//-               
//- Owner:        John Jakeman

#ifndef SHARED_REGRESS_ORTHOG_POLY_APPROX_DATA_HPP
#define SHARED_REGRESS_ORTHOG_POLY_APPROX_DATA_HPP

#include "SharedOrthogPolyApproxData.hpp"
#include "LinearSolver.hpp"

namespace Pecos {


/// Container class for various regression configuration options

/** The RegressionConfigOptions class provides a simple container class
    for regression configuration options related to ... */

class RegressionConfigOptions
{
  //
  //- Heading: Friends
  //

  //friend class PolynomialApproximation;

public:

  /// default constructor
  RegressionConfigOptions();
  /// constructor
  RegressionConfigOptions(bool cv, int seed, const RealVector& noise_tols,
			  Real l2_penalty, bool normalize_cv,
			  unsigned short init_lev, unsigned short growth_fact,
			  unsigned short num_advance);
  /// copy constructor
  RegressionConfigOptions(const RegressionConfigOptions& rc_options);
  /// destructor
  ~RegressionConfigOptions();

//private:

  /// flag for use of automatic cross-validation for parameter
  /// selection in regression approaches
  bool crossValidation;
  /// random seed for data fold definition within cross validation
  int randomSeed;

  /// noise tolerance(s) for compressed sensing algorithms; vector form
  /// used in cross-validation
  RealVector noiseTols;
  /// the L2 penalty parameter for LASSO (elastic net variant)
  Real l2Penalty;

  /// flag indicating whether to scale the CV error estimates by the size
  /// of the candidate basis expansion within adapted basis selection
  bool normalizeCV;

  /// initial sparse grid level that provides the starting point for basis
  /// adaptation in ADAPTED_BASIS_GENERALIZED mode
  unsigned short initSGLevel;
  /// a scalar growth factor for defining the tpMultiIndex contribution
  /// from a particular trial index set for basis adaptation in
  /// ADAPTED_BASIS_GENERALIZED mode
  unsigned short multiIndexGrowthFactor;

  /// number of front expansions per iteration for basis adaptation in
  /// ADAPTED_BASIS_EXPANDING_FRONT mode
  unsigned short numAdvancements;
};


inline RegressionConfigOptions::RegressionConfigOptions():
  crossValidation(false), randomSeed(0), l2Penalty(0.), normalizeCV(false),
  initSGLevel(0), multiIndexGrowthFactor(2), numAdvancements(3)
{ }


inline RegressionConfigOptions::
RegressionConfigOptions(bool cv, int seed, const RealVector& noise_tols,
			Real l2_penalty, bool normalize_cv,
			unsigned short init_lev, unsigned short growth_fact,
			unsigned short num_advance):
  crossValidation(cv), randomSeed(seed), noiseTols(noise_tols),
  l2Penalty(l2_penalty), normalizeCV(normalize_cv), initSGLevel(init_lev),
  multiIndexGrowthFactor(growth_fact), numAdvancements(num_advance)
{ }


inline RegressionConfigOptions::
RegressionConfigOptions(const RegressionConfigOptions& rc_options):
  crossValidation(rc_options.crossValidation),
  randomSeed(rc_options.randomSeed), noiseTols(rc_options.noiseTols),
  l2Penalty(rc_options.l2Penalty), normalizeCV(rc_options.normalizeCV),
  initSGLevel(rc_options.initSGLevel),
  multiIndexGrowthFactor(rc_options.multiIndexGrowthFactor),
  numAdvancements(rc_options.numAdvancements)
{ }


inline RegressionConfigOptions::~RegressionConfigOptions()
{ }


/// Derived approximation class for multivariate orthogonal polynomial
/// approximation with coefficient estimation via regression.

/** The SharedRegressOrthogPolyApproxData class provides a global
    approximation based on multivariate orthogonal polynomials, where
    the coefficients are computed using regression approaches such as
    least squares (L2) or compressed sensing (L1).  It is used
    primarily for polynomial chaos expansion aproaches to UQ. */

class SharedRegressOrthogPolyApproxData: public SharedOrthogPolyApproxData
{
  //
  //- Heading: Friends
  //

  friend class RegressOrthogPolyApproximation;

public:

  //
  //- Heading: Constructor and destructor
  //

  /// lightweight constructor
  SharedRegressOrthogPolyApproxData(short basis_type,
				    const UShortArray& approx_order,
				    size_t num_vars);
  /// full constructor
  SharedRegressOrthogPolyApproxData(short basis_type,
				    const UShortArray& approx_order,
				    size_t num_vars,
				    const ExpansionConfigOptions& ec_options,
				    const BasisConfigOptions& bc_options);
  /// destructor
  ~SharedRegressOrthogPolyApproxData();

  //
  //- Heading: Member functions
  //

  /// helper function for incrementing that is modular on trial set and
  /// multi-index
  void increment_trial_set(const UShortArray& trial_set,
			   UShort2DArray& aggregated_mi);

  /// update cvErrorRef and bestExpTerms after improvement in solution
  void update_reference(Real cv_error, const UShort2DArray& aggregated_mi);
  /// update cvErrorRef and bestExpTerms after improvement in solution
  void restore_best_solution(UShort2DArray& aggregated_mi);

  /// set regressConfigOptions
  void configuration_options(const RegressionConfigOptions& rc_options);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void allocate_data();
  //void increment_data();
  //void decrement_data();
  //void restore_data();
  //void finalize_data();

private:

  //
  //- Heading: Member functions
  //

  /// update shared approxOrder based on settings computed for one of the QoI
  void update_approx_order(unsigned short new_order);

  /// pack polynomial contributions to Psi matrix for regression
  void pack_polynomial_data(const RealVector& c_vars, const UShortArray& mi,
			    bool add_val,  double* pack_val,  size_t& pv_cntr,
			    bool add_grad, double* pack_grad, size_t& pg_cntr);
  /// pack response contributions to RHS for regression
  void pack_response_data(const SurrogateDataResp& sdr,
			  bool add_val,  double* pack_val,  size_t& pv_cntr,
			  bool add_grad, double* pack_grad, size_t& pg_cntr);

  //
  //- Heading: Data
  //

  /// container for regression configuration options
  RegressionConfigOptions regressConfigOptions;

  /// Wrapper class that is used to solve regression problems
  CompressedSensingTool CSTool;

  /// lower matrix factor in factorization
  RealMatrix lowerFactor;
  /// upper matrix factor in factorization
  RealMatrix upperFactor;
  /// pivoting history of block-LU factorization
  RealMatrix pivotHistory;
  /// pivoting vector
  IntVector pivotVect;

  /// sparse grid driver for adapting a CS candidate basis using greedy
  /// adaptation via the generalized sparse grid algorithm; it's state
  /// is reset for each response QoI
  CombinedSparseGridDriver csgDriver;
  /// the cross validation error reference point for adapting a CS
  /// candidate basis; it's state is reset for each response QoI
  Real cvErrorRef;
  /// size of expansion that corresponds to the best solution identified
  /// (may not be the last solution)
  size_t bestExpTerms;
};


inline SharedRegressOrthogPolyApproxData::
SharedRegressOrthogPolyApproxData(short basis_type,
				  const UShortArray& approx_order,
				  size_t num_vars):
  SharedOrthogPolyApproxData(basis_type, approx_order, num_vars)
{ }


inline SharedRegressOrthogPolyApproxData::
SharedRegressOrthogPolyApproxData(short basis_type,
				  const UShortArray& approx_order,
				  size_t num_vars,
				  const ExpansionConfigOptions& ec_options,
				  const BasisConfigOptions& bc_options):
  SharedOrthogPolyApproxData(basis_type, approx_order, num_vars,
			     ec_options, bc_options)
{ }


inline SharedRegressOrthogPolyApproxData::~SharedRegressOrthogPolyApproxData()
{ }


inline void SharedRegressOrthogPolyApproxData::
configuration_options(const RegressionConfigOptions& rc_options)
{ regressConfigOptions = rc_options; }


inline void SharedRegressOrthogPolyApproxData::
update_reference(Real cv_error, const UShort2DArray& aggregated_mi)
{
  cvErrorRef   = cv_error;
  bestExpTerms = aggregated_mi.size();
}

} // namespace Pecos

#endif

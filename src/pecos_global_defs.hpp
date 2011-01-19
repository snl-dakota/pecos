/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_GLOBAL_DEFS_H
#define PECOS_GLOBAL_DEFS_H

#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <cfloat>  // for DBL_MIN, DBL_MAX
#include <cmath>
#include <cstdlib>


namespace Pecos {

// --------------
// Special values
// --------------
/// the value for PI used in various numerical routines
const double PI = boost::math::constants::pi<double>();

/// special value returned by index() when entry not found
const size_t _NPOS = ~(size_t)0; // one's complement

/// used in ostream data output functions
const int WRITE_PRECISION = 10;

// define special values for ranVarTypesX/U
enum { DESIGN, STD_NORMAL, NORMAL, BOUNDED_NORMAL, LOGNORMAL, BOUNDED_LOGNORMAL,
       STD_UNIFORM, UNIFORM, LOGUNIFORM, TRIANGULAR, STD_EXPONENTIAL,
       EXPONENTIAL, STD_BETA, BETA, STD_GAMMA, GAMMA, GUMBEL, FRECHET, WEIBULL,
       HISTOGRAM_BIN, HISTOGRAM_POINT, STOCHASTIC_EXPANSION, INTERVAL, STATE };

// define special values for secondaryACVarMapTargets/secondaryADVarMapTargets
enum { NO_TARGET, CDV_LWR_BND, CDV_UPR_BND, DDRIV_LWR_BND, DDRIV_UPR_BND,
       N_MEAN, N_STD_DEV, N_LWR_BND, N_UPR_BND, LN_MEAN, LN_STD_DEV,
       LN_LAMBDA, LN_ZETA, LN_ERR_FACT, LN_LWR_BND, LN_UPR_BND,
       U_LWR_BND, U_UPR_BND, LU_LWR_BND, LU_UPR_BND,
       T_MODE, T_LWR_BND, T_UPR_BND, E_BETA,
       BE_ALPHA, BE_BETA, BE_LWR_BND, BE_UPR_BND, GA_ALPHA, GA_BETA,
       GU_ALPHA, GU_BETA, F_ALPHA, F_BETA, W_ALPHA, W_BETA,
       P_LAMBDA, BI_P_PER_TRIAL, BI_TRIALS, NBI_P_PER_TRIAL, NBI_TRIALS,
       GE_P_PER_TRIAL, HGE_TOT_POP, HGE_SEL_POP, HGE_FAILED,
       CSV_LWR_BND, CSV_UPR_BND, DSRIV_LWR_BND, DSRIV_UPR_BND };

/// derived basis polynomial types (orthogonal polynomial order follows
/// uncertain variable spec order of normal, uniform, exponential, beta, gamma)
enum { /* NO_POLY, */ HERMITE, LEGENDRE, LAGUERRE, JACOBI, GENERALIZED_LAGUERRE,
       CHEBYSHEV, NUMERICALLY_GENERATED, LAGRANGE };

/// integration rules within VPISparseGrid
enum { CLENSHAW_CURTIS=1, FEJER2, GAUSS_PATTERSON, GAUSS_LEGENDRE,
       GAUSS_HERMITE, GEN_GAUSS_HERMITE, GAUSS_LAGUERRE, GEN_GAUSS_LAGUERRE,
       GAUSS_JACOBI, GOLUB_WELSCH, GENZ_KEISTER };

/// growth rules within VPISparseGrid
enum { DEFAULT_GROWTH, SLOW_LINEAR, SLOW_LINEAR_ODD, MODERATE_LINEAR,
       SLOW_EXPONENTIAL, MODERATE_EXPONENTIAL, FULL_EXPONENTIAL };

/// options for synchronizing linear and exponential growth rule settings
enum { SLOW_RESTRICTED_GROWTH, MODERATE_RESTRICTED_GROWTH,
       UNRESTRICTED_GROWTH };

/// solution approaches for calculating the polynomial basis coefficients
/// (options for ConfigurationOptions::expCoeffsSolnApproach)
enum { QUADRATURE, CUBATURE, SPARSE_GRID, REGRESSION, SAMPLING };
/// options for ConfigurationOptions::nestingOverride (inactive)
enum { NO_OVERRIDE=0, NESTED, NON_NESTED };
/// options for ConfigurationOptions::refinementType (inactive)
enum { NO_REFINEMENT=0, UNIFORM_P_REFINEMENT, ADAPTIVE_P_REFINEMENT };
/// options for ConfigurationOptions::refinementControl
enum { DEFAULT_CONTROL=0, TOTAL_SOBOL, SPECTRAL_DECAY, GENERALIZED_SPARSE };
/// options for ConfigurationOptions::vbdSetting
enum { NO_VBD=0, UNIVARIATE_VBD, ALL_VBD };

///
enum { LINEAR_EQUIDISTANT, LINEAR, QUADRATIC_EQUIDISTANT, QUADRATIC,
       CUBIC_EQUIDISTANT, CUBIC };


// ----------------
// Standard streams
// ----------------
#define PCout std::cout
#define PCerr std::cerr


// --------------
// Global objects
// --------------
/// Dummy struct for overloading letter-envelope constructors.
/** BaseConstructor is used to overload the constructor for the base class
    portion of letter objects.  It avoids infinite recursion (Coplien p.139)
    in the letter-envelope idiom by preventing the letter from instantiating
    another envelope.  Putting this struct here avoids circular dependencies. */
struct BaseConstructor {
  BaseConstructor(int = 0) {} ///< C++ structs can have constructors
};


// ----------------
// Global functions
// ----------------

/// global function which handles serial or parallel aborts
void abort_handler(int code);


inline void abort_handler(int code)
{ std::exit(code); } // for now, prior to use of MPI

} // namespace Pecos

#endif // PECOS_GLOBAL_DEFS_H

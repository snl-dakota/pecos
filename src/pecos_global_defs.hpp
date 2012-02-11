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

// define special values for vector/matrix data copying modes
enum { DEFAULT_COPY=0, SHALLOW_COPY, DEEP_COPY };

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

/// derived basis approximation types
enum { NO_BASIS, FOURIER_BASIS, EIGEN_BASIS,
       GLOBAL_NODAL_INTERPOLATION_POLYNOMIAL,
       PIECEWISE_NODAL_INTERPOLATION_POLYNOMIAL,
       GLOBAL_HIERARCHICAL_INTERPOLATION_POLYNOMIAL,
       PIECEWISE_HIERARCHICAL_INTERPOLATION_POLYNOMIAL,
       GLOBAL_ORTHOGONAL_POLYNOMIAL, PIECEWISE_ORTHOGONAL_POLYNOMIAL };

/// derived basis polynomial types (orthogonal polynomial order follows
/// uncertain variable spec order of normal, uniform, exponential, beta, gamma)
enum { NO_POLY=0, HERMITE_ORTHOG, LEGENDRE_ORTHOG, LAGUERRE_ORTHOG,
       JACOBI_ORTHOG, GEN_LAGUERRE_ORTHOG, CHEBYSHEV_ORTHOG, NUM_GEN_ORTHOG,
       LAGRANGE_INTERP, HERMITE_INTERP, PIECEWISE_LINEAR_INTERP,
       PIECEWISE_QUADRATIC_INTERP, PIECEWISE_CUBIC_INTERP };

/// integration rules within VPISparseGrid (1-12: CC through User-closed)
/// and beyond (GOLUB_WELSCH, NEWTON_COTES)
enum { NO_RULE=0, CLENSHAW_CURTIS, FEJER2, GAUSS_PATTERSON, GAUSS_LEGENDRE,
       GAUSS_HERMITE, GEN_GAUSS_HERMITE, GAUSS_LAGUERRE, GEN_GAUSS_LAGUERRE,
       GAUSS_JACOBI, GENZ_KEISTER, /*USER_OPEN, USER_CLOSED,*/ GOLUB_WELSCH,
       NEWTON_COTES };

// growth rules within VPISparseGrid
//enum { DEFAULT_GROWTH, SLOW_LINEAR, SLOW_LINEAR_ODD, MODERATE_LINEAR,
//       SLOW_EXPONENTIAL, MODERATE_EXPONENTIAL, FULL_EXPONENTIAL };

/// options for synchronizing linear and exponential growth rule settings
/// (consistent with slow/moderate/full growth for new level_to_growth_*
/// functions in sandia_rules.cpp)
enum { SLOW_RESTRICTED_GROWTH, MODERATE_RESTRICTED_GROWTH,
       UNRESTRICTED_GROWTH };

/// solution approaches for calculating the polynomial basis coefficients
/// (options for ExpansionConfigOptions::expCoeffsSolnApproach)
enum { QUADRATURE, CUBATURE, COMBINED_SPARSE_GRID, HIERARCHICAL_SPARSE_GRID,
       LOCAL_REFINABLE, REGRESSION, SAMPLING };
/// options for BasisConfigOptions::nestingOverride (inactive)
enum { NO_NESTING_OVERRIDE=0, NESTED, NON_NESTED };
/// options for overriding the default growth restriction policy
enum { NO_GROWTH_OVERRIDE=0, RESTRICTED, UNRESTRICTED };
/// options for ExpansionConfigOptions::refinementType (inactive)
enum { NO_REFINEMENT=0, P_REFINEMENT, H_REFINEMENT };
/// options for ExpansionConfigOptions::refinementControl
enum { NO_CONTROL=0, UNIFORM_CONTROL, LOCAL_ADAPTIVE_CONTROL,
       DIMENSION_ADAPTIVE_CONTROL_SOBOL, DIMENSION_ADAPTIVE_CONTROL_DECAY,
       DIMENSION_ADAPTIVE_CONTROL_GENERALIZED };
/// options for ExpansionConfigOptions::vbdSetting
enum { NO_VBD=0, UNIVARIATE_VBD, ALL_VBD };

/// options for local basis functions within PiecewiseInterpPolynomial
enum { LINEAR_EQUIDISTANT, LINEAR, QUADRATIC_EQUIDISTANT, QUADRATIC,
       CUBIC_EQUIDISTANT, CUBIC };

/// special values for polynomial expansion combination
enum { NO_COMBINE=0,  ADD_COMBINE, MULT_COMBINE, ADD_MULT_COMBINE };

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


/** Templatized abort_handler_t method that allows for convenient return from 
    methods that otherwise have no sensible return from error clauses.  Usage:
    MyType& method() { return abort_handler<MyType&>(-1); } */
template <typename T>
T abort_handler_t(int code)
{
  abort_handler(code);
  throw code;
}

} // namespace Pecos

#endif // PECOS_GLOBAL_DEFS_H

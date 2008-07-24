/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_GLOBAL_DEFS_H
#define PECOS_GLOBAL_DEFS_H

#include "pecos_data_types.h"
#include <iostream>
#include <cstdlib>


namespace Pecos {

// --------------
// Special values
// --------------
/// special value returned by index() when entry not found
const size_t _NPOS = ~(size_t)0; // one's complement

// define special values for ranVarTypesX/U
enum { DESIGN, NORMAL, BOUNDED_NORMAL, LOGNORMAL, BOUNDED_LOGNORMAL, UNIFORM,
       LOGUNIFORM, TRIANGULAR, EXPONENTIAL, BETA, GAMMA, GUMBEL, FRECHET,
       WEIBULL, STATE };

// define special values for secondaryACVarMapTargets/secondaryADVarMapTargets
enum { NO_TARGET, CDV_LWR_BND, CDV_UPR_BND, DDV_LWR_BND, DDV_UPR_BND,
       N_MEAN, N_STD_DEV, N_LWR_BND, N_UPR_BND,
       LN_MEAN, LN_STD_DEV, LN_ERR_FACT, LN_LWR_BND, LN_UPR_BND,
       U_LWR_BND, U_UPR_BND, LU_LWR_BND, LU_UPR_BND,
       T_MODE, T_LWR_BND, T_UPR_BND, E_BETA,
       B_ALPHA, B_BETA, B_LWR_BND, B_UPR_BND, GA_ALPHA, GA_BETA,
       GU_ALPHA, GU_BETA, F_ALPHA, F_BETA, W_ALPHA, W_BETA,
       CSV_LWR_BND, CSV_UPR_BND, DSV_LWR_BND, DSV_UPR_BND };


// ----------------
// Standard streams
// ----------------
#define Cout std::cout
#define Cerr std::cerr


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

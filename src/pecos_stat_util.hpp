/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_STAT_UTIL_HPP
#define PECOS_STAT_UTIL_HPP

#include "pecos_data_types.hpp"
#include "BoundedNormalRandomVariable.hpp"
#include "BoundedLognormalRandomVariable.hpp"
#include "LoguniformRandomVariable.hpp"
#include "TriangularRandomVariable.hpp"
#include "BetaRandomVariable.hpp"
#include "GammaRandomVariable.hpp"
#include "GumbelRandomVariable.hpp"
#include "FrechetRandomVariable.hpp"
#include "WeibullRandomVariable.hpp"
#include "HistogramBinRandomVariable.hpp"
#include "InvGammaRandomVariable.hpp"
#include "PoissonRandomVariable.hpp"
#include "BinomialRandomVariable.hpp"
#include "NegBinomialRandomVariable.hpp"
#include "GeometricRandomVariable.hpp"
#include "HypergeometricRandomVariable.hpp"
#include "RangeVariable.hpp"
#include "SetVariable.hpp"
#include "IntervalRandomVariable.hpp"
#include "DiscreteSetRandomVariable.hpp"


namespace Pecos {

inline Real gamma_function(Real x)
{ return bmth::tgamma(x); }

} // namespace Pecos

#endif

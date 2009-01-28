/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_DATA_TYPES_H
#define PECOS_DATA_TYPES_H

#ifdef HAVE_CONFIG_H
#include "pecos_config.h"
#endif /* HAVE_CONFIG_H */
#include "pecos_global_defs.hpp"

#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_SerialSpdDenseSolver.hpp"

#ifdef HAVE_BOOST
// WJB - ToDo: investigate error in boost/math/tools/traits.hpp with SunProCC
//#include <boost/math/distributions.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/weibull.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/policies/policy.hpp>
namespace bmth = boost::math;
namespace bmp  = bmth::policies;
#endif

#include <complex>


namespace Pecos {

// avoid problems with circular dependencies by using fwd declarations
//class BasisFunction;


// -----------------------------------
// Aliases for fundamental data types:
// -----------------------------------
typedef double Real;

// --------
// Strings:
// --------
typedef std::string String;

// -----------------------------------
// Non-default boost math/policy types
// -----------------------------------
#ifdef HAVE_BOOST
typedef bmth::
  normal_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  normal_dist;
typedef bmth::
  gamma_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  gamma_dist;
typedef bmth::
  exponential_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  exponential_dist;
typedef bmth::
  beta_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  beta_dist;
typedef bmth::
  weibull_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  weibull_dist;
typedef bmth::
  chi_squared_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  chi_squared_dist;
typedef bmth::
  students_t_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  students_t_dist;
typedef bmth::
  fisher_f_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  fisher_f_dist;
#endif

// --------------------------------
// Numerical arrays (serial dense):
// --------------------------------
typedef Teuchos::SerialDenseVector<int, Real>                RealVector;
typedef Teuchos::SerialDenseVector<int, int>                 IntVector;
typedef Teuchos::SerialDenseVector<int, std::complex<Real> > ComplexVector;
typedef Teuchos::SerialDenseMatrix<int, Real>                RealMatrix;
typedef Teuchos::SerialSymDenseMatrix<int, Real>             RealSymMatrix;


// ---------------------------------
// Numerical solvers (serial dense):
// ---------------------------------
typedef Teuchos::SerialDenseSolver<int, Real>    RealSolver;
typedef Teuchos::SerialSpdDenseSolver<int, Real> RealSpdSolver;


// ---------------------------------------
// Admin/bookkeeping arrays (serial only):
// ---------------------------------------
typedef std::deque<bool>            BoolDeque; // See Meyers' Effective STL, #18
typedef std::vector<Real>           RealArray;
typedef std::vector<int>            IntArray;
typedef std::vector<unsigned int>   UIntArray;
typedef std::vector<short>          ShortArray;
typedef std::vector<unsigned short> UShortArray;
typedef std::vector<size_t>         SizetArray;
typedef std::vector<std::complex<Real> >    ComplexArray;
typedef std::vector<std::pair<Real, Real> > RealPairArray;
typedef std::vector<String>         StringArray;
typedef std::vector<RealArray>      Real2DArray;
typedef std::vector<RealVector>     RealVectorArray;
typedef std::vector<RealMatrix>     RealMatrixArray;
typedef std::vector<RealSymMatrix>  RealSymMatrixArray;
//typedef std::vector<BasisFunction>  BasisFunctionArray;

typedef std::set<int>               IntSet;
typedef std::set<Real>              RealSet;
typedef std::map<int, short>        IntShortMap;
typedef std::map<int, int>          IntIntMap;
typedef std::map<int, RealVector>   IntRealVectorMap;


// ---------
// Iterators
// ---------
typedef IntSet::iterator            ISIter;
typedef IntSet::const_iterator      ISCIter;
typedef IntShortMap::iterator       IntShMIter;
typedef IntIntMap::iterator         IntIntMIter;
typedef IntIntMap::const_iterator   IntIntMCIter;


} // namespace Pecos

#endif // PECOS_DATA_TYPES_H

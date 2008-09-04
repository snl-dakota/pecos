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

#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_SerialSpdDenseSolver.hpp"


namespace Pecos {

// -----------------------------------
// Aliases for fundamental data types:
// -----------------------------------
typedef double Real;


// --------
// Strings:
// --------
typedef std::string String;


// --------------------------------
// Numerical arrays (serial dense):
// --------------------------------
typedef Teuchos::SerialDenseVector<int, Real>    RealVector;
typedef Teuchos::SerialDenseVector<int, int>     IntVector;
typedef Teuchos::SerialDenseMatrix<int, Real>    RealMatrix;
typedef Teuchos::SerialSymDenseMatrix<int, Real> RealSymMatrix;


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
typedef std::vector<String>         StringArray;
typedef std::vector<RealVector>     RealVectorArray;
typedef std::vector<RealMatrix>     RealMatrixArray;
typedef std::vector<RealSymMatrix>  RealSymMatrixArray;

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

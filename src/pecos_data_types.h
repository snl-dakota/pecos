/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

#ifndef PECOS_DATA_TYPES_H
#define PECOS_DATA_TYPES_H

#include "Teuchos_SerialDenseVector.h"
#include "Teuchos_SerialDenseMatrix.h"


namespace Pecos {

// -----------------------------------
// Aliases for fundamental data types:
// -----------------------------------
typedef double Real;


// --------
// Strings:
// --------
typedef std::string String;


// -----------------------------------
// Numerical arrays (serial/parallel):
// -----------------------------------
typedef Teuchos::SerialDenseVector<size_t, Real>    RealVector;
typedef Teuchos::SerialDenseVector<size_t, int>     IntVector;
typedef Teuchos::SerialDenseMatrix<size_t, Real>    RealMatrix;
typedef Teuchos::SerialSymDenseMatrix<size_t, Real> RealSymMatrix;


// ---------------------------------------
// Admin/bookkeeping arrays (serial only):
// ---------------------------------------
typedef std::deque<bool>           BoolDeque; // See Effective STL (Meyers), #18
typedef std::vector<Real>          RealArray;
typedef std::vector<int>           IntArray;
typedef std::vector<short>         ShortArray;
typedef std::vector<size_t>        SizetArray;
typedef std::vector<String>        StringArray;
typedef std::vector<RealVector>    RealVectorArray;
typedef std::vector<RealMatrix>    RealMatrixArray;
typedef std::vector<RealSymMatrix> RealSymMatrixArray;

typedef std::set<int>              IntSet;
typedef std::set<Real>             RealSet;
typedef std::map<int, short>       IntShortMap;
typedef std::map<int, int>         IntIntMap;
typedef std::map<int, RealVector>  IntRealVectorMap;


// ---------
// Iterators
// ---------
typedef IntSet::iterator           ISIter;
typedef IntSet::const_iterator     ISCIter;
typedef IntShortMap::iterator      IntShMIter;
typedef IntIntMap::iterator        IntIntMIter;
typedef IntIntMap::const_iterator  IntIntMCIter;

} // namespace Pecos

#endif // DATA_TYPES_H

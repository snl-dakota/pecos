/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
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

#include "boost/multi_array.hpp"
#include <algorithm>  // for std::find
#include <complex>
#include <map>
#include <set>
#include <string>
#include <vector>


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
typedef std::vector<IntSet>         IntSetArray;
typedef std::vector<RealSet>        RealSetArray;
typedef std::map<int, short>        IntShortMap;
typedef std::map<int, int>          IntIntMap;
typedef std::map<int, RealVector>   IntRealVectorMap;

typedef boost::multi_array_types::index_range     idx_range;
typedef boost::multi_array<unsigned int, 1>       UIntMultiArray;
typedef UIntMultiArray::array_view<1>::type       UIntMultiArrayView;
typedef UIntMultiArray::const_array_view<1>::type UIntMultiArrayConstView;

// ---------
// Iterators
// ---------
typedef IntSet::iterator            ISIter;
typedef IntSet::const_iterator      ISCIter;
typedef RealSet::iterator           RSIter;
typedef RealSet::const_iterator     RSCIter;
typedef IntShortMap::iterator       IntShMIter;
typedef IntIntMap::iterator         IntIntMIter;
typedef IntIntMap::const_iterator   IntIntMCIter;


/// equality operator for UIntArray and UIntMultiArrayConstView
inline bool operator==(const UIntArray& ua, UIntMultiArrayConstView umav)
{
  // Check for equality in array lengths
  size_t len = ua.size();
  if ( umav.size() != len )
    return false;

  // Check each unsigned integer
  size_t i;
  for (i=0; i<len; i++)
    if ( umav[i] != ua[i] )
      return false;

  return true;
}


template <typename PecosContainerType>
inline typename PecosContainerType::difference_type
find_index(const PecosContainerType& v,
	   const typename PecosContainerType::value_type& val)
{
  typename PecosContainerType::const_iterator iter
    = std::find(v.begin(), v.end(), val);
  return iter != v.end() ? iter - v.begin() : _NPOS;
}


/// copy Teuchos::SerialDenseVector<OrdinalType, ScalarType> to
/// std::vector<ScalarType>
template <typename OrdinalType, typename ScalarType> 
void copy_data(const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& sdv,
	       std::vector<ScalarType>& v)
{
  OrdinalType size_sdv = sdv.length();
  if (size_sdv != v.size())
    v.resize(size_sdv);
  for (OrdinalType i=0; i<size_sdv; ++i)
    v[i] = sdv[i];
}


/// copy std::vector<ScalarType> to
/// Teuchos::SerialDenseVector<OrdinalType, ScalarType>
template <typename OrdinalType, typename ScalarType> 
void copy_data(const std::vector<ScalarType>& v,
	       Teuchos::SerialDenseVector<OrdinalType, ScalarType>& sdv)
{
  size_t size_v = v.size();
  if (sdv.length() != size_v)
    sdv.sizeUninitialized(size_v);
  for (OrdinalType i=0; i<size_v; ++i)
    sdv[i] = v[i];
}


/// global std::ostream insertion operator for std::vector
template <class T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& data)
{
  size_t i=0, len = data.size();
  for (i=0; i<len; ++i)
    s << "                     " << data[i] << '\n';
  return s;
}

} // namespace Pecos

#endif // PECOS_DATA_TYPES_H

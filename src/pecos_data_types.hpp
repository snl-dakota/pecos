/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_DATA_TYPES_H
#define PECOS_DATA_TYPES_H

#if defined(HAVE_CONFIG_H) && !defined(DISABLE_DAKOTA_CONFIG_H)
  // HAVE_CONFIG_H is STILL set in Dakota/src (EVEN IN THE CMAKE BUILD!) so
  // use a "disable config header" conditional to help manage the transition
  #include "pecos_config.h"
#endif // HAVE_CONFIG_H

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
class SurrogateDataVars;
class SurrogateDataResp;

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
typedef std::vector<RealArray>      Real2DArray;
typedef std::vector<Real2DArray>    Real3DArray;
typedef std::vector<int>            IntArray;
typedef std::vector<IntArray>       Int2DArray;
typedef std::vector<short>          ShortArray;
typedef std::vector<unsigned short> UShortArray;
typedef std::vector<UShortArray>    UShort2DArray;
typedef std::vector<UShort2DArray>  UShort3DArray;
typedef std::vector<size_t>         SizetArray;
typedef std::vector<SizetArray>     Sizet2DArray;
typedef std::list<size_t>           SizetList;
typedef std::vector<std::complex<Real> >    ComplexArray;
typedef std::vector<std::pair<Real, Real> > RealPairArray;
typedef std::vector<String>         StringArray;
typedef std::vector<RealVector>     RealVectorArray;
typedef std::vector<RealMatrix>     RealMatrixArray;
typedef std::vector<RealSymMatrix>  RealSymMatrixArray;

//typedef std::vector<BasisFunction>  BasisFunctionArray;
typedef std::vector<SurrogateDataVars> SDVArray;
typedef std::vector<SurrogateDataResp> SDRArray;
typedef std::vector<SDVArray>          SDV2DArray;
typedef std::vector<SDRArray>          SDR2DArray;

typedef std::set<int>                  IntSet;
typedef std::multiset<unsigned short>  UShortMultiSet;
typedef std::set<Real>                 RealSet;
typedef std::vector<IntSet>            IntSetArray;
typedef std::vector<RealSet>           RealSetArray;
typedef std::map<int, short>           IntShortMap;
typedef std::map<int, int>             IntIntMap;
typedef std::map<int, RealVector>      IntRealVectorMap;
typedef std::map<UShortMultiSet, Real> UShortMultiSetRealMap;

typedef boost::multi_array_types::index_range      idx_range;
typedef boost::multi_array<size_t, 1>              SizetMultiArray;
typedef SizetMultiArray::array_view<1>::type       SizetMultiArrayView;
typedef SizetMultiArray::const_array_view<1>::type SizetMultiArrayConstView;

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


/// equality operator for SizetArray and SizetMultiArrayConstView
inline bool operator==(const SizetArray& sa, SizetMultiArrayConstView smav)
{
  // Check for equality in array lengths
  size_t i, len = sa.size();
  if ( smav.size() != len )
    return false;

  // Check each size_t
  for (i=0; i<len; i++)
    if ( smav[i] != sa[i] )
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
  return (iter == v.end()) ? _NPOS : std::distance(v.begin(), iter);
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


/// copy Teuchos::SerialDenseVector<OrdinalType, ScalarType> to same
/// (used in place of operator= when a deep copy is required regardless
/// of Teuchos DataAccess mode)
template <typename OrdinalType, typename ScalarType> 
void copy_data(const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& sdv1,
	       Teuchos::SerialDenseVector<OrdinalType, ScalarType>& sdv2)
{
  OrdinalType size_sdv1 = sdv1.length();
  if (size_sdv1 != sdv2.length())
    sdv2.sizeUninitialized(size_sdv1);
  for (OrdinalType i=0; i<size_sdv1; ++i)
    sdv2[i] = sdv1[i];
}


/// copy Teuchos::SerialSymDenseMatrix<OrdinalType, ScalarType> to same
/// (used in place of operator= when a deep copy is required regardless
/// of Teuchos DataAccess mode)
template <typename OrdinalType, typename ScalarType> 
void copy_data(const
	       Teuchos::SerialSymDenseMatrix<OrdinalType, ScalarType>& ssdm1,
	       Teuchos::SerialSymDenseMatrix<OrdinalType, ScalarType>& ssdm2)
{
  OrdinalType size_ssdm1 = ssdm1.numRows();
  if (size_ssdm1 != ssdm2.numRows())
    ssdm2.shapeUninitialized(size_ssdm1);
  ssdm2.assign(ssdm1); // copies values
}


/// copy ScalarType* to ScalarType*
template <typename OrdinalType, typename ScalarType> 
void copy_data(const ScalarType* ptr1, const OrdinalType ptr_len,
	       ScalarType* ptr2)
{
  for (OrdinalType i=0; i<ptr_len; ++i)
    ptr2[i] = ptr1[i];
}


/// copy ScalarType* to Teuchos::SerialDenseVector<OrdinalType, ScalarType>
template <typename OrdinalType, typename ScalarType> 
void copy_data(const ScalarType* ptr, const OrdinalType ptr_len,
	       Teuchos::SerialDenseVector<OrdinalType, ScalarType>& sdv)
{
  if (sdv.length() != ptr_len)
    sdv.sizeUninitialized(ptr_len);
  for (OrdinalType i=0; i<ptr_len; ++i)
    sdv[i] = ptr[i];
}


/// copy ScalarType* to std::deque<ScalarType>
template <typename OrdinalType, typename ScalarType> 
void copy_data(const ScalarType* ptr, const OrdinalType ptr_len,
	       std::deque<ScalarType>& deq)
{
  if (deq.size() != ptr_len)
    deq.resize(ptr_len);
  for (OrdinalType i=0; i<ptr_len; ++i)
    deq[i] = ptr[i];
}


/// copy std::deque<ScalarType> to ScalarType*
template <typename OrdinalType, typename ScalarType> 
void copy_data(const std::deque<ScalarType>& deq, ScalarType* ptr,
	       const OrdinalType ptr_len)
{
  for (OrdinalType i=0; i<ptr_len; ++i)
    ptr[i] = deq[i];
}


/// copy Teuchos::SerialDenseVector<OrdinalType, ScalarType> to
/// ith row of Teuchos::SerialDenseMatrix<OrdinalType, ScalarType>
template <typename OrdinalType, typename ScalarType> 
void copy_row(const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& sdv,
	      Teuchos::SerialDenseMatrix<OrdinalType, ScalarType>& sdm,
	      OrdinalType row)
{
  OrdinalType i, len = sdv.length();
  //if (sdm.numCols() != len)
  //  PCerr << std::endl;
  for (i=0; i<len; ++i)
    sdm(row, i) = sdv[i];
}


/// copy ith row of Teuchos::SerialDenseMatrix<OrdinalType, ScalarType>
/// to Teuchos::SerialDenseVector<OrdinalType, ScalarType>
template <typename OrdinalType, typename ScalarType> 
void copy_row(const Teuchos::SerialDenseMatrix<OrdinalType, ScalarType>& sdm,
	      OrdinalType row,
	      const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& sdv)
{
  OrdinalType i, len = sdm.numCols();
  if (sdv.length() != len)
    sdv.sizeUninitialized(len);
  for (OrdinalType i=0; i<len; ++i)
    sdv[i] = sdm(row, i);
}


/// std::ostream write for Teuchos::SerialDenseVector
template <typename OrdinalType, typename ScalarType>
void write_data(std::ostream& s,
		const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& v)
{
  s.setf(std::ios::scientific);
  s << std::setprecision(WRITE_PRECISION);
  OrdinalType len = v.length();
  for (OrdinalType i=0; i<len; i++)
    s << "                     " << std::setw(WRITE_PRECISION+7)
      << v[i] << '\n';
}


/// formatted ostream insertion operator for SerialDenseMatrix
template <typename OrdinalType, typename ScalarType>
void write_data(std::ostream& s,
                const Teuchos::SerialDenseMatrix<OrdinalType, ScalarType>& m,
                bool brackets, bool row_rtn, bool final_rtn)
{
  OrdinalType i, j, nrows = m.numRows(), ncols = m.numCols();
  s.setf(std::ios::scientific); // formatting optimized for T = double
  s << std::setprecision(WRITE_PRECISION);
  if (brackets)  s << "[[ ";
  for (i=0; i<nrows; ++i) {
    for (j=0; j<ncols; ++j)
      s << std::setw(WRITE_PRECISION+7) << m(i,j) << ' ';
    if (row_rtn && i!=m.numRows()-1)
      s << "\n   ";
  }
  if (brackets)  s << "]] ";
  if (final_rtn) s << '\n';
}


/// formatted ostream insertion operator for SerialSymDenseMatrix
template <typename OrdinalType, typename ScalarType>
void write_data(std::ostream& s,
                const Teuchos::SerialSymDenseMatrix<OrdinalType, ScalarType>& m,
                bool brackets, bool row_rtn, bool final_rtn)
{
  OrdinalType i, j, nrows = m.numRows();
  s.setf(std::ios::scientific); // formatting optimized for T = double
  s << std::setprecision(WRITE_PRECISION);
  if (brackets)  s << "[[ ";
  for (i=0; i<nrows; ++i) {
    for (j=0; j<nrows; ++j)
      s << std::setw(WRITE_PRECISION+7) << m(i,j) << ' ';
    if (row_rtn && i!=m.numRows()-1)
      s << "\n   ";
  }
  if (brackets)  s << "]] ";
  if (final_rtn) s << '\n';
}


/// global std::ostream insertion operator for std::vector
template <class T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& data)
{
  s.setf(std::ios::scientific);
  s << std::setprecision(WRITE_PRECISION);
  size_t i=0, len = data.size();
  for (i=0; i<len; ++i)
    s << "                     " << std::setw(WRITE_PRECISION+7)
      << data[i] << '\n';
  return s;
}


/// global std::ostream insertion operator for std::set
template <class T>
std::ostream& operator<<(std::ostream& s, const std::set<T>& data)
{
  for (typename std::set<T>::const_iterator cit = data.begin();
       cit != data.end(); ++cit)
    s << "                     " << std::setw(WRITE_PRECISION+7)
      << *cit << '\n';
  return s;
}

} // namespace Pecos

#endif // PECOS_DATA_TYPES_H

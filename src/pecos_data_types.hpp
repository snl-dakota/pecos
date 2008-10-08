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
//typedef Teuchos::SerialDenseVector<int, std::complex<Real> > ComplexVector;
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


// --------------------
// Conversion templates
// --------------------
/// copy T* to std::vector<std::vector<T> >
template <class T>
void copy_data(const T* ptr, const int ptr_len, const String& ptr_type,
	       std::vector<std::vector<T> >& p2a, size_t num_vec,
	       size_t vec_len)
{
  if (num_vec && vec_len) { // both specified
    if (ptr_len != num_vec*vec_len) {
      PCerr << "Error: pointer allocation (" << ptr_len << ") does not equal "
	    << "num_vec*vec_len (" << num_vec << '*' << vec_len << ") in "
	    << "copy_data(T*, Pecos::std::vector<std::vector<T> >)."
	    << std::endl;
      abort_handler(-1);
    }
  }
  else if (num_vec) { // only num_vec is non-zero
    if (ptr_len%num_vec) {
      PCerr << "Error: pointer allocation (" << ptr_len << ") not evenly "
	    << "divisible by number of vectors (" << num_vec << ") in "
	    << "copy_data(T*, Pecos::std::vector<std::vector<T> >)."
	    << std::endl;
      abort_handler(-1);
    }
    vec_len = ptr_len/num_vec;
  }
  else if (vec_len) { // only vec_len is non-zero
    if (ptr_len%vec_len) {
      PCerr << "Error: pointer allocation (" << ptr_len << ") not evenly "
	    << "divisible by vector length (" << vec_len << ") in copy_data("
	    << "T*, Pecos::std::vector<std::vector<T> >)." << std::endl;
      abort_handler(-1);
    }
    num_vec = ptr_len/vec_len;
  }
  else { // neither specified
    PCerr << "Error: either num_vec or vec_len must be specified in "
	  << "copy_data(T*, Pecos::std::vector<std::vector<T> >)." << std::endl;
    abort_handler(-1);
  }
  size_t i, j;
  if (p2a.size() != num_vec)
    p2a.resize(num_vec);
  for (i=0; i<num_vec; i++)
    if (p2a[i].size() != vec_len)
      p2a[i].resize(vec_len);
  int cntr = 0;
  if (ptr_type == "c" || ptr_type == "C") {
    for (i=0; i<num_vec; i++)   // loop over rows
      for (j=0; j<vec_len; j++) // loop over columns
        p2a[i][j] = ptr[cntr++];
  }
  else {
    for (j=0; j<vec_len; j++)   // loop over columns
      for (i=0; i<num_vec; i++) // loop over rows
        p2a[i][j] = ptr[cntr++];
  }
}


/// copy std::vector<std::vector<T> > to T*
template <class T>
void copy_data(const std::vector<std::vector<T> >& p2a, T* ptr,
	       const int ptr_len, const String& ptr_type)
{
  bool c_type = (ptr_type == "c" || ptr_type == "C");
  size_t i, j, num_vec = p2a.size(), total_len = 0, max_vec_len = 0;
  for (i=0; i<num_vec; i++) { // loop over vectors in array
    size_t vec_len = p2a[i].size();
    total_len += vec_len;
    if (!c_type && vec_len > max_vec_len)
      max_vec_len = vec_len;
  }
  if (ptr_len != total_len) {
    PCerr << "Error: bad ptr_len in copy_data(Pecos::std::vector<std::vector<T>"
	  << " >, T* ptr)." << std::endl;
    abort_handler(-1);
  }
  int cntr = 0;
  if (c_type) {
    for (i=0; i<num_vec; i++) { // loop over rows
      size_t vec_len = p2a[i].size(); // allowed to vary
      for (j=0; j<vec_len; j++) // loop over columns
        ptr[cntr++] = p2a[i][j];
    }
  }
  else {
    for (j=0; j<max_vec_len; j++) // loop over longest column
      for (i=0; i<num_vec; i++) // loop over rows
	if (j < p2a[i].size())
	  ptr[cntr++] = p2a[i][j];
  }
}

} // namespace Pecos

#endif // PECOS_DATA_TYPES_H

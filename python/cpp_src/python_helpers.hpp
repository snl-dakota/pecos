#ifndef PYTHON_HELPERS_HPP
#define PYTHON_HELPERS_HPP

// NumPy include
// Python.h must be included before any headers to avoid 
// warning: "_POSIX_C_SOURCE" redefined. 
// Python.h is included in numpy_include.hpp
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
#include "Teuchos_SerialDenseVector.hpp"

#include "OptionsList.hpp"
#include <boost/shared_ptr.hpp>

namespace Surrogates{

/** \brief Copy a numpy ndarray into a Teuchos SerialDenseVector (SDV).
 * T is the scalar type of the SDV and S is the scalar type of the numpy array
 * Having these two template parameters instead of just one T allows the
 # conversion of both NPY_INT and NPY_LONG into A SDV<int,int>.
 */
template< typename T, typename S >
void copyNumPyToTeuchosVector(PyObject * pyArray,
			      Teuchos::SerialDenseVector<int,T > & tvec);

template< typename T >
PyObject * copyTeuchosVectorToNumPy(Teuchos::SerialDenseVector< int,T > &tvec);

bool setPythonParameter(OptionsList & opts_list,
			const std::string      & name,
			PyObject               * value);

bool updateOptionsListWithPyDict(PyObject    * dict,
				   OptionsList & opts_list);

OptionsList * pyDictToNewOptionsList(PyObject* dict);

PyObject * getPythonParameter(const OptionsList & plist,
			      const std::string & name);


bool updatePyDictWithOptionsList(PyObject          * dict,
				 const OptionsList & opts_list);

PyObject * optionsListToNewPyDict(const OptionsList & opts_list);

template< typename TYPE >
int NumPy_TypeCode();


} //namespace Surrogates
#endif // PYTHON_HELPERS_HPP

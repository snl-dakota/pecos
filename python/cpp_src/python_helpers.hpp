#ifndef PYTHON_HELPERS_HPP
#define PYTHON_HELPERS_HPP

// NumPy include
#define NO_IMPORT_ARRAY
// Python.h must be included before any headers to avoid 
// warning: "_POSIX_C_SOURCE" redefined. 
// Python.h is included in numpy_include.hpp
#include "numpy_include.hpp"


#include "OptionsList.hpp"
#include <boost/shared_ptr.hpp>

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
#endif // PYTHON_HELPERS_HPP

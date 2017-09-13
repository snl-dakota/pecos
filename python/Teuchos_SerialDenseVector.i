// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

// Teuchos_SerialDenseVector.i is a SWIG interface file that provides SWIG
// directives to handle Teuchos SerialDenseVector types.  This class is not
// wrapped, but instead typemaps are defined so that the python user
// can use NumPy arrays or Python sequences instead.
%{
#include "Teuchos_SerialDenseVector.hpp"
using Teuchos::SerialDenseVector;
%}
%import  "Teuchos_SerialDenseVector.hpp"

////////////////////////////////////////////////////////////////////////
// The philosophy is that wherever Teuchos SerialDenseVector classes are used in
// C++, NumPy arrays will be used in python.  Thus we need the NumPy
// SWIG directives.
%include "numpy.i"

////////////////////////////////////////////////////////////////////////
// Define a macro that takes a C++ data ordinal and scalar type
// (ORDINAL_TYPE,SCALAR_TYPE) and a
// corresponding NumPy typecode (TYPECODE) and define all of the
// typemaps needed to handle that TYPE array.
%define %teuchos_sdv_typemaps(ORDINAL_TYPE, SCALAR_TYPE, TYPECODE)

// If an SerialDenseVector argument has a template parameter argument that is
// a const TYPE, then we know that the argument is input only.
// Therefore we allow any type of sequence to be converted to a
// PyArrayObject and then extract the resulting data pointer to
// construct the SerialDenseVector.  If the conversion creates a new
// PyArrayObject, then we have to be sure to decrement its reference
// count once the SerialDenseVector has been used.

//////////////////////////////////////
// Teuchos::SerialDenseVector<const ORDINAL_TYPE,const SCALAR_TYPE> //
//////////////////////////////////////
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
  (Teuchos::SerialDenseVector<const ORDINAL_TYPE,const SCALAR_TYPE>)
{
  $1 = is_array($input) || PySequence_Check($input);
}

%typemap(in) Teuchos::SerialDenseVector<const ORDINAL_TYPE,const SCALAR_TYPE>
(int is_new = 0,
 PyArrayObject * npArray = NULL)
{
  npArray = obj_to_array_fortran_allow_conversion($input, TYPECODE, &is_new);
  if (!npArray) SWIG_fail;
  $1 = Teuchos::SerialDenseVector<const ORDINAL_TYPE,const SCALAR_TYPE>(
	       Teuchos::Copy,
	       (SCALAR_TYPE*) array_data(npArray),
	       array_size(npArray, 0));
}

%typemap(freearg) Teuchos::SerialDenseVector<const ORDINAL_TYPE,const SCALAR_TYPE>
{
  if (is_new$argnum) Py_DECREF(npArray$argnum);
}

%typemap(out) Teuchos::SerialDenseVector<const ORDINAL_TYPE,const SCALAR_TYPE>
{
  npy_intp dims[1] = { $1.length() };
  $result = PyArray_SimpleNewFromData(1, dims, TYPECODE, (void*) $1.values());
  if (!$result) SWIG_fail;
}

//////////////////////////////////////////////
// Teuchos::SerialDenseVector<const ORDINAL_TYPE,const SCALAR_TYPE> const & //
//////////////////////////////////////////////
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
  (Teuchos::SerialDenseVector<const ORDINAL_TYPE,const SCALAR_TYPE> const &)
{
  $1 = is_array($input) || PySequence_Check($input);
}

%typemap(in) Teuchos::SerialDenseVector<const ORDINAL_TYPE,const SCALAR_TYPE> const &
(int is_new = 0,
 PyArrayObject * npArray = NULL,
 Teuchos::SerialDenseVector<const ORDINAL_TYPE,const SCALAR_TYPE> temp)
{
  npArray = obj_to_array_fortran_allow_conversion($input, TYPECODE, &is_new);
  if (!npArray) SWIG_fail;
  temp = Teuchos::SerialDenseVector<const ORDINAL_TYPE,const SCALAR_TYPE>(
		 Teuchos::Copy,
		 (SCALAR_TYPE*) array_data(npArray),
		 array_size(npArray, 0));
  $1 = &temp;
}

%typemap(freearg) Teuchos::SerialDenseVector<const ORDINAL_TYPE,const SCALAR_TYPE> const &
{
  if (is_new$argnum) Py_DECREF(npArray$argnum);
}

// If an SerialDenseVector argument has template parameter argument that is a
// non-const TYPE, then the default behavior is to assume that the
// array is input/output.  Therefore the input python argument must be
// a NumPy array.

////////////////////////////////
// Teuchos::SerialDenseVector<ORDINAL_TYPE, SCALAR_TYPE> //
////////////////////////////////
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
  (Teuchos::SerialDenseVector<ORDINAL_TYPE, SCALAR_TYPE>)
{
  $1 = is_array($input) || PySequence_Check($input);
}

%typemap(in) Teuchos::SerialDenseVector<ORDINAL_TYPE, SCALAR_TYPE>
{
  int is_new = 0;
  PyArrayObject * npArray = obj_to_array_fortran_allow_conversion(
      $input, TYPECODE, &is_new);
  if (!npArray) SWIG_fail;
  $1 = Teuchos::SerialDenseVector<ORDINAL_TYPE, SCALAR_TYPE>(
		Teuchos::Copy,
		(SCALAR_TYPE*) array_data(npArray),
		array_size(npArray, 0));
}

%typemap(out) Teuchos::SerialDenseVector<ORDINAL_TYPE, SCALAR_TYPE>
{
  npy_intp dims[1] = { $1.length() };
  $result = PyArray_SimpleNewFromData(1, dims, TYPECODE, (void*) $1.values());
  if (!$result) SWIG_fail;
}

////////////////////////////////////////
// Teuchos::SerialDenseVector<ORDINAL_TYPE, SCALAR_TYPE> const & //
////////////////////////////////////////
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
  (Teuchos::SerialDenseVector<ORDINAL_TYPE, SCALAR_TYPE> const &)
{
  $1 = is_array($input) || PySequence_Check($input);
}

%typemap(in) Teuchos::SerialDenseVector<ORDINAL_TYPE, SCALAR_TYPE> const &
(Teuchos::SerialDenseVector<ORDINAL_TYPE, SCALAR_TYPE> temp)
{
  int is_new = 0;
  PyArrayObject * npArray = obj_to_array_fortran_allow_conversion(
      $input, TYPECODE, &is_new);
  if (!npArray) SWIG_fail;
  temp = Teuchos::SerialDenseVector<ORDINAL_TYPE, SCALAR_TYPE>(
		Teuchos::Copy,
		(SCALAR_TYPE*) array_data(npArray),
		array_size(npArray, 0));
  $1 = &temp;
}

%typemap(out) Teuchos::SerialDenseVector<ORDINAL_TYPE, SCALAR_TYPE> const &
{
  npy_intp dims[1] = { $1->length() };
  $result = PyArray_SimpleNewFromData(1, dims, TYPECODE, (void*) $1->values());
  if (!$result) SWIG_fail;
}

// If an SerialDenseVector is an output argument

////////////////////////////////
// Teuchos::SerialDenseVector<ORDINAL_TYPE, SCALAR_TYPE> & result//
////////////////////////////////

// Specify how to return result to python
%typemap(argout) Teuchos::SerialDenseVector<ORDINAL_TYPE,SCALAR_TYPE> &argout
{
  npy_intp dims[1] = { $1->length() };
  //PyObject* npArray = PyArray_SimpleNewFromData(1, dims, TYPECODE, (void*) $1->values());
  PyObject* npArray = PyArray_SimpleNew( 1, dims, TYPECODE );
  if (!$result) SWIG_fail;
  SCALAR_TYPE *array = (SCALAR_TYPE *)PyArray_DATA( npArray );
  for ( int i = 0; i < $1->length(); i++ )
    *array++ = $1->values()[i];
  $result = SWIG_Python_AppendOutput($result,npArray);
}

// Remove result from the python function call.
%typemap(in,numinputs=0) Teuchos::SerialDenseVector<ORDINAL_TYPE,SCALAR_TYPE> &argout
  ( Teuchos::SerialDenseVector<ORDINAL_TYPE,SCALAR_TYPE> TSDV )
{
  // Must declare the variable that needs to be passed to the c++ function
  TSDV = Teuchos::SerialDenseVector<ORDINAL_TYPE,SCALAR_TYPE>();
  $1 = &TSDV;
};

%typemap(freearg) Teuchos::SerialDenseVector<ORDINAL_TYPE,SCALAR_TYPE> &argout
{
  // Do nothing. Needed so that default free arg is not used
  // Check with Bill Spotz
};

%typemap(directorin) Teuchos::SerialDenseVector<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> &%{
  npy_intp dims$argnum[1] = { $1_name.length() };
  CPP_SCALAR_TYPE * data$argnum = $1_name.values();
  PyArray_Descr * dtype$argnum = PyArray_DescrFromType( NPYTYPE );
  PyArrayObject* py_array$argnum =
    (PyArrayObject*) PyArray_NewFromDescr( &PyArray_Type, dtype$argnum, 1,
					   dims$argnum, NULL, (void*)data$argnum,
					   NPY_FARRAY, NULL );

  // Cast to (PyObject*) is necessary
  $input = (PyObject*)py_array$argnum;
  %}

%typemap(directorargout) Teuchos::SerialDenseVector<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> &argout%{
  int is_new = 0;
  PyArrayObject * npArray$argnum = obj_to_array_fortran_allow_conversion($input, TYPECODE, &is_new);
  if (!npArray$argnum)
    throw( std::runtime_error( "Conversion to array failed" ) );
  int m$argnum = array_size(npArray$argnum, 0);
  Teuchos::SerialDenseVector<ORDINAL_TYPE, SCALAR_TYPE> temp$argnum(
     Teuchos::View, (SCALAR_TYPE*) array_data(npArray$argnum), m$argnum);
  $1.sizeUninitialized(m$argnum);
  $1.assign(temp$argnum);
%}
%enddef

////////////////////////////////////////////////////////////////////////
// Call the %teuchos_sdv_typemaps() macro for specific data types
// that are supported by NumPy
%teuchos_sdv_typemaps(int, signed char       , NPY_BYTE     )
%teuchos_sdv_typemaps(int, unsigned char     , NPY_UBYTE    )
%teuchos_sdv_typemaps(int, short             , NPY_SHORT    )
%teuchos_sdv_typemaps(int, unsigned short    , NPY_USHORT   )
%teuchos_sdv_typemaps(int, int               , NPY_INT      )
%teuchos_sdv_typemaps(int, unsigned int      , NPY_UINT     )
%teuchos_sdv_typemaps(int, long              , NPY_LONG     )
%teuchos_sdv_typemaps(int, unsigned long     , NPY_ULONG    )
%teuchos_sdv_typemaps(int, long long         , NPY_LONGLONG )
%teuchos_sdv_typemaps(int, unsigned long long, NPY_ULONGLONG)
%teuchos_sdv_typemaps(int, float             , NPY_FLOAT    )
%teuchos_sdv_typemaps(int, double            , NPY_DOUBLE   )

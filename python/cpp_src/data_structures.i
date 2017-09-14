%{
#include "numpy_include.hpp"//gets rid of deprecated warnings
#include "Teuchos_SerialDenseVector.hpp"
using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
%}
%import  "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

////////////////////////////////////////////////////////////////////////
// The philosophy is that wherever std::vector of Teuchos SerialDenseVector
// or Teuchos SerialDenseVector classes are used in C++,
// a list of NumPy arrays will be used in python.
// Thus we need the NumPy SWIG directives.
%include "numpy.i"

////////////////////////////////////////////////////////////////////////
// Define a macro that takes a C++ data ordinal and scalar type
// (ORDINAL_TYPE,SCALAR_TYPE) and a
// corresponding NumPy typecode (NPYTYPE) and define all of the
// typemaps needed to handle that TYPE array.

%define %teuchos_serial_dense_vector_list(CPP_ORDINAL_TYPE, CPP_SCALAR_TYPE,
					  NPYTYPE);

%typemap(in,numinputs=0) std::vector< Teuchos::SerialDenseVector<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> > &argout
   ( std::vector< Teuchos::SerialDenseVector<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> >  TSDVL )
{
  TSDVL = std::vector< Teuchos::SerialDenseVector<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> >();
  $1 = &TSDVL;
};

// Specify how to return result to python
%typemap(argout) std::vector< Teuchos::SerialDenseVector<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> > &argout
{
  Py_ssize_t len = $1->size();
  PyObject* list = PyList_New( len );
  for ( int j = 0; j < len; j++ )
    {
      // PyArray_SimpleNew will set refcnt to 1
      npy_intp dims[1] = { (*$1)[j].length() };
      PyObject* array = PyArray_SimpleNew( 1, dims, NPYTYPE );
      if (!array) SWIG_fail;
      CPP_SCALAR_TYPE *array_j = (CPP_SCALAR_TYPE *)PyArray_DATA( array );
      for ( int k = 0; k < dims[0]; k++ )
	*array_j++ = (*$1)[j].values()[k];
      PyList_SET_ITEM( list, j, array );
    }
  $result = SWIG_Python_AppendOutput($result,list);
}

%typemap(freearg) std::vector< Teuchos::SerialDenseVector<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> > &argout
{
  // do nothing
}

// Here is the default input typemap directive.
%typemap(in) std::vector< Teuchos::SerialDenseVector<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> > &
   ( std::vector< Teuchos::SerialDenseVector<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> > TSDVL, PyObject* seq = NULL, std::vector<bool> new_array_built )
{
  seq = PySequence_Fast( $input, "expected a sequence" );
  Py_ssize_t len = PySequence_Size( $input );
  TSDVL = std::vector< Teuchos::SerialDenseVector<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> >( (int)len );

  new_array_built.resize( (int)len );

  for ( int i = 0; i < (int)len; i++ )
    {
      PyObject* item; item = PySequence_Fast_GET_ITEM( seq, i );
      // If the object is already (a subclass of) the ndarray that satisfies
      // the requirements then a new reference is returned. Otherwise,
      // a new array is constructed. When the requirements are satisfied
      // and the new reference is returned we must call DECREF in the
      // corresponding freearg typemap on item.
      PyObject * py_array = PyArray_FROM_OTF( item , NPYTYPE, NPY_IN_ARRAY );
      if ( PyArray_NDIM( py_array ) != 1 )
	throw( std::runtime_error( "list element must be 1 dimensional" ) );
      if ( py_array != item ) new_array_built[i] = true;
      else new_array_built[i] = false;
      if ( py_array == NULL ) SWIG_fail;
      // Now we need to check that the NumPy array that we have is 1D.  If
      // it is higher dimension, we will raise a Python exception.
      if (PyArray_NDIM(py_array) > 1)
	{
	  PyErr_SetString(PyExc_ValueError,"Array data must be one dimensional");
	  SWIG_fail;
	}
      TSDVL[i] = Teuchos::SerialDenseVector<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE>(Teuchos::View, (CPP_SCALAR_TYPE*)PyArray_DATA(py_array), (CPP_ORDINAL_TYPE)PyArray_DIM(py_array,0));
    }
  $1 = &TSDVL;
};

%typemap(freearg) std::vector< Teuchos::SerialDenseVector<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> > &
{
  for ( int i = 0; i < (int)PySequence_Size( seq$argnum ); i++ )
    {
      if ( !new_array_built$argnum[i] )
	{
	  PyObject* item = PySequence_Fast_GET_ITEM( seq$argnum, i );
	  Py_XDECREF( item );
	}
    }

  // PySequence_Fast calls Py_INCREF so we must call Py_XDECREF here
  Py_XDECREF( seq$argnum );
};

%enddef

%define %teuchos_serial_dense_matrix_list(CPP_ORDINAL_TYPE, CPP_SCALAR_TYPE,
					  NPYTYPE);

%typemap(in,numinputs=0) std::vector< Teuchos::SerialDenseMatrix<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> > &argout
   ( std::vector< Teuchos::SerialDenseMatrix<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> >  TSDML )
{
  TSDML = std::vector< Teuchos::SerialDenseMatrix<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> >();
  $1 = &TSDML;
};

// Specify how to return result to python
%typemap(argout) std::vector< Teuchos::SerialDenseMatrix<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> > &argout
{
  Py_ssize_t len = $1->size();
  PyObject* list = PyList_New( len );
  for ( int k = 0; k < len; k++ )
    {
      // PyArray_SimpleNew will set refcnt to 1
      CPP_ORDINAL_TYPE m = (*$1)[k].numRows(), n = (*$1)[k].numCols();
      npy_intp dims[2] = { m, n };
      // Create a fortran ordered numpy array
      PyObject* array = PyArray_EMPTY( 2, dims, NPYTYPE, 1 );
      if (!array) SWIG_fail;
      CPP_SCALAR_TYPE *matrix = (CPP_SCALAR_TYPE *)PyArray_DATA( array );
      for ( int j = 0; j < n; j++ ){
	for ( int i = 0; i < m; i++ ){
	  matrix[j*m+i] = (*$1)[k](i,j);
	}
      }
      PyList_SET_ITEM( list, k, array );
    }
  $result = SWIG_Python_AppendOutput($result,list);
}

%typemap(freearg) std::vector< Teuchos::SerialDenseMatrix<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> > &argout
{
  // do nothing
}

// Here is the default input typemap directive.
%typemap(in) std::vector< Teuchos::SerialDenseMatrix<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> > &
   ( std::vector< Teuchos::SerialDenseMatrix<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> > TSDML, PyObject* seq = NULL, std::vector<bool> new_array_built )
{
  seq = PySequence_Fast( $input, "expected a sequence" );
  Py_ssize_t len = PySequence_Size( $input );
  TSDML = std::vector< Teuchos::SerialDenseMatrix<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> >( (int)len );

  new_array_built.resize( (int)len );

  for ( int i = 0; i < (int)len; i++ )
    {
      PyObject* item; item = PySequence_Fast_GET_ITEM( seq, i );
      // If the object is already (a subclass of) the ndarray that satisfies
      // the requirements then a new reference is returned. Otherwise,
      // a new array is constructed. When the requirements are satisfied
      // and the new reference is returned we must call DECREF in the
      // corresponding freearg typemap on item.
      PyObject *py_array = PyArray_FROM_OTF( item , NPYTYPE, NPY_IN_ARRAY );
      if ( PyArray_NDIM( py_array ) != 2 )
	throw( std::runtime_error( "list element must be 2 dimensional" ) );
      if ( py_array != item ) new_array_built[i] = true;
      else new_array_built[i] = false;
      if ( py_array == NULL ) SWIG_fail;
      // Now we need to check that the NumPy array that we have is 1D.  If
      // it is higher dimension, we will raise a Python exception.
      if (PyArray_NDIM(py_array) > 1)
	{
	  PyErr_SetString(PyExc_ValueError,"Array data must be one dimensional");
	  SWIG_fail;
	}
      CPP_ORDINAL_TYPE stride = (CPP_ORDINAL_TYPE)( PyArray_STRIDE(py_array,1) /
						    PyArray_ITEMSIZE(py_array) ),
	m = (CPP_ORDINAL_TYPE)PyArray_DIM(py_array,0),
	n = (CPP_ORDINAL_TYPE)PyArray_DIM(py_array,1);
      TSDML[i] = Teuchos::SerialDenseMatrix<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE>(Teuchos::View, (CPP_SCALAR_TYPE*)PyArray_DATA(py_array), stride, m, n );
    }
  $1 = &TSDML;
};

%typemap(freearg) std::vector< Teuchos::SerialDenseMatrix<CPP_ORDINAL_TYPE,CPP_SCALAR_TYPE> > &
{
  for ( int i = 0; i < (int)PySequence_Size( seq$argnum ); i++ )
    {
      if ( !new_array_built$argnum[i] )
	{
	  PyObject* item = PySequence_Fast_GET_ITEM( seq$argnum, i );
	  Py_XDECREF( item );
	}
    }

  // PySequence_Fast calls Py_INCREF so we must call Py_XDECREF here
  Py_XDECREF( seq$argnum );
};

%enddef

// Let's instantiate some concrete versions of the available typemaps
// CPP_ORDINAL_TYPE to match Python's native interger type (long)
%teuchos_serial_dense_vector_list( int, double, NPY_DOUBLE );
%teuchos_serial_dense_vector_list( int, int, NPY_INT );
%teuchos_serial_dense_vector_list( int, Complex, NPY_CDOUBLE );

// typemaps for std::vector<Teuchos_SerialDenseMatrix> setting
%teuchos_serial_dense_matrix_list( int, double, NPY_DOUBLE );
%teuchos_serial_dense_matrix_list( int, int, NPY_INT );
%teuchos_serial_dense_matrix_list( int, Complex, NPY_CDOUBLE );

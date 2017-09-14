/* univariate_polynomials.i */
%module(directors=1,package="PyDakota",autodoc=1) univariate_polynomials

/* %include <std_vector.i> */
/* %template() std::vector<double>; */
/* %template() std::vector<int>; */
%{
#include <Python.h> 
#include <numpy/arrayobject.h>
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"//gets rid of deprecated warnings

  //#include "numpy_include.hpp"
#include <vector>
#include "BasisPolynomial.hpp"
typedef double Real;
typedef std::vector<double>                RealArray;
using namespace Pecos;
%}

//%include "fundamentals.i"

//%include "pecos_data_types.hpp"
//%ignore type1_value(unsigned_short);

%include "numpy.i"
%fragment("NumPy_Fragments");
%init %{
  import_array();
%}

// Standard exception handling
%include "exception.i"
%exception
{
  try{                                          
    // use these print statements to debug
    //printf("Entering function : $name\n");
   $action
   if (PyErr_Occurred()) SWIG_fail;
   //else    printf("Exiting function : $name\n");
  }
  SWIG_CATCH_STDEXCEPT
    catch (...) {
    SWIG_exception(SWIG_UnknownError, "unknown C++ exception");
    throw(std::runtime_error(""));
  }
}
// If get
// SystemError: return without exception seterror
// this can mean swig_fail is being called but not working for some reason

%define %std_vector_typemaps(SCALAR_TYPE, TYPECODE)
%typemap(out) std::vector<SCALAR_TYPE> const &
{
  npy_intp dims[1] = { $1->size() };
  $result = PyArray_SimpleNew( 1, dims, TYPECODE );
  if (!$result) SWIG_fail;
  SCALAR_TYPE *array = (SCALAR_TYPE *)PyArray_DATA( $result );
  for (std::vector<SCALAR_TYPE>::iterator it=$1->begin() ; it!=$1->end(); ++it){
    *array++ = *it;
  }
}
%enddef
%std_vector_typemaps( double            , NPY_DOUBLE   )

%include "pecos_global_defs.hpp"
%include "BasisPolynomial.hpp"

typedef double Real;
typedef std::vector<double>                RealArray;

%pythoncode %{
    import numpy as np
%}

%extend Pecos::BasisPolynomial{  
%pythoncode
%{
    def values(self,samples,degree):
        """
        Evaluate the polynomial basis at a set of samples for all degrees
        d=0,....,degree.

        TODO: currently this wraps self.type1_value. In future create
        a new C++ function
        RealMatrix& values(const RealVector& samples, degree)
        and use swig to wrap this. rename values in .i file and use
        this function to deal with cases when samples is a scalar and a
        np.ndarray

        Parameters
        ----------
        samples : np.ndarray (nsamples) or double
            The samples at which to evaluate the polynomial.

        degree : integer
            The maximum degree.

        Returns
        -------
        values : np.ndarray (nsamples x nterms)
            The basis value, for each degree, at each sample.
        """
        if np.isscalar(samples):
            samples = np.asarray([samples])
        assert samples.ndim==1
        nsamples = samples.shape[0]; nterms=degree+1
        values = np.empty((nsamples, nterms),dtype=float)
        for i in range(nsamples):
            for j in range(nterms):
                values[i,j] = self.type1_value(samples[i], j)
        return values
%}
}

/* math_tools.i */
%module(directors=1,implicitconv="1", package="PyDakota",autodoc=1) math_tools
%{
  //#include <Python.h>
#include <numpy/arrayobject.h>
#include "numpy_include.hpp"//gets rid of deprecated warnings
  
// Local includes
#include "math_tools.hpp"
#include "linear_algebra.hpp"
%}

// %ignore and rename must be included before the function decleration, i.e.
// before the %include
%ignore *::operator[];
%ignore *::operator=;
%ignore *::print;

%ignore Surrogates::qr_solve( const RealMatrix &, const RealMatrix &, RealMatrix & );
%rename(tensor_product_indices_cpp) tensor_product_indices;

//%include "teuchos_data_types.hpp"

// include Teuchos enums, such as TRANS, NO_TRANS
%include "Teuchos_BLAS_types.hpp"
%include "fundamentals.i"
%include "math_tools.hpp"
%include "linear_algebra.hpp"

namespace Surrogates{
%template(cartesian_product_int) cartesian_product<int,int>;
%template(cartesian_product_double) cartesian_product<int,double>;
%template(outer_product_int) outer_product<int,int>;
%template(outer_product_double) outer_product<int,double>;
}


%pythoncode %{
import numpy
def remove_common_columns( A, B ):
    """
    Remove columns from A that are also in B
    """
    D = numpy.hstack( ( B, A ) )
    order = numpy.lexsort( D )
    C = D.copy()[:,order]
    I = numpy.hstack((numpy.arange(B.shape[1])+1,-numpy.arange(A.shape[1]) ))
    I = I[order]
    diff = numpy.diff( C, axis = 1 )
    ui = numpy.ones( C.shape[1], 'bool' )

    ui[1:] = ( numpy.absolute( diff ) >
               numpy.finfo( numpy.double ).eps ).any( axis = 0 )
    ui = ( I[ui] <= 0 )
    return A[:,ui]

def unique_matrix_rows(A):
    """
    Remove duplicate columns from A
    """
    return numpy.vstack(set(tuple(row) for row in A))

def unique_matrix_cols(A):
    return unique_matrix_cols(A.T).T

def cartesian_product(input_sets,elem_size=1):
    """Wrapper of cpp function that converts to accepted types."""
    if input_sets[0].dtype==float or input_sets[0].dtype==numpy.float64:
        return cartesian_product_double(input_sets,elem_size)
    elif (input_sets[0].dtype==numpy.int64 or input_sets[0].dtype==numpy.int32 or input_sets[0].dtype==int):
        for i in xrange(len(input_sets)):
            if input_sets[i].dtype!=numpy.int32:
                input_sets[i] = numpy.asarray(input_sets[i],dtype=numpy.int32)
        return cartesian_product_int(input_sets,elem_size)
    else:
        raise Exception, 'element type not supported'

def outer_product(input_sets,elem_size=1):
    """Wrapper of cpp function that converts to accepted types."""
    if input_sets[0].dtype==float or input_sets[0].dtype==numpy.float64:
        return outer_product_double(input_sets)
    elif (input_sets[0].dtype==numpy.int64 or input_sets[0].dtype==numpy.int32 or input_sets[0].dtype==int):
        for i in xrange(len(input_sets)):
            if input_sets[i].dtype!=numpy.int32:
                input_sets[i] = numpy.asarray(input_sets[i],dtype=numpy.int32)
        return outer_product_int(input_sets)
    else:
        raise Exception, 'element type not supported'

def tensor_product_indices(degrees):
    """Wrapper of cpp function that converts to accepted integer type."""
    assert degrees.dtype==numpy.int32 or degrees.dtype==numpy.int64
    if degrees.dtype!=numpy.int32:
        degrees = numpy.asarray(degrees,dtype=numpy.int32)
    return tensor_product_indices_cpp(degrees)
%}

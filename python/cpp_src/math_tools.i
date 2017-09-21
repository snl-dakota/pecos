/* math_tools.i */
%module(directors=1,implicitconv="1", package="PyDakota",autodoc=1) math_tools
%{
  #define SWIG_FILE_WITH_INIT
  
  // #define NO_IMPORT_ARRAY causes missing symbol error when loading
  // _math_tools module when used with %fragment("NumPy_Fragments");
  #include "numpy_include.hpp"
  
// Local includes
#include "math_tools.hpp"
#include "linear_algebra.hpp"
#include "Teuchos_SerialDenseVector.hpp"
%}

// Global swig features
%feature("autodoc", "1");
%feature("compactdefaultargs");

// %ignore and rename must be included before the function decleration, i.e.
// before the %include
%ignore *::operator[];
%ignore *::operator=;
%ignore *::print;

// We utilize very little of numpy.i, so some of the fragments
// do not get automatically instantiated.  This forces the issue.
%include "numpy.i"
%fragment("NumPy_Fragments"); //needed to include obj_to_array_fortran_allow_conversion

// include Teuchos enums, such as TRANS, NO_TRANS
%include "Teuchos_BLAS_types.hpp"
%include "Teuchos_SerialDenseVector.i"
%include "Teuchos_SerialDenseMatrix.i"
%include "data_structures.i"

%ignore Surrogates::qr_solve( const RealMatrix &, const RealMatrix &, RealMatrix & );
%rename(tensor_product_indices_cpp) tensor_product_indices;

// Must specify here to ensure that functions involving the function with
// parameters renamed using typedef can be wrapped.
// If no match is found using the above rules SWIG applies a typedef
// reduction to the type and repeats the typemap search for the reduced type
typedef double Real;
typedef Teuchos::SerialDenseVector<int,double> RealVector;
typedef Teuchos::SerialDenseVector<int,int> IntVector;
typedef Teuchos::SerialDenseVector<int,Complex> ComplexVector;
typedef Teuchos::SerialDenseMatrix<int,double> RealMatrix;
typedef Teuchos::SerialDenseMatrix<int,int> IntMatrix;
typedef Teuchos::SerialDenseMatrix<int,Complex> ComplexMatrix;

%apply IntVector &argout { IntVector &result }
%apply IntVector &argout { IntVector &result_0 }
%apply IntVector &argout { IntVector &result_1 }
%apply IntMatrix &argout { IntMatrix &result }
%apply IntMatrix &argout { IntMatrix &result_0 }
%apply IntMatrix &argout { IntMatrix &result_1 }

%apply RealVector &argout { RealVector &result }
%apply RealVector &argout { RealVector &result_0 }
%apply RealVector &argout { RealVector &result_1 }
%apply RealMatrix &argout { RealMatrix &result }
%apply RealMatrix &argout { RealMatrix &result_0 }
%apply RealMatrix &argout { RealMatrix &result_1 }

%apply std::vector<RealVector> &argout {std::vector<RealVector> &result}
%apply std::vector<RealVector> &argout {std::vector<RealVector> &result_0}
%apply std::vector<RealVector> &argout {std::vector<RealVector> &result_1}
%apply std::vector<RealMatrix>  &argout {std::vector<RealMatrix> &result}
%apply std::vector<RealMatrix>  &argout {std::vector<RealMatrix> &result_0}
%apply std::vector<RealMatrix>  &argout {std::vector<RealMatrix> &result_1}

%apply std::vector<IntVector> &argout {std::vector<IntVector> &result}
%apply std::vector<IntVector> &argout {std::vector<IntVector> &result_0}
%apply std::vector<IntVector> &argout {std::vector<IntVector> &result_1}
%apply std::vector<IntMatrix>  &argout {std::vector<IntMatrix> &result}
%apply std::vector<IntMatrix>  &argout {std::vector<IntMatrix> &result_0}
%apply std::vector<IntMatrix>  &argout {std::vector<IntMatrix> &result_1}

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

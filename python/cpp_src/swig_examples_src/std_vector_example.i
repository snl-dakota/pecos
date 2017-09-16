%module(implicitconv="1", autodoc="1",package="PyDakota.swig_examples") std_vector_example

%include <std_vector.i>
%template(DoubleStdVector) std::vector<double>;
%template(IntStdVector) std::vector<int>;
%template() std::vector<short>;
%{
#include <Python.h> 
#include <numpy/arrayobject.h>
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
#include <vector>
#include "std_vector_example.hpp"
%}
//%include "fundamentals.i"

%include "numpy.i"
%fragment("NumPy_Fragments");
%init %{
  import_array();
%}

%ignore *::operator=;

%include "stl.i"
%include "std_vector_example_type_defs.hpp"
%include "std_vector_example.hpp"

 // I cannot make following work. only can make work import typedefs file
 //typedef double Real; */
 //typedef std::vector<short> Other::ShortArray; */
 //typedef std::vector<Real> RealArray; */

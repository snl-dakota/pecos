%module(implicitconv="1", autodoc="1",package="PyDakota.swig_examples") options_list_interface
%{
#include "numpy_include.hpp"//gets rid of deprecated warnings
#include "options_list_interface.hpp"
#include "python_helpers.hpp"
%}
%include "numpy.i"
%fragment("NumPy_Fragments");
%init %{
  import_array();
%}

%import "OptionsList.i"
%include "options_list_interface.hpp"

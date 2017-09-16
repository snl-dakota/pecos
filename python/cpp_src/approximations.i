/* approximations.i */

%define %pydakota_docstring
"
PyDakota.surrogates is the python interface to the Dakota tools and
approximation packages:

    https://dakota.sandia.gov/
"
%enddef

%module(package      = "PyDakota",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %pydakota_docstring) approximation

%feature("director") Function;

// The following code is inserted directly into the wrapper file
%{
#include <Python.h> 
#include <numpy/arrayobject.h>
#include "numpy_include.hpp"//gets rid of deprecated warnings

// Approximation includes
#include "Function.hpp"
#include "CppFunction.hpp"
#include "Approximation.hpp"
#include "PolyApproximation.hpp"
#include "Monomial.hpp"
#include "PolynomialChaosExpansion.hpp"
#include "Variables.hpp"
#include "BoundedVariables.hpp"
#include "VariableTransformation.hpp"
#include "AffineVariableTransformation.hpp"
#include "polynomial_approximation_drivers.hpp"
#include "RegressionBuilder.hpp"
  using std::string;
  using namespace Pecos;
#include "typedefs_for_python_wrapper.hpp"
%}

%ignore *::operator[];

// How do I make math_tools a separate module of submodule?
%include "fundamentals.i"
%include "stl.i"
%template() std::vector<short>;
%include "typedefs_for_python_wrapper.hpp"

%shared_ptr(Surrogates::Function)
%shared_ptr(Surrogates::CppFunction)
%shared_ptr(Surrogates::Approximation)
%shared_ptr(Surrogates::PolyApproximation)
%shared_ptr(Surrogates::Monomial)
%shared_ptr(Surrogates::PolynomialChaosExpansion)
%shared_ptr(Surrogates::Variables)
%shared_ptr(Surrogates::BoundedVariables)
%shared_ptr(Surrogates::VariableTransformation)
%shared_ptr(Surrogates::AffineVariableTransformation)

%include "OptionsList.i"
%include "Function.hpp"
%include "CppFunction.hpp"
%include "Variables.hpp"
%include "BoundedVariables.hpp"
%include "VariableTransformation.hpp"
%include "AffineVariableTransformation.hpp"
// If classes uses other classes then the classes it use must
// be declared first. E.g. Approximation.hpp calls VariableTransformation.hpp
// Also if derived class does not have redefinition of virtual function in
// base class The base class must be %include so that python can see the
 // virtual function
%include "Approximation.hpp"
%include "PolyApproximation.hpp"
%include "Monomial.hpp"
%include "PolynomialChaosExpansion.hpp"
%include "SurrogateBuilder.hpp"
%include "RegressionBuilder.hpp"

%include "polynomial_approximation_drivers.hpp"

%pythoncode %{
import numpy
class PyFunction(Function):
    def __init__(self,target_function):
        """
        Parameters
        ----------
        target_function : callable function
            Calls to target funcation are assumed to follow
            vals = target_function(sample). Where vals
            is a 1D array and sample is a 1D array or scalar.
        """
        Function.__init__(self)
        self.target_function = target_function

    def value(self,samples):
        """
        Evaluate the function at a set of samples.

        The number of QoI of the vectored_valued function is determined
        by probing the target_function with the first sample in the set
        of samples.

        Parameters
        ----------
        samples : (num_vars x num_samples)
            The coordinates of the samples.
            If samples is a 1D array it will be assumed that
            samples is only one sample and be converted to
            a matrix (num_vars x 1)

        Returns
        -------
        values : (num_samples x num_qoi) matrix
            The vector-valued function value at the samples
        """
	if samples.ndim==1:
            samples = samples.reshape(samples.shape[0],1)
        num_samples = samples.shape[1]
        values = numpy.empty((num_samples),float)
        values_0 = self.target_function(samples[:,0])
        if numpy.isscalar(values_0):
            values_0 = numpy.array([values_0])
        assert values_0.ndim==1
        num_qoi = values_0.shape[0]
        values = numpy.empty((num_samples,num_qoi),float)
        values[0,:] = values_0
        for i in xrange(1,samples.shape[1]):
            values_i = self.target_function(samples[:,i])
            if numpy.isscalar(values_i):
                values_i = numpy.array([values_i])
            values[i,:] = values_i
        return values
 %}

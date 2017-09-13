/* regression.i */
%module(directors="1", implicitconv="1", autodoc="1",package="PyDakota") regression
%feature("director") LinearSystemSolver;
%{
// Required for interfacing with NumPy
#include <Python.h>
#include <numpy/arrayobject.h>

// Local includes
#include "linear_solvers.hpp"
#include "OMPSolver.hpp"
#include "LARSolver.hpp"
#include "LSQSolver.hpp"
#include "EqConstrainedLSQSolver.hpp"
#include "CrossValidatedSolver.hpp"
#include "CrossValidationIterator.hpp"
#include "LinearSystemCrossValidationIterator.hpp"
#include "LSQCrossValidationIterator.hpp"
#include "CrossValidatedSolver.hpp"

using namespace Surrogates;
%}

%rename(extract_values_cpp) extract_values;
// note following is applied to that class and all derived classes
%rename(solve_cpp) Surrogates::LinearSystemSolver::solve;
%rename(run_cpp) Surrogates::LinearSystemCrossValidationIterator::run;
%ignore *::operator[];

%include "fundamentals.i"
%import "math_tools.i" // If I create math_tools as a seperate module I start
 // getting errors when trying to import both regression and math_tools
 // and creating a ParameterList from a python dictionary.
 // The error looks like Warning! The following Teuchos::RCPNode objects
 // were created but have...
%include "OptionsList.i"

%shared_ptr(Surrogates::LinearSystemSolver)
%shared_ptr(Surrogates::SparseSolver)
%shared_ptr(Surrogates::OMPSolver)
%shared_ptr(Surrogates::LARSolver)
%shared_ptr(Surrogates::LSQSolver)
%shared_ptr(Surrogates::EqConstrainedLSQSolver)
%shared_ptr(Surrogates::CrossValidatedSolver)

%shared_ptr(Surrogates::CrossValidatedSolver)
%shared_ptr(Surrogates::CrossValidationIterator)
%shared_ptr(Surrogates::LinearSystemCrossValidationIteratorBase)
%shared_ptr(Surrogates::LinearSystemCrossValidationIterator)
%shared_ptr(Surrogates::LSQCrossValidationIterator)

// %%include of a base class needs to be called before %include for derived
// classes. LinearSolver base also needed to get enums. perhaps move enums to a type definitions file
%include "LinearSystemSolver.hpp" 
%include "OMPSolver.hpp"
%include "LARSolver.hpp"
%include "LSQSolver.hpp"
%include "EqConstrainedLSQSolver.hpp"

%include "CrossValidationIterator.hpp"
%include "LinearSystemCrossValidationIterator.hpp"
%include "LSQCrossValidationIterator.hpp"

%include "CrossValidatedSolver.hpp"
%include "linear_solvers.hpp"

%extend Surrogates::CrossValidationIterator {
  %pythoncode
    %{
    def extract_values( self, values, indices ):
        if ( values.ndim==1):
            values  = values.reshape(values.shape[0],1)
        return self.extract_values_cpp( values, indices ).squeeze()
    %}
}

%extend Surrogates::LinearSystemSolver {
  %pythoncode
    %{
     def solve( self, A, rhs, opts ):
        if (rhs.ndim==1):
            rhs = rhs.reshape(rhs.shape[0],1)
        return self.solve_cpp( A, rhs, opts )
    %}
}

%extend Surrogates::LinearSystemCrossValidationIterator {
  %pythoncode
    %{
     def run( self, A, rhs, opts ):
        if (rhs.ndim==1):
            rhs = rhs.reshape(rhs.shape[0],1)
        return self.run_cpp( A, rhs, opts )
    %}
}

%pythoncode
%{import numpy
def cast_linear_cv_iterator(cv_iterator, regression_type):
    if regression_type in [ORTHOG_MATCH_PURSUIT, LEAST_ANGLE_REGRESSION,
                          LASSO_REGRESSION]:
        return cast_to_linear_system_cross_validation_iterator(cv_iterator)
    else:
        return cast_to_least_squares_cross_validation_iterator(cv_iterator)

def cast_linear_system_solver(solver, regression_type):
    if regression_type==ORTHOG_MATCH_PURSUIT:
        return cast_linear_system_solver_to_ompsolver(solver)
    elif regression_type in [LEAST_ANGLE_REGRESSION,LASSO_REGRESSION]:
        return cast_linear_system_solver_to_larssolver(solver)
    elif regression_type == EQ_CONS_LEAST_SQ_REGRESSION:
        return cast_linear_system_solver_to_equalityconstrainedlsqsolver(
            solver)
    elif regression_type in [SVD_LEAST_SQ_REGRESSION,QR_LEAST_SQ_REGRESSION,
                             LU_LEAST_SQ_REGRESSION]:
        return cast_linear_system_solver_to_lsqsolver(solver)
    else:
        raise Exception, 'incorrect regression_type specified'

%}

#ifndef NESTA_HPP
#define NESTA_HPP

#include "linear_algebra.hpp"
#include "math_tools.hpp"

void create_orthogonal_matrix( const RealMatrix &Amatrix, const RealMatrix &rhs,
			       RealMatrix &scaled_matrix, 
			       RealMatrix &scaled_rhs);

void find_initial_solution( const RealMatrix &Amatrix, const RealMatrix &rhs,
			    RealMatrix & initial_solution );

void copy_matrix_col_to_vector( const RealMatrix &matrix, int col_num, 
				Real scalar, RealVector &vector);

void add( Real scalar1, const RealMatrix &matrix1, 
	  Real scalar2, const RealMatrix &matrix2, 
	  RealMatrix &result );

void add( Real scalar1, const RealMatrix &matrix1, 
	  Real scalar2, const RealMatrix &matrix2, 
	  Real scalar3, const RealMatrix &matrix3,
	  RealMatrix &result );

Real evaluate_l1_constraint( const RealMatrix &solution, Real mu, 
			     RealVector &objective_deriv );

void nesta( const RealMatrix &Amatrix, const RealMatrix &rhs, 
	    Real muf, Real delta, Real tol_var, 
	    int max_continuation_iters, int verbosity, 
	    const RealMatrix &warm_start_solution, RealMatrix &solution, 
	    RealMatrix &metrics);

int nesterov_l1_minimization( const RealMatrix &Amatrix, 
			      const RealMatrix &rhs, 
			      const RealMatrix &initial_solution,
			      Real mu, Real delta, 
			      int max_nesterov_iters, Real tol_var,
			      int verbosity,
			      RealMatrix &solution );

#endif // NESTA_HPP

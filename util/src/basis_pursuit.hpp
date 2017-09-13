#ifndef PECOS_BASIS_PURSUIT_H
#define PECOS_BASIS_PURSUIT_H

#include "teuchos_data_types.hpp"


namespace Teuchos {
  class ParameterList;
}


// BMA: These are to be removed.  

namespace Surrogates {

void basis_pursuit_solve( RealMatrix &A, const RealVector &b, RealMatrix &result_0, 
			  RealMatrix &result_1, Teuchos::ParameterList & params );

void bp_dn_solve( RealMatrix &A, const RealVector &b, RealMatrix &result_0, 
		  RealMatrix &result_1, Teuchos::ParameterList & params );


//! @name Interior-point optimization methods auxilary functions.
//@{ 

/**
 * \brief Compute the duality gap of the primal-dual interior point 
 * optimization algorithm. This is called by the function 
 * BasisPursuitPrimalDual().
 */
Real BP_surrogate_duality_gap( RealVector &primal_residual,
			       RealVector &fu1, RealVector &fu2, 
			       RealVector &lamu1, RealVector &lamu2, 
			       RealVector &AtV, Real mu, Real pdtol,
			       Real &tau, Real &slackness_norm );

/**
 * \brief Compute the central point on the central path computed by
 * BPDN_log_barrier_interior_point_method using newton's method.
 */
int BPDN_compute_central_point( RealMatrix &A, const RealVector &b, RealVector &x, 
				RealVector &u, RealMatrix &AtA, 
				Real epsilon, Real &tau, 
				Real newton_tol, int newton_maxiter,
				Real conjugate_gradients_tol = 1e-8,
				int verbosity = 0 );

//@}

//! @name Interior-point l1 minimization methods.
//@{ 

/**
 * \brief Compute the Basis Pursuit ( BP ) sparse solution to 
 \f[ \arg \! \min \|x\|_1\quad\text{subject to}\quad Ax=b\f] by 
 * solving the associated linear program using the primal-dual interior-point
 * method.
 *
 * \param A ( M x N ) matrix of the linear system AX=B
 *
 * \param B ( M x num_rhs ) matrix of the linear system AX=B
 *
 * \param X ( N x num_rhs ) solutions to AX=B
 * 
 * \param primal_dual_tol controls the accuracy of the internal optimzation 
 * alogrithm
 *
 * \param conjugate_gradients_tol Specfies the tolerance of the 
 * the conjugate gradients method used to solve the newton step. If 
 * conjugate_gradients_tol < 0 cholesky factorization will be used.
 * 
 * \param verbosity turn print statements on and off. 
 */
void BP_primal_dual_interior_point_method( RealMatrix &A, const RealVector &b, 
					   RealMatrix &result, 
					   Real primal_dual_tol,
					   Real conjugate_gradients_tol,
					   int verbosity );

/**
 * \brief Compute the Basis Pursuit Denoising sparse solution to 
 \f[ \arg \! \min \|x\|_1\quad\text{subject to}\quad\|Ax-b\|_2\le
 \varepsilon\f] using 
 * the log barrier method to solve the associated quadratic cone problem.
 *
 * \param A ( M x N ) matrix of the linear system AX=B
 *
 * \param B ( M x num_rhs ) matrix of the linear system AX=B
 *
 * \param X ( N x num_rhs ) solutions to AX=B
 *
 * \param epsilon defines the inequality constraint.
 *
 * \param log_barrier_tol controls the accuracy of the internal optimzation 
 * alogrithm   
 *
 * \param conjugate_gradients_tol Specfies the tolerance of the 
 * the conjugate gradients method used to solve the newton step. If 
 * conjugate_gradients_tol < 0 cholesky factorization will be used.
 *
 * \param verbosity turn print statements on and off. 
 * 0: off, 1: warnings on,  2: all print statements on.
 */
void BPDN_log_barrier_interior_point_method( RealMatrix &A, const RealVector &b, 
					     RealMatrix &result, 
					     Real epsilon, 
					     Real log_barrier_tol,
					     Real conjugate_gradients_tol,
					     int verbosity );
//@}

} // namespace Surrogates

#endif // include guard

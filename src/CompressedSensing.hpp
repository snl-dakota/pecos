/**
 * \file CompressedSensing.hpp
 * \author John D. Jakeman
 * \date 2 January 2012
 * \brief Methods used to solve various compressed sensing problems.
 */

#ifndef COMPRESSED_SENSING_HPP
#define COMPRESSED_SENSING_HPP

#include "LinearAlgebra.hpp"
#include "pecos_global_defs.hpp"

namespace Pecos {

//solverType solverTypeCast( int i );

/**
 * \brief Specify a set of options for using the CompressedSensingTool
 */
struct CompressedSensingOptions
{
  //solverType solver; //!< Specify which regression solver to use. See solverType
  short solver; //!< Specify which regression solver to use: see pecos_global_defs
  Real solverTolerance; //!< Specify the internal tolerance of the solver
  Real epsilon;         //!< Specify the residual tolerance of the solver
  Real delta;           //!< Specify the regularization parameter value
  int maxNumIterations; //!< Specify the maximum number of solver iterations
  bool standardizeInputs;  //!< Specify if the inputs need to be standardized
  bool storeHistory;       //<! Specify if the solution history should be stored 
  Real conjugateGradientsTolerance; //<! Specify wether to use conjugate gradients internally to solve newton step in BP and BPDN.  If < 0 cholesky factorization will be used.
  int verbosity;           //!< The verbosity level. 0: off, 1: warnings on,  2: all print statements on.
  int numFunctionSamples; //!< The number of function samples used to construct A and B. Used when A contains gradient information. If zero then numFunctionSamples = A.numRows()

  CompressedSensingOptions() : 
    solver( DEFAULT_LEAST_SQ_REGRESSION ), solverTolerance( -1. ),
    epsilon( 0.0 ), delta( 0.0 ),
    maxNumIterations( std::numeric_limits<int>::max() ), 
    standardizeInputs( false ), storeHistory( false ), 
    conjugateGradientsTolerance( -1 ), verbosity( 0 ), numFunctionSamples( 0 )
  {};

  void print()
  {
    std::cout << "Solver: " << solver << "\n";
    std::cout << "Solver Tolerance: " << solverTolerance << "\n";
    std::cout << "Epsilon: " << epsilon << "\n";
    std::cout << "Delta: " << delta << "\n";
    std::cout << "MaxNumIterations: " << maxNumIterations << "\n";
    std::cout << "StandardizeInputs: " << standardizeInputs << "\n";
    std::cout << "StoreHistory: " << storeHistory << "\n";
    std::cout << "Verbosity: " << verbosity << "\n";
  };
};

typedef std::vector< std::vector<CompressedSensingOptions> > CompressedSensingOptionsList;

/**
 * \class CompressedSensingTool
 * \brief Tool that implements a number of popular compressed sensing algorithms 
 */
class CompressedSensingTool
{
private:
  
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
  int BPDN_compute_central_point( RealMatrix &A, RealVector &b, RealVector &x, 
				  RealVector &u, RealMatrix &AtA, 
				  Real epsilon, Real &tau, 
				  Real newton_tol, int newton_maxiter,
				  Real conjugate_gradients_tol = 1e-8,
				  int verbosity = 0 );

  //@}

public:

  /// Default constructor
  CompressedSensingTool()
  {   
    std::cout.precision( std::numeric_limits<Real>::digits10 );
    std::cout.setf( std::ios::scientific );
  };

  // Deconstructor
  ~CompressedSensingTool(){};

  //! @name Interface tools.
  //@{ 

  /**
   * \brief Wrapper to call any of the compressed sensing methods
   *
   * \param A ( M x N ) matrix of the linear system AX=B
   *
   * \param B ( M x num_rhs ) matrix of the linear system AX=B
   *
   * \param solutions (output) vector containing multiple solutions
   * to AX=B. Each entry in solutions is a ( N x num_rhs ) matrix
   * opts.solver=LS will return only one solution whilst methods such
   * as OMP, LARS, and LASSO will return a history of solutions
   *
   * \param opts specifies the method options
   *
   * \param opts_list specifies the method options that can be used to
   * reproduce all the solutions found.
   */
  void solve( RealMatrix &A, 
	      RealMatrix &B, 
	      RealMatrixList &solutions,
	      CompressedSensingOptions &opts,
	      CompressedSensingOptionsList &opts_list );

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
  void BP_primal_dual_interior_point_method( RealMatrix &A, RealMatrix &B, 
					     RealMatrix &X, 
					     Real primal_dual_tol = 1.e-3,
					  Real conjugate_gradients_tol = 1e-8,
					     int verbosity = 0 );

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
  void BPDN_log_barrier_interior_point_method( RealMatrix &A, RealMatrix &B, 
					       RealMatrix &X, 
					       Real epsilon = 1.e-3, 
					       Real log_barrier_tol = 1.e-3,
					   Real conjugate_gradients_tol = 1e-8,
					       int verbosity  = 0 );
  //@}

  //! @name Greedy methods.
  //@{

  /**
   * \brief Compute a greedy approximation to 
   \f[ \arg \! \min \|x\|_1\quad\text{subject to}\quad\|Ax-b\|_2\le
   \varepsilon\f]
   *
   * If  max_num_iterations > min(M,N) then if the algorithm completes
   * the last solution will be the least squares solution. However sometimes
   * due to numerical roundoff error the algorithm may exit early and 
   * the least squares solution will not be computed.
   *
   * \param A ( M x N ) matrix of the linear system Ax=b
   *
   * \param b ( M x 1 ) vector of the linear system Ax=b
   *
   * \param solutions (output) On exit solution contains the solution
   * at each iteration of the algorithm.
   *
   * \param solution_metrics (output) Contains metrics about the solutions
   * contained in solutions. Specifically for each solution the residual
   * and number of non-zero-terms is stored.
   *
   * \param epsilon defines the inequality constraint.  If set to zero
   * the method will return all possible solutions with sparsity
   * increasing in unit increments from 1 to M
   *
   * \param max_num_iterations specify the maximum number of non-zeros
   * 
   * \param verbosity turn print statements on and off. 
   * 0: off, 1: warnings on,  2: all print statements on.
   */
  void orthogonal_matching_pursuit( RealMatrix &A, RealVector &b, 
				    RealMatrix &solutions,
				    RealMatrix &solution_metrics,
				    Real epsilon = 1.e-3,
		       int max_num_iterations = std::numeric_limits<int>::max(),
				    int verbosity = 0 );

  //@}

  //! @name Homotopy methods.
  //@{

  /**
   * \brief Compute the least angle regression ( and the lasso modification ) 
   * solution to \f[ \arg \! \min \|Ax-b\|^2_2 \quad\text{subject to}\quad
   \|x\|_1\le \tau \f]
   *
   * \param A ( M x N ) matrix of the linear system Ax=b
   *
   * \param b ( M x 1 ) vector of the linear system Ax=b
   *
   * \param epsilon controls the exit condition. If set to zero
   * the method will return the set of solutions for all \f$ \tau \f$
   * with sparsity ranging from 1 to M. If epsilon is non zero
   * then the method will terminate when the residual 
   * \f$ \|Ax-b\|_2<\varepsilon \f$
   *
   * \param solutions (output) On exit solution contains the solution
   * at each iteration of the homotopy algorithm.
   *
   * \param solution_metrics (output) Contains metrics about the solutions
   * contained in solutions. Specifically for each solution the residual
   * and number of non-zero-terms is stored.
   *
   * \param solver specify whether to compute the least angle regression (LARS)
   * or the lasso solution
   *
   * \param max_num_iterations specify the maximum number of iterations.
   * For LARS this is also the maximum number of non-zeros. For lasso
   * the numbre of non-zeros will likely be less than the number of iterations
   *
   * \param verbosity turn print statements on and off. 
   * 0: off, 1: warnings on,  2: all print statements on.
   */
  void least_angle_regression( RealMatrix &A, 
			       RealVector &b, 
			       RealMatrix &solutions,
			       RealMatrix &solution_metrics,
			       Real epsilon = 1e-3, 
			       //solverType solver = LEAST_ANGLE_REGRESSION,
			       short solver = LEAST_ANGLE_REGRESSION,
			       Real delta = 0.0,
	      int max_num_iterations = std::numeric_limits<int>::max(),
			       int verbosity = 0 );

  //@}

  /** \brief Standardize the components of a linear system so that the columns 
   * of A have unit length and the values B have mean zero.
   *
   * \param A the original A matrix
   *
   * \param B the original values (rhs)
   *
   * \param A_stand (output) the standardized A matrix
   *
   * \param B the (output) standardized values (rhs)
   *
   * \param A_column_norms the normalization factors used to standardize A
   *
   * \param B_means the means used to center B.
   *
   * \todo Need to make the columns of A_stand have mean zero. That is 
   * for each column of A I need to subtract its mean and use the
   * new centered column in A_stand. Also need to center the columns of 
   * B_stand. At the moment B_stand = B. Both these methods
   * require special treatment when Ax=B the first column of A is [1,...,1]^T
   * and when it is not.
   */
  void standardize_inputs( RealMatrix &A, RealMatrix &B, 
			   RealMatrix &A_stand, 
			   RealMatrix &B_stand,
			   RealVector &A_column_norms,
			   RealVector &A_column_means,
			   RealVector &B_means ); 
};

} // namespace Pecos

#endif



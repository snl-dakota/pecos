#ifndef SPGL1_HPP
#define SPGL1_HPP

#include "linear_algebra.hpp"
#include "math_tools.hpp"
#include "LinearSystemSolver.hpp" 
//move LinearSolver to new file LinearSolver.hpp"
// and use linear_solvers.hpp to include all types of solvers

namespace Surrogates {

class SPGL1Counters{
public:

  /// The number of matrix vector products using AMatrix
  int numAMatrixProd_;

  /// The number of matrix vector products using AMatrixTrans
  int numAMatrixTransProd_;

  /// The number of newton iterations used to find the BP/BPDN solution
  int newtonIters_;
  
  /// The total number of line iterations used to find the BP/BPDN solution
  int totalLineIters_;
  
  /// The total number of spectral projected gradient iterations used to 
  /// find the approximate Lasso solution
  int SPGIters_;
  
  /// The number of line errors
  int numLineErrors_;

  int numFixedNNZIters_;
  
  SPGL1Counters() : numAMatrixProd_(0), numAMatrixTransProd_(0){initialize();};

  ~SPGL1Counters(){};

  void initialize(){
    numAMatrixProd_ = 0;
    numAMatrixTransProd_ = 0;
    newtonIters_ = 0;
    totalLineIters_ = 0;
    SPGIters_ = 0;
    numLineErrors_ = 0;
    numFixedNNZIters_ = 0;
  } 
};

// Flags for spectral projected gradient iterations
enum { CONTINUE, EXIT_ROOT_FOUND, EXIT_BPSOL1_FOUND, 
       EXIT_BPSOL2_FOUND, EXIT_OPTIMAL, EXIT_ITERATIONS, 
       EXIT_LINE_ERROR, EXIT_SUBOPTIMAL_BP, EXIT_MATVEC_LIMIT,
       EXIT_NO_NNZ_CHANGE};

// Flags for line searches
enum{ LS_CONVERGED, LS_MAX_ITERATIONS_REACHED, LS_NOT_DESCENDING};

/// Compute the l1 norm of a vector vec
inline Real one_norm( const RealVector &vec ){
  int num_entries = vec.length();
  Real norm = 0.;
  for (int i=0; i<num_entries; i++)
    norm += std::abs( vec[i] );
  return norm;
};

/// Compute the dual norm of a vector vec
inline Real dual_norm( const RealVector &vec ){
  int num_entries = vec.length();
  Real dual_norm = 0.;
  for (int i = 0; i<num_entries; i++)
    dual_norm = std::max( dual_norm, std::abs( vec[i] ) );
  return dual_norm;
};

/// Returns the value of a with the sign of b.
/// This function replicates the functionality of Fortran sign function.
inline Real sign_fortran( Real a, Real b){
  return std::abs(a) * sgn(b);
};

/** \brief Sort (descending order) the elements of vec by magnitude and 
 * return the sorted absolute values and the indices that sort vec into
 * descending order
 *
 * \param[in] vec The vector to be sorted (raw values, not absolute values) 
 *
 * \param[out] indices The indices that sort vec into descending order
 * \param[out] abs_vec The sorted absolute values of vec
 */
void magnitude_sort( const RealVector &vec, IntVector &indices,
		     RealVector &abs_vec );

/*class MatrixOperator{
  /// Perform matrix vector multiply
  /// TODO: Add a memory efficient version for PCE which does not create
  /// Amatrix explicitly
  void multiply( const RealVector &vec, Teuchos::ETransp trans,
		 RealVector &result ){
    result.multiply( trans, Teuchos::NO_TRANS, 
		     -1.0, A_copy, vec, 1.0 );
  }
};*/

class SPGL1Solver : public LinearSolver
{
protected:

  
  /// Tolerance controlling size of duality gap at completion
  Real optTol_;
  
  /// Tolerance to which Basis Pursuit solution is computed
  Real bpTol_;

  /// The minimum gradient step size
  Real minGradStepSize_;

  /// The maximum gradient step size
  Real maxGradStepSize_;

  /// The maximum number of iterations for non-monotone line search
  int maxNMLinesearchIters_;
  
  /// The maximum number of iterations for projected line search
  int maxProjLinesearchIters_;

  /// The maximum number of linesearch errors allowed
  int maxLineErrors_;
  
  /// Sufficient descent tolerance for linesearch. 
  Real sufficientDescentTol_;

  /// Relative toleranace for change in objective before invoking another
  /// newton update
  Real newtonUpdateTol_;

  /// Number of previous objective functions that are stored for reference
  int numPrevObjStore_;

  /// The maxiumum number of iterations used to find each lasso solution
  int maxSPGIters_;

  /// The maximum number of spg iterations that can occur with no change in 
  /// the number of nonzero coefficients
  int maxNumFixedNNZIters_;

  /// Machine representation of infinity
  Real inf;

  /// Machine precision
  Real machEps;

  bool solveSingleLasso_;

  SPGL1Counters counters_;

  /// The matrix A of the linear system Ax=b
  RealMatrix AMatrix_;

  ///MatrixOperator AMatrixOp

  /// A vector the length of the solution that is used for projection
  // and adding of vectors
  RealVector solWorkVec_;

  /** \brief find the initial solution for the spgl1 algorithm
   *
   * \param[in] Amatrix The matrix A of the linear system Ax=b
   * \param[in] rhs The RHS b of the linear system Ax=b
   * 
   * \param[out] initial_solution The initial solution A*b
   */
  void find_initial_solution( const RealVector &rhs,
			      RealVector &initial_solution );

  /** \brief Compute the one norm projection of the vector vec onto the 
   * one-norm ball of radius tau
   *
   * \param[in] vec
   * \param[in] tau the radius of the one-norm ball
   *
   * \param[out] vec_proj The one-norm projection of vec onto the one ball 
   * of radius tau
   */
    void project( const RealVector &vec, Real tau, RealVector &vec_proj );

  /** \brief Find the current active set.
   *
   * \param[in] sol The current solution.
   * \param[in] grad The gradient of the current solution.
   * \param[in] nonzero_indices  A vector of primal/dual indicators.
   * \param[in] opt_tol 
   * 
   * \param[out] nnz_sol The no. of non-zero elements in the solution.
   * \param[out] nnz_new_indices The no. of elements in new_nonzero_indices.
   * \param[out] nnz_diff The no. of elements that changed in the support.
   * \param[out] new_nonzero_indices 

   * Only nnz_diff is used to determine exit conditions.
   * The other nnz_ variables are used for subspace minimization
   * which currently have not implemented
   */
  void num_active_variables( const RealVector &sol, 
			     const RealVector &grad, 
			     const IntVector &nonzero_indices,
			     Real opt_tol,
			     int &nnz_sol, int &nnz_new_indices, int &nnz_diff,
			     IntVector &new_nonzero_indices ) const;

  /** \brief Compute matrix vector product Amatrix * vec
   *  result = beta * result + alpha*Amatrix*vec
   *
   * \param[in] vec Vector involved in matrix vector product
   * \param[in] trans Flag specifying whether to compute A'v or Av
   * \param[in] alpha The scaling factor for A'*v
   * \param[in] beata The scaling factor for result
   *
   * \param[out] result The vector resulting from the matrix vector product. 
   *             result may not be empty. If it is empty memory 
   *             will be allocated
   */
  void matrix_multiply( const RealVector &vec, Teuchos::ETransp trans,
			Real alpha, Real beta, RealVector &result );

  /** \brief Projected backtracking linesearch.
   *
   * \param[in] rhs The rhs (b) of the linear system Ax=b
   * \param[in] descent_dir the descent direction of the line search
   * \param[in] max_prev_objs The maximum of all previous objectives of SPG
   * \param[in] lasso_param The lasso parameter tau
   *
   * \param[out] solution The solution found by the line search
   * \param[out] residual The vector b-Ax
   * \param[out] objective The value of the objective evaluated at solution
   * \param[out] proj_linesearch_iters The number of line search iterations taken
   * \param[out] proj_grad_step_len The size of the step found by the linesearch 
   * \param[out] linesearch_exit_code The exit code for the line search
   */
  void spg_projected_linesearch( const RealVector &rhs, 
				 const RealVector &descent_dir, 
				 Real max_prev_objs, Real lasso_param,
				 RealVector &solution,  RealVector &residual,
				 Real &objective, int &proj_linesearch_iters,
				 Real &proj_grad_step_len, 
				 int &linesearch_exit_code );

  /** \brief Nonmonotone line search that does not impose functional decrease 
   *  at every iteration
   *
   * \param[in] rhs The rhs (b) of the linear system Ax=b
   * \param[in] delta_sol The change in the current solution from the previous
   *                      solution
   * \param[in] delta_sol The dot product of the gradient of the objective
                          (gradient) with the change in the current solution  
                          from the previous solution (delta_sol)
   * \param[in] max_prev_objs The maximum of all previous objectives of SPG
   *
   * \param[out] solution The solution found by the line search
   * \param[out] residual The vector b-Ax
   * \param[out] nm_linesearch_iters The number of line search iterations taken
   * \param[out] linesearch_exit_code The exit code for the line search
   *
   * \param[in,out] objective On entry: the value of the objective at the current
   *                          solution. 
   *                          On exit: the value of the objective evaluated at 
   *			    the linesearch solution
   */
  void spg_non_monotone_linesearch( const RealVector &rhs, 
				    const RealVector &delta_sol, 
				    Real gradient_delta_sol, 
				    Real max_prev_objs, 
				    RealVector &solution, RealVector &residual,
				    Real &objective, int &nm_linesearch_iters,
				    int &linesearch_exit_code);

  /*\brief 
   *
   */
  void linesearch( const RealVector &rhs, Real lasso_param, 
		   const RealVector &prev_solution, const RealVector &gradient, 
		   Real gradient_step, const RealVector &prev_objs, 
		   RealVector &solution, Real &objective, RealVector &residual, 
		   Real &proj_grad_step_len, int &exit_code );

  /** \brief Compute the gradient of the new solution
   *
   * \param[in] solution The current solution
   * \param[in] prev_solution The previous solution
   * \param[in] prev_gradient The previous gradient
   * \param[in] residual The vector b-Ax
   *
   * \param[out] gradient The gradient of the current solution
   * \param[out] gradient_step The size of the gradient step
   */
  void update_gradient( const RealVector &solution, 
			const RealVector &prev_solution,
			const RealVector &prev_gradient, 
			const RealVector &residual, 
			RealVector &gradient, Real &gradient_step);

  /// Print cause of algorithm terminating
  void print_termination_message( int exit_code );

  void test_exit_conditions( Real objective, Real lasso_param, 
			     const RealVector &solution, 
			     const RealVector & gradient, 
			     const RealVector &rhs, const RealVector &residual, 
			     Real residual_norm, Real lagrange_multiplier, 
			     Real residual_tol, 
			     const IntVector &nonzero_indices, 
			     Real rhs_norm, Real prev_objective,
			     int &exit_code, Real proj_grad_step_len,
			     bool &update_lasso_param);

  void test_exit_conditions_for_multiple_lasso_param( Real duality_gap_rel, 
						      Real objective_rel_error, 
						      Real residual_rel_error, 
						      Real residual_norm, 
						      Real rhs_norm, 
						      Real lagrange_multiplier, 
						      Real residual_tol, 
						      Real objective, 
						      Real prev_objective, 
						      bool &update_lasso_param, 
						      int &exit_code);
    
public:
  SPGL1Solver() : optTol_(1e-4), bpTol_(1e-6), minGradStepSize_(1e-16),
		  maxGradStepSize_(1e5), maxNMLinesearchIters_(10),
		  maxProjLinesearchIters_(10), maxLineErrors_(10),
		  sufficientDescentTol_(1e-4),
		  newtonUpdateTol_(1e-6), numPrevObjStore_(3), maxSPGIters_(0),
		  maxNumFixedNNZIters_(std::numeric_limits<double>::max()),
		  inf(std::numeric_limits<double>::max()),
		  machEps(std::numeric_limits<double>::epsilon()),
		  solveSingleLasso_(true){};

  ~SPGL1Solver(){};

  Real spgl1( const RealMatrix &Amatrix, const RealVector &rhs, 
	      Real lasso_param, Real residual_tol, RealVector &result );

  /** \brief Solve the lasso problem (possibly approximately) using 
   * the spectral projected gradient method
   *
   * \param[in] rhs The rhs (b) of the linear system Ax=b

   * \param[in,out] residual The vector b-Ax
   * \param[in,out] gradient The gradient of the BPDN objective function
   * \param[in,out] solution On entry: thee  solution is an initial guess. 
                             On exit: the solution of the lasso problem. 
   * \param[in,out] nonzero_indices
   * \param[in,out] prev_objectives
   */
  void spectral_projected_gradient( const RealVector &rhs, RealVector &residual,
				    RealVector &gradient, Real lasso_param, 
				    Real &objective, Real &lagrange_multiplier,
				    RealVector &solution, Real residual_tol, 
				    IntVector &nonzero_indices, Real &rhs_norm,
				    Real &prev_objective, 
				    Real &proj_grad_step_len, int &exit_code,
				    Real &best_objective, Real &gradient_step, 
				    RealVector &prev_objs );

  /** \brief Apply defaults and user options and check they are consistent.
   *
   * \param[in] num_terms The number of columns in AMatrix_
   *
   * \param[in,out] lasso_param The regularization parameter of the lasso problem
   *                            If < 0 then apply default (0)
   * \param[in,out] residual_tol The regularization paramater of BPDN. If < 0 
   *                             then apply default (0)
   */
  void apply_user_options( int num_cols, Real &lasso_param, Real &residual_tol );

  void solve( RealMatrix &A, RealMatrix &B, RealMatrix &result_0, 
	      RealMatrix &result_1 );

  void set_max_iters( int max_iters );
};

} // namespace Surrogates

#endif // SPGL1_HPP

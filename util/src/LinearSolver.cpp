#include "LinearSolver.hpp"
#include "basis_pursuit.hpp"
#include "least_angle_regression.hpp"
#include "linear_solvers.hpp"
#include "orthogonal_matching_pursuit.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

namespace Surrogates {




Teuchos::ParameterList LinearSolver::pecos_default_params()
{
  Teuchos::ParameterList params;

  params.set<unsigned>("Linear Solver Type", SVD_LEAST_SQ_REGRESSION);
  params.set<Real>("Solver Tolerance", -1.0);
  // The control formerly known as epsilon: 
  params.set<Real>("Residual Tolerance", 0.0);
  // TODO: remove delta from LARS
  params.set<Real>("Delta", 0.0);
  params.set<int>("Max Iters", std::numeric_limits<int>::max());
  params.set<bool>("Normalize Inputs", false);
  params.set<int>("Verbosity", 0);
  params.set<bool>("Store History", false);
  // TODO: remove CG tol once BP* are retired
  params.set<Real>("CG Tolerance", -1.0);
  params.set<int>("Num Function Samples", 0);
  params.set<bool>("Normalize Choice", false);

  return params;
}

Teuchos::ParameterList LinearSolver::default_params()
{
  Teuchos::ParameterList params;

  params.set<unsigned>("Linear Solver Type", SVD_LEAST_SQ_REGRESSION);
  params.set<Real>("Solver Tolerance", -1.0);
  // The control formerly known as epsilon: 
  params.set<Real>("Residual Tolerance", 0.0);
  params.set<Real>("Delta", 0.0);
  params.set<int>("Max Iters", std::numeric_limits<int>::max());
  params.set<bool>("Normalize Inputs", true);
  params.set<int>("Verbosity", 0);
  params.set<bool>("Store History", false);
  params.set<Real>("CG Tolerance", -1.0);
  params.set<int>("Num Function Samples", 0);
  params.set<bool>("Normalize Choice", false);

  return params;
}


LinearSolver::LinearSolver()
{  /* empty ctor*/ }

LinearSolver::~LinearSolver()
{  /* empty dtor*/ }


// By default A is copied to avoid destorying, B is viewed, probably
// can be const, except for Teuchos::View semantics
void LinearSolver::solve(const RealMatrix &A, RealMatrix &B, RealMatrix &result_0, 
			 RealMatrix &result_1, Teuchos::ParameterList & params)
{
  // TODO: don't allow default, or warn?!?
  unsigned solver_type = params.get<unsigned>("Linear Solver Type");

  // options shared among solvers
  int verbosity               = params.get("Verbosity", 0);
  bool normalize_inputs       = params.get("Normalize Inputs", false); // BMA: default should be true

  // only used by OMP and LARS: put in switch?
  bool normalize_choice       = params.get("Normalize Choice", false);
  int max_iters               = params.get("Max Iters", std::numeric_limits<int>::max());
  Real residual_tol           = params.get("Residual Tolerance", 1.e-6);


  // always copy A and only support one RHS b for now (all sub-solvers
  // only support one for now)
  RealMatrix A_copy( Teuchos::Copy, A, A.numRows(), A.numCols() );
  if ( B.numCols() != 1 )
    throw( std::runtime_error(" OMPSolver::solve() B must be a vector") );
  RealVector b( Teuchos::View, B[0], B.numRows() );

  RealVector column_norms;
  if ( normalize_inputs )
    normalise_columns( A_copy, column_norms );

  switch(solver_type) {
    
    //case SVD_LEAST_SQ_REGRESSION: {
    // Real solver_tol             = params.get("Solver Tolerance", 1.e-6);
    //svd_lsq_solve(A_copy, b, result_0, result_1, solver_tol);
    //break;
    //}

    //case EQ_CONS_LEAST_SQ_REGRESSION: {
    //  int num_primary_eqs         = params.get<int>("Num Primary Equations");
    //eq_cons_lsq_solve(A_copy, b, result_0, result_1, num_primary_eqs);
    //break;
    //}

  case ORTHOG_MATCH_PURSUIT: {

    /**
     * \brief Find the solution min ||x||_0 such that |AX = B||_2 < eps
     */
    
    // Discuss with John .... ordering

    /// enforces a set of columns to be chosen first
    IntVector ordering = params.get("Column Ordering", IntVector());

    // Need to handle ordering - RWH
    orthogonal_matching_pursuit<Real>( A_copy, b, result_0, result_1,
				       residual_tol, max_iters, verbosity,
				       ordering, normalize_choice, true, 100 );

    break;
  }

  case LASSO_REGRESSION:
  case LEAST_ANGLE_REGRESSION: {
    /**
     * \brief Find the solution min ||x||_0 such that |AX = B||_2 < eps
     */

    Real delta = params.get("Delta", 0.0);

    int maxNNZ = params.get<int>("Max Nonzeros"); // no default means the parameter must be set

    least_angle_regression( A_copy, b, result_0, result_1,
			    residual_tol, solver_type, delta,
			    max_iters, maxNNZ,
			    verbosity, normalize_choice, false, true, 500 );
    break;
  }
  case BASIS_PURSUIT: {
    basis_pursuit_solve(A_copy, b, result_0, result_1, params);
    break;
  }
  case BASIS_PURSUIT_DENOISING: {
    bp_dn_solve(A_copy, b, result_0, result_1, params);
    break;
  }
  default: {
    std::cerr << "No default regression anymore!\n";
    break;
  }
  }

  if ( normalize_inputs )
    adjust_coefficients( column_norms, result_0 );

}

 
// BMA TODO: merge specialized multi-RHS solvers from csopts_toreview....

/** Default implementation of multiple RHS solve that just calls the
    derived solver for each col in B.  TODO: currently solve only
    supports 1 RHS... */
void LinearSolver::
multi_rhs_solve(const RealMatrix& A, const RealMatrix& B, 
		std::vector<RealMatrix>& solutions, 
		Teuchos::ParameterList & params)
{
  int num_cols = B.numCols();
  solutions.resize(num_cols);
  for (int rhs_ind = 0; rhs_ind<num_cols; ++rhs_ind) {
    // isolate the const_cast with a copy
    RealVector b = Teuchos::getCol(Teuchos::Copy,
				   const_cast<RealMatrix&>(B), rhs_ind);
    // TODO: Ask JDJ what to do with result_1 for say OMP...
    RealMatrix result_1;
    solve(A, b, solutions[rhs_ind], result_1, params);
  }
}



}

#ifndef LINEARSOLVER_HPP
#define LINEARSOLVER_HPP

#include "teuchos_data_types.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace Surrogates {


 /**
 * \brief Find the solution min ||x||_0 such that |AX = B||_2 < eps
 */
void svd_lsq_solve( const RealMatrix &A, const RealVector &B, RealMatrix &result_0,
                   RealMatrix &result_1, Real solver_tol );


//----------------------------------------------------------------

/**
 * \brief Find the solution min ||x||_0 such that |AX = B||_2 < eps
 */
void eq_cons_lsq_solve( const RealMatrix &A, const RealVector &B, 
                       RealMatrix &result_0, RealMatrix &result_1,
                       const int num_primary_eqs );
  

class LinearSolver
{
public:

  // BMA: consider naming this enum type and moving inside LinearSolver;
// eventually strongly typed with C++11
enum {
  SVD_LEAST_SQ_REGRESSION, EQ_CONS_LEAST_SQ_REGRESSION, ORTHOG_MATCH_PURSUIT,
  LASSO_REGRESSION, LEAST_ANGLE_REGRESSION, BASIS_PURSUIT, 
  BASIS_PURSUIT_DENOISING
};


  /// default constructor
  LinearSolver();

  /// virtual destructor as we intend this class to be specialized
  virtual ~LinearSolver();

  static Teuchos::ParameterList pecos_default_params();

  static Teuchos::ParameterList default_params();


  /**
   * \brief Find a regularized solution to AX = B
   */
  virtual void solve( const RealMatrix &A, RealMatrix &B, RealMatrix &result_0, 
		      RealMatrix &result_1, Teuchos::ParameterList & params );
  
  /** A multiple RHS solver that returns a compressed sensing solution
      matrix for each RHS.  TODO: call this just solve and derived
      classes will specialize this if they support multi-RHS, e.g.,
      SVD regression. */
  virtual void multi_rhs_solve(const RealMatrix& A, 
			       const RealMatrix& B, 
			       std::vector<RealMatrix>& solutions, 
			       Teuchos::ParameterList & params);

  //virtual void solve_using_points( RealMatrix &build_points, RealMatrix &B, 
  //      			   RealMatrix &result_0, 
  //      			   RealMatrix &result_1 )
  //{
  //  std::string msg = "solve_using_points() Has not been implemented for ";
  //  msg += "this class.";
  //  throw( std::runtime_error( msg ) );
  //};

  static void normalise_columns( RealMatrix &A, RealVector &result )
  {
    int M = A.numRows(), N = A.numCols();
    result.sizeUninitialized( N );
    for ( int i = 0; i < N; i++ )
      {
        RealVector col( Teuchos::View, A[i], M );
        result[i] = col.normFrobenius();
        col *= 1./result[i];
      }
  }

  static void adjust_coefficients( const RealVector &normalisation_factors, 
				   RealMatrix &coefficients )
  {
    int num_coeff = coefficients.numRows(), num_qoi = coefficients.numCols();
    for ( int i = 0; i < num_qoi; i++ )
      {
        for ( int j = 0; j < num_coeff; j++ )
          coefficients(j,i) /= normalisation_factors[j];
      } 
  }

};

} // namespace Surrogates

#endif //LINEARSOLVER_HPP

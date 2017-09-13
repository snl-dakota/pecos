#include "basis_pursuit.hpp"
#include "LinearSolver.hpp"
#include "linear_algebra.hpp"
#include "math_tools.hpp"  // for sum(RealVector)
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Surrogates {

/** High-level solver wrappers for BP approaches */

// BMA: These are to be removed, so just passing a PL for now

void basis_pursuit_solve( RealMatrix &A, const RealVector &b, RealMatrix &result_0, 
			  RealMatrix &result_1, Teuchos::ParameterList & params)
{
  int verbosity_              = params.get("Verbosity", 0);
  bool normaliseInputs_      = params.get("Normalize Inputs", false); // BMA: default should be true

  Real solverTol_ = params.get("Solver Tolerance", 1.0e-6);
  Real conjugateGradientTol_ = params.get("CG Tolerance", 1.0e-6);

  RealVector column_norms;
  if ( normaliseInputs_ )
    LinearSolver::normalise_columns( A, column_norms );

  BP_primal_dual_interior_point_method( A, b, 
					result_0, 
					solverTol_, 
					conjugateGradientTol_,
					verbosity_ );

  if ( normaliseInputs_ )
    LinearSolver::adjust_coefficients( column_norms, result_0 );

  result_1.shapeUninitialized( 2, 1 );
  result_1(0,0) = 0.;
  int num_non_zeros = 0;
  for ( int i = 0; i < result_0.numRows(); i++ )
    if ( std::abs( result_0(i,0) ) > std::numeric_limits<double>::epsilon() )
      num_non_zeros++;
  result_1(1,0) = num_non_zeros;
}

void bp_dn_solve( RealMatrix &A, const RealVector &b, RealMatrix &result_0, 
		  RealMatrix &result_1, Teuchos::ParameterList & params)
{
  int verbosity_              = params.get("Verbosity", 0);
  bool normaliseInputs_      = params.get("Normalize Inputs", false); // BMA: default should be true

  Real solverTol_ = params.get("Solver Tolerance", 1.0e-6);
  Real conjugateGradientTol_ = params.get("CG Tolerance", 1.0e-6);

  RealVector residualTols_ = params.get("Residual Tolerance List", RealVector());

  if ( residualTols_.length() <= 0 )
    throw( std::runtime_error(" BPDNSolver::solve() set residual tols") );

  RealVector column_norms;
  if ( normaliseInputs_ )
    LinearSolver::normalise_columns( A, column_norms );

  result_0.shapeUninitialized( A.numCols(), residualTols_.length() );
  result_1.shapeUninitialized( 2, residualTols_.length() );
  for ( int j = 0; j < residualTols_.length(); j++ )
    {
      RealMatrix x;
      BPDN_log_barrier_interior_point_method( A, b, 
					      x, 
					      residualTols_[j],
					      solverTol_, 
					      conjugateGradientTol_,
					      verbosity_ );
	
      if ( normaliseInputs_ )
	LinearSolver::adjust_coefficients( column_norms, x );
	
      result_1(0,j) = residualTols_[j];
      int num_non_zeros = 0;
      for ( int i = 0; i < result_0.numRows(); i++ )
	{
	  if ( std::abs( x(i,0) ) > std::numeric_limits<double>::epsilon() )
	    num_non_zeros++;
	  result_0(i,j) = x(i,0);
	}
      result_1(1,j) = num_non_zeros;
    }
}


/** Implementation functions for BP */

Real BP_surrogate_duality_gap( RealVector &primal_residual,
			       RealVector &f_1, 
			       RealVector &f_2, 
			       RealVector &lambda_1, 
			       RealVector &lambda_2, 
			       RealVector &Atv, 
			       Real mu, 
			       Real primal_dual_tol, 
			       Real &t, 
			       Real &slackness_norm )
{
  int M(primal_residual.numRows() ), N( f_1.numRows() );

  Real sdg( 0.0 );
  sdg = - f_1.dot( lambda_1 ) - f_2.dot( lambda_2 );
  t = mu * 2 * N / sdg;
       
  slackness_norm = 0.0;
  for ( int j = 0; j < N; j++ )
    {
      Real centrality_residual_j  = -lambda_1[j] * f_1[j] - 1.0 / t;
      Real centrality_residual_nj = -lambda_2[j] * f_2[j] - 1.0 / t;
      Real dual_residual_j = lambda_1[j] - lambda_2[j] + Atv[j];
      Real dual_residual_nj = 1.0 - lambda_1[j] - lambda_2[j];
      slackness_norm += centrality_residual_j * centrality_residual_j;
      slackness_norm += centrality_residual_nj * centrality_residual_nj;
      slackness_norm += dual_residual_j * dual_residual_j;
      slackness_norm += dual_residual_nj * dual_residual_nj;
    }
  for ( int i = 0; i < M; i++ )
    {
      slackness_norm += primal_residual[i] * primal_residual[i];
    }
  slackness_norm = std::sqrt( slackness_norm );

  return sdg;
};

void BP_primal_dual_interior_point_method( RealMatrix &A, 
					   const RealVector &b, 
					   RealMatrix &result,
					   Real primal_dual_tol, 
					   Real cg_tol, 
					   int verbosity )
{
  Teuchos::LAPACK<int, Real> la;
  
  // Extract matrix shapes
  int M( A.numRows() ), N( A.numCols() );


  // Initialise memory for the solutions
  result.shapeUninitialized( N, 1 );
  RealVector x( N, false );
 
  // Basis Pursuit algorithm parameters
  int max_iter( 30 );  // maximum number of primal dual iterations
  Real alpha( 1.e-2 ); // used in residual termination condition
  Real beta( 0.5 );    // used in back tracking procedure
  Real mu( 10.0 );     // used to update slackness condition

  // Allocate memory for the elements use to construct the Hessian
  // of the objective function
  RealVector f_1( N, false ), f_2( N, false ), lambda_1( N, false ), 
    lambda_2( N, false ), lambda_diff( N, false ), newton_step_rhs( M, false );
  RealMatrix newton_step_pos_def_matrix( M, M );
  
  // Allocate memory for storing and computing newton steps
  RealVector Atdv( N, false ), dx( N, false ), du( N, false ), 
    dlambda_1( N, false ), dlambda_2( N, false ), dv( M, false ), 
    Adx( M, false );

  // Allocate memory for the updated variables to be temporarily stored
  RealVector x_new( N, false ), u_new( N, false ), v_new( M, false ), 
    r_new( M, false ), f_1_new( N, false ), f_2_new( N, false ), 
    Atv_new( N, false ), lambda_1_new( N, false ), lambda_2_new( N, false );

  //----------------------------------------------------------//
  // Check whether starting point is in the feasiable region  //
  // If not use the least squares solution as new start point //
  //----------------------------------------------------------//
  RealVector r(Teuchos::Copy, b.values(), b.length());
  //Not necessary because X is set to zero( uninitialized ) above
  //r.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, x, -1.0 );

  Real frob_norm_of_b = b.normFrobenius();
  Real frob_norm_of_r = r.normFrobenius();
  
  Real initial_guess_tol = ( cg_tol < 0 ) ? 1e-8 : cg_tol;
  if ( ( frob_norm_of_r * frob_norm_of_r ) /
       ( frob_norm_of_b * frob_norm_of_b ) > initial_guess_tol )
    {
      if ( verbosity > 0 )
	{
	  std::cout << "BP_primal_dual_interior_point_method() ";
	  std::cout << "Initial guess is not feasiable. ";
	  std::cout << "Computing least squares estimate\n";
	}

      // Solve AX = b
      int rank;
      RealVector singular_values; // singular values
      svd_solve( A, b, x, singular_values, rank );
      // Compute reciprocal of condition number
      Real rcond = singular_values[singular_values.length()-1] / 
	singular_values[0];
      if ( rcond < 1e-14 )
	{
	  std::string msg = "BP_primal_dual_interior_point_method() ";
	  msg += "A is ill conditioned. Cannot find starting point";
	  throw( std::runtime_error( msg ) );
	}
      // Compute the primal_residual $r_\mathrm{pri} = Ax-b$;
      r.assign( b );
      r.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, x, -1.0 );
    }
 
  //-------------------------------------------------------------------//
  // Convert Basis Pursuit problem into traditional linear programming //
  // problem.                                                          //
  // Ensure the constraints $f_{u_1,i}$ and $f_{u_2,i}$ are always     //
  // negative by choosing $u_i=0.95x_i+0.10\max_i{x}$, $i=1,\ldots,N$  //
  //-------------------------------------------------------------------//
  RealVector u( N, false );

  Real x_max = x.normInf();
  for ( int j = 0; j < N; j++ )
    {
      u[j] = 0.95 * std::fabs( x[j] ) + 0.10 * x_max;
      f_1[j] = x[j] - u[j];
      f_2[j] = -x[j] - u[j];
      lambda_1[j] = -1.0 / f_1[j];
      lambda_2[j] = -1.0 / f_2[j]; 
      lambda_diff[j] = lambda_1[j] - lambda_2[j];
    }

  // Compute the dual variable v
  RealVector v( M, false ), Atv( N, false );
  v.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, -1.0 , A, lambda_diff, 
	      0.0 );
  Atv.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, A, v, 0.0 );

  // Compute surrogate duality gap (sdg) and update slackness condition (t)
  Real t, slackness_norm;
  Real sdg = BP_surrogate_duality_gap( r, f_1, f_2, lambda_1, lambda_2, 
				       Atv, mu, primal_dual_tol,
				       t, slackness_norm );   

  if ( verbosity > 1 )
    {
      std::cout << "At Initialisation:\n";
      std::cout << "\tObjective: " << sum( u ) << "\n";
      std::cout << "\tPrimal-dual gap: " << sdg << "\n";
      std::cout << "\tPrimal residual: ";
      std::cout << r.normFrobenius() << "\n";
    }

  //------------------------------------------//
  // Iterate though the primal-dual algorithm //
  //------------------------------------------//
  RealVector z_1( N, false ), z_2( N, false ), D_1( N, false ), 
    D_2( N, false ), D_3( N, false ), z_3( M, false ), tmp1( N, false );
  RealMatrix tmp2( N, M, false );
      
  int primal_dual_iter( 0 );
  bool done = ( ( sdg < primal_dual_tol ) || (primal_dual_iter >= max_iter));
  while ( !done )
    {
      primal_dual_iter++;
      
      //-------------------------------------------------------//
      // Set up linear system to compute newton step direction //
      //-------------------------------------------------------//
      for ( int j = 0; j < N; j++ )
	{
	  z_1[j] = - ( - 1.0 / f_1[j] + 1.0 / f_2[j] ) / t - Atv[j];
	  z_2[j] = - 1.0 - ( 1.0 / f_1[j] + 1.0 / f_2[j] ) / t;
	  D_1[j] = -lambda_1[j] / f_1[j] - lambda_2[j] / f_2[j];
	  D_2[j] =  lambda_1[j] / f_1[j] - lambda_2[j] / f_2[j];
	  D_3[j] =  D_1[j] - D_2[j] * D_2[j] / D_1[j];
	  tmp1[j] =  z_1[j] / D_3[j] - z_2[j] * D_2[j] / (D_3[j] * D_1[j]);
	}
      for ( int i = 0; i < M; i++ )
	{
	  z_3[i] = -r[i];
	  newton_step_rhs[i] = z_3[i];
	}

      newton_step_rhs.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
				1.0, A, tmp1, -1.0 );

      //-------------------------------------------------------//
      // Set up linear system to compute newton step direction //
      //-------------------------------------------------------//

      // Todo I neet to access elements of At. I do this by accessing
      // underlying data of A and utilising flat 1D structure. Can this 
      // be done a different way
      Real *A_matrix;
      A_matrix = A.values();
      for ( int j = 0; j < N; j++ )
	{
	  for ( int i = 0; i < M; i++ )
	    {
	      tmp2(j,i) = A_matrix[j*M+i] / D_3[j];
	    }
	}
      newton_step_pos_def_matrix.multiply( Teuchos::NO_TRANS, 
					   Teuchos::NO_TRANS, 
					   1.0, A, tmp2, 0.0 );

      int info;
      Real rcond( -1. );
      Real r_norm;
      if ( cg_tol < 0 )
	{
	  // Compute the direction of the newton step for v
	  info = cholesky_solve( newton_step_pos_def_matrix, 
				 newton_step_rhs, dv, rcond );
	  if ( info > 0 )
	    {
	      if ( verbosity > 0 )
		{
		  std::cout << "BP_primal_dual_interior_point_method() ";
		  std::cout << "The Hessian matrix is no longer positive ";
		  std::cout << "definite. If epsilon < 1e-8 this is most ";
		  std::cout << "likely due to roundoff error. ";
		  std::cout << "Try a larger epsilon. ";
		  std::cout << "Returning the last solution.\n";
		}
	      x_new.assign( x );
	      break;
	    }
	  if ( rcond < 1.e-14 )
	    {
	      if ( verbosity > 0 )
		{
		  std::cout << "BP_primal_dual_interior_point_method() ";
		  std::cout << "Matrix ill-conditioned. ";
		  std::cout << "Returning the last solution.\n";
		}
	      x_new.assign( x );
	      break;
	    }
	}
      else
	{
	  dv = 0.0;
          int iters;
	  info = conjugate_gradients_solve( newton_step_pos_def_matrix, 
					    newton_step_rhs, dv, r_norm, 
                                            iters, cg_tol, M, verbosity );
	  if ( ( r_norm > 0.5 ) || ( info > 1 ) )
	    {
	      if ( verbosity > 0 )
		{
		  std::cout << "BP_primal_dual_interior_point_method() ";
		  std::cout << "Matrix ill-conditioned. ";
		  std::cout << "Returning the last solution.\n";
		}
	      x_new.assign( x );
	      break;
	    }
	}


      //--------------------------//
      // Compute newton step size //
      //--------------------------//
      Atdv.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 
		     1.0, A, dv, 0.0 );
      Real step_size( 1.0 );
      for ( int j = 0; j < N; j++ )
	{
	  dx[j] = ( z_1[j] - z_2[j] * D_2[j] / D_1[j] - Atdv[j] ) / D_3[j];
	  du[j] = ( z_2[j] - D_2[j] * dx[j] ) / D_1[j];
	  dlambda_1[j] = lambda_1[j] / f_1[j] * ( -dx[j] + du[j] )
	    - lambda_1[j] - 1.0 / ( t * f_1[j] );
	  dlambda_2[j] = lambda_2[j] / f_2[j] * ( dx[j] + du[j] )	
	    - lambda_2[j] - 1.0 / ( t * f_2[j] );
	  
	  if ( dlambda_1[j] < 0 )
	    {
	      step_size = std::min( step_size, -lambda_1[j] / dlambda_1[j] );
	    }
	  if ( dlambda_2[j] < 0 )
	    {
	      step_size = std::min( step_size, -lambda_2[j] / dlambda_2[j] );
	    }
	  if (  ( dx[j] - du[j] ) > 0 )
	    {
	      step_size = std::min( step_size, -f_1[j] / ( dx[j] - du[j] ) );
	    }
	  if (  ( -dx[j] - du[j] ) > 0 )
	    {
	      step_size = std::min( step_size, -f_2[j] / ( -dx[j] - du[j] ));
	    }
	}
      step_size *= 0.99;

      Adx.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, dx, 0.0 );

      //-------------------------------------------------------------------//
      // Conduct line search to ensure the norm of the residuals decreases //
      // sufficiently                                                      //
      //-------------------------------------------------------------------//
      bool sufficient_decrease( false );
      // the iteration count of the backward line seach
      int line_search_iter( 0 ); 

      // Continue while the norm of the residuals has not decreased 
      // sufficiently for any of the quantities of interest.
      while ( !sufficient_decrease )
	{
	  sufficient_decrease = true;
	  Real new_slackness_norm( 0.0 );
	  for ( int j = 0; j < N; j++ )
	    {
	      x_new[j]     =  x[j] + step_size * dx[j];
	      u_new[j]     =  u[j] + step_size * du[j];
	      Atv_new[j]   =  Atv[j] + step_size * Atdv[j];
	      lambda_1_new[j] =  lambda_1[j] + step_size * dlambda_1[j];
	      lambda_2_new[j] =  lambda_2[j] + step_size * dlambda_2[j];
	      f_1_new[j]   =  x_new[j] - u_new[j];
	      f_2_new[j]   = -x_new[j] - u_new[j];

	      // 0:N-1
	      new_slackness_norm += ( lambda_1_new[j] - lambda_2_new[j] + Atv_new[j] ) * ( lambda_1_new[j] - lambda_2_new[j] + Atv_new[j] );
	      // N:2N-1
	      new_slackness_norm += ( 1.0 - lambda_1_new[j] - lambda_2_new[j] ) * ( 1.0 - lambda_1_new[j] - lambda_2_new[j] );
	      // 2N:3N-1
	      new_slackness_norm += ( -lambda_1_new[j] * f_1_new[j] - 1.0 / t ) * ( -lambda_1_new[j] * f_1_new[j] - 1.0 / t );
	      // 3N:4N-1
	      new_slackness_norm += ( -lambda_2_new[j] * f_2_new[j] - 1.0 / t ) * ( -lambda_2_new[j] * f_2_new[j] - 1.0 / t );
	    }
      
	  for ( int i = 0; i < M; i++ )
	    {
	      v_new[i] = v[i] + step_size * dv[i];
	      r_new[i] = r[i] + step_size * Adx[i];
	      new_slackness_norm += r_new[i] * r_new[i];
	    }

	  new_slackness_norm = std::sqrt( new_slackness_norm );
	  if ( new_slackness_norm > ( 1.0 - alpha * step_size ) * slackness_norm )
	    {
	      sufficient_decrease = false;
	    }
	  step_size *= beta;
  
	  line_search_iter++;
	  if ( line_search_iter > 30 )
	    {
	      if ( verbosity > 0 )
		{
		  std::cout << "BP_primal_dual_interior_point_method() ";
		  std::cout << "Line search failed. Returning the last ";
		  std::cout << "solution.\n";
		}
	      x_new.assign( x );
	      info = 3;
	      break;
	    }
	}
  
      for ( int j = 0; j < N; j++ )
	{
	  x[j] =  x_new[j];
	  u[j] =  u_new[j];
	  f_1[j] = f_1_new[j];
	  f_2[j] = f_2_new[j];
	  lambda_1[j] = lambda_1_new[j];
	  lambda_2[j] = lambda_2_new[j];
	  Atv[j] = Atv_new[j];
	}
      for ( int i = 0; i < M; i++ )
	{
	  r[i] = r_new[i];
	  v[i] = v_new[i];
	}

      // Update surrogate duality gap and the slackness condition
      Real sdg = BP_surrogate_duality_gap( r, f_1, f_2, lambda_1, 
					   lambda_2, Atv,
					   mu, primal_dual_tol, t,
					   slackness_norm );
  
  
      done = ( ( sdg < primal_dual_tol ) || ( primal_dual_iter >= max_iter));

      if ( verbosity > 1 )
	{
	  std::cout << "Newton iteration: " << primal_dual_iter << "\n";
	  std::cout << "\tObjective: " << sum( u ) << "\n";
	  std::cout << "\tPrimal-dual gap: " << sdg << "\n";
	  std::cout << "\tPrimal residual: ";
	  std::cout << r.normFrobenius() << "\n";
	  if ( cg_tol < 0 )
	    std::cout << "\tCondition number: " << rcond << "\n";
	  else
	    std::cout << "\tConjugate gradients residual: " << r_norm <<"\n";
	}
    }
  RealVector result_view( Teuchos::View, result[0], N );
  result_view.assign( x );
};

int BPDN_compute_central_point( RealMatrix &A, 
				const RealVector &b, 
				RealVector &x, 
				RealVector &u,
				RealMatrix &AtA, 
				Real epsilon, 
				Real &t, 
				Real newton_tol, 
				int newton_max_iter, 
				Real cg_tol, 
				int verbosity )
{
  int newton_info( 0 ); // No error

  // line search parameters
  Real alpha ( 0.01 );
  Real beta ( 0.5 );  

  // size of solution vector
  int N ( A.numCols() );

  // Number of pts in used to construct A
  int M  ( A.numRows() );

  // Compute the residual
  RealVector r(Teuchos::Copy, b.values(), b.length());
  r.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, x, -1.0 );

  // Set up variables necessary for initialising loop.
  // The objective function is
  // sum(u) + (-1/t) ( log(-f_1) + log(-f_2) + log(f3) )
  RealVector f_1( N, false ), f_2( N, false );
  Real objective( 0.0 ); 
  for ( int j = 0; j < N; j++ )
    {
      f_1[j] = x[j] - u[j];  // constraint x-u <= 0
      f_2[j] = -x[j] - u[j]; // constraint -x-u <= 0
      objective += u[j] - ( std::log( -f_1[j] ) + std::log( -f_2[j] ) ) / t;
    }
  //constraint ( ||Ax-b||_2^2-epsilon^2 ) <= 0
  Real f3 = 0.5 * ( r.dot( r ) - epsilon * epsilon );
  objective -= std::log( -f3 ) / t;

  // Allocate memory for the elements use to construct the Hessian
  // of the objective function
  RealVector D_1( N, false ), D_2( N, false ), D_3( N, false ), 
    z_1( N, false ), z_2( N, false );
  
  // Allocate memory for the gradient and Hessian of the 
  // objective function
  RealMatrix newton_step_pos_def_matrix( N, N, false );
  RealVector newton_step_rhs( N, false );

  // Allocate memory for the newton step
  RealVector dx( N, false ), du( N, false );

  // Allocate memory for the updated variables to be temporarily stored
  RealVector x_new( N, false ), u_new( N, false ), r_new( M, false ), 
    f_1_new( N, false ), f_2_new( N, false ), Atr( N, false ), Adx( M, false );

  int newton_iter( 0 );
  bool done ( false );
  while ( !done )
    {
      // Form the elements of the Jacobian
      Atr.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, A, r, 0.0 );
      newton_step_pos_def_matrix.multiply( Teuchos::NO_TRANS, Teuchos::TRANS, 
					   1.0, Atr, Atr, 0.0 );
      for ( int j = 0; j < N; j++ )
	{
	  Real f_1_inv = 1.0 / f_1[j],  f_2_inv = 1.0 / f_2[j];
	  
	  z_1[j] = f_1_inv - f_2_inv + Atr[j] / f3;
	  z_2[j] = -t - f_1_inv - f_2_inv;
	  
	  D_1[j] =  f_1_inv * f_1_inv + f_2_inv * f_2_inv;
	  D_2[j] = -f_1_inv * f_1_inv + f_2_inv * f_2_inv;
	  D_3[j] = D_1[j] - D_2[j] * D_2[j] / D_1[j]; 
	  newton_step_rhs[j]  = z_1[j]  - z_2[j] * D_2[j] / D_1[j]; 
	  
	  Real f3_inv = 1.0 / f3;
	  for ( int i = 0; i < N; i++ )
	    {
	      newton_step_pos_def_matrix(i,j) = f3_inv * f3_inv * 
		newton_step_pos_def_matrix(i,j) - f3_inv * AtA(i,j);
	    }
	  newton_step_pos_def_matrix(j,j) += D_3[j];
	}

      // Compute the direction of the newton step for x
      Real rcond( -1. );
      Real r_norm;
      int info;
      if ( cg_tol < 0 )
	{

	  info = cholesky_solve( newton_step_pos_def_matrix, 
				 newton_step_rhs, dx, rcond );	      
	  if ( info > 0 )
	    {
	      if ( verbosity > 0 )
		{
		  std::cout << "BPDN_compute_central_point() returning the ";
		  std::cout << "last iterate. ";
		  std::cout << "The Hessian matrix is no longer positive ";
		  std::cout << "definite. ";
		  std::cout << "If epsilon < 1e-8 this is most likely due to ";
		  std::cout << "roundoff error. Try a larger epsilon. ";
		  std::cout << "Returning the last solution\n";
		}
	      x_new.assign( x ); u_new.assign( u );
	      newton_info = 1; // Matrix is not positive definite
	      break;
	    }
	  if ( rcond < 1.e-14 )
	    {
	      if ( verbosity > 0 )
		{
		  std::cout << "BPDN_compute_central_point() ";
		  std::cout << "Matrix ill-conditioned. ";
		  std::cout << "Returning previous solution.\n";
		}
	      x_new.assign( x );  u_new.assign( u );
	      newton_info = 2; // Jacobian matrix is ill-conditioned.
	      break;
	    }
	}
      else
	{
	  dx = 0.0;
          int iters;
	  info = conjugate_gradients_solve( newton_step_pos_def_matrix, 
					    newton_step_rhs, dx, r_norm, 
					    iters, cg_tol, N, verbosity );

	  if ( ( r_norm > 0.5 ) || ( info > 1 ) )
	    {
	      if ( verbosity > 0 )
		{
		  std::cout << "BPDN_compute_central_point() ";
		  std::cout << "Matrix ill-conditioned. ";
		  std::cout << "Returning previous solution.\n";
		}
	      x_new.assign( x );  u_new.assign( u );
	      newton_info = 1; // Matrix is not positive definite
	      break;
	    }
	}

      Adx.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, dx, 0.0 );
      // Compute the direction of the newton step for u and 
      // the minimum step size that stays in the interior
      Real step_size = 1.0;
      for ( int j = 0; j < N; j++ )
	{
	  du[j] = z_2[j] / D_1[j] - dx[j] * D_2[j] / D_1[j];
	  if ( dx[j] - du[j] > 0.0 )
	    {
	      step_size = std::min( step_size, 
				    -f_1[j] / ( dx[j] - du[j] ) );
	    }
	  if ( -dx[j] - du[j] > 0.0 )
	    {
	      step_size = std::min( step_size,
				    -f_2[j] / ( -dx[j] - du[j] ) );
	    }
	}
      Real aqe = Adx.dot( Adx ), bqe = 2.0 * r.dot( Adx ),
	cqe = r.dot( r ) -  epsilon * epsilon;
	     
      step_size = std::min( step_size, ( -bqe + std::sqrt( bqe * bqe - 4.0 * aqe * cqe ) ) / ( 2.0 * aqe ) );
      step_size *= 0.99;
	
      // Backtracking line search 
      bool sufficient_decrease( false );
      int line_search_iter( 0 );
      Real objectivep( 0.0 ), f3p( 0.0 ), lambda2( 0.0 );
      while ( !sufficient_decrease )
	{
	  objectivep = 0.0;
	  for ( int j = 0; j < N; j++ )
	    {
	      x_new[j] = x[j] + step_size * dx[j];
	      u_new[j]  = u[j] + step_size * du[j];
	      f_1_new[j] =  x_new[j] - u_new[j];
	      f_2_new[j] = -x_new[j] - u_new[j];
	      objectivep += u_new[j] - ( std::log( -f_1_new[j] ) + 
					 std::log( -f_2_new[j] ) ) / t;
	    }
	  for ( int i = 0; i < M; i ++ )
	    {
	      r_new[i] = r[i] + step_size * Adx[i];
	    }
	  f3p = 0.5 * ( r_new.dot( r_new ) - epsilon * epsilon );
	  objectivep -= std::log( -f3p ) / t;
	  lambda2 =  0.0;
	  for ( int j = 0; j < N; j++ )
	    {
	      lambda2 += -( z_1[j] * dx[j] + z_2[j] * du[j] ) / t;
	    }
	  lambda2 *= -1.0;
	  Real flin = objective - lambda2 * alpha * step_size;	    
	  step_size *= beta;
	  if ( objectivep <= flin )
	    {
	      sufficient_decrease = true;
	    }
	  line_search_iter++;
	  if ( line_search_iter > 30 )
	    {
	      if ( verbosity > 0 )
		{
		  std::cout << "BPDN_compute_central_point() ";
		  std::cout << "Line search failed. Returning the last ";
		  std::cout << "solution.\n";
		}
	      x_new.assign( x );
	      newton_info = 3; // line search failed
	      break;
	    }
	}

      for ( int j = 0; j < N; j++ )
	{
	  x[j] =  x_new[j];
	  u[j] =  u_new[j];
	  f_1[j] = f_1_new[j];
	  f_2[j] = f_2_new[j];
	}
      r = r_new;
      f3 = f3p;
      objective = objectivep;

      if ( lambda2 * 0.5 < newton_tol )
	{
	  done = true;
	}
      else
	{
	  done = false;
	}

      newton_iter++;
      if ( newton_iter >= newton_max_iter )
	{
	  done = true;
	}
      if ( verbosity > 1 )
	{
	  Real dx_norm = dx.normFrobenius();
	  dx_norm *= dx_norm;
	  Real du_norm = du.normFrobenius();
	  du_norm *= du_norm;
	  Real norm = std::sqrt( dx_norm + du_norm ); 
	  std::cout << "Newton iteration: " << newton_iter << "\n";
	  std::cout << "\tObjective: " << objective << "\n";
	  std::cout << "\tStep size: " << step_size * norm << "\n";
	  if ( cg_tol < 0 )
	    std::cout << "\tCondition number: " << rcond << "\n";
	  else
	    std::cout << "\tConjugate gradients residual: " << r_norm << "\n";
	}
    }
  return newton_info;
};

void BPDN_log_barrier_interior_point_method( RealMatrix &A, const RealVector &b, 
					     RealMatrix &result, Real epsilon, 
					     Real log_barrier_tol, Real cg_tol,
					     int verbosity )
{
  // Extract matrix shapes
  int N( A.numCols() );

  // Initialise memory for the solutions
  result.shapeUninitialized( N, 1);

  // Basis Pursuit Denoising algorithm parameters
  Real mu ( 10 );
  Real newton_tol ( log_barrier_tol );
  int newton_max_iter ( 30 );
 
  // Compute A'A. this will be used each time BPDN_compute_central_point_method
  // is called
  RealMatrix AtA( N, N );
  AtA.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, A, A, 0.0 ); 
 
  //------------------------------------------------------------//
  // Check whether starting point is in the feasiable region    //
  // If not use the least squares solution as new start point X //
  //------------------------------------------------------------//
  RealVector x( N, false ); // solution
  RealVector r(Teuchos::Copy, b.values(), b.length());       // residual
  //Not necessary because X is set to zero( uninitialized) above
  //r.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, x, -1.0 );

  if ( r.normFrobenius() > epsilon )
    {
      // Starting point is not feasiable
      if ( verbosity > 0 )
	{
	  std::cout << "Initial guess is not feasiable. ";
	  std::cout << "Computing least squares estimate\n";
	}

      // Solve AX = b
      int rank;
      RealVector singular_values; // singular values
      svd_solve( A, b, x, singular_values, rank );
      // Compute reciprocal of condition number
      Real rcond = singular_values[singular_values.length()-1] / 
	singular_values[0];
      if ( rcond < 1e-14 )
	{
	  std::string msg = "BPDN_log_barrier_interior_point_method() ";
	  msg += "A is ill conditioned. Cannot find starting point";
	  throw( std::runtime_error( msg ) );
	}
    }

  // Default: choose initial value of t so that the duality gap after 
  // the first step will be about the original norm  
  Real t = std::max( ( 2.0 * (Real)N + 1.0 ) / x.normOne(), 1.0 );

  RealVector u( N, false );
  Real x_max = x.normInf();
  for ( int j = 0; j < N; j++ )
    {
      u[j] = 0.95 * std::fabs( x[j] ) + 0.10 * x_max;
    }

  int num_log_barrier_iter = std::ceil( ( std::log( 2. * (Real)N + 1. ) - 
					  std::log( log_barrier_tol ) - 
					  std::log( t ) ) / std::log( mu));
      
  if ( verbosity > 1 )
    {
      std::cout << "\nInitial l1 norm: " << x.normOne() << "\n";
      std::cout << "Initial objective: " << sum( u ) << "\n";
      std::cout << "Number of log-barrier iterations: ";
      std::cout << num_log_barrier_iter << "\n";
    }

  // Run newton steps
  for ( int iter = 0; iter < num_log_barrier_iter; iter++ )
    {
      if ( verbosity > 1 )
	std::cout << "\nLog-barrier iteration: " << iter + 1<< "\n";
	  
      int info = BPDN_compute_central_point( A, b, x, u, AtA, epsilon, t, 
					     newton_tol, newton_max_iter, 
					     cg_tol, verbosity );

      if ( verbosity > 1 )
	{

	  std::cout.precision( std::numeric_limits<Real>::digits10 );
	  std::cout.setf( std::ios::scientific );
	  std::cout << "l1 norm: " << x.normOne() << "\n";
	  std::cout << "t: " << t << std::endl;
	}
      if ( info > 0 ) break;
      t *= mu;
    }
      
  RealVector result_view( Teuchos::View, result[0], N );
  result_view.assign( x );
};

} // namespace Surrogates

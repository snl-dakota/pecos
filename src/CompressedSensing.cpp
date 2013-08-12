#include "CompressedSensing.hpp"
#include "MathTools.hpp"

namespace Pecos {

/*
solverType solverTypeCast( int i )
{
  solverType sType;
  switch( i )
    {
    case 5:
      {
	sType = DEFAULT_REGRESSION;
	break;
      }
    case 6:
      {
	sType = DEFAULT_LEAST_SQ_REGRESSION;
	break;
      }
    case 7:
      {
	sType = SVD_LEAST_SQ_REGRESSION;
	break;
      }
    case 8:
      {
	sType = EQ_CON_LEAST_SQ_REGRESSION;
	break;
      }
    case 9:
      {
	sType = BASIS_PURSUIT;
	break;
      }
    case 10:
      {
	sType = BASIS_PURSUIT_DENOISING;
	break;
      }
    case 11:
      {
	sType = ORTHOG_MATCH_PURSUIT;
	break;
      }
    case 12:
      {
	sType = LASSO_REGRESSION;
	break;
      }
    case 13:
      {
	sType = LEAST_ANGLE_REGRESSION;
	break;
      }
    default:
      {
	std::stringstream msg;
	msg << "solverTypeCast() out of range. Cannot cast " << i;
	msg << " to enum solverType";
	throw( std::runtime_error( msg.str() ) );
      }
    };
  return sType;
};
*/

Real CompressedSensingTool::BP_surrogate_duality_gap( 
						     RealVector &primal_residual,
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

void CompressedSensingTool::BP_primal_dual_interior_point_method( RealMatrix &A, 
								  RealMatrix &B, 
								  RealMatrix &X,
							 Real primal_dual_tol, 
								  Real cg_tol, 
								  int verbosity )
{
  Teuchos::LAPACK<int, Real> la;
  
  // Extract matrix shapes
  int M( A.numRows() ), N( A.numCols() ), num_rhs( B.numCols() );


  // Initialise memory for the solutions
  X.shapeUninitialized( N, num_rhs );
 
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

  for ( int k = 0; k < num_rhs; k++ )
    {
      //------------------------------------------------------------//
      // Check whether starting point is in the feasiable region    //
      // If not use the least squares solution as new start point X //
      //------------------------------------------------------------//
      RealVector x( Teuchos::View, X[k], N );
      RealVector b( Teuchos::View, B[k], M );
      RealVector r( b );
      //Not necessary because X is set to zero( uninitialized) above
      //r.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, x, -1.0 );

      Real frob_norm_of_b = b.normFrobenius();
      Real frob_norm_of_r = r.normFrobenius();

      Real initial_guess_tol = ( cg_tol < 0 ) ? 1e-8 : cg_tol;
      if ( ( frob_norm_of_r * frob_norm_of_r ) /
	   ( frob_norm_of_b * frob_norm_of_b ) > initial_guess_tol )
	{
	  if ( verbosity > 0 )
	    {
	      PCout << "BP_primal_dual_interior_point_method() ";
	      PCout << "Initial guess is not feasiable. ";
	      PCout << "Computing least squares estimate\n";
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
	  PCout << "At Initialisation:\n";
	  PCout << "\tObjective: " << sum( u ) << "\n";
	  PCout << "\tPrimal-dual gap: " << sdg << "\n";
	  PCout << "\tPrimal residual: ";
	  PCout << r.normFrobenius() << "\n";
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
		      PCout << "BP_primal_dual_interior_point_method() ";
		      PCout << "The Hessian matrix is no longer positive ";
		      PCout << "definite. If epsilon < 1e-8 this is most ";
		      PCout << "likely due to roundoff error. ";
		      PCout << "Try a larger epsilon. ";
		      PCout << "Returning the last solution.\n";
		    }
		  x_new.assign( x );
		  break;
		}
	      if ( rcond < 1.e-14 )
		{
		  if ( verbosity > 0 )
		    {
		      PCout << "BP_primal_dual_interior_point_method() ";
		      PCout << "Matrix ill-conditioned. ";
		      PCout << "Returning the last solution.\n";
		    }
		  x_new.assign( x );
		  break;
		}
	    }
	  else
	    {
	      dv = 0.0;
	      info = conjugate_gradients_solve( newton_step_pos_def_matrix, 
						newton_step_rhs, dv, r_norm, 
						cg_tol, M, verbosity );
	      if ( ( r_norm > 0.5 ) || ( info > 1 ) )
		{
		  if ( verbosity > 0 )
		    {
		      PCout << "BP_primal_dual_interior_point_method() ";
		      PCout << "Matrix ill-conditioned. ";
		      PCout << "Returning the last solution.\n";
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
		      PCout << "BP_primal_dual_interior_point_method() ";
		      PCout << "Line search failed. Returning the last ";
		      PCout << "solution.\n";
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
	      PCout << "Newton iteration: " << primal_dual_iter << "\n";
	      PCout << "\tObjective: " << sum( u ) << "\n";
	      PCout << "\tPrimal-dual gap: " << sdg << "\n";
	      PCout << "\tPrimal residual: ";
	      PCout << r.normFrobenius() << "\n";
	      if ( cg_tol < 0 )
		PCout << "\tCondition number: " << rcond << "\n";
	      else
		PCout << "\tConjugate gradients residual: " << r_norm <<"\n";
	    }
	}

      RealVector X_col( Teuchos::View, X[k], N );
      X_col.assign( x );
    }
};

int CompressedSensingTool::BPDN_compute_central_point( RealMatrix &A, 
						       RealVector &b, 
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
  RealVector r( b );
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
		  PCout << "BPDN_compute_central_point() returning the ";
		  PCout << "last iterate. ";
		  PCout << "The Hessian matrix is no longer positive ";
		  PCout << "definite. ";
		  PCout << "If epsilon < 1e-8 this is most likely due to ";
		  PCout << "roundoff error. Try a larger epsilon. ";
		  PCout << "Returning the last solution\n";
		}
	      x_new.assign( x ); u_new.assign( u );
	      newton_info = 1; // Matrix is not positive definite
	      break;
	    }
	  if ( rcond < 1.e-14 )
	    {
	      if ( verbosity > 0 )
		{
		  PCout << "BPDN_compute_central_point() ";
		  PCout << "Matrix ill-conditioned. ";
		  PCout << "Returning previous solution.\n";
		}
	      x_new.assign( x );  u_new.assign( u );
	      newton_info = 2; // Jacobian matrix is ill-conditioned.
	      break;
	    }
	}
      else
	{
	  dx = 0.0;
	  info = conjugate_gradients_solve( newton_step_pos_def_matrix, 
					    newton_step_rhs, dx, r_norm, 
					    cg_tol, N, verbosity );

	  if ( ( r_norm > 0.5 ) || ( info > 1 ) )
	    {
	      if ( verbosity > 0 )
		{
		  PCout << "BPDN_compute_central_point() ";
		  PCout << "Matrix ill-conditioned. ";
		  PCout << "Returning previous solution.\n";
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
		  PCout << "BPDN_compute_central_point() ";
		  PCout << "Line search failed. Returning the last ";
		  PCout << "solution.\n";
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
	  PCout << "Newton iteration: " << newton_iter << "\n";
	  PCout << "\tObjective: " << objective << "\n";
	  PCout << "\tStep size: " << step_size * norm << "\n";
	  if ( cg_tol < 0 )
	    PCout << "\tCondition number: " << rcond << "\n";
	  else
	    PCout << "\tConjugate gradients residual: " << r_norm << "\n";
	}
    }
  return newton_info;
};

void CompressedSensingTool::BPDN_log_barrier_interior_point_method( RealMatrix &A, RealMatrix &B, RealMatrix &X, Real epsilon, Real log_barrier_tol, Real cg_tol, int verbosity )
{
  // Extract matrix shapes
  int M( A.numRows() ), N( A.numCols() ), num_rhs( B.numCols() );

  // Initialise memory for the solutions
  X.shapeUninitialized( N, num_rhs );

  // Basis Pursuit Denoising algorithm parameters
  Real mu ( 10 );
  Real newton_tol ( log_barrier_tol );
  int newton_max_iter ( 30 );
 
  // Compute A'A. this will be used each time BPDN_compute_central_point_method
  // is called
  RealMatrix AtA( N, N );
  AtA.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, A, A, 0.0 ); 

  // Allocate memory for solution to one rhs vector
  RealVector  u( N, false ); 
 
  for ( int k = 0; k < num_rhs; k++ )
    {
      //------------------------------------------------------------//
      // Check whether starting point is in the feasiable region    //
      // If not use the least squares solution as new start point X //
      //------------------------------------------------------------//
      RealVector x( Teuchos::View, X[k], N ); // solution
      RealVector b( Teuchos::View, B[k], M ); // rhs
      RealVector r( b );                      // residual
      //Not necessary because X is set to zero( uninitialized) above
      //r.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, x, -1.0 );

      if ( r.normFrobenius() > epsilon )
	{
	  // Starting point is not feasiable
	  if ( verbosity > 0 )
	    {
	      PCout << "Initial guess is not feasiable. ";
	      PCout << "Computing least squares estimate\n";
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
	  PCout << "\nInitial l1 norm: " << x.normOne() << "\n";
	  PCout << "Initial objective: " << sum( u ) << "\n";
	  PCout << "Number of log-barrier iterations: ";
	  PCout << num_log_barrier_iter << "\n";
	}

      // Run newton steps
      for ( int iter = 0; iter < num_log_barrier_iter; iter++ )
	{
	  if ( verbosity > 1 )
	    PCout << "\nLog-barrier iteration: " << iter + 1<< "\n";
	  
	  int info = BPDN_compute_central_point( A, b, x, u, AtA, epsilon, t, 
						 newton_tol, newton_max_iter, 
						 cg_tol, verbosity );

	  if ( verbosity > 1 )
	    {

	      PCout.precision( std::numeric_limits<Real>::digits10 );
	      PCout.setf( std::ios::scientific );
	      PCout << "l1 norm: " << x.normOne() << "\n";
	      PCout << "t: " << t << std::endl;
	    }
	  if ( info > 0 ) break;
	  t *= mu;
	}
      
      RealVector X_col( Teuchos::View, X[k], N );
      X_col.assign( x );
    }
};

void CompressedSensingTool::orthogonal_matching_pursuit( RealMatrix &A, 
							 RealVector &b, 
							 RealMatrix &solutions,
						    RealMatrix &solution_metrics,
							 Real epsilon, 
						   int max_num_non_zero_entries,
							 int verbosity )
{
  Teuchos::BLAS<int, Real> blas;

  int M( A.numRows() ), N( A.numCols() );

  // Determine the maximum number of iterations
  //if ( max_num_non_zero_entries == 0 ) max_num_non_zero_entries = M;
  int max_num_indices( std::min( M, max_num_non_zero_entries ) );
  max_num_indices = std::min( N, max_num_indices );

  // Vector to store non-zero solution entries
  RealVector x_sparse;	

  int memory_chunk_size = std::min( N, M );
  // if I use min( tmp, M ) where tmp < M then seq faults will occur.
  // Not sure why has something to do with resizing memory
  int initial_N = std::min( memory_chunk_size, N );

  // Initialise entries of all solutions to zero
  solutions.shape( N, initial_N );

  // Allocate memory to store solution metrics
  solution_metrics.shapeUninitialized( 2, initial_N );

  // Matrix to store Q and R factorization
  RealMatrix Q( M, initial_N ), R( initial_N, initial_N );

  // Compute residual
  RealVector residual( b );

  // Compute correlation of columns with residual
  RealMatrix Atb( N, 1, false );
  Atb.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 
		1.0, A, residual, 0.0 );
  RealMatrix correlation( Teuchos::Copy, Atb, N, 1 );
  
  // Compute norm of residual
  Real b_norm_sq( b.dot( b ) ), residual_norm( std::sqrt( b_norm_sq ) );

  // Matrix to store the full rank matrix associated with x_sparse
  RealMatrix A_sparse_memory( M, initial_N, false );
  RealMatrix AtA_sparse_memory( N, initial_N, false );
  RealMatrix Atb_sparse_memory( N, 1, false );
  bool residual_computed = false;

  verbosity = 2;
  if ( verbosity > 1 )
    {
      PCout << "Orthogonal Matching Pursuit" << std::endl;
      std::printf( "Iter\tAdded\tResidual\tl1 norm of x\n" );
    }

  int num_active_indices( 0 );
  IntVector active_index_set( max_num_indices );
  bool done = false;
  while ( !done )
    {
      // Find the column that has the largest inner product with
      // the residual
      // Warning IAMX returns the index of the element with the 
      // largest magnitude but IAMAX assumes indexing 1,..,N not 0,...,N-1
      int active_index = blas.IAMAX( correlation.numRows(), 
				     correlation[0], 1 ) - 1;

      // todo define active_index_set as std::set and use find function
      for ( int i = 0; i < num_active_indices; i++ )
	{
	  if ( active_index_set[i] == active_index )
	    {
	      if ( verbosity > 1 ){
		PCout << "Exiting: New active index has already been added. ";
		PCout << "This has likely occured because all correlations ";
		PCout << "are roughly the same size. ";
		PCout << "This means problem has been solved to roughly ";
		PCout << "machine precision. Check this before continuing.";
		PCout << std::endl;
		}
	      done = true;
	    }
	}
      if ( done ) break;

      if ( Q.numCols() <= num_active_indices )
	{
	  Q.reshape( Q.numRows(), Q.numCols() + memory_chunk_size );
	  R.reshape( R.numRows() + memory_chunk_size, 
		     R.numCols() + memory_chunk_size );
	  AtA_sparse_memory.reshape( N, 
				     AtA_sparse_memory.numCols() + memory_chunk_size );
	  A_sparse_memory.reshape( N, A_sparse_memory.numCols()+memory_chunk_size);
	  solutions.reshape( N, solutions.numCols() + memory_chunk_size );
	  solution_metrics.reshape( 2, solution_metrics.numCols() + memory_chunk_size );
	}
      // Update the QR factorisation.	 
      RealMatrix A_col( Teuchos::View, A, M, 1, 0, active_index );
      int colinear = qr_factorization_update_insert_column( Q, R, A_col, 
							    num_active_indices );

      if ( !colinear )
	{
	  active_index_set[num_active_indices] = active_index;

	  RealMatrix AtA_sparse( Teuchos::View, AtA_sparse_memory, N, 
				 num_active_indices + 1, 0, 0 ),
	    Atb_sparse( Teuchos::View, Atb_sparse_memory, 
			num_active_indices + 1, 1, 0, 0 );
	  
	  int index( active_index_set[num_active_indices] );
	  Atb_sparse(num_active_indices,0) = Atb(index,0);
	  for ( int n = 0; n < N; n++ )
	    {
	      //AtA_sparse(n,num_active_indices) = AtA(n,index);
	      RealVector a_n( Teuchos::View, A[n], A.numRows() ),
		a_j( Teuchos::View, A[index], A.numRows() );
	      AtA_sparse(n,num_active_indices) = a_n.dot( a_j );
	    }

	  /*
	  // Create subset representations of Atb and AtA
	  RealMatrix Atb_sparse( num_active_indices + 1, 1, false ), 
	    AtA_sparse( N, num_active_indices + 1, false );
	  for ( int i = 0; i < num_active_indices + 1; i++ )
	    {
	      int index( active_index_set[i] );
	      Atb_sparse(i,0) = Atb(index,0);
	      for ( int n = 0; n < N; n++ )
		{
		  AtA_sparse(n,i) = AtA(n,index);
		}
		}*/

	  //Solve R'z = A'b via back substitution
	  RealMatrix z;
	  RealMatrix R_new( Teuchos::View, R, num_active_indices+1, 
			    num_active_indices+1, 0, 0 );
	  
	  substitution_solve( R_new, Atb_sparse, z, Teuchos::TRANS );

	  //Solve Rx = z via back substitution to obtain signal
	  substitution_solve( R_new, z, x_sparse );

	  correlation.assign( Atb );
	  correlation.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
				-1.0, AtA_sparse, x_sparse, 1.0 );
	  Real z_norm = z.normFrobenius();
	  //residual_norm = std::sqrt( b_norm_sq - z_norm * z_norm );
	  // the above suffers from numerical precion problems so
	  // compute residual exactly. This is more expensive though
	  // so only do when necessary
	  if ( b_norm_sq - z_norm * z_norm  < 1e-8 )
	    {

	      // Create subset representations of A
	      RealMatrix A_sparse( Teuchos::View, A_sparse_memory, 
				   M, num_active_indices+1, 0, 0 );
	      for ( int m = 0; m < M; m++ )
		A_sparse(m,num_active_indices) = A(m,active_index);

	      if ( !residual_computed )
		{
		  // This is the first time A_sparse has been needed, so 
		  // update A_sparse_memory will all missing information
		  for ( int i = 0; i < num_active_indices; i++ )
		    {
		      int index( active_index_set[i] );
		      for ( int m = 0; m < M; m++ )
			A_sparse(m,i) = A(m,index);
		    }
		  residual_computed = true;
		}
	      
	      /*RealMatrix A_sparse( M, num_active_indices + 1, false );
	      for ( int i = 0; i < num_active_indices + 1; i++ )
		{
		  int index( active_index_set[i] );
		  for ( int m = 0; m < M; m++ )
		    {
		      A_sparse(m,i) = A(m,index);
		    }
		    }*/
	      RealVector residual( b );
	      residual.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
				 -1.0, A_sparse, x_sparse, 1.0 );
	      residual_norm = residual.normFrobenius();
	    }
	  else
	    {
	      residual_norm = std::sqrt( b_norm_sq - z_norm * z_norm );
	    }
	    
	  num_active_indices++;   
	}
      else
	{
	  //New column was co linear so ignore
	  if ( verbosity > 0 )
	    {
	      std::stringstream msg;
	      msg << "No variable added. Column " << active_index;
	      msg << " was colinear " << std::endl;
	      PCout << msg.str();
	    }
	}
      for ( int n = 0; n < num_active_indices; n++ )
	{
	  solutions(active_index_set[n],num_active_indices-1) = x_sparse[n];
	}
      solution_metrics(0,num_active_indices-1) = residual_norm;
      solution_metrics(1,num_active_indices-1) = active_index;     

      if ( verbosity > 1 )
	std::printf( "%d\t%d\t%1.5e\t%1.5e\n", num_active_indices, 
		     active_index, residual_norm, 
		     x_sparse.normOne() );
 
      if ( residual_norm <= epsilon )
	{
	  if ( verbosity > 1 )
	    PCout << "Exiting: residual norm lower than tolerance" << std::endl;
	  done = true;
	}
      
      if ( num_active_indices >= max_num_indices )
	{
	  if ( verbosity > 1 )
	     PCout<<"Exiting: maximum number of covariates reached"<<std::endl;;
	  done = true;
	}

      if ( colinear )
	{
	  // Usually occurs when A is a vandermonde matrix and the inputs 
	  // have been standardized. The first column will be colinear
	  // This condition should only occur when 
	  // num_covariates = max_num_covariates - 1
	   if ( verbosity > 1 )
	     PCout << "Exiting: attempted to add colinear vector" << std::endl;
	  done = true;
	}
    }
  // remove unused memory
  solutions.reshape( N, num_active_indices );
  solution_metrics.reshape( 2, num_active_indices );
}

// Need to add check for linear dependence to both lars and fast omp
// LARS only adds one variable at a time. At the moment an error 
// is thrown if two or more variables are selected. I need
// to implement a choice if this happens
// Check that when num_covariates == N that the least squares solution is returned
// return solution at each iteration
void CompressedSensingTool::least_angle_regression( RealMatrix &A, 
						    RealVector &b, 
						    RealMatrix &solutions,
						    RealMatrix &solution_metrics,
						    Real epsilon, 
						    //solverType solver,
						    short solver,
						    Real delta,
						    int max_num_iterations,
						    int verbosity )
{
  Teuchos::BLAS<int, Real> blas;

  int M( A.numRows() ), N( A.numCols() );
  
  int max_num_covariates = std::min( M, N );
  if ( delta > std::numeric_limits<Real>::epsilon() )
    max_num_covariates = std::min( N, max_num_iterations );
  
  int max_num_iter = std::min( 2*M, max_num_iterations );

  // Lasso will usually require more iterations than max covariates
  // as variables are added and deleted during the algorithm. However
  // to allow cross validation I restrict the number of iterations 
  // max_num_iter to max_num_covariates. If I do not do this then
  // LASSO can use different iteration numbers for different training
  // data but the same cross validation model options
  //if ( solver == LASSO )
  // max_num_iter *= 10;

  int memory_chunk_size = std::min( M, 500 );

  // Allocate memory to store solution metrics
  PCout << max_num_iter << std::endl;
  PCout << max_num_iterations << std::endl;
  PCout << M << std::endl;
  PCout << N << std::endl;
  solution_metrics.shapeUninitialized( 2, max_num_iter );

  // Memory for active indices
  std::vector<int>::iterator vec_iter;
  std::vector<int> active_indices;
  
  // Store inactive indices
  std::set<int>::iterator inactive_index_iter;
  std::set<int> inactive_indices;
  for ( int n = 0;  n < N; n++) inactive_indices.insert( n );
  
  // Matrix to store the full rank matrix associated with x_sparse
  RealMatrix A_sparse;
  
  // Matrix to store cholesky factorization of A'A and initialize to zero
  int size = std::min( N, memory_chunk_size );
  RealMatrix U( size, size );
  //RealMatrix U( N, N );

  // Initialise all entries of x to zero
  solutions.shape( N, memory_chunk_size );
  //solutions.shape( N, max_num_iter );
  
  // Compute residual
  RealVector residual( b );

  // Compute correlation of columns with residual
  RealMatrix Atb( N, 1, false );
  Atb.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 
		1.0, A, residual, 0.0 );
  RealMatrix correlation( Teuchos::Copy, Atb, N, 1 );

  if ( verbosity > 1 )
    {
      PCout << "LASSO/LARS ( delta = " << delta << " )\n";
      std::printf( "Iter\tAdded\tDropped\t\tSparsity\tC\t\tResidual\tl1 norm of x\n" );
    }
      
  bool done = false;

  // if solver == LARS sign_condition_violated will always be true. 
  // if solver = LASSO we must alwways add a variable during the first iteration.
  // The value of sign_condition_violated will there after will be depend on
  // whether the LASSO algorithm requires a variable to be added or removed.  
  bool sign_condition_violated = false; 

  // The index of the column of A we are adding to the active set
  int index_to_add;

  // The index of the element in the active_indices array which
  // specifies how to find which column of A violates the sign condition.
  // i.e. active_indices[index_to_drop] is the column of A which
  // violates the sign condition
  int index_to_drop( -1 );

  // The lasso approximation of b. Initialize to zero
  RealMatrix b_hat( M, 1 );

  int homotopy_iter( 0 );
  while ( !done )
    {     
      int num_covariates( active_indices.size() );
      // needs to be negative because if mean of b data is zero
      // then the 0 variable will not be added because its abs_correlation
      // will be 0 and a segfault will occur.
      Real max_abs_correlation( -1.0 );
      int prev_iter( std::max( 0, homotopy_iter -1 ) );
      for ( inactive_index_iter = inactive_indices.begin(); 
	    inactive_index_iter != inactive_indices.end(); 
	    inactive_index_iter++ )
	{
	  int n = *inactive_index_iter;
	  Real correlation_n = correlation(n,0);
	  Real x_n = solutions(n,prev_iter);
	  Real abs_correlation_n = std::abs( correlation_n-delta*x_n );
	  if ( abs_correlation_n > max_abs_correlation )
	    {
	      max_abs_correlation = abs_correlation_n;
	      index_to_add = n;
	    }
	}


      if ( U.numRows() <= num_covariates )
	{
	  U.reshape( U.numRows() + memory_chunk_size, 
		     U.numCols() + memory_chunk_size );
	}
      if ( solutions.numCols() <= homotopy_iter )
	{
	  solutions.reshape( N, solutions.numCols() + memory_chunk_size );
	}

      // Add the new index ( if it exists ) to the active set and 
      // update the Cholesky factorization
      int colinear = 0;
      if ( !sign_condition_violated )
	{
	  RealMatrix A_col( Teuchos::View, A, M, 1, 0, 
			    index_to_add );
	  colinear = cholesky_factorization_update_insert_column( A_sparse, 
								  U, 
								  A_col, 
								  num_covariates,
								  delta );
	  if ( colinear == 1 )
	    {
	      if ( verbosity > 0  )
		{
		  // Usually occurs when A is a vandermonde matrix and the inputs 
		  // have been standardized. The first column will be colinear
		  // This condition should only occur when 
		  // num_covariates = max_num_covariates - 1
		  std::stringstream msg;
		  msg << "Exiting: attempted to add colinear vector\n";
		  PCout << msg.str();
		}
	      break;
	    }

	  if ( verbosity > 1 )
	    std::printf( "%d\t%d\t\t\t%d\t\t", homotopy_iter, 
			 index_to_add, (int)active_indices.size()+1  );
	  column_append( A_col, A_sparse );
	  active_indices.push_back( index_to_add );
	  inactive_index_iter = inactive_indices.find( index_to_add );
	  inactive_indices.erase( inactive_index_iter );
	  num_covariates = (int)active_indices.size();

	  // store which variable was added to the active index set
	  solution_metrics(1,homotopy_iter) = index_to_add;
	}
	
      // Get the signs of the correlations
      RealMatrix s_sparse( num_covariates, 1, false );
      for ( int active_index = 0; active_index < num_covariates; active_index++ )
	{
	  s_sparse(active_index,0) = 
	    sgn( correlation(active_indices[active_index],0) );
	}

      // normalisation_factor = 1 / sqrt( s'*inv(A'*A)*s )
      // inv(A'A)*s  = > solve A'*A*z \ s => U'U \ s 
      // so solve two problems w = U' \ s then z = U \ w
      RealMatrix w, z;	      
      RealMatrix U_old( Teuchos::View, U, num_covariates, num_covariates, 0, 0 );
      substitution_solve( U_old, s_sparse, w, Teuchos::TRANS );
      substitution_solve( U_old, w, z, Teuchos::NO_TRANS );
      Real normalisation_factor = 1.0 / std::sqrt( blas.DOT( num_covariates, 
							     s_sparse[0], 
							     1, z[0], 1 ) );


      // Compute unit vector making equal angles, less than 90^o , 
      // with the columns of A_sparse
      RealMatrix w_sparse( z );
      w_sparse *= normalisation_factor;

      // Compute the equiangular vector
      RealMatrix u_sparse( M, 1, false );
      u_sparse.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
			 1.0, A_sparse, w_sparse, 0.0 );

      // Compute angles betwen A_j and u_sparse
      RealMatrix angles( N, 1, false );
      angles.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 
		       1.0, A, u_sparse, 0.0 ); 

      // Compute the lars step size. (2.13) (Efron 2004)
      // When thea active_indices contains all covariates, gamma_hat
      // is not defined. By convention in this situation the algorithm takes 
      // gamma_ =  corelationMax / normalisation_factor , 
      // making the solution x equal to the ordinary least squares 
      // estimate for the full set of N covariates.
      Real gamma_hat = max_abs_correlation / normalisation_factor;
      if ( num_covariates < N )
	{
	  for ( inactive_index_iter = inactive_indices.begin(); 
		inactive_index_iter != inactive_indices.end(); 
		inactive_index_iter++ )
	    {
	      int n = *inactive_index_iter;
	      Real gamma_hat1 = ( max_abs_correlation - correlation(n,0) ) / 
		( normalisation_factor - angles(n,0) );
	      
	      gamma_hat = ( ( gamma_hat1 < gamma_hat ) && ( gamma_hat1 > 0 ) )
		? gamma_hat1 : gamma_hat;
	      
	      Real gamma_hat2 =  ( max_abs_correlation + correlation(n,0) ) /
 		( normalisation_factor + angles(n,0) );

	      gamma_hat = ( ( gamma_hat2 < gamma_hat ) && ( gamma_hat2 > 0 ) )
		? gamma_hat2 : gamma_hat;
	    }
	}

      // Find the first occurance of a solution element changing sign
      // (3.5) (Efron 2004)
      Real gamma_tilde = std::numeric_limits<Real>::max();
      index_to_drop = -1;
      if ( solver == LASSO_REGRESSION )
	{
	  for ( int n = 0; n < num_covariates; n++ )
	    {
	      Real gamma = -solutions(active_indices[n],prev_iter) / 
		w_sparse(n,0);
	      if ( ( gamma > 0 ) && ( gamma < gamma_tilde ) )
		{
		  index_to_drop = n;
		  gamma_tilde = gamma;
		}
	    }
	}

      Real gamma_min ( gamma_hat );
      sign_condition_violated = false;
      if ( gamma_tilde < gamma_hat )
	{
	  // The lasso sign restriction is violated. Must remove the index
	  // active_indices[index_to_drop] from the 
	  // active index set
	  gamma_min = gamma_tilde;
	  sign_condition_violated = true;
	}

      // Update the solution. Note that the currrent solution is stored in
      // the previous column of solutions.
      for ( int n = 0; n < num_covariates; n++ )
	{
	  solutions(active_indices[n],homotopy_iter) = 
	    solutions(active_indices[n],prev_iter) + gamma_min * w_sparse(n,0);
	}

      // Update the residual. 
      // b_new = b_old + gamma_min * u_sparse (2.12) (Efron 2004)
      // => r_new = b - b_new = b - ( b_old + gamma_min * u_sparse )
      //          = r_old - gamma_min * u_sparse
      for ( int m = 0; m < M; m++ )
	{
	  b_hat(m,0) += gamma_min * u_sparse(m,0);
	  residual[m] -= gamma_min * u_sparse(m,0);
	}
      
      // Update the correlation (2.15) (Efron 2004)
      for ( int n = 0; n < N; n++ )
	{
	  correlation(n,0) -= ( gamma_min * angles(n,0) );
	}

      // Compute residual and corrrelation directly.
      if (false)
	{
	  RealVector x( Teuchos::View, solutions[homotopy_iter], N );
	  residual.assign( b );
	  residual.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
			     -1.0, A, x, 1.0 );
	  correlation.multiply( Teuchos::TRANS, Teuchos::NO_TRANS,
				1.0, A, residual, 0.0 );
	}

      Real residual_norm = residual.normFrobenius();
      RealVector x( Teuchos::View, solutions[homotopy_iter], N );

      if ( sign_condition_violated )
	{
	  if ( index_to_drop < 0 )
	    {
	      std::string msg = "least_angle_regression() No ";
	      msg += "index needed to be deleted.";
	      throw( std::runtime_error( msg ) );
	    }

	  if ( verbosity > 1 )
	    {
	      std::printf( "%1.5e\t%1.5e\t%1.5e\n", max_abs_correlation, 
			   residual_norm, x.normOne() );
	      std::printf( "%d\t\t%d\t\t%d\t\t",homotopy_iter+1, 
			   active_indices[index_to_drop],
			   (int)active_indices.size()-1 );
	    }

	  cholesky_factorization_update_delete_column( U, 
						       index_to_drop, 
						       num_covariates );
	  
	  
	  // store which variable was removed from active index set
	  // minus sign indicates variable was removed
	  solution_metrics(1,homotopy_iter) =  -active_indices[index_to_drop];

	  // delete column from A_sparse and resize
	  delete_column( index_to_drop, A_sparse );
	  inactive_indices.insert( active_indices[index_to_drop] );
	  std::vector<int>::iterator it;
	  it =  active_indices.begin() + index_to_drop;
	  active_indices.erase( active_indices.begin() + index_to_drop );
	  num_covariates = (int)active_indices.size();
	}

      if ( ( verbosity > 1 ) &&( !sign_condition_violated ) )
      	std::printf( "%1.5e\t%1.5e\t%1.5e\n", max_abs_correlation, 
		     residual_norm, x.normOne() );

      if ( ( homotopy_iter > 0 ) && 
	   ( residual_norm > solution_metrics(0,homotopy_iter-1 ) ) )
	{
	   if ( verbosity > 1 )
	     PCout << "Exiting: Residues are small and started to increase due to numerical errors\n";
	  done = true;
	}
      else
	{
	  solution_metrics(0,homotopy_iter) = residual_norm;
	  homotopy_iter++;
	}
 
      if ( residual_norm <= epsilon )
	{
	  if ( verbosity > 1 )
	    PCout << "Exiting: residual norm lower than tolerance\n";
	  done = true;
	}
      
      if ( num_covariates >= max_num_covariates )
	{
	  if ( verbosity > 1 )
	    PCout << "Exiting: maximum number of covariates reached\n";
	  done = true;
	}

      if ( colinear )
	{
	  // Usually occurs when A is a vandermonde matrix and the inputs 
	  // have been standardized. The first column will be colinear
	  // This condition should only occur when 
	  // num_covariates = max_num_covariates - 1
	   if ( verbosity > 1 )
	     PCout << "Exiting: attempted to add colinear vector\n";
	  done = true;
	}
     
      if ( homotopy_iter == max_num_iter )
	{
	  if ( verbosity > 1 )
	    PCout << "Exiting: maximum number of iterations reached\n";
	  done = true;
	}

    }
  // remove unused memory
  solutions.reshape( N, homotopy_iter );
  solution_metrics.reshape( 2, homotopy_iter );

  // The current solutions are the naive elastic net esimates
  // Rescale to avoid Real shrinkage (12) (Zou 2005)
  if ( delta > 0 )
    for ( int iter = 0; iter < homotopy_iter; iter++ )
      {
	for ( int n = 0; n < N; n++ )
	  {
	    solutions(n,iter) *= ( 1.0 + delta );
	  }
      }  
};
//Check out for elastic nets
//http://www2.imm.dtu.dk/pubdb/views/publication_details.php?id=3897
//http://www.mlpack.org/doxygen.php?doc=classmlpack_1_1regression_1_1LARS.html#_details

void CompressedSensingTool::standardize_inputs( RealMatrix &A, RealMatrix &B, 
						RealMatrix &A_stand, 
						RealMatrix &B_stand,
						RealVector &A_column_norms,
						RealVector &A_column_means,
						RealVector &B_means )
{

  int M( A.numRows() ), N( A.numCols() ), num_rhs( B.numCols() );
  A_stand.shapeUninitialized( M, N );
  B_stand.shapeUninitialized( M, num_rhs );

  RealVector A_column_var;
  A_column_var.sizeUninitialized( N );

  A_column_means.sizeUninitialized( N );
  A_column_norms.sizeUninitialized( N );
  for ( int n = 0; n < N; n ++ )
    {
      RealVector A_col( Teuchos::View, A[n], M );
      A_column_means[n] = mean( M, A_col.values() );
      A_column_var[n] = variance( M, A_col.values() );
      // If column is a constant vector do not center it
      if ( A_column_var[n] < 10*std::numeric_limits<Real>::epsilon() )
	A_column_means[n] = 0.0;
      for ( int m = 0; m < M; m++ )
	{
	  A_stand(m,n) = A(m,n) - A_column_means[n];
	}
    }

  for ( int n = 0; n < N; n ++ )
    {
      RealVector A_stand_col( Teuchos::View, A_stand[n], M );
      A_column_norms[n] = A_stand_col.normFrobenius();
    }

  for ( int n = 0; n < N; n ++ )
    {
      for ( int m = 0; m < M; m++ )
	{
	  A_stand(m,n) = A_stand(m,n) / A_column_norms[n];
	}
    }

  // If the first column is a vector of ones then center the rhs (B) values
  // That is only center values if A is a vandermonde matrix built by
  // a linear model
  if ( ( A_column_var[0] < 10*std::numeric_limits<Real>::epsilon() ) )
    {
      B_means.sizeUninitialized( num_rhs );
      for ( int k = 0; k < num_rhs; k++ )
	{
	  B_means[k] = 0.0;
	  for ( int m = 0; m < M; m++ )
	    {
	      B_means[k] += B(m,k);
	    }
	  B_means[k] /= M;

	  for ( int m = 0; m < M; m++ )
	    {
	      B_stand(m,k) = B(m,k) - B_means[k];
	    } 
	}
    }
  else
    {
      B_stand = B; 
    }
};

void CompressedSensingTool::solve( RealMatrix &A, RealMatrix &B, 
				   RealMatrixList &solutions, 
				   CompressedSensingOptions &opts,
				   CompressedSensingOptionsList &opts_list )
{
  int M( A.numRows() ), N( A.numCols() ), num_rhs( B.numCols() );

  RealMatrix A_stand, B_stand;
  RealVector A_column_norms, A_column_means, B_means;
  if ( opts.standardizeInputs )
    {
      standardize_inputs( A, B, A_stand, B_stand, A_column_norms, 
			  A_column_means, B_means );
    }
  else
    {
      A_stand = A;
      B_stand = B;
    }

  //solverType solver( opts.solver );
  short solver( opts.solver );
  Real solver_tolerance( opts.solverTolerance );
  if ( ( M >= N ) &&  (solver != LASSO_REGRESSION ) &&
       ( solver != LEAST_ANGLE_REGRESSION ) &&
       ( solver != ORTHOG_MATCH_PURSUIT ) && 
       ( solver !=  EQ_CON_LEAST_SQ_REGRESSION ) )
    {
      if ( opts.verbosity > 0 )
	{
	  std::string msg = "CompressedSensingTool::solve() A ";
	  msg += " Matrix is over-determined computing least squares solution.";
	  PCout <<  msg << std::endl;
	}
      solver = SVD_LEAST_SQ_REGRESSION;
      solver_tolerance = -1.0;
    }

  solutions.resize( num_rhs );
  opts_list.resize( num_rhs );
  if ( solver == SVD_LEAST_SQ_REGRESSION || solver == BASIS_PURSUIT ||
       solver == BASIS_PURSUIT_DENOISING )
    {
      // These methods only return one solution for each rhs
      RealMatrix single_solution;
      switch( solver )
	{
	case SVD_LEAST_SQ_REGRESSION:
	  {
	    PCout << "Using SVD least squares regression" << std::endl;
	    RealVector singular_values;
	    int rank(0);
	    svd_solve( A_stand, B_stand, single_solution, singular_values, rank, 
		      solver_tolerance );
#ifdef DEBUG
	    // For SVD, condition number can be extracted from singular values
	    PCout << "\nCondition number for LLS using LAPACK SVD (GELSS) is "
		  << singular_values[0]/singular_values[N-1] << '\n';
#endif // DEBUG
	    break;
	  }
	case BASIS_PURSUIT:
	  {
	    PCout << "Using basis pursuit" << std::endl;
	    BP_primal_dual_interior_point_method( A_stand, B_stand, 
						  single_solution, 
						  solver_tolerance, 
					  opts.conjugateGradientsTolerance,
						  opts.verbosity );
	    break;
	  }
	case BASIS_PURSUIT_DENOISING:
	  {
	    PCout << "Using basis pursuit denoising" << std::endl;
	    BPDN_log_barrier_interior_point_method( A_stand, B_stand, 
						    single_solution, 
						    opts.epsilon, 
						    solver_tolerance, 
					    opts.conjugateGradientsTolerance,
						    opts.verbosity );
	    break;
	  }
	default:
	  {
	    std::string msg = "CompressedSensingTool::solve() ";
	    msg += " incorrect solver specified";
	    throw( std::runtime_error( msg ) );
	  }
	}
      // Pack solutions into a  list of length (num_rhs) containing 
      // ( N x num_solutions=1 ) matrix
      for ( int k = 0; k < num_rhs; k++ )
	{
	  reshape( solutions[k], N, 1 );
	  RealMatrix solution( Teuchos::View, single_solution, N, 1, 0, k );
	  solutions[k].assign( solution );

	  // Store solution metrics in CompressedSensingOptions format.
	  opts_list[k].resize( 1 );
	  // First copy the original options
	  opts_list[k][0] = opts;
	  // There is no need to update the needed metrics as they will not
	  // effect the solution produced if used again.
	}
    }
  else
    {
      // These methods return multiple solutions for each rhs
      // solutions[k] will be a ( N x num_solutions ) matrix 
      // containing multiple solutions x of Ax = b
      for ( int k = 0; k < num_rhs; k++ )
	{
	  // ( M x 1 ) rhs vector of Ax = b
	  RealVector b( Teuchos::View, B_stand[k], M );
	  switch( solver )
	    {
	    case ORTHOG_MATCH_PURSUIT:
	      {
		PCout << "Using Orthogonal Matching Pursuit" << std::endl;
		RealMatrix solution_metrics;
		orthogonal_matching_pursuit( A_stand, b, solutions[k], 
					     solution_metrics,
					     opts.epsilon,
					     opts.maxNumIterations,
					     opts.verbosity );
		// Store solution metrics in CompressedSensingOptions format
		int num_solutions( solutions[k].numCols() );

		opts_list[k].resize( num_solutions );
		for ( int l = 0; l < num_solutions; l++ )
		  {
		    // First copy the original options
		    opts_list[k][l] = opts;
		    // Then update the needed metrics
		    opts_list[k][l].epsilon = solution_metrics(0,l);
		    opts_list[k][l].maxNumIterations = 
		      solution_metrics(1,l);
		    
		  }
		break;
	      }
	    case LASSO_REGRESSION:
	      {
		if ( opts.delta > 0  )
		  PCout << "Using elastic net regression" << std::endl;
		else
		  PCout << "Using LASSO regression" << std::endl;
		RealMatrix solution_metrics;
		least_angle_regression( A_stand, b, solutions[k], 
					solution_metrics, opts.epsilon,
					LASSO_REGRESSION, opts.delta,
					opts.maxNumIterations, 
					opts.verbosity );
		int num_solutions( solutions[k].numCols() );
		opts_list[k].resize( num_solutions );
		for ( int l = 0; l < num_solutions; l++ )
		  {
		    // First copy the original options
		    opts_list[k][l] = opts;
		    // Then update the needed metrics
		    opts_list[k][l].epsilon = solution_metrics(0,l);
		    opts_list[k][l].maxNumIterations = 
		      solution_metrics(1,l);
		  }
		break;
	      };
	    case LEAST_ANGLE_REGRESSION:
	      {
		PCout << "Using least angle regression" << std::endl;
		RealMatrix solution_metrics;
		least_angle_regression( A_stand, b, solutions[k], 
					solution_metrics, opts.epsilon,
					LEAST_ANGLE_REGRESSION, opts.delta,
					opts.maxNumIterations,
					opts.verbosity );
		int num_solutions( solutions[k].numCols() );
		opts_list[k].resize( num_solutions );
		for ( int l = 0; l < num_solutions; l++ )
		  {
		    // First copy the original options
		    opts_list[k][l] = opts;
		    // Then update the needed metrics
		    opts_list[k][l].epsilon = solution_metrics(0,l);
		    opts_list[k][l].maxNumIterations = 
		      solution_metrics(1,l);
		  }
		break;
	      };
	    case EQ_CON_LEAST_SQ_REGRESSION:
	      {
		PCout << "Using equality constained least squares regression\n";
		if ( opts.numFunctionSamples == 0 )
		  {
		    std::stringstream msg;
		    msg << "opts.numFunctionSamples is set to zero. ";
		    msg << "Please specify a positive value.";
		    throw( std::runtime_error( msg.str() ) );
		  }
		RealMatrix C_eq( Teuchos::View, A_stand, opts.numFunctionSamples,
				 A_stand.numCols(),  0, 0 );
		RealMatrix A_eq( Teuchos::View, A_stand, 
				 A_stand.numRows()-opts.numFunctionSamples,
				 A_stand.numCols(), opts.numFunctionSamples, 0);
		RealVector d_eq( Teuchos::View, B_stand[k], 
				 opts.numFunctionSamples );
		RealVector b_eq( Teuchos::View, 
				 B_stand[k] + opts.numFunctionSamples, 
				 B_stand.numRows()-opts.numFunctionSamples );
 		equality_constrained_least_squares_solve( A_eq, b_eq,
							  C_eq, d_eq,
							  solutions[k] );
		opts_list[k].resize( 1 );
		opts_list[k][0] = opts;
		break;
	      };
	    default:
	      {
		std::string msg = "CompressedSensingTool::solve() ";
		msg += " incorrect solver specified";
		throw( std::runtime_error( msg ) );
	      }
	    };
	  
	  if ( !opts.storeHistory )
	    {
	      // Only return final solution 
	      int num_solutions( solutions[k].numCols() );
	      RealMatrix final_solution( Teuchos::Copy, solutions[k], 
					 N, 1, 0, num_solutions - 1 );
	      solutions[k].reshape( N, 1 );
	      solutions[k] = final_solution;
	      // Only return the options that produced the final solution
	      opts_list[k][0] = opts_list[k][num_solutions-1];
	      opts_list[k].resize( 1 );
	    }
	}
    }
  
  if ( opts.standardizeInputs )
    {
      for (int k = 0; k < num_rhs; k++ )
	{
	  int num_solutions( solutions[k].numCols() );
	  for ( int l = 0; l < num_solutions; l++ )
	    {
	      // Re-adjust solution to account for normalisation of A
	      for (int j = 0; j < N; j++ )
		{
		  solutions[k](j,l) /= A_column_norms[j];
		}
	      RealVector solution_kl( Teuchos::View, solutions[k][l], N );
	      // Readjust solution to account for centering of A
	      solutions[k](0,l) -= A_column_means.dot( solution_kl );
	      // Readjust solution to account for centering of B
	      if ( B_means.length() == num_rhs )
		solutions[k](0,l) += B_means[k];
	    }
	}
    }
};

} // namespace Pecos

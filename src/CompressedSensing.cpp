/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Description:  Functions for solving undertermined linear systems using 
// compressed sensing

#include "CompressedSensing.hpp"

namespace Pecos {


  void CompressedSensing::write_matrix_to_file( RealMatrix& M, 
						std::string filename )
  {
    std::ofstream fid;
    fid.precision(16);
    fid.open(filename.c_str());
    fid << M;
    fid.close();
  }

  void CompressedSensing:: copy_data_to_matrix( Real *A_matrix, int m, int n, 
						RealMatrix& A )
  {
    if ( ( m != A.numRows() ) || ( n != A.numCols() ) )
      A.shapeUninitialized(m,n);
    int k = 0;
    for (int j = 0; j < n; j++)
      {
	for (int i = 0; i < m; i++)
	  {
	    A(i,j) = A_matrix[k];
	    k++;
	  }
      }	   
  };

  void CompressedSensing::copy_data_from_matrix( RealMatrix& A, Real *A_matrix  )
  {
    int m( A.numRows() ), n( A.numCols() );

    int k = 0;
    for (int j = 0; j < n; j++)
      {
	for (int i = 0; i < m; i++)
	  {
	    A_matrix[k] = A(i,j);
	    k++;
	  }
      }	   
  };

  Real CompressedSensing::DotProduct ( RealVector &v1, RealVector& v2 )
  {
    if ( v1.length() != v2.length() )
      PCout << "DotProduct error: vectors do not have the same length\n";
    Real sum = 0.0;
    for ( int i = 0; i < v1.length(); i++ )
      {
	sum += v1[i] * v2[i];
      }
    return sum;
  };

  int CompressedSensing::CholeskySolve( RealMatrix& A, RealMatrix& B,
					RealMatrix& X )
  {
    bool lapack_err = false;
    Teuchos::LAPACK<int, Real> la;

    int m( A.numRows() ), num_rhs( B.numCols() );
    Real *A_matrix, *B_matrix;
    A_matrix = new Real [m*m];
    B_matrix = new Real [m*num_rhs];
    copy_data_from_matrix( A, A_matrix );
    copy_data_from_matrix( B, B_matrix );
    int info;
    char uplo( 'L' );
    // Compute the Cholesky lower triangular factorisation of A
    la.POTRF( uplo, m, A_matrix, m, &info );
    if ( info > 0 ) 
      {
	PCout << "The matrix A is not positive definite"<< std::endl;
	lapack_err = true;
	return lapack_err;
      }
    if ( info < 0 ) 
      {
	lapack_err = true;
	return lapack_err;
      }

    // Solves the system of linear equations A*X = B with a symmetric
    // positive definite matrix A=LL' using the Cholesky factorization
    la.POTRS( uplo, m, num_rhs, A_matrix, m, B_matrix, m, &info );
    copy_data_to_matrix( B_matrix, m, num_rhs, X );

    delete [] A_matrix;
    delete [] B_matrix;
  
    if ( info < 0 ) lapack_err = true;

    return lapack_err;
  };

  double CompressedSensing::GetSurrogateDualityGapAndSlackness( 
					     RealMatrix& primal_residual,
					     RealMatrix& fu1, 
					     RealMatrix& fu2, 
					     RealMatrix& lamu1, 
					     RealMatrix& lamu2, 
					     RealMatrix& AtV,
					     double mu,
					     double tol,
					     RealVector& sdg,
					     RealVector& tau,
					     RealVector& slackness_norm,
					     IntVector& complete)
  {
    int m = primal_residual.numRows();
    int n = fu1.numRows();
    int num_rhs = fu1.numCols();

    double sdg_max = 0.0;
    for ( int k = 0; k < num_rhs; k++ )
      {
	RealVector fu1_column_k( Teuchos::View, fu1[k], n );
	RealVector fu2_column_k( Teuchos::View, fu2[k], n );
	RealVector lamu1_column_k( Teuchos::View, lamu1[k], n );
	RealVector lamu2_column_k( Teuchos::View, lamu2[k], n );
	sdg[k] = - DotProduct( fu1_column_k , lamu1_column_k ) -	\
	  DotProduct( fu2_column_k , lamu2_column_k );
	sdg_max = std::max( sdg_max, sdg[k] );
	tau[k] = mu * 2 * n / sdg[k];
       
	slackness_norm[k] = 0.0;
	for ( int j = 0; j < n; j++ )
	  {
	    double centrality_residual_jk = -lamu1(j,k) * fu1(j,k) - \
	      1.0 / tau[k];
	    double centrality_residual_njk = -lamu2(j,k) * fu2(j,k) - \
	      1.0 / tau[k];
	    double dual_residual_jk = lamu1(j,k) - lamu2(j,k) + AtV(j,k);
	    double dual_residual_njk = 1.0 - lamu1(j,k) - lamu2(j,k);
	    slackness_norm[k] += centrality_residual_jk * \
	      centrality_residual_jk;
	    slackness_norm[k] += centrality_residual_njk *	\
	      centrality_residual_njk;
	    slackness_norm[k] += dual_residual_jk * \
	      dual_residual_jk;
	    slackness_norm[k] += dual_residual_njk * \
	      dual_residual_njk;
	  }
	for ( int i = 0; i < m; i++ )
	  {
	    slackness_norm[k] += primal_residual(i,k) * primal_residual(i,k);
	  }
	slackness_norm[k] = \
	  std::sqrt(slackness_norm[k]);
	//PCout << sdg[k] << "," << tol << "," << complete[k] << "\n";
	complete[k] = ( ( complete[k] ) || ( sdg[k] < tol ) );
      }
    return sdg_max;
  };

  void CompressedSensing::BasisPursuit( RealMatrix& A, RealMatrix& B, 
					RealMatrix& X )
  {
    Teuchos::LAPACK<int, Real> la;

    // Basis Pursuit algorithm parameters
    Real tol( 1.e-3 );   // primal-duality gap toleranace
    int max_iter( 50 );  // maximum number of primal dual iterations
    Real alpha( 1.e-2 ); // used in residual termination condition
    Real beta( 0.5 );    // used in back tracking procedure
    Real mu( 10.0 );     // used to update slackness condition
    Real cgtol( 1.e-8 ); // tolerance of the conjugate gradient algorithm

    // Extract matrix shapes
    int m( A.numRows() ), n( A.numCols() ), num_rhs( B.numCols() );

    //-----------------------------------//
    // Estimate solution vector X = A'*B //
    //-----------------------------------//
    X.reshape( n, num_rhs );
    X.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, A, B, 0.0 );

    //------------------------------------------------------------//
    // Check whether starting point is in the feasiable region    //
    // If not use the least squares solution as new start point X //
    //------------------------------------------------------------//
    RealMatrix primal_residual( B );
    primal_residual.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 
			      1.0, A, X, -1.0 );
  
    bool compute_least_squares = false;
    for ( int k = 0; k < num_rhs; k++ )
      {
	RealVector b_column( Teuchos::View, B[k], n );
	RealVector primal_residual_column_k( Teuchos::View, 
					     primal_residual[k], n );
	Real frob_norm_of_b_k = b_column.normFrobenius();
	Real frob_norm_of_primal_residual_k = \
	  primal_residual_column_k.normFrobenius();
	if ( ( frob_norm_of_primal_residual_k * 
	       frob_norm_of_primal_residual_k ) /
	     ( frob_norm_of_b_k * frob_norm_of_b_k ) > cgtol )
	  {
	    // Starting point is not feasiable for at least one of the 
	    // columns of b.
	    compute_least_squares = true;
	    break;
	  }
      }
  
    if ( compute_least_squares )
      {
	// Solve (A'A)X = A'b where A*A' is positive definite symmetric matrix;
	// Speicifically calcualte X =  A'inv(AA')*b

	// TODO: replace AAt.multiply by SYRK blas routine which is currently
	// not in teuchos
	RealMatrix AAt( m, m ), W( n, num_rhs );
	AAt.multiply( Teuchos::NO_TRANS, Teuchos::TRANS, 1.0, A, A, 0.0 );
	CholeskySolve( AAt, B, W );
	X.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, A, W, 0.0 );
      }

    //-------------------------------------------------------------------//
    // Convert Basis Pursuit problem into traditional linear programming //
    // problem.                                                          //
    // Ensure the constraints $f_{u_1,i}$ and $f_{u_2,i}$ are always     //
    // negative by choosing $u_i=0.95x_i+0.10\max_i{x}$, $i=1,\ldots,N$  //
    //-------------------------------------------------------------------//

    // TODO:: Perhaps write linear programming function seperate
    // from basis pursuit function. merge fu1 and fu2 and lamu1 and lamu2
    // into fu and lamu respectively which have shape ( 2*n , num_rhs )

    RealMatrix U( n, num_rhs ), fu1( n, num_rhs ), fu2( n, num_rhs );
    RealMatrix lamu1( n, num_rhs ), lamu2( n, num_rhs );
    RealMatrix lamudiff( n, num_rhs );
    for ( int k = 0; k < num_rhs ; k++ )
      {
	RealVector x_column( Teuchos::View, X[k], n );
	Real x_column_k_max = x_column.normInf();
	for ( int j = 0; j < n; j++ )
	  {
	    U(j,k) = 0.95 * fabs( X(j,k) ) + 0.10 * x_column_k_max;
	    fu1(j,k) = X(j,k) - U(j,k);
	    fu2(j,k) = -X(j,k) - U(j,k);
	    lamu1(j,k) = -1.0 / fu1(j,k);
	    lamu2(j,k) = -1.0 / fu2(j,k); 
	    lamudiff(j,k) = lamu1(j,k) - lamu2(j,k);
	  }
      }

    // Compute the dual variable v
    RealMatrix V( m, num_rhs ), AtV( n, num_rhs );
    V.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, -1.0 , A, lamudiff, 0.0 );
    AtV.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, A, V, 0.0 );
    primal_residual = B;

    // Compute the primal_residual $r_\mathrm{pri} = AX-B$;
    primal_residual.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
			     1.0, A, X, -1.0);
  
    // Compute surrogate duality gap (sdg) and update slackness condition (tau)
    RealVector sdg( num_rhs ), tau( num_rhs ), slackness_norm( num_rhs );
    IntVector complete( num_rhs );
    Real sdg_max;
    sdg_max = GetSurrogateDualityGapAndSlackness( primal_residual,
						  fu1, 
						  fu2, 
						  lamu1, 
						  lamu2, 
						  AtV,
						  mu,
						  tol,
						  sdg,
						  tau,
						  slackness_norm,
						  complete);

    //PCout << "X_matrix:\n"; write_data(PCout, X, false, true, true);
    //PCout << "V_matrix:\n"; write_data(PCout, V, false, true, true);
    //PCout << "AtV_matrix:\n"; write_data(PCout, AtV, false, true, true);
    //PCout << "pimal_residual_matrix:\n"; write_data(PCout, primal_residual,false, true, true);
    //PCout << "sdg_vector:\n"; write_data(PCout, sdg, false, true, true);
    //PCout << "tau_vector:\n"; write_data(PCout, tau, false, true, true);
    //PCout << "slackness_vector:\n"; write_data(PCout, slackness_norm, false,true, true);

    //------------------------------------------//
    // Iterate though the primal-dual algorithm //
    //------------------------------------------//
    int primal_dual_iter( 0 );
    bool done = ( ( sdg_max < tol ) || ( primal_dual_iter >= max_iter ) );
    RealMatrix w1( n, num_rhs ), w2( n, num_rhs );
    RealMatrix sig1( n, num_rhs ), sig2( n, num_rhs ), sigx( n, num_rhs );
    RealMatrix  w3( m, num_rhs ), w1p( m, num_rhs ), tmp1( n, num_rhs ), \
      tmp2( n, m );
    for ( int k = 0; k < num_rhs; k++)
      complete[k] = false;
    while ( !done )
      {
	primal_dual_iter++;
      
	//-------------------------------------------------------//
	// Set up linear system to compute newton step direction //
	//-------------------------------------------------------//
	for ( int k = 0; k < num_rhs ; k++ )
	  {
	    for ( int j = 0; j < n; j++ )
	      {
		w1(j,k) = - ( - 1.0 / fu1(j,k) + 1.0 / fu2(j,k) ) / tau[k] \
		  - AtV(j,k);
		w2(j,k) = - 1.0 - ( 1.0 / fu1(j,k) + 1.0 / fu2(j,k) ) / tau[k];
		sig1(j,k) = -lamu1(j,k) / fu1(j,k) - lamu2(j,k) / fu2(j,k);
		sig2(j,k) =  lamu1(j,k) / fu1(j,k) - lamu2(j,k) / fu2(j,k);
		sigx(j,k) =  sig1(j,k) - sig2(j,k) * sig2(j,k) / sig1(j,k);
		tmp1(j,k)  =  w1(j,k) / sigx(j,k) - \
		  w2(j,k) * sig2(j,k) / ( sigx(j,k) * sig1(j,k) );
	      }
	    for ( int i = 0; i < m; i++ )
	      {
		w3(i,k) = primal_residual(i,k);
		w1p(i,k)  = w3(i,k);
	      }
	  }
	w1p.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
		      1.0, A, tmp1, -1.0 );
	//PCout << "\n Iteration: " << primal_dual_iter << "\n";
	//PCout << "w1_matrix:\n"; write_data(PCout, w1, false, true, true);
	//PCout << "w2_matrix:\n"; write_data(PCout, w2, false, true, true);
	//PCout << "w3_matrix:\n"; write_data(PCout, w3, false, true, true);
	//PCout << "sig1_matrix:\n"; write_data(PCout, sig1, false, true, true);
	//PCout << "sig2_matrix:\n"; write_data(PCout, sig2, false, true, true);
	//PCout << "sigx_matrix:\n"; write_data(PCout, sigx, false, true, true);
	//PCout << "w1p_matrix:\n"; write_data(PCout, w1p, false, true, true);

	//--------------------------------------------------------//
	// Compute newton step direction by solving linear system //
	//--------------------------------------------------------//
	RealMatrix H11p( m, m );
	RealMatrix dV( m, num_rhs );

	// Todo I neet to access elements of At. I do this by accessing
	// underlying data of A and utilising flat 1D structure. Can this 
	// be done a different way
	Real *A_matrix;
	A_matrix = A.values();
	for ( int k = 0; k < num_rhs ; k++ )
	  {
	    if ( !complete[k] )
	      {
		for ( int i = 0; i < m; i++ )
		  {
		    for ( int j = 0; j < n; j++ )
		      {
			tmp2(j,i) = 1.0 / sigx(j,k) * A_matrix[j*m+i];
		      }
		  }
		H11p.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
			       1.0, A, tmp2, 0.0 );
		//PCout << "tmp2_matrix:\n"; write_data(PCout, tmp2, false, true, true);
		//PCout << "H11p_matrix:\n"; write_data(PCout, H11p, false, true,true);
		
		RealVector w1p_column_k( Teuchos::View, w1p[k], m );
		RealVector dv_column_k( m, num_rhs );
		CholeskySolve( H11p, w1p_column_k, dv_column_k );
		// Copy dv_column_k to dV
		for ( int i = 0 ; i < m; i++ )
		  {
		    dV(i,k) = dv_column_k[i];
		  }
	      }
	  }
	//--------------------------//
	// Compute newton step size //
	//--------------------------//
	RealMatrix AtdV( n, num_rhs );
	RealMatrix dX( n, num_rhs ), dU( n, num_rhs ), dlamu1( n, num_rhs ),\
	  dlamu2( n, num_rhs );
	RealVector step_size( num_rhs );
	AtdV.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 
		       1.0, A, dV, 0.0);

	for ( int k = 0; k < num_rhs ; k++ )
	  {
	    if ( !complete[k] )
	      {
		step_size[k] = 1.0;
		for ( int j = 0; j < n; j++ )
		  {
		    dX(j,k) = ( w1(j,k) - w2(j,k) * sig2(j,k) / sig1(j,k) - \
				AtdV(j,k) ) / sigx(j,k);
		    dU(j,k) = ( w2(j,k) - sig2(j,k) * dX(j,k) ) / sig1(j,k);
		    dlamu1(j,k) = lamu1(j,k) / fu1(j,k) * ( -dX(j,k) + dU(j,k) )\
		      - lamu1(j,k) - 1.0 / ( tau[k] * fu1(j,k) );
		    dlamu2(j,k) = lamu2(j,k) / fu2(j,k) * ( dX(j,k) + dU(j,k) )\
		      - lamu2(j,k) - 1.0 / ( tau[k] * fu2(j,k) );

		    if ( dlamu1(j,k) < 0 )
		      {
			step_size[k] = std::min( step_size[k], 
						 -lamu1(j,k) / dlamu1(j,k) );
		      }
		    if ( dlamu2(j,k) < 0 )
		      {
			step_size[k] = std::min( step_size[k], 
						 -lamu2(j,k) / dlamu2(j,k) );
		      }
		    if (  ( dX(j,k) - dU(j,k) ) > 0 )
		      {
			step_size[k] = std::min( step_size[k], -fu1(j,k) / \
						 ( dX(j,k) - dU(j,k) ) );
		      }
		    if (  ( -dX(j,k) - dU(j,k) ) > 0 )
		      {
			step_size[k] = std::min( step_size[k], -fu2(j,k) / \
						 ( -dX(j,k) - dU(j,k) ) );
		      }
		  }
		step_size[k] *= 0.99;
	      }
	  }
	RealMatrix AdX( m, num_rhs );
	AdX.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
		      1.0, A, dX, 0.0);

	//PCout << "dV_matrix:\n"; write_data(PCout, dV, false, true, true);
	//PCout << "dX_matrix:\n"; write_data(PCout, dX, false, true, true);
	//PCout << "AdX_matrix:\n"; write_data(PCout, AdX, false, true, true);
	//PCout << "AtdV_matrix:\n"; write_data(PCout, AtdV, false, true, true);
	//PCout << "step_size_matrix:\n"; write_data(PCout, step_size, false,true, true);

	//-------------------------------------------------------------------//
	// Conduct line search to ensure the norm of the residuals decreases //
	// sufficiently                                                      //
	//-------------------------------------------------------------------//
	bool sufficient_decrease = false;
	int backiter = 0; // the iteration count of the backward line seach

	// Continue while the norm of the residuals has not decreased 
	// sufficiently for any of the quantities of interest.
	RealMatrix X_p( n, num_rhs ), U_p( n, num_rhs );
	RealMatrix lamu1_p( n, num_rhs ), lamu2_p( n, num_rhs );
	RealMatrix fu1_p( n, num_rhs ), fu2_p( n, num_rhs );
	RealMatrix AtV_p( n, num_rhs );
	RealMatrix V_p( m, num_rhs ), primal_residual_p( m, num_rhs );
	RealVector slackness_p_norm( num_rhs );
	while ( !sufficient_decrease )
	  {
	    sufficient_decrease = true;
	    for ( int k = 0; k < num_rhs ; k++ )
	      {
		if ( !complete[k] )
		  {
		    slackness_p_norm[k] = 0.0;
		    for ( int j = 0; j < n; j++ )
		      {
			X_p(j,k)     =  X(j,k) + step_size[k] * dX(j,k);
			U_p(j,k)     =  U(j,k) + step_size[k] * dU(j,k);
			AtV_p(j,k)   =  AtV(j,k) + step_size[k] * AtdV(j,k);
			lamu1_p(j,k) =  lamu1(j,k) + step_size[k] * dlamu1(j,k);
			lamu2_p(j,k) =  lamu2(j,k) + step_size[k] * dlamu2(j,k);
			fu1_p(j,k)   =  X_p(j,k) - U_p(j,k);
			fu2_p(j,k)   = -X_p(j,k) - U_p(j,k);

			// 0:N-1
			slackness_p_norm[k] += \
			  ( lamu1_p(j,k) - lamu2_p(j,k) + AtV_p(j,k) ) * \
			  ( lamu1_p(j,k) - lamu2_p(j,k) + AtV_p(j,k) );
			// N:2N-1
			slackness_p_norm[k] += \
			  ( 1.0 - lamu1_p(j,k) - lamu2_p(j,k) ) * \
			  ( 1.0 - lamu1_p(j,k) - lamu2_p(j,k) );
			// 2N:3N-1
			slackness_p_norm[k] += \
			  ( -lamu1_p(j,k) * fu1_p(j,k) - 1.0 / tau[k] ) * \
			  ( -lamu1_p(j,k) * fu1_p(j,k) - 1.0 / tau[k] );
			// 3N:4N-1
			slackness_p_norm[k] += \
			  ( -lamu2_p(j,k) * fu2_p(j,k) - 1.0 / tau[k] ) * \
			  ( -lamu2_p(j,k) * fu2_p(j,k) - 1.0 / tau[k] );
		      }

		    for ( int i = 0; i < m; i++ )
		      {
			V_p(i,k)  = V(i,k) + step_size[k] * dV(i,k);
			primal_residual_p(i,k) = primal_residual(i,k) +	\
			  step_size[k] * AdX(i,k);
			slackness_p_norm[k] += primal_residual_p(i,k) *	\
			  primal_residual_p(i,k);
		      }

		    slackness_p_norm[k] =		\
		      sqrt( slackness_p_norm[k] );
		    if ( slackness_p_norm[k] >	\
			 ( 1.0 - alpha * step_size[k] ) * slackness_norm[k] )
		      {
			sufficient_decrease = false;
		      }
		    step_size[k] *= beta;
		  }
		else
		  {
		    for ( int j = 0; j < n; j++ )
		      {
			X_p(j,k) = X(j,k);
		      }
		  }
	      }
	  
	    backiter++;
	    if ( backiter > 32 )
	      {
		PCout << "Stuck backtracking, returning last iterate.\n";
		X_p = X;
		break;
	      }
	  }
	X = X_p; U = U_p;
	V = V_p; AtV = AtV_p;
	lamu1 = lamu1_p; lamu2 = lamu2_p;
	fu1 = fu1_p; fu2 = fu2_p;
	primal_residual = primal_residual_p;
	// Update surrogate duality gap and the slackness condition
	sdg_max = GetSurrogateDualityGapAndSlackness( primal_residual,
						      fu1, 
						      fu2, 
						      lamu1, 
						      lamu2, 
						      AtV,
						      mu,
						      tol,
						      sdg,
						      tau,
						      slackness_norm,
						      complete);

	//PCout << "sdg_vector:\n"; write_data(PCout, sdg, false, true, true);
	//PCout << "tau_vector:\n"; write_data(PCout, tau, false, true, true);
	//PCout << "slackness_vector:\n"; write_data(PCout, slackness_norm, false,true, true);
	
	done = ( ( sdg_max < tol ) || ( primal_dual_iter >= max_iter ) );      
      }
    //PCout << "X_matrix:\n"; write_data(PCout, X, false, true, true);
};

  int CompressedSensing::ArgMaxMagnitude(RealVector& v)
  {
    //The following blas call does not work
    //Teuchos::BLAS<int, double> blas;
    //Real *v_raw;
    //v_raw = v.values();
    //return blas.IAMAX( v.length(), v_raw, 1 );
    int argMax = 0;
    double max = fabs(v[0]);
    for ( int i = 0; i < v.length(); i++ )
      {
	//if ( fabs(v[i]) - max > std::numeric_limits<float>::epsilon() )
	if ( fabs(v[i]) > max )
	  {
	    max = fabs(v[i]);
	    argMax = i;
	  }
      }
    return argMax;
  };

  void CompressedSensing::Product(RealMatrix& A, RealVector& x, RealVector& y,
				  Teuchos::ETransp trans, Real alpha, Real beta)
  {
    Teuchos::BLAS<int, double> blas;
    int m( A.numRows() ), n( A.numCols() );
    Real *A_raw, *x_raw, *y_raw;
    A_raw = A.values(); x_raw = x.values(); y_raw = y.values();
    blas.GEMV( trans, m, n, alpha, A_raw, m, x_raw, 1, beta, y_raw, 1);
  };

  int CompressedSensing::SVDSolve( RealMatrix& A, RealMatrix& B,
				   RealMatrix& X )
  {
    bool lapack_err = false;
    Teuchos::LAPACK<int, Real> la;

    int m( A.numRows() ), n( A.numCols() ), num_rhs( B.numCols() );
    Real *A_matrix, *b_vectors, *X_matrix;
    A_matrix = new Real [m*m];
    b_vectors = new Real [m*num_rhs];
    X_matrix = new Real [m*num_rhs];
    copy_data_from_matrix( A, A_matrix );
    copy_data_from_matrix( B, b_vectors );

    double* s_vector = new double [n]; // output: singular values  
    double rcond   = -1.; // input:  use macheps to rank singular vals of A
    int    rank    =  0;  // output: effective rank of matrix A
    int    info    =  0;  // LAPACK output flag

    // Get the optimal work array size
    int lwork = -1; // special code for workspace query
    double *work  = new double [1]; // temporary work array
    la.GELSS( m, n, num_rhs, A_matrix, m, b_vectors,
	      m, s_vector, rcond, &rank, work, lwork, &info );
    lwork = (int)work[0]; // optimal work array size returned by query
    delete [] work;
    work  = new double [lwork]; // Optimal work array

    la.GELSS( m, n, num_rhs, A_matrix, m, b_vectors,
	      m, s_vector, rcond, &rank, work, lwork, &info );

    copy_data_to_matrix( b_vectors, m, num_rhs, X );

    delete [] A_matrix;
    delete [] b_vectors;
    delete [] X_matrix;
    delete [] s_vector;
    delete [] work;
    if ( info ) lapack_err = true;

    return lapack_err;
  };
  
  // Use the modified gram-schmidt algorithm to add a new column to a
  // matrix A and ensure all columns remain orthogonal.
  // See A Wavelet Tour of Signal Processing By Stephane Mallat p 428
  void CompressedSensing::UpdateOrthogonalMatrix(RealMatrix& Q, int n, 
						 RealVector& new_column)
  {
    Teuchos::BLAS<int, double> blas;

    int m( Q.numRows() ); 
    Real *new_column_raw;
    new_column_raw = new_column.values();
    Real *Q_column_j_raw;
    for ( int j = 0; j < n; j++ )
      {
	RealVector Q_column_j( Teuchos::View, Q[j], m );
	Q_column_j_raw = Q[j];
	double scale = DotProduct( Q_column_j, new_column);
	// Subtract projection of new column onto previous column
	blas.AXPY( m, -scale, Q_column_j_raw, 1, new_column_raw, 1 ); 
      }
    // Add new orthogonalised column to Q
    double new_column_norm = new_column.normFrobenius();
    for ( int i = 0; i < m; i++ )
      {
	Q[n][i] = new_column[i] / new_column_norm;
      }
  };

  void CompressedSensing::OrthogonalMatchingPursuit( RealMatrix& A, 
						     RealMatrix& B, 
						     RealMatrix& X )
  {
    double tol = 1.e-8;
 
    int m( A.numRows() ), n( A.numCols() ), num_rhs( B.numCols() );

    // Initialise all entries of x to zero
    X.reshape( n, num_rhs );

    for ( int k = 0; k < num_rhs; k++ )
      {
 
	int max_iter( m );
      
	RealMatrix Q( m, max_iter );

	// Compute residual
	RealVector b( Teuchos::Copy, B[k], m ), residual( b ), x( n );

	// Compute correlation of columns with residual
	RealVector AtR( n );
	Product( A, residual, AtR, Teuchos::TRANS, 1.0, 0.0 );
  
	// Compute norm of residual
	Real residual_norm( DotProduct( residual, residual ) );
	Real b_norm( DotProduct( b, b ) );

	int iter( 0 );
	IntVector active_index_set( n );
	Real *AtR_vector;
	bool done = false;
	RealVector x_orthogonal( n );
	while ( !done )
	  {
	    int arg_max_AtR = ArgMaxMagnitude(AtR);
	    active_index_set[iter] = arg_max_AtR;	  
	    RealVector new_column( Teuchos::Copy, A[arg_max_AtR], m );
	    UpdateOrthogonalMatrix(Q, iter, new_column);
	    active_index_set[iter] = arg_max_AtR;
	    //PCout << "Adding variable " << arg_max_AtR << "\n";

	    Product( Q , b , x_orthogonal, Teuchos::TRANS, 1.0, 0.0 );

	    // Compute residual. Note that here we are using the
	    // orthogonalised values of x and thus this is not the true
	    // residual. The residual is not the same but 
	    // <R^mf, u_m> = <R^mf, g_m> u is orthonormal basis and g is
	    // original basis. Plot how two residuals compare.

	    // For Matrix operator= to work Teuchos::copy=true must be 
	    // set during creation of the  matrix being copied.
	    residual = b;
	    Product( Q, x_orthogonal, residual, Teuchos::NO_TRANS, -1.0, 1.0 );
	    residual_norm = DotProduct ( residual, residual );
      
	    iter++;      
	    if ( (residual_norm <= tol*b_norm) || ( iter >= max_iter-1 ) )
	      done = true;
	    else
	      {
		Product( A, residual, AtR, Teuchos::TRANS, 1.0, 0.0 );
	      }      
	  }

  
	RealMatrix A_trunc( m, iter );
	for (int j = 0; j < iter; j++ )
	  {
	    for ( int i = 0; i < m; i++ )
	      A_trunc(i,j) = A(i,active_index_set[j]);
	  }
	RealVector x_trunc( iter );
	SVDSolve( A_trunc, b, x_trunc );

	for (int j = 0; j < iter; j++ )
	  {
	    X(active_index_set[j],k) = x_trunc[j];
	  }
      }
  };

} // namespace pecos

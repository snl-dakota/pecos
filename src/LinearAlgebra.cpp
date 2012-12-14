#include "LinearAlgebra.hpp"

namespace Pecos {
int cholesky_solve( RealMatrix& A, RealMatrix& B, RealMatrix& X, Real &rcond )
{
  Teuchos::LAPACK<int, Real> la;

  int m( A.numRows() ), num_rhs( B.numCols() );
  RealMatrix A_copy( A );
  // reshape must be called before assign
  X.reshape( m, num_rhs );
  X.assign( B );
  int info;
  char uplo( 'L' );
  // Compute the Cholesky lower triangular factorisation of A
  la.POTRF( uplo, m, A_copy.values(), A_copy.stride(), &info );
  if ( info > 0 ) 
    {
      std::cout << "cholesky_solve() The matrix A is not positive definite\n";
      return info;
    }
  if ( info < 0 ) 
    {
      std::cout << "cholesky_solve() Incorrect arguments specified to POTRF()\n";
      return info;
    }

  // Compute the reciprocal of the condition number of A 
  // from its cholesky decompostion computed by POTRF. 
  // The cholesky decomposition is stored in  A_copy.
  if ( rcond < 0 )
    {
      Real *work = new Real [3*m]; // workspace array 
      int *iwork = new int [m];  // workspace array
      la.POCON( uplo, m, A_copy.values(), A_copy.stride(), A.normOne(), 
		&rcond, work, iwork, &info );
      delete [] work;
      delete [] iwork;
      if ( info < 0 ) 
	{
	  std::cout << "cholesky_solve() Incorrect arguments specified to ";
	  std::cout << "POCON()\n";
	  return info;
	}
    }

  // Solves the system of linear equations A*X = B with a symmetric
  // positive definite matrix A=LL' using the Cholesky factorization
  if ( info == 0 )
    {
      la.POTRS( uplo, m, num_rhs, A_copy.values(), A_copy.stride(), X.values(), 
		X.stride(), &info );
    }
  
  return info;
};

void qr_solve( RealMatrix &A, RealMatrix &B, RealMatrix &X, 
	       Teuchos::ETransp trans )
{
  Teuchos::LAPACK<int, Real> la;

  RealMatrix A_copy( A );
  int M( A.numRows() ), N( A.numCols() ), num_rhs( B.numCols() );
  X.reshape( N, num_rhs );
  X.assign( B );

  //---------------------------------//
  // Get the optimal work array size //
  //---------------------------------//
  
  int lwork;     // Size of Teuchos::LAPACK work array
  Real *work;  // Teuchos::LAPACK work array
  int info;      // Teuchos::LAPACK output flag 
  int lda = A_copy.stride();
  int ldb = X.stride();

  lwork = -1;             // special code for workspace query
  work  = new Real [1]; // temporary work array
  la.GELS( Teuchos::ETranspChar[trans], M, N, num_rhs, A_copy.values(), 
	   lda, X.values(), ldb, work, lwork, &info );
  // Note la.GELS does not work because line 1358 in Teuchos_Teuchos::LAPACK.hpp 
  // last argument to DGELSS_F77 uses a & which should not be there
  lwork = (int)work[0];  // optimal work array size returned by query
  delete [] work;
  work  = new Real [lwork]; // Optimal work array

  //---------------------------------//
  // Solve Ax = b                    //
  //---------------------------------//

  la.GELS( Teuchos::ETranspChar[trans], M, N, num_rhs, A_copy.values(), lda, 
	   X.values(), ldb, work, lwork, &info );
  if ( info < 0 )
    {
      std::stringstream msg;
      msg << "qr_solve() dgels failed. ";
      msg << "The " << std::abs( info ) << "-th argument had an ";
      msg << "illegal value";
      throw( std::runtime_error( msg.str() ) );
    }
  if ( info > 0 )
    {
      std::stringstream msg;
      msg << "QR Solve dgels failed. ";
      msg << "The " << info << "-th diagonal element of the ";
      msg << "triangular factor of A is zero, so that A does not have";
      msg << "full rank; the least squares solution could not be computed.";
      throw( std::runtime_error( msg.str() ) );
    }
  delete [] work;
};


void svd_solve( RealMatrix &A, RealMatrix &B, RealMatrix &X,
		RealVector &S, int &rank, Real rcond )
{
    Teuchos::LAPACK<int, Real> la;

  //-----------------//
  // Allocate memory //
  //-----------------//

  int M( A.numRows() ),  N( A.numCols() ), num_rhs( B.numCols() );
  RealMatrix A_copy( A );
  S.sizeUninitialized( std::min( M, N ) );

  //---------------------------------//
  // Get the optimal work array size //
  //---------------------------------//
  
  int lwork;     // Size of Teuchos::LAPACK work array
  Real *work;  // Teuchos::LAPACK work array
  int info;      // Teuchos::LAPACK output flag
  int lda = A_copy.stride();
  int ldb = std::max( std::max( B.stride(), lda ), N );

  X.shapeUninitialized( M, num_rhs );
  X.assign( B );
  X.reshape( ldb, num_rhs );
;

  lwork = -1;             // special code for workspace query
  work  = new Real [1]; // temporary work array
  la.GELSS( M, N, num_rhs, A_copy.values(), lda, X.values(), ldb, 
	    S.values(), rcond, &rank, work, lwork, &info );
  lwork = (int)work[0];  // optimal work array size returned by query

  delete [] work;
  work  = new Real [lwork]; // Optimal work array

  //---------------------------------//
  // Solve Ax = b                    //
  //---------------------------------//
  la.GELSS( M, N, num_rhs, A_copy.values(), lda, X.values(), ldb, 
	    S.values(), rcond, &rank, work, lwork, &info );
  X.reshape( N, num_rhs );

  delete [] work;
};

void backward_substitution_solve( RealMatrix &A, 
				  RealMatrix &B, 
				  RealMatrix &X,
				  Teuchos::ETransp trans,
				  Teuchos::EUplo uplo,
				  Teuchos::EDiag diag )
{
  int M( A.numRows() ), num_rhs( B.numCols() );

  Teuchos::LAPACK<int, Real> lapack;
  X.reshape( M, num_rhs );
  X.assign( B );

  int info;
  lapack.TRTRS( Teuchos::EUploChar[uplo], Teuchos::ETranspChar[trans], 
		Teuchos::EDiagChar[diag], 
		M, num_rhs, A.values(), A.stride(), 
		X.values(), X.stride(), &info );

  if ( info < 0 )
    {
      std::stringstream msg;
      msg << "backwardsSubstitutionSolve() dtrts failed. ";
      msg << "The " << std::abs( info ) << "-th argument had an ";
      msg << "illegal value";
      throw( std::runtime_error( msg.str() ) );
    }
  if ( info > 0 )
    {
      std::stringstream msg;
      msg << "backwardsSubstitutionSolve() dtrts failed. ";
      msg << "The " << info << "-th diagonal element of A is zero ";
      msg << "indicating that the matrix is singular and the solutions ";
      msg << "X have not been computed.";
      throw( std::runtime_error( msg.str() ) );
    }
   
}

int qr_factorization_update_insert_column( RealMatrix &Q, RealMatrix &R, 
					   RealMatrix &col, int iter )
{
  int info( 0 );
  int M( col.numRows() );
  Real col_norm = col.normFrobenius();

  if ( iter == 0 )
    {
      R(0,0) = col_norm;
      for ( int m = 0; m < M; m++ )
	Q(m,0) = col(m,0) / col_norm;
    }
  else
    {
      RealMatrix Q_old( Teuchos::View, Q, M, iter, 0, 0 );
      RealMatrix w( iter, 1 ); 
      w.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 
		  1.0, Q_old, col, 0.0 );
      Real w_norm = w.normFrobenius();

      if ( col_norm * col_norm - w_norm*w_norm <= 
	   std::numeric_limits<Real>::epsilon() )
	{
	  // New column is colinear. That is, it is in the span of the active
	  // set
	  info = 1;
	}
      else
	{
	  R(iter,iter) = std::sqrt( col_norm * col_norm - w_norm * w_norm );
	  RealMatrix R_col( Teuchos::View, R, iter, 1, 0, iter );
	  // must use assign below because operator= will not work
	  // because it will call deleteArrays when w is a copy 
	  // ( which it is here )
	  R_col.assign( w ); 
  
	  RealMatrix Qw( M, 1, false );
	  Qw.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
		       1.0, Q_old, w, 0.0 );

	  for ( int m = 0; m < M; m ++ )
	    Q(m,iter) = ( col(m,0) - Qw(m,0) ) / R(iter,iter);
	}
    }
  return info;
};

int cholesky_factorization_update_insert_column( RealMatrix &A, RealMatrix &U, 
						 RealMatrix &col, int iter,
						 Real delta )
{
  int info( 0 );
  Real col_norm = col.normFrobenius();

  if ( iter == 0 )
    {
      U(0,0) = std::sqrt( col_norm * col_norm + delta );
    }
  else
    {
      RealMatrix w( iter, 1 );
      RealMatrix U_old( Teuchos::View, U, iter, iter, 0, 0 );
      // compute column k in gramian matrix A'A
      RealMatrix gramian_col( iter, 1 );
      gramian_col.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 
			    1.0, A, col, 0.0 );
      backward_substitution_solve( U_old, gramian_col, w, Teuchos::TRANS,
				   Teuchos::UPPER_TRI );
      Real w_norm = w.normFrobenius();
      
      if ( col_norm * col_norm + delta - w_norm*w_norm <= 
	   std::numeric_limits<Real>::epsilon() )
	{
	  // New column is colinear. That is, it is in the span of the active
	  // set
	  info = 1;
	}
      else
	{
	  U(iter,iter) = std::sqrt(( col_norm*col_norm + delta )-w_norm*w_norm );
	  RealMatrix U_col( Teuchos::View, U, iter, 1, 0, iter );
	  U_col.assign( w );
	}
    }
  return info;
};

void givens_rotation( RealVector &x, RealVector &x_rot, RealMatrix &givens_matrix )
{
  givens_matrix.reshape( 2, 2 );
  x_rot.sizeUninitialized( x.length() );

  if ( x[1] == 0 )
    {
      //givens_matrix = I
      givens_matrix(0,0) = 1.0; givens_matrix(1,1) = 1.0;
      x_rot.assign( x );
    }
  else
    {
      Real x_norm = x.normFrobenius();

      givens_matrix(0,0) = x[0] / x_norm;
      givens_matrix(0,1) = x[1] / x_norm;
      givens_matrix(1,0) = -x[1] / x_norm;
      givens_matrix(1,1) = x[0] / x_norm;

      x_rot[0] = x_norm;
      x_rot[1] = 0.0;
    }
};

void cholesky_factorization_update_delete_column( RealMatrix &U, 
						  int col_index,
						  int N )
{
  if ( col_index != N - 1 )
    {
      // delete column but do not resize the matrix U
      delete_column( col_index, U, false );
    };
  
  RealVector x( 2, false );
  for ( int n = col_index; n < N-1; n++ )
    {
      RealMatrix givens_matrix;
      RealVector x_rot;
      x[0] = U(n,n); x[1] = U(n+1,n);
      givens_rotation( x, x_rot, givens_matrix );
      U(n,n) = x_rot[0]; U(n+1,n) = x_rot[1];
      if ( n < N - 2 )
	{
	  RealMatrix U_sub( Teuchos::View, U, 2, N - n - 1, n, n + 1 );
	  RealMatrix U_sub_rot( U_sub.numRows(), U_sub.numCols() );
	  U_sub_rot.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
			      1.0, givens_matrix, U_sub, 0.0 );
	  U_sub.assign( U_sub_rot );
	}
    }

  // Zero out last row and column of U
  for ( int m = 0; m < N; m++ ) U(m,N-1) = 0.0;
  for ( int n = 0; n < N; n++ ) U(N-1,n) = 0.0;
};

int conjugate_gradients_solve( RealMatrix &A, RealVector &b, RealVector &x, 
			       Real &relative_residual_norm,
			       Real cg_tol, int max_iter,
			       int verbosity )
{
  int info( 0 );

  int M( A.numRows() ), N( A.numCols() );

  if ( x.length() != N )
    // Initialize to zero
    x.size( N );

  RealVector current_x( x );

  Real b_norm = b.normFrobenius();
  RealVector residual( b );
  residual.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
		     -1.0, A, current_x, 1.0 );

  if ( b_norm < std::numeric_limits<Real>::epsilon() ) 
    b_norm = 1.0;

  relative_residual_norm = residual.normFrobenius() / b_norm ;

  if ( Teuchos::ScalarTraits<Real>::isnaninf( relative_residual_norm ) )
    {
      // At least one of the matrix inputs contains nan of inf
      if ( verbosity > 2 )
	{
	  std::stringstream msg;
	  msg << "conjugate_gradient_solve() Warning: at least one of the ";
	  msg << "matrix inputs contains nan and/or inf.\n";
	  std::cout << msg.str();
	}
      info = 3;
      return info;
    }
  
  if ( relative_residual_norm <= cg_tol )
    {
      info = 0;
      return info;
    }

  if ( verbosity > 2 )
    {
      std::cout << "CG iteration: " << 0 << ", residual: ";
      std::cout << relative_residual_norm << "\n";
    }

  RealVector p( residual );
  RealVector Ap( A.numRows() );
  
  int iter( 0 );
  bool done( false );
  Real rtr_old( residual.dot( residual ) );

  while ( !done )
    {
      Ap.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
		   1.0, A, p, 0.0 );

      Real ptAp = p.dot( Ap );

      Real alpha = rtr_old / ptAp;

      if ( Teuchos::ScalarTraits<Real>::isnaninf( alpha ) || ( ptAp <= 0 ) )
	{
	  if ( verbosity > 2 )
	    {
	      // The matrix is not positive definite
	      std::stringstream msg;
	      msg << "conjugate_gradient_solve() Warning: A is not postive ";
	      msg << "definite.\n";
	      std::cout << msg.str();
	    }
	  info = 2;
	  return info;
	}

      for ( int n = 0; n < N; n++ )
	current_x[n] += alpha * p[n];

      if ( (iter+1)%50 == 0 )
	{
	  // Avoid numerical drift due to rounding error
	  residual.assign( b );
	  residual.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
			     -1.0, A, current_x, 1.0 );
	}
      else
	{
	  for ( int m = 0; m < M; m++ )
	    residual[m] -= alpha * Ap[m];
	}
      
      Real rtr = residual.dot( residual );

      for ( int m = 0; m < M; m++ )
	p[m] = residual[m] + rtr / rtr_old * p[m];

      rtr_old = rtr;

      Real current_relative_residual_norm = std::sqrt( rtr_old ) / b_norm;

      if ( current_relative_residual_norm < relative_residual_norm )
	{
	  relative_residual_norm = current_relative_residual_norm;
	  for ( int n = 0; n < N; n++ )
	    x[n] = current_x[n];
	}

      if ( verbosity > 2 )
	{
	  std::cout << "CG iteration: " << iter + 1<< ", residual: ";
	  std::cout <<  current_relative_residual_norm << "\n";
	}

      if ( current_relative_residual_norm < cg_tol )
	{
	  if ( verbosity > 2 )
	    {
	      std::stringstream msg;
	      msg << "conjugate_gradient_solve() Exiting residual below ";
	      msg << "tolerance.\n";
	      std::cout << msg.str();
	    }
	  done = true;
	  info = 0;
	}

      iter++;
      if ( iter == max_iter )
	{
	  if ( verbosity > 2 )
	    {
	      std::stringstream msg;
	      msg << "conjugate_gradient_solve() Exiting maximum number of ";
	      msg << "iterations reached.\n";
	      std::cout << msg.str();
	    }
	  done = true;
	  info = 1;
	}
    }    
  return info;
};

void equality_constrained_least_squares_solve( RealMatrix &A, 
					       RealVector &b,
					       RealMatrix &C, 
					       RealVector &d,
					       RealMatrix &x, 
					       int verbosity )
{
  RealMatrix A_copy( A ), C_copy( C );
  RealVector b_copy( b ), d_copy( d );

  int M( A_copy.numRows() ), N( A_copy.numCols() ), lda( A_copy.stride() ), 
    ldc( C_copy.stride() );

  x.shapeUninitialized( N, 1 );

  Teuchos::LAPACK<int, Real> la;
  double* work;    // LAPACK work array
  int info( 0 );   // LAPACK output flag
  int lwork; // size of LAPACK work array
  int num_cons( C_copy.numRows() ); // number of equality constraints

  // Get the optimal work array size
  lwork = -1; // special code for workspace query
  work  = new double [1]; // temporary work array
  la.GGLSE( M, N, num_cons, A_copy.values(), lda, 
	    C_copy.values(), ldc, b_copy.values(), d_copy.values(), x.values(), 
	    work, lwork, &info );
  lwork = (int)work[0]; // optimal work array size returned by query
  delete [] work;
  work  = new double [lwork]; // Optimal work array

  // Least squares computation using LAPACK's DGGLSE subroutine which uses
  // a GRQ factorization method for solving the eq-constrained LLS problem
  info = 0;
  la.GGLSE( M, N, num_cons, A_copy.values(), lda, C_copy.values(),
	    ldc, b_copy.values(), d_copy.values(), x.values(), 
	    work, lwork, &info );

  if ( info < 0 )
    {
      std::stringstream msg;
      msg << "equality_constrained_least_squares() dgglse failed. ";
      msg << "The " << std::abs( info ) << "-th argument had an ";
      msg << "illegal value";
      throw( std::runtime_error( msg.str() ) );
    }
  if ( info == 1 )
    {
      std::stringstream msg;
      msg << "the upper triangular factor R associated with C in the ";
      msg << "generalized RQ factorization of the pair (C, A) is ";
      msg << "singular, so that rank(C) < num_cons; the least squares ";
      msg << "solution could not be computed.";
      throw( std::runtime_error( msg.str() ) );
    }
  if ( info == 2 )
    {
      std::stringstream msg;
      msg << "the (N-P) by (N-P) part of the upper trapezoidal factor ";
      msg << "T associated with A in the generalized RQ factorization ";
      msg << "of the pair (C, A) is singular, so that\n";
      msg << "rank( (A) ) < N; the least squares solution could not\n";
      msg << "    ( (C) )\n";
      msg << "be computed.";
      throw( std::runtime_error( msg.str() ) );
    }
};


} // namespace Pecos

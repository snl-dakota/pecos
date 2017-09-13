#include "nesta.hpp"

void create_orthogonal_matrix( const RealMatrix &Amatrix, const RealMatrix &rhs,
			       RealMatrix &scaled_matrix, 
			       RealMatrix &scaled_rhs){
  // QR = A' => A=R'Q' 
  // Ax = b  => R'Q'x = b => Q'x = inv(R')b
  // This only truly correct for the noiseless case.
  // It works for delta > 0, but the noise model is now
  // slightly different, so the constraints may be off.
  // instead use compute_svd_factorization()
  RealMatrix Qfactor, Rfactor;
  RealMatrix Amatrix_trans( Amatrix, Teuchos::TRANS );
  // Must use a qr factorizaztion such that if m>n, R has only n rows.
  qr_factorization( Amatrix_trans, Qfactor, Rfactor );
  // cannot use substitution solve as R is usuall not square
  qr_solve(Rfactor, rhs, scaled_rhs, Teuchos::TRANS );
  
  transpose(Qfactor,scaled_matrix);
}

void find_initial_solution( const RealMatrix &Amatrix, const RealMatrix &rhs,
			    RealMatrix & initial_solution ){
  initial_solution.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 
			     1.0, Amatrix, rhs, 0.0 );

}

void copy_matrix_col_to_vector( const RealMatrix &matrix, int col_num, 
				Real scalar, RealVector &vector){
  for (int i=0; i<matrix.numRows(); i++)
    vector[i] = matrix(i,col_num) * scalar;
}

void add( Real scalar1, const RealMatrix &matrix1, 
	  Real scalar2, const RealMatrix &matrix2, 
	  RealMatrix &result ){
  for (int j=0; j<matrix1.numCols(); j++)
    for (int i=0; i<matrix1.numRows(); i++)
      result(i,j) = scalar1 * matrix1(i,j) + scalar2 * matrix2(i,j);
}

void add( Real scalar1, const RealMatrix &matrix1, 
	  Real scalar2, const RealMatrix &matrix2, 
	  Real scalar3, const RealMatrix &matrix3,
	  RealMatrix &result ){
  for (int j=0; j<matrix1.numCols(); j++)
    for (int i=0; i<matrix1.numRows(); i++)
      result(i,j) = scalar1 * matrix1(i,j) + scalar2 * matrix2(i,j) + 
	scalar3 * matrix3(i,j);
}

// replace mu with hubber_func_tol
Real evaluate_l1_constraint( const RealMatrix &solution, Real mu, 
			     RealVector &objective_deriv ){
  // in python code
  // obj_deriv = sol.copy() / numpy.maximum( mu, abs( sol ) );
  // objective = ( numpy.dot( obj_deriv , sol ) - mu*dual_prox_func ); 

  Teuchos::BLAS<int, Real> blas;
  int len_sol = solution.numRows();
  for (int i=0;i<len_sol;i++)
    objective_deriv[i] = solution(i,0) / std::max( mu,std::abs(solution(i,0)) );
        
  // Compute dual prox function p_d(dual_sol)
  Real objective_deriv_norm = objective_deriv.normFrobenius(); 
  Real dual_prox_func = 0.5*objective_deriv_norm*objective_deriv_norm;
  // Compute objective which is based on the huber function
  Real objective = blas.DOT( len_sol, objective_deriv.values(), 1, 
			      solution.values(), 1 ) - mu*dual_prox_func;
  return objective;
}

void nesta( const RealMatrix &Amatrix, const RealMatrix &rhs, 
	    Real muf, Real delta, Real tol_var, 
	    int max_continuation_iters, int verbosity, 
	    const RealMatrix &warm_start_solution, RealMatrix &solution, 
	    RealMatrix &metrics ){

  if ( Amatrix.numRows() != rhs.numRows() )
    throw( std::runtime_error("nesta: Matrix and RHS are inconistent") );

  int num_cols = Amatrix.numCols();
  solution.shapeUninitialized( num_cols, 1 );

  if ( warm_start_solution.numRows() == 0 )
    find_initial_solution( Amatrix, rhs, solution );
  else{
    solution.assign( warm_start_solution );
  }

  int max_nesterov_iters = 1e4;

  int sol_abs_max_idx = magnitude_argmax( num_cols, solution.values() );
  Real mu0 = 0.9 * std::abs( solution(sol_abs_max_idx,0) );
  Real gamma = std::pow(muf/mu0,1./(Real)max_continuation_iters);
  Real gammat = std::pow(10.*tol_var,1./(Real)max_continuation_iters);
  tol_var = 0.1;
  Real mu = mu0;
  int num_total_iters = 0;
  RealMatrix initial_solution( solution.numRows(), 1, false );
  for ( int continuation_it=0; continuation_it<max_continuation_iters; continuation_it++ ){
    mu *= gamma;
    tol_var *= gammat;
    initial_solution.assign( solution );
    if ( verbosity > 0 )
      std::cout << "\tBeginning l1 Minimization; mu = " << mu << std::endl;
    int num_nesterov_iters = nesterov_l1_minimization( Amatrix, rhs, 
						       initial_solution, mu, 
						       delta, 
						       max_nesterov_iters,
						       tol_var, verbosity,
						       solution );
    num_total_iters += num_nesterov_iters;
  }
}

int nesterov_l1_minimization( const RealMatrix &Amatrix, 
			      const RealMatrix &rhs, 
			      const RealMatrix &initial_solution,
			      Real mu, Real delta, 
			      int max_nesterov_iters,
			      Real tol_var, int verbosity,
			      RealMatrix &solution ){
  if ( delta < 0. ){
    std::string msg = "nesterov_l1_minimization - ensure delta >= 0";
    throw( std::runtime_error( msg ) );
  }

  int num_rows = Amatrix.numRows(), num_cols = Amatrix.numCols();

  RealVector objective_deriv( num_cols, false );
  RealMatrix diff( num_rows, 1, false );
  RealMatrix cp( num_cols, 1, false );
  RealMatrix Acp( num_rows, 1, false );
  RealMatrix AtAcp( num_cols, 1, false );
  RealMatrix Ayk( num_rows, 1, false );
  RealMatrix Azk( num_rows, 1, false );
  RealMatrix Axk( num_rows, 1, false );
  RealMatrix solution_prev( num_cols, 1, false );

  RealMatrix Atb( num_cols, 1, false );
  Atb.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 
		1.0, Amatrix, rhs, 0.0 );
  RealMatrix residual( Teuchos::Copy, rhs, num_rows, 1 );
    
  solution.assign( initial_solution );
  RealVector wk( num_cols ); // initialize to zero
  RealMatrix sol_prev( Teuchos::Copy, solution, num_rows, 1 );
        
  Real Ak = 0;
  Real Lmu = 1. / mu;
  RealVector yk( Teuchos::Copy, solution.values(), num_cols );
  RealVector zk( Teuchos::Copy, solution.values(), num_cols );
  // smallest positive normalized floating-point number in IEEE 
  RealVector objective_history( 10, false );
  objective_history[0] = std::numeric_limits<double>::min();
  bool done = false;

  Real Lmu_inv = 1./Lmu;
  int it = 0;
  int num_obj_history = 0;
  if ( verbosity > 0 )
    std::printf(" Iter\t   fmu\t\tVar of fmu\t Residual\n"); 
  for ( it=0; it<max_nesterov_iters; it++ ){

    //////////////////
    // Dual problem //
    //////////////////
    Real objective = evaluate_l1_constraint( solution, mu, objective_deriv );

    ////////////////////
    // Primal problem //
    ////////////////////
    
    //update u_k

    // the following line is "q" in eq. (3.7) in the paper
    add( 1., solution, -Lmu_inv, objective_deriv, cp ); 
    Acp.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
		  1.0, Amatrix, cp, 0.0 );
    AtAcp.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 
		    1.0, Amatrix, Acp, 0.0 );

    if ( delta > 0 ){
      add( 1., rhs, -1., Acp, diff );
      Real Lambda = std::max(0., Lmu*( diff.normFrobenius()/delta-1.));
      Real gamma = Lambda / ( Lambda + Lmu );
      add( Lambda / Lmu * (1.-gamma), Atb, 1., cp, -gamma, AtAcp, yk );
      // for calculating the residual, we'll avoid calling A()
      // by storing A(yk) here (using A'*A = I):
      add(  Lambda / Lmu * (1.-gamma), rhs, 1., Acp, -gamma, Acp, Ayk );
    }else{
      // this assumes that A is a projection

      // TODO: Check what part of this is the projection step
      add( 1., Atb, 1., cp, -1., AtAcp, yk ); 
      Ayk.assign( rhs );
    }

    ///////////////////////////////////
    // Evaluate termination criteria //
    ///////////////////////////////////

    // only take mean over last 10 iterations
    Real mean_objective_history = mean( std::min( num_obj_history+1, 10 ),
					objective_history.values() );
    Real qp = std::numeric_limits<double>::max();
    if ( it != 0 )
      qp = std::abs( objective - mean_objective_history )/mean_objective_history;

    if ( (qp <= tol_var) && (done) ){ 
      if ( verbosity > 0 )
	std::cout << "Exiting: change in objective was below tolerance\n";
      break;
    }

    if ( (qp <= tol_var) && (!done) )
      done = true;

    if ( objective < std::numeric_limits<double>::epsilon() ){
      if ( verbosity > 0 )
	std::cout << "Exiting: Objective was smaller than machine precision\n";
      done = true;
      break;
    }

    // We only want to take the mean over last 10 iterations
    // so overwrite correct entry in history
    objective_history[(num_obj_history+1)%10] = objective;
    num_obj_history++;

    // update zk
    Real apk = 0.5*((Real)it+1.);
    Ak += apk; 
    Real tauk = 2./((Real)it+3.); 
           
    add( apk, objective_deriv, 1., wk, wk );          
    add( 1, initial_solution, -Lmu_inv, wk, cp );
    Acp.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
		  1.0, Amatrix, cp, 0.0 );
    AtAcp.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 
		    1.0, Amatrix, Acp, 0.0 );

    if ( delta > 0 ){
      add( 1., rhs, -1., Acp, diff );
      Real Lambda = std::max( 0.,Lmu * (diff.normFrobenius()/delta - 1.));
      Real gamma = Lambda / ( Lambda + Lmu );
      add( Lambda / Lmu*(1.-gamma), Atb, 1., cp, -gamma, AtAcp, zk );
      // for calculating the residual, we'll avoid calling A()
      // by storing A(zk) here (using A'*A = I):
      add( Lambda / Lmu*(1.-gamma), rhs, 1., Acp,  - gamma, Acp, Azk );
    }else{
      // this assumes A is a projection
      add( 1., Atb, 1., cp, -1., AtAcp, zk );
      Azk.assign( rhs );
    }

    solution_prev.assign( solution );
    add( tauk, zk, (1.-tauk), yk, solution );
    add( tauk, Azk,  (1.-tauk), Ayk, Axk );

    // prevent loss of precision
    if ( it%10 )
      Axk.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
		    1.0, Amatrix, solution, 0.0 );

    if ( verbosity > 1 ){
      // compute the residual
      add( 1., rhs, -1., Axk, residual );
      Real residual_norm = residual.normFrobenius();
      std::printf("%3d\t %.3e\t %.2e\t %.2e\n",it+1,objective,qp,residual_norm); 
    }
  }
  return it+1;
}

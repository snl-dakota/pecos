#include "spgl1.hpp"
#include "nesta.hpp" // needed to import addition of two matrices. Consider moving this utility to higher level

// verbosity levels
// 0 SILENT
// 1 SUMMARY STATEMENT
// 2 LASSO STEP INFO
// 3 LINESEARCH INFO

void filter( RealVector &vec ){
  int num_entries = vec.length();
  for (int i=0; i<num_entries; i++)
    vec[i] = ( (std::abs( vec[i] ) > 100*std::numeric_limits<double>::epsilon())
	       ? vec[i] : 0 );
}

void SPGL1Solver::
find_initial_solution( const RealVector &rhs, RealVector &initial_solution ) {
  //matrix_multiply( rhs, Teuchos::TRANS, 0.0, 1.0, initial_solution );
  initial_solution = 0.;
};

void SPGL1Solver::
project( const RealVector &vec, Real lasso_param, RealVector &vec_proj ){
  Real eps = std::numeric_limits<double>::epsilon();

  if ( vec_proj.length() != vec.length() )
    vec_proj.sizeUninitialized( vec.length() );

  int num_entries = vec.length();
  Real vec_norm = one_norm( vec );

  if ( lasso_param >= vec_norm ){
    vec_proj = vec;
    return;
  }
  else if ( lasso_param < eps ){
    vec_proj = 0.;
    return;
  }

  for (int i=0; i<num_entries; i++)
    vec_proj[i] = std::abs( vec[i] );
      
  IntVector indices;
  RealVector vec_abs;
  magnitude_sort( vec_proj, indices, vec_abs );
  
  Real infeasiable_mu = -lasso_param;
  Real gamma_final = 0.;
  Real gamma = 0.;
  for (int i=0; i<num_entries; i++){
    // accumuate the infeasiability Eq (4.2)
    infeasiable_mu = infeasiable_mu + vec_abs[i];
    //define the current solution of Eq (4.3)
    gamma = infeasiable_mu / (double) (i+1);
    if ( gamma >= vec_abs[i] ) 
      // The remaining iterations all satisfy Eq (4.5). I.e.
      // solution_abs[k+1:] = 0.
      break;
    gamma_final = gamma;
  }

  for (int i=0; i<num_entries; i++){
    vec_proj[indices[i]] = sign_fortran( std::max( 0., vec_abs[i]-gamma_final ),
					 vec[indices[i]] );
  }
}

void SPGL1Solver::
num_active_variables( const RealVector &sol, 
		      const RealVector &grad, 
		      const IntVector &non_zero_indices,
		      Real opt_tol,
		      int &nnz_sol, int &nnz_new_indices, int &nnz_diff,
		      IntVector &new_non_zero_indices ) const{

  int num_entries = sol.length();
  new_non_zero_indices.sizeUninitialized(num_entries);

  Real sol_tol = std::min( .1, 10*opt_tol );
  Real grad_tol = std::min( .1, 10*opt_tol );
  Real grad_norm = dual_norm( grad );

  // The nnz counts are based on simple primal indicator.
  nnz_sol = 0;
  nnz_new_indices = 0;
  nnz_diff = 0;
  for (int i=0; i<num_entries; i++){
    // Primal/dual based indicators
    if ( ( ( sol[i] >  sol_tol ) && ( grad[i]+grad_norm < grad_tol ) ) 
	 || ( ( sol[i] < -sol_tol ) && ( grad_norm-grad[i] < grad_tol ) ) )
      new_non_zero_indices[i] = 1;
    else
      new_non_zero_indices[i] = 0;
       
    nnz_new_indices += new_non_zero_indices[i];
    if ( std::abs( sol[i] ) >= sol_tol ) nnz_sol++;
    if ( new_non_zero_indices[i] != non_zero_indices[i] ) nnz_diff++;
  }
  if ( nnz_new_indices == 0 ) nnz_diff = -inf;
};

void SPGL1Solver::
matrix_multiply( const RealVector &vec, Teuchos::ETransp trans, Real alpha,
		 Real beta, RealVector &result ){

  if ( (trans == Teuchos::NO_TRANS) && (result.length() != AMatrix_.numRows()) )
    result.size(  AMatrix_.numRows() );
  if ( (trans == Teuchos::TRANS) && (result.length() != AMatrix_.numCols()) )
    result.size(  AMatrix_.numCols() );

  result.multiply( trans, Teuchos::NO_TRANS, alpha, AMatrix_, vec, beta );

  if (trans==Teuchos::NO_TRANS)
    counters_.numAMatrixProd_++;
  else
    // A'*x is only called when updating the gradient.
    // Thus the number of products with A' is number of iterations is 1 or 2
    // +1 is due to finding initial solution + 1 is due to resetting
    // to best solution (bpdn only)
    counters_.numAMatrixTransProd_++;
};

Real SPGL1Solver::
spgl1( const RealMatrix &Amatrix, const RealVector &rhs, Real lasso_param,
       Real residual_tol, RealVector &solution ){
  if ( Amatrix.numRows() != rhs.numRows() )
    throw( std::runtime_error("spgl1: Matrix and RHS are inconistent") );

  if ( verbosity_ > 1 ){
    std::printf("optTol: %1.1e; bpTol : %1.1e;",solverTol_,solverSecondaryTol_);
    if ( !solveSingleLasso_ )
      std::printf(" newtonTol: %1.1e",newtonUpdateTol_);
    std::printf("\n");
  }
  
  bpTol_ = solverSecondaryTol_;

  // Store copy of Amatrix to allow for generic matrix multiplication
  // This will eventually allow for passing matrix vector operators 
  // must be a deep copy.
  AMatrix_.shapeUninitialized( Amatrix.numRows(), Amatrix.numCols() ) ;
  AMatrix_.assign( Amatrix );

  int num_cols = AMatrix_.numCols();

  counters_.initialize(); 
  apply_user_options( num_cols, lasso_param, residual_tol );

  solution.sizeUninitialized( num_cols );
  // Must be called after counters are initialized
  find_initial_solution( rhs, solution );

  IntVector nonzero_indices;
  nonzero_indices.size( num_cols ); // initialize to zero
  
  RealVector prev_objs( numPrevObjStore_, false );
  prev_objs = -inf;

  int exit_code = CONTINUE;

  Real proj_grad_step_len = 1.;
  Real rhs_norm = rhs.normFrobenius();
  if ( rhs_norm <= residual_tol ){
    if ( verbosity_ > 1 )
      std::cout <<  "residual_tol >= ||rhs||. The exact solution is x = 0\n";
    solution = 0.;
    return rhs_norm;
  }

  solWorkVec_.size( solution.length() ); // initialize to zero
  project( solWorkVec_, lasso_param, solution );
  RealVector residual( rhs );
  RealVector gradient( solution.length(), false );
  matrix_multiply( solution, Teuchos::NO_TRANS, -1.0, 1.0, residual );
  matrix_multiply( residual, Teuchos::TRANS, -1.0, 0.0, gradient );

  //Real residual_norm = residual.normFrobenius();
  //Real objective = residual_norm*residual_norm / 2.;
  Real objective = residual.dot( residual ) / 2.;
  Real residual_norm = std::sqrt( objective * 2. );
  // print initial objective and gradient norm
  //std::printf( "%1.16e\t%1.16e\n",residual_norm,dual_norm( gradient ) );

  prev_objs[0] = objective;
  Real best_objective = objective;
  Real prev_objective  = objective;

  RealVector delta_sol;
  add( 1.0, solution, -1.0, gradient, solWorkVec_ );
  project( solWorkVec_, lasso_param, delta_sol );
  solWorkVec_ = delta_sol;
  add( 1.0, solWorkVec_, -1.0, solution, delta_sol );
  Real delta_sol_norm = delta_sol.normInf(); 

  Real gradient_step = 0.;
  if ( delta_sol_norm < ( 1. / maxGradStepSize_ ) )
    gradient_step = maxGradStepSize_;
  else
    gradient_step = std::min(  maxGradStepSize_, 
			       std::max( minGradStepSize_, 
					 1. / delta_sol_norm ) );

  if ( verbosity_ > 1 )
    std::printf("\t %5s  %13s  %13s  %9s  %9s  %6s  %6s  %6s\n",
                "Iter","Objective","Relative Gap","Rel Error","Grad Norm",
		"Grad Step","NNZ Sol","NNZ Grad");

  while (true){
    residual_norm = residual.normFrobenius();
    Real lagrange_multiplier = dual_norm( gradient );// should be -gradient but dual norm does not care about sign
    Real residual_abs_error = residual_norm - residual_tol;

    if ( !solveSingleLasso_ ){
      Real prev_lasso_param = lasso_param;
      lasso_param = std::max(0.,lasso_param+(residual_norm*residual_abs_error)/
			     lagrange_multiplier );
      // --------------------------------------------
      // Take next newton step
      // --------------------------------------------
      counters_.newtonIters_++;
      if ( lasso_param < prev_lasso_param ){
	// The one-norm ball has decreased.  Need to make sure 
	// that the next iterate if feasible, which we do by 
	// projecting it.
	solWorkVec_ = solution;
	project( solWorkVec_, lasso_param, solution );
      }
      
      if ( verbosity_ > 1 )
	std::printf("\n Newton Step %d:\t Lasso_Param=%1.15e\n", 
		    counters_.newtonIters_, lasso_param );
      
      // else:
      // No need for newton iterations as we are solving the lasso problem
    }else{
      if ( verbosity_ > 1 )
	std::printf("\n Solving the Lasso:\t Lasso_Param=%1.15e\n", lasso_param);
    }

    spectral_projected_gradient( rhs, residual, gradient, lasso_param,
				 objective, lagrange_multiplier, solution, 
				 residual_tol, nonzero_indices, rhs_norm,
				 prev_objective, proj_grad_step_len, exit_code,
				 best_objective, gradient_step, prev_objs );

    if ( exit_code ) break;
  }
  print_termination_message( exit_code );

  return residual.normFrobenius();
}

void SPGL1Solver::
spectral_projected_gradient( const RealVector &rhs, RealVector &residual,
			     RealVector &gradient, Real lasso_param, 
			     Real &objective, Real &lagrange_multiplier,
			     RealVector &solution, Real residual_tol, 
			     IntVector &nonzero_indices, Real &rhs_norm,
			     Real &prev_objective, 
			     Real &proj_grad_step_len, int &exit_code,
			     Real &best_objective, Real &gradient_step, 
			     RealVector &prev_objs ){

  // set optTol from generic solverTol of LinearSolver base class
  optTol_ = solverTol_;
  // set optTol from generic solverTertiaryTol of LinearSolver base class
  newtonUpdateTol_ = solverTertiaryTol_;

  RealVector prev_solution;
  RealVector prev_gradient;
  RealVector prev_residual;
  RealVector best_solution;
  Real residual_norm = 0.;
  while ( true ){
    counters_.SPGIters_++;

    prev_solution = solution;
    Real prev_objective = objective; 
    prev_gradient = gradient; 
    prev_residual = residual;
    // ---------------------------------------------------------------
    // Projected gradient step and linesearch.
    // ---------------------------------------------------------------
    linesearch( rhs, lasso_param, prev_solution, gradient, gradient_step,
		prev_objs, solution, objective, residual, proj_grad_step_len,
		exit_code );

    if ( one_norm( solution ) > lasso_param + optTol_ ){
      std::string msg = "spectral_projected_gradient: primal norm infeasiable";
      throw( std::runtime_error( msg ) );
    }

    // ---------------------------------------------------------------
    // Update gradient and compute new Barzilai-Borwein scaling.
    // ---------------------------------------------------------------
    update_gradient( solution, prev_solution, prev_gradient, residual, 
		     gradient, gradient_step );
    
    // ---------------------------------------------------------------
    // Update function history.
    // ---------------------------------------------------------------
    if ( objective > residual_tol*residual_tol / 2. ){
      // Don't update if superoptimal.
      prev_objs[counters_.SPGIters_%numPrevObjStore_] = objective;
      if ( best_objective > objective ){
	best_objective = objective;
	best_solution = solution;
      }
    }
    
    // TODO: Check if I need these two lines    
    Real lagrange_multiplier = dual_norm( gradient );
    residual_norm = residual.normFrobenius();

    bool update_lasso_param = false;
    test_exit_conditions( objective, lasso_param, solution, gradient, rhs, 
			  residual, residual_norm, lagrange_multiplier, 
			  residual_tol, nonzero_indices, rhs_norm,
			  prev_objective, exit_code, proj_grad_step_len,
			  update_lasso_param );

    if ( ( exit_code ) || ( update_lasso_param ) ) break;

    //residualNorms_.push_back(residual_norm );

  }

  // ---------------------------------------------------------------
  // Restore best solution (only if solving single problem).
  // ---------------------------------------------------------------
  if ( ( !solveSingleLasso_ ) && ( objective > best_objective ) && ( exit_code ) ){

    if ( verbosity_ > 2 ){
      residual_norm = std::sqrt( 2. * best_objective );
      std::printf( "Restoring best iterate to objective %13.7e\n",residual_norm);
    }

    solution = best_solution;
    residual =  rhs;
    matrix_multiply( solution, Teuchos::NO_TRANS, -1.0, 1.0, residual );
    matrix_multiply( residual, Teuchos::TRANS, -1.0, 0.0, gradient );
    
    lagrange_multiplier = dual_norm( gradient );
    residual_norm = residual.normFrobenius();
  }  
}



void SPGL1Solver::
print_termination_message( int exit_code ){
  if ( verbosity_ > 0 ){
    if ( exit_code == EXIT_ITERATIONS ) 
      std::cout << "exiting: maximum iterations reached\n";
    else if (exit_code==EXIT_ROOT_FOUND)
      std::cout << "exiting: root found\n"; 
    else if ( (exit_code==EXIT_BPSOL1_FOUND) || (exit_code==EXIT_BPSOL2_FOUND))
      std::cout << "exiting: found a BP solution\n";
    else if (exit_code==EXIT_OPTIMAL) 
      std::cout << "exiting: optimal solution found\n";
    else if (exit_code== EXIT_LINE_ERROR) 
      std::cout << "exiting: line search error\n";
    else if (exit_code==EXIT_SUBOPTIMAL_BP) 
      std::cout << "exiting: suboptimal bp solution\n";
    else if (exit_code==EXIT_MATVEC_LIMIT) 
      std::cout << "exiting: maimum mat-vec operations reached\n";
    if (exit_code==EXIT_NO_NNZ_CHANGE)
      std::cout << "exiting: found a possible active set\n";
    // else incorrect exit code given

    std::printf("Products with A: %d\n", counters_.numAMatrixProd_);
    std::printf("Products with trans(A): %d\n",counters_.numAMatrixTransProd_);
    std::printf("Newton iterations: %d\n",counters_.newtonIters_ );
    std::printf("Line search iterations: %d\n", counters_.totalLineIters_ );
    //std::printf("Computation time: %1.2f\n"; ,time.time()-t0);
  }
}

void SPGL1Solver::
spg_projected_linesearch( const RealVector &rhs, const RealVector &descent_dir, 
			  Real max_prev_objs, Real lasso_param,
			  RealVector &solution,  RealVector &residual,
			  Real &objective, int &proj_linesearch_iter,
			  Real &proj_grad_step_len, int &linesearch_exit_code ){

  int num_terms = solution.length();
  
  // initial stepsize from current solution of linesearch
  Real step_size = 1.;
  // scale to safeguard projection against large step sizes
  Real safeguard_scale = 1. ;
  // norm of change in solution of current linesearch interate
  Real norm_delta_sol = 0.;
  // the number of line search iterations
  proj_linesearch_iter = 0;
  // the number of iterations that must be safe guarded
  int safegaurd_iter = 0;

  if ( verbosity_ > 2 )
    std::printf( "\t\t %5s  %13s  %13s  %13s  %8s\n", "iters","obj","step size",
		 "detla obj","scale" );  
  RealVector new_solution( solution.length(), false );
  RealVector new_residual( rhs.length(), false );
  RealVector delta_sol( solution.length(), false );
  Real new_objective = 0.;
  while ( true ){

    // find the candidate point of the line search
    // x_{k+1} = P( x_{k}-lambda_k*sigma*d_{k} )
    add( 1.0, solution, -step_size*safeguard_scale, descent_dir, solWorkVec_ );
    project( solWorkVec_, lasso_param, new_solution );
    
    // update the residual
    // r_{k+1} = b - A*x_{k+1}
    new_residual.assign( rhs );
    matrix_multiply( new_solution, Teuchos::NO_TRANS, -1.0, 1.0, new_residual );

    // compute new objective value
    ///Real residual_norm = new_residual.normFrobenius();
    //new_objective = residual_norm*residual_norm / 2.;
    new_objective = new_residual.dot( new_residual ) / 2.;
  
    // compute change in solution
    // s_{k} = x_{k+1} - x_{k}
    add( 1.0, new_solution, -1.0, solution, delta_sol );

    // compute change in objective
    Real delta_obj = safeguard_scale * descent_dir.dot( delta_sol );

    //if ( numpy.any( delta_obj >= 0. ) ){
    if ( delta_obj >= 0. ){
      // exit: the line seach is not descending
      linesearch_exit_code = LS_NOT_DESCENDING;
      break;
    }

    if ( verbosity_ > 2 )
      std::printf( "\t\t %2d  %13.16e  %13.7e  %13.16e  %8.1e\n", 
		   proj_linesearch_iter, new_objective, step_size, delta_obj,
		   safeguard_scale );

    // check if there has been a sufficient decrease in the objective
    // This is the nonmonotone Armijo condition
    if (new_objective<max_prev_objs+sufficientDescentTol_*step_size*delta_obj){
      // exit: the line search worked
      linesearch_exit_code = LS_CONVERGED;
      break;
    }

    if ( proj_linesearch_iter >= maxProjLinesearchIters_ ){
      // exit: too many line search iterations
      linesearch_exit_code = LS_MAX_ITERATIONS_REACHED;
      break;
    }

    proj_linesearch_iter++;

    // decrease step length
    step_size /= 2.;

    // safeguard: If step_size is huge, then even damped search
    // directions can give exactly the same point after projection.  If
    // we observe this in adjacent iterations, we drastically damp the
    // next search direction.
    Real prev_norm_delta_sol  = norm_delta_sol;
    norm_delta_sol = delta_sol.normFrobenius() / std::sqrt(num_terms);
    if ( std::abs( norm_delta_sol - prev_norm_delta_sol ) 
	 <= 1e-6 * norm_delta_sol ){
      Real descent_dir_norm = descent_dir.normFrobenius() / std::sqrt(num_terms);
      safeguard_scale = norm_delta_sol / descent_dir_norm / 
	std::pow(2.,safegaurd_iter);
      safegaurd_iter++;
    }
  }
  solution = new_solution;
  objective = new_objective;
  residual = new_residual;
  proj_grad_step_len = step_size;
  
}

void SPGL1Solver::
spg_non_monotone_linesearch( const RealVector &rhs, 
			     const RealVector &delta_sol, 
			     Real gradient_delta_sol, 
			     Real max_prev_objs, 
			     RealVector &solution, RealVector &residual,
			     Real &objective, int &nm_linesearch_iters,
			     int &linesearch_exit_code){

  // initial stepsize from current solution of linesearch
  Real step_size = 1.;
  // the number of line search iterations
  nm_linesearch_iters   = 0;

  gradient_delta_sol = -std::abs(gradient_delta_sol); //if complex numbers

  if ( verbosity_ > 2 )
    std::printf( "\t\t %5s  %13s  %13s  %13s\n", "iters","obj","step size",
		 "detla obj" );  
  RealVector new_solution( solution.length(), true );
  RealVector new_residual;
  Real new_objective = 0.;
  while ( true ){
    
    // Evaluate trial point and function value.
    // x_{k+1} = x_{k} + alpha*dk
    add( 1.0, solution, step_size, delta_sol, new_solution );

    // update the residual
    // r_{k+1} = b - A*x_{k+1}
    new_residual.assign( rhs );
    matrix_multiply( new_solution, Teuchos::NO_TRANS, -1.0, 1.0, new_residual );

    // compute new objective value
    //Real residual_norm = new_residual.normFrobenius();
    //Real new_objective = residual_norm*residual_norm / 2.;
    new_objective = new_residual.dot( new_residual ) / 2.;

    if ( verbosity_ > 2 )
      std::printf( "\t\t %2d  %13.16e  %13.7e  %13.16e\n", 
		   nm_linesearch_iters, new_objective, step_size, 
		   gradient_delta_sol );


    // check exit conditions
    if ( new_objective < max_prev_objs+
	 sufficientDescentTol_*step_size*gradient_delta_sol ){
      // sufficient descent condition
      linesearch_exit_code = LS_CONVERGED;
      break;
    }else if ( nm_linesearch_iters >= maxNMLinesearchIters_ ){
      // too many line search iterations
      linesearch_exit_code = LS_MAX_ITERATIONS_REACHED;
      break;
    }

    nm_linesearch_iters++;

    // safeguarded quadratic interpolation.
    Real tmp = 0.;
    if ( step_size <= 0.1 )
      step_size  = step_size / 2.;
    else
      tmp = ( (-gradient_delta_sol*step_size*step_size) / 
	      (2.*(new_objective-objective-step_size*gradient_delta_sol ) ) );

    if ( ( tmp < 0.1 ) || ( tmp > 0.9*step_size ) ||( isnan( tmp ) ) )
      tmp = step_size / 2.;

    step_size = tmp;
  }
  solution = new_solution;
  objective = new_objective;
  residual = new_residual;
}

void SPGL1Solver::
linesearch( const RealVector &rhs, Real lasso_param, 
	    const RealVector &prev_solution, 
	    const RealVector &gradient, 
	    Real gradient_step,
	    const RealVector &prev_objs, 
	    RealVector &solution, 
	    Real &objective,
	    RealVector &residual,
	    Real &proj_grad_step_len,
	    int &exit_code ){
  RealVector delta_sol( solution.length(), false );
  RealVector descent_dir( gradient );
  descent_dir *= gradient_step;

  int  proj_linesearch_iters = 0;
  int linesearch_exit_code = LS_NOT_DESCENDING;
  int argmax_prev_objs = argmax( prev_objs.length(), prev_objs.values() );
  Real max_prev_objs = prev_objs[argmax_prev_objs];
  // Typical spectral projected gradient  method does not include this step.
  // I have included it to be consistent with spgl1.m
  spg_projected_linesearch( rhs, descent_dir, 
			    max_prev_objs, lasso_param,
			    solution,  residual,
			    objective, proj_linesearch_iters,
			    proj_grad_step_len, 
			    linesearch_exit_code );
  
  counters_.totalLineIters_ += proj_linesearch_iters;

  if ( linesearch_exit_code ){
    if ( verbosity_ > 2 )
      std::cout << "projected backtrack failed.\n"
		<< "Retrying with non monotone line search\n";
    solution = prev_solution;
    
    // compute x_{k}-lambda_{k}*g_{k}
    add( 1.0, solution, -gradient_step, gradient, solWorkVec_ );
    // compute P(x_{k}-lambda_{k}*g_{k})
    project( solWorkVec_, lasso_param, delta_sol );

    // compute the search direction dk = P(x_{k}-lambda_{k}*g_{k})-x_{k}
    solWorkVec_ = delta_sol;
    add( 1.0, solWorkVec_, -1.0, solution, delta_sol );
    // todo (above) remove extra copy by using extra memory.
    
    Real gradient_delta_sol = gradient.dot( delta_sol );
    int nm_linesearch_iters = 0;
    spg_non_monotone_linesearch( rhs, delta_sol, gradient_delta_sol, 
				 max_prev_objs, solution, residual,
				 objective, nm_linesearch_iters,
				 linesearch_exit_code );
    counters_.totalLineIters_ += nm_linesearch_iters;
  }
  
  exit_code = linesearch_exit_code;
  if ( linesearch_exit_code ){
    // second line search failed. Revert to previous iterate and 
    // damp max BB step
    if ( counters_.numLineErrors_ >= maxLineErrors_  )
      exit_code = EXIT_LINE_ERROR;
    else{
      maxGradStepSize_ /= 10.;
      if ( verbosity_ > 2 )
	std::printf( "Linesearch failed with error %d. Damping max BB scaling to %6.1e.",linesearch_exit_code,maxGradStepSize_);
      counters_.numLineErrors_++;
      exit_code = 0;
    }
  }
}

void SPGL1Solver::
update_gradient( const RealVector &solution, const RealVector &prev_solution,
		 const RealVector &prev_gradient, const RealVector &residual, 
		 RealVector &gradient, Real &gradient_step){

  // compute gradient g_{k+1}
  matrix_multiply( residual, Teuchos::TRANS, -1.0, 0.0, gradient );

  // -------------------------------------------- //
  // Algorithm 2.1 Step 3 (p5 Birgin_MR_JSS_2014) //
  // -------------------------------------------- //

  // compute change in solution
  // s_{k} = x_{k+1} - x_{k}
  // todo try not to reallocate this vector every time
  RealVector delta_sol( solution.length(), false ); 
  add( 1.0, solution, -1.0, prev_solution, delta_sol );

  // Compute change in gradient 
  // y_{k} = g_{k+1} - g_{k}
  RealVector delta_grad( gradient.length(), false );
  add( 1.0, gradient, -1.0, prev_gradient, delta_grad );

  // compute ||s_{k}||^2
  //Real delta_sol_norm = delta_sol.normFrobenius();
  //Real delta_sol_norm_sq = delta_sol_norm*delta_sol_norm;
  Real delta_sol_norm_sq = delta_sol.dot( delta_sol );

  // compute s_{k}^T*y_{k}
  Real delta_sol_delta_grad = delta_sol.dot( delta_grad );

  if ( delta_sol_delta_grad <= 0. ){
    gradient_step = maxGradStepSize_;
    //throw(std::runtime_error("I have not found a situation where this happens"));
  }else
    gradient_step = std::min( maxGradStepSize_, 
			      std::max(minGradStepSize_, 
				       delta_sol_norm_sq/delta_sol_delta_grad));
  // -------------------------------------------- //
  // End Step 3
  // -------------------------------------------- //
  
}

void SPGL1Solver::
test_exit_conditions( Real objective, Real lasso_param, 
		      const RealVector &solution, const RealVector & gradient, 
		      const RealVector &rhs, const RealVector &residual, 
		      Real residual_norm, Real lagrange_multiplier, 
		      Real residual_tol, const IntVector &nonzero_indices, 
		      Real rhs_norm, Real prev_objective,
		      int &exit_code, Real proj_grad_step_len,
		      bool &update_lasso_param ){

  RealVector res_work_vec( residual.length(), false );
  res_work_vec = residual;
  res_work_vec -= rhs;
  Real duality_gap= residual.dot(res_work_vec)+lasso_param*lagrange_multiplier;
  Real duality_gap_rel = std::abs( duality_gap ) / std::max( 1., objective );

  Real residual_abs_error = residual_norm - residual_tol;
  Real objective_abs_error = objective - residual_tol * residual_tol / 2.;
  Real residual_rel_error = std::abs(residual_abs_error)/
    std::max(1.,residual_norm);
  Real objective_rel_error = ( std::abs( objective_abs_error ) / 
			       std::max( 1., objective ) );

  // Only nnz_change is used to determine exit conditions.
  // The other nnz_ variables are used for subspace minimization
  // which currently I have turned off
  int nnz_sol = 0, nnz_new_indices = 0, nnz_diff = 0;
  IntVector new_nonzero_indices;
  num_active_variables( solution, gradient, nonzero_indices,
			optTol_, nnz_sol, nnz_new_indices, nnz_diff,
			new_nonzero_indices );

  if ( nnz_diff )
    // number of iterations with fixed set of non-zeros
    counters_.numFixedNNZIters_= 0; 
  else{
    counters_.numFixedNNZIters_++;
    if (counters_.numFixedNNZIters_ >= maxNumFixedNNZIters_ )
      exit_code = EXIT_NO_NNZ_CHANGE;
  }
  
  if ( solveSingleLasso_ ){
    if ( ( duality_gap_rel < optTol_ ) || ( residual_norm < optTol_*rhs_norm ) )
      exit_code = EXIT_OPTIMAL;
  }else
    // Instead of using a fixed threshold for solving the current
    // lasso problem we compare the duality gap to the relative error
    // in the objective
    test_exit_conditions_for_multiple_lasso_param(duality_gap_rel, 
						  objective_rel_error, 
						  residual_rel_error, 
						  residual_norm, 
						  rhs_norm, lagrange_multiplier,
						  residual_tol, objective, 
						  prev_objective,
						  update_lasso_param,
						  exit_code );
        
  if ( (!exit_code) && (counters_.SPGIters_ >= maxSPGIters_) )
    exit_code = EXIT_ITERATIONS;

  if ( verbosity_ > 1 )
    std::printf( "\t %5i  %13.16e  %13.16e  %9.2e  %9.16e  %6.1f  %6i  %6i\n",
		 counters_.SPGIters_, residual_norm, duality_gap_rel,
		 residual_rel_error, lagrange_multiplier,
		 std::log10( proj_grad_step_len ), nnz_sol,
		 nnz_new_indices );
  //TODO linesearch prints objective but I print residual_norm here
  // make consistent. residual_norm = sqrt(objective*2)

}

void SPGL1Solver::
apply_user_options( int num_cols, Real &lasso_param, Real &residual_tol ){
  if ( maxSPGIters_ <= 0 )
    maxSPGIters_ = 10*num_cols;
  
  if ( ( lasso_param < 0 ) && ( residual_tol < 0 ) ){
    // Solve the BP problem
    solveSingleLasso_ = false;
    lasso_param = 0.;
    residual_tol = 0.;
  }else if ( residual_tol < 0 )
    // Solve a single lasso problem
    solveSingleLasso_ = true;
  else{
    // Solve the BPDN problem
    if ( lasso_param < 0 )
      lasso_param = 0.;
    solveSingleLasso_ = false;
  }
}

void SPGL1Solver::
test_exit_conditions_for_multiple_lasso_param( Real duality_gap_rel, 
					       Real objective_rel_error, 
					       Real residual_rel_error, 
					       Real residual_norm, 
					       Real rhs_norm, 
					       Real lagrange_multiplier, 
					       Real residual_tol, 
					       Real objective, 
					       Real prev_objective, 
					       bool &update_lasso_param, 
					       int &exit_code){
  if ( ( duality_gap_rel <= std::max(optTol_,objective_rel_error) ) ||
       ( residual_rel_error <= optTol_ ) ){
    // The problem is nearly is solved to the accuracy for the current
    // lasso parameter. Check optimality of the current root for the 
    // Newton problem.
    if ( residual_norm <= bpTol_ * rhs_norm )
      exit_code = EXIT_BPSOL1_FOUND;
    if (lagrange_multiplier <= bpTol_ * rhs_norm )
      exit_code = EXIT_BPSOL2_FOUND;
    if ( residual_rel_error <= optTol_ )
      exit_code = EXIT_ROOT_FOUND;
    if (residual_norm <= residual_tol)
      exit_code = EXIT_SUBOPTIMAL_BP;
  }
  
  bool test_rel_change_1 = (std::abs(objective-prev_objective)<=
			    newtonUpdateTol_*objective);
  bool test_rel_change_2=
    (std::abs(objective-prev_objective)<=
     newtonUpdateTol_*objective*(std::abs(residual_norm-residual_tol)));
     //1.e-1*objective*(std::abs(residual_norm-residual_tol)));
  // The uncommented line above was used in spgl1.m but I think 1e-1 is a bug

  update_lasso_param =
    (((test_rel_change_1 && (residual_norm > 2.*residual_tol ) ) ||
      (test_rel_change_2 && ( residual_norm <= 2.*residual_tol ) ) )
     && ( !exit_code ) );
}

void SPGL1Solver::
solve( RealMatrix &A, RealMatrix &B, RealMatrix &result_0, 
       RealMatrix &result_1 ){
  if ( B.numCols() != 1 )
    throw( std::runtime_error("SPGL1Solver::solve() B must be a vector") ); 
  
  // A and b must be copies or  = operator will make shallow copy
  // and cause havoc
  RealVector b( Teuchos::Copy, B[0], B.numRows() );
  RealMatrix A_copy( Teuchos::Copy, A, A.numRows(), A.numCols() );

  if ( A.numRows() > A.numCols() ){
    std::string msg = "SPGL1Solver::solve() A must be underdetermined";
    throw(std::runtime_error(msg)); 
  }

  RealVector column_norms;
  if ( normaliseInputs_ )
    normalise_columns( A_copy, column_norms );

  RealVector solution;
  result_0.shapeUninitialized( A.numCols(), residualTols_.length() );
  result_1.shapeUninitialized( 1, residualTols_.length() );
  for (int j=0; j<residualTols_.length(); j++){
    Real residual_norm = spgl1( A_copy, b, 0., residualTols_[j], solution );
    for (int i=0; i<solution.length(); i++)
      result_0(i,j) = solution[i];
    result_1(0,j) = residual_norm;
  }

  if ( normaliseInputs_ )
    adjust_coefficients( column_norms, result_0 );
};


void SPGL1Solver::set_max_iters( int max_iters )
{
  LinearSolver::set_max_iters( max_iters );
  maxSPGIters_ = max_iters;
};

void magnitude_sort( const RealVector &vec, IntVector &indices,
		     RealVector &abs_vec ){
  int num_entries = vec.length();
  std::vector<int> work_indices( num_entries );
  for ( int i = 0; i < num_entries; i++ )
    work_indices[i] = i;
  
  std::sort( work_indices.begin(), work_indices.end(), 
   magnitude_index_sorter<Teuchos::SerialDenseVector<int,double> >(vec) );
	     
  abs_vec.sizeUninitialized( num_entries );
  indices.sizeUninitialized( num_entries );
  for ( int i = 0; i < num_entries; i++ ){
    indices[i] = work_indices[i];
    abs_vec[i] = std::abs( vec[indices[i]] );
  }
}

// In exact Newton method for solving bpdn is highly sensitive to newtonUpdateTol. I found one example where rubbish when using 1e-4,1e-5 but when I used 1e-6 or smaller I got the correct result. A good test is use spgl1 to solve lasso by specifying the tau parameter. Then compute corresponding residual. Then make sure that spgl1 taking that residual reproduces the tau that generated it.

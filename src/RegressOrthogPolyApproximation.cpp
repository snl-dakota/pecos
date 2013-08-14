/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        RegressOrthogPolyApproximation
//- Description:  Implementation code for RegressOrthogPolyApproximation class
//-               
//- Owner:        John Jakeman

#include "RegressOrthogPolyApproximation.hpp"
#include "pecos_global_defs.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

// headers necessary for cross validation
#include "MathTools.hpp"
#include "CrossValidationIterator.hpp"
#include "LinearModelModules.hpp"
#include "CrossValidationModules.hpp"

//#define DEBUG


namespace Pecos {

int RegressOrthogPolyApproximation::min_coefficients() const
{
  // return the minimum number of data instances required to build the 
  // surface of multiple dimensions
  if (expConfigOptions.expansionCoeffFlag ||
      expConfigOptions.expansionCoeffGradFlag)
    // Now that L1-regression has been implemented. There is no longer a need 
    // to enforce a lower bound on the number of data instances.
    return 1;
    // numExpansionTerms is computed by the allocate_arrays() call in
    // compute_coefficients(), which is too late for use of this fn by
    // ApproximationInterface::minimum_samples() in DataFitSurrModel::
    // build_global(), so numExpansionTerms must be calculated.
    //return total_order_terms(approxOrder);
  else
    return 0;
}


void RegressOrthogPolyApproximation::compute_coefficients()
{
  if (!expConfigOptions.expansionCoeffFlag &&
      !expConfigOptions.expansionCoeffGradFlag) {
    PCerr << "Warning: neither expansion coefficients nor expansion "
	  << "coefficient gradients\n         are active in "
	  << "RegressOrthogPolyApproximation::compute_coefficients().\n"
	  << "         Bypassing approximation construction." << std::endl;
    return;
  }

  // For testing of anchor point logic:
  //size_t index = surrData.size() - 1;
  //surrData.anchor_point(surrData.variables_data()[index],
  //                      surrData.response_data()[index]);
  //surrData.pop(1);

  // anchor point, if present, is handled differently for different
  // expCoeffsSolnApproach settings:
  //   SAMPLING:   treat it as another data point
  //   QUADRATURE/CUBATURE/COMBINED_SPARSE_GRID: error
  //   LEAST_SQ_REGRESSION: use equality-constrained least squares
  size_t i, j, num_total_pts = surrData.size();
  if (surrData.anchor())
    ++num_total_pts;
  if (!num_total_pts) {
    PCerr << "Error: nonzero number of sample points required in RegressOrthog"
	  << "PolyApproximation::compute_coefficients()." << std::endl;
    abort_handler(-1);
  }

  // Array sizing can be divided into two parts:
  // > data used in all cases (size in allocate_arrays())
  // > data not used in expansion import case (size here)
  allocate_arrays();
#ifdef DEBUG
  gradient_check();
#endif // DEBUG

  // calculate polynomial chaos coefficients
  surrData.data_checks();
  regression();

  computedMean = computedVariance = 0;
}


/** In this case, regression is used in place of spectral projection.  That
    is, instead of calculating the PCE coefficients using inner products, 
    linear least squares is used to estimate the PCE coefficients which
    best match a set of response samples.  The least squares estimation is
    performed using DGELSS (SVD) or DGGLSE (equality-constrained) from
    LAPACK, based on anchor point and derivative data availability. */
void RegressOrthogPolyApproximation::regression()
{
  set_fault_info();
  
  CSOpts.solver = expConfigOptions.expCoeffsSolnApproach;
  bool fn_constrained_lls = (basisConfigOptions.useDerivs && 
			     faultInfo.constr_eqns &&
			     faultInfo.constr_eqns < numExpansionTerms);
  if (CSOpts.solver==DEFAULT_REGRESSION)
    if ((fn_constrained_lls || faultInfo.anchor_fn || faultInfo.anchor_grad) 
	&& (!faultInfo.under_determined))
      CSOpts.solver=EQ_CON_LEAST_SQ_REGRESSION;
    else if (!faultInfo.under_determined)
      CSOpts.solver=SVD_LEAST_SQ_REGRESSION;
    else 
      CSOpts.solver=LASSO_REGRESSION;

  if (CSOpts.solver==DEFAULT_LEAST_SQ_REGRESSION)
    if ((fn_constrained_lls || faultInfo.anchor_fn || faultInfo.anchor_grad) 
	&& (!faultInfo.under_determined))
      CSOpts.solver=EQ_CON_LEAST_SQ_REGRESSION;
    else
      CSOpts.solver=SVD_LEAST_SQ_REGRESSION;

  // Set solver parameters
  if ( CSOpts.solver == LASSO_REGRESSION )
    CSOpts.delta = l2Penalty;
  if ( noiseTols.length() > 0 )
    CSOpts.epsilon = noiseTols[0];
  else
    {
      noiseTols.size( 1 );
      noiseTols[0] = CSOpts.epsilon;
      if ( CSOpts.solver == BASIS_PURSUIT_DENOISING ) noiseTols[0] = 1e-3;
    }
  if ( CSOpts.solver != SVD_LEAST_SQ_REGRESSION )
      CSOpts.solverTolerance = expConfigOptions.convergenceTol;
  else
    CSOpts.solverTolerance = -1.0;
  CSOpts.verbosity = std::max(0, expConfigOptions.outputLevel - 1);
  if ( expConfigOptions.maxIterations > 0 )
    CSOpts.maxNumIterations = expConfigOptions.maxIterations;

  // Solve the regression problem using L1 or L2 minimization approaches
  bool regression_err = 0;
  if (CSOpts.solver==EQ_CON_LEAST_SQ_REGRESSION && !crossValidation){
    if ((fn_constrained_lls || faultInfo.anchor_fn || faultInfo.anchor_grad) 
	&& (!faultInfo.under_determined))
      {
	CSOpts.numFunctionSamples = surrData.size();
	run_regression();
      }
    else{
      PCout << "Could not perform equality constrained least-squares. ";
      if (faultInfo.under_determined){
	CSOpts.solver = LASSO_REGRESSION;
	PCout << "Using LASSO regression instead\n";
      }
      else
	{
	  CSOpts.solver = SVD_LEAST_SQ_REGRESSION;
	  PCout << "Using SVD least squares regression instead\n";
	}
      //regression_err = L2_regression(num_data_pts_fn, num_data_pts_grad, reuse_solver_data);
      run_regression();
    }
  }
  else{
    //regression_err = L2_regression(num_data_pts_fn, num_data_pts_grad, reuse_solver_data);
    run_regression();
  }

  if (regression_err) { // if numerical problems in LLS, abort with error
    PCerr << "Error: nonzero return code from least squares solution in "
	  << "RegressOrthogPolyApproximation::regression()" << std::endl;
    abort_handler(-1);
  }
}


void RegressOrthogPolyApproximation::set_fault_info()
{
  size_t constr_eqns, anchor_fn, anchor_grad, num_data_pts_fn,
    num_data_pts_grad, total_eqns, num_surr_data_pts;
  bool under_determined = false, reuse_solver_data;

  // compute order of data contained within surrData
  short data_order = (expConfigOptions.expansionCoeffFlag) ? 1 : 0;
  if (surrData.anchor()) {
    if (!surrData.anchor_gradient().empty())     data_order |= 2;
    //if (!surrData.anchor_hessian().empty())    data_order |= 4;
  }
  else {
    if (!surrData.response_gradient(0).empty())  data_order |= 2;
    //if (!surrData.response_hessian(0).empty()) data_order |= 4;
  }
  // verify support for basisConfigOptions.useDerivs, which indicates usage of
  // derivative data with respect to expansion variables (aleatory or combined)
  // within the expansion coefficient solution process, which must be
  // distinguished from usage of derivative data with respect to non-expansion
  // variables (the expansionCoeffGradFlag case).
  bool config_err = false;
  if (basisConfigOptions.useDerivs) {
    if (!(data_order & 2)) {
      PCerr << "Error: useDerivs configuration option lacks data support in "
	    << "RegressOrthogPolyApproximation::regression()" << std::endl;
      config_err = true;
    }
    if (expConfigOptions.expansionCoeffGradFlag) {
      PCerr << "Error: useDerivs configuration option conflicts with gradient "
	    << "expansion request in RegressOrthogPolyApproximation::"
	    << "regression()" << std::endl;
      config_err = true;
    }
    //if (data_order & 4)
    //  PCerr << "Warning: useDerivs configuration option does not yet support "
    //	      << "Hessian data in RegressOrthogPolyApproximation::regression()"
    //	      << std::endl;
  }
  if (config_err)
    abort_handler(-1);

  // compute data counts
  const SizetShortMap& failed_resp_data = surrData.failed_response_data();
  size_t num_failed_surr_fn = 0, num_failed_surr_grad = 0;
  SizetShortMap::const_iterator fit; bool faults_differ = false;
  for (fit=failed_resp_data.begin(); fit!=failed_resp_data.end(); ++fit) {
    short fail_bits = fit->second;
    if (fail_bits & 1) ++num_failed_surr_fn;
    if (fail_bits & 2) ++num_failed_surr_grad;
    // if failure omissions are not consistent, manage differing Psi matrices
    if ( (fail_bits & data_order) != data_order ) faults_differ = true;
  }
  num_surr_data_pts = surrData.size();
  num_data_pts_fn   = num_surr_data_pts - num_failed_surr_fn;
  num_data_pts_grad = num_surr_data_pts - num_failed_surr_grad;
  anchor_fn = false;
  anchor_grad = false;
  if (surrData.anchor()) {
    short failed_anchor_data = surrData.failed_anchor_data();
    if ((data_order & 1) && !(failed_anchor_data & 1)) anchor_fn   = true;
    if ((data_order & 2) && !(failed_anchor_data & 2)) anchor_grad = true;
  }

  // detect underdetermined system of equations (following fault omissions)
  // for either expansion coeffs or coeff grads (switch logic treats together)
  reuse_solver_data
    = (expConfigOptions.expansionCoeffFlag &&
       expConfigOptions.expansionCoeffGradFlag && !faults_differ);
  constr_eqns = 0;
  if (expConfigOptions.expansionCoeffFlag) {
    constr_eqns = num_data_pts_fn;
    if (anchor_fn)   constr_eqns += 1;
    if (anchor_grad) constr_eqns += numVars;
    total_eqns = (basisConfigOptions.useDerivs) ?
      constr_eqns + num_data_pts_grad*numVars : constr_eqns;
    if (total_eqns < numExpansionTerms) under_determined = true;
  }
  if (expConfigOptions.expansionCoeffGradFlag) {
    total_eqns = (anchor_grad) ? num_data_pts_grad+1 : num_data_pts_grad;
    if (total_eqns < numExpansionTerms) under_determined = true;
  }

  faultInfo.set_info( constr_eqns, anchor_fn, anchor_grad,
		      under_determined, num_data_pts_fn, num_data_pts_grad,
		      reuse_solver_data, total_eqns, num_surr_data_pts,
		      numVars, basisConfigOptions.useDerivs,
		      expansionCoeffGrads.numRows() );
};


void RegressOrthogPolyApproximation::run_regression()
{
  // Assume all function values are stored in top block of matrix in rows
  // 0 to num_surr_data_pts-1. Gradient information will be stored
  // in the bottom block of the matrix in rows num_surr_data_pts to
  // num_surr_data_pts + num_data_pts_grad * numVars. All the gradient 
  // information of point 0 will be stored consecutively then all the gradient
  // data of point 1, and so on.

  // Currently nothing is done  to modify the regression linear system matrices
  // A and B if surrData.anchor() is true, as currently surrData.anchor()
  // is always false. If in the future surrData.anchor() is enabled then
  // A must be adjusted to include the extra constraint information associated
  // with the anchor data. That is, if using EQ_CON_LEAST_SQUARES C matrix 
  // (top block of A ) must contain the fn and grad data of anchor point.
  // This will violate the first assumption discussed above and effect cross
  // validation. For this reason no modification is made as yet.

  size_t i, j, a_cntr = 0, b_cntr = 0, num_surr_data_pts = surrData.size(),
    num_deriv_vars = expansionCoeffGrads.numRows();
  int num_rows_A =  0, // number of rows in matrix A
      num_cols_A = numExpansionTerms, // number of columns in matrix A
      num_coeff_rhs, num_grad_rhs = num_deriv_vars, num_rhs;
  bool add_val, add_grad;

  RealMatrix A, B;
  
  bool multiple_rhs
    = (expConfigOptions.expansionCoeffFlag &&
       expConfigOptions.expansionCoeffGradFlag);

  bool anchor_fn = false, anchor_grad = false;
  if (surrData.anchor()) {
    anchor_fn = true;
    anchor_grad = true;
  }

  int num_data_pts_fn = num_surr_data_pts, 
    num_data_pts_grad = num_surr_data_pts;

  size_t a_grad_cntr = 0, b_grad_cntr = 0;

  CompressedSensingOptionsList opts_list;
  RealMatrixList solutions;
  CSOpts.standardizeInputs = false; // false essential when using derivatives
  
  if (expConfigOptions.expansionCoeffFlag) {

    // matrix/vector sizing
    num_rows_A = (basisConfigOptions.useDerivs) ?
      num_data_pts_fn + num_data_pts_grad * numVars : num_data_pts_fn;
    num_coeff_rhs = 1;
    num_rhs = (multiple_rhs) ? num_coeff_rhs + num_grad_rhs : num_coeff_rhs;
    PCout << "Applying regression to compute " << numExpansionTerms
	  << " chaos coefficients using " << num_rows_A << " equations.\n";

    A.shapeUninitialized(num_rows_A,num_cols_A);
    B.shapeUninitialized(num_rows_A,num_rhs);
    Real *A_matrix = A.values(), *b_vectors = B.values();
    // The "A" matrix is a contiguous block of memory packed in column-major
    // ordering as required by F77 for the GELSS subroutine from LAPACK.  For
    // example, the 6 elements of A(2,3) are stored in the order A(1,1),
    // A(2,1), A(1,2), A(2,2), A(1,3), A(2,3).
    for (i=0; i<numExpansionTerms; ++i) {
      a_cntr = num_rows_A*i;
      a_grad_cntr = a_cntr + num_data_pts_fn;
      const UShortArray& mi = multiIndex[i];
      for (j=0;j<num_surr_data_pts; ++j) {
	add_val = true; add_grad = basisConfigOptions.useDerivs;
	pack_polynomial_data(surrData.continuous_variables(j), mi, add_val,
			     A_matrix, a_cntr, add_grad, A_matrix, a_grad_cntr);
      }
    }
    
    // response data (values/gradients) define the multiple RHS which are
    // matched in the LS soln.  b_vectors is num_data_pts (rows) x num_rhs
    // (cols), arranged in column-major order.
    b_cntr = 0;
    b_grad_cntr = num_data_pts_fn;
    const SDRArray& sdr_array = surrData.response_data();
    for (i=0; i<num_surr_data_pts; ++i) {
      add_val = true; add_grad = basisConfigOptions.useDerivs;
      pack_response_data(sdr_array[i], add_val, b_vectors, b_cntr, add_grad,
			 b_vectors, b_grad_cntr);
    }

    // if no RHS augmentation, then solve for coeffs now
    if (!multiple_rhs) {

      // Perform cross validation loop over degrees here.
      // Current cross validation will not work for equality 
      // constrained least squares
      if ( crossValidation )
	// run cross validation
	run_cross_validation( A, B, num_data_pts_fn );
      else
	{
	  
	  IntVector index_mapping; 
	  remove_faulty_data( A, B, index_mapping,
			      faultInfo,
			      surrData.failed_response_data() );
	  CSTool.solve( A, B, solutions, CSOpts, opts_list );
	  
	  copy_data(solutions[0][0], numExpansionTerms, expansionCoeffs);
	}
    }
  }

  if (expConfigOptions.expansionCoeffGradFlag) {
    if (!multiple_rhs) {
      num_rows_A = num_data_pts_grad;
      num_rhs    = num_grad_rhs; num_coeff_rhs = 0;
      A.shapeUninitialized(num_rows_A,num_cols_A);
      B.shapeUninitialized(num_rows_A,num_rhs);
      Real *A_matrix   = A.values();

      // repack "A" matrix with different Psi omissions
      a_cntr = 0;
      for (i=0; i<numExpansionTerms; ++i) {
	const UShortArray& mi = multiIndex[i];
	for (j=0; j<num_surr_data_pts; ++j) {
	  add_val = false; add_grad = true;
	  if (add_grad) {
	    A_matrix[a_cntr] = multivariate_polynomial(
						       surrData.continuous_variables(j), mi);
	    ++a_cntr;
	  }
	}
      }
    }
    
    PCout << "Applying regression to compute gradients of " << numExpansionTerms
	  << " chaos coefficients using " << num_rows_A << " equations.\n";

    // response data (values/gradients) define the multiple RHS which are
    // matched in the LS soln.  b_vectors is num_data_pts (rows) x num_rhs
    // (cols), arranged in column-major order.
    Real *b_vectors  = B.values();
    b_cntr = 0;
    for (i=0; i<num_surr_data_pts; ++i) {
      add_val = false; add_grad = true;
      if (add_grad) {
	const RealVector& resp_grad = surrData.response_gradient(i);
	for (j=0; j<num_grad_rhs; ++j) // i-th point, j-th grad component
	  b_vectors[(j+num_coeff_rhs)*num_data_pts_grad+b_cntr] = resp_grad[j];
	++b_cntr;
      }
    }

    // solve
    IntVector index_mapping; 
    remove_faulty_data( A, B, index_mapping,
			faultInfo, surrData.failed_response_data());
    CSTool.solve( A, B, solutions, CSOpts, opts_list );
    
    if (multiple_rhs)
      {
	for ( int j = 0; j < numExpansionTerms; j++ )
	  expansionCoeffs[j] = solutions[0](j,0);
      }
    
    for (i=0; i<numExpansionTerms; ++i)
      for (j=0; j<num_grad_rhs; ++j)
	expansionCoeffGrads(j,i)
	  = solutions[j+num_coeff_rhs](i,0);

  }
}

void RegressOrthogPolyApproximation::
run_cross_validation( RealMatrix &A, RealMatrix &B, size_t num_data_pts_fn )
{
  RealMatrix A_copy( Teuchos::Copy, A, A.numRows(), A.numCols() );
  RealMatrix B_copy( Teuchos::Copy, B, B.numRows(), B.numCols() );
  int num_rhs = B.numCols(), num_dims( approxOrder.size() );
  // Do cross validation for varing polynomial orders up to 
  // a maximum order defined by approxOrder[0]
  int min_order = 1;
  if ( min_order > approxOrder[0] ) min_order = approxOrder[0];

  /// The options used to create the best PCE for each QOI
  std::vector<CompressedSensingOptions> bestCompressedSensingOpts_;

  /// The options of the best predictors of the predictors produced by each item
  /// in predictorOptionsList_. Information is stored for each PCE degree
  std::vector<RealMatrixList> predictorOptionsHistory_;

  /// The best predictors of the predictors produced by each item 
  /// in predictorOptionsList_. Information is stored for each PCE degree
  std::vector<RealMatrixList> predictorIndicatorsHistory_;

  /// The indicators of each partition for the best predictors of the 
  /// predictors produced by each item in predictorOptionsList_. 
  /// Information is stored for each PCE degree
  std::vector<RealMatrixList> predictorPartitionIndicatorsHistory_;
  bestCompressedSensingOpts_;

  std::vector<CompressedSensingOptions> best_cs_opts( num_rhs );
  
  RealVector min_best_predictor_indicators( num_rhs );
  min_best_predictor_indicators = std::numeric_limits<Real>::max();
  bestCompressedSensingOpts_.resize( num_rhs );
  IntVector best_cross_validation_orders( num_rhs );
  predictorOptionsHistory_.resize( approxOrder[0] - min_order + 1 );
  predictorIndicatorsHistory_.resize( approxOrder[0] - min_order + 1 );
  predictorPartitionIndicatorsHistory_.resize( approxOrder[0] - min_order + 1 );
  int cnt( 0 );
  for ( int order = min_order; order <= approxOrder[0]; order++ )
    {
      if (expConfigOptions.outputLevel >= NORMAL_OUTPUT)
	PCout << "Testing PCE order " << order << std::endl;
      int num_basis_terms = nchoosek( num_dims + order, order );
      RealMatrix vandermonde_submatrix( Teuchos::View, 
					A_copy,
					A_copy.numRows(),
					num_basis_terms, 0, 0 );

      RealVector best_predictor_indicators;
      estimate_compressed_sensing_options_via_cross_validation( 
							       vandermonde_submatrix, 
							       B_copy, 
							       best_cs_opts,
							       best_predictor_indicators,
							       predictorOptionsHistory_[cnt], 
							       predictorIndicatorsHistory_[cnt],  
							       predictorPartitionIndicatorsHistory_[cnt],
							       num_data_pts_fn );

      // Only execute on master processor
      //      if ( is_master() )
      if ( true )
	{
	  for ( int k = 0; k < num_rhs; k++ )
	    {
	      if ( best_predictor_indicators[k] < 
		   min_best_predictor_indicators[k] )
		{
		  min_best_predictor_indicators[k] = 
		    best_predictor_indicators[k];
		  best_cross_validation_orders[k] = order;
		  bestCompressedSensingOpts_[k] = best_cs_opts[k];
		}
	      if (expConfigOptions.outputLevel >= NORMAL_OUTPUT)
		{
		  PCout<<"Cross validation error for rhs "<<k<<" and degree ";
		  PCout << order << ": " <<  best_predictor_indicators[k]<< "\n";
		}
	    }
	}
      cnt++;
    }
  bestApproxOrder = best_cross_validation_orders;
  int num_basis_terms = nchoosek( num_dims + bestApproxOrder[0], 
				  bestApproxOrder[0] );
  if (expConfigOptions.outputLevel >= NORMAL_OUTPUT)
    {
      PCout << "Best approximation order: " << bestApproxOrder[0]<< "\n";
      PCout << "Best cross validation error: ";
      PCout <<  min_best_predictor_indicators[0]<< "\n";
    }
  // set CSOpts so that best PCE can be built. We are assuming num_rhs=1
  RealMatrix vandermonde_submatrix( Teuchos::View, 
				    A_copy,
				    A_copy.numRows(),
				    num_basis_terms, 0, 0 );
  CompressedSensingOptionsList opts_list;
  RealMatrixList solutions;
  bestCompressedSensingOpts_[0].storeHistory = false;
  bestCompressedSensingOpts_[0].print();
  IntVector index_mapping;
  remove_faulty_data( vandermonde_submatrix, B_copy, index_mapping,
		      faultInfo, surrData.failed_response_data() );
  CSTool.solve( vandermonde_submatrix, B_copy, solutions, 
		bestCompressedSensingOpts_[0], opts_list );

  // Resize solutions so that it can be used with large vandermonde.
  if (expansionCoeffs.length()!=numExpansionTerms)
    expansionCoeffs.sizeUninitialized(numExpansionTerms);
  for ( int i=0; i<num_basis_terms; ++i)
    expansionCoeffs[i] = solutions[0](i,0);
  for ( int i=num_basis_terms; i < numExpansionTerms; ++i)
    expansionCoeffs[i] = 0.0;
};

void RegressOrthogPolyApproximation::gridSearchFunction( RealMatrix &opts,
						  int M, int N, 
						  int num_function_samples )
{
  // Setup a grid based search
  bool is_under_determined = M < N;
  
  // Define the 1D grids for under and over-determined LARS, LASSO, OMP, BP and 
  // LS
  std::vector<RealVector> opts1D( 9 );
  opts1D[0].size( 1 ); // solver type
  opts1D[0][0] = CSOpts.solver;
  opts1D[1].size( 1 ); // Solver tolerance. 
  opts1D[1][0] = CSOpts.solverTolerance;
  opts1D[2] = noiseTols; // epsilon.
  opts1D[3].size( 1 ); // delta
  opts1D[3] = CSOpts.delta;
  opts1D[4].size( 1 ); // max_number of non_zeros
  opts1D[4] = CSOpts.maxNumIterations; 
  opts1D[5].size( 1 );  // standardizeInputs
  opts1D[5] = false;
  opts1D[6].size( 1 );  // storeHistory
  opts1D[6] = true;  
  opts1D[7].size( 1 );  // Verbosity. Warnings on
  opts1D[7] =  std::max(0, expConfigOptions.outputLevel - 1);
  opts1D[8].size( 1 );  // num function samples
  opts1D[8] = num_function_samples;
      
  // Form the multi-dimensional grid
  cartesian_product( opts1D, opts );
};

void RegressOrthogPolyApproximation::
estimate_compressed_sensing_options_via_cross_validation( RealMatrix &vandermonde_matrix, RealMatrix &rhs, std::vector<CompressedSensingOptions> &best_cs_opts, RealVector &best_predictor_indicators, RealMatrixList &predictor_options_history, RealMatrixList &predictor_indicators_history, RealMatrixList &predictor_partition_indicators_history, size_t num_data_pts_fn ){
  // Initialise the cross validation iterator
  CrossValidationIterator CV;
  CV.mpi_communicator( MPI_COMM_WORLD );
  CV.verbosity(  std::max(0, expConfigOptions.outputLevel - 1) );
  // Set and partition the data
  CV.set_data( vandermonde_matrix, rhs, num_data_pts_fn );
  int num_folds( 10 );
  // Keep copy of state
  CompressedSensingOptions cs_opts_copy = CSOpts;
  
  if ( ( ( num_data_pts_fn / num_folds == 1 ) && 
	 ( num_data_pts_fn - 3 < vandermonde_matrix.numCols() ) )  || 
       ( num_data_pts_fn / num_folds == 0 ) )
    // use one at a time cross validation
    num_folds = num_data_pts_fn;
  if ( num_data_pts_fn == vandermonde_matrix.numCols() )
    {
      PCout << "Warning: The cross validation results will not be consistent. ";
      PCout << "The number of function samples = the number of basis terms, ";
      PCout << "thus only underdetermined matrices will be generated during ";
      PCout << "cross validation even though the system is fully determined.\n";
    }

  if ( ( CSOpts.solver == EQ_CON_LEAST_SQ_REGRESSION ) &&
       ( num_folds = num_data_pts_fn ) && 
       ( vandermonde_matrix.numRows() * ( num_data_pts_fn - 1 ) / num_data_pts_fn  <= vandermonde_matrix.numCols() ) )
    // includes exactly determined case
    {
      PCout << "EQ_CON_LEAST_SQ_REGRESSION could not be used. ";
      PCout << "The cross validation training vandermonde matrix is ";
      PCout << "under-determined\n";
      CSOpts.solver = LASSO_REGRESSION;
    }
  if ( ( CSOpts.solver == EQ_CON_LEAST_SQ_REGRESSION ) && ( num_folds = num_data_pts_fn ) &&  ( num_data_pts_fn - 1 >= vandermonde_matrix.numCols() ) )
    {
      PCout << "EQ_CON_LEAST_SQ_REGRESSION could not be used. ";
      PCout << "The cross validation training vandermonde matrix is ";
      PCout << "over-determined\n";
      CSOpts.solver = DEFAULT_REGRESSION;
    }
    
  CV.setup_k_folds( num_folds );
  
  // Tell the cross validation iterator what options to test
  RealMatrix opts;
  gridSearchFunction( opts, vandermonde_matrix.numRows(),
		      vandermonde_matrix.numCols(), num_data_pts_fn );
  CV.predictor_options_list( opts );

  // Perform cross validation
  CV.run( &rmse_indicator, &linear_predictor_analyser, 
	  &normalised_mean_selector,
	  &linear_predictor_best_options_extractor,
	  faultInfo, surrData.failed_response_data() );

  // Get results of cross validation
  RealMatrix best_predictor_options;
  CV.get_best_predictor_info( best_predictor_options, 
			      best_predictor_indicators );

  CV.get_history_data( predictor_options_history, 
		       predictor_indicators_history,
		       predictor_partition_indicators_history );

  //if ( CV.is_master() )
  if ( true )
    {
      int len_opts(  best_predictor_options.numRows() ), 
	num_rhs( rhs.numCols() );
      for ( int k = 0; k < num_rhs; k++ )
	{
	  RealVector col( Teuchos::View,  best_predictor_options[k], len_opts );
	  set_linear_predictor_options( col, best_cs_opts[k] );
	}
    }

  //restore state
  CSOpts = cs_opts_copy;
};

} // namespace Pecos

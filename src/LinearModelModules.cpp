#include "LinearModelModules.hpp"
#include "MathTools.hpp"

namespace Pecos {

void set_linear_predictor_options( RealVector &predictor_opts,  
				   CompressedSensingOptions &cs_opts )
{
  cs_opts.solver = predictor_opts[0];
  cs_opts.solverTolerance = predictor_opts[1];
  cs_opts.epsilon = predictor_opts[2];
  cs_opts.delta = predictor_opts[3];
  cs_opts.maxNumIterations = (int)predictor_opts[4];
  cs_opts.standardizeInputs = (bool)predictor_opts[5];
  cs_opts.storeHistory = (bool)predictor_opts[6];
  cs_opts.verbosity = (int)predictor_opts[7];
  cs_opts.numFunctionSamples = (int)predictor_opts[8];
};

void extract_linear_predictor_options(CompressedSensingOptionsList &cs_opts_list,
				      RealMatrixList &predictor_options_list )
{
  int num_qoi( cs_opts_list.size() );
  predictor_options_list.resize( num_qoi );
  for ( int k = 0; k < num_qoi; k++ )
    {
      int num_predictors( cs_opts_list[k].size() );
      predictor_options_list[k].shapeUninitialized( 9, num_predictors );
      for ( int j = 0; j < num_predictors; j++ )
	{
	  predictor_options_list[k](0,j) = cs_opts_list[k][j].solver;
	  predictor_options_list[k](1,j) = cs_opts_list[k][j].solverTolerance;
	  predictor_options_list[k](2,j) = cs_opts_list[k][j].epsilon;
	  predictor_options_list[k](3,j) = cs_opts_list[k][j].delta;
	  predictor_options_list[k](4,j) = cs_opts_list[k][j].maxNumIterations;
	  predictor_options_list[k](5,j) = cs_opts_list[k][j].standardizeInputs;
	  predictor_options_list[k](6,j) = cs_opts_list[k][j].storeHistory;
	  predictor_options_list[k](7,j) = cs_opts_list[k][j].verbosity;
	  predictor_options_list[k](8,j) = cs_opts_list[k][j].numFunctionSamples;
	}
    }
};

void linear_predictor_analyser( RealMatrix &A_training, 
				RealMatrix &B_training, 
				RealMatrix &A_validation,
				RealMatrix &B_validation,
				RealVector &predictor_opts,
				IndicatorFunction *indicator_function,
				RealMatrixList &indicators_list,
				RealMatrixList &predictor_options_list,
				FaultInfo &fault_info,
				const SizetShortMap& failed_resp_data,
				IntVector &training_indices,
				IntVector &validation_indices )
{
 
  // Extract the predictor options and store in the format needed
  CompressedSensingOptions cs_opts;
  set_linear_predictor_options( predictor_opts, cs_opts );

  // Compute a set of coefficients that solve the linear system
  // Ax = b, where A = training_samples and b = training_values;
  CompressedSensingTool cs_tool;
  RealMatrixList coefficient_sets;
  CompressedSensingOptionsList cs_opts_list;
  remove_faulty_data( A_training, B_training, training_indices,
  		      fault_info, failed_resp_data );

  remove_faulty_data( A_validation, B_validation, validation_indices,
		      fault_info, failed_resp_data );
  
  cs_tool.solve( A_training, B_training, coefficient_sets,
		 cs_opts, cs_opts_list );

  // Convert cs_opts_list into predictor_options_list
  extract_linear_predictor_options( cs_opts_list, predictor_options_list );

  // Evaluate the predictor on the validation data
  int num_qoi( coefficient_sets.size() );
  int num_validation_samples( B_validation.numRows() );

  indicators_list.resize( num_qoi );
  for ( int k = 0; k < num_qoi; k++ )
    {
      RealMatrix B_prediction( A_validation.numRows(), 
			       coefficient_sets[k].numCols() );
      // Use the linear model to predict at the validation points
      B_prediction.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0,
			     A_validation, coefficient_sets[k], 0.0 );

      // Get view of validation data
      RealMatrix b_validation( Teuchos::View, B_validation, 
			       num_validation_samples, 1,
			       0, k );
      
      // The number of coefficients sets can be different for each QOI
      int num_coefficient_sets( coefficient_sets[k].numCols() );
      indicators_list[k].shapeUninitialized( num_coefficient_sets, 1 );
      for ( int j = 0; j < num_coefficient_sets; j++ )
	{
	  // Get view of prediction data for the jth predictor created
	  RealMatrix b_prediction( Teuchos::View, B_prediction, 
				   num_validation_samples, 1,
				   0, j );
	  // Get view of indicator elements to avoid unnecessary copying
	  RealMatrix indicators( Teuchos::View, indicators_list[k], 
				 1, 1, j, 0 );
	  indicator_function( b_validation, b_prediction, indicators );
	}
    }
};

void linear_predictor_best_options_extractor( std::vector<RealMatrixList> &partition_options, IntVector &best_predictor_indices, int num_training_samples, int num_samples, RealMatrix &best_predictor_options )
{
  int num_partitions = partition_options.size(), 
    num_qoi = best_predictor_indices.length();
  
  for ( int k = 0; k < num_qoi; k++ )
    {
      std::vector<Real> epsilons( num_partitions );

      // Get the average number of max_iterations and the average residual
      // for the best predictors
      Real ave_epsilon( 0.0 );
      int best_max_num_iterations( 0 );
      int argmin_k = best_predictor_indices[k];
      for ( int i = 0; i < num_partitions; i++ )
	{
	  epsilons[i] = partition_options[i][k](2,argmin_k);
	  ave_epsilon += epsilons[i];
	  //best_max_num_iterations += partition_options[i][k](4,argmin_k);
	  best_max_num_iterations += argmin_k;
	}
      double sample_size_ratio = (Real)num_samples / 
	(Real)num_training_samples;
      ave_epsilon /= (Real)num_partitions;
      double median_epsilon = median( epsilons );
      // using ave_epsilon does not work well for rosenbrock
      // so set to best_epsilon to median_epsilon
      double best_epsilon = median_epsilon;
      // sqrt used because epsilon is the sqrt( r'r );
      best_epsilon *= std::sqrt( sample_size_ratio );
      best_max_num_iterations /= num_partitions;
      best_max_num_iterations = std::ceil( (Real)best_max_num_iterations *
					   sample_size_ratio )+1;
      // Copy the partition_options of the best predictor for the
      // first parition. Then adjust the values that can vary with the
      // partition with their averages
      RealVector best_opts_k( Teuchos::View, 
			      partition_options[0][k][argmin_k],
			      partition_options[0][k].numRows() );
      append_column( best_opts_k, best_predictor_options );
      best_predictor_options(2,k) = best_epsilon;
      best_predictor_options(4,k) = best_max_num_iterations;
    }
}

} //namespace Pecos

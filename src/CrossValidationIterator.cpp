#include "CrossValidationIterator.hpp"
#include "MathTools.hpp"
 
namespace Pecos {


void CrossValidationIterator::setup_k_folds( int num_partitions )
{
  numPartitions_ = num_partitions;
  
  k_folds_check_inputs();
    
  // Randomly permute the row indices of the training points and values.
  // This is used to partition the training data into k sets.
  IntMatrix permutations;
  get_permutations( permutations, numFunctionSamples_, 1, seed_ );
  
  // Generate the k sets for cross validation
  int num_samples_per_fold = numFunctionSamples_ / 
    numPartitions_;

  int num_small_folds = num_partitions * ( num_samples_per_fold + 1 ) - 
    numFunctionSamples_;
  int num_large_folds = numFunctionSamples_ - num_partitions * num_samples_per_fold;

  int max_num_training_samples = num_large_folds * ( num_samples_per_fold + 1 )+
    (num_small_folds - 1) * num_samples_per_fold;
  int max_num_validation_samples = num_samples_per_fold + 1;

  RealMatrix training_indices, validation_indices;
  training_indices.reshape( max_num_training_samples, numPartitions_ );
  validation_indices.reshape( max_num_validation_samples, numPartitions_ );

  RealVector  training_indices_sizes( numPartitions_ ),
    validation_indices_sizes( numPartitions_ );

  int partition_id( 0 );
  // Create the first partitions with the validation sets containing a large
  // fold
  int num_validation_samples = num_samples_per_fold + 1;
  for ( int j = 0; j < num_large_folds; j++ )
    {

      // Store validation indices
      for ( int i = 0; i < num_validation_samples ; i++ )
	{
	  validation_indices(i,partition_id) = 
	    permutations(j*num_validation_samples+i,0);
	}
      validation_indices_sizes[partition_id] = num_validation_samples;

      // Store training indices
      for ( int i = 0; i < numFunctionSamples_ - num_validation_samples ; i++ )
	{
	  training_indices(i,partition_id) = 
	    permutations(((j+1)*num_validation_samples+i)%numFunctionSamples_,0);
	}
      training_indices(max_num_training_samples-1,partition_id) = -1;
      training_indices_sizes[partition_id] = max_num_training_samples-1;
      partition_id++;
    }

  // The remaining partitions will consist of validation sets containing a small
  // fold
  num_validation_samples = num_samples_per_fold;
  int k = num_large_folds * ( num_samples_per_fold + 1 );
  for ( int j = 0; j < num_small_folds; j++ )
    {

      // Store validation indices
      for ( int i = 0; i < num_validation_samples ; i++ )
	{
	  validation_indices(i,partition_id) = 
	    permutations(k + j*num_validation_samples+i,0);
	}
      validation_indices(max_num_validation_samples-1,partition_id) = -1;
      validation_indices_sizes[partition_id] = num_validation_samples;

      // Store training indices
      for ( int i = 0; i < numFunctionSamples_ - num_validation_samples ; i++ )
	{
	  training_indices(i,partition_id) = 
	    permutations((k+(j+1)*num_validation_samples+i)%numFunctionSamples_,0);
	}
      training_indices_sizes[partition_id] = max_num_training_samples;
      partition_id++;
    }

  set_partition_indices( training_indices, validation_indices,
			 training_indices_sizes, validation_indices_sizes );
};

void CrossValidationIterator::partition_data( RealMatrix &training_samples, 
					      RealMatrix &training_values, 
					      RealMatrix &validation_samples,
					      RealMatrix &validation_values,
					      int fold )
{
  int num_training_samples( trainingIndicesSizes_[fold] ),
    num_covariates( samples_.numCols() ), num_qoi( values_.numCols() ),
    num_validation_samples( validationIndicesSizes_[fold] );

  int num_data_per_function_sample = samples_.numRows()/numFunctionSamples_;
 
  reshape( training_samples, 
	   num_data_per_function_sample * num_training_samples, 
	   num_covariates );
  reshape( training_values, 
	   num_data_per_function_sample * num_training_samples, 
	   num_qoi );
  reshape( validation_samples, 
	   num_data_per_function_sample * num_validation_samples, 
	   num_covariates );
  reshape( validation_values, 
	   num_data_per_function_sample * num_validation_samples, 
	   num_qoi );

  // Partition the sample data
  for ( int n = 0; n < num_covariates; n++ )
    {
      // Partition training samples
      for ( int i = 0; i < num_training_samples ; i++ )
	{
	  if ( trainingIndices_(i,fold) < 0 )
	    break;
	  for ( int j = 0; j < num_data_per_function_sample; j++ )
	    {
	      training_samples(j*num_training_samples+i,n) = 
		samples_(j*numFunctionSamples_+trainingIndices_(i,fold),n);
	    }
	}
      // Partition validation samples
      for ( int i = 0; i < num_validation_samples ; i++ )
	{
	  if ( validationIndices_(i,fold) < 0 )
	    break;
	  for ( int j = 0; j < num_data_per_function_sample; j++ )
	    {
	      validation_samples(j*num_validation_samples+i,n) = 
		samples_(j*numFunctionSamples_+validationIndices_(i,fold),n);
	    }
	}
    }

  // Partition the values data
  for ( int k = 0; k < num_qoi; k++ )
    {
      // Partition training function values
      for ( int i = 0; i < num_training_samples ; i++ )
	{
	  if ( trainingIndices_(i,fold) < 0 )
	    break;
	  for ( int j = 0; j < num_data_per_function_sample; j++ )
	    {
	      training_values(j*num_training_samples+i,k) = 
		values_(j*numFunctionSamples_+trainingIndices_(i,fold),k);
	    }
	}
      // Partition validation function values
      for ( int i = 0; i < num_validation_samples ; i++ )
	{
	  if ( validationIndices_(i,fold) < 0 )
	    break;
	  for ( int j = 0; j < num_data_per_function_sample; j++ )
	    {
	      validation_values(j*num_validation_samples+i,k) = 
	      values_(j*numFunctionSamples_+validationIndices_(i,fold),k);
	    }
	}
    }
};

/** \brief Run cross validation
 *
 * \param indicator_function The cross validation error indicator
 *
 * \param analyser Function that builds a predictor with one set of predictor
 * options and uses indicator_function to determine the quantify the 
 * cross validation accuracy of the predictor that was built
 */
void CrossValidationIterator::run( IndicatorFunction *indicator_function, 
				   Analyser *analyser, Selector *selector,
				   BestOptionsExtractor *best_options_extractor )
{
  if ( verbosity_ > 0 )
    {
      std::cout << "Processor " << processorId_ << " starting cross ";
      std::cout << "validation run\n";
    }
  RealMatrix training_samples, training_values, validation_samples,
    validation_values;
  int num_predictor_options_list( predictorOptionsList_.numCols() );
  for ( int opts_iter = 0; opts_iter < num_predictor_options_list; opts_iter++ )
    {
      RealVector opts( Teuchos::View, predictorOptionsList_[opts_iter], 
		       predictorOptionsList_.numRows() );
      std::vector<RealMatrixList> all_partition_indicators( numPartitions_ );
      std::vector<RealMatrixList> all_partition_options( numPartitions_ );
      RealMatrixList predictor_options_list;
      for ( int j = 0; j < numPartitions_; j++ )
	{ 
	  partition_data( training_samples, training_values,
			  validation_samples, validation_values, j );
	  RealMatrixList indicators_list;
	  analyser( training_samples, training_values,
		    validation_samples, validation_values, 
		    opts,
		    indicator_function,
		    indicators_list,
		    predictor_options_list );
	  all_partition_indicators[j] = indicators_list;
	  all_partition_options[j] = predictor_options_list;
	}

      // Determine  the "best" predictor options
      IntVector best_predictor_indices;
      RealVector best_predictor_indicators;
      RealMatrix best_predictor_partition_indicators;
      selector( all_partition_indicators,
		best_predictor_indices,
		best_predictor_indicators,
		best_predictor_partition_indicators );

      RealMatrix best_predictor_options;
      best_options_extractor( all_partition_options, best_predictor_indices,
			      training_samples.numRows(), samples_.numRows(),
			      best_predictor_options );

      update_best_predictor_info( best_predictor_indicators,
				  best_predictor_options,
				  best_predictor_partition_indicators  );

    }
  if ( verbosity_ > 0 )
    {
      std::cout << "Processor " << processorId_ << " finished cross ";
      std::cout << "validation run\n";
    }
};

void CrossValidationIterator::update_best_predictor_info( RealVector &predictor_indicators, RealMatrix &predictor_options, RealMatrix &predictor_partition_indicators)
{
  // Update history
  predictorOptionsHistory_.push_back( predictor_options );
  predictorIndicatorsHistory_.push_back( predictor_indicators );
  predictorPartitionIndicatorsHistory_.push_back( predictor_partition_indicators );
  // Update best predictor info
  int num_qoi(  predictor_options.numCols() ), 
    len_options( predictor_options.numRows() );
  
  if ( bestPredictorIndicators_.numRows() == 0 )
    {
      bestPredictorIndicators_.sizeUninitialized( num_qoi );
      bestPredictorIndicators_ = std::numeric_limits<Real>::max();
      bestPredictorOptions_.shapeUninitialized( len_options, num_qoi );
    }
  for ( int k = 0; k < num_qoi; k++ )
    {
      if (  predictor_indicators[k] < bestPredictorIndicators_[k] )
	{
	  bestPredictorIndicators_[k] =  predictor_indicators[k];
	  RealVector opts_k( Teuchos::View, predictor_options[k],
				  predictor_options.numRows() );
	  fill_column( k, opts_k, bestPredictorOptions_ );
	}
    }
};

void CrossValidationIterator::get_best_predictor_info( RealMatrix &best_predictor_options, RealVector &best_predictor_indicators )
{
#ifdef ENABLE_LIBHEAT_MPI
  if ( is_master() )
    {
      int num_qoi(  bestPredictorOptions_.numCols() ),
	len_options( bestPredictorOptions_.numRows() );
      best_predictor_options = bestPredictorOptions_;
      best_predictor_indicators = bestPredictorIndicators_;

      RealMatrix slave_best_predictor_options;
      RealVector slave_best_predictor_indicators;
      for ( int proc_id = 1; proc_id < num_processors(); proc_id++ )
	{
	  receive( slave_best_predictor_options, proc_id, MPICommunicator_ );
	  receive( slave_best_predictor_indicators, proc_id, MPICommunicator_ );

	  if ( verbosity_ > 0 )
	    {
	      std::cout << "Master received best predictor data from slave ";
	      std::cout << proc_id << "\n";
	    }

	  if ( ( proc_id == 1 ) && ( predictorOptionsList_.numCols() == 0 ) )
	    // No work was assigned to master. So use the results from the 
	    // first slave process to define inital values of 
	    // best_predictor_options and best_predictor_indicators
	    {
	      best_predictor_options = slave_best_predictor_options;
	      best_predictor_indicators = slave_best_predictor_indicators;
	      num_qoi = best_predictor_options.numCols();
	      len_options = best_predictor_options.numRows();
	    }
	  else
	    // Update best predictor info
	    {
	      if ( slave_best_predictor_indicators.length() > 0 )
		{
		  for ( int k = 0; k < num_qoi; k++ )
		    {
		      if (  slave_best_predictor_indicators[k] < 
			    best_predictor_indicators[k] )
			{
			  best_predictor_indicators[k] = 
			    slave_best_predictor_indicators[k];
			  RealVector slave_best_opts_k( Teuchos::View, 
						slave_best_predictor_options[k],
							len_options );
			  fill_column( k, slave_best_opts_k, 
				       best_predictor_options );
			}
		    }
		}
	    }
	}
    }
  else
    {
      send( bestPredictorOptions_, master_processor_id(), MPICommunicator_ );
      send( bestPredictorIndicators_, master_processor_id(), MPICommunicator_ );

      if ( verbosity_ > 0 )
	{
	  std::cout << "Slave " << processorId_ << " sent best predictor data ";
	  std::cout << "to master\n";
	}
    };

#else

  best_predictor_options = bestPredictorOptions_;
  best_predictor_indicators =  bestPredictorIndicators_;

#endif
};

void CrossValidationIterator::get_history_data( RealMatrixList &predictor_options_history, RealMatrixList &predictor_indicators_history, RealMatrixList &predictor_partition_indicators_history )
{
#ifdef ENABLE_LIBHEAT_MPI
  if ( is_master() )
    {
      predictor_options_history.clear();
      predictor_indicators_history.clear();
      predictor_partition_indicators_history.clear();

      RealMatrixList slave_predictor_options_history,
	slave_predictor_indicators_history,
	slave_predictor_partition_indicators_history;
      for ( int proc_id = 1; proc_id < num_processors(); proc_id++ )
	{
	  receive( slave_predictor_options_history, proc_id, MPICommunicator_ );
	  receive( slave_predictor_indicators_history, proc_id, 
		   MPICommunicator_ );
	  receive( slave_predictor_partition_indicators_history, proc_id, 
		   MPICommunicator_ );
	  
	  if ( verbosity_ > 0 )
	    {
	      std::cout << "Master received predictor history from slave ";
	      std::cout << processor_id() << "\n";
	    }
	  predictor_options_history.insert( predictor_options_history.end(), 
					slave_predictor_options_history.begin(),
					slave_predictor_options_history.end() );
	  predictor_indicators_history.insert(predictor_indicators_history.end(),
				      slave_predictor_indicators_history.begin(),
				      slave_predictor_indicators_history.end() );
	  predictor_partition_indicators_history.insert(
			   predictor_partition_indicators_history.end(),
			   slave_predictor_partition_indicators_history.begin(),
	  		   slave_predictor_partition_indicators_history.end() );
	  
	}
      
      predictor_options_history.insert( predictor_options_history.end(),
					predictorOptionsHistory_.begin(),
					predictorOptionsHistory_.end() );
      predictor_indicators_history.insert( predictor_indicators_history.end(),
					   predictorIndicatorsHistory_.begin(),
					   predictorIndicatorsHistory_.end() );
      predictor_partition_indicators_history.insert( 
				  predictor_partition_indicators_history.end(),
				  predictorPartitionIndicatorsHistory_.begin(),
				  predictorPartitionIndicatorsHistory_.end() );
     
      reshape_history_data( predictor_options_history, 
			    predictor_indicators_history,
			    predictor_partition_indicators_history );
    }
  else
    {
      send( predictorOptionsHistory_, master_processor_id(), MPICommunicator_ );
      send( predictorIndicatorsHistory_, master_processor_id(), MPICommunicator_ );
      send( predictorPartitionIndicatorsHistory_, master_processor_id(), MPICommunicator_ );

      if ( verbosity_ > 0 )
	{
	  std::cout << "Slave " << processorId_ << " sent best predictor data ";
	  std::cout << "to master\n";
	}
    }

#else

      predictor_options_history = predictorOptionsHistory_;
      predictor_indicators_history = predictorIndicatorsHistory_;
      predictor_partition_indicators_history = 
	predictorPartitionIndicatorsHistory_;

      reshape_history_data( predictor_options_history, 
			    predictor_indicators_history,
			    predictor_partition_indicators_history );

#endif
};

void CrossValidationIterator::reshape_history_data( RealMatrixList &predictor_options_history, RealMatrixList &predictor_indicators_history, RealMatrixList &predictor_partition_indicators_history )
{
  int num_opts( predictor_partition_indicators_history.size() );
  int num_rhs( predictor_partition_indicators_history[0].numCols() );
  RealMatrixList partition_indicators_history( num_rhs );
  RealMatrixList indicators_history( num_rhs );
  RealMatrixList options_history( num_rhs );
  for ( int k = 0; k < num_rhs; k++ ) 
    {
      partition_indicators_history[k].shapeUninitialized( numPartitions_,
							  num_opts );
      indicators_history[k].shapeUninitialized( num_opts, 1 );
      int len_opts( predictor_options_history[0].numRows() );
      options_history[k].shapeUninitialized( len_opts, num_opts );
      for ( int i = 0; i < num_opts; i++)
	{
	  indicators_history[k](i,0) = predictor_indicators_history[i](k,0);
	  for ( int j = 0; j < numPartitions_; j++ ) 
	    {
	      partition_indicators_history[k](j,i) = 
		predictor_partition_indicators_history[i](j,k);
	    }
	  for ( int j = 0; j < len_opts; j++ ) 
	    {
	      options_history[k](j,i) =  predictor_options_history[i](j,k); 
	    };
	}
    }

  predictor_options_history = options_history;
  predictor_indicators_history = indicators_history;
  predictor_partition_indicators_history = partition_indicators_history;

};

void CrossValidationIterator::set_data( RealMatrix &samples, RealMatrix &values,
					int num_function_samples )
{
  // pass in and set useGradients_ and numFunctionSamples_ here

#ifdef ENABLE_LIBHEAT_MPI

  if ( is_master() )
    {
      for ( int proc_id = 1; proc_id < num_processors(); proc_id++ )
	{
	  send( samples, proc_id, MPICommunicator_ );
	  send( values, proc_id, MPICommunicator_ );
	  if ( verbosity_ > 0 )
	    {
	      std::cout << "Master sent data to slave ";
	      std::cout << proc_id << "\n";
	    }
	}
      samples_ = samples;
      values_ = values;
      if ( verbosity_ > 0 )
	std::cout << "Master has set data\n";
    }
  else
    {
      receive( samples_, master_processor_id(), MPICommunicator_ );
      receive( values_, master_processor_id(), MPICommunicator_ );

      if ( verbosity_ > 0 )
	std::cout << "Slave " << processorId_ << " has set data\n";
    };

#else

  samples_ = samples;
  values_ = values;

#endif

  if ( num_function_samples == 0 )
    {
      numFunctionSamples_ = samples_.numRows();
      useGradients_ = false;
    }
  else
    {
      if ( ( num_function_samples != samples_.numRows() ) && 
	   ( samples_.numRows() % num_function_samples != 0 ) )
	{
	  std::stringstream msg;
	  msg << "num_function_samples is inconistent with the data.\n";
	  msg << "If data contains gradients, set num_function_samples = ";
	  msg << "num_dims * num_total_samples / ( num_dims + 1 )\n";
	  msg << "If data does not contain gradients, set ";
	  msg << "num_function_samples = num_total_samples;";
	  throw( std::runtime_error( msg.str() ) );
	}
      numFunctionSamples_ = num_function_samples;
    }
  if ( numFunctionSamples_ != samples_.numRows() )
    useGradients_= true;
  
};

void CrossValidationIterator::set_partition_indices( RealMatrix &training_indices, RealMatrix &validation_indices, RealVector &training_indices_sizes, RealVector &validation_indices_sizes )
{
  
#ifdef ENABLE_LIBHEAT_MPI

  if ( is_master() )
    {
      for ( int proc_id = 1; proc_id < num_processors(); proc_id++ )
	{
	  send( training_indices, proc_id, MPICommunicator_ );
	  send( validation_indices, proc_id, MPICommunicator_ );
	  send( training_indices_sizes, proc_id, MPICommunicator_ );
	  send( validation_indices_sizes, proc_id, MPICommunicator_ );
	  if ( verbosity_ > 0 )
	    {
	      std::cout << "Master sent partition indices to slave ";
	      std::cout << proc_id << "\n";
	    }
	}
      trainingIndices_ = training_indices;
      validationIndices_ = validation_indices;
      trainingIndicesSizes_ = training_indices_sizes;
      validationIndicesSizes_ = validation_indices_sizes;
      if ( verbosity_ > 0 )
	std::cout << "Master has set partition indices\n";
    }
  else
    {
      receive( trainingIndices_, master_processor_id(), MPICommunicator_ );
      receive( validationIndices_, master_processor_id(), MPICommunicator_ );
      receive( trainingIndicesSizes_, master_processor_id(), MPICommunicator_ );
      receive( validationIndicesSizes_, master_processor_id(), MPICommunicator_);

      if ( verbosity_ > 0 )
	std::cout << "Slave " << processorId_ << " has set partition indices\n";
    };
  
#else
  
  trainingIndices_ = training_indices;
  validationIndices_ = validation_indices;
  trainingIndicesSizes_ = training_indices_sizes;
  validationIndicesSizes_ = validation_indices_sizes;
  
#endif 
};

void CrossValidationIterator::predictor_options_list( RealMatrix &predictor_options_list_in )
{

#ifdef ENABLE_LIBHEAT_MPI
  
  if ( is_master() )
    {
      int len_options( predictor_options_list_in.numRows() );
      int num_options(  predictor_options_list_in.numCols() );
      // Determine the number of options to pass to each processor
      int load( std::max( 1, num_options / num_processors() ) );
      for ( int proc_id = 1; proc_id < num_processors(); proc_id++ )
	{
	  int proc_load = load, start_col( (proc_id-1)*load ),
	    num_rows( len_options );
	  if ( start_col >= num_options )
	    {
	      proc_load = 0;
	      start_col = 0;
	      num_rows = 0;
	    }
	  RealMatrix predictor_options_subset( Teuchos::View, 
					       predictor_options_list_in,
					       num_rows, proc_load, 0, 
					       start_col );
	  send( predictor_options_subset, proc_id, MPICommunicator_ );
	  if ( verbosity_ > 0 )
	    {
	      std::cout << "Master sent predictor options to slave ";
	      std::cout << proc_id << "\n";
	    }
	}

      // Assign the remaining options to the master
      int master_load( std::max(0,num_options - (num_processors()-1) * load ) );
      RealMatrix predictor_options_subset( Teuchos::View, 
					   predictor_options_list_in,
					   len_options, master_load, 
					   0, 
					   load*(num_processors()-1));
      predictorOptionsList_ = predictor_options_subset;
      if ( verbosity_ > 0 )
	std::cout << "Master has set predictor options\n";
    }
  else
    {
      RealMatrix predictor_options_subset;
      receive( predictor_options_subset, master_processor_id(), 
	       MPICommunicator_ );
      predictorOptionsList_ = predictor_options_subset;
      if ( verbosity_ > 0 )
	std::cout << "Slave " << processorId_ << " has set predictor options\n";
    }

#else

  predictorOptionsList_ = predictor_options_list_in;

#endif
};

void CrossValidationIterator::seed( int seed_in )
{
  seed_ = seed_in;
};

void CrossValidationIterator::verbosity( int verbosity_in )
{
  verbosity_ = verbosity_in;
};

void CrossValidationIterator::k_folds_check_inputs()
{
  if ( numPartitions_ <= 0 )
    {
      std::stringstream msg;
      msg << "CrossValidationIterator:setup_k_folds() The number of folds must";
      msg << "be greater than 0";
      throw( std::runtime_error( msg.str() ) );
    }
  
  if ( numPartitions_ > numFunctionSamples_ )
    {
      std::stringstream msg;
      msg << "CrossValidationIterator:setup_k_folds() The number of folds: ";
      msg << numPartitions_ << " must be less than or equal the number of ";
      msg << "samples: " << numFunctionSamples_;
      throw( std::runtime_error( msg.str() ) );
    }
};

int CrossValidationIterator::num_partitions()
{
  return numPartitions_;
};

void CrossValidationIterator::get_partition_indices( RealMatrix &training_indices, RealMatrix &validation_indices )
{
   training_indices = trainingIndices_;
   validation_indices = validationIndices_;
};

} //namespace Pecos

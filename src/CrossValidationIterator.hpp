#ifndef CROSS_VALIDATION_ITERATOR_HPP
#define CROSS_VALIDATION_ITERATOR_HPP

#include "LinearAlgebra.hpp"
#include "ParallelObject.hpp"
#include "FaultTolerance.hpp"

namespace Pecos {

/** \brief Compute the cross validation error indicator.
 *
 * \param validation_values the 'truth' values.
 *
 * \param prediction_values the values predicted by the predictor
 *
 * \param indicators. The error indicator for each quantity of interest,
 * i.e. for each column in prediction_values
 */
typedef void ( IndicatorFunction )( RealMatrix &validation_values,
				    RealMatrix &prediction_values,
				    RealMatrix &indicators );

/** \brief Given a realization of the predictor options, build
 * a set of predictors from training data and return
 * the cross validation error indicators of each predictor built.
 *
 * Although only one realization of the predictor options is given the 
 * analyser may produce a set of predictors. For example, if the analyser
 * builds a set of linear models using LARS. LARS will return a predictor
 * corresponding to each step of the algorithm.
 *
 * \param indicator_function specifies how to calculate the error indicators
 *
 * \param indicators (output) the error indicators for each predictor built.
 *
 * \param predictor_options_list (output) A list of matricies that contatin 
 * the matrix of the predictor options associated with each predictor built.
 * Each list entry corresponds to a quantity of interest, i.e. a column
 * of the values matrix. Each realisation of the predictor options is 
 * stored in a column of the  matrix. Example: if the analyser
 * builds a set of linear models using LARS then the predictor_options will
 * contain the predictor_opts but will also set the num_non_zeros.
 */
typedef void ( Analyser )( RealMatrix &training_samples, 
			   RealMatrix &training_values, 
			   RealMatrix &validation_samples,
			   RealMatrix &validation_values,
			   RealVector &predictor_opts,
			   IndicatorFunction *indicator_function,
			   RealMatrixList &indicators_list,
			   RealMatrixList &predictor_options_list,
			   FaultInfo &fault_info,
			   const SizetShortMap& failed_resp_data,
			   IntVector &training_indices );

/** \brief Select the 'best' predictor from a set of indicators.
 *
 * \param partition_indicators  
 * ( ( ( num_predictors x 1 ) x num_qoi ) x num_partitions )
 * 4D object containing the cross validation error indicators for a set of 
 * predictors built
 * on each partitioned data set created during cross validation.
 *
 * \param best_indices (ouput) ( num_qoi x 1 ) vector containing
 * the index of the 'best' predictor for each quantity of interest. 
 *
 * \param best_predictor_indicators (ouput) ( num_qoi x 1 ) vector containing
 * the error indicator of the 'best' predictor for 
 * each quantity of interest. It is a function of the partition indicators,
 * e.g. mean
 *
 * \param best_predictor_partition_indicators (output) 
 * ( num_partitions x num_qoi ) matrix containing the indicators of the 'best' 
 * predictor for each quantity of interest and partition.
 */
typedef void ( Selector )( std::vector<RealMatrixList> &partition_indicators, 
			   IntVector &best_predictor_indices,
			   RealVector &best_predictor_indicators,
			   RealMatrix &best_predictor_partition_indicators );

/**
 * \brief Extract the options that will create the best predictors on the
 * full data set. 
 *
 * Simply using the options used to build the best predictor on the 
 * cross validation trainging data will not suffice. For example, cross 
 * validation for compressed sensing attempts to find the max_num_iterations
 * and/or the residual ( epsilon ) that produces the best predictor.  
 * If LARS is used and max_num_iterations = the size of the training set
 * this would imply that all variables are important, but setting 
 * max_num_iterations = the size of the training set for the full
 * data set would likely stop the iteration to early
 *
 * \param partition_options  
 * ( ( ( len opts x num_predictors  ) x num_qoi ) x num_partitions )
 * 4D object containing the cross validation error indicators for a set of 
 * predictors built
 * on each partitioned data set created during cross validation.
 *
 * \param best_predictor_indices. The indices that allow us to find the best
 * indicators in partition_options
 *
 * \param num_training_samples the size of the cross validation training sample.
 * For k-folds cross validation not all training samples will be the same size
 * so use the size of the training set of the last partition.
 *
 * \param num_samples the size of the entire data.
 *
 * \param best_predictor_options (output) ( len_opts x num_qoi ) matrix 
 * containg the options that will be used to build the best predictors on
 * the full data set.
 */
typedef void ( BestOptionsExtractor )( std::vector<RealMatrixList> &partition_options, IntVector &best_predictor_indices, int num_training_samples, int num_samples, RealMatrix &best_predictor_options );

/**
 * \class CrossValidationIterator
 * \brief Implements the generic components of cross validation
 */
class CrossValidationIterator : public ParallelObject
{
protected:
  
  /// The number of data partitions used in the cross validation procedure
  int numPartitions_;

  /// The seed used to control the random number generator
  int seed_;

  /// Samples used to build the predictors
  RealMatrix samples_;

  /// Values used to build the predictors
  RealMatrix values_;

  /// Indices partitioning the data into k-1 folds used to build
  /// the predictors
  IntMatrix trainingIndices_;

  /// Indices partitioning the data into 1 folds for computing cross validation
  /// score of the predictors
  IntMatrix validationIndices_;

  /// The size of the training set for each partition
  IntVector trainingIndicesSizes_;

    /// The size of the validation set for each partition
  IntVector validationIndicesSizes_;

  /// A set of options that will be used to create instances of the predictor
  RealMatrix predictorOptionsList_;

  /// The error indicator of the best predictor for each QOI
  RealVector bestPredictorIndicators_;

  /// The predictor options the produce the best predictor for each QOI
  RealMatrix bestPredictorOptions_;

  /// The options of the best predictors of the predictors produced by each item 
  /// in predictorOptionsList_.
  RealMatrixList predictorOptionsHistory_;

  /// The best predictors of the predictors produced by each item 
  /// in predictorOptionsList_.
  RealMatrixList predictorIndicatorsHistory_;

  /// The indicators of each partition for the best predictors of the 
  /// predictors produced by each item in predictorOptionsList_.
  RealMatrixList predictorPartitionIndicatorsHistory_;
  
  /// Controls the amount of information written to standard I/O
  int verbosity_;

  /// Specify whether values_ contains gradient information
  bool useGradients_;
  
  /// Specify the number of function values in samples_
  int numFunctionSamples_;


  /// Set the indices that partition the data into training and validation sets
  void set_partition_indices( IntMatrix &training_indices, 
			      IntMatrix &validation_indices,
			      IntVector &trainining_indices_sizes,
			      IntVector &validation_indices_sizes );


  /** \brief Update the best predictor history and the best predictor 
   *  information, that will allows us to build the best predictor.
   */
  void update_best_predictor_info( RealVector &best_predictor_indicators, 
				   RealMatrix &best_predictor_options, 
			       RealMatrix &best_predictor_partition_indicators );

  /** \brief Check that the inputs to k_folds cross validation have been set
   *
   * Specifically check that the number of samples is consistent with the number
   * of folds;
   */
  void k_folds_check_inputs();

  /** \brief reshape history data to a list of size num_rhs with each list item
      containing data realated to each predictorOptionsList_;
   */
  void reshape_history_data( RealMatrixList &predictor_options_history,
			     RealMatrixList &predictor_indicators_history, 
			RealMatrixList &predictor_partition_indicators_history );

public:

  /// Default constructor
  CrossValidationIterator()
    : numPartitions_( 10 ), seed_( 0 ), verbosity_( 0 ), useGradients_( false ),
      numFunctionSamples_( 0 )
  {};
  
  /// Deconstructor
  ~CrossValidationIterator(){};

  /** \brief Return the indices that partition the data into 
   * training and validation sets using K-fold cross validation splitting.
   *
   * Each partition is split into a training and validation set. The training
   * set consists of k-1 folds and the validation set contains the remaing fold.
   * The k folds will likely be of different size. Here we ensure the size of 
   * each fold differs by no more than one. We name these folds small and large
   * folds. The later containg one more point than the former. Define 
   * num_samples_per_fold = floor( num_samples / num_partitions ) then of
   * the k partitions (num_samples - num_partitions * num_samples_per_fold)
   * will contain (num_samples_per_fold + 1) samples. These folds are refered 
   * to as large folds. The remaining
   * (num_partitions * ( num_samples_per_fold + 1 ) - num_samples)
   * folds will contain num_samples_per_fold. These folds are refered 
   * to as small folds.
   * Due to the varying size of the folds, the size of the validation and 
   * training data sets can vary between partitions. When a validation or 
   * training set of a partition is less than the maximum possible size. 
   * The remaining entries of the index vector are set to -1. 
   * For example if num_samples_per_fold = 2
   * a validation set can consist of 2 or 3 entries. When the jth validation fold
   * contains only two entries the third entry of validation_indices[j], i.e.
   * validation_indices(2,j) = -1
   */
  void setup_k_folds( int num_partitions = 10 );


  /** \brief Create the training and validation data for a partition
   *
   * If gradient information is used, this function
   * assumes that the function data is stored in the first (0-num_pts) 
   * rows of values_. Then the partial derivative with respect to the first
   * dimension is stored in (num_pts-2*num_pts) rows of the matrix and so on.
   */
  void partition_data( RealMatrix &training_samples, 
		       RealMatrix &training_values, 
		       RealMatrix &validation_samples,
		       RealMatrix &validation_values, 
		       int fold );


  /** \brief Run cross validation
   *
   * \param indicator_function The cross validation error indicator
   *
   * \param analyser Function that builds a predictor with one set of predictor
   * options and uses indicator_function to determine the quantify the 
   * cross validation accuracy of the predictor that was built
   *
   * \param selector Function that selects the 'best' predictor from 
   * a set of indicators
   *
   * \param best_options_extractor Function that extracts the options 
   * that will create the best predictors on the full data set. 
   */
  void run( IndicatorFunction *indicator_function,
	    Analyser* analyser, Selector *selector,
	    BestOptionsExtractor *best_options_extractor,
	    FaultInfo &fault_info,
	    const SizetShortMap& failed_resp_data );

  /** \brief Return the options that create the best predictor for each QOI.
   * Also return the error indicators of each best predictor.
   * 
   * \param best_predictor_options The predictor options the produce the best 
   * predictor for each QOI
   * 
   * \param best_predictor_indicators The error indicator of the best predictor 
   * for each QOI
   */
  void get_best_predictor_info( RealMatrix &best_predictor_options,
				RealVector &best_predictor_indicators );

  /** \brief Return the predictor error indicators, the predictor partition
   * error indicators and the predictor options of all the predictors
   * created duing the cross validation. 
   *
   * Note predictor_options_history
   * will unlikely be the same as predictorOptionsList_ as the entries in
   * the former represent the best predictor for each of the entries in the 
   * later. 
   *
   * \param predictor_options_history The options of the best predictors of
   * the predictors produced by each item in predictorOptionsList_.
   *
   * \param predictor_indicators_history The best predictors of the predictors 
   * produced by each item in predictorOptionsList_.
   *
   * \param predictor_partition_indicators_history The indicators of each 
   * partition for the best predictors of the predictors produced by each item 
   * in predictorOptionsList_
   */ 
  void get_history_data( RealMatrixList &predictor_options_history,
			 RealMatrixList &predictor_indicators_history,
			 RealMatrixList &predictor_partition_indicators_history);
  
  
  /** \brief Set the data to be used in the cross validation study
   * \param samples the sample data
   * \param values the value data
   * \param num_function_samples the number of function samples in samples.
   * If gradients is used then set
   * num_function_samples = num_dims * num_total_samples / ( num_dims + 1 )
   * if gradients are not used then set
   * num_function_samples = num_total_samples or 0. 
   * num_function_samples = 0 will cause 
   * num_function_samples = num_total_samples to be set internally
   */
  void set_data( RealMatrix &samples, RealMatrix &values,
		 int num_function_samples = 0 );


  /// Set the predictor options that will be tested during cross validation
  void predictor_options_list( RealMatrix &predictor_options_list_in );

  /// Set the seed used to control the random number generator
  void seed( int seed_in );

  /// Set the verbosity level
  void verbosity( int verbosity_in );

  /// Return the number of partitions used in the cross validation
  int num_partitions();

  /// Get the indices used to partition the data into training and
  /// validation sets
  void get_partition_indices( IntMatrix &training_indices, 
			      IntMatrix &validation_indices );
  
};

} // namespace Pecos
#endif

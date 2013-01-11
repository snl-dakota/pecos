#ifndef LINEAR_MODEL_MODULES_HPP
#define LINEAR_MODEL_MODULES_HPP

#include "CrossValidationIterator.hpp"
#include "CompressedSensing.hpp"

namespace Pecos {

void set_linear_predictor_options( RealVector &predictor_opts,  
				   CompressedSensingOptions &cs_opts );

void extract_linear_predictor_options(CompressedSensingOptionsList &cs_opts_list,
				      RealMatrixList &predictor_options_list );

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
				IntVector &validation_indices );

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
void linear_predictor_best_options_extractor( std::vector<RealMatrixList> &partition_options, IntVector &best_predictor_indices, int num_training_samples, int num_samples, RealMatrix &best_predictor_options );

} // namespace Pecos

#endif

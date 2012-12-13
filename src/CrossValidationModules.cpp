#include "CrossValidationModules.hpp"
#include "MathTools.hpp"

namespace Pecos {

void rmse_indicator( RealMatrix &validation_values,
		     RealMatrix &prediction_values,
		     RealMatrix &indicators )
{
  Teuchos::BLAS<int, Real> blas;
  RealMatrix diff( validation_values );
  diff -= prediction_values;

  int M( diff.numRows() ), num_qoi( diff.numCols() );
  reshape( indicators, num_qoi, 1 );
  for ( int k = 0; k < num_qoi; k++ )
    {
      Real *diffCol = diff[k];
      //Sum the absolute values of the entries of col. 
      indicators(k,0) = blas.NRM2( M, diffCol, 1 ) / sqrt( (Real)M );
    }
};

void normalised_mean_selector( std::vector<RealMatrixList> &partition_indicators,
			       IntVector &best_predictor_indices,
			       RealVector &best_predictor_indicators,
			       RealMatrix &best_predictor_partition_indicators )
{
  int num_partitions( partition_indicators.size() ), 
    num_qoi( partition_indicators[0].size() );

  best_predictor_indices.resize( num_qoi );
  best_predictor_indicators.resize( num_qoi );
  best_predictor_partition_indicators.shapeUninitialized( num_partitions, 
							  num_qoi );

  for ( int k = 0; k < num_qoi; k++ )
    {

      // Some predictor types will produce different number of predictors
      // for the same predictor option but on different partitions. We need
      // to allocate memory to allow for this situation
      int max_num_predictors( 0 );
      for ( int i = 0; i < num_partitions; i++ )
	{
	  int num_predictors( partition_indicators[i][k].numRows() );
	  max_num_predictors = std::max( max_num_predictors, num_predictors );
	}
      
      RealMatrix values( num_partitions, max_num_predictors, false );      
      for ( int i = 0; i < num_partitions; i++ )
	{
	  int num_predictors( partition_indicators[i][k].numRows() );
	  for ( int j = 0; j < num_predictors; j++ )
	    { 
	      values(i,j) = partition_indicators[i][k](j,0);
	    }
	  for ( int j = num_predictors; j < max_num_predictors; j++ )
	    {
	      values(i,j) = std::numeric_limits<Real>::max();
	    }
	}

      int argmin( 0 );
      Real min_objective( std::numeric_limits<Real>::max() );
      for ( int j = 0; j < values.numCols(); j++ )
	{
	  // Compute the mean and variance of the error 
	  // of each cross validation partition;
	  Real cv_mean = mean(  num_partitions, values[j] );
	  Real cv_var = variance( num_partitions, values[j] );
	  Real objective = std::abs( cv_mean ) * cv_var;

	  if (  objective < min_objective )
	    {
	      argmin = j;
	      min_objective = objective;
	    }
	}
      best_predictor_indices[k] = argmin;
      best_predictor_indicators[k] = min_objective;
      for ( int i = 0; i < num_partitions; i++ )
	{
	  best_predictor_partition_indicators(i,k) = 
	    partition_indicators[i][k](argmin,0);
	}
    }
};
  
} // namespace Pecos

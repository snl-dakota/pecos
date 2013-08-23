#ifndef CROSS_VALIDATION_MODULES_HPP
#define CROSS_VALIDATION_MODULES_HPP

#include "LinearAlgebra.hpp"

namespace Pecos {

void rmse_indicator( RealMatrix &validation_values,
		     RealMatrix &prediction_values,
		     RealMatrix &indicators );

void normalised_mean_selector( RealMatrix2DArray &partition_indicators,
			       IntVector &best_predictor_indices,
			       RealVector &best_predictor_indicators,
			       RealMatrix &best_predictor_partition_indicators );
} // namespace Pecos
#endif


#include "RandomVariables.hpp"

using namespace Surrogates;

void
RandomVariables::generate_samples(int num_vars, 
                                  int num_samples,
      		                  const Teuchos::ParameterList & params, 
                                  IntMatrix & samples)
{ 
  samples.shapeUninitialized(num_vars, num_samples);
  return RandomNumberGenerator::generate(samples, num_vars, num_samples, params);
}


void
RandomVariables::generate_samples(int num_vars, 
                                  int num_samples,
      		                  const Teuchos::ParameterList & params, 
                                  RealMatrix & samples)
{ 
  samples.shapeUninitialized(num_vars, num_samples);
  return RandomNumberGenerator::generate(samples, num_vars, num_samples, params);
}

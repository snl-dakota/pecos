/**
 * \file MonteCarlo.hpp
 * \brief A Monte Carlo sampler that includes methods to calculate 
 * basic statistics
 * \author John D. Jakeman 
 */

// http://www.stat.cmu.edu/~abrock/beowulf/mpiandc.hpptml

#ifndef MONTE_CARLO_HPP
#define MONTE_CARLO_HPP

#include "LinearAlgebra.hpp"
#include "Model.hpp"
#include "RandomNumberGenerator.hpp"

enum distribution
{
  uniform,
  beta,
  gaussian
};

/**
 * \class MonteCarlo
 * \brief A Monte Carlo sampler that includes methods to calculate 
 * basic statistics 
 */
class MonteCarlo
{
protected:
  RandomNumberGenerator RNG_;
 
 public:

  MonteCarlo();

  ~MonteCarlo();

  void generateSamples( RealMatrix &samples, int nSamples, distribution dist,
			RealVector &distParameters, RealVector &domain,
			unsigned int seed = 0 );

  void latin_hypercube_design( int num_dims, int nSamples, RealMatrix &samples, 
			     distribution dist, RealVector &distParameters, 
			     RealVector &domain, unsigned int seed );

  /**
   * \brief compute the sample mean of each row of a matrix
   */
  void mean( RealMatrix &values, RealVector& mean );

  /**
   * \brief compute the sample variance of each row of a matrix
   */
  void variance( RealMatrix &values, RealVector& variance );
};
#endif

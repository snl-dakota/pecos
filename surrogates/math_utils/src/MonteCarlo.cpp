#include "MonteCarlo.hpp"
#include "Printing.hpp"
#include "MathTools.hpp"
#include <boost/math/distributions/normal.hpp>

MonteCarlo::MonteCarlo(){};

MonteCarlo::~MonteCarlo(){};

void MonteCarlo::generateSamples( RealMatrix& samples, int num_samples,
				  distribution dist, RealVector &dist_parameters,
				  RealVector &domain, unsigned int seed )
{
  int num_dims = domain.length() / 2;
  switch( dist )
    {
    case uniform:
      {
	RNG_.uniform( samples, num_dims, num_samples, seed );
	rescale( samples, domain, 0 );
	break;
      }
    case beta:
      {
	if ( dist_parameters.length() != 2 )
	  {
	    // use uniform distribution as default
	    dist_parameters.resize( 2 );
	    dist_parameters[0] = 1.0; 
	    dist_parameters[1] = 1.0;
	  }
	Real alpha( dist_parameters[0] ), beta( dist_parameters[1] );
	RNG_.beta( samples, num_dims, num_samples, alpha, beta, seed );
	rescale( samples, domain, 0 );
	break;
      }
    case gaussian:
      {
	if ( dist_parameters.length() != 2 )
	  {
	    // use standard normal as default
	    dist_parameters.resize( 2 );
	    dist_parameters[0] = 0.0; 
	    dist_parameters[1] = 1.0;
	  }
	Real mean( dist_parameters[0] ), variance( dist_parameters[1] );
	RNG_.gaussian( samples, num_dims, num_samples, mean, variance, seed );
	break;
      }
    default:
      {
	std::string msg = "MonteCarlo:generateSamples() ";
	msg += "Incorrect distribution specified";
	throw( std::runtime_error( msg ) );
      }
    }
};

void MonteCarlo::mean( RealMatrix& values, RealVector& mu )
{
  int num_samples( values.numRows() ), num_qoi( values.numCols() );
  mu.resize( num_qoi );
  for ( int i = 0; i < num_qoi; i++ )
    {
      mu[i] = 0.0;
      for ( int j = 0; j < num_samples; j++ )
	{
	  mu[i] += values(j,i);
	}
      mu[i] /= (Real)num_samples;
    } 
};

void MonteCarlo::variance( RealMatrix& values, RealVector& var )
{ 
  int num_samples( values.numRows() ), num_qoi( values.numCols() );
  RealVector mean( num_qoi );
  var.resize( num_qoi );
  for ( int i = 0; i < num_qoi; i++ )
    {
      var[i] = 0.0;
      mean[i] = 0.0;
      for ( int j = 0; j < num_samples; j++ )
	{
	  mean[i] += values(j,i);
	  var[i] += values(j,i)*values(j,i);
	}
      mean[i] /= (Real)num_samples;
      var[i] = ( var[i]  - (Real)num_samples * mean[i] * mean[i] ) / 
	( (Real)num_samples  - 1.0 );
    } 
};

void MonteCarlo::latin_hypercube_design( int num_dims, int num_samples, 
				       RealMatrix &samples, distribution dist,
				       RealVector &dist_parameters,
				       RealVector &domain, unsigned int seed )
{
  RealMatrix lhs;
  ::latin_hypercube_design( num_samples, num_dims, lhs, seed );
  switch( dist )
    {
    case uniform:
      {
	// return base hypercube
	samples = lhs;
	break;
      }
    case beta:
      {
	if ( dist_parameters.length() != 2 )
	  {
	    // use uniform distribution as default
	    dist_parameters.resize( 2 );
	    dist_parameters[0] = 1.0; 
	    dist_parameters[1] = 1.0;
	  }
	Real alpha( dist_parameters[0] ), beta( dist_parameters[1] );
	
	boost::math::beta_distribution<Real> beta_distribution( alpha, beta );
	
	samples.reshape( num_dims, num_samples );
	
	for ( int i = 0; i < num_dims; i++)
	  {
	    for ( int j = 0; j < num_samples; j++)
	      {
		// Evaluate Beta CDF at u;
		samples(i,j) = boost::math::quantile( beta_distribution, 
						      lhs(i,j) );
	      }
	  }

	rescale( samples, domain, 0 );
	break;
      }
    case gaussian:
      {
	if ( dist_parameters.length() != 2 )
	  {
	    // use standard normal as default
	    dist_parameters.resize( 2 );
	    dist_parameters[0] = 0.0; 
	    dist_parameters[1] = 1.0;
	  }
	Real mu( dist_parameters[0] ), sigma2( dist_parameters[1] );

	boost::math::normal normal_distribution( mu , sigma2 ); 
       	
	samples.reshape( num_dims, num_samples );
	
	for ( int i = 0; i < num_dims; i++)
	  {
	    for ( int j = 0; j < num_samples; j++)
	      {
		// Evaluate Beta CDF at u;
		samples(i,j) = boost::math::quantile( normal_distribution, 
						      lhs(i,j) );
	      }
	  }
      }
    default:
      {
	std::string msg = "MonteCarlo:latin_hypercube_design() ";
	msg += "Incorrect distribution specified";
	throw( std::runtime_error( msg ) );
      }
    };
};

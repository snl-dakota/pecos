
#include <string>

#include <Teuchos_UnitTestHarness.hpp> 
#include <Teuchos_ParameterList.hpp> 

#include <RandomVariables.hpp>

namespace {

  template<typename T>
  void compute_stats(const T& matrix_values, RealMatrix& stat_vals)
  {
    stat_vals.shapeUninitialized(matrix_values.numRows(),3);
    Real avg = 0.0, min = 1.e15, max = -1.e15;
    for( int i=0; i<matrix_values.numRows(); ++i ) {
      avg = 0.0;
      min = 1.e15;
      max = -1.e15;
      for( int j=0; j<matrix_values.numCols(); ++j ) {
        avg += (Real)matrix_values(i,j);
        if( min > matrix_values(i,j) )
          min = matrix_values(i,j);
        if( max < matrix_values(i,j) )
          max = matrix_values(i,j);
      }
      stat_vals(i,0) = min;
      stat_vals(i,1) = max;
      stat_vals(i,2) = avg/(Real)matrix_values.numCols();
    }
  }

}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(math_tools, random_vars_syntax)
{
  int num_vars = 1, num_samples = 1;
  RealMatrix samples;

  Teuchos::ParameterList params;
  // Test for required behavior with incorrect syntax, eg missing required "distribution" setting
  TEST_THROW( Surrogates::RandomVariables::generate_samples(num_vars, num_samples, params, samples), std::exception );

  params.set("distribution", "uniform");
  // Test for required behavior with incorrect syntax, eg missing required "seed" setting
  TEST_THROW( Surrogates::RandomVariables::generate_samples(num_vars, num_samples, params, samples), std::exception );

  params.set("seed", 1337);
  // Should work now
  TEST_NOTHROW( Surrogates::RandomVariables::generate_samples(num_vars, num_samples, params, samples) );

  // Test for required behavior with incorrect syntax, eg incorrect "distribution" setting
  params.set("distribution", "not_there");
  TEST_THROW( Surrogates::RandomVariables::generate_samples(num_vars, num_samples, params, samples), std::runtime_error );

  params.set("distribution", "normal");
  // Test for required behavior with incorrect syntax, eg missing required "mean" and "variance" settings
  TEST_THROW( Surrogates::RandomVariables::generate_samples(num_vars, num_samples, params, samples), std::exception );
  params.set("mean", 0.0);
  TEST_THROW( Surrogates::RandomVariables::generate_samples(num_vars, num_samples, params, samples), std::exception );
  params.set("variance", 1.0);
  TEST_NOTHROW( Surrogates::RandomVariables::generate_samples(num_vars, num_samples, params, samples) );

  // Now test beta distribution params
  Teuchos::ParameterList params2;
  params2.set("distribution", "beta");
  params2.set("seed", 1337);
  // Test for required behavior with incorrect syntax, eg missing required "alpha" and "beta" settings
  TEST_THROW( Surrogates::RandomVariables::generate_samples(num_vars, num_samples, params2, samples), std::exception );
  params2.set("alpha", 1.0);
  TEST_THROW( Surrogates::RandomVariables::generate_samples(num_vars, num_samples, params2, samples), std::exception );
  params2.set("beta", 1.0);
  TEST_NOTHROW( Surrogates::RandomVariables::generate_samples(num_vars, num_samples, params2, samples) );
}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(math_tools, random_vars)
{
  Teuchos::ParameterList params;
  params.set("seed", 1337);
  int num_vars = 3, num_samples = 200;

  RealMatrix samples;
  // What we are testing
  // ----------------------------
  params.set("distribution", "uniform");
  Surrogates::RandomVariables::generate_samples(num_vars, num_samples, params, samples);
  // ----------------------------

  // Test structure of returned samples
  TEST_EQUALITY( num_vars, samples.numRows() );
  TEST_EQUALITY( num_samples, samples.numCols() );

  // Test statistics of samples - what is reasonable here? - RWH
  RealMatrix stats;
  compute_stats(samples, stats);
  Real diff = 0.0;
  //stats.print(std::cout);
  for( int i=0; i<stats.numRows(); ++i ) {
    TEST_ASSERT( (stats(i,0) > 0.0) && (stats(i,0) < 0.01) ); // test min value within 1%
    diff = std::abs(stats(i,1) - 1.0);
    TEST_ASSERT( (stats(i,1) <= 1.0) && (diff < 0.02) ); // test max value within 1%
    diff = std::abs(stats(i,2) - 0.5);
    TEST_ASSERT( diff < 0.03 ); // test average value within 3%
  }


  // Now test int version of uniform distribution
  IntMatrix int_samples;
  // What we are testing
  // ----------------------------
  params.set("distribution", "uniform");
  Surrogates::RandomVariables::generate_samples(num_vars, num_samples, params, int_samples);
  // ----------------------------

  // Test structure of returned samples
  TEST_EQUALITY( num_vars, int_samples.numRows() );
  TEST_EQUALITY( num_samples, int_samples.numCols() );

  // Test statistics of samples - what is reasonable here? - RWH
  compute_stats(int_samples, stats);
  diff = 0.0;
  //stats.print(std::cout);
  for( int i=0; i<stats.numRows(); ++i ) {
    TEST_EQUALITY( stats(i,0), 0 ); // test min value
    TEST_EQUALITY( stats(i,1), 1 ); // test max value
    diff = std::abs(stats(i,2) - 0.5);
    TEST_ASSERT( diff < 0.06 ); // test average value within 6%
  }


  // Now test normal distribution
  samples.shape(1,1);;
  params.set("seed", 11234);
  num_vars = 2, num_samples = 500;
  // What we are testing
  // ----------------------------
  params.set("distribution", "normal");
  params.set("mean", 2.0);
  params.set("variance", 1.5);
  Surrogates::RandomVariables::generate_samples(num_vars, num_samples, params, samples);
  // ----------------------------

  // Test structure of returned samples
  TEST_EQUALITY( num_vars, samples.numRows() );
  TEST_EQUALITY( num_samples, samples.numCols() );

  // Test statistics of samples - what is reasonable here? - RWH
  compute_stats(samples, stats);
  diff = 0.0;
  //stats.print(std::cout);
  for( int i=0; i<stats.numRows(); ++i ) {
    diff = std::abs(stats(i,2) - 2.0);
    TEST_ASSERT( diff < 0.09 ); // test average value within 9%
  }


  // Now test beta distribution
  samples.shape(1,1);;
  Teuchos::ParameterList params2;
  params2.set("seed", 1337);
  num_vars = 3, num_samples = 300;
  // What we are testing
  // ----------------------------
  params.set("distribution", "beta");
  Real alpha = 2.0, beta = 5.0;
  params.set("alpha", alpha);
  params.set("beta", beta);
  Surrogates::RandomVariables::generate_samples(num_vars, num_samples, params, samples);
  // ----------------------------

  // Test structure of returned samples
  TEST_EQUALITY( num_vars, samples.numRows() );
  TEST_EQUALITY( num_samples, samples.numCols() );

  // Test statistics of samples - what is reasonable here? - RWH
  compute_stats(samples, stats);
  diff = 0.0;
  //stats.print(std::cout);
  for( int i=0; i<stats.numRows(); ++i ) {
    diff = std::abs(stats(i,2) - alpha/(alpha+beta));
    TEST_ASSERT( diff < 0.09 ); // test average value within 9%
  }
}

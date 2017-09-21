#include <ctype.h>
#include <string>
#include <cmath>
#include <iostream>

#include <Teuchos_UnitTestHarness.hpp>

#include <teuchos_data_types.hpp>
#include <CppFunction.hpp>
#include <Monomial.hpp>
#include <BoundedVariables.hpp>
#include <AffineVariableTransformation.hpp>
#include <RegressionBuilder.hpp>

using namespace Surrogates;

namespace {

//--------------------------------------
// Helper functions


/**
\brief Evaluate a quadratic polynomial with no interaction terms, i.e.
*  \[q_i = (i+1) + \sum_j (i+1)*x_j^2+(i+2)x_j\]
*/
void additive_quadratic_function(const Real* sample, Real* func_vals, const OptionsList &opts){
   const int num_vars = opts.get<int>("num_vars");
   const int num_qoi =  opts.get<int>("num_qoi");

   for (int i=0; i<num_qoi; i++){
     func_vals[i] = (Real)(i+1);
     for (int j=0; j<num_vars; j++){
       func_vals[i] += (Real)(i+1)*sample[j]*sample[j]+(Real)(i+2)*sample[j];
     }
   }
}

bool test_polynomial_regression_builder(RegressionType regression_type, Real atol){
  int num_vars = 2, num_qoi = 2, num_samples = 20;
  int degree = 3;

  // Generate realizations of function variables

  // Define the function variables
  RealVector ranges;
  define_homogeneous_ranges(num_vars, 0., 1., ranges);
  boost::shared_ptr<Variables> variables(new BoundedVariables);
  boost::dynamic_pointer_cast<BoundedVariables>(variables)->set_ranges(ranges);

  // Define the variable transformation. Need to decide if a seperate
  // transformation should exist for approximation and for mapping samples
  // generated. For example approximation accepts samples in user x-space
  // but operates in a standardized u-space. But sample generator produces
  // samples in (possibly another standardized space) and must map these to
  // x-space, or even to approximation u-space.
  boost::shared_ptr<VariableTransformation>
    var_transform(new AffineVariableTransformation);
  var_transform->set_variables(variables);

  // Define the function to approximate
  RealMatrix func_vals;
  CppFunction model;
  OptionsList model_opts;
  model_opts.set("num_vars",num_vars);//total number of variables
  model_opts.set("num_qoi",num_qoi);//total number of quantities of interest
  model.set_options(model_opts);
  model.set_function(&additive_quadratic_function);

  // Initialize the approximation
  OptionsList monomial_opts;
  monomial_opts.set("max_total_degree",degree);//max degree of total-degree polynomial space
  monomial_opts.set("num_qoi", num_qoi);//number of quantites of interest
  Monomial monomial;
  monomial.set_options(monomial_opts);
  monomial.set_variable_transformation(var_transform);
  IntMatrix basis_indices;
  compute_hyperbolic_indices(num_vars, degree, 1., basis_indices);
  monomial.set_basis_indices(basis_indices);

  // Generate the approximation coefficients using a regression based method
  //     generate at set of samples to build approximation
  OptionsList regression_opts;
  regression_opts.set("regression_type",SVD_LEAST_SQ_REGRESSION);
  regression_opts.set("num_samples",num_samples);
  regression_opts.set("sample_type","probabilistic_MC");
  RegressionBuilder builder;
  builder.set_target_function(model);
  builder.build(regression_opts, monomial);
  

  // Check that approximation is an interpolant
  // TODO samples and function values (in a data object) should be return
  // in status opts so that they can be accessed outide driver like here.
  // requires wrapping or serial dense matrix for parameter list
  RealMatrix approx_vals;
  //monomial.value(samples,approx_vals);
  if (!allclose(approx_vals,func_vals,atol))
    return false;

  // Check that approximation is exact everywhere
  RealMatrix test_samples;
  int num_test_samples = 100;
  int seed = 1337;
  generate_uniform_samples(num_vars, num_test_samples, seed, *var_transform,
		      test_samples);

  model.value(test_samples, func_vals);
  monomial.value(test_samples,approx_vals);
  if (!allclose(approx_vals,func_vals,atol))
    return false;

  // Check the monomial coefficients
  RealMatrix coeffs;
  monomial.get_coefficients(coeffs);
  Real true_coeffs_array[] = {1, 2, 2, 1, 1, 0, 0, 0, 0, 0,
			      2, 3, 3, 2, 2, 0, 0, 0, 0, 0, };
  RealMatrix true_coeffs(Teuchos::Copy,true_coeffs_array,10,10,2);
  //coeffs.print(std::cout);
  //true_coeffs.print(std::cout);
  return allclose(coeffs,true_coeffs,atol);

}

//----------------------------------------------------------------
// Unit tests

/**
\brief Generate a monomial interpolant using qr factorization to
solve a determined linear system Ac=f, where A is the vandermonde-type
interpolation matrix and f(z) is the function values at a set of samples z,
and c is the vector of monomial coefficients
*/
TEUCHOS_UNIT_TEST(polynomial_approximation, monomial)
{
  TEST_ASSERT(test_polynomial_regression_builder(ORTHOG_MATCH_PURSUIT, 1e-8));

  TEST_ASSERT(test_polynomial_regression_builder(LEAST_ANGLE_REGRESSION, 1e-8));

  TEST_ASSERT(test_polynomial_regression_builder(LASSO_REGRESSION, 1e-8));

  TEST_ASSERT(test_polynomial_regression_builder(SVD_LEAST_SQ_REGRESSION, 1e-8));
}

}
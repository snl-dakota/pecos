#ifndef POLYNOMIAL_APPROXIMATION_DRIVERS_HPP
#define POLYNOMIAL_APPROXIMATION_DRIVERS_HPP

#include <PolynomialApproximation.hpp>
#include <math_tools.hpp>
// TODO: to be removed
#include <boost/random.hpp>
#include <boost/random/uniform_real.hpp>

#ifdef NEVER

namespace Surrogates {

  enum SolverTypes{QR_LSQ,SVD_LSQ};
  

void regression_solve(const RealMatrix &matrix, const RealMatrix &rhs,
		      const OptionsList &opts,
		      RealMatrix &coeffs, OptionsList &status){
  switch (opts.get<int>("solver_type")){
    case QR_LSQ:
      //resize_if_needed(coeffs, matrix.numRows(), rhs.numRows());
      qr_solve(matrix,rhs,coeffs,Teuchos::NO_TRANS);
      break;
    default:
      error("Incorrect regression solver specified.");
      break;
    }
}


  
  void regression_driver(Function &target_function, PolynomialApproximation &approx, OptionsList &opts ){

    // Generate samples to build approximation
    int num_samples = opts.get<int>("num_samples");
    std::string sample_type = opts.get<std::string>("sample_type");

    // Create mc sampler and pass in sample type. For now hack support for
    // uniform mc sampler
    RealMatrix samples;
    int seed = 1337;
    boost::shared_ptr<VariableTransformation> var_transform;
    approx.get_variable_transformation(var_transform);
    gen_uniform_samples(approx.num_vars(), num_samples, seed,
			var_transform, samples);
    samples.print(std::cout);

    // Evaluate the function at the build samples
    RealMatrix values;
    target_function.value(samples, values);
    values.print(std::cout);

    // Generate matrix of the linear sytem to be used in regression solve
    RealMatrix basis_matrix;
    approx.generate_basis_matrix(samples,basis_matrix);
    basis_matrix.print(std::cout);

    // Use regressiong to find coefficients
    OptionsList regression_status;
    RealMatrix coeffs;
    regression_solve(basis_matrix, values, opts, coeffs,
		     regression_status);
    coeffs.print(std::cout);

    // Set the approximation coefficients
    approx.set_coefficients(coeffs);

   // figure out how to return class using swig argout
   // Then consider returning a status parameterList
  }

} // namespace Surrogates
#endif

#endif //POLYNOMIAL_APPROXIMATION_DRIVERS_HPP

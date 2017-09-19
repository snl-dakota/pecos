import unittest, numpy
from PyDakota.approximation import *
from PyDakota.regression import *
from PyDakota.math_tools import *
class TestMonomialApproximation(unittest.TestCase):
    def setUp(self):
        pass

    def test_pyfunction(self):
        num_vars = 2; num_samples = 10
        additive_quadratic_function = \
          lambda x: numpy.sum(x**2)*numpy.arange(1,3) + \
          numpy.sum(x)*numpy.arange(2,4) + numpy.arange(1,3)
        function = PyFunction(additive_quadratic_function)
        samples = numpy.random.uniform(-1,1,(num_vars,num_samples))
        values = function.value(samples)

        # Check function.value produces the same result as the function
        # it is wrapping
        true_values = numpy.empty((num_samples,2),float)
        for i in xrange(num_samples):
            true_values[i,:] = additive_quadratic_function(samples[:,i])
        assert numpy.allclose(true_values, values)

        # Check that function.value can be called with a 1D array
        assert numpy.allclose(
            function.value(samples[:,0]), true_values[0,:])

        # check that function.value can be used with a function that
        # returns a scalar
        additive_quadratic_scalar_valued_function = \
          lambda x: numpy.sum(x**2) + numpy.sum(x)*2 +1.
        function = PyFunction(additive_quadratic_scalar_valued_function)
        values = function.value(samples)
        true_values = numpy.empty((num_samples,1),float)
        for i in xrange(num_samples):
            true_values[i,0] = additive_quadratic_scalar_valued_function(
                samples[:,i])
        assert numpy.allclose(true_values, values)

    def test_define_homogeneous_ranges(self):
        num_vars = 3
        ranges = define_homogeneous_ranges(num_vars, 0., 1.);
        true_ranges = numpy.ones((2*num_vars),float)
        true_ranges[::2]=0.
        assert numpy.allclose(ranges,true_ranges)

        ranges = define_homogeneous_ranges(num_vars, -3., 2.);
        true_ranges = numpy.ones((2*num_vars),float)*2.
        true_ranges[::2]=-3.
        assert numpy.allclose(ranges,true_ranges)

    def test_compute_hyperbolic_indices(self):
        num_vars = 2; degree = 3
        basis_indices = compute_hyperbolic_indices(num_vars, degree, 1.)
        num_indices = cardinality_of_total_degree_polynomial(
            num_vars,degree)
        assert basis_indices.shape==(num_vars,num_indices)
        true_indices = numpy.array(
            [[0,0],[1,0],[0,1],[2,0],[1,1],
             [0,2],[3,0],[2,1],[1,2],[3,0]]).T
        for i in xrange(num_indices):
            found = False
            for j in xrange(num_indices):
                if numpy.allclose(basis_indices[:,j],true_indices[:,j]):
                    found = True
                    break
            assert found
    
    def initialize_homogeneous_polynomial_01(self, num_vars, degree,
                                             use_tensor_product_indices,
                                             poly_type):

        # Define the function variables
        variables = BoundedVariables()
        ranges = define_homogeneous_ranges(num_vars, 0., 1.);
        variables.set_ranges(ranges)
        
        # Define the variable transformation used to covert data in
        # the user define space into the space native to the approximation
        var_transform = AffineVariableTransformation()
        var_transform.set_variables(variables)

        # TODO replace this seitch with a c++ factory
        if poly_type=='monomial':
            approx = Monomial()
        else:
            from PyDakota.univariate_polynomials import LEGENDRE_ORTHOG
            approx = PolynomialChaosExpansion()
            basis_types = [LEGENDRE_ORTHOG]*num_vars
            approx.initialize_polynomial_basis_from_basis_types(
                basis_types)

        if use_tensor_product_indices:
            degrees = numpy.array([degree]*num_vars,dtype=numpy.int32)
            basis_indices = tensor_product_indices(degrees)
        else:
            basis_indices=compute_hyperbolic_indices(num_vars, degree, 1.)
        # Check exception is thrown in basis_indices is set before
        # variable transformation
        self.assertRaises(
            RuntimeError,approx.set_basis_indices,basis_indices)
        approx.set_variable_transformation(var_transform)
        approx.set_basis_indices(basis_indices)

        return approx

    def test_generate_monomial_basis_matrix(self):
        num_vars = 2; degree = 3; use_tensor_product_indices=True
        approx = self.initialize_homogeneous_polynomial_01(
            num_vars, degree, use_tensor_product_indices, 'monomial')
        num_samples = 10

        index_sets_1d = [numpy.arange(degree+1)]*num_vars
        tensor_product_indices = cartesian_product(index_sets_1d, 1)
        assert numpy.allclose(
            tensor_product_indices,approx.get_basis_indices())

        samples = numpy.random.uniform(0,1,(num_vars,num_samples))
        basis_matrix = approx.generate_basis_matrix(samples)

        for i in range(num_samples):
            basis_vals_1d = [(2*samples[0,i]-1)**numpy.arange(degree+1)]
            for d in range(1,num_vars):
                basis_vals_1d+=[(2*samples[d,i]-1)**numpy.arange(degree+1)]
            assert numpy.allclose(
                basis_matrix[i,:],outer_product(basis_vals_1d,1))
        

    def test_build_polynomial(self):
        num_vars = 2; degree = 3; use_tensor_product_indices=False
        approx = self.initialize_homogeneous_polynomial_01(
            num_vars, degree, use_tensor_product_indices, poly_type='pce')

        # Define the function to approximate
        additive_quadratic_function = \
          lambda x: numpy.sum(x**2)*numpy.arange(1,3) + \
          numpy.sum(x)*numpy.arange(2,4) + numpy.arange(1,3)
        function = PyFunction(additive_quadratic_function)

        # Check that if value is called before coefficients is set an
        # exception is thrown
        samples = numpy.random.uniform(-1,1,(num_vars,10))
        self.assertRaises(RuntimeError,approx.value,samples)

        # Construct the approximation
        num_build_samples = 2*approx.num_terms()
        build_samples = numpy.random.uniform(0,1,(num_vars,num_build_samples))
        build_function_vals = function.value(build_samples)
        basis_matrix = approx.generate_basis_matrix(build_samples)
        #coeffs = qr_solve(basis_matrix,build_function_vals,NO_TRANS)
        opts_dict = {'regression_type':SVD_LEAST_SQ_REGRESSION}
        linsys_opts = OptionsList(opts_dict)
        solver = regression_solver_factory(linsys_opts);
        solver.solve(
            basis_matrix, build_function_vals, linsys_opts);
        coeffs = solver.get_final_solutions();
        approx.set_coefficients(coeffs)

        # Check that approximation is an interpolant
        approx_vals = approx.value(build_samples)
        assert numpy.allclose(approx_vals,build_function_vals)

        # Check that approximation is exact everywhere
        num_samples = 100
        samples = numpy.random.uniform(0,1,(num_vars,num_samples))
        approx_vals = approx.value(samples)
        function_vals = function.value(samples)
        assert numpy.allclose(approx_vals,function_vals)
        opts_dict = {'num_samples':num_build_samples,
                     'regression_type':SVD_LEAST_SQ_REGRESSION,
                     'sample_type':'probabilistic_MC'}
        regression_opts = OptionsList(opts_dict)
        builder = RegressionBuilder()
        builder.set_target_function(function)
        builder.build(regression_opts, approx)

        # Check that approximation is an interpolant
        approx_vals = approx.value(build_samples)
        assert numpy.allclose(approx_vals,build_function_vals)

        # Check that approximation is exact everywhere
        num_samples = 100
        samples = numpy.random.uniform(-1,1,(num_vars,num_samples))
        approx_vals = approx.value(samples)
        function_vals = function.value(samples)
        assert numpy.allclose(approx_vals,function_vals)

if __name__ == '__main__':
    unittest.main()

/**
 * \file RandomNumberGenerator.hpp
 * \brief Easy to use wrapper to the Boost psuedo random number generator
 * \author John D. Jakeman
 */


#ifndef RAND_NUMBER_GENERATOR_HPP
#define RAND_NUMBER_GENERATOR_HPP

#include <boost/random.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/distributions/beta.hpp>
#include <sys/time.h>

#include <Teuchos_ParameterList.hpp>

#include "math_tools.hpp"

namespace Surrogates {

/**
 * \class RandomNumberGenerator
 * \brief Easy to use wrapper to the Boost psuedo random number generator
 */
class RandomNumberGenerator
{
  public:
    /**
     * Default constructor
     */
    RandomNumberGenerator() { }

    /**
     * Deconstructor
     */
    ~RandomNumberGenerator() { }

    /**
     * Public interface
     */
    static void generate( IntMatrix& A, int M, int N, const Teuchos::ParameterList & params );
    static void generate( RealMatrix& A, int M, int N, const Teuchos::ParameterList & params );

  protected:

    static void uniform( IntMatrix& A, int M, int N, unsigned int seed );

    /**
     * \brief Generate a matrix of values drawn from the uniform distribution
     * on [0,1].
     */
    static void uniform( RealMatrix& A, int M, int N, unsigned int seed );

    /**
     * \brief Generate a matrix of values drawn from the Gaussian distribution
     *  N(mu,sigma).
     */
    static void gaussian( RealMatrix& A, int M, int N, Real mean, Real variance,
        unsigned int seed );

    /**
     * \brief Generate a matrix of values drawn from the Beta distribution
     *  B(alpha,beat).
     */
    static void beta( RealMatrix& A, int M, int N, Real alpha, Real beta,
        unsigned int seed );

    /**
     * \brief Generate N random permutations of the vector [1,2,...,M].
     */
    static void permutation( IntMatrix &result, int M, int N, unsigned int seed );
};

}; // namespace Surrogates

#endif

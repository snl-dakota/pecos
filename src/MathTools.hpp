/**
 * @file MathTools.hpp
 * @author John D. Jakeman
 * @date 31 October 2011
 * @brief Miscelaneous math functions.
 */

#ifndef MATH_TOOLS_HPP
#define MATH_TOOLS_HPP

#include "LinearAlgebra.hpp"

namespace Pecos {

/**
 * \brief Map a linear index of a d-dimensional array to the equivalent
 * d-dimensional index
 *
 * Example:
 \f[
 \begin{bmatrix}
  1 & 4 & 7\\
  2 & 5 & 8\\
  3 & 6 & 9
  \end{bmatrix}
  \rightarrow
  \begin{bmatrix}
  1,1 & 1,2 & 1,3\\
  2,1 & 2,2 & 2,3\\
  3,1 & 3,2 & 3,3
  \end{bmatrix}
 \f]
 *
 * @param sizes 
 * @param index
 * @param numElems
 * @return multiIndex the d-dimensional index
 */
void ind2sub( const IntVector &sizes, int index, int numElems,
	      IntVector &multiIndex );


/**
 * \brief Compute the cartesian product of an arbitray number of sets. 
 * 
 * These sets can consist of numbers of be themselves sets of vectors
 * @param inputSets the sets to be used in the cartesian product
 * @param elemSize the size of the vectors within each set.
 * @return the cartesian product
 */
template<typename O, typename T>
void cartesian_product( const std::vector< Teuchos::SerialDenseVector<O,T> > &inputSets, Teuchos::SerialDenseMatrix<O,T> &result, int elemSize  = 1)
{
  #ifdef DEBUG
  if ( inputSets.size() == 0 ) 
    {
      throw(std::invalid_argument("CartesianProduct() there must be at least one input set"));
    }
  #endif
  int numElems = 1;
  int nSets = inputSets.size();
  IntVector sizes;
  sizes.resize( nSets );
  IntVector multiIndex;
  for ( int i = 0; i < nSets; i++ )
    {
      sizes[i] = inputSets[i].length() / elemSize;
      numElems *= sizes[i];
    }
  result.reshape( nSets*elemSize, numElems );
  for ( int i = 0; i < numElems; i++ )
    {
      ind2sub( sizes, i, numElems, multiIndex );
      for ( int j = 0; j < nSets; j++ )
	{
	  for ( int k = 0; k < elemSize; k++ )
	    {
	      result( j*elemSize+k, i ) = inputSets[j][multiIndex[j]*elemSize+k];
	    }
	}
    }
};

/**
 * \brief Construct the outer product of an arbitray number of sets.
 *
 * Example: 
 \f[ \{1,2\}\times\{3,4\}=\{1\times3, 2\times3, 1\times4, 2\times4\} = 
 \{3, 6, 4, 8\} \f]
 * @param inputSets "Vector.hpp" of sets to be used in the outer product
 * @return result the outer product
*/
template<typename O, typename T>
void outerProduct( const std::vector< Teuchos::SerialDenseVector<O,T> > &inputSets, Teuchos::SerialDenseVector<O,T> &result )
{
  int numElems = 1;
  int nSets = inputSets.size();
  IntVector sizes;
  IntVector multiIndex;
  sizes.resize(nSets);
  for ( int i = 0; i < nSets; i++ )
    {
      sizes[i] = inputSets[i].length();
      numElems *= sizes[i];
    }
  result.resize( numElems );
  for ( int i = 0; i < numElems; i++ )
    {
      result[i] = 1.0;
      ind2sub( sizes, i, numElems, multiIndex );
      for ( int j = 0; j < nSets; j++ )
	{
	  result[i] *= inputSets[j][multiIndex[j]];
	}
    }
};

/**
 * \brief Build a vector containing all integers in [p,p+n) seperated by incr
 *
 * E.g Range( 10, 1, 2 ) -> [1,3,5,7,9]
 */
template<typename O, typename T>
void range( Teuchos::SerialDenseVector<O,T> &v, O m, O n, T incr  = 1 )
{
  int i = m, j = 0;
  v.size( ( n - m + 1 ) / (O)incr );
  while( i < n )
    {
      v[j] = i;
      i += incr;
      j++;
    }
  v.resize( j );
};

/** 
 * \brief Generate a vector whose n elements are equidistantly spaced in the
 * interval [a,b].
 *
 * @param [a,b] the range
 * @param n the number of times the interval [a,b] is divided
 * @return vec n numbers equidistant is [a,b]
 */
void linspace( RealVector &vec, Real a, Real b, int n );

/**
 * \brief Round a Real to the nearest integer.
 */
Real round( Real r );

/*
 * \brief Function to allow for searching for maximum absolute 
 * value (Real) in a containter
 *
 * @param a first value
 * @param b second value
 * @return true if a < b, false otherwise
 */
bool absCompare( Real a, Real b );

/** 
 * \brief Return the nearest integer to x
 * @param x the value
 * return int the neartest integer to x
 */
int nearest_integer( Real x );

/**
 * \brief Calulate n! 
 *
 * @param n the factorial of interest
 * @return n!
 */
Real factorial( int n );

/**
 * \brief Return the number of possible combinations of k elements that can be
 * chosen from n elements.
 *
 *\f[ { n \choose k} = \frac{n!}{k!(n-k)!}\f]
 * @param n number of elements
 * @param k the number of elements to be chosen. k<=n. 
 * @return the total number of combinations.
 */
int nchoosek( int n, int k );

/**
 * \brief Transform a point x in the interval [a,b] onto the interval [0,1] and
 * vice versa.
 *
 * @param x the value to be transformed
 * @param [a,b] the new/old interval
 * @param dir specifies whether to map to [a,b] or from [a,b]. That is 
 * dir=0  returns [0,1] -> [a,b]
 * dir!=0 returns [a,b]  -> [0,1]
 * @return the rescaled value
*/
Real rescale( Real x, Real a, Real b, int dir );

/**
 * \brief Transform a point x in the interval 
 * \f$[a_1,b_1] \times \cdots \times [a_d,b_d] \f$ onto 
 * the hypercube interval \f$[0,1]^d\f$ and vice versa.
 *
 * @param x the value to be transformed
 * @param [a,b] the new/old interval
 * @param dir specifies whether to map to [a,b] or from [a,b]. That is 
 * dir=0  returns [0,1] -> [a,b]
 * dir!=0 returns [a,b]  -> [0,1]
 * @return the rescaled value
*/
void rescale( RealMatrix &x, const RealVector &domain, int dir );

/**
 * \brief Construct a d-dimensional mesh on a hyper-rectangle. The meshes are 
 * equidistant with respect to each coordinate direction.
 *
 * @param numPts1D array specifying the number of meshpoints for each dimension
 * @param domain the min and max value of each dimension. 
 *        That is \f$ [a_1,b_1,...,a_d,b_d]\f$
 * @return the multi-dimensional mesh coordinates. ( num_dims x nPts ) array. 
 * nPts = numPts1D[0]*...*numPts1D[d-1]
 */
void mesh_grid( const IntVector &numPts1D, 
	       const RealVector &domain, 
	       RealMatrix &points );

/**
 * Get all the multi-dimensional indices of a multi-dimensional polynomial
 * basis with a specified degree sum. A cleaner interface is provided by
 * GetMultiDimensionalPolynomial_indices()
 * @return B the multi-dimensional indices
 */
void fill_seq( IntMatrix &B, int* elems, int len_elems, int* pos, int len_pos, 
	       int choices_made, int first_pos, int order, int &row );

/**
 * \brief Get all the multi-dimensional indices of a multi-dimensional polynomial
 * basis with a specified degree sum.
 *
 * @param num_dims the dimensionailty
 * @param degree the degree sum of the polynomial indices wanted
 * @return B the multi-dimensional polynomial_indices
 */
void get_multi_dimensional_polynomial_indices( int num_dims, int degree, IntMatrix &B );


/**
 * \brief Get the multi-dimensional hypercube \f$[a,b]^d\f$
 *
 * @param num_dims the dimension of the computational domain.
 * @param \f$[a,b]\f$ the one-dimensional bounds of the hypercube.
 * @return the multi-dimensional bounds of the hypercube.
 */
void set_hypercube_domain( RealVector &domain, int num_dims, Real a, Real b );

/// Perturb the columns of a matrix
void get_permutations( IntMatrix &permutations, 
		       int M , int N, unsigned int seed );

enum lp_norm
  {
    l1_norm,
    l2_norm,
    linf_norm
  };
  
/**
 * \brief Return the index of the element of x with the minimum value.
 *
 *\f[
     \arg \! \min_i x_i\quad i=1,\ldots,n
  \f] 
 */
template<typename T>
int argmin( int n, T* x )
{
  int amin( 0 );
  T min = x[0];
  for ( int i = 1; i < n; i++ )
    {
      if ( x[i] < min )
	{
	  min = x[i];
	  amin = i;
	}
    }
  return amin;
};

/**
 * \brief Return the index of the element of x with the maximum value.
 *
 *\f[
     \arg \! \max_i x_i \quad i=1,\ldots,n
  \f] 
 */
template<typename T>
int argmax( int n, T* x )
{
  int amax( 0 );
  T max = x[0];
  for ( int i = 1; i < n; i++ )
    {
      if ( x[i] > max )
	{
	  max = x[i];
	  amax = i;
	}
    }
  return amax;
};

/**
 * \brief Return the sum of the elements in the array x. 
 *
 *\f[
     \sum_{i=1}^N x_i
  \f]
 */
template<typename T>
T sum( int n, T* x )
{
  T sum = x[0];
  for ( int i = 1; i < n; i++ )
    sum += x[i];
  return sum;
};

/**
 * \brief Return the sum of the elements in the vector x. 
 *
 *\f[
     \sum_{i=1}^N x_i
  \f]
  *
  * This will not return the same result as sum( int n, T* x ) if
  * the vector is a subvector of another vector. That is stride does not equal
  * the number of rows.
 */
template < typename O, typename T >
T sum( Teuchos::SerialDenseVector<O,T> &v )
{
  T tmp = Teuchos::ScalarTraits<T>::zero();
  for ( O i = 0; i < v.length(); i++ )
    tmp += v[i];
  return tmp;
};

/**
 * \brief Return the mean of the array x. 
 *
 * \f[
    \bar{x} = \frac{1}{N}\sum_{i=1}^N x_i
  \f]
 */
Real mean( int n, Real *x );


/**
 * \brief Return the sample variance of the array x. 
 *
 * \f[
    \sigma^2 = \frac{1}{N-\text{ddof}}\sum_{i=1}^N ( x_i-\bar{x} )^2
  \f]
  * 
  * \param ddof Delta Degrees of Freedom (ddof).  
  * The divisor used in calculations is \f$N - \text{ddof}\f$, 
  * where \$fN\$f represents the number of elements.
  * By default ddof is one.
 */
Real variance( int n, Real *x, int ddof = 1 );


/**
 * \brief Calculate one of several different types of \f$\ell_p\f$ error norms
 * of a set of vectors. 
 *
 * Each column of the referenceValues matrix is considered 
 * independently with the correpsonding column in the approximationValues.
 * The error between the each set of two columns tested is returned.
 * \param referenceValues a matrix containing the reference values
 * \param approximateValues a matrix containig the approximate values
 * \param activeColumns specifies which set of columns to consider. If empty
 *        all columns are considered
 * \pram normalise if true the error norms are normalised. The l_inf norm
 *                 is normalised by the half the range of the data in the 
 *                 referenceValues column. The L_2 norm is normalised by the
 *                 standard deviation of the data in the reference values column.
 */
void lp_error( RealMatrix &referenceValues, 
	      RealMatrix &approximateValues,
	      std::vector<lp_norm> error_norms, RealMatrix &error, 
	      IntVector &activeColumns, bool normalise = false );

/**
 *\brief Construct a Latin Hyperube Design (LHD)
 */
void latin_hypercube_design( int nPts, int num_dims, RealMatrix &points, int seed );


/**
 * \breif Return  the sign of x.
 *
 * \return 1 if the corresponding element of x is greater than zero;
 * 0 if the corresponding element of X equals zero;
 *-1 if the corresponding element of X is less than zero
 */
template <typename T> 
int sgn( T x )
{
  return ( T(0) < x ) - ( x < T(0) );
}

template <typename T> 
int num_non_zeros( T *data, int n )
{
  int num_non_zero_count( 0 );
  for ( int i = 0; i < n; i++ )
    {
      if ( std::abs( data[i] ) > 0 ) 
	num_non_zero_count++;
    }
  return num_non_zero_count;
}

} // namespace Pecos

#endif

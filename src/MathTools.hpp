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
 * \brief Map a scalar index of a flat 1D array to the equivalent
 *        d-dimensional index
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
 * @param sizes the number of elems in each dimension. For a 2D index
 * sizes = [numRows, numCols]
 * @param index the scalar index
 * @param num_elems the total number of elements in the d-dimensional matrix
 * @return result the d-dimensional index
 */
void ind2sub( const IntVector &sizes, int index, int num_elems,
	      IntVector &result );

/**
 * \brief Map a d-dimensional index to the scalar index of the equivalent flat 
 * 1D array
 * 
 * Example:
 \f[
 \begin{bmatrix}
 1,1 & 1,2 & 1,3\\
 2,1 & 2,2 & 2,3\\
 3,1 & 3,2 & 3,3
 \end{bmatrix}
 \rightarrow
 \begin{bmatrix}
 1 & 4 & 7\\
 2 & 5 & 8\\
 3 & 6 & 9
 \end{bmatrix}
 \f]
 *
 * @param sizes the number of elems in each dimension. For a 2D index
 * sizes = [numRows, numCols]
 * @param multi_index the d-dimensional index
 * @return scalar_index the scalar index of the flat array
 */
int sub2ind( const IntVector &sizes, const IntVector &multi_index );

/**
 * \brief Compute the cartesian product of an arbitray number of sets. 
 * 
 * These sets can consist of numbers of be themselves sets of vectors
 * @param inputSets the sets to be used in the cartesian product
 * @param elem_size the size of the vectors within each set.
 * @return the cartesian product
 */
template<typename O, typename T>
void cartesian_product( const std::vector< Teuchos::SerialDenseVector<O,T> > &input_sets, Teuchos::SerialDenseMatrix<O,T> &result, int elem_size )
{
#ifdef DEBUG
  if ( input_sets.size() == 0 ) 
    {
      throw(std::invalid_argument("CartesianProduct() there must be at least one input set"));
    }
#endif
  int num_elems = 1;
  int num_sets = input_sets.size();
  IntVector sizes;
  sizes.resize( num_sets );
  IntVector multi_index;
  for ( int i = 0; i < num_sets; i++ )
    {
      sizes[i] = input_sets[i].length() / elem_size;
      num_elems *= sizes[i];
    }
  result.reshape( num_sets*elem_size, num_elems );
  for ( int i = 0; i < num_elems; i++ )
    {
      ind2sub( sizes, i, num_elems, multi_index );
      for ( int j = 0; j < num_sets; j++ )
	{
	  for ( int k = 0; k < elem_size; k++ )
	    {
	      result( j*elem_size+k, i ) = 
		input_sets[j][multi_index[j]*elem_size+k];
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
 * @param input_sets "Vector.hpp" of sets to be used in the outer product
 * @return result the outer product
 */
template<typename O, typename T>
void outer_product( const std::vector< Teuchos::SerialDenseVector<O,T> > &input_sets, Teuchos::SerialDenseVector<O,T> &result )
{
  int num_elems = 1;
  int num_sets = input_sets.size();
  IntVector sizes( num_sets, false );
  for ( int i = 0; i < num_sets; i++ )
    {
      sizes[i] = input_sets[i].numRows();
      num_elems *= sizes[i];
    }
  IntVector multi_index( num_sets, false );
  result.sizeUninitialized( num_elems );
  for ( int i = 0; i < num_elems; i++ )
    {
      result[i] = 1.0;
      ind2sub( sizes, i, num_elems, multi_index );
      for ( int j = 0; j < num_sets; j++ )
	result[i] *= input_sets[j][multi_index[j]];
    }
}

/**
 * \brief Build a vector containing all integers in [m,n) seperated by incr
 *
 * E.g Range( 1, 10, 2 ) -> [1,3,5,7,9]
 */
template<typename O, typename T>
void range( Teuchos::SerialDenseVector<O,T> &result, T m, T n, T incr = 1 )
{
  T i = m;
  int j = 0;
  int len = (int)( n - m ) / (int)incr;
  if ( ( (int)( n - m ) % (int)incr ) != 0 ) len++;
  result.size( len );
  while( i < n )
    {
      result[j] = i;
      i += incr;
      j++;
    }
};

/** 
 * \brief Generate a vector whose n elements are equidistantly spaced in the
 * interval [a,b].
 *
 * @param [a,b] the range
 * @param n the number of times the interval [a,b] is divided
 * @return vec n numbers equidistant is [a,b]
 */
void linspace( RealVector &result, Real a, Real b, int n );

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
bool abs_compare( Real a, Real b );

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
 *\brief transform points in domain_in to points in domain_out
 */
void hypercube_map( RealMatrix &x, const RealVector &domain_in,  
		    const RealVector &domain_out,
		    RealMatrix &result );

/**
 * \brief Construct a d-dimensional mesh on a hyper-rectangle. The meshes are 
 * equidistant with respect to each coordinate direction.
 *
 * @param num_pts_1d array specifying the number of meshpoints for each dimension
 * @param domain the min and max value of each dimension. 
 *        That is \f$ [a_1,b_1,...,a_d,b_d]\f$
 * @return the multi-dimensional mesh coordinates. ( num_dims x num_pts ) array. 
 * num_pts = num_pts_1d[0]*...*num_pts_1d[d-1]
 */
void mesh_grid( const IntVector &num_pts_1d, 
		const RealVector &domain, 
		RealMatrix &result );

/**
 * Get all the multi-dimensional indices of a multi-dimensional polynomial
 * basis with a specified degree sum. A cleaner interface is provided by
 * GetMultiDimensionalPolynomial_indices()
 * @return B the multi-dimensional indices
 */
void get_multi_dimensional_polynomial_subspace_indices( IntMatrix &B, 
							int* elems, 
							int len_elems, 
							int* pos, 
							int len_pos, 
							int choices_made, 
							int first_pos, 
							int order, int &row );

/**
 * \brief Get all the multi-dimensional indices of a multi-dimensional polynomial
 * basis with a specified degree sum.
 *
 * @param num_dims the dimensionailty
 * @param degree the degree sum of the polynomial indices wanted
 * @return B the multi-dimensional polynomial_indices
 */
void get_multi_dimensional_polynomial_indices( int num_dims, int degree, 
					       IntMatrix &result );


/**
 * \brief Get the multi-dimensional hypercube \f$[a,b]^d\f$
 *
 * @param num_dims the dimension of the computational domain.
 * @param \f$[a,b]\f$ the one-dimensional bounds of the hypercube.
 * @return result the multi-dimensional bounds of the hypercube.
 */
void set_hypercube_domain( RealVector &result, int num_dims, Real a, Real b );

/// Perturb the columns of a matrix
void get_permutations( IntMatrix &permutations, 
		       int M , int N, unsigned int seed );

void compute_hyperbolic_level_subdim_indices( int num_dims, int level, 
					      int num_active_dims, 
					      Real p,
					      IntMatrix &result );

void compute_hyperbolic_level_indices( int num_dims, int level, Real p, 
				       IntMatrix &result );

void compute_hyperbolic_indices( int num_dims, int level, Real p, 
				 IntMatrix &result );

// weights ( num_dims x 1 ) vector with values in [0,1]
void compute_anisotropic_hyperbolic_indices( int num_dims, int level, Real p, 
					     RealVector &weights,
					     IntMatrix &result );

enum lp_norm
  {
    l1_norm,
    l2_norm,
    linf_norm,
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
 * Each column of the reference_values matrix is considered 
 * independently with the correpsonding column in the approximation_values.
 * The error between the each set of two columns tested is returned.
 * \param reference_values a matrix containing the reference values
 * \param approximate_values a matrix containig the approximate values
 * \param active_columns specifies which set of columns to consider. If empty
 *        all columns are considered
 * \pram normalise if true the error norms are normalised. The l_inf norm
 *                 is normalised by the half the range of the data in the 
 *                 reference_values column. The L_2 norm is normalised by the
 *                 standard deviation of the data in the reference values column.
 */
void lp_error( RealMatrix &reference_values, 
	       RealMatrix &approximate_values,
	       std::vector<lp_norm> error_norms, RealMatrix &error, 
	       IntVector &active_columns, bool normalise = false );

/**
 *\brief Construct a Latin Hyperube Design (LHD)
 */
void latin_hypercube_design( int num_pts, int num_dims, RealMatrix &result, 
			     int seed );


void compute_next_combination( int num_dims, int level, IntVector &index, 
			       bool &extend, int &h, int &t );

void compute_combinations( int num_dims, int level, IntMatrix &result );


/**
 * \breif Return  the sign of x.
 *
 * \return 1 if the corresponding element of x is greater than zero;
 * 0 if the corresponding element of X equals zero;
 *-1 if the corresponding element of X is less than zero
 */
template<typename T> 
int sgn( T x )
{
  return ( T(0) < x ) - ( x < T(0) );
}

template<typename T> 
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

/**
 *\brief return the median of a std::vector
 */
template<typename T>
double median( std::vector<T> &v )
{
  size_t n = v.size();
  typename std::vector<T>::iterator target = v.begin() + n / 2;
  std::nth_element( v.begin(), target, v.end() );
  if( n % 2 != 0 )
    {
      return (double)*target;
    }
  else
    {
      std::vector<Real>::iterator target_neighbour = 
	std::max_element( v.begin(), target );
      return (double)( *target + *target_neighbour ) / 2.0;
    }
};

/**
 *\brief return the median of a RealVector.
 */
template<typename O, typename T>
double median( Teuchos::SerialDenseVector<O,T> &v )
{
  std::vector<T> tmp( v.length() );
  for ( int i = 0; i < v.length(); i++ )
    tmp[i] = v[i];
  return median( tmp );
}

template<typename O, typename T>
double pnorm( Teuchos::SerialDenseVector<O,T> &v, double p )
{
  double sum = 0.;
  for ( O i = 0; i < v.length(); i++ )
    {
      sum += std::pow( std::abs( (double)v[i] ), p );
    }
  return std::pow( sum, 1./p );
}

template<typename O, typename T>
O num_nonzeros( Teuchos::SerialDenseVector<O,T> &v )
{
  O num_nonzeros = Teuchos::ScalarTraits<O>::zero();
  for ( int d = 0; d < v.length(); d++ )
    {
      if ( v[d] != Teuchos::ScalarTraits<T>::zero() ) num_nonzeros++;
    }
  return num_nonzeros;
}

template<typename O, typename T>
void nonzero( Teuchos::SerialDenseVector<O,T> &v, 
	      IntVector &result )
{
  int num_nonzeros = 0;
  result.sizeUninitialized( v.length() );
  for ( int d = 0; d < v.length(); d++ )
    {
      if ( v[d] != Teuchos::ScalarTraits<T>::zero() )  
	{
	  result[num_nonzeros] = d;
	  num_nonzeros++;
	}
    }
  // chop off unwanted memory
  result.resize( num_nonzeros );
}

void lhs_indices( int num_dims, int num_pts, int duplication, int seed, 
		  IntMatrix &result );

/// Construct the improved distributed hypercube sample
void ilhs( int num_dims, int num_pts, int duplication, int seed, 
	   RealMatrix &result );

/// Construct basic random hypercube sample
void lhs( int num_dims, int num_pts, int seed, 
	  RealMatrix &result );

template<typename T>
class index_sorter {
public:
  index_sorter(T &values){ values_ = values; }
  bool operator()( int lhs, int rhs) const {
    return values_[lhs] < values_[rhs];
  }
private:
  T values_;
};

template<typename T>
class magnitude_index_sorter {
public:
  magnitude_index_sorter(T &values){ values_ = values; }
  bool operator()( int lhs, int rhs) const {
    return std::abs( values_[lhs] ) > std::abs( values_[rhs] );
  }
private:
  T values_;
};

// sorts in ascending order
template<typename O, typename T>
void argsort( Teuchos::SerialDenseVector<O,T> &values, 
	      IntVector &result )
{
  std::vector<O> indices( values.length() );
  for ( O i = 0; i < values.length(); i++ )
    indices[i] = i;
  
  std::sort( indices.begin(), indices.end(), index_sorter< Teuchos::SerialDenseVector<O,T> >( values ) );
	     
  result.sizeUninitialized( values.length() );
  for ( O i = 0; i < values.length(); i++ )
    result[i] = indices[i];
}

template<typename O, typename T>
void argsort( Teuchos::SerialDenseVector<O,T> &v,
	      Teuchos::SerialDenseVector<O,T> &result_0,
	      IntVector &result_1 )
{
  argsort( v, result_1 );
  result_0.sizeUninitialized( v.length() );
  for ( O i = 0; i < v.length(); i++ )
    result_0[i] = v[result_1[i]];
}

// sorts in descending order absolute value of entries
template<typename O, typename T>
void magnitude_argsort( Teuchos::SerialDenseVector<O,T> &values, 
			IntVector &result )
{
  std::vector<O> indices( values.length() );
  for ( O i = 0; i < values.length(); i++ )
    indices[i] = i;
  
  std::sort( indices.begin(), indices.end(), magnitude_index_sorter< Teuchos::SerialDenseVector<O,T> >( values ) );
	     
  result.sizeUninitialized( values.length() );
  for ( O i = 0; i < values.length(); i++ )
    result[i] = indices[i];
}

/**
 * \brief find the interval containing a target value. Assumes data is 
 * in ascending order
 */
template<typename O, typename T>
int binary_search( T target, Teuchos::SerialDenseVector<O,T> &data )
{
  O low = 0, high = data.length()-1, mid;
  while ( low <= high )
    {
      mid = low + ( high - low ) / 2;
      if ( target < data[mid] ) high = mid - 1;
      else if ( target > data[mid] ) low = mid + 1;
      else return mid;
    }
  if ( high < 0 ) return 0;
  else if ( high < low ) return high;
  else return low;
}

class LinearInterpolant1D
{
private:
  RealVector pts_, vals_;

public:
  LinearInterpolant1D( RealVector &pts, RealVector &vals )
  {
    pts_ = pts; vals_ = vals;
  };

  ~LinearInterpolant1D(){};
  
  void interpolate( RealVector &pts, RealVector &result )
  {
    int num_pts = pts.length();
    if ( result.length() != num_pts )
      result.sizeUninitialized( num_pts );
    for ( int i = 0; i < num_pts; i++ )
      {
	// enforce constant interpolation when interpolation is outside the
	// range of pts_
	if ( pts[i] <= pts_[0] )
	  result[i] = vals_[0];
	else if ( pts[i] >= pts_[pts_.length()-1] )
	  result[i] = vals_[pts_.length()-1];
	else
	  {
	    // assumes binary search returns index of the closest point in pts_ 
	    // to the left of pts[i]
	    int index = binary_search( pts[i], pts_ );
	    result[i] = vals_[index] + ( ( vals_[index+1] - vals_[index] ) / 
					 ( pts_[index+1] - pts_[index] ) ) 
	      * ( pts[i] - pts_[index] );
	  }
      }
  };
};

/**
 * \brief Reverse the contents of a vector. 
 *
 * Useful for changing an ordered array to/from ascending/descending order
 */
template<typename O, typename T>
void reverse( Teuchos::SerialDenseVector<O,T> &v )
{
  O n = v.length();
  Teuchos::SerialDenseVector<O,T> tmp( v );
  for ( O i = 0; i < n; i++ )
    v[i] = tmp[n-i-1];
}

}  // namespace Pecos

#endif

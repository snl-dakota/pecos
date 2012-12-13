#include "MathTools.hpp"

namespace Pecos {

void ind2sub( const IntVector &sizes, int index, int numElems,
	      IntVector &multiIndex )
/* Map a linear index of a d-dimensional array to the equivalent
   d-dimensional index
   Example:
   | 1  4  7 |      | 1,1  1,2  1,3 |
   | 2  5  8 |  --> | 2,1  2,2  2,3 |
   | 3  6  9 |      | 3,1  3,2  3,3 |
*/
{

  int denom = numElems;
  int nSets = sizes.length();
  multiIndex.resize(nSets);
  for ( int i = nSets-1; i >= 0; i-- )
    {
      denom /= sizes[i];
      multiIndex[i] = index / denom;
      index = index % denom;
    }
};

void linspace( RealVector &vec, Real a, Real b, int n )
/* Generate a vector whose n elements are equdistantly distributed in the
   interval [a,b].
*/
{
  vec.resize( n );
  Real h = ( b - a ) / (Real)(n-1);
  for ( int i = 0; i < n ; i++ )
    {
      vec[i] = a + (Real)i * h;
    }
};

Real round( Real r ) 
{
  return (r > 0.0) ? std::floor(r + 0.5) : std::ceil(r - 0.5);
};

bool absCompare( Real a, Real b )
{
  return ( std::fabs( a ) < std::fabs( b ) );
};

int nearest_integer( Real x )
{
  int value;
  
  if ( x < 0.0 )
    {
      value = - ( int ) ( std::fabs ( x ) + 0.5 );
    }
  else
    {
      value =   ( int ) ( std::fabs ( x ) + 0.5 );
    }

  return value;
};

Real factorial(int n)
/* Calulate n! */
{
  Real temp = 1.0;
  for ( int i = 1; i <= n; i++ )
    {
      temp *= (Real)(i);
    }
  return temp;
};

int nchoosek(int n, int k)
/* Determine the number of possible combinations of k elements that can be
   chosen from n elements. */
{
  Real value = 1;
  for ( int i = 0; i < n-k; i++ )
    {
      value *= (Real)(n-i) / (Real)(n-k-i);
    }
  return (int)round( value );
  //return (int)(Factorial(n)/(Factorial(n-k)*Factorial(k)));
};

Real rescale( Real x, Real a, Real b, int dir )
{
  if ( dir == 0 )
    {
      return a + ( b-a ) * x; // [0,1] -> [a,b]
    }
  else
    {
      return ( x - a ) / ( b - a); // [a,b] -> [0,1]
    }
};

void rescale( RealMatrix &x, const RealVector &domain, int dir )
{
  for ( int i = 0; i < x.numCols(); i++ )
    {
      for ( int d = 0; d < x.numRows(); d++ )
	{
	  x(d,i) = rescale( x(d,i), domain[2*d], domain[2*d+1], dir );
	}
    }
};

void mesh_grid( const IntVector &numPts1D, 
	       const RealVector &domain, 
	       RealMatrix &points )
{
#ifdef DEBUG
  if ( numPts1D.length() % 2 != 0 )
    {
      
      throw(std::invalid_argument("Ensure numPts1D.length() %2 == 0."));
    }
#endif

  int num_dims( domain.length() / 2 );

  std::vector< RealVector > pointSets1D;
  pointSets1D.resize( num_dims );  
  for ( int d = 0; d < num_dims; d++ )
    {
      linspace( pointSets1D[d], domain[2*d], domain[2*d+1], numPts1D[d] );
    }
  cartesian_product( pointSets1D, points );
};

void fill_seq( IntMatrix &B, int* elems, int len_elems, int* pos, int len_pos, 
	      int choices_made, int first_pos, int order, int &col )
{
  int start,temp;
  // Have we selected the number of required elements?
  if ( choices_made >= len_pos )
    {
      start = 0;
      temp = 0;
      for ( int j = 0; j < len_pos; j++ )
	{
	  if ( pos[j]-start == 0 )
	    {
	      B(j,col) = 0;
	    }
	  else
	    {
	      B(j,col) = pos[j] - start;
	      temp = temp + ( pos[j] - start );
	    }
	  start = pos[j] + 1;
	}
      B(len_pos,col) = order - temp;
      col++;
      return;
    }
  
  // Are there enough remaining elements to be selected?
  if ( ( len_elems - first_pos ) < ( len_pos - choices_made ) )
    {
      return;
    }

  // Try to select new elements in the same position or
  // to the right of the last selected one.
  for ( int j = first_pos; j < len_elems; j++ )
    {
      pos[choices_made] = j;
      fill_seq( B, elems, len_elems, pos, len_pos, choices_made+1, j+1, order, 
	       col);
      //change to ii if items can be selected repeatedly: Do Not for gPC
    }
  return;
};

void get_multi_dimensional_polynomial_indices( int num_dims, int degree, IntMatrix &B )
{
  // m: degree of polynomial m=0,...,degree
  // i: index of basis function i=0,...,P-1
  // j: index of random dimension j=0,...,num_dims-1
  
  int len_elems, len_pos, choices_made,first_pos;
  int *elems,*pos;
  int P = nchoosek( num_dims + degree - 1, num_dims - 1 );
  if ( ( B.numRows() != num_dims ) || ( B.numCols() != P ) )
    {
      std::string msg = "get_multi_dimensional_polynomial_indices() ";
      msg += "The size of B must be set before entry to function.";
      throw( std::runtime_error( msg ) );
    }

  int col = 0;
  len_pos = num_dims - 1;
  pos = new int [len_pos];
  
  // Generate all multimonomial indices of degree degree
  len_elems = ( num_dims - 1 ) + degree;
  elems = new int [len_elems];

  for ( int i = 0; i < len_elems; i++ )
    {
      elems[i] = i;
      if ( i < len_pos )
	{
	  pos[i] = 0;
	}
    }
  choices_made = 0;
  first_pos = 0;
  fill_seq( B, elems, len_elems, pos, len_pos, choices_made, first_pos, degree, 
	   col );
  delete[] elems;
  delete[] pos;
};

void set_hypercube_domain( RealVector &domain, int num_dims, Real a, Real b )
{
  domain.resize( 2*num_dims );
  for ( int d = 0; d < num_dims; d++ )
    {
      domain[2*d] = a;
      domain[2*d+1] = b;
    }
};

void lp_error( RealMatrix &referenceValues, RealMatrix &approximateValues,
	      std::vector<lp_norm> error_norms, RealMatrix &error, 
	      IntVector &activeColumns, bool normalise )
{
  Teuchos::BLAS<int, Real> blas;

  if ( referenceValues.numRows() != approximateValues.numRows() )
    throw( std::runtime_error( "lp_error() Matrix sizes do not match." ) );

  if ( activeColumns.length() == 0 )
    range( activeColumns, 0, approximateValues.numCols() );
 
  RealMatrix diff( referenceValues );
  diff -= approximateValues;

  int M( diff.numRows() );

  reshape( error, activeColumns.length(), (int)error_norms.size() );

  for ( int j = 0; j < (int)error_norms.size(); j++ )
    {
      switch ( error_norms[j] )
	{
	case linf_norm:
	  {
	    // Infinity norm
	    for ( int i = 0; i < activeColumns.length(); i++ )
	      {
		// Access the active columns
		int colNumber( activeColumns[i] );
		Real* diffCol = diff[colNumber];
		// Find the index of the element of the col with the 
		// maximum magnitude. 
		int amax = blas.IAMAX( M, diffCol, 1 ) - 1;
		error(i,j) = std::abs(diffCol[amax]);
		if ( normalise )
		  {
		    // Normalise by the half the range of the data in
		    // col
		    Real minima, maxima;
		    Real *refCol = referenceValues[colNumber];
		    minima = refCol[argmin( M, refCol )];
		    maxima = refCol[argmax( M, refCol )];
		    error(i,j) /= ( 0.5 * ( maxima - minima ) );
		  }
	      }
	    break;
	  }
	case l1_norm:
	  {
	    // L1 norm
	    for ( int i = 0; i < activeColumns.length(); i++ )
	      {
		// Access the active column
		int colNumber( activeColumns[i] );
		Real* diffCol = diff[colNumber];
		//Sum the absolute values of the entries of col. 
		error(i,j) = blas.ASUM( M, diffCol, 1 ) / (Real)M;
	      }
	    break;
	  }
	case l2_norm:
	  {
	    // L2 norm
	    for ( int i = 0; i < activeColumns.length(); i++ )
	      {
		// Access the active column
		int colNumber( activeColumns[i] );
		Real *diffCol = diff[colNumber];
		//Sum the absolute values of the entries of col. 
		error(i,j) = blas.NRM2( M, diffCol, 1 ) / std::sqrt( (Real)M );
		if ( normalise )
		  {
		    // Normalise by the standard deviation of the data in
		    // col
		    Real *refCol = referenceValues[colNumber];
		    Real stddev = std::sqrt( variance( M, refCol ) );
		    error(i,j) /= stddev;
		  }
	      }
	    break;
	  }
	}
    }
};

Real mean( int n, Real *x )
{
  return sum( n, x ) / (Real)n;
};

Real variance( int n, Real *x, int dof )
{
  Real mu = mean( n, x );

  Real variance( 0.0 );
  for ( int i = 0; i < n; i++ )
    {
      variance += ( x[i] - mu ) * ( x[i] - mu );
    }
  return variance / (Real)( n - dof );
}

void get_permutations( IntMatrix &permutations, 
		       int M , int N, unsigned int seed )
{
  std::srand( seed );
  
  permutations.reshape( M, N );
  IntMatrix numbers;
  for ( int j = 0; j < N; j++ )
    {
      std::vector<int> random( M );
      for ( int i = 0; i < M; i++ ) 
	{
	  random[i] = i;
	}

      std::random_shuffle( random.begin(), random.end() );

      for ( int i = 0; i < M; i++ )
	{
	  permutations(i,j) = random[i];
	}
    }
};

} // namespace Pecos

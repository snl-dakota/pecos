#include "tensor_product_basis.hpp"

TensorProductBasis::TensorProductBasis() :
  numDims_(0)
{};
  
TensorProductBasis::~TensorProductBasis()
{clear();};

void TensorProductBasis::clear()
{
  bases1D_.clear();
  numDims_ = 0;
  dimensionList_.clear();
};

void TensorProductBasis::set_num_dims( int num_dims )
{
  numDims_ = num_dims; 
};

/**
 * \brief Initialise a set of one-dimensional bases.
 * dimension list has num unique bases as rows
 * each row has variable length. The entries in each row are the dimensions
 * that use that unique basis. E.g.
 * Legendre 2 4
 * Jacobi 0 1 3
 * means that pce is [Jacobi, Jacobi, Legendre, Jacobi, Legendre].
 * this can produce significant speedups when evaluating the tensor
 * product basis
 */
void TensorProductBasis::set_bases( const std::vector< OrthogPoly_ptr >& bases, 
		const std::vector< IntVector >& dimension_list, 
		int set_size )
{
  if ( bases.size() != dimension_list.size() )
    {
      std::string msg = "TensorProductBasis::set_bases() ";
      msg += "Ensure each dimenion is assigned a basis";
      throw( std::runtime_error( msg ) );
    }
  bases1D_.resize( set_size );
  int count = 0;
  for ( int i = 0; i < (int)dimension_list.size(); i++ )
    {
      for ( int j = 0; j < (int)dimension_list[i].length(); j++ )
	{
	  bases1D_[dimension_list[i][j]] = bases[i];
	  count++;
	}
    }
  if ( count != set_size )
    { 
      std::string msg = "TensorProductBasis::set_bases() ";
      msg += "The dimension list does not match the set size";
      throw( std::runtime_error( msg ) );
    }
  set_num_dims( set_size );
  dimensionList_ = dimension_list;
};

Real TensorProductBasis::value( const RealVector &pt,
				PolyIndex_ptr index ) const
{
  RealVector result;
  value_set( pt, index, result );
  return result[0];
}

void TensorProductBasis::value_set( const RealMatrix &pts,
				    PolyIndex_ptr index,
				    RealVector &result ) const
{
  int num_pts = pts.numCols();
  RealMatrix pts_trans( pts, Teuchos::TRANS );
  result.sizeUninitialized( num_pts );
  result = 1.;
  for ( int d = 0; d < index->effective_dimension(); d++ )
    {
      int dim = index->dimension( d );
      int degree = index->level( d );
      RealVector x( Teuchos::View, pts_trans[dim], num_pts );
      RealVector basis_1d;
      bases1D_[dim]->value_set( x, degree, basis_1d );
      for ( int j = 0; j < num_pts; j++ )
	result[j] *= basis_1d[j];
    }
};

void TensorProductBasis::
compute_max_degree_for_each_dimension( const std::vector<PolyIndex_ptr>& indices,
				       IntVector &max_degree_1d ) const{
  int num_indices = indices.size();
  max_degree_1d.size( numDims_ ); // initialize to zero
  for ( int i = 0; i < num_indices; i++ )
    {
      for ( int d = 0; d < indices[i]->effective_dimension(); d++ )
	{
	  int dim = indices[i]->dimension( d );
	  int degree = indices[i]->level( d );
	  max_degree_1d[dim] = std::max( max_degree_1d[dim], degree );
	}
    }
}

void TensorProductBasis::
compute_1d_basis_values_for_each_dimension( 
    const RealMatrix &pts,
    const IntVector &max_degree_1d, 
    std::vector< std::vector<RealVector> > &basis_values_1d ) const{
  int num_pts = pts.numCols();
  basis_values_1d.resize( numDims_ );
      
  for ( int dim = 0; dim < numDims_; dim++ ){
    basis_values_1d[dim].resize( num_pts );
    //RealVector x( Teuchos::View, pts_trans[dim], num_pts );
    for ( int k = 0; k < num_pts; k++ ){
      //bases1D_[dim]->value( x[k], max_degree_1d[dim], 
      bases1D_[dim]->value( pts(dim,k), max_degree_1d[dim], 
			    basis_values_1d[dim][k] );
    }
  }
};

void TensorProductBasis::
compute_1d_basis_derivatives_for_single_dimension( 
    const RealMatrix &pts,
    int derivative_dim,
    int max_degree,
    std::vector<RealVector> &basis_derivatives_1d ) const {
  int num_pts = pts.numCols();
  basis_derivatives_1d.resize( num_pts );
  //RealVector x( Teuchos::View, pts_trans[dim], num_pts );
  for ( int k = 0; k < num_pts; k++ ){
    //bases1D_[dim]->derivative( x[k], max_degree_1d[dim], 
    bases1D_[derivative_dim]->derivative( pts(derivative_dim,k), 
					  max_degree, 
					  basis_derivatives_1d[k] );
  }
}

void TensorProductBasis::value_set( const RealMatrix &pts,
				    std::vector<PolyIndex_ptr>& indices,
				    RealMatrix &result ) const
{
  int num_pts = pts.numCols();

  int num_indices = indices.size();
  result.shapeUninitialized( num_pts, num_indices );
  result = 1.;

  IntVector max_degree_1d;
  compute_max_degree_for_each_dimension( indices, max_degree_1d );

  std::vector< std::vector<RealVector> > basis_values_1d;
  compute_1d_basis_values_for_each_dimension( pts, max_degree_1d, 
					      basis_values_1d );

  for ( int i = 0; i < num_indices; i++ ){
    for ( int d = 0; d < indices[i]->effective_dimension(); d++ ){
      int dim = indices[i]->dimension( d );
      int degree = indices[i]->level( d );
      for ( int j = 0; j < num_pts; j++ )
	result(j,i) *= basis_values_1d[dim][j][degree];
    }
  }
};

Real TensorProductBasis::derivative( const RealVector &pt,
				     const PolyIndex_ptr index,
				     int derivative_dim ) const
{
  RealVector result;
  derivative_set( pt, index, derivative_dim, result );
  return result[0];
}

void TensorProductBasis::derivative_set( const RealMatrix &pts,
					 PolyIndex_ptr index,
					 int derivative_dim,
					 RealVector &result ) const
{
  int num_pts = pts.numCols();
  RealMatrix pts_trans( pts, Teuchos::TRANS );
  result.sizeUninitialized( num_pts );
  result = 1.;
  for ( int d = 0; d < index->effective_dimension(); d++ )
    {
      int dim = index->dimension( d );
      if ( dim != derivative_dim )
	{
	  int degree = index->level( d );
	  RealVector x( Teuchos::View, pts_trans[dim], num_pts );
	  RealVector basis_1d;
	  bases1D_[dim]->value_set( x, degree, basis_1d );
	  for ( int j = 0; j < num_pts; j++ )
	    result[j] *= basis_1d[j];
	}
    }
  
  int dim = derivative_dim;
  IntVector dimension_data;
  index->get( dim, dimension_data );
  int degree = dimension_data[0];
  RealVector x( Teuchos::View, pts_trans[dim], num_pts );
  RealVector basis_1d;
  bases1D_[dim]->derivative_set( x, degree, basis_1d );
  for ( int j = 0; j < num_pts; j++ )
    result[j] *= basis_1d[j];
};

void TensorProductBasis::gradient_set(const RealMatrix &pts,
				      const std::vector<PolyIndex_ptr>& indices,
				      RealMatrix &result ) const{
  int num_pts = pts.numCols();
  int num_indices = indices.size();

  IntVector max_degree_1d;
  compute_max_degree_for_each_dimension( indices, max_degree_1d );

  std::vector< std::vector<RealVector> > basis_values_1d;
  compute_1d_basis_values_for_each_dimension( pts, max_degree_1d, 
					      basis_values_1d );
  
  // derivatives are returned so that the first num_pts rows correspond to 
  // the 1st derivative for each point, the next set of rows are the 2nd 
  // derivative for each point and so on
  result.shapeUninitialized( numDims_*num_pts, num_indices );
  int start_row = 0;
  for ( int dim=0; dim< numDims_; dim++ ){
    RealMatrix dim_result( Teuchos::View, result, num_pts, 
			   num_indices, start_row, 0 );
    derivative_set( pts, indices, dim, max_degree_1d, basis_values_1d, 
		    dim_result );
    start_row += num_pts;
  }
}

void TensorProductBasis::
derivative_set( const RealMatrix &pts,
		const std::vector<PolyIndex_ptr>& indices,
		int derivative_dim,
		RealMatrix &result ) const{

  // This is slow but necessary if a derivative only in a single direction
  // is required
  IntVector max_degree_1d;
  compute_max_degree_for_each_dimension( indices, max_degree_1d );

  std::vector< std::vector<RealVector> > basis_values_1d;
  compute_1d_basis_values_for_each_dimension( pts, max_degree_1d, 
					      basis_values_1d );

  derivative_set( pts, indices, derivative_dim, max_degree_1d, 
		  basis_values_1d, result );
}

void TensorProductBasis::
derivative_set( const RealMatrix &pts,
		const std::vector<PolyIndex_ptr>& indices,
		int derivative_dim,
		const IntVector &max_degree_1d,
		const std::vector<std::vector<RealVector> > &basis_values_1d,
		RealMatrix &result ) const{

  int num_pts = pts.numCols();
  int num_indices = indices.size();

  std::vector<RealVector> basis_derivatives_1d;
  compute_1d_basis_derivatives_for_single_dimension(pts, derivative_dim,
						   max_degree_1d[derivative_dim],
						   basis_derivatives_1d );

  // allow result to be a view by not allocating new memory if result 
  // already has the correct shape
  if ( ( result.numRows() != num_pts ) || ( result.numCols() != num_indices ) )
    result.shapeUninitialized( num_pts, num_indices );

  result = 1.;
  for ( int i = 0; i < num_indices; i++ ){
    bool derivative_dim_active = false;
    for ( int d = 0; d < indices[i]->effective_dimension(); d++ ){
      int dim = indices[i]->dimension( d );
      int degree = indices[i]->level( d );
      if ( dim != derivative_dim ){
	for ( int j = 0; j < num_pts; j++ )
	  result(j,i) *= basis_values_1d[dim][j][degree];
      }else{
	derivative_dim_active = true;
	for ( int j = 0; j < num_pts; j++ )
	  result(j,i) *= basis_derivatives_1d[j][degree];
      }
    }
    // if the poly is a consant in the direction of the derivative then
    // the derivative will be zero
    if ( !derivative_dim_active ){
      for ( int j = 0; j < num_pts; j++ )
	result(j,i) = 0.;
    }
  }
}

void TensorProductBasis::
derivative_set_deprecated( const RealMatrix &pts,
			   std::vector<PolyIndex_ptr>& indices,
			   int derivative_dim,
			   RealMatrix &result ) const
{
  int num_pts = pts.numCols();
  RealMatrix pts_trans( pts, Teuchos::TRANS );

  int num_indices = indices.size();
  result.shapeUninitialized( num_pts, num_indices );
  result = 1.;

  IntVector max_degree_1d( numDims_ );
  std::vector< std::set<int> > degrees_1d( numDims_ );
  std::vector<bool> derivative_dim_active( num_indices, false );
  for ( int i = 0; i < num_indices; i++ )
    {
      for ( int d = 0; d < indices[i]->effective_dimension(); d++ )
	{
	  int dim = indices[i]->dimension( d );
	  int degree = indices[i]->level( d );
	  max_degree_1d[dim] = std::max( max_degree_1d[dim], degree );
	  degrees_1d[dim].insert( degree );
	  if ( dim == derivative_dim )
	    derivative_dim_active[i] = true;
	}
    }
    
  // compute 1d bases one pt at a time. Unlike value_set this should save
  // memory by not computing each basis for all indices up front
  std::set<int>::iterator it;
  std::vector< RealVector > basis_values_1d( numDims_ );
  for ( int k = 0; k < num_pts; k++ )
    {
      for ( int dim = 0; dim < numDims_; dim++ )
	{
	  Real x =  pts(dim,k);
	  if ( dim != derivative_dim )
	    {
	      basis_values_1d[dim].sizeUninitialized( max_degree_1d[dim]+1 );
	      for ( it = degrees_1d[dim].begin(); it != degrees_1d[dim].end(); 
		    ++it )
		{
		  if ( (*it) > 0  )
		    basis_values_1d[dim][*it]=bases1D_[dim]->value( x, (*it) );
		}
	    }
	  else
	    {
	      basis_values_1d[dim].sizeUninitialized( max_degree_1d[dim]+1 );
	      for ( it = degrees_1d[dim].begin(); it != degrees_1d[dim].end(); 
		    ++it )
		{
		  basis_values_1d[dim][*it]=bases1D_[dim]->derivative( x, (*it) );
		}
	    }
	}

      for ( int i = 0; i < num_indices; i++ )
	{
	  if( derivative_dim_active[i] )
	    {
	      for ( int d = 0; d < indices[i]->effective_dimension(); d++ )
		{
		  int dim = indices[i]->dimension( d );
		  int degree = indices[i]->level( d );
		  result(k,i) *= basis_values_1d[dim][degree];
		}
	    }
	  else result(k,i) = 0.;
	}
    }
};

OrthogPoly_ptr TensorProductBasis::get_basis_1d( int dim )
{
  return bases1D_[dim];
}


#ifdef PRECOMPUTE_BASES_OLD
void TensorProductBasis::value_set( const RealMatrix &pts,
				    std::vector<PolyIndex_ptr>& indices,
				    RealMatrix &result ) const
{
  int num_pts = pts.numCols();
  RealMatrix pts_trans( pts, Teuchos::TRANS );

  int num_indices = indices.size();
  result.shapeUninitialized( num_pts, num_indices );
  result = 1.;

  IntVector max_degree_1d( numDims_ );
  std::vector< std::set<int> > degrees_1d( numDims_ );
  for ( int i = 0; i < num_indices; i++ )
    {
      for ( int d = 0; d < indices[i]->effective_dimension(); d++ )
	{
	  int dim = indices[i]->dimension( d );
	  int degree = indices[i]->level( d );
	  max_degree_1d[dim] = std::max( max_degree_1d[dim], degree );
	  degrees_1d[dim].insert( degree );
	}
    }
    
  // compute 1d bases
  std::set<int>::iterator it;
  std::vector< std::vector<RealVector> > basis_values_1d( numDims_ );
  for ( int dim = 0; dim < numDims_; dim++ )
    {
      basis_values_1d[dim].resize( max_degree_1d[dim]+1 );
      RealVector x( Teuchos::View, pts_trans[dim], num_pts );
      for ( it = degrees_1d[dim].begin(); it != degrees_1d[dim].end(); ++it )
	{
	  if ( (*it) > 0  )
	    {
	      bases1D_[dim]->value_set( x, (*it), 
					basis_values_1d[dim][*it] );
	    }
	}
    }

  for ( int i = 0; i < num_indices; i++ )
    {
      for ( int d = 0; d < indices[i]->effective_dimension(); d++ )
	{
	  int dim = indices[i]->dimension( d );
	  int degree = indices[i]->level( d );
	  for ( int j = 0; j < num_pts; j++ )
	    result(j,i) *= basis_values_1d[dim][degree][j];
	}
    }
};
#endif //PRECOMPUTE_BASES_OLD

#ifdef PRECOMPUTE_BASES_SINGLE

void TensorProductBasis::value_set( const RealMatrix &pts,
				    std::vector<PolyIndex_ptr>& indices,
				    RealMatrix &result ) const
{
  
  int num_pts = pts.numCols();
  int num_indices = indices.size();
  result.shapeUninitialized( num_pts, num_indices );
  result = 1.;

  IntVector max_degree_1d( numDims_ );
  std::vector< std::set<int> > degrees_1d( numDims_ );
  for ( int i = 0; i < num_indices; i++ )
    {
      for ( int d = 0; d < indices[i]->effective_dimension(); d++ )
	{
	  int dim = indices[i]->dimension( d );
	  int degree = indices[i]->level( d );
	  max_degree_1d[dim] = std::max( max_degree_1d[dim], degree );
	  degrees_1d[dim].insert( degree );
	}
    }

  std::set<int>::iterator it;
  std::vector< RealVector > basis_values_1d( numDims_ );

  for ( int j = 0; j < num_pts; j++ )
    {
      // compute 1d bases
      for ( int dim = 0; dim < numDims_; dim++ )
	{
	  basis_values_1d[dim].resize( max_degree_1d[dim]+1 );
	  for ( it = degrees_1d[dim].begin(); it != degrees_1d[dim].end(); ++it )
	    {
	      if ( (*it) > 0  )
		{
		  basis_values_1d[dim][*it] = bases1D_[dim]->value( pts(dim,j), (*it) );
		}
	    }
	}
      
      for ( int i = 0; i < num_indices; i++ )
	{
	  for ( int d = 0; d < indices[i]->effective_dimension(); d++ )
	    {
	      int dim = indices[i]->dimension( d );
	      int degree = indices[i]->level( d );
	      result(j,i) *= basis_values_1d[dim][degree];
	    }
	}
    }
};
#endif //PRECOMPUTE_BASES_SINGLE

#ifdef NAIVE_COMPUTE

void TensorProductBasis::value_set( const RealMatrix &pts,
				    std::vector<PolyIndex_ptr>& indices,
				    RealMatrix &result ) const
{
  int num_pts = pts.numCols();
  int num_indices = indices.size();
  result.shapeUninitialized( num_pts, num_indices );
  result = 1.;

  for ( int i = 0; i < num_indices; i++ )
    {
      for ( int d = 0; d < indices[i]->effective_dimension(); d++ )
	{
	  int dim = indices[i]->dimension( d );
	  int degree = indices[i]->level( d );
	  for ( int j = 0; j < num_pts; j++ )
	    {
	      Real value = bases1D_[dim]->value( pts(dim,j), degree );
	      result(j,i) *= value;
	    }
	}
    }
};
#endif // NAIVE_COMPUTE

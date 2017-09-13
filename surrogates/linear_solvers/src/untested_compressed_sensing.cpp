#include "untested_compressed_sensing.hpp"
#include "MathTools.hpp"
#include "Printing.hpp"


void get_forward_neighbours( PolyIndex_ptr index, 
		      int num_dims,
		      bool expansion_type,
		      std::set<int> &non_zero_array_indices,
		      IndexHashSet<PolynomialIndex>::hashSet &indices,
		      IndexHashSet<PolynomialIndex>::hashSet &neighbours )
{
  IndexHashSet<PolynomialIndex>::const_iterator iter;
  for ( int dim = 0; dim < num_dims; dim++ )
    {
      PolyIndex_ptr forward_neighbour( new PolynomialIndex() );
      index->get_forward_neighbour( dim, forward_neighbour );
      iter = neighbours.find( forward_neighbour );
      if ( iter == neighbours.end() )
	{
	  bool admissible = true;
	  for ( int d = 0; d < forward_neighbour->effective_dimension(); d++ )
	    {
	      int back_dim = forward_neighbour->dimension( d );
	      PolyIndex_ptr back_neighbour( new PolynomialIndex() );
	      forward_neighbour->get_backward_neighbour( back_dim,
							 back_neighbour );
	      iter = indices.find( back_neighbour );
	      if ( iter == indices.end() )
		{
		  admissible = false;
		  break;
		}
	      else if ( expansion_type == 1 )
		{
		  std::set<int>::const_iterator it = 
		    non_zero_array_indices.find( (*iter)->get_array_index() );
		  if ( it == non_zero_array_indices.end() )
		    admissible = false;
		}
	    }
	  if ( admissible )
	    {
	      neighbours.insert( forward_neighbour );
	    }
	}
    }
}

void get_ancestors( PolyIndex_ptr index, int num_dims, 
		    std::set<int> &non_zero_array_indices,
		    IndexHashSet<PolynomialIndex>::hashSet &indices,
		    IndexHashSet<PolynomialIndex>::hashSet &ancestors )
{
  IndexHashSet<PolynomialIndex>::const_iterator iter;
  for ( int d = 0; d < index->effective_dimension(); d++ )
    {
      int back_dim = index->dimension( d );
      PolyIndex_ptr ancestor( new PolynomialIndex() );
      index->get_backward_neighbour( back_dim, ancestor );
      iter = indices.find( ancestor );
      if ( iter == indices.end() )
	{
	  throw( std::runtime_error("Ancestor must exist") );
	}
      std::set<int>::const_iterator it = 
	non_zero_array_indices.find( (*iter)->get_array_index() );
      if ( it == non_zero_array_indices.end() )
	{
	  /*std::cout << "The ancestor: " << ancestor->to_string( num_dims )
		    << " of the index " << index->to_string( num_dims ) 
		    << " was missing\n";*/
	  ancestors.insert( *iter );
	  get_ancestors( ancestor, num_dims, non_zero_array_indices, 
			 indices, ancestors );
	}
    }
}

void get_children( PolyIndex_ptr index, int num_dims, int depth, 
		   int expansion_type,
		   std::set<int> &non_zero_array_indices,
		   IndexHashSet<PolynomialIndex>::hashSet &indices,
		   IndexHashSet<PolynomialIndex>::hashSet &children )
{
  
  //std::cout << "-------------------\n";
  //std::cout << index->to_string( num_dims ) << std::endl;

  IndexHashSet<PolynomialIndex>::const_iterator iter, it;
  IndexHashSet<PolynomialIndex>::hashSet neighbours, expanded_indices;
  std::vector< PolyIndex_ptr > active_indices( 1 );
  active_indices[0] = index;
  //shallow copy is necessary for checking admissibility criterion
  expanded_indices = indices; 

  for ( int i = 0; i < depth; i++ )
    {
      int expansion_type_i = 0;
      if ( i == 0 ) expansion_type_i = expansion_type;
      for ( int j = 0; j < (int)active_indices.size(); j++ )
	{
	  get_forward_neighbours( active_indices[j], num_dims, 
				  expansion_type_i,
				  non_zero_array_indices,
				  expanded_indices, 
				  neighbours );
	}
      active_indices.resize( neighbours.size() );
      int k = 0;
      for ( iter = neighbours.begin(); iter != neighbours.end(); ++iter )
	{
	  active_indices[k] = *iter;
	  // below is necessary for checking admissibility criterion
	  it = indices.find( *iter );
	  if ( it == indices.end() )
	    {
	      children.insert( *iter );
	      expanded_indices.insert( *iter );
	    }
	  k++;
	}
      neighbours.clear();
    }
}

void expand_basis( PolyIndex_ptr index,
		   std::vector<PolyIndex_ptr> &indices, 
		   IndexHashSet<PolynomialIndex>::hashSet &indices_set,
		   int num_dims, int depth,
		   RealMatrix &build_points,
		   RealMatrix &A, RealVector &b,
		   TensorProductBasis_ptr basis,
		   bool normalise,
		   RealVector &column_norms,
		   RealMatrix &Atb,
		   std::set<int> &non_zero_array_indices, 
		   int expansion_type,
		   IndexHashSet<PolynomialIndex>::hashSet &ancestors )
{
  //std::cout << "expand" << index->to_string( num_dims ) << std::endl;
  int M = A.numRows();
  Teuchos::BLAS<int, Real> blas;
  // add new basis elements
  IndexHashSet<PolynomialIndex>::const_iterator iter, it;
  int num_old_indices = indices.size();
	  
  //get_ancestors( index, num_dims, non_zero_array_indices,
  //indices_set, ancestors );
  IndexHashSet<PolynomialIndex>::hashSet children;
  get_children( index, num_dims, depth, 
		expansion_type,
		non_zero_array_indices,
		indices_set, 
		children );

  int num_new_indices = ancestors.size() + children.size();
  std::vector<PolyIndex_ptr> new_indices( num_new_indices );
  int i = 0;
  for ( iter = children.begin(); iter != children.end(); ++iter )
    {
      bool add_candidate = false;
      if ( ancestors.empty() )
	add_candidate = true;
      else
	{
	  IndexHashSet<PolynomialIndex>::const_iterator it = 
	     ancestors.find( *iter );
	  if ( it == ancestors.end() )
	    add_candidate = true;
	}
      if ( add_candidate )
	{ 
	  //std::cout << "c" << (*iter)->to_string( num_dims )
	  //	    << "," << (*iter)->get_array_index()
	  //	    << "," << num_old_indices + i << std::endl;
	  indices_set.insert( *iter );
	  (*iter)->set_array_index( num_old_indices + i );
	  new_indices[i] = *iter;
	  indices.push_back( *iter );
	  i++;
	}
    }
  new_indices.resize( i );
  num_new_indices = new_indices.size();

  RealMatrix A_ext;
  basis->value_set( build_points, new_indices, A_ext );
  if ( normalise )
    {
      RealVector new_column_norms;
      normalise_columns( A_ext, new_column_norms );
      column_norms.resize( num_old_indices + num_new_indices );
      for ( i = 0; i < num_new_indices; i++ )
	column_norms[num_old_indices+i] = new_column_norms[i];
    }
  // append A_ext to A
  column_append( A_ext, A );
  // update Atb
  Atb.reshape( num_old_indices + num_new_indices, 1 );
  for ( i = 0; i < num_new_indices; i++ )
    Atb(num_old_indices+i,0) = blas.DOT( M, A_ext[i], 1, b.values(), 1 );

}


void tomp_update_memory( RealMatrix &Q, RealMatrix &R, 
			 RealMatrix &AtA_sparse_memory, 
			 RealMatrix &A_sparse_memory,
			 RealMatrix &solutions, RealMatrix &solution_metrics,
			 int new_size,
			 int memory_chunk_size )
{
  if ( Q.numCols() <= new_size )
    {
      //disp( "updating memory" );
      Q.reshape( Q.numRows(), Q.numCols() + memory_chunk_size );
      R.reshape( R.numRows() + memory_chunk_size, 
		 R.numCols() + memory_chunk_size );
      AtA_sparse_memory.reshape( AtA_sparse_memory.numRows(), 
				 AtA_sparse_memory.numCols() + 
				 memory_chunk_size );
      A_sparse_memory.reshape( A_sparse_memory.numRows(), 
			       A_sparse_memory.numCols() + memory_chunk_size );
      solutions.reshape( solutions.numRows(), 
			 solutions.numCols() + memory_chunk_size );
      solution_metrics.reshape( solution_metrics.numRows(),
				solution_metrics.numCols() + memory_chunk_size );
    }  
}

bool index_already_added( std::set<int> &non_zero_array_indices, 
			  int active_index, int verbosity )
{
  // todo define active_index_set as std::set and use find function
  std::set<int>::iterator iter = non_zero_array_indices.find( active_index );
  if ( iter != non_zero_array_indices.end() )
    {
      if ( verbosity > 1 ){
	std::cout << "Exiting: New active index " <<  active_index
		  << " has already been added. "
		  << "This has likely occured because all correlations "
		  << "are roughly the same size. "
		  << "This means problem has been solved to roughly "
		  << "machine precision. Check this before continuing. ";
      }
      return true;
    }
  return false;
}

#ifdef QR_UPDATED_TO_USE_NEW_QR_FACTORIZTION_UPDATE
// Cannot use because qr_factorization update changed. Only small change so easy to fix
void compute_residual_using_subtree(
PolynomialChaosExpansion *pce,			    
RealMatrix &A,
RealVector &b,
RealMatrix &A_sparse_memory, 
RealMatrix &AtA_sparse_memory,
RealMatrix &Atb,
RealMatrix &Atb_sparse_memory,
RealMatrix &Q,
RealMatrix &R,
RealMatrix &solutions,
RealMatrix &solution_metrics,
RealVector &residual,
RealVector &x_sparse,
int active_index,
std::vector<PolyIndex_ptr> &indices, 
IndexHashSet<PolynomialIndex>::hashSet &indices_set,
std::set<int> &non_zero_array_indices,
bool update,
int memory_chunk_size,
IntVector &active_indices,
bool add_ancestors )


{
  int M = Q.numRows();
  Teuchos::BLAS<int, Real> blas;

  IndexHashSet<PolynomialIndex>::hashSet ancestors;
  if ( add_ancestors )
    get_ancestors( indices[active_index], pce->dimension(), 
		   non_zero_array_indices,
		   indices_set, ancestors );

  ancestors.insert( indices[active_index] );
  tomp_update_memory( Q, R, AtA_sparse_memory, A_sparse_memory,
		      solutions, solution_metrics, 
		      non_zero_array_indices.size() + ancestors.size(),
		      memory_chunk_size );


  IndexHashSet<PolynomialIndex>::const_iterator iter;
  int num_active_indices = non_zero_array_indices.size();
  for ( iter = ancestors.begin(); iter != ancestors.end(); ++iter )
    {
      int index = (*iter)->get_array_index();
      if ( update )
	{
	  non_zero_array_indices.insert( index );
	  active_indices[num_active_indices] = index;
	}
      
      RealMatrix A_col( Teuchos::View, A, M, 1, 0, index );
      int colinear = qr_factorization_update_insert_column( Q, R, A_col, 
							    num_active_indices );
      if ( colinear ) 
	{
	  //return std::numeric_limits<double>::max();
	  //throw( std::runtime_error("Colinear") );
	  std::cout << "warning colinear" << std::endl;
	  break;
	}

      RealMatrix Atb_sparse( Teuchos::View, Atb_sparse_memory, 
			     num_active_indices + 1, 1, 0, 0 );
	  
      Atb_sparse(num_active_indices,0) = Atb(index,0);

      //Solve R'z = A'b via back substitution
      RealMatrix R_new( Teuchos::View, R, num_active_indices+1, 
			num_active_indices+1, 0, 0 );

      RealMatrix A_sparse( Teuchos::View, A_sparse_memory, 
			   M, num_active_indices+1, 0, 0 );
      for ( int m = 0; m < M; m++ )
	A_sparse(m,num_active_indices) = A(m,index);

      num_active_indices++;
    }
  RealMatrix A_sparse( Teuchos::View, A_sparse_memory, 
		       M, num_active_indices, 0, 0 );
  RealMatrix Atb_sparse( Teuchos::View, Atb_sparse_memory, 
			 num_active_indices, 1, 0, 0 );

  RealMatrix z;
  RealMatrix R_new( Teuchos::View, R, num_active_indices, 
		    num_active_indices, 0, 0 );
  substitution_solve( R_new, Atb_sparse, z, Teuchos::TRANS );
  //Solve Rx = z via back substitution to obtain signal
  substitution_solve( R_new, z, x_sparse );
  residual.assign( b );
  blas.GEMV( Teuchos::NO_TRANS, M, num_active_indices, -1., 
	     A_sparse.values(), A_sparse.stride(), 
	     x_sparse.values(), 1, 1., residual.values(), 1 );
}

void tree_orthogonal_matching_pursuit( RealMatrix &A,
				       RealVector &b,
				       PolynomialChaosExpansion *pce,
				       RealMatrix &solutions,
				       RealMatrix &solution_metrics,
				       Real epsilon,
				       Real ratio,
				       int max_num_non_zero_entries,
				       int verbosity,
				       bool add_ancestors )
{
  Teuchos::BLAS<int, Real> blas;
  int M( A.numRows() ), N( A.numCols() );

  // Determine the maximum number of iterations
  int max_num_indices( std::min( M, max_num_non_zero_entries ) );
  max_num_indices = std::min( N, max_num_indices );

  // Vector to store non-zero solution entries
  RealVector x_sparse;	

  int memory_chunk_size = std::min( N, M );
  // if I use min( tmp, M ) where tmp < M then seq faults will occur.
  // Not sure why has something to do with resizing memory
  int initial_N = std::min( memory_chunk_size, N );

  // Initialise entries of all solutions to zero
  solutions.shape( N, initial_N );

  // Allocate memory to store solution metrics
  solution_metrics.shapeUninitialized( 2, initial_N );

  // Matrix to store Q and R factorization
  RealMatrix Q( M, initial_N ), R( initial_N, initial_N );

  // Compute residual
  RealVector residual( b );

  // Compute correlation of columns with residual
  RealMatrix Atb( N, 1, false );
  Atb.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 
		1.0, A, residual, 0.0 );
  RealMatrix correlation( Teuchos::Copy, Atb, N, 1 );
  
  // Compute norm of residual
  Real b_norm_sq( b.dot( b ) ), residual_norm( std::sqrt( b_norm_sq ) );

  // Matrix to store the full rank matrix associated with x_sparse
  RealMatrix A_sparse_memory( M, initial_N, false );
  RealMatrix AtA_sparse_memory( N, initial_N, false );
  RealMatrix Atb_sparse_memory( N, 1, false );

  if ( verbosity > 1 )
    {
      std::cout << "Tree Orthogonal Matching Pursuit\n";
      std::printf( "Iter\tNNZ\tResidual\tCorrelation\tl1 norm of x\n" );
    }

  std::vector<PolyIndex_ptr> indices;
  pce->get_basis_indices( indices );
  IndexHashSet<PolynomialIndex>::hashSet indices_set;
  for ( int i = 0; i < (int)indices.size(); i++ )
    {
      indices_set.insert( indices[i] );
      indices[i]->set_array_index( i );
    }
  std::set<int> non_zero_array_indices;
  IntVector active_indices( N, false );
  active_indices = -1;

  int iter = 0;
  bool done = false;
  while ( !done )
    {
      
      Real r_norm_min = std::numeric_limits<double>::max();
      int active_index = -1;
      // sort correlation
      IntVector sort_indices;
      RealVector corr_vector( Teuchos::View, correlation.values(), N );
      magnitude_argsort( corr_vector, sort_indices );
      Real max_correlation = std::abs( correlation(sort_indices[0],0) );
      for ( int i = 0; i < N - (int)non_zero_array_indices.size(); i++ )
	{
	  int index = sort_indices[i];
	  if ( std::abs( correlation(index,0) ) < ratio*max_correlation ) break;
	  compute_residual_using_subtree( pce, A, b,
					  A_sparse_memory, AtA_sparse_memory,
					  Atb, Atb_sparse_memory,
					  Q, R, solutions, solution_metrics,
					  residual, x_sparse,
					  index, indices, 
					  indices_set, non_zero_array_indices,
					  false, memory_chunk_size,
					  active_indices,
					  add_ancestors );
	  Real r_norm_i = residual.normFrobenius();
	  if ( r_norm_i < r_norm_min )
	    {
	      r_norm_min = r_norm_i;
	      active_index = index;
	    }
	}
      if ( index_already_added( non_zero_array_indices, active_index,
				verbosity ) )
	break;

      Real chosen_correlation = correlation(active_index,0);

      // Update solution with the best set of indices
      compute_residual_using_subtree( pce, A, b,
				      A_sparse_memory, AtA_sparse_memory,
				      Atb, Atb_sparse_memory,
				      Q, R, solutions, solution_metrics,
				      residual, x_sparse,
				      active_index, indices, 
				      indices_set, non_zero_array_indices,
				      true, memory_chunk_size,
				      active_indices,
				      add_ancestors );
      

      // Compute new correlations
      // correlation = A' * residual: O(MN) N >= M

      blas.GEMV( Teuchos::TRANS, M, N, 1., A.values(), A.stride(), 
		 residual.values(), 1, 0., correlation.values(), 1 );

      // Update solution history
      residual_norm = residual.normFrobenius();
      int num_active_indices = non_zero_array_indices.size();
      
      for ( int n = 0; n < num_active_indices; n++ )
	solutions(active_indices[n],iter) = x_sparse[n];

      solution_metrics(0,iter) = residual_norm;
      solution_metrics(1,iter) = num_active_indices;     

      if ( verbosity > 1 )
	std::printf( "%d\t%d\t%1.5e\t%1.5e\t%1.5e\n", iter, num_active_indices, 
		     residual_norm, chosen_correlation,
		     x_sparse.normOne() );
 
      if ( residual_norm <= epsilon )
	{
	  if ( verbosity > 1 )
	    std::cout << "Exiting: residual norm lower than tolerance\n";
	  done = true;
	}

      if ( num_active_indices >= max_num_indices )
	{
	  if ( verbosity > 1 )
	    std::cout << "Exiting: maximum number of non-zero terms obtained\n";
	  done = true;
	}
      iter++;
    }
  // remove unused memory
  solutions.reshape( N, iter );
  solution_metrics.reshape( 2, iter );
}
#endif

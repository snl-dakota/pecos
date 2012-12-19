#include "FaultTolerance.hpp"

namespace Pecos {

void fail_booleans(SizetShortMap::const_iterator& fit, size_t j,
		   bool& add_val, bool& add_grad,
		   const SizetShortMap& failed_response_data )
{
  if (fit != failed_response_data.end() && fit->first == j) {
    short fail_asv = fit->second;
    if (fail_asv & 1) add_val  = false;
    if (fail_asv & 2) add_grad = false;
    ++fit;
  }
}

void remove_faulty_data( RealMatrix &A, RealMatrix &B, 
		    bool expansion_coeff_grad_flag,
		    IntVector &index_mapping,
		    FaultInfo fault_info,
		    const SizetShortMap& failed_resp_data )
{
  size_t  num_rows_A = A.numRows(), num_cols_A  = A.numCols(), 
    num_rhs( B.numCols() );
    
  IntVector sorted_indices, I;
  if ( index_mapping.length() == 0 )
    {
      if ( num_rows_A % fault_info.num_surr_data_pts != 0  ){
	std::string msg = "remove_faulty_data() A matrix and index_mapping are ";
	msg += "inconsistent";
	throw( std::runtime_error( msg ) );
      }

      range(index_mapping,0,(int)fault_info.num_surr_data_pts);  
    }
  // compute the number of rows in a that are function values. The remainder
  // are gradient values
  size_t num_data_pts_fn_in_A = num_rows_A / fault_info.num_surr_data_pts;
  if ( num_rows_A % fault_info.num_surr_data_pts != 0 ) num_data_pts_fn_in_A++;
  num_data_pts_fn_in_A = num_rows_A / num_data_pts_fn_in_A;

  size_t m;
  // When fault tolerance is applied with cross validation some of the 
  // entries in the index_mapping will be -1. These values need to be removed.
  // -1's will always be at the end of the mapping
  for ( m=index_mapping.length()-1; m>=0; --m){
    if (index_mapping[m]>=0) break;
  }
  IntVector index_map( Teuchos::Copy, index_mapping.values(), (int)m+1 );

  // sorted indices are needed to use fail_booleans
  // I is needed to access the correct row in A.
  argsort( index_map, sorted_indices, I );
  
  SizetShortMap::const_iterator fit;

  if ( !expansion_coeff_grad_flag )
    {
      // A  = [A_fn; A_grad]
      size_t i, j, k, l, a_fn_cntr, a_grad_cntr;
      // Initialise memory of A_fn and A_grad to largest possible sizes.
      // cross validation will result in not all entries being filled.
      // thus we must resize after the entries are added.
      RealMatrix A_fn( fault_info.num_data_pts_fn, num_cols_A, false ), 
	A_grad( fault_info.num_data_pts_grad*fault_info.num_vars, num_cols_A, 
		false );
      for (i=0; i<num_cols_A; ++i) {
	a_fn_cntr = 0; a_grad_cntr = 0;
	for ( j=0, l=0, fit=failed_resp_data.begin(); 
	      j<fault_info.num_surr_data_pts;  j++){
	  int index = I[l];
	  bool add_val = true, add_grad = fault_info.use_derivatives;
	  fail_booleans(fit, sorted_indices[l], add_val, add_grad, 
			failed_resp_data );
	  if ( sorted_indices[l] != j ) add_val = add_grad = false;
	  if ( add_val ){
	    A_fn(a_fn_cntr,i) = A(index,i);
	    a_fn_cntr++;
	  }
	  if ( add_grad ){
	    for ( k=0; k<fault_info.num_vars; ++k, ++a_grad_cntr ){
	      A_grad(a_grad_cntr,i) = 
		A(num_data_pts_fn_in_A+index*fault_info.num_vars+k,i);
	    }
	  }
	  if ( l < index_map.length() && sorted_indices[l] == j ) l++;
	  if ( l >= index_map.length() ) break;
	}
      }
      // Resize matrix to include only non-faulty data.
      // To efficiently copy a submatrix first take view of submatrix
      // then use assign ( do not use = operator as this will only provide 
      // a view
      RealMatrix A_fn_tmp( Teuchos::View, A_fn, a_fn_cntr, num_cols_A );
      A.reshape( a_fn_cntr, num_cols_A );
      A.assign( A_fn_tmp );
      if ( fault_info.use_derivatives )
	{
	  RealMatrix A_grad_tmp( Teuchos::View, A_grad, a_grad_cntr, 
				 num_cols_A );
	  row_append( A_grad_tmp, A );
	}

      size_t b_fn_cntr = 0, b_grad_cntr = 0;
      RealMatrix B_fn( fault_info.num_data_pts_fn, num_rhs, false ), 
	B_grad( fault_info.num_data_pts_grad*fault_info.num_vars, num_rhs, 
		false );
      for (j=0, l=0, fit=failed_resp_data.begin(); 
	   j<fault_info.num_surr_data_pts; ++j) {
	int index = I[l];
	bool add_val = true, add_grad = fault_info.use_derivatives;
	fail_booleans(fit, sorted_indices[l], add_val, add_grad,
		      failed_resp_data);
	if ( sorted_indices[l] != j ) add_val = add_grad = false;
	if ( add_val )
	  { B_fn(b_fn_cntr,0) = B(index,0); ++b_fn_cntr; }
	if ( add_grad ) {
	  for (k=0; k<fault_info.num_vars; ++k, ++b_grad_cntr){
	    B_grad(b_grad_cntr,0) = 
	      B(num_data_pts_fn_in_A+index*fault_info.num_vars+k,0);
	  }
	}
	if ( l < index_map.length() && sorted_indices[l] == j ) l++;
	if ( l >= index_map.length() ) break;
      }
      // Resize matrix to include only non-faulty data
      RealMatrix B_fn_tmp( Teuchos::View, B_fn, b_fn_cntr, num_rhs );
      B.reshape( b_fn_cntr, num_rhs );
      B.assign( B_fn_tmp );
      if ( fault_info.use_derivatives )
	{
	  RealMatrix B_grad_tmp( Teuchos::View, B_grad, b_grad_cntr, num_rhs );
	  row_append( B_grad_tmp, B );
	}
    }
  /*else
    {
    //when computing coef of grad expansion
    }
  */
};

} // namespace Pecos

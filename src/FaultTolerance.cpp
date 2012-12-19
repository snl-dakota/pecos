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
  //PCout << A.numRows()  << "," << A.numCols() << " A\n";

  argsort( index_mapping, sorted_indices, I );
  //index_mapping.print(std::cout);
  //I.print(std::cout);
  //sorted_indices.print(std::cout);
  
  SizetShortMap::const_iterator fit;

  if ( !expansion_coeff_grad_flag )
    {
      size_t i, j, k, l, a_fn_cntr, a_grad_cntr;
      RealMatrix A_fn( fault_info.num_data_pts_fn, num_cols_A, false ), 
	A_grad( fault_info.num_data_pts_grad*fault_info.num_vars, num_cols_A, false );
      for (i=0; i<num_cols_A; ++i) {
	a_fn_cntr = 0; a_grad_cntr = 0;
	for ( j=0, l=0, fit=failed_resp_data.begin(); 
	      j<fault_info.num_surr_data_pts; j++){
	  int index = index_mapping[I[l]];
	  //PCout << index << " i\n";
	  bool add_val = true, add_grad = fault_info.use_derivatives;
	  //PCout << sorted_indices[l] << "," << fit->first<<" s\n";
	  fail_booleans(fit, sorted_indices[l], add_val, add_grad, 
			failed_resp_data );
	  if ( add_val ){
	    A_fn(a_fn_cntr,i) = A(index,i);
	    a_fn_cntr++;
	  }
	  if ( add_grad ){
	    for ( k=0; k<fault_info.num_vars; ++k, ++a_grad_cntr ){
	      //PCout << a_grad_cntr << "," << num_surr_data_pts+index*fault_info.num_vars+k << " \n";
	      A_grad(a_grad_cntr,i) = 
		A(fault_info.num_surr_data_pts+index*fault_info.num_vars+k,i);
	    }
	  }
	  if ( l < index_mapping.length() && sorted_indices[l] == j ) l++;
	}
      }

      A = A_fn;
      if ( fault_info.use_derivatives )
	row_append( A_grad, A );

      if ( A.numRows() != fault_info.total_eqns )
	throw( std::runtime_error(""));

      size_t b_fn_cntr = 0, b_grad_cntr = 0;
      RealMatrix B_fn( fault_info.num_data_pts_fn, num_rhs, false ), 
	B_grad( fault_info.num_data_pts_grad*fault_info.num_vars, num_rhs, false );
      for (i=0, l=0, fit=failed_resp_data.begin(); 
	   i<fault_info.num_surr_data_pts; ++i) {
	int index = index_mapping[I[l]];
	bool add_val = true, add_grad = fault_info.use_derivatives;
	fail_booleans(fit, sorted_indices[l], add_val, add_grad,
		      failed_resp_data);
	if ( add_val )
	  { B_fn(b_fn_cntr,0) = B(index,0); ++b_fn_cntr; }
	if ( add_grad ) {
	  for (k=0; k<fault_info.num_vars; ++k, ++b_grad_cntr){
	      B_grad(b_grad_cntr,0) = 
		B(fault_info.num_surr_data_pts+index*fault_info.num_vars+k,0);
	  }
	}
	if ( l < index_mapping.length() && sorted_indices[l] == i ) l++;
      }

      B = B_fn;
      if ( fault_info.use_derivatives )
	row_append( B_grad, B );
    }
  /*else
    {
    //when computing coef of grad expansion
    }
  */
};

} // namespace Pecos

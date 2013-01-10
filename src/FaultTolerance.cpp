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

// \todo check if system becomes under-determined. perhaps return flag
// specifying this. Then if flag is true a different solver can still be
// used
void remove_faulty_data( RealMatrix &A, RealMatrix &B, 
			 IntVector &index_mapping,
			 FaultInfo fault_info,
			 const SizetShortMap& failed_resp_data )
{
  RealMatrix A_new, B_new;
  
  B.print(std::cout);
  index_mapping.print(std::cout);
  
  // tmp remove
  PCout << "Using " << fault_info.num_data_pts_fn << " equations to solve for ";
  PCout << A.numCols() << " coefficients\n";

  PCout << "num rows " << A.numRows() << std::endl;

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
  I.print(std::cout);
  sorted_indices.print(std::cout);
  
  SizetShortMap::const_iterator fit;

  size_t num_grad_rhs = fault_info.num_grad_rhs;
  size_t num_coef_rhs = B.numCols() - num_grad_rhs;
  if ( B.numCols() < num_grad_rhs ) { num_coef_rhs = 1; num_grad_rhs = 0; }
  if ( num_grad_rhs > A.numCols() ){
    num_coef_rhs = 1;
    num_grad_rhs = 0;
  }

  bool multiple_rhs = false;
  if ( ( num_coef_rhs > 0 ) && ( num_grad_rhs > 0 ) ) multiple_rhs = true;

  PCout <<  num_coef_rhs << "," << num_grad_rhs << "," << num_rhs << "," << B.numCols() << "\n";
  PCout <<  fault_info.num_data_pts_fn << "," << fault_info.num_data_pts_grad << "," << fault_info.num_vars << "\n";

  if ( num_rhs != num_grad_rhs + num_coef_rhs ){
    std::string msg = "remove_faulty_data. Fault_info is inconsistent with A";
    throw(std::runtime_error( msg ) );
  }
  RealMatrix B_grad;
  if ( ( multiple_rhs ) || ( num_coef_rhs == 1 ) )
    {
      PCout << "#############\n";
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
	  //PCout << index << "," << j << "," << l << "," << sorted_indices[l] << "," << fit->first << std::endl;
	  //fail_booleans(fit, sorted_indices[l], add_val, add_grad, 
	  //failed_resp_data );
	  fail_booleans(fit, j, add_val, add_grad, //new
			failed_resp_data );
	  if ( sorted_indices[l] != j ) {
	    add_val = add_grad = false;
	  }
	  if ( add_val ){
	    //PCout << A(index,i) << "," << B(index,0)<<std::endl;
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
      A_new.shapeUninitialized( a_fn_cntr, num_cols_A );
      A_new.assign( A_fn_tmp );
      if ( fault_info.use_derivatives )
	{
	  RealMatrix A_grad_tmp( Teuchos::View, A_grad, a_grad_cntr, 
				 num_cols_A );
	  row_append( A_grad_tmp, A_new );
	}

      PCout << "$$$$$$$$$$$$$\n";
      size_t b_fn_cntr = 0, b_grad_cntr = 0;
      RealMatrix B_fn( fault_info.num_data_pts_fn, num_rhs, true );
      B_grad.shapeUninitialized( fault_info.num_data_pts_grad * 
				 fault_info.num_vars, num_rhs );
      for (j=0, l=0, fit=failed_resp_data.begin(); 
	   j<fault_info.num_surr_data_pts; ++j) {
	int index = I[l];
	bool add_val = true, add_grad = fault_info.use_derivatives;
	fail_booleans(fit, j, add_val, add_grad, //new
		      failed_resp_data );
	//fail_booleans(fit, sorted_indices[l], add_val, add_grad,
	//	      failed_resp_data);
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
      B_new.shapeUninitialized( b_fn_cntr, num_rhs );
      PCout << "MATRICES:" << std::endl;
      B_new.assign( B_fn_tmp );
      A_new.print(std::cout);
      B_new.print(std::cout);
      if ( fault_info.use_derivatives )
	{
	  PCout << "@@@@@@@@@@@@@@@@@@@@@\n";
	  RealMatrix B_grad_tmp( Teuchos::View, B_grad, b_grad_cntr, 
				 num_rhs );
	  row_append( B_grad_tmp, B_new );
	}
      PCout << "Using " << A_new.numRows() << " equations to solve for ";
      PCout << A_new.numCols() << " coefficients\n";
    }

  if ( num_grad_rhs > 0 )
    {
      PCout << "!!!!!!!!!!!!!!\n";
      // used when computing coef of grad expansion

      if ( !multiple_rhs )
	{
	  PCout << "&&&&&&&&&&&&&&&&&\n";
	  size_t i, j, k, l, a_grad_cntr;
	  RealMatrix A_grad( fault_info.num_data_pts_grad, num_cols_A, false );
	  for (i=0; i<num_cols_A; ++i) {
	    a_grad_cntr = 0;
	    for ( j=0, l=0, fit=failed_resp_data.begin(); 
		  j<fault_info.num_surr_data_pts;  j++){
	      int index = I[l];
	      bool add_val = false, add_grad = true;
	      fail_booleans(fit, sorted_indices[l], add_val, add_grad, 
			    failed_resp_data );
	      if ( sorted_indices[l] != j ) add_val = add_grad = false;
	      if ( add_grad ){
		A_grad(a_grad_cntr,i) =  A(index,i);
		a_grad_cntr++;
	      }
	      if ( l < index_map.length() && sorted_indices[l] == j ) l++;
	      if ( l >= index_map.length() ) break;
	    }
	  }
	  RealMatrix A_grad_tmp( Teuchos::View, A_grad, a_grad_cntr, 
				 num_cols_A );

	  A_new.reshape( a_grad_cntr, num_cols_A );
	  A_new.assign( A_grad_tmp );
	  B_new.shapeUninitialized( fault_info.num_data_pts_grad, 
				    num_rhs );
	}
      //PCout << A_new.numRows() << "," << A_new.numCols() << std::endl;
      //PCout << B_new.numRows() << "," << B_new.numCols() << std::endl;
      size_t i, j, k, l, b_grad_cntr = 0;
      for (j=0, l=0, fit=failed_resp_data.begin(); 
	   j<fault_info.num_surr_data_pts; ++j) {
	int index = I[l];
	bool add_val = false, add_grad = true;
	fail_booleans(fit, sorted_indices[l], add_val, add_grad,
		      failed_resp_data);
	if ( sorted_indices[l] != j ) add_val = add_grad = false;
	if ( add_grad ) {
	  for (k=0; k<num_grad_rhs; ++k ){
	    //PCout << j << "," << b_grad_cntr <<  "," << k << "," << num_coef_rhs+k << "," << index << "\n";
	    B_new(b_grad_cntr,num_coef_rhs+k) = 
	      B(index,num_coef_rhs+k);
	  }
	  ++b_grad_cntr;
	}
	if ( l < index_map.length() && sorted_indices[l] == j ) l++;
	if ( l >= index_map.length() ) break;
      }
    }

  //RealMatrix Z( A_new );
  //Z -= A;
  //PCout << Z.normInf() << std::endl;
  //Z = B_new;
  //Z -= B;
  //PCout <<  Z.normInf() << std::endl;

  A = A_new;
  B = B_new;
  
};

} // namespace Pecos

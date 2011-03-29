/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/** \file pecos_hierarchical_driver.cpp
    \brief A test program for HierarchicalLinearBasis class. */


#include "../src/HierarchicalLinearBasis.hpp"

using namespace Pecos;
int main(int argc, char** argv)
{
  
  HierarchicalLinearBasis *a, *b, *c, *d;
  //Allocate Objects using default parameters and non default.
  a = new HierarchicalLinearBasis();
  b = new HierarchicalLinearBasis(-6.5,4.7);

  //Allocate with left >= right.  Should throw an exception.
  bool caughtException1 = false;
  bool caughtException2 = false;
  try {
    c = new HierarchicalLinearBasis(10.3,9.1);
  } catch (const std::exception& exception) {
    caughtException1 = true;
  }
  //Allocate with left==right.  Should throw an exception.
  try {
    d = new HierarchicalLinearBasis(10.3,10.3);
  } catch (const std::exception& exception) {
    caughtException2 = true;
  }

  if ( !( caughtException1 && caughtException2 ) ){
    std::cout << "Test fail: Constructor should have thrown"
	      << " an exception.";
    return EXIT_FAILURE;
  }

  delete a,b,c,d;

  a = new HierarchicalLinearBasis();
  b = new HierarchicalLinearBasis(-1.2,5.6);

  if(a->get_lower_interp_bound() != 0.0){
    std::cout << "Test fail: get_lower_interp_bound()";
    return EXIT_FAILURE;
  }
  if(a->get_lower_interp_bound() != 0.0){
    std::cout << "Test fail: get_lower_interp_bound()";
    return EXIT_FAILURE;
  }
  if(a->get_lower_interp_bound() != 0.0){
    std::cout << "Test fail: get_lower_interp_bound()";
    return EXIT_FAILURE;
  }
  if(a->get_lower_interp_bound() != 0.0){
    std::cout << "Test fail: get_lower_interp_bound()";
    return EXIT_FAILURE;
  }

  delete a,b;

  a = new HierarchicalLinearBasis(-1,1);
  IntArray basis_index(2);
  //Test the highest basis function.  Should be constant 1
  // on the interval and zero off of it.
  basis_index[0] = 1;
  basis_index[1] = 0;
  if( !((a->get_type1_value(0.0,basis_index) == 1.0) &&
        (a->get_type1_value(-0.5,basis_index) == 1.0) &&
        (a->get_type1_value(0.5,basis_index) == 1.0) &&
        (a->get_type1_value(-1.0,basis_index) == 1.0) && 
        (a->get_type1_value(1.0,basis_index) == 1.0) &&
        (a->get_type1_value(-1.0001,basis_index) == 0.0) &&
	(a->get_type1_value(1.0001,basis_index) == 0.0) ) ){
    std::cout << "Test fail: first basis function value";
    return EXIT_FAILURE;
  }
  
  //Test first functions gradients.  Should be zero everywhere
  if( !((a->get_type1_gradient(0.0,basis_index) == 0.0) &&
        (a->get_type1_gradient(-0.5,basis_index) == 0.0) &&
	(a->get_type1_gradient(0.5,basis_index) == 0.0) &&
	(a->get_type1_gradient(-1.0,basis_index) == 0.0) && 
	(a->get_type1_gradient(1.0,basis_index) == 0.0) &&
	(a->get_type1_gradient(-1.0001,basis_index) == 0.0) &&
	(a->get_type1_gradient(1.0001,basis_index) == 0.0) ) ){
    std::cout <<"Test fail: first basis function gradient";
    return EXIT_FAILURE;
  }

  // Test the left boundary basis function.
  basis_index[0] = 2;
  basis_index[1] = 0;
  if( !( (a->get_type1_value(0.0,basis_index) == 0.0) &&
	 (a->get_type1_value(-0.5,basis_index) == 0.5) &&
	 (a->get_type1_value(0.5,basis_index) == 0.0) &&
	 (a->get_type1_value(-1.0,basis_index) == 1.0) && 
	 (a->get_type1_value(1.0,basis_index) == 0.0) &&
	 (a->get_type1_value(-1.0001,basis_index) == 0.0) &&
	 (a->get_type1_value(1.0001,basis_index) == 0.0) ) ){
    std::cout << "Test fail: left boundary basis function value";
    return EXIT_FAILURE;
  }
  //Test functions gradients.
  if( !( (a->get_type1_gradient(0.0,basis_index) == 0.0) &&
	 (a->get_type1_gradient(-0.5,basis_index) == -1.0) &&
	 (a->get_type1_gradient(-0.2,basis_index) == -1.0) &&
	 (a->get_type1_gradient(0.5,basis_index) == 0.0) &&
	 (a->get_type1_gradient(-1.0,basis_index) == 0.0) && 
	 (a->get_type1_gradient(1.0,basis_index) == 0.0) &&
	 (a->get_type1_gradient(-1.0001,basis_index) == 0.0) &&
	 (a->get_type1_gradient(1.0001,basis_index) == 0.0) ) ){
    std::cout << "Test fail: left boundary basis function gradient";
    return EXIT_FAILURE;
  }
      
  //Testing Right basis element
  basis_index[0] = 2;
  basis_index[1] = 1;
  if ( !( (a->get_type1_value(0.0,basis_index) == 0.0) &&
	  (a->get_type1_value(-0.5,basis_index) == 0.0) &&
	  (a->get_type1_value(0.5,basis_index) == 0.5) &&
	  (a->get_type1_value(-1.0,basis_index) == 0.0) && 
	  (a->get_type1_value(1.0,basis_index) == 1.0) &&
	  (a->get_type1_value(-1.0001,basis_index) == 0.0) && 
	  (a->get_type1_value(1.0001,basis_index) == 0.0) ) ){
    std::cout << "Test fail: right boundary basis function value";
    return EXIT_FAILURE;
  }
  
  //Test first functions gradients.  should be 1 on the suppport
  // and zero off.
  if ( !( (a->get_type1_gradient(0.0,basis_index) == 0.0) &&
	  (a->get_type1_gradient(-0.5,basis_index) == 0.0) &&
	  (a->get_type1_gradient(0.5,basis_index) == 1.0) &&
	  (a->get_type1_gradient(-1.0,basis_index) == 0.0) && 
	  (a->get_type1_gradient(1.0,basis_index) == 0.0) &&
	  (a->get_type1_gradient(-1.0001,basis_index) == 0.0) &&
	  (a->get_type1_gradient(1.0001,basis_index) == 0.0) ) ){
    std::cout << "Test fail: right boundary basis function grad";
    return EXIT_FAILURE;
  }
  //Test the 3rd function on level 4. Support is [0,.5] with peak at
  //0.25
  basis_index[0] = 4;
  basis_index[1] = 2;
  if( !( (a->get_type1_value(0.0,basis_index) == 0.0) &&
	 (a->get_type1_value(0.5,basis_index) == 0.0) &&
	 (a->get_type1_value(0.25,basis_index) == 1.0) &&
	 (a->get_type1_value(0.25 - 0.25/2,basis_index) == 0.5) && 
	 (a->get_type1_value(0.25 + 0.25/2,basis_index) == 0.5) &&
	 (a->get_type1_value(-0.0001,basis_index) == 0.0) &&
	 (a->get_type1_value(0.50001,basis_index) == 0.0) ) ){
    std::cout << "Test fail: [4,3] basis function value";
    return EXIT_FAILURE;
  } 
  
  if( !( (a->get_type1_gradient(0.0,basis_index) == 0.0) &&
	 (a->get_type1_gradient(0.5,basis_index) == 0.0) &&
	 (a->get_type1_gradient(0.25,basis_index) == 0.0) && 
	 (a->get_type1_gradient(0.25 - 0.25/2,basis_index) == 4.0) && 
	 (a->get_type1_gradient(0.25 + 0.25/2,basis_index) == -4.0) &&
	 (a->get_type1_gradient(-0.0001,basis_index) == 0.0) &&
	 (a->get_type1_gradient(0.5001,basis_index) == 0.0) ) ){
    std::cout << "Test fail: [4,3] basis function grad";
    return EXIT_FAILURE;
  }

  
  //Test that nonsense index values throw an exception
  bool caughtException = false;
  //First try an invalid level value
  basis_index[0] = 0;
  basis_index[1] = 0;
  try{
    a->get_type1_value(0.0,basis_index);
  } catch( const std::exception& exception ) {
    caughtException = true;
  }
  if(!caughtException){
    std::cout << "Test fail: evaluation of [0,0] basis should"
	      << "throw exception.";
    return EXIT_FAILURE;
  }
  caughtException = false;
  
  //Now try an invalid index
  basis_index[0] = 4;
  basis_index[1] = -1;
  try{
    a->get_type1_value(0.0,basis_index);
  } catch( const std::exception& exception ) {
    caughtException = true;
  }
  if(!caughtException){
    std::cout << "Test fail: evaluation of [4,-1] basis should"
	      << "throw exception.";
    return EXIT_FAILURE;
  }
  caughtException = false;

  //Try another invalid index.  This time too large
  basis_index[0] = 4;
  basis_index[1] = 4;
  try{
    a->get_type1_value(0.0,basis_index);
  } catch( const std::exception& exception ) {
    caughtException = true;
  }
  if(!caughtException){
    std::cout << "Test fail: evaluation of [4,4] basis should"
	      << "throw exception.";
    return EXIT_FAILURE;
  }

  delete a;

  //getInterpPoint()
  a = new HierarchicalLinearBasis(-2,6);
  basis_index[0] = 5;
  basis_index[1] = 0;
  if( !( a->get_interp_point(basis_index) == -1.5) ){
    std::cout << "test fail: first getInterpPoint test.";
    return EXIT_FAILURE;
  }
  basis_index[0] = 5;
  basis_index[1] = 7;
  if( !( a->get_interp_point(basis_index) == 5.5) ){
    std::cout << "test fail: second getInterpPoint test.";
    return EXIT_FAILURE;
  }
  basis_index[0] = 5;
  basis_index[1] = 4;
  if( !( a->get_interp_point(basis_index) == 2.5) ){
    std::cout << "test fail: third getInterpPoint test.";
    return EXIT_FAILURE;
  }

  delete a;

  //conversion utilities test
  basis_index[0] = 1;
  basis_index[1] = 0;
  if( !( HierarchicalLinearBasis::intArray_to_int(basis_index) == 0) ){
    std::cout << "test fail: first conversion test.";
    return EXIT_FAILURE;
  }
  if( !(basis_index == HierarchicalLinearBasis::int_to_intArray(0) ) ){
    std::cout << "test fail: second conversion test.";
    return EXIT_FAILURE;
  }

  basis_index[0] = 2;
  basis_index[1] = 0;
  if( !(HierarchicalLinearBasis::intArray_to_int(basis_index) == 1)){
    std::cout << "test fail: third conversion test.";
    return EXIT_FAILURE;
  }
  if(!(basis_index == HierarchicalLinearBasis::int_to_intArray(1))){
    std::cout << "test fail: fourth conversion test.";
    return EXIT_FAILURE;
  }
  basis_index[0] = 2;
  basis_index[1] = 1;
  if( !(HierarchicalLinearBasis::intArray_to_int(basis_index) == 2) ){
    std::cout << "test fail: fifth conversion test.";
    return EXIT_FAILURE;
  }
  if( !(basis_index == HierarchicalLinearBasis::int_to_intArray(2))){
    std::cout << "test fail: sixth conversion test.";
    return EXIT_FAILURE;
  }
  
  basis_index[0] = 5;
  basis_index[1] = 3;
  if(!(HierarchicalLinearBasis::intArray_to_int(basis_index) == 12)){
    std::cout << "test fail: seventh conversion test.";
    return EXIT_FAILURE;
  }
  if(!(basis_index == HierarchicalLinearBasis::int_to_intArray(12))){
    std::cout << "test fail: eigth conversion test.";
    return EXIT_FAILURE;
  }

  basis_index[0] = 6;
  basis_index[1] = 0;
  if(!(HierarchicalLinearBasis::intArray_to_int(basis_index) == 17)){
    std::cout << "test fail: ninth conversion test.";
    return EXIT_FAILURE;
  }
  if(!(basis_index == HierarchicalLinearBasis::int_to_intArray(17))){
    std::cout << "test fail: tenth conversion test.";
    return EXIT_FAILURE;
  }

  basis_index[0] = 6;
  basis_index[1] = 5;
  if(!(HierarchicalLinearBasis::intArray_to_int(basis_index) == 22)){
    std::cout << "test fail: 11th conversion test.";
    return EXIT_FAILURE;
  }
  if(!(basis_index == HierarchicalLinearBasis::int_to_intArray(22))){
    std::cout << "test fail: 12th conversion test.";
    return EXIT_FAILURE;
  }
  
  basis_index[0] = 6;
  basis_index[1] = 11;
  if(!(HierarchicalLinearBasis::intArray_to_int(basis_index) == 28)){
    std::cout << "test fail: 13th conversion test.";
    return EXIT_FAILURE;
  }
  if(!(basis_index == HierarchicalLinearBasis::int_to_intArray(28))){
    std::cout << "test fail: 14th conversion test.";
    return EXIT_FAILURE;
  }

  basis_index[0] = 6;
  basis_index[1] = 15;
  if(!(HierarchicalLinearBasis::intArray_to_int(basis_index) == 32)){
    std::cout << "test fail: 15th conversion test.";
    return EXIT_FAILURE;
  }
  if(!(basis_index == HierarchicalLinearBasis::int_to_intArray(32))){
    std::cout << "test fail: 16th conversion test.";
    return EXIT_FAILURE;
  }

  //test log
  if(log(pow(2,512))/log(2) != 512){
    std::cout << "log broke down due to rounding." << std::endl;
    return EXIT_FAILURE;
  }

  //invalid Index test
  IntArray basis_index_too_long(3);
  IntArray basis_index_too_short(1);

  // Test for exceptions on invalidly lengthed
  // indices.
  caughtException = false;
  try{
    HierarchicalLinearBasis::check_valid_index(basis_index_too_short);
  } catch( const std::length_error& exception ) {
    caughtException = true;
  }
  if(!caughtException){
    std::cout << "test fail: missed exception when vector too short.";
    return EXIT_FAILURE;
  }
  caughtException = false;

  try{
    HierarchicalLinearBasis::check_valid_index(basis_index_too_long);
  } catch( const std::length_error& exception ) {
    caughtException = true;
  }
  if(!caughtException){
    std::cout << "test fail: missed exception when vector too long.";
    return EXIT_FAILURE;
  }
  caughtException = false;

  basis_index[0] = 0;
  basis_index[1] = 0;
  try{
    HierarchicalLinearBasis::check_valid_index(basis_index);
  } catch( const std::invalid_argument& exception ) {
    caughtException = true;
  }
  if(!caughtException){
    std::cout << "test fail: missed exception for 1st invalid idx.";
    return EXIT_FAILURE;
  }
  caughtException = false;

  basis_index[0] = 2;
  basis_index[1] = 2;
  try{
    HierarchicalLinearBasis::check_valid_index(basis_index);
  } catch( const std::invalid_argument& exception ) {
    caughtException = true;
  }
  if(!caughtException){
    std::cout << "test fail: missed exception for 2nd invalid idx";
    return EXIT_FAILURE;
  }



  return EXIT_SUCCESS;
  
}

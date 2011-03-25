/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HierarchicalLinearBasisTestClass
//- Description:  Implements Unit Test Class for 1-D hierarchical piecewise linear basis
//-               
//- Owner:        Christopher Miller: University of Maryland at College Park
//- Contact:      cmiller@math.umd.edu

#include "HierarchicalLinearBasisTestClass.hpp"
#include <cppunit/TestCaller.h>
#include <cppunit/TestAssert.h>

using namespace Pecos;

void HierarchicalLinearBasisTestClass::setup()
{
}

void HierarchicalLinearBasisTestClass::tearDown()
{
}

void HierarchicalLinearBasisTestClass::creationTest()
{
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
  CPPUNIT_ASSERT(caughtException1);
  CPPUNIT_ASSERT(caughtException2);
  CPPUNIT_ASSERT(c == NULL);
  CPPUNIT_ASSERT(d == NULL);
}

void HierarchicalLinearBasisTestClass::destructionTest()
{
  delete a,b,c,d;
}

void HierarchicalLinearBasisTestClass::getAttributeTest()
{
  a = new HierarchicalLinearBasis();
  b = new HierarchicalLinearBasis(-1.2,5.6);

  CPPUNIT_ASSERT(a->get_lower_interp_bound() == 0.0);
  CPPUNIT_ASSERT(a->get_upper_interp_bound() == 1.0);
  CPPUNIT_ASSERT(b->get_lower_interp_bound() == -1.2);
  CPPUNIT_ASSERT(b->get_upper_interp_bound() == 5.6);

  delete a,b;

}

void HierarchicalLinearBasisTestClass::evaluationTest()
{
  a = new HierarchicalLinearBasis(-1,1);
  IntArray basis_index(2);
  //Test the highest basis function.  Should be constant 1
  // on the interval and zero off of it.
  basis_index[0] = 1;
  basis_index[1] = 0;
  CPPUNIT_ASSERT((a->get_type1_value(0.0,basis_index)) == 1.0);
  CPPUNIT_ASSERT(a->get_type1_value(-0.5,basis_index) == 1.0);
  CPPUNIT_ASSERT(a->get_type1_value(0.5,basis_index) == 1.0);
  CPPUNIT_ASSERT(a->get_type1_value(-1.0,basis_index) == 1.0 && 
		 a->get_type1_value(1.0,basis_index) == 1.0);
  CPPUNIT_ASSERT(a->get_type1_value(-1.0001,basis_index) == 0.0 &&
		 a->get_type1_value(1.0001,basis_index) == 0.0);
  
  //Test first functions gradients.  Should be zero everywhere
  CPPUNIT_ASSERT(a->get_type1_gradient(0.0,basis_index) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-0.5,basis_index) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(0.5,basis_index) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-1.0,basis_index) == 0.0 && 
		 a->get_type1_gradient(1.0,basis_index) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-1.0001,basis_index) == 0.0 &&
		 a->get_type1_gradient(1.0001,basis_index) == 0.0);

  //Same test except using the single int index method
  unsigned short basis_int_idx = 0;
  CPPUNIT_ASSERT(a->get_type1_value(0.0, basis_int_idx) == 1.0);
  CPPUNIT_ASSERT(a->get_type1_value(-0.5, basis_int_idx) == 1.0);
  CPPUNIT_ASSERT(a->get_type1_value(0.5, basis_int_idx) == 1.0);
  CPPUNIT_ASSERT(a->get_type1_value(-1.0, basis_int_idx) == 1.0 && 
		 a->get_type1_value(1.0, basis_int_idx) == 1.0);
  CPPUNIT_ASSERT(a->get_type1_value(-1.0001, basis_int_idx) == 0.0 &&
		 a->get_type1_value(1.0001, basis_int_idx) == 0.0);
  
  CPPUNIT_ASSERT(a->get_type1_gradient(0.0,basis_int_idx ) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-0.5, basis_int_idx) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(0.5, basis_int_idx) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-1.0, basis_int_idx) == 0.0 && 
		 a->get_type1_gradient(1.0, basis_int_idx) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-1.0001, basis_int_idx) == 0.0 &&
		 a->get_type1_gradient(1.0001, basis_int_idx) == 0.0);

  // Test the left boundary basis function.
  basis_index[0] = 2;
  basis_index[1] = 0;
  CPPUNIT_ASSERT(a->get_type1_value(0.0,basis_index) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_value(-0.5,basis_index) == 0.5);
  CPPUNIT_ASSERT(a->get_type1_value(0.5,basis_index) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_value(-1.0,basis_index) == 1.0 && 
		 a->get_type1_value(1.0,basis_index) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_value(-1.0001,basis_index) == 0.0 &&
		 a->get_type1_value(1.0001,basis_index) == 0.0);
  
  //Test functions gradients.
  CPPUNIT_ASSERT(a->get_type1_gradient(0.0,basis_index) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-0.5,basis_index) == -1.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-0.2,basis_index) == -1.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(0.5,basis_index) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-1.0,basis_index) == 0.0 && 
		 a->get_type1_gradient(1.0,basis_index) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-1.0001,basis_index) == 0.0 &&
		 a->get_type1_gradient(1.0001,basis_index) == 0.0);

  //Same test except using the single int index method
  basis_int_idx = 1;
  CPPUNIT_ASSERT(a->get_type1_value(0.0, basis_int_idx) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_value(-0.5, basis_int_idx) == 0.5);
  CPPUNIT_ASSERT(a->get_type1_value(0.5, basis_int_idx) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_value(-1.0, basis_int_idx) == 1.0 && 
		 a->get_type1_value(1.0, basis_int_idx) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_value(-1.0001, basis_int_idx) == 0.0 &&
		 a->get_type1_value(1.0001, basis_int_idx) == 0.0);
  
  CPPUNIT_ASSERT(a->get_type1_gradient(0.0, basis_int_idx) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-0.5, basis_int_idx) == -1.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-0.2,basis_int_idx) == -1.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(0.5, basis_int_idx) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-1.0, basis_int_idx) == 0.0 && 
		 a->get_type1_gradient(1.0, basis_int_idx) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-1.0001, basis_int_idx) == 0.0 &&
		 a->get_type1_gradient(1.0001, basis_int_idx) == 0.0);

  //Testing Right basis element
  basis_index[0] = 2;
  basis_index[1] = 1;
  CPPUNIT_ASSERT(a->get_type1_value(0.0,basis_index) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_value(-0.5,basis_index) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_value(0.5,basis_index) == 0.5);
  CPPUNIT_ASSERT(a->get_type1_value(-1.0,basis_index) == 0.0 && 
		 a->get_type1_value(1.0,basis_index) == 1.0);
  CPPUNIT_ASSERT(a->get_type1_value(-1.0001,basis_index) == 0.0 &&
		 a->get_type1_value(1.0001,basis_index) == 0.0);
  
  //Test first functions gradients.  should be 1 on the suppport
  // and zero off.
  CPPUNIT_ASSERT(a->get_type1_gradient(0.0,basis_index) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-0.5,basis_index) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(0.5,basis_index) == 1.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-1.0,basis_index) == 0.0 && 
		 a->get_type1_gradient(1.0,basis_index) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-1.0001,basis_index) == 0.0 &&
		 a->get_type1_gradient(1.0001,basis_index) == 0.0);

  //Same test except using the single int index method
  basis_int_idx = 2;
  CPPUNIT_ASSERT(a->get_type1_value(0.0, basis_int_idx) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_value(-0.5, basis_int_idx) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_value(0.5, basis_int_idx) == 0.5);
  CPPUNIT_ASSERT(a->get_type1_value(-1.0, basis_int_idx) == 0.0 && 
		 a->get_type1_value(1.0, basis_int_idx) == 1.0);
  CPPUNIT_ASSERT(a->get_type1_value(-1.0001, basis_int_idx) == 0.0 &&
		 a->get_type1_value(1.0001, basis_int_idx) == 0.0);
  
  CPPUNIT_ASSERT(a->get_type1_gradient(0.0, basis_int_idx) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-0.5, basis_int_idx) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(0.5, basis_int_idx) == 1.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-1.0, basis_int_idx) == 0.0 && 
		 a->get_type1_gradient(1.0, basis_int_idx) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-1.0001, basis_int_idx) == 0.0 &&
		 a->get_type1_gradient(1.0001, basis_int_idx) == 0.0);


  //Test the 3rd function on level 4. Support is [0,.5] with peak at
  //0.25
  basis_index[0] = 4;
  basis_index[1] = 2;
  CPPUNIT_ASSERT(a->get_type1_value(0.0,basis_index) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_value(0.5,basis_index) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_value(0.25,basis_index) == 1.0);
  CPPUNIT_ASSERT(a->get_type1_value(0.25 - 0.25/2,basis_index) == 0.5 && 
		 a->get_type1_value(0.25 + 0.25/2,basis_index) == 0.5);
  CPPUNIT_ASSERT(a->get_type1_value(-0.0001,basis_index) == 0.0 &&
		 a->get_type1_value(0.50001,basis_index) == 0.0);
  
  CPPUNIT_ASSERT(a->get_type1_gradient(0.0,basis_index) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(0.5,basis_index) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(0.25,basis_index) == 0.0);
  
  CPPUNIT_ASSERT(a->get_type1_gradient(0.25 - 0.25/2,basis_index) == 4.0 && 
		 a->get_type1_gradient(0.25 + 0.25/2,basis_index) == -4.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-0.0001,basis_index) == 0.0 &&
		 a->get_type1_gradient(0.5001,basis_index) == 0.0);

  //Same test except using the single int index method
  basis_int_idx = 7;
  CPPUNIT_ASSERT(a->get_type1_value(0.0,basis_int_idx) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_value(0.5,basis_int_idx) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_value(0.25,basis_int_idx) == 1.0);
  CPPUNIT_ASSERT(a->get_type1_value(0.25 - 0.25/2,basis_int_idx) == 0.5 && 
		 a->get_type1_value(0.25 + 0.25/2,basis_int_idx) == 0.5);
  CPPUNIT_ASSERT(a->get_type1_value(-0.0001,basis_int_idx) == 0.0 &&
		 a->get_type1_value(0.50001,basis_int_idx) == 0.0);
  
  //Test first functions gradients.  should be -1 on the suppport
  // and zero off.
  CPPUNIT_ASSERT(a->get_type1_gradient(0.0,basis_int_idx) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(0.5,basis_int_idx) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(0.25,basis_int_idx) == 0.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(0.25 - 0.25/2,basis_int_idx) == 4.0 && 
		 a->get_type1_gradient(0.25 + 0.25/2,basis_int_idx) == -4.0);
  CPPUNIT_ASSERT(a->get_type1_gradient(-0.0001,basis_int_idx) == 0.0 &&
		 a->get_type1_gradient(0.5001,basis_int_idx) == 0.0);

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
  CPPUNIT_ASSERT(caughtException);
  caughtException = false;
  
  //Now try an invalid index
  basis_index[0] = 4;
  basis_index[1] = -1;
  try{
    a->get_type1_value(0.0,basis_index);
  } catch( const std::exception& exception ) {
    caughtException = true;
  }
  CPPUNIT_ASSERT(caughtException);
  caughtException = false;

  //Try another invalid index.  This time too large
  basis_index[0] = 4;
  basis_index[1] = 4;
  try{
    a->get_type1_value(0.0,basis_index);
  } catch( const std::exception& exception ) {
    caughtException = true;
  }
  CPPUNIT_ASSERT(caughtException);

  delete a;
}

void HierarchicalLinearBasisTestClass::getInterpPointTest()
{
  a = new HierarchicalLinearBasis(0,1);
  IntArray basis_index(2);
  basis_index[0] = 1;
  basis_index[1] = 0;
  CPPUNIT_ASSERT(a->get_interp_point(basis_index) == 0.5);
  basis_index[0] = 2;
  basis_index[1] = 0;
  CPPUNIT_ASSERT(a->get_interp_point(basis_index) == 0.0);
  basis_index[0] = 2;
  basis_index[1] = 1;
  CPPUNIT_ASSERT(a->get_interp_point(basis_index) == 1.0);
  basis_index[0] = 4;
  basis_index[1] = 2;
  CPPUNIT_ASSERT(a->get_interp_point(basis_index) == (0.5 + 0.75) / 2);
  
  unsigned short int_basis_idx = 0;
  CPPUNIT_ASSERT(a->get_interp_point(int_basis_idx) == 0.5);
  int_basis_idx = 1;
  CPPUNIT_ASSERT(a->get_interp_point(int_basis_idx) == 0.0);
  int_basis_idx = 2;
  CPPUNIT_ASSERT(a->get_interp_point(int_basis_idx) == 1.0);
  int_basis_idx = 7;
  CPPUNIT_ASSERT(a->get_interp_point(basis_index) == (0.5 + 0.75) / 2);
  delete a;
  
  // Try a non unit non positive example
  a = new HierarchicalLinearBasis(-2,6);
  basis_index[0] = 5;
  basis_index[1] = 0;
  CPPUNIT_ASSERT(a->get_interp_point(basis_index) == -1.5);
  basis_index[0] = 5;
  basis_index[1] = 7;
  CPPUNIT_ASSERT(a->get_interp_point(basis_index) == 5.5);
  basis_index[0] = 5;
  basis_index[1] = 4;
  CPPUNIT_ASSERT(a->get_interp_point(basis_index) == 2.5);

  //Test for exceptions on invalid indexes.
  bool caughtException = false;
  basis_index[0] = 0;
  basis_index[1] = 0;
  try{
    a->get_interp_point(basis_index);
  } catch( const std::exception& exception ) {
    caughtException = true;
  }
  CPPUNIT_ASSERT(caughtException);
  caughtException = false;

  basis_index[0] = 2;
  basis_index[1] = 2;
  try{
    a->get_interp_point(basis_index);
  } catch( const std::exception& exception ) {
    caughtException = true;
  }
  CPPUNIT_ASSERT(caughtException);
 
  delete a;
}

void HierarchicalLinearBasisTestClass::conversionUtilitiesTest()
{
  IntArray basis_index(2);
  basis_index[0] = 1;
  basis_index[1] = 0;
  CPPUNIT_ASSERT(HierarchicalLinearBasis::intArray_to_int(basis_index) == 0);
  CPPUNIT_ASSERT(basis_index == HierarchicalLinearBasis::int_to_intArray(0));

  basis_index[0] = 2;
  basis_index[1] = 0;
  CPPUNIT_ASSERT(HierarchicalLinearBasis::intArray_to_int(basis_index) == 1);
  CPPUNIT_ASSERT(basis_index == HierarchicalLinearBasis::int_to_intArray(1));

  basis_index[0] = 2;
  basis_index[1] = 1;
  CPPUNIT_ASSERT(HierarchicalLinearBasis::intArray_to_int(basis_index) == 2);
  CPPUNIT_ASSERT(basis_index == HierarchicalLinearBasis::int_to_intArray(2));

  basis_index[0] = 5;
  basis_index[1] = 3;
  CPPUNIT_ASSERT(HierarchicalLinearBasis::intArray_to_int(basis_index) == 12);
  CPPUNIT_ASSERT(basis_index == HierarchicalLinearBasis::int_to_intArray(12));

  basis_index[0] = 6;
  basis_index[1] = 0;
  CPPUNIT_ASSERT(HierarchicalLinearBasis::intArray_to_int(basis_index) == 17);
  CPPUNIT_ASSERT(basis_index == HierarchicalLinearBasis::int_to_intArray(17));

  basis_index[0] = 6;
  basis_index[1] = 5;
  CPPUNIT_ASSERT(HierarchicalLinearBasis::intArray_to_int(basis_index) == 22);
  CPPUNIT_ASSERT(basis_index == HierarchicalLinearBasis::int_to_intArray(22));

  basis_index[0] = 6;
  basis_index[1] = 11;
  CPPUNIT_ASSERT(HierarchicalLinearBasis::intArray_to_int(basis_index) == 28);
  CPPUNIT_ASSERT(basis_index == HierarchicalLinearBasis::int_to_intArray(28));

  basis_index[0] = 6;
  basis_index[1] = 15;
  CPPUNIT_ASSERT(HierarchicalLinearBasis::intArray_to_int(basis_index) == 32);
  CPPUNIT_ASSERT(basis_index == HierarchicalLinearBasis::int_to_intArray(32));

}

void HierarchicalLinearBasisTestClass::invalidIndexTest()
{
  IntArray basis_index(2);
  IntArray basis_index_too_long(3);
  IntArray basis_index_too_short(1);
  
  // Test for exceptions on invalidly lengthed
  // indices.
  bool caughtException = false;
  try{
    HierarchicalLinearBasis::check_valid_index(basis_index_too_short);
  } catch( const std::length_error& exception ) {
    caughtException = true;
  }
  CPPUNIT_ASSERT(caughtException);
  caughtException = false;

  try{
    HierarchicalLinearBasis::check_valid_index(basis_index_too_long);
  } catch( const std::length_error& exception ) {
    caughtException = true;
  }
  CPPUNIT_ASSERT(caughtException);
  caughtException = false;

  basis_index[0] = 0;
  basis_index[1] = 0;
  try{
    HierarchicalLinearBasis::check_valid_index(basis_index);
  } catch( const std::invalid_argument& exception ) {
    caughtException = true;
  }
  CPPUNIT_ASSERT(caughtException);
  caughtException = false;

  basis_index[0] = 2;
  basis_index[1] = 2;
  try{
    HierarchicalLinearBasis::check_valid_index(basis_index);
  } catch( const std::invalid_argument& exception ) {
    caughtException = true;
  }
  CPPUNIT_ASSERT(caughtException);

}

void HierarchicalLinearBasisTestClass::copyAndAssignmentTest()
{
  // Test copy
  HierarchicalLinearBasis basis1(1.789,2.657);
  HierarchicalLinearBasis basis2(basis1);

  CPPUNIT_ASSERT(basis2.get_lower_interp_bound() == 1.789 &&
		 basis2.get_upper_interp_bound() == 2.657);

  // Test assignment
  HierarchicalLinearBasis basis3;
  basis3 = basis2;
  CPPUNIT_ASSERT(basis3.get_lower_interp_bound() == 1.789 &&
		 basis3.get_upper_interp_bound() == 2.657);

}

CppUnit::Test* HierarchicalLinearBasisTestClass::suite()
{
  CppUnit::TestSuite* suiteOfTests = 
    new CppUnit::TestSuite("HierarchicalLinearBasisTestClass");
  suiteOfTests->addTest(new CppUnit::
			TestCaller<HierarchicalLinearBasisTestClass>
			("creationTest", &HierarchicalLinearBasisTestClass::
			 creationTest) );
  suiteOfTests->addTest(new CppUnit::
			TestCaller<HierarchicalLinearBasisTestClass>
			("getAttribute", &HierarchicalLinearBasisTestClass::
			 getAttributeTest) );
  suiteOfTests->addTest(new CppUnit::
			TestCaller<HierarchicalLinearBasisTestClass>
			("evaluationTest", &HierarchicalLinearBasisTestClass::
			 evaluationTest) );
  suiteOfTests->addTest(new CppUnit::
			TestCaller<HierarchicalLinearBasisTestClass>
			("getInterpPointTest", 
			 &HierarchicalLinearBasisTestClass::
			 getInterpPointTest) );
  suiteOfTests->addTest(new CppUnit::
			TestCaller<HierarchicalLinearBasisTestClass>
			("copyAndAssignmentTest",
			 &HierarchicalLinearBasisTestClass::
			 copyAndAssignmentTest) );
  suiteOfTests->addTest(new CppUnit::
			TestCaller<HierarchicalLinearBasisTestClass>
			("conversionUtilitiesTest",
			 &HierarchicalLinearBasisTestClass::
			 conversionUtilitiesTest) );
  suiteOfTests->addTest(new CppUnit::
			TestCaller<HierarchicalLinearBasisTestClass>
			("destructionTest", &HierarchicalLinearBasisTestClass::
			 destructionTest) );
  suiteOfTests->addTest(new CppUnit::
			TestCaller<HierarchicalLinearBasisTestClass>
			("invalidIndexTest", &HierarchicalLinearBasisTestClass::
			 invalidIndexTest) );
  return suiteOfTests;
}

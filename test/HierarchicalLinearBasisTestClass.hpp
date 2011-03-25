/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HierarchicalLinearBasisTestClass
//- Description:  Unit Test Class for 1-D hierarchical piecewise linear basis
//-               
//- Owner:        Christopher Miller: University of Maryland at College Park
//- Contact:      cmiller@math.umd.edu


#ifndef HIERARCHICALTEST_HPP
#define HIERARCHICALTEST_HPP

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>
#include <cppunit/Test.h>
#include "HierarchicalLinearBasis.hpp"

using namespace Pecos;
class HierarchicalLinearBasisTestClass : public CppUnit :: TestFixture
{


public:
  void setup ();
  void tearDown();
  
  void creationTest();
  void getAttributeTest();
  void evaluationTest();
  void getInterpPointTest();
  void copyAndAssignmentTest();
  void destructionTest();
  void conversionUtilitiesTest();
  void invalidIndexTest();

  static CppUnit::Test* suite();

  protected:
  
  private:
  HierarchicalLinearBasis *a, *b, *c, *d;

};

#endif

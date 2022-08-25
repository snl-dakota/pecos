/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */


//#include <ctype.h>
#include <limits>

#include <Teuchos_UnitTestHarness.hpp> 

#include "pecos_data_types.hpp"
#include "pecos_global_defs.hpp"

using namespace Pecos;

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(pecos_utils, is_small)
{
  // Behavior without reference value
  TEST_ASSERT( is_small(0.0) );

  Real test_val = SMALL_NUMBER;
  TEST_ASSERT( is_small(test_val) );

  test_val += std::numeric_limits<Real>::epsilon()*SMALL_NUMBER;
  TEST_ASSERT( !is_small(test_val) );

  test_val = (1.0-std::numeric_limits<Real>::epsilon())*SMALL_NUMBER;
  TEST_ASSERT( is_small(test_val) );
  TEST_ASSERT( is_small(-test_val) );


  // Behavior with reference value (both positive and negative)
  test_val = SMALL_NUMBER;
  TEST_ASSERT( !is_small( test_val, 0.0) ); // value is not < reference value
  TEST_ASSERT( !is_small(-test_val, 0.0) ); // abs(value) is not < reference value

  test_val = (1.0-std::numeric_limits<Real>::epsilon())*SMALL_NUMBER;
  TEST_ASSERT( is_small( test_val,  SMALL_NUMBER) );
  TEST_ASSERT( is_small(-test_val,  SMALL_NUMBER) );
  TEST_ASSERT( is_small( test_val, -SMALL_NUMBER) );
  TEST_ASSERT( is_small(-test_val, -SMALL_NUMBER) );

  Real ref_val = (1.0+std::numeric_limits<Real>::epsilon())*SMALL_NUMBER;
  TEST_ASSERT( !is_small( test_val,  ref_val) );
  TEST_ASSERT( !is_small(-test_val,  ref_val) );
  TEST_ASSERT( !is_small( test_val, -ref_val) );
  TEST_ASSERT( !is_small(-test_val, -ref_val) );

  // Not quite small enough when using very small reference value
  test_val = 10.0*SMALL_NUMBER_SQ;
  TEST_ASSERT( !is_small( test_val,  ref_val) );
  TEST_ASSERT( !is_small(-test_val,  ref_val) );
  TEST_ASSERT( !is_small( test_val, -ref_val) );
  TEST_ASSERT( !is_small(-test_val, -ref_val) );

  // ... but now the values are small enough 
  test_val = (1.0-std::numeric_limits<Real>::epsilon())*SMALL_NUMBER_SQ;
  TEST_ASSERT( is_small( test_val,  ref_val) );
  TEST_ASSERT( is_small(-test_val,  ref_val) );
  TEST_ASSERT( is_small( test_val, -ref_val) );
  TEST_ASSERT( is_small(-test_val, -ref_val) );
}

/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */


//#include <ctype.h>
#include <limits>

#define BOOST_TEST_MODULE pecos_utils
#include <boost/test/included/unit_test.hpp>

#include "pecos_data_types.hpp"
#include "pecos_global_defs.hpp"

using namespace Pecos;

//----------------------------------------------------------------

BOOST_AUTO_TEST_CASE(test_pecos_utils_is_small)
{
  // Test zero without reference value
  BOOST_CHECK( is_small(0.0) );

  // ... and with zero reference value
  BOOST_CHECK( is_small(0.0, 0.0) );
  Real test_val = SMALL_NUMBER;
  BOOST_CHECK( is_small(test_val, 0.0) );
  test_val += std::numeric_limits<Real>::epsilon()*SMALL_NUMBER;
  BOOST_CHECK( !is_small(test_val, 0.0) );

  // Test without reference value
  test_val = SMALL_NUMBER;
  BOOST_CHECK( is_small(test_val) );
  test_val += std::numeric_limits<Real>::epsilon()*SMALL_NUMBER;
  BOOST_CHECK( !is_small(test_val) );

  test_val = (1.0-std::numeric_limits<Real>::epsilon())*SMALL_NUMBER;
  BOOST_CHECK( is_small(test_val) );
  BOOST_CHECK( is_small(-test_val) );


  // Behavior with reference value (both positive and negative)
  test_val = SMALL_NUMBER;
  BOOST_CHECK( is_small( test_val, 0.0) ); // NOTE that this allows     value  > reference value
  BOOST_CHECK( is_small(-test_val, 0.0) ); // NOTE that this allows abs(value) > reference value

  test_val = (1.0-std::numeric_limits<Real>::epsilon())*SMALL_NUMBER;
  BOOST_CHECK( is_small( test_val,  SMALL_NUMBER) );
  BOOST_CHECK( is_small(-test_val,  SMALL_NUMBER) );
  BOOST_CHECK( is_small( test_val, -SMALL_NUMBER) );
  BOOST_CHECK( is_small(-test_val, -SMALL_NUMBER) );

  Real ref_val = (1.0+std::numeric_limits<Real>::epsilon())*SMALL_NUMBER;
  BOOST_CHECK( !is_small( test_val,  ref_val) );
  BOOST_CHECK( !is_small(-test_val,  ref_val) );
  BOOST_CHECK( !is_small( test_val, -ref_val) );
  BOOST_CHECK( !is_small(-test_val, -ref_val) );

  // Not quite small enough when using very small reference value
  test_val = 10.0*SMALL_NUMBER_SQ;
  BOOST_CHECK( !is_small( test_val,  ref_val) );
  BOOST_CHECK( !is_small(-test_val,  ref_val) );
  BOOST_CHECK( !is_small( test_val, -ref_val) );
  BOOST_CHECK( !is_small(-test_val, -ref_val) );

  // ... but now the values are small enough 
  test_val = (1.0-std::numeric_limits<Real>::epsilon())*SMALL_NUMBER_SQ;
  BOOST_CHECK( is_small( test_val,  ref_val) );
  BOOST_CHECK( is_small(-test_val,  ref_val) );
  BOOST_CHECK( is_small( test_val, -ref_val) );
  BOOST_CHECK( is_small(-test_val, -ref_val) );
}

//----------------------------------------------------------------

BOOST_AUTO_TEST_CASE(test_pecos_utils_is_small_sq)
{
  // Behavior without reference value
  BOOST_CHECK( is_small_sq(0.0) );

  // ... and with zero reference value
  BOOST_CHECK( is_small_sq(0.0, 0.0) );
}

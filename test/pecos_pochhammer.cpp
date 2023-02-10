/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */


#include <ctype.h>
#include <string>

#define BOOST_TEST_MODULE pecos_pochhammer
#include <boost/test/included/unit_test.hpp>

#include "pecos_data_types.hpp"
#include "BasisPolynomial.hpp"

using namespace Pecos;

//----------------------------------------------------------------

BOOST_AUTO_TEST_CASE(test_pochhammer_simple1)
{

  Real testval = -5.0;
  std::vector<Real> values;
  for( unsigned short n = 0; n < 10; ++n ) {
    values.push_back(BasisPolynomial::pochhammer(testval, n));
    //std::cout << n << "\t" << values[n] << std::endl;
  }

  BOOST_CHECK(    1 ==  values[0] );
  BOOST_CHECK(   -5 ==  values[1] );
  BOOST_CHECK(   20 ==  values[2] );
  BOOST_CHECK(  -60 ==  values[3] );
  BOOST_CHECK(  120 ==  values[4] );
  BOOST_CHECK( -120 ==  values[5] );
  BOOST_CHECK(    0 ==  values[6] );
  BOOST_CHECK(    0 ==  values[7] );
  BOOST_CHECK(    0 ==  values[8] );
  BOOST_CHECK(    0 ==  values[9] );
}

//----------------------------------------------------------------

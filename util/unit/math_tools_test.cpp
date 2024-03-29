/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */


#include <string>

#define BOOST_TEST_MODULE pecos_math_tools
#include <boost/test/included/unit_test.hpp>

#include "math_tools.hpp"

using namespace Pecos;
using namespace Pecos::util;

namespace {
  const int MAX_DEGREE = 10;
  const int MAX_DIM = 10;
}

//----------------------------------------------------------------

BOOST_AUTO_TEST_CASE(test_math_tools_multi_index)
{
  IntMatrix indices;

  for(int dim = 1; dim < MAX_DIM+1; ++dim)
  {
    IntMatrix ones(1,dim);
    ones.putScalar(1);
    for(int degree = 1; degree < MAX_DEGREE+1; ++degree)
    {
      // What we are testing
      // ----------------------------
      Pecos::util::get_multi_dimensional_polynomial_indices( dim, degree, indices );
      // ----------------------------

      // Test structure of multi-index array
      int num_indices = Pecos::util::nchoosek((dim+degree-1),(dim-1));
      BOOST_CHECK( dim == indices.numRows() );
      BOOST_CHECK( num_indices == indices.numCols() );

      // Test contents of multi-index array
      IntMatrix test_sums(1,num_indices);
      test_sums.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, ones, indices, 0);
      for( int i=0; i<num_indices; ++i )
        BOOST_CHECK( degree == test_sums(0,i) );
    }
  }
}

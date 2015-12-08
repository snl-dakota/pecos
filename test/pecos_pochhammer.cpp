
#include <ctype.h>
#include <string>

#include <Teuchos_UnitTestHarness.hpp> 

#include "pecos_data_types.hpp"
#include "BasisPolynomial.hpp"

using namespace Pecos;

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(pochhammer, simple1)
{

  Real testval = 1.23;
  for( unsigned short n = 0; n < 10; ++n )
    std::cout << n << "\t" << BasisPolynomial::pochhammer(testval, n) << std::endl;

  TEST_ASSERT( true );
}

//----------------------------------------------------------------

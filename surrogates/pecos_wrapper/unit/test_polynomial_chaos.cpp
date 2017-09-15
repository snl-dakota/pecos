#include "Teuchos_UnitTestHarness.hpp"
#include "PolynomialChaosExpansionWrap.hpp"

using namespace Surrogates;

namespace {

bool test_polynomial_approximation(short basis_type, Real tol){
  PolynomialChaosExpansionWrap poly;
  return true;
}

TEUCHOS_UNIT_TEST(pecos_polynomial_approximation, pce)
{
  TEST_ASSERT(test_polynomial_approximation(Pecos::LEGENDRE_ORTHOG, 1e-8));
}

} // namespace

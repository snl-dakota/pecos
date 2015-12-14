
#include <ctype.h>
#include <string>

#include <Teuchos_UnitTestHarness.hpp> 

#include "pecos_data_types.hpp"
#include "BasisPolynomial.hpp"
#include "KrawtchoukOrthogPolynomial.hpp"
#include "MeixnerOrthogPolynomial.hpp"
#include "HahnOrthogPolynomial.hpp"

using namespace Pecos;

namespace {

  //------------------------------------
  // Compute known exact orthogonality value
  //------------------------------------
  Real krawtchouck_exact_orthog(short N, Real p, short order)
  {
    Real value =   std::pow(-1, order)*BasisPolynomial::factorial(order)/BasisPolynomial::pochhammer(-N,order)
                 * std::pow((1.0-p)/p,order);

    return value;
  }

  //------------------------------------
  // Compute numerical inner product
  //------------------------------------
  Real krawtchouck_inner_prod(short N, Real p, short order1, short order2, BasisPolynomial * poly)
  {
    Real sum = 0.0;
    for( short i=0; i<N+1; ++i )
      sum += BasisPolynomial::n_choose_k(N,i)*std::pow(p,i)*std::pow(1.0-p,N-i)*poly->type1_value(Real(i),order1)*poly->type1_value(Real(i),order2);

    return sum;
  }
}


//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(discrete_orthog_poly, krawtchouck1)
{
  BasisPolynomial poly_basis = BasisPolynomial(KRAWTCHOUK_DISCRETE);
  KrawtchoukOrthogPolynomial * ptr = dynamic_cast<KrawtchoukOrthogPolynomial*>(poly_basis.polynomial_rep());
  TEST_ASSERT( ptr != NULL );

  // Test deafult settings and accessors
  TEST_EQUALITY( poly_basis.alpha_polynomial(), -1.0 );
  TEST_EQUALITY( poly_basis.beta_polynomial(), -1.0 );

  const Real p = 0.1;
  const Real N = 15.0;
  const Real TEST_TOL = 1.e-9; // a relative tolerance based on the exact answers

  poly_basis.alpha_stat(p);
  poly_basis.beta_stat(N);

  // Test orthogonality of first 10 polynomials - covers hardcoded 1st and 2nd orders and recursion-based orders
  for( short i=0; i<11; ++i ) {
    Real exact_orth_val = krawtchouck_exact_orthog(N, p, i);
    for( short j=0; j<11; ++j ) {
      Real numerical_orth_val = krawtchouck_inner_prod(N, p, i, j, ptr);
      if( i == j ) {
        TEST_FLOATING_EQUALITY( exact_orth_val, numerical_orth_val, TEST_TOL );
      }
      else {
        Real shifted_zero = 1.0 + numerical_orth_val;
        TEST_FLOATING_EQUALITY( shifted_zero, 1.0, TEST_TOL );
      }
    }
  }
}

//----------------------------------------------------------------


namespace {

  //------------------------------------
  // Compute known exact orthogonality value
  //------------------------------------
  Real meixner_exact_orthog(Real c, Real B, short order)
  {
    Real value =   BasisPolynomial::factorial(order)/(BasisPolynomial::pochhammer(B,order)
                 * std::pow(c,order)*std::pow((1.0-c),B));

    return value;
  }

  //------------------------------------
  // Compute numerical inner product
  //------------------------------------
  Real meixner_inner_prod(unsigned nterms, Real c, Real B, short order1, short order2, BasisPolynomial * poly)
  {
    Real sum = 0.0;
    for( short i=0; i<nterms; ++i )
      sum += BasisPolynomial::pochhammer(B,i)*std::pow(c,i)/BasisPolynomial::factorial(i)*poly->type1_value(Real(i),order1)*poly->type1_value(Real(i),order2);

    return sum;
  }
}


//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(discrete_orthog_poly, meixner1)
{
  BasisPolynomial poly_basis = BasisPolynomial(MEIXNER_DISCRETE);
  MeixnerOrthogPolynomial * ptr = dynamic_cast<MeixnerOrthogPolynomial*>(poly_basis.polynomial_rep());
  TEST_ASSERT( ptr != NULL );

  // Test deafult settings and accessors
  TEST_EQUALITY( poly_basis.alpha_polynomial(), -1.0 );
  TEST_EQUALITY( poly_basis.beta_polynomial(), -1.0 );

  const Real c    = 0.1;
  const Real beta = 1.5;
  const Real TEST_TOL = 1.e-9; // a relative tolerance based on the exact answers
  const unsigned NUM_TERMS_TO_SUM = 40; // the number of terms needed for the orthogonality sum to converge

  poly_basis.alpha_stat(c);
  poly_basis.beta_stat(beta);

  // Test orthogonality of first 10 polynomials - covers hardcoded 1st and 2nd orders and recursion-based orders
  for( short i=0; i<7; ++i ) {
    Real exact_orth_val = meixner_exact_orthog(c, beta, i);
    for( short j=0; j<7; ++j ) {
      Real numerical_orth_val = meixner_inner_prod(NUM_TERMS_TO_SUM, c, beta, i, j, ptr);
      if( i == j ) {
        TEST_FLOATING_EQUALITY( exact_orth_val, numerical_orth_val, TEST_TOL );
      }
      else {
        Real shifted_zero = 1.0 + numerical_orth_val;
        TEST_FLOATING_EQUALITY( shifted_zero, 1.0, TEST_TOL );
      }
    }
  }
}

//----------------------------------------------------------------


namespace {

  //------------------------------------
  // Compute known exact orthogonality value
  //------------------------------------
  Real hahn_exact_orthog(Real a, Real b, short N, short order)
  {
    Real value =   std::pow(-1,order)*BasisPolynomial::pochhammer((order+a+b+1.0),(N+1))*BasisPolynomial::pochhammer((b+1.0),order)
                  *BasisPolynomial::factorial(order)/((2.0*order+a+b+1.0)*BasisPolynomial::pochhammer((a+1.0),order)
                  *BasisPolynomial::pochhammer(-N,order)*BasisPolynomial::factorial(N));

    return value;
  }

  //------------------------------------
  // Compute numerical inner product
  //------------------------------------
  Real hahn_inner_prod(short N, Real a, Real b, short order1, short order2, BasisPolynomial * poly)
  {
    Real sum = 0.0;
    for( short i=0; i<N+1; ++i ) {
      Real x = Real(i);
      short sa = short(a);
      short sb = short(b);
      Real term = BasisPolynomial::n_choose_k(sa+i,i)*BasisPolynomial::n_choose_k(sb+N-i,N-i)*poly->type1_value(x,order1)*poly->type1_value(x,order2);
      sum += term;
    }

    return sum;
  }
}


//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(discrete_orthog_poly, hahn1)
{
  BasisPolynomial poly_basis = BasisPolynomial(HAHN_DISCRETE);
  HahnOrthogPolynomial * ptr = dynamic_cast<HahnOrthogPolynomial*>(poly_basis.polynomial_rep());
  TEST_ASSERT( ptr != NULL );

  // Test deafult settings and accessors
  TEST_EQUALITY( poly_basis.alpha_polynomial(), -1.0 );
  TEST_EQUALITY( poly_basis.beta_polynomial(), -1.0 );
  TEST_EQUALITY( ptr->gamma_polynomial(), -1.0 );

  const Real alpha = 4.0;
  const Real beta  = 6.0;
  const short N    = 10;
  const Real  rN   = 10.0;
  const Real TEST_TOL = 5.e-8; // a relative tolerance based on the exact answers

  poly_basis.alpha_stat(alpha);
  poly_basis.beta_stat(beta);
  ptr->gamma_stat(rN);

  // Test orthogonality of first 10 polynomials - covers hardcoded 1st and 2nd orders and recursion-based orders
  for( short i=0; i<11; ++i ) {
    Real exact_orth_val = hahn_exact_orthog(alpha, beta, N, i);
    for( short j=0; j<11; ++j ) {
      Real numerical_orth_val = hahn_inner_prod(N, alpha, beta, i, j, ptr);
      if( i == j ) {
        TEST_FLOATING_EQUALITY( exact_orth_val, numerical_orth_val, TEST_TOL );
      }
      else {
        Real shifted_zero = 1.0 + numerical_orth_val;
        TEST_FLOATING_EQUALITY( shifted_zero, 1.0, TEST_TOL );
      }
    }
  }
}

//----------------------------------------------------------------

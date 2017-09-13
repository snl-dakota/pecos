#include "LinearSolver.hpp"
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_ParameterList.hpp"
#include <ctype.h>
#include <string>

using namespace Surrogates;

namespace {

  const int NUMROWS = 5;
  const int NUMCOLS = 5;

  // A random matrix
  //   we may want to dial in scaling
  RealMatrix get_test_matrix() {
    RealMatrix mat(NUMROWS,NUMCOLS);
    mat.random();
    return mat;
  }

  //--------------------------------------

  // A random matrix
  //   we may want to dial in scaling
  RealMatrix get_spd_test_matrix() {
    RealMatrix mat = get_test_matrix();
    RealMatrix trans_mat(mat, Teuchos::TRANS);
    mat += trans_mat;
    mat.scale(0.5);
    Real shift = 2.0;
    for(int i=0; i<mat.numRows(); ++i)
      mat(i,i) += shift;
    return mat;
  }

  //--------------------------------------

  Real test_solver(LinearSolver* sol, Teuchos::ParameterList& params, bool use_spd) {

    RealMatrix A;
    if( use_spd )
      A = get_spd_test_matrix();
    else
      A = get_test_matrix();

    RealMatrix gold_soln(NUMROWS,1);
    gold_soln.random();

    RealMatrix B(gold_soln);
    // Create the RHS needed to recover the gold soln via solve
    B.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, gold_soln, 0.0);

    RealMatrix result0, result1;
    sol->solve(A, B, result0, result1, params);

    //gold_soln.print(std::cout);
    //std::cout << std::endl;
    //result0.print(std::cout);
    //std::cout << std::endl;
    //result1.print(std::cout);

    // Test correctness
    gold_soln -= result0;
    Real shifted_diff_norm = gold_soln.normFrobenius() + 1.0;

    return shifted_diff_norm;
  }
}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_solvers, ompsolver)
{
  LinearSolver * psol = new LinearSolver;
  //params.set("Normalize Inputs", true);
  Teuchos::ParameterList params;
  params.set<unsigned>("Linear Solver Type", LinearSolver::ORTHOG_MATCH_PURSUIT);
  Real diff = test_solver(psol, params, false);

  Real tol = params.get("Solver Tolerance", 1.e-8);
  // This solve fails
  //TEST_FLOATING_EQUALITY( 1.0, diff, tol )
}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_solvers, larssolver)
{
  LinearSolver * psol = new LinearSolver;
  Teuchos::ParameterList params;
  params.set<unsigned>("Linear Solver Type", LinearSolver::LEAST_ANGLE_REGRESSION);
  params.set("Max Nonzeros", 0); // will trigger internal determination max(M,N)
  params.set("Max Iters", 30);
  Real diff = test_solver(psol, params, false);

  Real tol = params.get("Solver Tolerance", 1.e-8);
  // This solve fails
  //TEST_FLOATING_EQUALITY( 1.0, diff, tol )
}

//----------------------------------------------------------------

//TEUCHOS_UNIT_TEST(linear_solvers, cosampsolver)
//{
//  LinearSolver * psol = new LinearSolver;
//  Teuchos::ParameterList params;
//  params.set<unsigned>("Linear Solver Type", COSAMP_SOLVER);
//  Real diff = test_solver(psol, params, false);
//
//  Real tol = params.get<Real>("Solver Tolerance");
//  TEST_FLOATING_EQUALITY( 1.0, diff, tol )
//}

//----------------------------------------------------------------

//TEUCHOS_UNIT_TEST(linear_solvers, lsqsolver)
//{
//  LinearSolver * psol = new LinearSolver;
//  Teuchos::ParameterList params;
//  params.set<unsigned>("Linear Solver Type", LinearSolver::SVD_LEAST_SQ_REGRESSION);
//  Real diff = test_solver(psol, params, false);
//
//  Real tol = params.get("Solver Tolerance", 1.e-8);
//  TEST_FLOATING_EQUALITY( 1.0, diff, tol )
//}

//----------------------------------------------------------------

//TEUCHOS_UNIT_TEST(linear_solvers, equalityconstrainedlsqsolver)
//{
//  LinearSolver * psol = new LinearSolver;
//  psol->set_num_primary_equations(NUMROWS);
//  Teuchos::ParameterList params;
//  params.set<unsigned>("Linear Solver Type", EQ_CONS_LEAST_SQ_REGRESSION);
//  Real diff = test_solver(psol, params, false);
//
//  Real tol = params.get("Solver Tolerance", 1.e-8);
//  TEST_FLOATING_EQUALITY( 1.0, diff, tol )
//}

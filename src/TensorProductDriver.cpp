/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 TensorProductDriver
//- Description: Implementation code for TensorProductDriver class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "TensorProductDriver.hpp"
#include "PolynomialApproximation.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: TensorProductDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


void TensorProductDriver::
initialize_grid(const ShortArray& u_types,
		const Pecos::BasisConfigOptions& bc_options)
{
  initialize_rules(u_types, bc_options);
  quadOrder.resize(numVars); levelIndex.resize(numVars);
}


void TensorProductDriver::
initialize_grid(const std::vector<BasisPolynomial>& poly_basis)
{
  initialize_rules(poly_basis);
  quadOrder.resize(numVars); levelIndex.resize(numVars);
}


/** This function selects the smallest nested rule order that meets the
    integrand precision of a corresponding Gauss rule.  It is similar to
    the moderate exponential growth option in sparse grids. */
void TensorProductDriver::
integrand_goal_to_nested_quadrature_order(size_t i,
					  unsigned short integrand_goal,
					  unsigned short& nested_quad_order)
{
  switch (collocRules[i]) {
  case CLENSHAW_CURTIS: case NEWTON_COTES: { // closed rules
    nested_quad_order = 1; // in case while loop falls through
    unsigned short level = 0, integrand_actual = 1;
    while (integrand_actual < integrand_goal) {
      ++level;
      nested_quad_order = (unsigned short)std::pow(2., (int)level) + 1;
      // integrand exactness for CC; also used for NC
      integrand_actual = (nested_quad_order % 2) ?
	nested_quad_order : nested_quad_order - 1;
    }
    break;
  }
  case FEJER2: { // open rule (integrand exactness for CC also used for F2)
    nested_quad_order = 1; // in case while loop falls through
    unsigned short level = 0, integrand_actual = 1;
    while (integrand_actual < integrand_goal) {
      ++level;
      nested_quad_order = (unsigned short)std::pow(2., (int)level+1) - 1;
      integrand_actual = (nested_quad_order % 2) ?
	nested_quad_order : nested_quad_order - 1;
    }
    break;
  }
  case GAUSS_PATTERSON: { // open rule
    nested_quad_order = 1; // in case while loop falls through
    unsigned short level = 0, integrand_actual = 1,
      previous = nested_quad_order;
    while (integrand_actual < integrand_goal) {
      ++level;
      // exponential growth
      nested_quad_order = (unsigned short)std::pow(2., (int)level+1) - 1;
      integrand_actual = 2*nested_quad_order - previous;//2m-1 - constraints + 1
      previous = nested_quad_order;
    }
    break;
  }
  case GENZ_KEISTER: { // open rule with lookup
    unsigned short level = 0, max_level = 5;
    while (level <= max_level && precGenzKeister[level] < integrand_goal)
      ++level;
    nested_quad_order = orderGenzKeister[level];
    /*
    nested_quad_order = 1;
    unsigned short integrand_goal = 2*ref_quad_order - 1, level = 0,
      integrand_actual = 1, previous = nested_quad_order, i_rule = GENZ_KEISTER,
      g_rule = FULL_EXPONENTIAL; // map l->o without restriction
    while (integrand_actual < integrand_goal) {
      ++level;
      webbur::level_growth_to_order_new(1, &level, &i_rule, &g_rule,
	                                nested_quad_order);
      integrand_actual = 2*nested_quad_order - previous;//2m-1 - constraints + 1
      previous = nested_quad_order;
    }
    */
    break;
  }
  default: { // Gauss rules
    nested_quad_order = 1;
    unsigned short integrand_actual = 1;
    while (integrand_actual < integrand_goal) {
      // allow even quad order (won't happen for current goal definitions)
      ++nested_quad_order;
      //nested_quad_order += 2; // moderate linear growth: odd for weakly nested
      integrand_actual = 2*nested_quad_order - 1;
    }
    break;
  }
  }
}


/** This function selects the smallest nested rule order that meets the
    quadrature order goal. */
void TensorProductDriver::
quadrature_goal_to_nested_quadrature_order(size_t i, unsigned short quad_goal,
					   unsigned short& nested_quad_order)
{
  switch (collocRules[i]) {
  case CLENSHAW_CURTIS: case NEWTON_COTES: { // closed nested rules
    nested_quad_order = 1; unsigned short level = 0;
    while (nested_quad_order < quad_goal) {
      ++level;
      nested_quad_order = (unsigned short)std::pow(2., (int)level) + 1;
    }
    break;
  }
  case FEJER2: case GAUSS_PATTERSON:{ // open nested rules
    nested_quad_order = 1; unsigned short level = 0;
    while (nested_quad_order < quad_goal) {
      ++level;
      nested_quad_order = (unsigned short)std::pow(2., (int)level+1) - 1;
    }
    break;
  }
  case GENZ_KEISTER: { // open nested rule with lookup
    nested_quad_order = 1; unsigned short level = 0, max_level = 5;
    while (level <= max_level && nested_quad_order < quad_goal) {
      ++level;
      nested_quad_order = orderGenzKeister[level];
    }
    break;
  }
  default: // open weakly/non-nested Gauss rules
    nested_quad_order = quad_goal; break;
  }
}


void TensorProductDriver::compute_grid(RealMatrix& variable_sets)
{
#ifdef DEBUG
  // -----------------------------------
  // Output number of collocation points
  // -----------------------------------
  PCout << "Total number of tensor-product quadrature points: "
	<< grid_size() << '\n';
#endif // DEBUG

  // -------------------------------------------------------------------
  // Get collocation points and integration weights and update 1D arrays
  // -------------------------------------------------------------------
  compute_tensor_grid(quadOrder, levelIndex, variable_sets, type1WeightSets,
		      type2WeightSets, collocKey);
}

} // namespace Pecos

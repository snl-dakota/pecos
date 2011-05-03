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
initialize_grid(const ShortArray& u_types, bool  nested_rules,
		bool  equidistant_rules,   short growth_rate,
		short nested_uniform_rule)
{
  numVars = u_types.size();
  quadOrder.resize(numVars);
  initialize_rules(u_types, nested_rules, equidistant_rules, growth_rate,
		   nested_uniform_rule, integrationRules, growthRules);
}


void TensorProductDriver::
initialize_grid(const std::vector<BasisPolynomial>& poly_basis,
		const UShortArray& quad_order, short growth_rate)
{
  numVars         = poly_basis.size();
  polynomialBasis = poly_basis; // shallow copy
  quadOrder       = quad_order;
  initialize_rules(poly_basis, growth_rate, integrationRules, growthRules);
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

  // ----------------------------------------------
  // Get collocation points and integration weights
  // ----------------------------------------------
  compute_tensor_grid(quadOrder, variable_sets, weightSets, collocKey,
		      collocPts1D, collocWts1D);
}

} // namespace Pecos

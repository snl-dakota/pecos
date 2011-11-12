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

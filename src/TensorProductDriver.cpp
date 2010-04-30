/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
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


void TensorProductDriver::anisotropic_weights(const RealVector& aniso_wts)
{
  if (aniso_wts.empty())
    dimIsotropic = true;
  else {
    if (aniso_wts.length() != numVars) {
      PCerr << "Error: length of anisotropic weights specification is "
	    << "inconsistent with\n       number of variables in "
	    << "TensorProductDriver::anisotropic_weights()." << std::endl;
      abort_handler(-1);
    }

    dimIsotropic = true;
    Real wt0 = aniso_wts[0];
    for (size_t i=1; i<numVars; ++i)
      if (std::abs(aniso_wts[i] - wt0) > DBL_EPSILON)
	{ dimIsotropic = false; break; }

    // TO DO: scaling
    anisoLevelWts = aniso_wts;
  }
}


void TensorProductDriver::
initialize_grid(const ShortArray& u_types, bool nested_rules,
		short growth_rate, short nested_uniform_rule)
{
  numVars = u_types.size();
  quadOrder.resize(numVars);
  initialize_rules(u_types, nested_rules, growth_rate, nested_uniform_rule,
		   integrationRules, growthRules);
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


void TensorProductDriver::compute_grid()
{
  // --------------------------------
  // Get number of collocation points
  // --------------------------------
  int num_colloc_pts = grid_size();
#ifdef DEBUG
  PCout << "Total number of tensor-product quadrature points: "
	<< num_colloc_pts << '\n';
#endif // DEBUG

  // ----------------------------------------------
  // Get collocation points and integration weights
  // ----------------------------------------------
  size_t i, j;
  if (gaussPts1D.empty() || gaussWts1D.empty()) {
    gaussPts1D.resize(numVars); gaussWts1D.resize(numVars);
    for (i=0; i<numVars; ++i)
      { gaussPts1D[i].resize(1); gaussWts1D[i].resize(1); }
  }
  for (i=0; i<numVars; ++i) {
    gaussPts1D[i][0] = polynomialBasis[i].gauss_points(quadOrder[i]);
    gaussWts1D[i][0] = polynomialBasis[i].gauss_weights(quadOrder[i]);
  }
  // Tensor-product quadrature: Integral of f approximated by
  // Sum_i1 Sum_i2 ... Sum_in (w_i1 w_i2 ... w_in) f(x_i1, x_i2, ..., x_in)
  // > project 1-D gauss point arrays (of potentially different type and order)
  //   into an n-dimensional stencil.
  // > compute and store products of 1-D Gauss weights at each point in stencil.
  weightSets.sizeUninitialized(num_colloc_pts);
  variableSets.shapeUninitialized(numVars, num_colloc_pts);// Teuchos: col major
  UShortArray gauss_indices(numVars, 0);
  for (i=0; i<num_colloc_pts; i++) {
    Real& wt_prod_i = weightSets[i];
    Real* vars_i  = variableSets[i];
    wt_prod_i = 1.0;
    for (j=0; j<numVars; j++) {
      vars_i[j]  = gaussPts1D[j][0][gauss_indices[j]];
      wt_prod_i *= gaussWts1D[j][0][gauss_indices[j]];
    }
    // increment the n-dimensional gauss point index set
    PolynomialApproximation::increment_indices(gauss_indices, quadOrder, true);
  }
}

} // namespace Pecos

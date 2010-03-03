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
#include "OrthogPolyApproximation.hpp"
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
initialize_grid(const ShortArray& u_types, bool nested_rules)
{
  numVars = u_types.size();
  quadOrder.resize(numVars);
  integrationRules.resize(numVars);

  short exp_growth = MODERATE_GROWTH;//SLOW_GROWTH,STANDARD_GROWTH;
  short nested_uniform_base_rule = GAUSS_PATTERSON;//CLENSHAW_CURTIS,FEJER2
  for (size_t i=0; i<numVars; i++) {
    switch (u_types[i]) {
    case STD_NORMAL:
      integrationRules[i] = //(exp_growth == SLOW_GROWTH) ? GAUSS_HERMITE_SLOW :
	GAUSS_HERMITE; break;
    case STD_UNIFORM:
      if (nested_rules)
	switch (nested_uniform_base_rule) {
	case GAUSS_PATTERSON: // closed fully nested
	  switch (exp_growth) {
	  case SLOW_GROWTH:	// valid up to level==127 (max m=255, m=2*l+1)
	    integrationRules[i] = GAUSS_PATTERSON_SLOW;     break;
	  case MODERATE_GROWTH:	// valid up to level==63 (max m=255, m=4*l+1)
	    integrationRules[i] = GAUSS_PATTERSON_MODERATE; break;
	  case STANDARD_GROWTH:	// valid up to level==7 (max m=2^{l+1}-1 = 255)
	    integrationRules[i] = GAUSS_PATTERSON;          break;
	  } break;
	case CLENSHAW_CURTIS:
	  switch (exp_growth) {
	  case SLOW_GROWTH:
	    integrationRules[i] = CLENSHAW_CURTIS_SLOW;     break;
	  case MODERATE_GROWTH:
	    integrationRules[i] = CLENSHAW_CURTIS_MODERATE; break;
	  case STANDARD_GROWTH:
	    integrationRules[i] = CLENSHAW_CURTIS;          break;
	  } break;
	case FEJER2:
	  switch (exp_growth) {
	  case SLOW_GROWTH:
	    integrationRules[i] = FEJER2_SLOW;     break;
	  case MODERATE_GROWTH:
	    integrationRules[i] = FEJER2_MODERATE; break;
	  case STANDARD_GROWTH:
	    integrationRules[i] = FEJER2;          break;
	  } break;
	}
      else
	integrationRules[i] = GAUSS_LEGENDRE; // no SLOW if not nested_rules
      break;
    case STD_EXPONENTIAL:
      integrationRules[i] =//(exp_growth == SLOW_GROWTH) ? GAUSS_LAGUERRE_SLOW :
	GAUSS_LAGUERRE;     break;
    case STD_BETA:
      integrationRules[i] =//(exp_growth == SLOW_GROWTH) ? GAUSS_JACOBI_SLOW :
	GAUSS_JACOBI;       break;
    case STD_GAMMA:
      integrationRules[i] =
	//(exp_growth == SLOW_GROWTH) ? GEN_GAUSS_LAGUERRE_SLOW :
	GEN_GAUSS_LAGUERRE; break;
    case BOUNDED_NORMAL: case LOGNORMAL:  case BOUNDED_LOGNORMAL:
    case LOGUNIFORM:     case TRIANGULAR: case GUMBEL:
    case FRECHET:        case WEIBULL:    case HISTOGRAM_BIN:
      integrationRules[i] = //(exp_growth == SLOW_GROWTH) ? GOLUB_WELSCH_SLOW :
	GOLUB_WELSCH;       break;
    default:
      PCerr << "Error: unsupported distribution type in TensorProductDriver "
	    << "for u_type = " << u_types[i] << std::endl;
      abort_handler(-1);
      break;
    }
  }

  ShortArray basis_types, gauss_modes;
  OrthogPolyApproximation::distribution_types(u_types, integrationRules,
					      basis_types, gauss_modes);
  OrthogPolyApproximation::distribution_basis(basis_types, gauss_modes,
					      polynomialBasis);
}


void TensorProductDriver::
initialize_grid_parameters(const ShortArray& u_types, 
			   const DistributionParams& dist_params)
{
  OrthogPolyApproximation::distribution_parameters(u_types, dist_params,
						   polynomialBasis);
}


void TensorProductDriver::compute_grid()
{
  // --------------------------------
  // Get number of collocation points
  // --------------------------------
  int num_colloc_pts = grid_size();
  PCout << "Total number of tensor-product quadrature points: "
	<< num_colloc_pts << '\n';

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

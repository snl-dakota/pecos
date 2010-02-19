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
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: TensorProductDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


void TensorProductDriver::dimension_preference(const RealVector& dim_pref)
{
  RealVector aniso_wts;
  if (!dim_pref.empty()) {
    size_t num_pref = dim_pref.length();
    aniso_wts.sizeUninitialized(num_pref);
    for (size_t i=0; i<num_pref; ++i)
      aniso_wts[i] = (dim_pref[i] == 0.) ? 0. : 1./dim_pref[i];
    // scaling occurs in anisotropic_weights() below
  }
  anisotropic_weights(aniso_wts);
}


void TensorProductDriver::anisotropic_weights(const RealVector& aniso_wts)
{
  if (aniso_wts.empty())
    isotropicTPQ = true;
  else {
    if (aniso_wts.length() != numVars) {
      PCerr << "Error: length of anisotropic weights specification is "
	    << "inconsistent with\n       number of variables in "
	    << "TensorProductDriver::anisotropic_weights()." << std::endl;
      abort_handler(-1);
    }

    isotropicTPQ = true;
    Real wt0 = aniso_wts[0];
    for (size_t i=1; i<numVars; ++i)
      if (std::abs(aniso_wts[i] - wt0) > DBL_EPSILON)
	{ isotropicTPQ = false; break; }

    // TO DO: scaling
    tpqAnisoLevelWts = aniso_wts;
  }
}


void TensorProductDriver::initialize_grid_parameters(const ShortArray& u_types)
{
  integrationRules.resize(numVars);

  // isotropic uniform).  Ultimately, just want consistency: use nested rules
  // with slow exponential growth where available and non-nested rules with
  // linear growth where nested is not available.
  bool nested_rules = true;
  size_t i;
  for (i=0; i<numVars; ++i)
    if (u_types[i] != STD_UNIFORM)
      { nested_rules = false; break; }

  for (i=0; i<numVars; i++) {
    switch (u_types[i]) {
    case STD_NORMAL:
      integrationRules[i] = GAUSS_HERMITE; break;
    case STD_UNIFORM:
      // GAUSS_PATTERSON valid only up to ssgLevel==7 (max m=2^{l+1}-1 = 255)
      // GAUSS_PATTERSON_SLOW valid up to ssgLevel==127 (max m=2*l+1   = 255)
      integrationRules[i] = (nested_rules) ?
	GAUSS_PATTERSON_SLOW : GAUSS_LEGENDRE;  break;
    case STD_EXPONENTIAL:
      integrationRules[i] = GAUSS_LAGUERRE;     break;
    case STD_BETA:
      integrationRules[i] = GAUSS_JACOBI;       break;
    case STD_GAMMA:
      integrationRules[i] = GEN_GAUSS_LAGUERRE; break;
    case BOUNDED_NORMAL: case LOGNORMAL:  case BOUNDED_LOGNORMAL:
    case LOGUNIFORM:     case TRIANGULAR: case GUMBEL:
    case FRECHET:        case WEIBULL:    case HISTOGRAM_BIN:
      integrationRules[i] = GOLUB_WELSCH;       break;
    default:
      PCerr << "Error: unsupported distribution type in TensorProductDriver "
	    << "for u_type = " << u_types[i] << std::endl;
      abort_handler(-1);
      break;
    }
  }
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
  weightSets.sizeUninitialized(num_colloc_pts);
  variableSets.shapeUninitialized(numVars, num_colloc_pts);// Teuchos: col major

  // TO DO
}

} // namespace Pecos

/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 CubatureDriver
//- Description: Implementation code for CubatureDriver class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "CubatureDriver.hpp"
#include "sandia_cubature.H"
#include "OrthogPolyApproximation.hpp"
#include "DistributionParams.hpp"
#include "HermiteOrthogPolynomial.hpp"
#include "LegendreOrthogPolynomial.hpp"
#include "JacobiOrthogPolynomial.hpp"
#include "GenLaguerreOrthogPolynomial.hpp"
//#include "NumericGenOrthogPolynomial.hpp"
//#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: CubatureDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


void CubatureDriver::
initialize_grid(const ShortArray& u_types, unsigned short order,
		unsigned short rule)
{
  numVars = u_types.size();
  integrand_order(order);
  integration_rule(rule);

  ShortArray basis_types, gauss_modes;
  IntArray int_rules(numVars, (int)rule);
  OrthogPolyApproximation::distribution_types(u_types, int_rules, basis_types,
					      gauss_modes);
  OrthogPolyApproximation::distribution_basis(basis_types, gauss_modes,
					      polynomialBasis);
}


void CubatureDriver::
initialize_grid_parameters(const ShortArray& u_types,
			   const DistributionParams& dist_params)
{
  // extract alpha, beta stats: values must be homogeneous
  bool err_flag = false;
  switch (integrationRule) {
  case GAUSS_JACOBI:
    alphaStat = dist_params.beta_alpha(0);
    betaStat  = dist_params.beta_beta(0);
    // verify homogeneity
    for (size_t i=1; i<numVars; ++i)
      if (dist_params.beta_alpha(i) != alphaStat ||
	  dist_params.beta_beta(i)  != betaStat)
	err_flag = true;
    break;
  case GEN_GAUSS_LAGUERRE:
    alphaStat = dist_params.gamma_alpha(0); break;
    // verify homogeneity
    for (size_t i=1; i<numVars; ++i)
      if (dist_params.gamma_alpha(i) != alphaStat)
	err_flag = true;
  }

  if (err_flag) {
    PCerr << "Error: inhomogeneous distribution parameters in "
	  << "CubatureDriver::initialize_grid_parameters()." << std::endl;
    abort_handler(-1);
  }

  OrthogPolyApproximation::distribution_parameters(u_types, dist_params,
						   polynomialBasis);
}


int CubatureDriver::grid_size()
{
  switch(integrationRule) {
  case GAUSS_HERMITE:
    switch (integrandOrder) {
    case 1: return webbur::en_her_01_1_size(numVars);    break; // 1
    case 2: return webbur::en_her_02_xiu_size(numVars);  break; // n+1
  //case 3: return webbur::en_her_03_1_size(numVars);    break; // 2n
    case 3: return webbur::en_her_03_xiu_size(numVars);  break; // 2n
    case 5: return (numVars >=2 && numVars <= 7) ?
	webbur::en_her_05_1_size(numVars) :                     // n^2+n+2
	webbur::en_her_05_2_size(numVars);               break; // 2n^2+1
    } break;
  case GAUSS_LEGENDRE:
    switch (integrandOrder) {
    case 1: return webbur::cn_leg_01_1_size(numVars);    break; // 1
    case 2: return webbur::cn_leg_02_xiu_size(numVars);  break; // n+1
  //case 3: return webbur::cn_leg_03_1_size(numVars);    break; // 2n
    case 3: return webbur::cn_leg_03_xiu_size(numVars);  break; // 2n
    case 5: return (numVars >=4 && numVars <= 6) ?
      webbur::cn_leg_05_1_size(numVars) :                       // n^2+n+2
      webbur::cn_leg_05_2_size(numVars);                 break; // 2n^2+1
    } break;
  case GAUSS_LAGUERRE:
    switch (integrandOrder) {
    case 1: return webbur::epn_lag_01_1_size(numVars);   break;
    case 2: return webbur::epn_lag_02_xiu_size(numVars); break;
    } break;
  case GAUSS_JACOBI:
    switch (integrandOrder) {
    case 1: return webbur::cn_jac_01_1_size(numVars,   alphaStat, betaStat);
      break;
    case 2: return webbur::cn_jac_02_xiu_size(numVars, alphaStat, betaStat);
      break;
    } break;
  case GEN_GAUSS_LAGUERRE:
    switch (integrandOrder) {
    case 1: return webbur::epn_glg_01_1_size(numVars,   alphaStat); break;
    case 2: return webbur::epn_glg_02_xiu_size(numVars, alphaStat); break;
    } break;
  case GOLUB_WELSCH:
    switch (integrandOrder) {
    case 2: return webbur::gw_02_xiu_size(numVars); break;
    } break; 
  }

  PCerr << "Error: unsupported rule in CubatureDriver::grid_size()."
	<< std::endl;
  abort_handler(-1);
  return 0;
}


void CubatureDriver::compute_grid()
{
  // --------------------------------
  // Get number of collocation points
  // --------------------------------
  int num_pts = grid_size();
#ifdef DEBUG
  PCout << "Total number of cubature integration points: " << num_pts << '\n';
#endif // DEBUG

  // ----------------------------------------------
  // Get collocation points and integration weights
  // ----------------------------------------------
  weightSets.sizeUninitialized(num_pts);
  variableSets.shapeUninitialized(numVars, num_pts); // Teuchos: col major
  double *pts = variableSets.values(), *wts = weightSets.values();
  bool err_flag = false, pt_scaling = false, wt_scaling = false;
  double pt_factor, wt_factor;
  switch(integrationRule) {
  case GAUSS_HERMITE: {
    switch (integrandOrder) {
    case 1: webbur::en_her_01_1(numVars,    num_pts, pts, wts); break;
    case 2: webbur::en_her_02_xiu(numVars,  num_pts, pts, wts); break;
  //case 3: webbur::en_her_03_1(numVars,    num_pts, pts, wts); break;
    case 3: webbur::en_her_03_xiu(numVars,  num_pts, pts, wts); break;
    case 5:
      if (numVars >=2 && numVars <= 7) {
	int option = 1; // two options for n=3,5,6
	webbur::en_her_05_1(numVars, option, num_pts, pts, wts);
      }
      else
	webbur::en_her_05_2(numVars, num_pts, pts, wts);        break;
    default: err_flag = true;                                   break;
    }
    HermiteOrthogPolynomial poly;
    pt_scaling = true; pt_factor = poly.point_factor();
    wt_scaling = true; wt_factor = std::pow(poly.weight_factor(), (int)numVars);
    break;
  }
  case GAUSS_LEGENDRE: {
    switch (integrandOrder) {
    case 1: webbur::cn_leg_01_1(numVars,    num_pts, pts, wts); break;
    case 2: webbur::cn_leg_02_xiu(numVars,  num_pts, pts, wts); break;
  //case 3: webbur::cn_leg_03_1(numVars,    num_pts, pts, wts); break;
    case 3: webbur::cn_leg_03_xiu(numVars,  num_pts, pts, wts); break;
    case 5:
      if (numVars >=4 && numVars <= 6) {
	int option = 1; // two options for n=5,6
	webbur::cn_leg_05_1(numVars, option, num_pts, pts, wts);
      }
      else
	webbur::cn_leg_05_2(numVars, num_pts, pts, wts);        break;
    default: err_flag = true;                                   break;
    }
    LegendreOrthogPolynomial poly;
    wt_scaling = true; wt_factor = std::pow(poly.weight_factor(), (int)numVars);
    break;
  }
  case GAUSS_LAGUERRE:
    switch (integrandOrder) {
    case 1: webbur::epn_lag_01_1(numVars,   num_pts, pts, wts); break;
    case 2: webbur::epn_lag_02_xiu(numVars, num_pts, pts, wts); break;
    default: err_flag = true;                                   break;
    } break;
  case GAUSS_JACOBI: {
    switch (integrandOrder) {
    case 1:
      webbur::cn_jac_01_1(numVars,   alphaStat, betaStat, num_pts, pts, wts);
      break;
    case 2:
      webbur::cn_jac_02_xiu(numVars, alphaStat, betaStat, num_pts, pts, wts);
      break;
    default: err_flag = true;                            break;
    }
    JacobiOrthogPolynomial poly(alphaStat, betaStat);
    wt_scaling = true; wt_factor = std::pow(poly.weight_factor(), (int)numVars);
    break;
  }
  case GEN_GAUSS_LAGUERRE: {
    switch (integrandOrder) {
    case 1: webbur::epn_glg_01_1(numVars,   alphaStat, num_pts, pts, wts);
      break;
    case 2: webbur::epn_glg_02_xiu(numVars, alphaStat, num_pts, pts, wts);
      break;
    default: err_flag = true;                            break;
    }
    GenLaguerreOrthogPolynomial poly(alphaStat);
    wt_scaling = true; wt_factor = std::pow(poly.weight_factor(), (int)numVars);
    break;
  }
  case GOLUB_WELSCH:
    switch (integrandOrder) {
    case 2: {
      double gamma0, delta0, c1, vol_1d; // TO DO
      webbur::gw_02_xiu(numVars, num_pts, gamma0, delta0, c1, vol_1d, pts, wts);
      break;
    }
    default: err_flag = true; break;
    } break;
  default:
    err_flag = true; break;
  }

  if (err_flag) {
    PCerr << "Error: unsupported rule in CubatureDriver::compute_grid()."
	  << std::endl;
    abort_handler(-1);
  }

  // scale points and weights
  if (pt_scaling)
    variableSets.scale(pt_factor);
  if (wt_scaling)
    weightSets.scale(wt_factor);
}

} // namespace Pecos

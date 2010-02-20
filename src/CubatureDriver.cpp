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
#include "NumericGenOrthogPolynomial.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: CubatureDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


int CubatureDriver::grid_size()
{
  Real alpha, beta; // TO DO
  switch(integrationRule) {
  case GAUSS_HERMITE:
    switch (integrandPrec) {
    case 1: return webbur::en_her_01_1_size(numVars);    break; // 1
    case 2: return webbur::en_her_02_xiu_size(numVars);  break; // n+1
  //case 3: return webbur::en_her_03_1_size(numVars);    break; // 2n
    case 3: return webbur::en_her_03_xiu_size(numVars);  break; // 2n
    case 5: return //(numVars >=2 && numVars <= 7) ?
	//webbur::en_her_05_1_size(numVars) :     // n^2+n+2
	webbur::en_her_05_2_size(numVars); break; // 2n^2+1
    } break;
  case GAUSS_LEGENDRE:
    switch (integrandPrec) {
    case 1: return webbur::cn_leg_01_1_size(numVars);    break; // 1
    case 2: return webbur::cn_leg_02_xiu_size(numVars);  break; // n+1
  //case 3: return webbur::cn_leg_03_1_size(numVars);    break; // 2n
    case 3: return webbur::cn_leg_03_xiu_size(numVars);  break; // 2n
  //case 5: return (numVars >=4 && numVars <= 6) ?
      //webbur::cn_leg_05_1_size(numVars) :       // n^2+n+2
      //webbur::cn_leg_05_2_size(numVars); break; // 2n^2+1
    } break;
  case GAUSS_LAGUERRE:
    switch (integrandPrec) {
    case 1: return webbur::epn_lag_01_1_size(numVars);   break;
    case 2: return webbur::epn_lag_02_xiu_size(numVars); break;
    } break;
  case GAUSS_JACOBI:
    switch (integrandPrec) {
    case 1: return webbur::cn_jac_01_1_size(numVars,   alpha, beta); break;
    case 2: return webbur::cn_jac_02_xiu_size(numVars, alpha, beta); break;
    } break;
  case GEN_GAUSS_LAGUERRE:
    switch (integrandPrec) {
    case 1: return webbur::epn_glg_01_1_size(numVars,   alpha); break;
    case 2: return webbur::epn_glg_02_xiu_size(numVars, alpha); break;
    } break;
  //case GOLUB_WELSCH:
    // TO DO
  default:
    PCerr << "Error: unsupported rule in CubatureDriver::grid_size()."
	  << std::endl;
    abort_handler(-1);
  }
}


void CubatureDriver::compute_grid()
{
  // --------------------------------
  // Get number of collocation points
  // --------------------------------
  int num_pts = grid_size();//
  PCout << "Total number of cubature integration points: " << num_pts
	<< '\n';

  // ----------------------------------------------//
  // Get collocation points and integration weights
  // ----------------------------------------------
  weightSets.sizeUninitialized(num_pts);
  variableSets.shapeUninitialized(numVars, num_pts); // Teuchos: col major
  double *pts = variableSets.values(), *wts = weightSets.values();
  Real alpha, beta; // TO DO
  switch(integrationRule) {
  case GAUSS_HERMITE:
    switch (integrandPrec) {
    case 1: webbur::en_her_01_1(numVars,    num_pts, pts, wts); break;
    case 2: webbur::en_her_02_xiu(numVars,  num_pts, pts, wts); break;
    case 3: webbur::en_her_03_xiu(numVars,  num_pts, pts, wts); break;
    case 5: webbur::en_her_05_2(numVars,    num_pts, pts, wts); break;
    } break;
  case GAUSS_LEGENDRE:
    switch (integrandPrec) {
    case 1: webbur::cn_leg_01_1(numVars,    num_pts, pts, wts); break;
    case 2: webbur::cn_leg_02_xiu(numVars,  num_pts, pts, wts); break;
    case 3: webbur::cn_leg_03_xiu(numVars,  num_pts, pts, wts); break;
    } break;
  case GAUSS_LAGUERRE:
    switch (integrandPrec) {
    case 1: webbur::epn_lag_01_1(numVars,   num_pts, pts, wts); break;
    case 2: webbur::epn_lag_02_xiu(numVars, num_pts, pts, wts); break;
    } break;
  case GAUSS_JACOBI:
    switch (integrandPrec) {
    case 1:
      webbur::cn_jac_01_1(numVars,   alpha, beta, num_pts, pts, wts); break;
    case 2:
      webbur::cn_jac_02_xiu(numVars, alpha, beta, num_pts, pts, wts); break;
    } break;
  case GEN_GAUSS_LAGUERRE:
    switch (integrandPrec) {
    case 1: webbur::epn_glg_01_1(numVars,   alpha, num_pts, pts, wts); break;
    case 2: webbur::epn_glg_02_xiu(numVars, alpha, num_pts, pts, wts); break;
    } break;
  //case GOLUB_WELSCH:
    // TO DO
  default:
    PCerr << "Error: unsupported rule in CubatureDriver::compute_grid()."
	  << std::endl;
    abort_handler(-1);
  }
}

} // namespace Pecos

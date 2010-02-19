/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 SparseGridDriver
//- Description: Implementation code for SparseGridDriver class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "SparseGridDriver.hpp"
#include "sandia_rules.H"
#include "sparse_grid_mixed_growth.H"
#include "sgmga.H"
#include "NumericGenOrthogPolynomial.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: SparseGridDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


void SparseGridDriver::dimension_preference(const RealVector& dim_pref)
{
  RealVector aniso_wts;
  if (!dim_pref.empty()) {
    size_t num_pref = dim_pref.length();
    aniso_wts.sizeUninitialized(num_pref);
    webbur::sgmga_importance_to_aniso(num_pref, dim_pref.values(),
				      aniso_wts.values());
#ifdef DEBUG
    PCout << "dimension preference:\n"; write_data(PCout, dim_pref);
    PCout << "anisotropic weights after sgmga_importance_to_aniso():\n";
    write_data(PCout, aniso_wts);
#endif
  }
  anisotropic_weights(aniso_wts);
}


void SparseGridDriver::anisotropic_weights(const RealVector& aniso_wts)
{
  if (aniso_wts.empty())
    isotropicSSG = true;
  else {
    if (aniso_wts.length() != numVars) {
      PCerr << "Error: length of sparse grid anisotropic weights "
	    << "specification is inconsistent with\n       number of variables "
	    << "in SparseGridDriver::anisotropic_weights()." << std::endl;
      abort_handler(-1);
    }

    isotropicSSG = true;
    Real wt0 = aniso_wts[0];
    for (size_t i=1; i<numVars; ++i)
      if (std::abs(aniso_wts[i] - wt0) > DBL_EPSILON)
	{ isotropicSSG = false; break; }

    if (!isotropicSSG) {
      int option = 1; // weights scaled so that minimum nonzero entry is 1
      ssgAnisoLevelWts = aniso_wts; // copy
      webbur::sgmga_aniso_normalize(option, numVars, ssgAnisoLevelWts.values());
#ifdef DEBUG
      PCout << "ssgAnisoLevelWts after sgmga_aniso_normalize():\n";
      write_data(PCout, ssgAnisoLevelWts);
#endif
      // enforce axis lower bounds, if present, for current ssgLevel.  An axis
      // lower bound defines a weight upper bound based on the current ssgLevel:
      // LB_i = level*wt_min/wt_i --> wt_i = level*wt_min/LB_i and wt_min=1.
      // Catch special case of dim_pref_i = 0 --> wt_i = LB_i = 0.
      if (!axisLowerBounds.empty()) {
	for (size_t i=0; i<numVars; ++i)
	  if (std::abs(axisLowerBounds[i]) > 1.e-10)
	    ssgAnisoLevelWts[i] = std::min((Real)ssgLevel/axisLowerBounds[i],
					   ssgAnisoLevelWts[i]);
#ifdef DEBUG
	PCout << "ssgAnisoLevelWts after axisLowerBounds enforcement:\n";
	write_data(PCout, ssgAnisoLevelWts);
#endif
      }
    }
  }
}


void SparseGridDriver::update_axis_lower_bounds()
{
  if (axisLowerBounds.empty())
    axisLowerBounds.sizeUninitialized(numVars);
  // An axisLowerBound is the maximum index coverage achieved on a coordinate
  // axis (when all other indices are zero); it defines a constraint for
  // minimum coordinate coverage in future refinements.  The linear index set
  // constraint is level*wt_min-|wt| < j.wt <= level*wt_min, which becomes
  // level-|wt| < j_i w_i <= level for wt_min=1 and all other indices=0.
  // The max feasible j_i is then level/w_i (except for special case w_i=0).
  if (isotropicSSG)
    axisLowerBounds = (Real)ssgLevel; // all weights = 1
  else // min nonzero weight scaled to 1 --> just catch special case w_i=0
    for (size_t i=0; i<numVars; ++i)
      axisLowerBounds[i] = (std::abs(ssgAnisoLevelWts[i]) > 1.e-10) ?
	(Real)ssgLevel/ssgAnisoLevelWts[i] : 0.;
}


void SparseGridDriver::
initialize_grid_parameters(const ShortArray& u_types,
  const RealVector& nuv_means,      const RealVector& nuv_std_devs,
  const RealVector& nuv_l_bnds,     const RealVector& nuv_u_bnds,
  const RealVector& lnuv_means,     const RealVector& lnuv_std_devs,
  const RealVector& lnuv_lambdas,   const RealVector& lnuv_zetas,
  const RealVector& lnuv_err_facts, const RealVector& lnuv_l_bnds,
  const RealVector& lnuv_u_bnds,    const RealVector& luuv_l_bnds,
  const RealVector& luuv_u_bnds,    const RealVector& tuv_modes,
  const RealVector& tuv_l_bnds,     const RealVector& tuv_u_bnds,
  const RealVector& buv_alphas,     const RealVector& buv_betas,
  const RealVector& gauv_alphas,    const RealVector& guuv_alphas,
  const RealVector& guuv_betas,     const RealVector& fuv_alphas,
  const RealVector& fuv_betas,      const RealVector& wuv_alphas,
  const RealVector& wuv_betas,      const RealVectorArray& hbuv_bin_prs)
{
  integrationRules.resize(numVars);
  compute1DPoints.resize(numVars);
  compute1DWeights.resize(numVars);
  numPolyParams.assign(numVars, 0);
  size_t i, num_total_params = 0, hbuv_cntr = 0;
  for (i=0; i<numVars; i++) {
    switch (u_types[i]) {
    case STD_GAMMA:
      numPolyParams[i] = 1; break;
    case STD_BETA: case LOGNORMAL: case LOGUNIFORM: case GUMBEL: case FRECHET:
    case WEIBULL:
      numPolyParams[i] = 2; break;
    case TRIANGULAR:
      numPolyParams[i] = 3; break;
    case BOUNDED_NORMAL: case BOUNDED_LOGNORMAL:
      numPolyParams[i] = 4; break;
    case HISTOGRAM_BIN:
      numPolyParams[i] = hbuv_bin_prs[hbuv_cntr].length();
      ++hbuv_cntr;          break;
    }
    num_total_params += numPolyParams[i];
  }
  polyParams.resize(num_total_params);

  // For now, this logic is more restrictive than it needs to be (nested iff
  // isotropic uniform).  Ultimately, just want consistency: use nested rules
  // with slow exponential growth where available and non-nested rules with
  // linear growth where nested is not available.
  bool nested_rules = true;
  for (i=0; i<numVars; ++i)
    if (u_types[i] != STD_UNIFORM)
      { nested_rules = false; break; }

  size_t nuv_cntr = 0, lnuv_cntr = 0, luuv_cntr = 0, tuv_cntr = 0, buv_cntr = 0,
    gauv_cntr = 0, guuv_cntr = 0, fuv_cntr = 0, wuv_cntr = 0, pp_cntr = 0;
  hbuv_cntr = 0;
  for (i=0; i<numVars; i++) {
    switch (u_types[i]) {
    case STD_NORMAL:
      compute1DPoints[i]  = webbur::hermite_compute_points_np;
      compute1DWeights[i] = webbur::hermite_compute_weights_np;
      integrationRules[i] = GAUSS_HERMITE; // open weakly nested
      ++nuv_cntr; break;
    case BOUNDED_NORMAL:
      compute1DPoints[i]  = bounded_normal_gauss_points;
      compute1DWeights[i] = bounded_normal_gauss_weights;
      integrationRules[i] = GOLUB_WELSCH; // open non-nested
      polyParams[pp_cntr] = nuv_means[nuv_cntr];    ++pp_cntr;
      polyParams[pp_cntr] = nuv_std_devs[nuv_cntr]; ++pp_cntr;
      polyParams[pp_cntr] = nuv_l_bnds[nuv_cntr];   ++pp_cntr;
      polyParams[pp_cntr] = nuv_u_bnds[nuv_cntr];   ++pp_cntr;
      ++nuv_cntr; break;
    case LOGNORMAL:
      compute1DPoints[i]  = lognormal_gauss_points;
      compute1DWeights[i] = lognormal_gauss_weights;
      integrationRules[i] = GOLUB_WELSCH; // open non-nested
      moments_from_lognormal_spec(lnuv_means, lnuv_std_devs, lnuv_lambdas,
	lnuv_zetas, lnuv_err_facts, lnuv_cntr, polyParams[pp_cntr],
	polyParams[pp_cntr+1]);
      pp_cntr += 2;
      ++lnuv_cntr; break;
    case BOUNDED_LOGNORMAL:
      compute1DPoints[i]  = bounded_lognormal_gauss_points;
      compute1DWeights[i] = bounded_lognormal_gauss_weights;
      integrationRules[i] = GOLUB_WELSCH; // open non-nested
      moments_from_lognormal_spec(lnuv_means, lnuv_std_devs, lnuv_lambdas,
	lnuv_zetas, lnuv_err_facts, lnuv_cntr, polyParams[pp_cntr],
	polyParams[pp_cntr+1]);
      pp_cntr += 2;
      polyParams[pp_cntr] = lnuv_l_bnds[lnuv_cntr]; ++pp_cntr;
      polyParams[pp_cntr] = lnuv_u_bnds[lnuv_cntr]; ++pp_cntr;
      ++lnuv_cntr; break;
    case STD_UNIFORM:
      // For tensor-product quadrature, Gauss-Legendre is used due to greater
      // polynomial exactness since nesting is not a concern.  For nested sparse
      // grids, Clenshaw-Curtis or Gauss-Patterson can be better selections.
      // However, sparse grids that are isotropic in level but anisotropic in
      // rule become skewed when mixing Gauss rules with CC.  For this reason,
      // CC is selected only if isotropic in rule (for now).
      if (nested_rules) {
	// GAUSS_PATTERSON valid only up to ssgLevel==7 (max m=2^{l+1}-1 = 255)
	// GAUSS_PATTERSON_SLOW valid up to ssgLevel==127 (max m=2*l+1   = 255)
	compute1DPoints[i]  = webbur::patterson_lookup_points_np;
	compute1DWeights[i] = webbur::patterson_lookup_weights_np;
	integrationRules[i] = GAUSS_PATTERSON_SLOW; // closed fully nested

	//compute1DPoints[i]  = webbur::clenshaw_curtis_compute_points_np;
	//compute1DWeights[i] = webbur::clenshaw_curtis_compute_weights_np;
	//integrationRules[i] = CLENSHAW_CURTIS_SLOW; // closed fully nested

	//compute1DPoints[i]  = webbur::fejer2_compute_points_np;
	//compute1DWeights[i] = webbur::fejer2_compute_weights_np;
	//integrationRules[i] = FEJER2_SLOW; // closed fully nested
      }
      else {
	compute1DPoints[i]  = webbur::legendre_compute_points_np;
	compute1DWeights[i] = webbur::legendre_compute_weights_np;
	integrationRules[i] = GAUSS_LEGENDRE; // open weakly nested
      }
      break;
    case LOGUNIFORM:
      compute1DPoints[i]  = loguniform_gauss_points;
      compute1DWeights[i] = loguniform_gauss_weights;
      integrationRules[i] = GOLUB_WELSCH; // open non-nested
      polyParams[pp_cntr] = luuv_l_bnds[luuv_cntr]; ++pp_cntr;
      polyParams[pp_cntr] = luuv_u_bnds[luuv_cntr]; ++pp_cntr;
      ++luuv_cntr; break;
    case TRIANGULAR:
      compute1DPoints[i]  = triangular_gauss_points;
      compute1DWeights[i] = triangular_gauss_weights;
      integrationRules[i] = GOLUB_WELSCH; // open non-nested
      polyParams[pp_cntr] = tuv_modes[tuv_cntr];  ++pp_cntr;
      polyParams[pp_cntr] = tuv_l_bnds[tuv_cntr]; ++pp_cntr;
      polyParams[pp_cntr] = tuv_u_bnds[tuv_cntr]; ++pp_cntr;
      ++tuv_cntr; break;
    case STD_EXPONENTIAL:
      compute1DPoints[i]  = webbur::laguerre_compute_points_np;
      compute1DWeights[i] = webbur::laguerre_compute_weights_np;
      integrationRules[i] = GAUSS_LAGUERRE; break; // open non-nested
    case STD_BETA:
      compute1DPoints[i]  = webbur::jacobi_compute_points_np;
      compute1DWeights[i] = webbur::jacobi_compute_weights_np;
      integrationRules[i] = GAUSS_JACOBI; // open non-nested
      // convert stat to poly and assume consistent ordering via b_cntr:
      polyParams[pp_cntr] = buv_betas[buv_cntr]  - 1.; ++pp_cntr;
      polyParams[pp_cntr] = buv_alphas[buv_cntr] - 1.; ++pp_cntr;
      ++buv_cntr; break;
    case STD_GAMMA:
      compute1DPoints[i]  = webbur::gen_laguerre_compute_points_np;
      compute1DWeights[i] = webbur::gen_laguerre_compute_weights_np;
      integrationRules[i] = GEN_GAUSS_LAGUERRE; // open non-nested
      // convert stat to poly and assume consistent ordering via g_cntr:
      polyParams[pp_cntr] = gauv_alphas[gauv_cntr] - 1.; ++pp_cntr;
      ++gauv_cntr; break;
    case GUMBEL:
      compute1DPoints[i]  = gumbel_gauss_points;
      compute1DWeights[i] = gumbel_gauss_weights;
      integrationRules[i] = GOLUB_WELSCH; // open non-nested
      polyParams[pp_cntr] = guuv_alphas[guuv_cntr]; ++pp_cntr;
      polyParams[pp_cntr] = guuv_betas[guuv_cntr];  ++pp_cntr;
      ++guuv_cntr; break;
    case FRECHET:
      compute1DPoints[i]  = frechet_gauss_points;
      compute1DWeights[i] = frechet_gauss_weights;
      integrationRules[i] = GOLUB_WELSCH; // open non-nested
      polyParams[pp_cntr] = fuv_alphas[fuv_cntr]; ++pp_cntr;
      polyParams[pp_cntr] = fuv_betas[fuv_cntr];  ++pp_cntr;
      ++fuv_cntr; break;
    case WEIBULL:
      compute1DPoints[i]  = weibull_gauss_points;
      compute1DWeights[i] = weibull_gauss_weights;
      integrationRules[i] = GOLUB_WELSCH; // open non-nested
      polyParams[pp_cntr] = wuv_alphas[wuv_cntr]; ++pp_cntr;
      polyParams[pp_cntr] = wuv_betas[wuv_cntr];  ++pp_cntr;
      ++wuv_cntr; break;
    case HISTOGRAM_BIN: {
      compute1DPoints[i]  = histogram_bin_gauss_points;
      compute1DWeights[i] = histogram_bin_gauss_weights;
      integrationRules[i] = GOLUB_WELSCH; // open non-nested
      const RealVector& hbuv_bin_prs_cntr = hbuv_bin_prs[hbuv_cntr];
      size_t j, num_hbuv = hbuv_bin_prs_cntr.length();
      for (j=0; j<num_hbuv; ++j, ++pp_cntr)
	polyParams[pp_cntr] = hbuv_bin_prs_cntr[j];
      ++hbuv_cntr; break;
    }
    default:
      PCerr << "Error: unsupported distribution type in SparseGridDriver for "
	    << "u_type = " << u_types[i] << std::endl;
      abort_handler(-1);
      break;
    }
  }
}


int SparseGridDriver::grid_size()
{
  return (isotropicSSG) ?
    webbur::sparse_grid_mixed_growth_size(numVars, ssgLevel,
      &integrationRules[0], &numPolyParams[0], &polyParams[0],
      &compute1DPoints[0], duplicateTol, webbur::level_to_order_default) :
    webbur::sgmga_size(numVars, ssgAnisoLevelWts.values(), ssgLevel,
      &integrationRules[0], &numPolyParams[0], &polyParams[0],
      &compute1DPoints[0], duplicateTol, webbur::level_to_order_default);
}


int SparseGridDriver::grid_size_total()
{
  return (isotropicSSG) ?
    webbur::sparse_grid_mixed_growth_size_total(numVars, ssgLevel,
      &integrationRules[0], webbur::level_to_order_default) :
    webbur::sgmga_size_total(numVars, ssgAnisoLevelWts.values(), ssgLevel,
      &integrationRules[0], webbur::level_to_order_default);
}


void SparseGridDriver::compute_grid()
{
  // --------------------------------
  // Get number of collocation points
  // --------------------------------
  int num_total_pts = grid_size_total(), num_colloc_pts = grid_size();
  PCout << "Total number of sparse grid integration points: "
	<< num_colloc_pts << '\n';

  // ----------------------------------------------
  // Get collocation points and integration weights
  // ----------------------------------------------
  weightSets.sizeUninitialized(num_colloc_pts);
  variableSets.shapeUninitialized(numVars, num_colloc_pts);// Teuchos: col major
  uniqueIndexMapping.resize(num_total_pts);

  int* sparse_order  = new int [num_colloc_pts*numVars];
  int* sparse_index  = new int [num_colloc_pts*numVars];
  if (isotropicSSG) {
    webbur::sparse_grid_mixed_growth_unique_index(numVars, ssgLevel,
      &integrationRules[0], &numPolyParams[0], &polyParams[0],
      &compute1DPoints[0], duplicateTol, num_colloc_pts, num_total_pts,
      webbur::level_to_order_default, &uniqueIndexMapping[0]);
    webbur::sparse_grid_mixed_growth_index(numVars, ssgLevel,
      &integrationRules[0], num_colloc_pts, num_total_pts,
      &uniqueIndexMapping[0], webbur::level_to_order_default, sparse_order,
      sparse_index);
    webbur::sparse_grid_mixed_growth_weight(numVars, ssgLevel,
      &integrationRules[0], &numPolyParams[0], &polyParams[0],
      &compute1DWeights[0], num_colloc_pts, num_total_pts,
      &uniqueIndexMapping[0], webbur::level_to_order_default,
      weightSets.values());
    webbur::sparse_grid_mixed_growth_point(numVars, ssgLevel,
      &integrationRules[0], &numPolyParams[0], &polyParams[0],
      &compute1DPoints[0], num_colloc_pts, sparse_order, sparse_index,
      webbur::level_to_order_default, variableSets.values());
  }
  else {
    webbur::sgmga_unique_index(numVars, ssgAnisoLevelWts.values(), ssgLevel,
      &integrationRules[0], &numPolyParams[0], &polyParams[0],
      &compute1DPoints[0], duplicateTol, num_colloc_pts, num_total_pts,
      webbur::level_to_order_default, &uniqueIndexMapping[0]);
    webbur::sgmga_index(numVars, ssgAnisoLevelWts.values(), ssgLevel,
      &integrationRules[0], num_colloc_pts, num_total_pts,
      &uniqueIndexMapping[0], webbur::level_to_order_default, sparse_order,
      sparse_index);
    webbur::sgmga_weight(numVars, ssgAnisoLevelWts.values(), ssgLevel,
      &integrationRules[0], &numPolyParams[0], &polyParams[0],
      &compute1DWeights[0], num_colloc_pts, num_total_pts,
      &uniqueIndexMapping[0], webbur::level_to_order_default,
      weightSets.values());
    webbur::sgmga_point(numVars, ssgAnisoLevelWts.values(), ssgLevel,
      &integrationRules[0], &numPolyParams[0], &polyParams[0],
      &compute1DPoints[0], num_colloc_pts, sparse_order, sparse_index,
      webbur::level_to_order_default, variableSets.values());
  }
  delete [] sparse_order;
  delete [] sparse_index;
#ifdef DEBUG
  PCout << "uniqueIndexMapping:\n" << uniqueIndexMapping << '\n';
#endif

  // compute scale factors
  bool construct_h  = false, construct_l = false, construct_j = false,
       construct_gl = false, construct_c = false;
  size_t i, j;
  for (i=0; i<numVars; ++i) {
    switch (integrationRules[i]) {
    case GAUSS_HERMITE:
      construct_h = true;  break;
    case GAUSS_LEGENDRE: case GAUSS_PATTERSON: case GAUSS_PATTERSON_SLOW:
      construct_l = true;  break;
    case GAUSS_JACOBI:
      construct_j = true;  break;
    case GEN_GAUSS_LAGUERRE:
      construct_gl = true; break;
    case CLENSHAW_CURTIS: case CLENSHAW_CURTIS_SLOW:
    case FEJER2:          case FEJER2_SLOW:
      construct_c = true;  break;
    }
  }
  BasisPolynomial hermite_poly, legendre_poly, jacobi_poly, gen_laguerre_poly,
    chebyshev_poly;;
  if (construct_h)  hermite_poly      = BasisPolynomial(HERMITE);
  if (construct_l)  legendre_poly     = BasisPolynomial(LEGENDRE); //, mode);
  if (construct_j)  jacobi_poly       = BasisPolynomial(JACOBI);
  if (construct_gl) gen_laguerre_poly = BasisPolynomial(GENERALIZED_LAGUERRE);
  if (construct_c)  chebyshev_poly    = BasisPolynomial(CHEBYSHEV);//, mode);
  Real wt_factor = 1.;
  RealVector pt_factor(numVars, false); pt_factor = 1.;
  size_t pp_cntr = 0;
  for (i=0; i<numVars; ++i) {
    switch (integrationRules[i]) {
    case GAUSS_HERMITE: // Gauss-Hermite open weakly nested
      pt_factor[i] = hermite_poly.point_factor();
      wt_factor   *= hermite_poly.weight_factor();    break;
    case GAUSS_LEGENDRE: case GAUSS_PATTERSON: case GAUSS_PATTERSON_SLOW:
      legendre_poly.gauss_mode(integrationRules[i]);
      wt_factor *= legendre_poly.weight_factor();     break;
    case GAUSS_JACOBI:
      jacobi_poly.beta_stat(polyParams[pp_cntr]+1.);       // convert poly->stat
      jacobi_poly.alpha_stat(polyParams[pp_cntr+1]+1.);    // convert poly->stat
      wt_factor *= jacobi_poly.weight_factor();       break;
    case GEN_GAUSS_LAGUERRE:
      gen_laguerre_poly.alpha_stat(polyParams[pp_cntr]+1.);// convert poly->stat
      wt_factor *= gen_laguerre_poly.weight_factor(); break;
    case CLENSHAW_CURTIS: case CLENSHAW_CURTIS_SLOW:
    case FEJER2:          case FEJER2_SLOW:
      chebyshev_poly.gauss_mode(integrationRules[i]);
      wt_factor *= chebyshev_poly.weight_factor();    break;
    //case GAUSS_LAGUERRE: case GOLUB_WELSCH: // scaling is OK
    }
    pp_cntr += numPolyParams[i];
  }

  // perform scaling
  for (i=0; i<num_colloc_pts; ++i) {
    weightSets[i] *= wt_factor;
    for (j=0; j<numVars; ++j)
      variableSets[i][j] *= pt_factor[j]; // ith column, jth row
  }

  /*
  // -----------------------------------
  // Get sparse grid index/base mappings
  // -----------------------------------
  size_t size = numVars*num_colloc_pts, cntr = 0;
  int* indices = new int [size];
  int* bases   = new int [size];

  webbur::sparse_grid_mixed_growth_index(numVars, ssgLevel,
    integrationRules, num_colloc_pts, num_total_pts, uniqueIndexMapping,
    webbur::level_to_order_default, bases, indices);

  IntArray key(2*numVars);
  unsigned short closed_order_max;
  level_to_order_closed_exponential(ssgLevel, closed_order_max);
  for (i=0; i<num_colloc_pts; i++) {
    for (j=0; j<numVars; j++, cntr++) {
      switch (integrationRules[j]) {
      case GAUSS_HERMITE: case GAUSS_LEGENDRE:
	key[j] = 2 * bases[cntr] + 1;                 // map to quad order
	key[j+numVars] = indices[cntr] + bases[cntr]; // 0-based index
	break;
      case CLENSHAW_CURTIS: case CLENSHAW_CURTIS_SLOW:
	key[j] = closed_order_max;      // promotion to highest grid
	key[j+numVars] = indices[cntr]; // already 0-based
	break;
      case GAUSS_LAGUERRE:
	key[j] = bases[cntr];               // already quad order
	key[j+numVars] = indices[cntr] - 1; // map to 0-based
	break;
      }
    }
    ssgIndexMap[key] = i;
#ifdef DEBUG
    PCout << "i = " << i << " key =\n" << key << std::endl;
#endif // DEBUG
  }
  delete [] indices;
  delete [] bases;
  */
}


void SparseGridDriver::
anisotropic_multi_index(Int2DArray& multi_index, RealArray& coeffs) const
{
  multi_index.clear();
  coeffs.clear();
  // Utilize webbur::sgmga_vcn_{ordered,coef} for 0-based index sets
  // (w*alpha_min-|alpha| < |alpha . j| <= w*alpha_min).
  // With scaling alpha_min = 1: w-|alpha| < |alpha . j| <= w.
  // In the isotropic case, reduces to w-N < |j| <= w, which is the same as
  // w-N+1 <= |j| <= w.
  IntArray x(numVars), x_max(numVars); //x_max = ssgLevel;
  Real wt_sum = 0., q_max = ssgLevel;
  for (size_t i=0; i<numVars; ++i) {
    const Real& wt_i = ssgAnisoLevelWts[i];
    wt_sum += wt_i;
    // minimum nonzero weight is scaled to 1, so just catch special case of 0
    x_max[i] = (std::abs(wt_i) > 1.e-10) ? (int)std::ceil(q_max/wt_i) : 0;
  }
  Real q_min = ssgLevel - wt_sum;
#ifdef DEBUG
  PCout << "q_min = " << q_min << " q_max = " << q_max;
#endif // DEBUG

  bool more = false;
  webbur::sgmga_vcn_ordered(numVars, ssgAnisoLevelWts.values(), &x_max[0],
			    &x[0], q_min, q_max, &more);
  while (more) {
    Real coeff = webbur::sgmga_vcn_coef(numVars, ssgAnisoLevelWts.values(),
					&x_max[0], &x[0], q_min, q_max);
    if (std::abs(coeff) > 1.e-10) {
      multi_index.push_back(x);
      coeffs.push_back(coeff);
    }
    webbur::sgmga_vcn_ordered(numVars, ssgAnisoLevelWts.values(),
			      &x_max[0], &x[0], q_min, q_max, &more);
  }
}


// TO DO: avoid recalculating existing Gauss pts within following functions
//
// Burkardt: instead of params/data, pass index and retain array of NGPoly
//
// Then add interface to pass in array of polynomials in initialize() instead
// of current initialize_grid_parameters().
//
// KDE: look at figTree --> if not too dense; otherwise, simple box kernel +
// Gaussian kernel
//
// add/activate STOCHASTIC_EXPANSION allowing moments or KDE
void SparseGridDriver::
bounded_normal_gauss_points(int order, int num_params, double* params,
			    double* data)
{
  if (num_params != 4) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::bounded_normal_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.bounded_normal_distribution(params[0], params[1], params[2], params[3]);
                                //(mean,      stdev,     lwr,       upr);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::
bounded_normal_gauss_weights(int order, int num_params, double* params,
			     double* data)
{
  if (num_params != 4) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::bounded_normal_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.bounded_normal_distribution(params[0], params[1], params[2], params[3]);
                                //(mean,      stdev,     lwr,       upr);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGridDriver::
lognormal_gauss_points(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::lognormal_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.lognormal_distribution(params[0], params[1]); //(mean, stdev);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::
lognormal_gauss_weights(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::lognormal_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.lognormal_distribution(params[0], params[1]); //(mean, stdev);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGridDriver::
bounded_lognormal_gauss_points(int order, int num_params, double* params,
			       double* data)
{
  if (num_params != 4) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::bounded_lognormal_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.bounded_lognormal_distribution(params[0], params[1], params[2],
				      params[3]); //(mean, stdev, lwr, upr);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::
bounded_lognormal_gauss_weights(int order, int num_params, double* params,
				double* data)
{
  if (num_params != 4) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::bounded_lognormal_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.bounded_lognormal_distribution(params[0], params[1], params[2],
				      params[3]); //(mean, stdev, lwr, upr);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGridDriver::
loguniform_gauss_points(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::loguniform_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.loguniform_distribution(params[0], params[1]); //(lwr, upr);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::
loguniform_gauss_weights(int order, int num_params, double* params,
			 double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::loguniform_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.loguniform_distribution(params[0], params[1]); //(lwr, upr);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGridDriver::
triangular_gauss_points(int order, int num_params, double* params, double* data)
{
  if (num_params != 3) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::triangular_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.triangular_distribution(params[0], params[1], params[2]);
                            //(mode,      lwr,       upr);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::
triangular_gauss_weights(int order, int num_params, double* params,
			 double* data)
{
  if (num_params != 3) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::triangular_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.triangular_distribution(params[0], params[1], params[2]);
                            //(mode,      lwr,       upr);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGridDriver::
gumbel_gauss_points(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::gumbel_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.gumbel_distribution(params[0], params[1]); //(alpha, beta);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::
gumbel_gauss_weights(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::gumbel_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.gumbel_distribution(params[0], params[1]); //(alpha, beta);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGridDriver::
frechet_gauss_points(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::frechet_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.frechet_distribution(params[0], params[1]); //(alpha, beta);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::
frechet_gauss_weights(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::frechet_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.frechet_distribution(params[0], params[1]); //(alpha, beta);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGridDriver::
weibull_gauss_points(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::weibull_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.weibull_distribution(params[0], params[1]); //(alpha, beta);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::
weibull_gauss_weights(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::weibull_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.weibull_distribution(params[0], params[1]); //(alpha, beta);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGridDriver::
histogram_bin_gauss_points(int order, int num_params, double* params,
			   double* data)
{
  if (num_params < 4 || num_params % 2) { // need at least 2 (x,c) pairs
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::histogram_bin_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  RealVector param_vec;
  copy_data(params, num_params, param_vec);
  ngop.histogram_bin_distribution(param_vec);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::
histogram_bin_gauss_weights(int order, int num_params, double* params,
			    double* data)
{
  if (num_params < 4 || num_params % 2) { // need at least 2 (x,c) pairs
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::histogram_bin_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  RealVector param_vec;
  copy_data(params, num_params, param_vec);
  ngop.histogram_bin_distribution(param_vec);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}

} // namespace Pecos

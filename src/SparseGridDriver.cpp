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
#include "sandia_sgmg.H"
#include "sandia_sgmga.H"
#include "DistributionParams.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: SparseGridDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {

/// initialize static member pointer to active driver instance
SparseGridDriver* SparseGridDriver::sgdInstance(NULL);


void SparseGridDriver::dimension_preference(const RealVector& dim_pref)
{
  RealVector aniso_wts;
  if (!dim_pref.empty()) {
    size_t num_pref = dim_pref.length();
    aniso_wts.sizeUninitialized(num_pref);
    webbur::sandia_sgmga_importance_to_aniso(num_pref, dim_pref.values(),
					     aniso_wts.values());
#ifdef DEBUG
    PCout << "dimension preference:\n"; write_data(PCout, dim_pref);
    PCout << "anisotropic weights after sandia_sgmga_importance_to_aniso():\n";
    write_data(PCout, aniso_wts);
#endif
  }
  anisotropic_weights(aniso_wts);
}


void SparseGridDriver::anisotropic_weights(const RealVector& aniso_wts)
{
  if (aniso_wts.empty())
    dimIsotropic = true;
  else {
    if (aniso_wts.length() != numVars) {
      PCerr << "Error: length of sparse grid anisotropic weights "
	    << "specification is inconsistent with\n       number of variables "
	    << "in SparseGridDriver::anisotropic_weights()." << std::endl;
      abort_handler(-1);
    }

    dimIsotropic = true;
    anisoLevelWts.resize(numVars);
    // truncate any negative values
    size_t i;
    for (i=0; i<numVars; ++i)
      anisoLevelWts[i] = (aniso_wts[i] < 0.) ? 0. : aniso_wts[i];
    // detect anisotropy
    Real wt0 = anisoLevelWts[0];
    for (i=1; i<numVars; ++i)
      if (std::abs(anisoLevelWts[i] - wt0) > DBL_EPSILON)
	{ dimIsotropic = false; break; }
    // normalize and enforce axis lower bounds/weight upper bounds
    if (!dimIsotropic) {
      int option = 1; // weights scaled so that minimum nonzero entry is 1
      webbur::sandia_sgmga_aniso_normalize(option, numVars,
					   anisoLevelWts.values());
#ifdef DEBUG
      PCout << "anisoLevelWts after sandia_sgmga_aniso_normalize():\n";
      write_data(PCout, anisoLevelWts);
#endif
      // enforce axis lower bounds, if present, for current ssgLevel.  An axis
      // lower bound defines a weight upper bound based on the current ssgLevel:
      // LB_i = level*wt_min/wt_i --> wt_i = level*wt_min/LB_i and wt_min=1.
      // Catch special case of dim_pref_i = 0 --> wt_i = LB_i = 0.
      if (!axisLowerBounds.empty()) {
	for (i=0; i<numVars; ++i)
	  if (axisLowerBounds[i] > 1.e-10) {                 // nonzero LB
	    Real wt_u_bnd = (Real)ssgLevel/axisLowerBounds[i];
	    anisoLevelWts[i] = (anisoLevelWts[i] > 1.e-10) ? // nonzero wt
	      std::min(wt_u_bnd, anisoLevelWts[i]) : wt_u_bnd;
	  }
#ifdef DEBUG
	PCout << "anisoLevelWts after axisLowerBounds enforcement:\n";
	write_data(PCout, anisoLevelWts);
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
  if (dimIsotropic)
    axisLowerBounds = (Real)ssgLevel; // all weights = 1
  else // min nonzero weight scaled to 1 --> just catch special case w_i=0
    for (size_t i=0; i<numVars; ++i)
      axisLowerBounds[i] = (anisoLevelWts[i] > 1.e-10) ? // nonzero wt
	(Real)ssgLevel/anisoLevelWts[i] : 0.;
}


void SparseGridDriver::
initialize_grid(const ShortArray& u_types,  size_t ssg_level,
		const RealVector& dim_pref, bool store_1d_gauss,
		short growth_rate, short nested_uniform_rule)
{
  numVars = u_types.size();
  compute1DPoints.resize(numVars);
  compute1DWeights.resize(numVars);

  store1DGauss = store_1d_gauss;
  level(ssg_level);
  dimension_preference(dim_pref);

  // For STANDARD exponential growth, use of nested rules is restricted to
  // isotropic uniform in order to enforce consistent growth rates:
  bool nested_rules = true;
  if (growth_rate == UNRESTRICTED_GROWTH)
    for (size_t i=0; i<numVars; ++i)
      if (u_types[i] != STD_UNIFORM)
	{ nested_rules = false; break; }
  // For MODERATE exponential growth, nested rules can be used heterogeneously
  // and synchronized with STANDARD Gaussian linear growth.
  // For SLOW exponential growth, nested rules can be used heterogeneously
  // and synchronized with SLOW Gaussian linear growth (not yet available).

  bool cheby_poly = false;
  for (size_t i=0; i<numVars; i++) {
    // set compute1DPoints/compute1DWeights
    if (u_types[i] == STD_UNIFORM && nested_rules && 
	nested_uniform_rule != GAUSS_PATTERSON) {
      compute1DPoints[i]  = chebyshev_points;
      compute1DWeights[i] = chebyshev_weights;
      cheby_poly = true;
    }
    else {
      compute1DPoints[i]  = basis_gauss_points;
      compute1DWeights[i] = basis_gauss_weights;
    }
  }
  if (cheby_poly) // gauss_mode set within loops
    chebyPolyPtr = new BasisPolynomial(CHEBYSHEV);

  initialize_rules(u_types, nested_rules, growth_rate, nested_uniform_rule,
		   integrationRules, growthRules);
}


void SparseGridDriver::
initialize_grid(const std::vector<BasisPolynomial>& poly_basis,
		size_t ssg_level, const RealVector& dim_pref,
		bool store_1d_gauss, short growth_rate)
{
  numVars         = poly_basis.size();
  polynomialBasis = poly_basis; // shallow copy
  store1DGauss    = store_1d_gauss;

  level(ssg_level);
  dimension_preference(dim_pref);

  compute1DPoints.resize(numVars);
  compute1DWeights.resize(numVars);
  for (size_t i=0; i<numVars; i++) {
    compute1DPoints[i]  = basis_gauss_points;
    compute1DWeights[i] = basis_gauss_weights;
  }

  initialize_rules(poly_basis, growth_rate, integrationRules, growthRules);
}


int SparseGridDriver::grid_size()
{
  // do this here (called at beginning of compute_grid()) since sgdInstance
  // required within compute1DPoints below
  sgdInstance = this;

  return (dimIsotropic) ?
    webbur::sgmg_size(numVars, ssgLevel, &integrationRules[0],
      &compute1DPoints[0], duplicateTol, &growthRules[0]) :
    webbur::sandia_sgmga_size(numVars, anisoLevelWts.values(), ssgLevel,
      &integrationRules[0], &compute1DPoints[0], duplicateTol, &growthRules[0]);
}


int SparseGridDriver::grid_size_total()
{
  return (dimIsotropic) ?
    webbur::sgmg_size_total(numVars, ssgLevel, &integrationRules[0],
      &growthRules[0]) :
    webbur::sandia_sgmga_size_total(numVars, anisoLevelWts.values(), ssgLevel,
      &integrationRules[0], &growthRules[0]);
}


void SparseGridDriver::compute_grid()
{
  // --------------------------------
  // Get number of collocation points
  // --------------------------------
  // Note:  grid_size() sets sgdInstance
  // TO DO: order calls for data reuse
  int num_colloc_pts = grid_size(), num_total_pts = grid_size_total();
#ifdef DEBUG
  PCout << "Total number of sparse grid integration points: "
	<< num_colloc_pts << '\n';
#endif // DEBUG

  // ----------------------------------------------
  // Get collocation points and integration weights
  // ----------------------------------------------
  weightSets.sizeUninitialized(num_colloc_pts);
  variableSets.shapeUninitialized(numVars, num_colloc_pts);// Teuchos: col major
  uniqueIndexMapping.resize(num_total_pts);

  int* sparse_order = new int [num_colloc_pts*numVars];
  int* sparse_index = new int [num_colloc_pts*numVars];
  if (dimIsotropic) {
    webbur::sgmg_unique_index(numVars, ssgLevel, &integrationRules[0],
      &compute1DPoints[0], duplicateTol, num_colloc_pts, num_total_pts,
      &growthRules[0], &uniqueIndexMapping[0]);
    webbur::sgmg_index(numVars, ssgLevel, &integrationRules[0], num_colloc_pts,
      num_total_pts, &uniqueIndexMapping[0], &growthRules[0], sparse_order,
      sparse_index);
    webbur::sgmg_weight(numVars, ssgLevel, &integrationRules[0],
      &compute1DWeights[0], num_colloc_pts, num_total_pts,
      &uniqueIndexMapping[0], &growthRules[0], weightSets.values());
    webbur::sgmg_point(numVars, ssgLevel, &integrationRules[0],
      &compute1DPoints[0], num_colloc_pts, sparse_order, sparse_index,
      &growthRules[0], variableSets.values());
  }
  else {
    webbur::sandia_sgmga_unique_index(numVars, anisoLevelWts.values(), ssgLevel,
      &integrationRules[0], &compute1DPoints[0], duplicateTol, num_colloc_pts,
      num_total_pts, &growthRules[0], &uniqueIndexMapping[0]);
    webbur::sandia_sgmga_index(numVars, anisoLevelWts.values(), ssgLevel,
      &integrationRules[0], num_colloc_pts, num_total_pts,
      &uniqueIndexMapping[0], &growthRules[0], sparse_order, sparse_index);
    webbur::sandia_sgmga_weight(numVars, anisoLevelWts.values(), ssgLevel,
      &integrationRules[0], &compute1DWeights[0], num_colloc_pts, num_total_pts,
      &uniqueIndexMapping[0], &growthRules[0], weightSets.values());
    webbur::sandia_sgmga_point(numVars, anisoLevelWts.values(), ssgLevel,
      &integrationRules[0], &compute1DPoints[0], num_colloc_pts, sparse_order,
      sparse_index, &growthRules[0], variableSets.values());
  }
  delete [] sparse_order;
  delete [] sparse_index;
#ifdef DEBUG
  PCout << "uniqueIndexMapping:\n" << uniqueIndexMapping << '\n';
#endif

  if (store1DGauss) { // 1D arrays not needed for PCE
    // ----------------------------
    // Define 1-D point/weight sets
    // ----------------------------
    if (gaussPts1D.empty())
      gaussPts1D.resize(numVars);
    if (gaussWts1D.empty())
      gaussWts1D.resize(numVars);
    // level_index (j indexing) range is 0:w, level (i indexing) range is 1:w+1
    unsigned short level_index, order;
    for (size_t i=0; i<numVars; i++) {
      gaussPts1D[i].resize(ssgLevel + 1); gaussWts1D[i].resize(ssgLevel + 1);
      switch (integrationRules[i]) {
      case CLENSHAW_CURTIS: case FEJER2:
	chebyPolyPtr->gauss_mode(integrationRules[i]); // integration mode
	for (level_index=0; level_index<=ssgLevel; level_index++) {
	  level_to_order(i, level_index, order);
	  gaussPts1D[i][level_index] = chebyPolyPtr->gauss_points(order);
	  gaussWts1D[i][level_index] = chebyPolyPtr->gauss_weights(order);
	}
	break;
      default: // Gaussian rules
	for (level_index=0; level_index<=ssgLevel; level_index++) {
	  level_to_order(i, level_index, order);
	  gaussPts1D[i][level_index] = polynomialBasis[i].gauss_points(order);
	  gaussWts1D[i][level_index] = polynomialBasis[i].gauss_weights(order);
	}
	break;
      }
    }
  }

  /*
  // -----------------------------------
  // Get sparse grid index/base mappings
  // -----------------------------------
  size_t size = numVars*num_colloc_pts, cntr = 0;
  int* indices = new int [size];
  int* bases   = new int [size];

  webbur::sgmg_index(numVars, ssgLevel, integrationRules, num_colloc_pts,
    num_total_pts, uniqueIndexMapping, &growthRules[0], bases, indices);

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
      case CLENSHAW_CURTIS:
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
  // Utilize webbur::sandia_sgmga_vcn_{ordered,coef} for 0-based index sets
  // (w*alpha_min-|alpha| < |alpha . j| <= w*alpha_min).
  // With scaling alpha_min = 1: w-|alpha| < |alpha . j| <= w.
  // In the isotropic case, reduces to w-N < |j| <= w, which is the same as
  // w-N+1 <= |j| <= w.
  IntArray x(numVars), x_max(numVars); //x_max = ssgLevel;
  Real wt_sum = 0., q_max = ssgLevel;
  for (size_t i=0; i<numVars; ++i) {
    const Real& wt_i = anisoLevelWts[i];
    wt_sum += wt_i;
    // minimum nonzero weight is scaled to 1, so just catch special case of 0
    x_max[i] = (wt_i > 1.e-10) ? (int)std::ceil(q_max/wt_i) : 0;
  }
  Real q_min = ssgLevel - wt_sum;
#ifdef DEBUG
  PCout << "q_min = " << q_min << " q_max = " << q_max;
#endif // DEBUG

  bool more = false;
  webbur::sandia_sgmga_vcn_ordered(numVars, anisoLevelWts.values(), &x_max[0],
				   &x[0], q_min, q_max, &more);
  while (more) {
    Real coeff = webbur::sandia_sgmga_vcn_coef(numVars, anisoLevelWts.values(),
					       &x_max[0], &x[0], q_min, q_max);
    if (std::abs(coeff) > 1.e-10) {
      multi_index.push_back(x);
      coeffs.push_back(coeff);
    }
    webbur::sandia_sgmga_vcn_ordered(numVars, anisoLevelWts.values(), &x_max[0],
				     &x[0], q_min, q_max, &more);
  }
}


// TO DO:
// KDE: look at figTree --> if lightweight; else, code simple box kernel and
// Gaussian kernel
//
// sampling on stoch expansion within \xi --> kernel added for each sample -->
// approximate PDF for response --> accept a KDE PDF within NumGenOrthogPoly
//
// add/activate STOCHASTIC_EXPANSION allowing moments or KDE
void SparseGridDriver::basis_gauss_points(int order, int index, double* data)
{
  const RealArray& gauss_pts
    = sgdInstance->polynomialBasis[index].gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::basis_gauss_weights(int order, int index, double* data)
{
  const RealArray& gauss_wts
    = sgdInstance->polynomialBasis[index].gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGridDriver::chebyshev_points(int order, int index, double* data)
{
  sgdInstance->chebyPolyPtr->gauss_mode(sgdInstance->integrationRules[index]);
  const RealArray& gauss_pts = sgdInstance->chebyPolyPtr->gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::chebyshev_weights(int order, int index, double* data)
{
  sgdInstance->chebyPolyPtr->gauss_mode(sgdInstance->integrationRules[index]);
  const RealArray& gauss_wts = sgdInstance->chebyPolyPtr->gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}

} // namespace Pecos

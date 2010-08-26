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
#include "PolynomialApproximation.hpp"
#include "sandia_rules.H"
#include "sandia_sgmg.H"
#include "sandia_sgmga.H"
#include "sandia_sgmgg.H"
#include "DistributionParams.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: SparseGridDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {

/// initialize static member pointer to active driver instance
SparseGridDriver* SparseGridDriver::sgdInstance(NULL);


void SparseGridDriver::
allocate_smolyak_arrays(UShort2DArray& multi_index, IntArray& coeffs)
{
  // Populate smolyakMultiIndex and smolyakCoeffs.  Identifies
  // use of polynomialBasis[variable][index] based on index 0:num_levels-1.
  // w = q - N = dimension-independent level.  For isotropic,
  //   w + 1 <= |i| <= w + N for i starts at 1 (used for index set defn.)
  //   w - N + 1 <= |j| <= w for j = i - 1 starts at 0 (used for generation)
  // For anisotropic, a weighted linear index set constraint is used.

  size_t i;
  if (dimIsotropic) { // initialize multi_index
    UShortArray levels(numVars, ssgLevel);
    PolynomialApproximation::total_order_multi_index(levels, multi_index,
						     numVars-1);
    size_t num_terms = multi_index.size();
    // initialize coeffs
    coeffs.resize(num_terms);
    for (i=0; i<num_terms; i++) {
      int wpNmi = ssgLevel - index_norm(multi_index[i]); // w+N-|i| = w-|j|
      coeffs[i] = (int)std::pow(-1., wpNmi)
	* (int)BasisPolynomial::n_choose_k(numVars - 1, wpNmi);
    }
  }
  else { // utilize Pecos wrapper to sgmga_vcn_{ordered,coef}
    multi_index.clear();
    coeffs.clear();
    // Utilize webbur::sandia_sgmga_vcn_{ordered,coef} for 0-based index sets
    // (w*alpha_min-|alpha| < |alpha . j| <= w*alpha_min).
    // With scaling alpha_min = 1: w-|alpha| < |alpha . j| <= w.
    // In the isotropic case, reduces to w-N < |j| <= w, which is the same as
    // w-N+1 <= |j| <= w.
    IntArray x(numVars), x_max(numVars); //x_max = ssgLevel;
    UShortArray index_set(numVars);
    Real wt_sum = 0., q_max = ssgLevel;
    for (i=0; i<numVars; ++i) {
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
    Real *aniso_wts = anisoLevelWts.values();
    int  *x0 = &x[0], *xm0 = &x_max[0], coeff;
    webbur::sandia_sgmga_vcn_ordered(numVars, aniso_wts, xm0, x0,
				     q_min, q_max, &more);
    while (more) {
      coeff = (int)webbur::sandia_sgmga_vcn_coef(numVars, aniso_wts, x0, q_max);
      if (coeff) {
	coeffs.push_back(coeff);
	for (i=0; i<numVars; ++i)
	  index_set[i] = (unsigned short)x[i];
	multi_index.push_back(index_set);
      }
      webbur::sandia_sgmga_vcn_ordered(numVars, aniso_wts, xm0, x0,
				       q_min, q_max, &more);
    }
  }

#ifdef DEBUG
  size_t num_terms = coeffs.size();
  PCout << "\nnum Smolyak terms = " << num_terms << '\n';
  for (i=0; i<num_terms; i++)
    PCout << "multi_index[" << i << "]:\n" << multi_index[i]
	  << "coeffs[" << i << "] = " << coeffs[i] << "\n\n";
#endif // DEBUG
}


void SparseGridDriver::
allocate_generalized_coefficients(const UShort2DArray& multi_index,
				  IntArray& coeffs) const
{
  size_t i, j, cntr = 0, num_sets = multi_index.size();
  if (coeffs.size() != num_sets)
    coeffs.resize(num_sets);
  int* mi = new int [numVars*num_sets];
  //copy_data(multi_index, mi); // UShort2DArray -> int*
  for (i=0; i<num_sets; ++i)
    for (j=0; j<numVars; ++j, ++cntr)
      mi[cntr] = multi_index[i][j]; // sgmgg packs by variable groups
  webbur2::sgmgg_coef_naive(numVars, num_sets, mi, &coeffs[0]);
  delete [] mi;
}


void SparseGridDriver::allocate_collocation_arrays()
{
  // define mapping from 1:numCollocPts to set of 1d interpolation indices
  size_t num_smolyak_indices = smolyakMultiIndex.size();
  collocKey.resize(num_smolyak_indices);
  expansionCoeffIndices.resize(num_smolyak_indices);
  UShortArray quad_order(numVars); //, gauss_indices(numVars);
  //IntArray key(2*numVars);
  //unsigned short closed_order_max;
  //level_to_order_closed_exponential(ssg_level, closed_order_max);
  size_t i, j, cntr = 0;
  for (i=0; i<num_smolyak_indices; ++i) {
    UShort2DArray& key_i       = collocKey[i];
    SizetArray&    coeff_map_i = expansionCoeffIndices[i];

    level_to_order(smolyakMultiIndex[i], quad_order);
    PolynomialApproximation::tensor_product_multi_index(quad_order, key_i);
    size_t num_tp_pts = key_i.size();
    coeff_map_i.resize(num_tp_pts);
    //gauss_indices = 0;
    for (j=0; j<num_tp_pts; ++j, ++cntr) {
      coeff_map_i[j] = uniqueIndexMapping[cntr];

      /*
      // update expansionCoeffIndices: sparse grid nesting of points may be
      // tracked by demoting each point to its lowest order representation
      // or by promoting each point to the highest order grid.
      for (k=0; k<numVars; k++) {
	switch (integrationRules[k]) {
	// For open weakly-nested rules (Gauss-Hermite, Gauss-Legendre),
	// demote to order=1,index=0 if a center pt for a particular order>1
	case GAUSS_HERMITE: case GAUSS_LEGENDRE:
	  if (numVars > 1 && quad_order[k] > 1 &&
	      gauss_indices[k] == (quad_order[k]-1)/2) // demoted base/index
	    { key[k] = 1;        key[k+numVars] = 0; }
	  else                                    // original base/index
	    { key[k] = quad_order[k]; key[k+numVars] = gauss_indices[k]; }
	  break;
	// For closed nested rules (Clenshaw-Curtis), base is a dummy and
	// index is promoted to the highest order grid
	case CLENSHAW_CURTIS:
	  key[k] = closed_order_max; // promoted base
	  if (sm_index_i[k] == 0)
	    key[k+numVars] = (closed_order_max-1)/2;      // promoted index
	  else {
	    key[k+numVars] = gauss_indices[k];
	    unsigned short delta_w = ssg_level - sm_index_i[k];
	    if (delta_w)
	      key[k+numVars] *= (size_t)pow(2., delta_w); // promoted index
	  }
	  break;
	// For open non-nested rules (Gauss-Laguerre), no modification reqd
        // Note: might need to check for symmetric cases of Gauss-Jacobi
	default: // original base/index
	  key[k] = quad_order[k]; key[k+numVars] = gauss_indices[k];
	  break;
	}
      }
#ifdef DEBUG
      PCout << "lookup key:\n" << key << std::endl;
#endif // DEBUG
      IntArraySizetMap::const_iterator cit = ssgIndexMap.find(key);
      if (cit == ssgIndexMap.end()) {
	PCerr << "Error: lookup on sparse grid index map failed in "
	      << "InterpPolyApproximation::find_coefficients()"
	      << std::endl;
	abort_handler(-1);
      }
      else
	coeff_map_i[j] = cit->second;

      // increment the n-dimensional gauss point index set
      if (j != num_tp_pts - 1)
        increment_indices(gauss_indices, quad_order, true);
      */
#ifdef DEBUG
      PCout << "collocKey[" << i << "][" << j << "]:\n" << key_i[j]
	    << "expansionCoeffIndices[" << i << "][" << j << "] = "
	    << coeff_map_i[j] << '\n';
#endif // DEBUG
    }
  }
}


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
initialize_grid(const ShortArray& u_types,  unsigned short ssg_level,
		const RealVector& dim_pref, bool store_1d_gauss,
		bool nested_rules, short growth_rate, short nested_uniform_rule)
{
  numVars = u_types.size();
  compute1DPoints.resize(numVars);
  compute1DWeights.resize(numVars);

  store1DGauss = store_1d_gauss;
  level(ssg_level);
  dimension_preference(dim_pref);

  // For unrestricted exponential growth, use of nested rules is restricted
  // to uniform/normal in order to enforce similar growth rates:
  if (nested_rules && growth_rate == UNRESTRICTED_GROWTH)
    for (size_t i=0; i<numVars; ++i)
      if (u_types[i] != STD_UNIFORM && u_types[i] != STD_NORMAL)
	{ nested_rules = false; break; }
  // For MODERATE and SLOW restricted exponential growth, nested rules can be
  // used heterogeneously and synchronized with STANDARD and SLOW Gaussian
  // linear growth, respectively.

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
  if (cheby_poly && !chebyPolyPtr) // gauss_mode set within loops
    chebyPolyPtr = new BasisPolynomial(CHEBYSHEV);

  initialize_rules(u_types, nested_rules, growth_rate, nested_uniform_rule,
		   integrationRules, growthRules);
}


void SparseGridDriver::
initialize_grid(const std::vector<BasisPolynomial>& poly_basis,
		unsigned short ssg_level, const RealVector& dim_pref,
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
  // TO DO: improve efficiency of these calls through data reuse
  numCollocPts = grid_size();
  int num_total_pts = grid_size_total();
#ifdef DEBUG
  PCout << "Total number of sparse grid integration points: "
	<< numCollocPts << '\n';
#endif // DEBUG

  // ----------------------------------------------
  // Get collocation points and integration weights
  // ----------------------------------------------
  weightSets.sizeUninitialized(numCollocPts);
  variableSets.shapeUninitialized(numVars, numCollocPts);// Teuchos: col major
  uniqueIndexMapping.resize(num_total_pts);

  int* sparse_order = new int [numCollocPts*numVars];
  int* sparse_index = new int [numCollocPts*numVars];
  if (dimIsotropic) {
    webbur::sgmg_unique_index(numVars, ssgLevel, &integrationRules[0],
      &compute1DPoints[0], duplicateTol, numCollocPts, num_total_pts,
      &growthRules[0], &uniqueIndexMapping[0]);
    webbur::sgmg_index(numVars, ssgLevel, &integrationRules[0], numCollocPts,
      num_total_pts, &uniqueIndexMapping[0], &growthRules[0], sparse_order,
      sparse_index);
    webbur::sgmg_weight(numVars, ssgLevel, &integrationRules[0],
      &compute1DWeights[0], numCollocPts, num_total_pts,
      &uniqueIndexMapping[0], &growthRules[0], weightSets.values());
    webbur::sgmg_point(numVars, ssgLevel, &integrationRules[0],
      &compute1DPoints[0], numCollocPts, sparse_order, sparse_index,
      &growthRules[0], variableSets.values());
  }
  else {
    webbur::sandia_sgmga_unique_index(numVars, anisoLevelWts.values(), ssgLevel,
      &integrationRules[0], &compute1DPoints[0], duplicateTol, numCollocPts,
      num_total_pts, &growthRules[0], &uniqueIndexMapping[0]);
    webbur::sandia_sgmga_index(numVars, anisoLevelWts.values(), ssgLevel,
      &integrationRules[0], numCollocPts, num_total_pts,
      &uniqueIndexMapping[0], &growthRules[0], sparse_order, sparse_index);
    webbur::sandia_sgmga_weight(numVars, anisoLevelWts.values(), ssgLevel,
      &integrationRules[0], &compute1DWeights[0], numCollocPts, num_total_pts,
      &uniqueIndexMapping[0], &growthRules[0], weightSets.values());
    webbur::sandia_sgmga_point(numVars, anisoLevelWts.values(), ssgLevel,
      &integrationRules[0], &compute1DPoints[0], numCollocPts, sparse_order,
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
    size_t i, num_levels = ssgLevel + 1;
    if (gaussPts1D.size() != num_levels || gaussWts1D.size() != num_levels) {
      gaussPts1D.resize(num_levels); gaussWts1D.resize(num_levels);
      for (i=0; i<num_levels; ++i)
	{ gaussPts1D[i].resize(numVars); gaussWts1D[i].resize(numVars); }
    }
    // level_index (j indexing) range is 0:w, level (i indexing) range is 1:w+1
    unsigned short level_index, order;
    for (i=0; i<numVars; i++) {
      switch (integrationRules[i]) {
      case CLENSHAW_CURTIS: case FEJER2:
	chebyPolyPtr->gauss_mode(integrationRules[i]); // integration mode
	for (level_index=0; level_index<num_levels; ++level_index) {
	  level_to_order(i, level_index, order);
	  gaussPts1D[level_index][i] = chebyPolyPtr->gauss_points(order);
	  gaussWts1D[level_index][i] = chebyPolyPtr->gauss_weights(order);
	}
	break;
      default: // Gaussian rules
	for (level_index=0; level_index<num_levels; ++level_index) {
	  level_to_order(i, level_index, order);
	  gaussPts1D[level_index][i] = polynomialBasis[i].gauss_points(order);
	  gaussWts1D[level_index][i] = polynomialBasis[i].gauss_weights(order);
	}
	break;
      }
    }
  }

  /*
  // -----------------------------------
  // Get sparse grid index/base mappings
  // -----------------------------------
  size_t size = numVars*numCollocPts, cntr = 0;
  int* indices = new int [size];
  int* bases   = new int [size];

  webbur::sgmg_index(numVars, ssgLevel, integrationRules, numCollocPts,
    num_total_pts, uniqueIndexMapping, &growthRules[0], bases, indices);

  IntArray key(2*numVars);
  unsigned short closed_order_max;
  level_to_order_closed_exponential(ssgLevel, closed_order_max);
  for (i=0; i<numCollocPts; i++) {
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


void SparseGridDriver::initialize_sets()
{
  // define set O (old) from smolyakMultiIndex and smolyakCoeffs:
  //oldMultiIndex = smolyakMultiIndex;
  oldMultiIndex.clear();
  oldMultiIndex.insert(smolyakMultiIndex.begin(), smolyakMultiIndex.end());

  // compute initial set A (active) by applying add_active_neighbors()
  // to the frontier of smolyakMultiIndex:
  size_t i, num_old_sets = smolyakCoeffs.size();
  // anisotropic test on coeff==1 is necessary but not sufficient for presence
  // on index set frontier, requiring an additional logic test within
  // add_active_neighbors().  This is currently the best we can do since the
  // weighted norm of the index set may differ from the level.
  for (i=0; i<num_old_sets; ++i)
    if ( smolyakCoeffs[i] == 1 && ( !dimIsotropic || // imperfect for aniso
	 ( dimIsotropic && index_norm(smolyakMultiIndex[i]) == ssgLevel ) ) )
	add_active_neighbors(smolyakMultiIndex[i]);

#ifdef DEBUG
  PCout << "SparseGridDriver::initialize_sets():\nold sets:\n" << oldMultiIndex
	<< "active sets:\n" << activeMultiIndex << std::endl;
#endif // DEBUG
}


void SparseGridDriver::push_trial_set(const UShortArray& set)
{
  trialIndexSet = set;
  smolyakMultiIndex.push_back(set);
  // update smolyakCoeffs from smolyakMultiIndex
  allocate_generalized_coefficients(smolyakMultiIndex, smolyakCoeffs);
}


void SparseGridDriver::pop_trial_set()
{
  UShort3DArray::iterator it = --collocKey.end();
  numCollocPts -= it->size(); // subtract number of trial points
  uniqueIndexMapping.resize(numCollocPts); // prune trial set from end

  //trialIndexSet.clear();
  smolyakMultiIndex.pop_back();
  // no need to update smolyakCoeffs as this will be updated on next push
  collocKey.pop_back();
  expansionCoeffIndices.pop_back();
}


void SparseGridDriver::compute_trial_grid()
{
  // compute variableSets/weightSets and update collocKey
  UShortArray quad_order(numVars);
  UShort2DArray new_key; collocKey.push_back(new_key);
  UShort3DArray::iterator it = --collocKey.end();
  Real2DArray pts_1d(numVars), wts_1d(numVars);
  level_to_order(trialIndexSet, quad_order);
  compute_tensor_grid(quad_order, *it, pts_1d, wts_1d);

  // if needed, update 3D with new 2D gauss pts/wts (in correct location)
  size_t i, num_levels = gaussPts1D.size(), max_level = 0;
  for (i=0; i<numVars; ++i)
    if (trialIndexSet[i] > max_level)
      max_level = trialIndexSet[i];
  if (max_level >= num_levels) {
    gaussPts1D.resize(max_level+1); gaussWts1D.resize(max_level+1);
    for (i=num_levels; i<=max_level; ++i)
      { gaussPts1D[i].resize(numVars); gaussWts1D[i].resize(numVars); }
    for (i=0; i<numVars; ++i) {
      unsigned short trial_index = trialIndexSet[i];
      if (trial_index >= num_levels) {
	gaussPts1D[trial_index][i] = pts_1d[i];
	gaussWts1D[trial_index][i] = wts_1d[i];
      }
    }
  }

  // update expansionCoeffIndices, uniqueIndexMapping
  //update_collocation_arrays();
  //
  // prior to having updated uniqueIndexMapping available from sandia_sgmgg,
  // manually increment the point sets (ignoring the merging of duplicates).
  // The right aggregated weights are picked up for PCE, but SC requires
  // consolidated pts/wts in moment calculcations!!
  size_t index, num_tp_pts = it->size();
  SizetArray new_coeff_map(num_tp_pts);
  for (i=0; i<num_tp_pts; ++i) {
    // can't search collocKey for duplicates since indices are relative to order
    index = numCollocPts + i; // temp hack, prior to sandia_sgmgg management
    uniqueIndexMapping.push_back(index);
    new_coeff_map[i] = index;
  }
  expansionCoeffIndices.push_back(new_coeff_map);
  numCollocPts += num_tp_pts; // for now, prior to sgmgg duplicate management
}


void SparseGridDriver::finalize_sets()
{
  // for final answer, push all evaluated sets into old and clear active
  //smolyakMultiIndex = oldMultiIndex; // not needed if all trials are popped
  smolyakMultiIndex.insert(smolyakMultiIndex.end(), activeMultiIndex.begin(),
			   activeMultiIndex.end());
  activeMultiIndex.clear();
}


void SparseGridDriver::update_sets(const UShortArray& set_star)
{
  // update set O by adding set_star to oldMultiIndex:
  oldMultiIndex.insert(set_star);

  // remove set_star from set A by erasing from activeMultiIndex: 
  activeMultiIndex.erase(set_star);
  // update activeMultiIndex (set A) based on neighbors of set_star: 
  add_active_neighbors(set_star);

  // update evaluation set smolyakMultiIndex:
  // this push is permanent (will not be popped)
  push_trial_set(set_star);

  // TO DO: prune irrelevant sets that have Coeff = 0
}


void SparseGridDriver::add_active_neighbors(const UShortArray& set)
{
  UShortArray trial_set = set;
  std::set<UShortArray>::const_iterator cit;
  size_t i, j;
  for (i=0; i<numVars; ++i) {
    // i^{th} candidate for set A (active) computed from forward neighbor:
    // increment by 1 in dimension i
    unsigned short& trial_set_i = trial_set[i];
    trial_set_i += 1;
    // anisotropic initialize_sets() candidates could be in oldMultiIndex
    // since smolyakCoeffs[i]==1 test is necessary but not sufficient
    if (dimIsotropic || oldMultiIndex.find(trial_set) == oldMultiIndex.end()) {
      // test all backwards neighbors for membership in set O (old)
      bool backward_old = true;
      for (j=0; j<numVars; ++j) {
	unsigned short& trial_set_j = trial_set[j];
	if (trial_set_j) { // if 0, then admissible by default
	  trial_set_j -= 1;
	  cit = oldMultiIndex.find(trial_set);
	  trial_set_j += 1; // restore
	  if (cit == oldMultiIndex.end())
	    { backward_old = false; break; }
	}
      }
      if (backward_old) // std::set<> will discard any active duplicates
	activeMultiIndex.insert(trial_set);
    }
    trial_set_i -= 1; // restore
  }
}

} // namespace Pecos

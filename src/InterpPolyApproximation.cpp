/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        InterpPolyApproximation
//- Description:  Implementation code for InterpPolyApproximation class
//-               
//- Owner:        Mike Eldred

#include "InterpPolyApproximation.hpp"
#include "TensorProductDriver.hpp"
#include "SparseGridDriver.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

//#define DEBUG
//#define INTERPOLATION_TEST


namespace Pecos {

int InterpPolyApproximation::min_coefficients() const
{
  // return the minimum number of data instances required to build the 
  // surface of multiple dimensions
  return (configOptions.expansionCoeffFlag ||
	  configOptions.expansionCoeffGradFlag) ? 1 : 0;
}


void InterpPolyApproximation::
distribution_types(short& poly_type_1d, short& rule, bool& push_pts)
{
  switch (basisType) {
  case PIECEWISE_INTERPOLATION_POLYNOMIAL:
    poly_type_1d = (configOptions.useDerivs) ?
      PIECEWISE_CUBIC_INTERP : PIECEWISE_LINEAR_INTERP;
    rule = NEWTON_COTES; push_pts = false;                   break;
  case GLOBAL_INTERPOLATION_POLYNOMIAL:
    poly_type_1d = (configOptions.useDerivs) ? HERMITE_INTERP : LAGRANGE_INTERP;
    rule = NO_RULE; push_pts = true;                         break;
  default:
    poly_type_1d = NO_POLY; rule = NO_RULE; push_pts = true; break;
  }
}


void InterpPolyApproximation::allocate_arrays()
{
  allocate_component_effects();
  allocate_total_effects();

  if (configOptions.expansionCoeffFlag &&
      expansionCoeffs.length() != numCollocPts)
    expansionCoeffs.sizeUninitialized(numCollocPts);
  if (configOptions.expansionCoeffGradFlag) {
    const SurrogateDataPoint& sdp
      = (anchorPoint.is_null()) ? dataPoints[0] : anchorPoint;
    size_t num_deriv_vars = sdp.response_gradient().length();
    if (expansionCoeffGrads.numRows() != num_deriv_vars ||
	expansionCoeffGrads.numCols() != numCollocPts)
      expansionCoeffGrads.shapeUninitialized(num_deriv_vars, numCollocPts);
  }

  // checking numCollocPts is insufficient due to anisotropy --> changes in
  // anisotropic weights could move points around without changing the total.
  //bool update_exp_form =
  //  ( (configOptions.expansionCoeffFlag &&
  //     expansionCoeffs.length()      != numCollocPts) ||
  //    (configOptions.expansionCoeffGradFlag &&
  //     expansionCoeffGrads.numCols() != numCollocPts ) );

  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    const UShortArray&   quad_order = tpq_driver->quadrature_order();

    // verify total number of collocation pts (should not include anchorPoint)
    size_t i, j, num_colloc_pts = 1;
    for (i=0; i<numVars; ++i)
      num_colloc_pts *= quad_order[i];
    if (num_colloc_pts != numCollocPts) {
      PCerr << "Error: inconsistent total collocation point count in "
	    << "InterpPolyApproximation::allocate_arrays()" << std::endl;
      abort_handler(-1);
    }

    //bool update_exp_form = (quad_order != quadOrderPrev);
    // *** TO DO: capture updates to parameterized/numerical polynomials
    //if (update_exp_form) {
    //}

    // can't use quad_order > quadOrderPrev logic since only 1 pt set is stored
    bool update_basis_form = (quad_order != quadOrderPrev);
    // *** TO DO: carefully evaluate interdependence between exp_form/basis_form
    if (update_basis_form) {
      // size and initialize polynomialBasis, one interpolant per variable
      if (polynomialBasis.empty())
	{ polynomialBasis.resize(1); polynomialBasis[0].resize(numVars); }
      const Real2DArray& colloc_pts_1d = tpq_driver->collocation_points_array();
      std::vector<BasisPolynomial>& poly_basis_0 = polynomialBasis[0];
      short poly_type_1d; short rule; bool push_pts, found;
      distribution_types(poly_type_1d, rule, push_pts);
      for (i=0; i<numVars; ++i) {
	found = false;
	for (j=0; j<i; ++j)
	  if (colloc_pts_1d[i] == colloc_pts_1d[j]) // vector equality in pts
	    { found = true; break; }
	if (found) // reuse previous instance via shared representation
	  poly_basis_0[i] = poly_basis_0[j];
	else { // instantiate a new unique instance
	  poly_basis_0[i] = BasisPolynomial(poly_type_1d, rule);
	  //if (push_pts)
	    poly_basis_0[i].interpolation_points(colloc_pts_1d[i]);
	}
      }
    }

    quadOrderPrev = quad_order;
    break;
  }
  case SPARSE_GRID: {
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    unsigned short    ssg_level  = ssg_driver->level();
    const RealVector& aniso_wts  = ssg_driver->anisotropic_weights();

    //bool update_exp_form
    //  = (ssg_level != ssgLevelPrev || aniso_wts != ssgAnisoWtsPrev);
    // *** TO DO: capture updates to parameterized/numerical polynomials
    //if (update_exp_form) {
    //}

    // Ignore weights since they only reduce the interpolation depth from the
    // level and the basis update uses a coarse increment based on level.  This
    // matches isotropic sparse grids, but forces fewer and larger updates in
    // the case of anisotropic or generalized grids.
    bool update_basis_form
      = (ssgLevelPrev == USHRT_MAX || ssg_level > ssgLevelPrev);
    if (update_basis_form)
      update_sparse_interpolation_basis(ssg_level);

    ssgLevelPrev = ssg_level; ssgAnisoWtsPrev = aniso_wts;
    break;
  }
  }
}


void InterpPolyApproximation::compute_coefficients()
{
  if (!configOptions.expansionCoeffFlag &&
      !configOptions.expansionCoeffGradFlag) {
    PCerr << "Warning: neither expansion coefficients nor expansion "
	  << "coefficient gradients\n         are active in "
	  << "InterpPolyApproximation::compute_coefficients().\n         "
	  << "Bypassing approximation construction." << std::endl;
    return;
  }

  // For testing of anchorPoint logic:
  //anchorPoint = dataPoints.front();
  //dataPoints.pop_front();

  size_t i, offset = 0;
  numCollocPts = dataPoints.size();
  // anchorPoint, if present, is the first expansionSample.
  if (!anchorPoint.is_null()) {
    offset        = 1;
    numCollocPts += offset;
  }
  if (!numCollocPts) {
    PCerr << "Error: nonzero number of sample points required in "
	  << "InterpPolyApproximation::compute_coefficients()." << std::endl;
    abort_handler(-1);
  }

  allocate_arrays();

  if (!anchorPoint.is_null()) {
    if (configOptions.expansionCoeffFlag)
      expansionCoeffs[0] = anchorPoint.response_function();
    if (configOptions.expansionCoeffGradFlag)
      Teuchos::setCol(anchorPoint.response_gradient(), 0, expansionCoeffGrads);
  }

  std::vector<SurrogateDataPoint>::iterator it = dataPoints.begin();
  for (i=offset; i<numCollocPts; ++i, ++it) {
    if (configOptions.expansionCoeffFlag)
      expansionCoeffs[i] = it->response_function();
    if (configOptions.expansionCoeffGradFlag)
      Teuchos::setCol(it->response_gradient(), (int)i, expansionCoeffGrads);
  }

#ifdef INTERPOLATION_TEST
  // SC should accurately interpolate the collocation data for TPQ and
  // SSG with fully nested rules, but will exhibit interpolation error
  // for SSG with other rules.
  it = dataPoints.begin();
  for (i=offset; i<numCollocPts; ++i, ++it)
    PCout << "Colloc pt " << i+1 << ": coeff = " << expansionCoeffs[i]
	  << " interpolation error = " << std::abs(expansionCoeffs[i] -
	     get_value(it->continuous_variables())) << '\n';
#endif // INTERPOLATION_TEST
}


void InterpPolyApproximation::increment_coefficients()
{
  bool err_flag = false;
  switch (configOptions.expCoeffsSolnApproach) {
  case SPARSE_GRID: {
    // As for allocate_arrays(), increments are performed in coarser steps
    // than may be strictly necessary: all increments are filled in for all
    // vars for a step in level (ignoring anisotropy or generalized indices).
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    const UShortArray& trial_set = ssg_driver->trial_index_set();
    unsigned short max_trial_index = 0;
    for (size_t i=0; i<numVars; ++i)
      if (trial_set[i] > max_trial_index)
	max_trial_index = trial_set[i];
    update_sparse_interpolation_basis(max_trial_index);
    break;
  }
  default:
    err_flag = true; break;
  }
  if (err_flag) {
    PCerr << "Error: unsupported grid definition in InterpPolyApproximation::"
	  << "increment_coefficients()" << std::endl;
    abort_handler(-1);
  }

  restore_expansion_coefficients();
}


void InterpPolyApproximation::decrement_coefficients()
{
  // leave polynomialBasis as is

  switch (configOptions.expCoeffsSolnApproach) {
  case SPARSE_GRID: {
    // move previous expansion data to current expansion
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    savedSmolyakMultiIndex.push_back(ssg_driver->trial_index_set());
    break;
  }
  }

  // not necessary to prune; next increment/restore/finalize takes care of this
  //if (configOptions.expansionCoeffFlag)
  //  expansionCoeffs.resize(numCollocPts);
  //if (configOptions.expansionCoeffGradFlag) {
  //  size_t num_deriv_vars = expansionCoeffGrads.numRows();
  //  expansionCoeffGrads.reshape(num_deriv_vars, numCollocPts);
  //}

  numCollocPts = dataPoints.size(); // data already decremented
  if (!anchorPoint.is_null())
    numCollocPts += 1;
}


void InterpPolyApproximation::restore_coefficients()
{
  // leave polynomialBasis as is (a previous increment is being restored)

  switch (configOptions.expCoeffsSolnApproach) {
  case SPARSE_GRID: {
    // move previous expansion data to current expansion
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    std::deque<UShortArray>::iterator sit
      = std::find(savedSmolyakMultiIndex.begin(), savedSmolyakMultiIndex.end(),
		  ssg_driver->trial_index_set());
    if (sit != savedSmolyakMultiIndex.end())
      savedSmolyakMultiIndex.erase(sit);
    break;
  }
  }

  restore_expansion_coefficients();
}


void InterpPolyApproximation::finalize_coefficients()
{
  // leave polynomialBasis as is (all previous increments are being restored)

  switch (configOptions.expCoeffsSolnApproach) {
  case SPARSE_GRID:
    // move previous expansion data to current expansion
    savedSmolyakMultiIndex.clear();
    break;
  }

  restore_expansion_coefficients();
}


void InterpPolyApproximation::restore_expansion_coefficients()
{
  size_t i, new_colloc_pts = dataPoints.size();
  if (!anchorPoint.is_null())
    new_colloc_pts += 1;

  if (configOptions.expansionCoeffFlag)
    expansionCoeffs.resize(new_colloc_pts);
  if (configOptions.expansionCoeffGradFlag) {
    size_t num_deriv_vars = expansionCoeffGrads.numRows();
    expansionCoeffGrads.reshape(num_deriv_vars, new_colloc_pts);
  }

  std::vector<SurrogateDataPoint>::iterator it = dataPoints.begin();
  std::advance(it, numCollocPts);
  for (i=numCollocPts; i<new_colloc_pts; ++i, ++it) {
    if (configOptions.expansionCoeffFlag)
      expansionCoeffs[i] = it->response_function();
    if (configOptions.expansionCoeffGradFlag)
      Teuchos::setCol(it->response_gradient(), (int)i, expansionCoeffGrads);
  }

  numCollocPts = new_colloc_pts;
}


void InterpPolyApproximation::
update_sparse_interpolation_basis(unsigned short max_level)
{
  SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
  const Real3DArray& colloc_pts_1d = ssg_driver->collocation_points_array();

  // resize if needed (leaving previous levels unmodified)
  size_t i, j, k, basis_size = polynomialBasis.size();
  // j range is 0:w inclusive; i range is 1:w+1 inclusive
  unsigned short num_levels = max_level + 1;
  if (num_levels > basis_size) {
    polynomialBasis.resize(num_levels);
    for (i=basis_size; i<num_levels; ++i)
      polynomialBasis[i].resize(numVars);
  }

  // fill gaps that may exist within any level
  short poly_type_1d; short rule; bool push_pts, found;
  distribution_types(poly_type_1d, rule, push_pts);
  for (i=0; i<num_levels; ++i) { // i -> 0:num_levels-1 -> 0:ssg_level
    for (j=0; j<numVars; ++j) {
      const RealArray& colloc_pts_1d_ij =   colloc_pts_1d[i][j];
      BasisPolynomial&    poly_basis_ij = polynomialBasis[i][j];
      if (poly_basis_ij.is_null() && !colloc_pts_1d_ij.empty()) {
	found = false;
	for (k=0; k<j; ++k)
	  if (colloc_pts_1d_ij == colloc_pts_1d[i][k] &&  // vector equality
	      !polynomialBasis[i][k].is_null())
	    { found = true; break; }
	if (found) // reuse previous instances via shared representations
	  poly_basis_ij = polynomialBasis[i][k]; // shared rep
	else { // instantiate new unique instances
	  poly_basis_ij = BasisPolynomial(poly_type_1d, rule);
	  //if (push_pts)
	    poly_basis_ij.interpolation_points(colloc_pts_1d_ij);
	}
      }
    }
  }
}


void InterpPolyApproximation::compute_component_effects()
{
  // perform subset sort
  constituentSets.resize(sobolIndices.length());
  get_subsets();

  // initialize partialVariance
  if (partialVariance.empty())
    partialVariance.sizeUninitialized(sobolIndices.length());
  partialVariance = 0.;

  Real& mean           = numericalMoments[0];
  Real& total_variance = numericalMoments[1];
  partialVariance[0]   = mean*mean; // init with mean sq

  // Solve for partial variance
  for (IntIntMIter map_iter=sobolIndexMap.begin();
       map_iter!=sobolIndexMap.end(); ++map_iter) {
    // partialVariance[0] stores the mean; it is not a component function
    // and does not follow the procedures for obtaining variance 
    if (map_iter->first) {
      partial_variance(map_iter->first);
      sobolIndices[map_iter->second]
	= partialVariance[map_iter->second]/total_variance;
      // total indices simply identify the membership of the sobolIndices 
      // and adds it to the appropriate bin
    }
  }
#ifdef DEBUG
  PCout << "In InterpPolyApproximation::compute_component_effects(), "
	<< "sobolIndices =\n"; write_data(PCout, sobolIndices);
#endif // DEBUG
}


void InterpPolyApproximation::compute_total_effects()
{
  // iterate through existing indices if all component indices are available
  totalSobolIndices = 0; // init total indices
  if (configOptions.vbdControl == ALL_VBD) {
    for (IntIntMIter itr=sobolIndexMap.begin(); itr!=sobolIndexMap.end(); ++itr)
      for (int k=0; k<numVars; ++k) {
        if (itr->first & (1 << k))
          totalSobolIndices[k] += sobolIndices[itr->second];
        totalSobolIndices[k] = std::abs(totalSobolIndices[k]);
      }
  }
  // If not available, compute total indices independently.
  // This approach parallels partial_variance_integral where the algorithm is 
  // separated by integration approach.
  else {
    Real& total_variance = numericalMoments[1];
    int j, set_value;
    switch (configOptions.expCoeffsSolnApproach) {
      case QUADRATURE: {
        for (j=0; j<numVars; ++j) {
          // define set_value that includes all but index of interest
          set_value = (int)std::pow(2.,int(numVars)) - (int)std::pow(2.,j) - 1;
          totalSobolIndices[j]
	    = std::abs(1 - (total_effects_integral(set_value)/total_variance));
        }
        break;
      }
      case SPARSE_GRID: {
        SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
        const IntArray&   sm_coeffs  = ssg_driver->smolyak_coefficients();
        // Smolyak recursion of anisotropic tensor products
        size_t i, num_smolyak_indices = sm_coeffs.size();
        // iterate each variable 
        for (j=0; j<numVars; ++j) {
          set_value = (int)std::pow(2.,int(numVars)) - (int)std::pow(2.,j) - 1; 
          for (i=0; i<num_smolyak_indices; ++i)
	    if (sm_coeffs[i])
	      totalSobolIndices[j]
		+= sm_coeffs[i]*total_effects_integral(set_value,i);
          totalSobolIndices[j]
	    = std::abs(1 - (totalSobolIndices[j]/total_variance));
        }
        break;
      }
    }
  }
#ifdef DEBUG
  PCout << "In InterpPolyApproximation::compute_total_effects(), "
	<< "totalSobolIndices =\n"; write_data(PCout, totalSobolIndices);
#endif // DEBUG
}


Real InterpPolyApproximation::total_effects_integral(int set_value)
{
  // Some other routine here
  TensorProductDriver* tpq_driver    = (TensorProductDriver*)driverRep;
  const UShort2DArray& key           = tpq_driver->collocation_key();
  const Real2DArray&   colloc_wts_1d = tpq_driver->collocation_weights_array();
  const UShortArray&   quad_order    = tpq_driver->quadrature_order();;

  // Distinguish between non-members and members of the given set, set_value
  BoolDeque nonmember_vars(numVars,true); 
  int num_mem_exp_coeffs = 1; // number of expansion coeffs in
                                    // member-variable-only expansion 
  IntVector indexing_factor; // factors indexing member variables 
  indexing_factor.sizeUninitialized(numVars);

  indexing_factor = 1;
        
  for (int k=0; k<numVars; ++k) {
    // if subset contains variable k, set key for variable k to true
    if (set_value & (1 << k)) {
      nonmember_vars[k] = false;	
      // information to properly index mem_exp_coeffs
      indexing_factor[k] = num_mem_exp_coeffs;
      num_mem_exp_coeffs *= quad_order[k];
    }	
  }
        
  // Create vector to store new coefficients
  RealVector mem_exp_coeffs(num_mem_exp_coeffs),
    mem_weights(num_mem_exp_coeffs);
 
  // Perform integration over non-member variables and store indices
  // of new expansion
  size_t i, j;
  for (i=0; i<numCollocPts; ++i) {
    const UShortArray& key_i = key[i];
    size_t mem_exp_coeffs_index = 0;	
    Real prod_i_nonmembers = 1, prod_i_members = 1;
    for (j=0; j<numVars; ++j) {
      // Convert key to corresponding index on mem_exp_coeffs
      mem_exp_coeffs_index += (nonmember_vars[j]) ? 
        0 : key_i[j]*indexing_factor[j];
      // Save the product of the weights of the member and non-member variables 
      if (nonmember_vars[j])
        prod_i_nonmembers *= colloc_wts_1d[j][key_i[j]];
      else
        prod_i_members *= colloc_wts_1d[j][key_i[j]];
    }
 
    // mem_weights is performed more time than necessary here, but it
    // seems to be the simplest place to put it
    mem_weights[mem_exp_coeffs_index] = prod_i_members;
    // sort coefficients by the "signature" of the member variables
    // (i.e. mem_exp_coeffs_index)
    mem_exp_coeffs[mem_exp_coeffs_index]
      += expansionCoeffs[i]*prod_i_nonmembers;
  }
 
  // Now integrate over the remaining variables	
  Real& total_mean = numericalMoments[0];
  Real  integral   = 0;
  for (i=0; i<num_mem_exp_coeffs; ++i)
    integral += std::pow((mem_exp_coeffs[i]-total_mean),2.0)*mem_weights[i];
  return integral;
}


Real InterpPolyApproximation::
total_effects_integral(int set_value, size_t tp_index)
{
  SparseGridDriver*    ssg_driver = (SparseGridDriver*)driverRep;
  const UShortArray&     sm_index = ssg_driver->smolyak_multi_index()[tp_index];
  const UShort2DArray&        key = ssg_driver->collocation_key()[tp_index];
  const SizetArray&  colloc_index = ssg_driver->collocation_indices()[tp_index];
  const Real3DArray& colloc_wts_1d = ssg_driver->collocation_weights_array();

  // Distinguish between non-members and members of the given set, set_value
  BoolDeque nonmember_vars(numVars,true); 
  int num_mem_exp_coeffs = 1; // number of expansion coeffs in
                              // member-variable-only expansion 
  IntVector indexing_factor; // factors indexing member variables 
  indexing_factor.sizeUninitialized(numVars);

  // create member variable key and get number of expansion coeffs in
  // member-variable-only expansion
  indexing_factor = 1;
	
  UShortArray quad_order;
  ssg_driver->level_to_order(sm_index, quad_order);
  for (int k=0; k<numVars; ++k) {
    // if subset contains variable k, set key for variable k to true
    if (set_value & (1 << k)) {
      nonmember_vars[k] = false;	
      // information to properly index mem_exp_coeffs
      indexing_factor[k] = num_mem_exp_coeffs;
      num_mem_exp_coeffs *= quad_order[k];
    }	
  }
        
  // Create vector to store new coefficients
  RealVector mem_exp_coeffs(num_mem_exp_coeffs),
    mem_weights(num_mem_exp_coeffs);
 
  // Perform integration over non-member variables and store indices
  // of new expansion
  size_t i, j, num_colloc_pts = key.size();
  for (i=0; i <num_colloc_pts; ++i) {
    const UShortArray& key_i = key[i];
    size_t mem_exp_coeffs_index = 0;	
    Real prod_i_nonmembers = 1, prod_i_members = 1;
    for (j=0; j<numVars; ++j) {
      // Convert key to corresponding index on mem_exp_coeffs
      mem_exp_coeffs_index += (nonmember_vars[j]) ? 0 :
	key_i[j]*indexing_factor[j];
      // Save the product of the weights of the member and non-member variables 
      if (nonmember_vars[j])
        prod_i_nonmembers *= colloc_wts_1d[sm_index[j]][j][key_i[j]];
      else
        prod_i_members *= colloc_wts_1d[sm_index[j]][j][key_i[j]];
    }
 
    // mem_weights is performed more time than necessary here, but it
    // seems to be the simplest place to put it
    mem_weights[mem_exp_coeffs_index] = prod_i_members;
    // sort coefficients by the "signature" of the member variables
    // (i.e. mem_exp_coeffs_index)
    mem_exp_coeffs[mem_exp_coeffs_index]
      += expansionCoeffs[colloc_index[i]]*prod_i_nonmembers;
  }
 
  // Now integrate over the remaining variables	
  Real  integral   = 0;
  Real& total_mean = numericalMoments[0];
  for (i=0; i<num_mem_exp_coeffs; ++i)
    integral += std::pow(mem_exp_coeffs[i] - total_mean, 2.)*mem_weights[i];
  return integral;
}


/** Find constituent subsets. */
void InterpPolyApproximation::get_subsets()
{
  // includes the "zero" set
  //int num_subsets = sobolIndices.length(); 

  // Here we want to utilize the integer representation of the subset
  // but we want to store it in a size appropriate container
  // so finding lower sets is given the argument of the integer rep (->first)
  // and stored in constituentSets in size-appropriate-index-map (->second)
  for (IntIntMIter map_iter=sobolIndexMap.begin();
       map_iter!=sobolIndexMap.end(); ++map_iter) {
    lower_sets(map_iter->first, constituentSets[map_iter->second]);
    constituentSets[map_iter->second].erase(map_iter->first);
  }
}


/** For input set, recursively finds constituent subsets with one
    fewer element */
void InterpPolyApproximation::
lower_sets(int plus_one_set, IntSet& top_level_set)
{
  // if this set has been stored before, stop
  if (top_level_set.count(plus_one_set))
    return;
  // otherwise store current set
  else
    top_level_set.insert(plus_one_set);
  // and find lower level sets
  for (int k=0; k<numVars; ++k)
    // this performs a bitwise comparison by shifting 1 by k spaces 
    // and comparing that to a binary form of plus_one_set; this allows 
    // the variable membership using integers instead of a d-array of bools
    if (plus_one_set & (1 << k)) 
      // if subset i contains variable k, remove that variable from the set 
      // by converting the bit-form of (1<<k) to an integer and subtract from
      // the plus_one_set
      lower_sets(plus_one_set-(int)std::pow(2.0,k),top_level_set);
}


/** Forms an interpolant over variables that are members of the given set.
    Finds the variance of the interpolant w.r.t. the variables in the set.
    Overloaded version supporting tensor-product quadrature. */
Real InterpPolyApproximation::partial_variance_integral(int set_value)
{
  TensorProductDriver* tpq_driver    = (TensorProductDriver*)driverRep;
  const UShort2DArray& key           = tpq_driver->collocation_key();
  const Real2DArray&   colloc_wts_1d = tpq_driver->collocation_weights_array();
  const UShortArray&   quad_order    = tpq_driver->quadrature_order();;

  // Distinguish between non-members and members of the given set, set_value
  BoolDeque nonmember_vars(numVars,true); 
  int num_mem_exp_coeffs = 1; // number of expansion coeffs in
                              // member-variable-only expansion 
  IntVector indexing_factor; // factors indexing member variables 
  indexing_factor.sizeUninitialized(numVars);

  // create member variable key and get number of expansion coeffs in
  // member-variable-only expansion
  indexing_factor = 1;
	
  for (int k=0; k<numVars; ++k) {
    // if subset contains variable k, set key for variable k to true
    if (set_value & (1 << k)) {
      nonmember_vars[k] = false;	
      // information to properly index mem_exp_coeffs
      indexing_factor[k] = num_mem_exp_coeffs;
      num_mem_exp_coeffs *= quad_order[k];
    }	
  }
	
  // Create vector to store new coefficients
  RealVector mem_exp_coeffs(num_mem_exp_coeffs),
    mem_weights(num_mem_exp_coeffs);

  // Perform integration over non-member variables and store indices
  // of new expansion
  size_t i, j;
  for (i=0; i<numCollocPts; ++i) {
    const UShortArray& key_i = key[i];
    size_t mem_exp_coeffs_index = 0;	
    Real prod_i_nonmembers = 1, prod_i_members = 1;
    for (j=0; j<numVars; ++j) {
      // Convert key to corresponding index on mem_exp_coeffs
      mem_exp_coeffs_index += (nonmember_vars[j]) ? 0 :
	key_i[j]*indexing_factor[j];
      // Save the product of the weights of the member and non-member variables 
      if (nonmember_vars[j])
	prod_i_nonmembers *= colloc_wts_1d[j][key_i[j]];
      else
	prod_i_members *= colloc_wts_1d[j][key_i[j]];
    }

    // mem_weights is performed more time than necessary here, but it
    // seems to be the simplest place to put it
    mem_weights[mem_exp_coeffs_index] = prod_i_members;
    // sort coefficients by the "signature" of the member variables
    // (i.e. mem_exp_coeffs_index)
    mem_exp_coeffs[mem_exp_coeffs_index]
      += expansionCoeffs[i]*prod_i_nonmembers;
  }

  // Now integrate over the remaining variables	
  Real integral = 0;
  for (i=0; i<num_mem_exp_coeffs; ++i)
    integral += std::pow(mem_exp_coeffs[i], 2.0)*mem_weights[i];
  return integral;	
}


/** Forms an interpolant over variables that are members of the given set.
    Finds the variance of the interpolant w.r.t. the variables in the set.
    Overloaded version supporting Smolyak sparse grids. */
Real InterpPolyApproximation::
partial_variance_integral(int set_value, size_t tp_index)
{
  SparseGridDriver*    ssg_driver = (SparseGridDriver*)driverRep;
  const UShortArray&     sm_index = ssg_driver->smolyak_multi_index()[tp_index];
  const UShort2DArray&        key = ssg_driver->collocation_key()[tp_index];
  const SizetArray&  colloc_index = ssg_driver->collocation_indices()[tp_index];
  const Real3DArray& colloc_wts_1d = ssg_driver->collocation_weights_array();

  // Distinguish between non-members and members of the given set, set_value
  BoolDeque nonmember_vars(numVars,true); 
  int num_mem_exp_coeffs = 1; // number of expansion coeffs in
                              // member-variable-only expansion 
  IntVector indexing_factor; // factors indexing member variables 
  indexing_factor.sizeUninitialized(numVars);

  // create member variable key and get number of expansion coeffs in
  // member-variable-only expansion
  indexing_factor = 1;
	
  UShortArray quad_order;
  ssg_driver->level_to_order(sm_index, quad_order);
  for (int k=0; k<numVars; ++k) {
    // if subset contains variable k, set key for variable k to true
    if (set_value & (1 << k)) {
      nonmember_vars[k] = false;	
      // information to properly index mem_exp_coeffs
      indexing_factor[k] = num_mem_exp_coeffs;
      num_mem_exp_coeffs *= quad_order[k];
    }	
  }
	
  // Create vector to store new coefficients
  RealVector mem_exp_coeffs(num_mem_exp_coeffs),
    mem_weights(num_mem_exp_coeffs);

  // Perform integration over non-member variables and store indices
  // of new expansion
  size_t i, j, num_colloc_pts = key.size();
  for (i=0; i <num_colloc_pts; ++i) {
    const UShortArray& key_i = key[i];
    size_t mem_exp_coeffs_index = 0;	
    Real prod_i_nonmembers = 1, prod_i_members = 1;
    for (j=0; j<numVars; ++j) {
      // Convert key to corresponding index on mem_exp_coeffs
      mem_exp_coeffs_index += (nonmember_vars[j]) ? 
	0 : key_i[j]*indexing_factor[j];
      // Save the product of the weights of the member and non-member variables 
      if (nonmember_vars[j])
	prod_i_nonmembers *= colloc_wts_1d[sm_index[j]][j][key_i[j]];
      else
	prod_i_members *= colloc_wts_1d[sm_index[j]][j][key_i[j]];
    }

    // mem_weights is performed more time than necessary here, but it
    // seems to be the simplest place to put it
    mem_weights[mem_exp_coeffs_index] = prod_i_members;
    // sort coefficients by the "signature" of the member variables
    // (i.e. mem_exp_coeffs_index)
    mem_exp_coeffs[mem_exp_coeffs_index]
      += expansionCoeffs[colloc_index[i]]*prod_i_nonmembers;
  }

  // Now integrate over the remaining variables	
  Real integral = 0;
  for (i=0; i<num_mem_exp_coeffs; ++i)
    integral += std::pow(mem_exp_coeffs[i], 2.0)*mem_weights[i];
  return integral;	
}


/** Computes the variance of component functions. Assumes that all
    subsets of set_value have been computed in advance which will be
    true so long as the partial_variance is called following
    appropriate enumeration of set value  */
void InterpPolyApproximation::partial_variance(int set_value)
{
  // Computes the integral first
  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE:
    partialVariance[sobolIndexMap[set_value]]
      = partial_variance_integral(set_value);
    break;
  case SPARSE_GRID: {
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    const IntArray&   sm_coeffs  = ssg_driver->smolyak_coefficients();
    // Smolyak recursion of anisotropic tensor products
    size_t num_smolyak_indices = sm_coeffs.size();
    for (size_t i=0; i<num_smolyak_indices; ++i)
      if (sm_coeffs[i])
	partialVariance[sobolIndexMap[set_value]] += sm_coeffs[i]
	  * partial_variance_integral(set_value, i);
    break;
  }
  }
  
  // Now subtract the contributions from constituent subsets
  IntSet::iterator itr;
  for (itr  = constituentSets[sobolIndexMap[set_value]].begin();
       itr != constituentSets[sobolIndexMap[set_value]].end(); ++itr) 
    partialVariance[sobolIndexMap[set_value]]
      -= partialVariance[sobolIndexMap[*itr]];
}

} // namespace Pecos

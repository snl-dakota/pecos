/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        SharedOrthogPolyApproxData
//- Description:  Class for Multivariate Orthogonal Polynomial Approximations
//-               
//- Owner:        Mike Eldred

#ifndef SHARED_ORTHOG_POLY_APPROX_DATA_HPP
#define SHARED_ORTHOG_POLY_APPROX_DATA_HPP

#include "SharedPolyApproxData.hpp"
#include "NumericGenOrthogPolynomial.hpp"
#include "CombinedSparseGridDriver.hpp"

namespace Pecos {


/// Derived approximation class for orthogonal polynomials (global
/// approximation).

/** The SharedOrthogPolyApproxData class provides a global approximation
    based on orthogonal polynomials.  It is used primarily for polynomial
    chaos expansions (for stochastic finite element approaches to
    uncertainty quantification). */

class SharedOrthogPolyApproxData: public SharedPolyApproxData
{
  //
  //- Heading: Friends
  //

  friend class OrthogPolyApproximation;

public:

  //
  //- Heading: Constructor and destructor
  //

  /// standard constructor
  SharedOrthogPolyApproxData(short basis_type, const UShortArray& approx_order,
			     size_t num_vars);
  /// alternate constructor
  SharedOrthogPolyApproxData(short basis_type, const UShortArray& approx_order,
			     size_t num_vars,
			     const ExpansionConfigOptions& ec_options,
			     const BasisConfigOptions&     bc_options);
  /// destructor
  ~SharedOrthogPolyApproxData();

  //
  //- Heading: Member functions
  //

  /// alternative form for setting multiIndex (expansion import) and
  /// allocating related arrays
  void allocate_data(const UShort2DArray& mi);
  /// get multiIndex
  const UShort2DArray& multi_index() const;

  /// retrieve size of multiIndex
  size_t expansion_terms() const;

  /// get approxOrder
  const UShortArray& expansion_order() const;
  /// set approxOrder
  void expansion_order(const UShortArray& order);
  /// uniformly increment approxOrder
  void increment_order();

  /// invoke initialize_orthog_poly_basis_types(),
  /// initialize_polynomial_basis(), and, if needed,
  /// update_basis_distribution_parameters()
  void construct_basis(const ShortArray& u_types,
		       const AleatoryDistParams& adp);

  /// set orthogPolyTypes
  void orthog_poly_basis_types(const ShortArray& opb_types);
  /// get orthogPolyTypes
  const ShortArray& orthog_poly_basis_types() const;

  /// get polynomialBasis
  const std::vector<BasisPolynomial>& polynomial_basis() const;
  /// set polynomialBasis
  void polynomial_basis(const std::vector<BasisPolynomial>& poly_basis);

  /// set NumericGenOrthogPolynomial::coeffsNormsFlag
  void coefficients_norms_flag(bool flag);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void allocate_data();

  void store_data();
  void pre_combine_data(short combine_type);
  void post_combine_data(short combine_type);

  //
  //- Heading: Member functions
  //

  /// initialize orthogonal polynomial basis types from u_types
  bool initialize_orthog_poly_basis_types(const ShortArray& u_types,
					  ShortArray& opb_types);

  /// convert a sparse grid index set and a growth setting to an integrand_order
  void sparse_grid_level_to_expansion_order(
    CombinedSparseGridDriver* csg_driver, const UShortArray& levels,
    UShortArray& exp_order);//, short growth_rate = UNRESTRICTED_GROWTH);
  /// convert quadrature orders to integrand orders using rigorous mappings
  void quadrature_order_to_integrand_order(IntegrationDriver* int_driver,
					   const UShortArray& quad_order,
					   UShortArray& int_order);
  /// convert integrand orders to expansion orders using rigorous mappings
  void integrand_order_to_expansion_order(const UShortArray& int_order,
					  UShortArray& exp_order);

  /// helper function for incrementing that is modular on sparse grid driver
  /// and multi-index
  void increment_trial_set(CombinedSparseGridDriver* csg_driver,
			   UShort2DArray& aggregated_mi);
  /// helper function for incrementing that is modular on trial set and
  /// multi-index
  void increment_trial_set(const UShortArray& trial_set,
			   UShort2DArray& aggregated_mi);
  /// helper function for decrementing that is modular on trial set
  /// and multi-index
  void decrement_trial_set(const UShortArray& trial_set,
			   UShort2DArray& aggregated_mi);
  /// helper function for restoring that is modular on trial set and multi-index
  void pre_restore_trial_set(const UShortArray& trial_set,
			     UShort2DArray& aggregated_mi);
  /// helper function for restoring that is modular on trial set and multi-index
  void post_restore_trial_set(const UShortArray& trial_set,
			      UShort2DArray& aggregated_mi);
  /// helper function for restoring that is modular on trial set and multi-index
  void restore_trial_set(const UShortArray& trial_set,
			 UShort2DArray& aggregated_mi);
  /// update cvErrorRef and bestExpTerms after improvement in solution
  void update_reference(Real cv_error, const UShort2DArray& aggregated_mi);
  /// update cvErrorRef and bestExpTerms after improvement in solution
  void restore_best_solution(UShort2DArray& aggregated_mi);

  /// update approxOrder and total-order multiIndex
  void allocate_total_order();
  /// allocate sobolIndexMap from multi_index
  void allocate_component_sobol(const UShort2DArray& multi_index);
  /// update sobolIndexMap using new multi_index terms (from multifidelity
  /// overlay or a new QoI in orthogonal least interpolation)
  void update_component_sobol(const UShort2DArray& multi_index);

  /// append multi-indices from append_mi that do not already appear
  /// in combined_mi
  void append_multi_index(const UShort2DArray& append_mi,
			  UShort2DArray& combined_mi);
  /// append multi-indices from append_mi that do not already appear
  /// in combined_mi
  void append_multi_index(const UShortArraySet& append_mi,
			  UShort2DArray& combined_mi);
  /// append multi-indices from append_mi that do not already appear
  /// in combined_mi; define append_mi_map and append_mi_map_ref
  void append_multi_index(const UShort2DArray& append_mi,
			  UShort2DArray& combined_mi, SizetArray& append_mi_map,
			  size_t& append_mi_map_ref);
  /// append multi-indices from append_mi that do not already appear
  /// in combined_mi; define append_mi_map and append_mi_map_ref
  void append_multi_index(const UShort2DArray& append_mi,
			  UShort2DArray& combined_mi, SizetSet& append_mi_map,
			  size_t& append_mi_map_ref);
  /// append multi-indices from append_mi that do not already appear in
  /// combined_mi (consistent ordering assumed); define append_mi_map
  /// and append_mi_map_ref
  void append_leading_multi_index(const UShort2DArray& append_mi,
				  UShort2DArray& combined_mi,
				  SizetSet& append_mi_map,
				  size_t& append_mi_map_ref);
  /// append multi-indices from append_mi that do not already appear
  /// in combined_mi, using previously defined append_mi_map and
  /// append_mi_map_ref for mapping
  void append_multi_index(const UShort2DArray& append_mi,
			  SizetArray& append_mi_map, size_t& append_mi_map_ref,
			  UShort2DArray& combined_mi);
  /// append multi-indices from append_mi that do not already appear in
  /// combined_mi, updating sparse_indices, exp_coeffs, and exp_coeff_grads
  void append_multi_index(SizetSet& sparse_indices,
			  const UShort2DArray& append_mi,
			  UShort2DArray& combined_mi, RealVector& exp_coeffs,
			  RealMatrix& exp_coeff_grads);

  /// returns the norm-squared of a particular multivariate polynomial,
  /// treating all variables as probabilistic
  Real norm_squared(const UShortArray& indices);
  /// returns the norm-squared of a particular multivariate polynomial,
  /// treating a subset of the variables as probabilistic
  Real norm_squared(const UShortArray& indices, const SizetList& rand_indices);

  /// calculate a particular multivariate orthogonal polynomial value
  /// evaluated at a particular parameter set
  Real multivariate_polynomial(const RealVector& x,const UShortArray& indices);
  /// calculate a particular multivariate orthogonal polynomial value over
  /// the nonrandom variable subset evaluated at a particular parameter set
  Real multivariate_polynomial(const RealVector& x,const UShortArray& indices,
			       const SizetList& non_rand_indices);

  /// compute multivariate orthogonal polynomial gradient for term
  /// corresponding to deriv_index, evaluated at x
  Real multivariate_polynomial_gradient(const RealVector& x, size_t deriv_index,
    const UShortArray& indices);
  /// compute multivariate orthogonal polynomial gradient for term corresponding
  /// to deriv_index, for nonrandom variable subset, and evaluated at x
  Real multivariate_polynomial_gradient(const RealVector& x, size_t deriv_index,
    const UShortArray& indices, const SizetList& non_rand_indices);

  /// calculate multivariate orthogonal polynomial gradient vector
  /// evaluated at a particular parameter set
  const RealVector& multivariate_polynomial_gradient_vector(const RealVector& x,
    const UShortArray& indices);
  /// calculate multivariate orthogonal polynomial gradient vector with
  /// respect to specified dvv and evaluated at a particular parameter set
  const RealVector& multivariate_polynomial_gradient_vector(const RealVector& x,
    const UShortArray& indices, const SizetArray& dvv);

  /// update the combined Pareto set with a new multi_index
  void update_pareto(const UShort2DArray& multi_index,
		     UShort2DArray& combined_pareto);
  /// update the combined Pareto set with a new multi_index
  void update_pareto(const UShort2DArray& multi_index,
		     UShortArraySet& combined_pareto);
  /// update the combined Pareto set with a new multi_index term
  void update_pareto(const UShortArray& mi_i, UShortArraySet& combined_pareto);
  // assess whether new_pareto is dominated by total_pareto
  //bool assess_dominance(const UShort2DArray& pareto,
  //			  const UShort2DArray& combined_pareto);
  /// assess bi-directional dominance for a new polynomial index set 
  /// against an incumbent polynomial index set
  void assess_dominance(const UShortArray& new_order,
			const UShortArray& existing_order,
			bool& new_dominated, bool& existing_dominated);

  /// test for nonzero indices in random variable subset
  bool zero_random(const UShortArray& indices) const;

  /// Generate the coefficient tag for variable j of given expansion term order
  void get_tag(char* tag, size_t j, unsigned short order) const;

  /// tests 1D gradient computations (active in DEBUG compile mode)
  void gradient_check();

  //
  //- Heading: Data
  //

  /// array of basis types for each one-dimensional orthogonal polynomial:
  /// HERMITE_ORTHOG, LEGENDRE_ORTHOG, LAGUERRE_ORTHOG, JACOBI_ORTHOG,
  /// GEN_LAGUERRE_ORTHOG, CHEBYSHEV_ORTHOG, or NUM_GEN_ORTHOG
  ShortArray orthogPolyTypes;

  /// array of one-dimensional basis polynomial objects which are used in
  /// constructing the multivariate orthogonal/interpolation polynomials
  std::vector<BasisPolynomial> polynomialBasis;

  /// order of orthogonal polynomial expansion
  UShortArray approxOrder;
  /// previous value of approxOrder; used for detecting when a multiIndex
  /// update is needed
  UShortArray approxOrderPrev;
  /// copy of approxOrder stored in store_coefficients() for use in
  /// combine_coefficients()
  UShortArray storedApproxOrder;

  /// number of exp terms-by-number of vars array for identifying the orders
  /// of the one-dimensional orthogonal polynomials contributing to each
  /// of the multivariate orthogonal polynomials
  UShort2DArray multiIndex;
  /// copy of multiIndex (aggregated expansion) stored in store_coefficients()
  /// for use in combine_coefficients()
  UShort2DArray storedMultiIndex;
  /// mapping of terms when aggregating storedMultiIndex with multiIndex in
  /// pre_combine_data()
  SizetArray storedMultiIndexMap;
  /// multi-index that is the result of expansion combination
  UShort2DArray combinedMultiIndex;

  /// numSmolyakIndices-by-numTensorProductPts-by-numVars array for
  /// identifying the orders of the one-dimensional orthogonal polynomials
  /// contributing to each of the multivariate orthogonal polynomials.
  /** For nested rules (GP, CC, or GK), the integration driver's collocKey
      is insufficient and we must track expansion orders separately. */
  UShort3DArray tpMultiIndex;
  /// sparse grid bookkeeping: mapping from num tensor-products by 
  /// tensor-product multi-indices into aggregated multiIndex
  Sizet2DArray tpMultiIndexMap;
  /// sparse grid bookkeeping: reference points for tpMultiIndexMap
  SizetArray tpMultiIndexMapRef;

  /// saved instances of tpMultiIndex that were computed but not selected
  std::deque<UShort2DArray> savedTPMultiIndex;
  /// saved instances of tpMultiIndexMap that were computed but not selected
  std::deque<SizetArray> savedTPMultiIndexMap;
  /// saved instances of tpMultiIndexMapRef that were computed but not selected
  std::deque<size_t> savedTPMultiIndexMapRef;

  /// index into saved sets of data to be restored (stored in this
  /// class for used by each ProjectOrthogPolyApproximation)
  size_t restoreIndex;

  /// sparse grid driver for adapting a CS candidate basis using greedy
  /// adaptation via the generalized sparse grid algorithm; it's state
  /// is reset for each response QoI
  CombinedSparseGridDriver csgDriver;
  /// the cross validation error reference point for adapting a CS
  /// candidate basis; it's state is reset for each response QoI
  Real cvErrorRef;
  /// size of expansion that corresponds to the best solution identified
  /// (may not be the last solution)
  size_t bestExpTerms;
  /// a scalar growth factor for defining the tpMultiIndex contribution
  /// from a particular trial index set in adapted basis mode
  unsigned short multiIndexGrowthFactor;
  /// initial sparse grid level that provides the starting point for basis
  /// adaptation
  unsigned short referenceSGLevel;
  /// flag indicating whether to apply normalization to CV error estimates
  /// within select_best_active()
  bool normalizeCV;
  /// number of consecutive cycles for which convergence criterion
  /// must be met prior to termination
  unsigned short softConvLimit;
  /// number of consecutive cycles for which convergence criterion
  /// must be met prior to termination
  unsigned short numAdvancements;

  /// Data vector for storing the gradients of individual expansion term
  /// polynomials (see multivariate_polynomial_gradient_vector())
  RealVector mvpGradient;

private:

  //
  //- Heading: Member functions
  //

};


inline SharedOrthogPolyApproxData::
SharedOrthogPolyApproxData(short basis_type, const UShortArray& approx_order,
			   size_t num_vars):
  SharedPolyApproxData(basis_type, num_vars), approxOrder(approx_order),
  multiIndexGrowthFactor(2),
  // paper consistency:
  referenceSGLevel(2), normalizeCV(false), softConvLimit(3), numAdvancements(3)
  // normal settings:
  //referenceSGLevel(0), normalizeCV(true), softConvLimit(2), numAdvancements(2)
{ }


inline SharedOrthogPolyApproxData::
SharedOrthogPolyApproxData(short basis_type, const UShortArray& approx_order,
			   size_t num_vars,
			   const ExpansionConfigOptions& ec_options,
			   const BasisConfigOptions&     bc_options):
  SharedPolyApproxData(basis_type, num_vars, ec_options, bc_options),
  approxOrder(approx_order), multiIndexGrowthFactor(2),
  // paper consistency:
  referenceSGLevel(2), normalizeCV(false), softConvLimit(3), numAdvancements(3)
  // normal settings:
  //referenceSGLevel(0), normalizeCV(true), softConvLimit(2), numAdvancements(2)
{ }


inline SharedOrthogPolyApproxData::~SharedOrthogPolyApproxData()
{ }


inline const UShort2DArray& SharedOrthogPolyApproxData::multi_index() const
{ return multiIndex; }


inline size_t SharedOrthogPolyApproxData::expansion_terms() const
{ return multiIndex.size(); }


inline const UShortArray& SharedOrthogPolyApproxData::expansion_order() const
{ return approxOrder; }


inline void SharedOrthogPolyApproxData::
expansion_order(const UShortArray& order)
{ approxOrder = order; } // multiIndex updated in allocate_arrays()


inline void SharedOrthogPolyApproxData::increment_order()
{
  // increment approxOrder (multiIndex updated in allocate_arrays())
  for (size_t i=0; i<numVars; ++i)
    approxOrder[i] += 1;
}


/** This function is invoked to create orthogPolyTypes and polynomialBasis
    for cases where they have not already been created by an
    IntegrationDriver (i.e., expansion_samples or regression). */
inline void SharedOrthogPolyApproxData::
construct_basis(const ShortArray& u_types, const AleatoryDistParams& adp)
{
  bool dist_params
    = initialize_orthog_poly_basis_types(u_types, orthogPolyTypes);
  ShortArray colloc_rules;
  initialize_collocation_rules(u_types, basisConfigOptions, colloc_rules);
  initialize_polynomial_basis(orthogPolyTypes, colloc_rules, polynomialBasis);
  if (dist_params)
    update_basis_distribution_parameters(u_types, adp, polynomialBasis);
}


inline void SharedOrthogPolyApproxData::
orthog_poly_basis_types(const ShortArray& opb_types)
{ orthogPolyTypes = opb_types; }


inline const ShortArray& SharedOrthogPolyApproxData::
orthog_poly_basis_types() const
{ return orthogPolyTypes; }


inline const std::vector<BasisPolynomial>& SharedOrthogPolyApproxData::
polynomial_basis() const
{ return polynomialBasis; }


inline void SharedOrthogPolyApproxData::
polynomial_basis(const std::vector<BasisPolynomial>& poly_basis)
{ polynomialBasis = poly_basis; }


inline void SharedOrthogPolyApproxData::coefficients_norms_flag(bool flag)
{
  size_t i, num_basis = orthogPolyTypes.size();
  for (i=0; i<num_basis; ++i)
    if (orthogPolyTypes[i] == NUM_GEN_ORTHOG)
      ((NumericGenOrthogPolynomial*)polynomialBasis[i].polynomial_rep())
	->coefficients_norms_flag(flag);
}


inline void SharedOrthogPolyApproxData::
restore_trial_set(const UShortArray& trial_set, UShort2DArray& aggregated_mi)
{
  pre_restore_trial_set(trial_set, aggregated_mi);
  post_restore_trial_set(trial_set, aggregated_mi);
}


inline void SharedOrthogPolyApproxData::
update_reference(Real cv_error, const UShort2DArray& aggregated_mi)
{
  cvErrorRef   = cv_error;
  bestExpTerms = aggregated_mi.size();
}


inline Real SharedOrthogPolyApproxData::
multivariate_polynomial(const RealVector& x, const UShortArray& indices)
{
  Real mvp = 1.; unsigned short order_1d;
  for (size_t i=0; i<numVars; ++i) {
    order_1d = indices[i];
    if (order_1d)
      mvp *= polynomialBasis[i].type1_value(x[i], order_1d);
  }
  return mvp;
}


/** All variables version. */
inline Real SharedOrthogPolyApproxData::
multivariate_polynomial(const RealVector& x, const UShortArray& indices,
			const SizetList& non_rand_indices)
{
  Real mvp = 1.; SizetList::const_iterator cit;
  unsigned short order_1d; size_t i;
  for (cit=non_rand_indices.begin(); cit!=non_rand_indices.end(); ++cit) {
    i = *cit; order_1d = indices[i];
    if (order_1d)
      mvp *= polynomialBasis[i].type1_value(x[i], order_1d);
  }
  return mvp;
}


inline Real SharedOrthogPolyApproxData::
multivariate_polynomial_gradient(const RealVector& x, size_t deriv_index,
				 const UShortArray& indices)
{
  Real mvp_grad = 1.;
  for (size_t k=0; k<numVars; ++k)
    mvp_grad *= (k == deriv_index) ?
      polynomialBasis[k].type1_gradient(x[k], indices[k]) :
      polynomialBasis[k].type1_value(x[k],    indices[k]);
  return mvp_grad;
}


/** All variables version. */
inline Real SharedOrthogPolyApproxData::
multivariate_polynomial_gradient(const RealVector& x, size_t deriv_index,
				 const UShortArray& indices,
				 const SizetList& non_rand_indices)
{
  Real mvp_grad = 1.; SizetList::const_iterator cit; size_t k;
  for (cit=non_rand_indices.begin(); cit!=non_rand_indices.end(); ++cit) {
    k = *cit;
    mvp_grad *= (k == deriv_index) ?
      polynomialBasis[k].type1_gradient(x[k], indices[k]) :
      polynomialBasis[k].type1_value(x[k],    indices[k]);
  }
  return mvp_grad;
}


inline const RealVector& SharedOrthogPolyApproxData::
multivariate_polynomial_gradient_vector(const RealVector& x,
					const UShortArray& indices)
{
  if (mvpGradient.length() != numVars)
    mvpGradient.sizeUninitialized(numVars);
  for (size_t i=0; i<numVars; ++i)
    mvpGradient[i] = multivariate_polynomial_gradient(x, i, indices);
  return mvpGradient;
}


inline const RealVector& SharedOrthogPolyApproxData::
multivariate_polynomial_gradient_vector(const RealVector& x,
					const UShortArray& indices,
					const SizetArray& dvv)
{
  size_t i, j, deriv_index, num_deriv_vars = dvv.size();
  if (mvpGradient.length() != num_deriv_vars)
    mvpGradient.sizeUninitialized(num_deriv_vars);
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // *** requires an "All" view
    mvpGradient[i] = multivariate_polynomial_gradient(x, deriv_index, indices);
  }
  return mvpGradient;
}


inline bool SharedOrthogPolyApproxData::
zero_random(const UShortArray& indices) const
{
  SizetList::const_iterator cit;
  for (cit=randomIndices.begin(); cit!=randomIndices.end(); ++cit)
    if (indices[*cit])
      return false;
  return true;
}


inline Real SharedOrthogPolyApproxData::norm_squared(const UShortArray& indices)
{
  // the norm squared of a particular multivariate polynomial is the product of
  // the norms squared of the numVars univariate polynomials that comprise it.
  Real norm_sq = 1.; unsigned short order_1d;
  for (size_t i=0; i<numVars; ++i) {
    order_1d = indices[i];
    if (order_1d)
      norm_sq *= polynomialBasis[i].norm_squared(order_1d);
  }
  return norm_sq;
}


/** All variables version. */
inline Real SharedOrthogPolyApproxData::
norm_squared(const UShortArray& indices, const SizetList& rand_indices)
{
  // the norm squared of a particular multivariate polynomial is the product of
  // the norms squared of the numVars univariate polynomials that comprise it.
  Real norm_sq = 1.; SizetList::const_iterator cit;
  unsigned short order_1d; size_t i;
  for (cit=rand_indices.begin(); cit!=rand_indices.end(); ++cit) {
    i = *cit; order_1d = indices[i];
    if (order_1d)
      norm_sq *= polynomialBasis[i].norm_squared(order_1d);
  }
  return norm_sq;
}

} // namespace Pecos

#endif

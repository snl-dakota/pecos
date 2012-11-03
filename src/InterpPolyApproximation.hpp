/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        InterpPolyApproximation
//- Description:  Class for Lagrange Interpolation Polynomial Approximation
//-               
//- Owner:        Mike Eldred

#ifndef INTERP_POLY_APPROXIMATION_HPP
#define INTERP_POLY_APPROXIMATION_HPP

#include "PolynomialApproximation.hpp"
#include "BasisPolynomial.hpp"

namespace Pecos {


/// Derived approximation class for interpolation polynomials (global
/// approximation).

/** The InterpPolyApproximation class provides a global approximation
    based on interpolation polynomials.  It is used primarily for
    stochastic collocation approaches to uncertainty quantification. */

class InterpPolyApproximation: public PolynomialApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  InterpPolyApproximation(short basis_type, size_t num_vars, bool use_derivs);
  /// destructor
  ~InterpPolyApproximation();

  //
  //- Heading: member functions
  //

  /// define n-D basis types from u_types and mode booleans for use in
  /// defining quadrature rules used by integration drivers
  /** These basis types may include orthogonal polynomials for
      purposes of computing their Gauss points and weights within
      integration drivers; thus they differ in general from the
      interpolation polynomial basis used for approximation. */
  static bool initialize_integration_basis_types(const ShortArray& u_types,
    const BasisConfigOptions& bc_options, ShortArray& basis_types);
  /// initialize basis types and collocation rules, construct the polynomial
  /// basis, and initialize the distribution parameters within this basis
  static void construct_basis(const ShortArray& u_types,
			      const DistributionParams& dp,
			      const BasisConfigOptions& bc_options,
			      std::vector<Pecos::BasisPolynomial>& poly_basis);

  //
  //- Heading: Virtual function redefinitions
  //

  /// compute the coefficients for the expansion of multivariate
  /// interpolation polynomials
  void compute_coefficients();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  int min_coefficients() const;

  /// update the coefficients for the expansion of multivariate Lagrange
  /// interpolation polynomials
  void increment_coefficients();
  /// restore the coefficients to their previous state prior to last increment
  void decrement_coefficients();
  /// restore the coefficients to a previously incremented state as
  /// identified by the current increment to the Smolyak multi index
  void restore_coefficients();
  /// finalize the coefficients by applying all previously evaluated increments
  void finalize_coefficients();

  /// size expansionType{1,2}Coeffs and expansionType1CoeffGrads
  void allocate_arrays();

  /// computes component (main and interaction) effect Sobol' indices
  void compute_component_effects();
  /// computes total effect Sobol' indices
  void compute_total_effects();

  /// compute numerical moments to order 4
  void compute_moments();
  /// compute numerical moments in all-variables mode to order 2
  void compute_moments(const RealVector& x);
  /// return numericalMoments
  const RealVector& moments() const;

  //
  //- Heading: New virtual functions
  //

  /// derived portion of allocate_arrays()
  virtual void allocate_expansion_coefficients() = 0;
  /// derived portion of compute_coefficients()
  virtual void compute_expansion_coefficients() = 0;
  /// increment expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  /// within increment_coefficients()
  virtual void increment_expansion_coefficients() = 0;
  /// decrement expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  /// within decrement_coefficients()
  virtual void decrement_expansion_coefficients() = 0;
  /// restore expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  /// within restore_coefficients()
  virtual void restore_expansion_coefficients() = 0;
  /// finalize expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  /// within finalize_coefficients()
  virtual void finalize_expansion_coefficients() = 0;

  /// compute moments of response using numerical integration
  virtual void compute_numerical_response_moments(size_t num_moments) = 0;
  /// compute moments of expansion using numerical integration
  virtual void compute_numerical_expansion_moments(size_t num_moments) = 0;

  virtual void compute_total_sobol_indices() = 0;
  virtual void compute_partial_variance(const BitArray& set_value);

  //
  //- Heading: Convenience functions
  //

  /// return value of type 1 interpolation polynomial using all dimensions
  Real type1_interpolant_value(const RealVector& x, const UShortArray& key,
			       const UShortArray& basis_index);
  /// return value of type 1 interpolation polynomial using interpolated
  /// (non-integrated) dimension subset
  Real type1_interpolant_value(const RealVector& x, const UShortArray& key,
			       const UShortArray& basis_index,
			       const SizetList& non_rand_indices);
  /// return gradient of type 1 interpolation polynomial using all dimensions
  Real type1_interpolant_gradient(const RealVector& x, size_t deriv_index,
				  const UShortArray& key,
				  const UShortArray& basis_index);
  /// return gradient of type 1 interpolation polynomial using interpolated
  /// (non-integrated) dimension subset
  Real type1_interpolant_gradient(const RealVector& x, size_t deriv_index,
				  const UShortArray& key,
				  const UShortArray& basis_index,
				  const SizetList& non_rand_indices);

  /// return value of type 2 interpolation polynomial using all dimensions
  Real type2_interpolant_value(const RealVector& x, size_t interp_index,
			       const UShortArray& key,
			       const UShortArray& basis_index);
  /// return value of type 2 interpolation polynomial using interpolated
  /// (non-integrated) dimension subset
  Real type2_interpolant_value(const RealVector& x, size_t interp_index,
			       const UShortArray& key,
			       const UShortArray& basis_index,
			       const SizetList& non_rand_indices);
  /// return gradient of type 2 interpolation polynomial using all dimensions
  Real type2_interpolant_gradient(const RealVector& x, size_t deriv_index,
				  size_t interp_index, const UShortArray& key,
				  const UShortArray& basis_index);
  /// return gradient of type 2 interpolation polynomial using interpolated
  /// (non-integrated) dimension subset
  Real type2_interpolant_gradient(const RealVector& x, size_t deriv_index,
				  size_t interp_index, const UShortArray& key,
				  const UShortArray& basis_index,
				  const SizetList& non_rand_indices);

  /// return type 1 product weight from integration of type 1 interpolation
  /// polynomials using integrated dimension subset
  Real type1_weight(const UShortArray& key, const UShortArray& basis_index, 
		    const SizetList& rand_indices);
  /// return type 2 product weight from integration of type 1/2 interpolation
  /// polynomials using integrated dimension subset
  Real type2_weight(size_t interp_index, const UShortArray& key,
		    const UShortArray& basis_index,
		    const SizetList& rand_indices);

  /// compute the value of a tensor interpolant on a tensor grid;
  /// contributes to value(x)
  Real tensor_product_value(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs,
    const UShortArray& basis_index,  const UShort2DArray& key,
    const SizetArray&  colloc_index);

  /// compute the gradient of a tensor interpolant on a tensor grid
  /// with respect to variables that are included in the polynomial
  /// basis; contributes to gradient_basis_variables(x)
  const RealVector& tensor_product_gradient_basis_variables(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs,
    const UShortArray& basis_index,  const UShort2DArray& key,
    const SizetArray&  colloc_index);
  /// compute the gradient of a tensor interpolant on a tensor grid
  /// with respect to variables that are included in the polynomial
  /// basis for given DVV; contributes to gradient_basis_variables(x, dvv)
  const RealVector& tensor_product_gradient_basis_variables(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs,
    const UShortArray& basis_index,  const UShort2DArray& key,
    const SizetArray& colloc_index,  const SizetArray& dvv);
  /// compute the gradient of a tensor interpolant on a tensor grid
  /// with respect to variables that are not included in the
  /// polynomial basis; contributes to gradient_nonbasis_variables(x)
  const RealVector& tensor_product_gradient_nonbasis_variables(
    const RealVector& x,             const RealMatrix& exp_t1_coeff_grads,
    const UShortArray& basis_index,  const UShort2DArray& key,
    const SizetArray& colloc_index);

  /// for a particular level, find index of basis v2 that matches basis v1
  bool find_basis(unsigned short level, size_t v1, size_t& v2);
  /// for a particular level, test for equality between basis v2 and basis v1
  bool same_basis(unsigned short level, size_t v1, size_t v2);

  /// compute total Sobol effects for an index within a sparse grid
  Real total_effects_integral(const BitArray& set_value,
			      const UShortArray& quad_order,
			      const UShortArray& lev_index,
			      const UShort2DArray& key,
			      const SizetArray& colloc_index);
  /// finds variance of sparse interpolant with respect to variables in the set
  Real partial_variance_integral(const BitArray& set_value,
				 const UShortArray& quad_order,
				 const UShortArray& lev_index,
				 const UShort2DArray& key,
				 const SizetArray& colloc_index);
  /// defines member_coeffs and member_wts for a particular set_value
  void member_coefficients_weights(const BitArray& set_value,
    const UShortArray& quad_order, const UShortArray& lev_index,
    const UShort2DArray& key, const SizetArray& colloc_index,
    RealVector& member_coeffs, RealVector& member_wts);

  //
  //- Heading: Data
  //

  /// 2D array of one-dimensional basis polynomial objects used in
  /// constructing the multivariate orthogonal/interpolation polynomials.
  /** Each variable (inner array size = numVars) has multiple
      integration orders associated with it (outer array size = max
      quadrature order nfor TPQ or sparse grid level + 1 for SSG). */
  std::vector<std::vector<BasisPolynomial> > polynomialBasis;

  /// total number of collocation points = number of type 1 terms in
  /// interpolation expansion (length of expansionType1Coeffs)
  int numCollocPts;

  /// the gradient of a tensor-product interpolant; a contributor to
  /// approxGradient
  RealVector tpGradient;
  /// the gradient of the mean of a tensor-product interpolant; a
  /// contributor to meanGradient
  RealVector tpMeanGrad;
  /// the gradient of the variance of a tensor-product interpolant; a
  /// contributor to varianceGradient
  RealVector tpVarianceGrad;

  /// the partialVariances of subset functions f_alpha
  RealVector partialVariance;

  /// {GLOBAL,PIECEWISE}_{NODAL,HIERARCHICAL}_INTERPOLATION_POLYNOMIAL or
  /// {GLOBAL,PIECEWISE}_ORTHOGONAL_POLYNOMIAL
  short basisType;

private:

  //
  //- Heading: Convenience functions
  //

  /// define the 1D basis type and collocation rule
  void initialize_polynomial_basis_type(short& poly_type_1d, short& rule);

  /// update polynomialBasis after a change in quadrature order
  void update_tensor_interpolation_basis();
  /// update polynomialBasis after a change in sparse grid level
  void update_sparse_interpolation_basis(unsigned short max_level);

  /// recursively identifies constituent subsets that are children of
  /// a parent set
  void proper_subsets(const BitArray& parent_set, BitArraySet& children);

  //
  //- Heading: Data
  //

};


inline InterpPolyApproximation::
InterpPolyApproximation(short basis_type, size_t num_vars, bool use_derivs):
  PolynomialApproximation(num_vars, use_derivs), numCollocPts(0),
  basisType(basis_type)
{ }


inline InterpPolyApproximation::~InterpPolyApproximation()
{ }


inline void InterpPolyApproximation::
construct_basis(const ShortArray& u_types, const DistributionParams& dp,
		const BasisConfigOptions& bc_options,
		std::vector<Pecos::BasisPolynomial>& poly_basis)
{
  ShortArray basis_types, colloc_rules;
  bool dist_params
    = initialize_integration_basis_types(u_types, bc_options, basis_types);
  initialize_collocation_rules(u_types, bc_options, colloc_rules);
  initialize_polynomial_basis(basis_types, colloc_rules, poly_basis);
  if (dist_params)
    update_basis_distribution_parameters(u_types, dp, poly_basis);
}


inline void InterpPolyApproximation::compute_moments()
{
  // standard variables mode supports four moments
  compute_numerical_response_moments(4);
  //if (expConfigOptions.outputLevel >= VERBOSE_OUTPUT)
    compute_numerical_expansion_moments(4);
  // numericalMoments and expansionMoments should be the same for faithful
  // interpolants (TPQ and SSG with nested rules), but will differ for other
  // cases (SSG with non-nested rules)
}


inline void InterpPolyApproximation::compute_moments(const RealVector& x)
{
  // all variables mode only supports first two moments
  mean(x); variance(x);
  //standardize_moments(numericalMoments);
  //compute_numerical_expansion_moments(4, x);

  // Note: it would be feasible to implement an all_variables version of
  // compute_numerical_expansion_moments() by evaluating the combined
  // expansion at {design/epistemic=initialPtU,aleatory=Gauss points}
  // --> cannot do this for compute_numerical_response_moments()
  //     (lacking required response data)
  // --> would require generation of new TPQ/SSG grid only over aleatory vars
  // --> could possibly retire redundant all_vars functions
}


inline const RealVector& InterpPolyApproximation::moments() const
{ return numericalMoments; }


inline Real InterpPolyApproximation::
type1_interpolant_value(const RealVector& x, const UShortArray& key,
			const UShortArray& basis_index)
{
  Real L1 = 1.;
  for (size_t j=0; j<numVars; ++j)
    L1 *= polynomialBasis[basis_index[j]][j].type1_value(x[j], key[j]);
  return L1;
}


/** All variables version. */
inline Real InterpPolyApproximation::
type1_interpolant_value(const RealVector& x, const UShortArray& key,
			const UShortArray& basis_index,
			const SizetList& non_rand_indices)
{
  Real L1 = 1.; SizetList::const_iterator cit; size_t j;
  for (cit=non_rand_indices.begin(); cit!=non_rand_indices.end(); ++cit) {
    j   = *cit;
    L1 *= polynomialBasis[basis_index[j]][j].type1_value(x[j], key[j]);
  }
  return L1;
}


inline Real InterpPolyApproximation::
type1_interpolant_gradient(const RealVector& x, size_t deriv_index,
			   const UShortArray& key,
			   const UShortArray& basis_index)
{
  Real L1_grad = 1.;
  for (size_t k=0; k<numVars; ++k)
    L1_grad *= (k == deriv_index) ?
      polynomialBasis[basis_index[k]][k].type1_gradient(x[k], key[k]) :
      polynomialBasis[basis_index[k]][k].type1_value(x[k],    key[k]);
  return L1_grad;
}


/** All variables version. */
inline Real InterpPolyApproximation::
type1_interpolant_gradient(const RealVector& x, size_t deriv_index,
			   const UShortArray& key,
			   const UShortArray& basis_index,
			   const SizetList& non_rand_indices)
{
  Real L1_grad = 1.; SizetList::const_iterator cit; size_t k;
  for (cit=non_rand_indices.begin(); cit!=non_rand_indices.end(); ++cit) {
    k        = *cit;
    L1_grad *= (k == deriv_index) ?
      polynomialBasis[basis_index[k]][k].type1_gradient(x[k], key[k]) :
      polynomialBasis[basis_index[k]][k].type1_value(x[k], key[k]);
  }
  return L1_grad;
}


inline Real InterpPolyApproximation::
type2_interpolant_value(const RealVector& x,    size_t interp_index,
			const UShortArray& key, const UShortArray& basis_index)
{
  Real L2 = 1.;
  for (size_t k=0; k<numVars; ++k)
    L2 *= (interp_index == k) ?
      polynomialBasis[basis_index[k]][k].type2_value(x[k], key[k]) :
      polynomialBasis[basis_index[k]][k].type1_value(x[k], key[k]);
  return L2;
}


/** All variables version. */
inline Real InterpPolyApproximation::
type2_interpolant_value(const RealVector& x,    size_t interp_index,
			const UShortArray& key, const UShortArray& basis_index,
			const SizetList& non_rand_indices)
{
  Real L2 = 1.; SizetList::const_iterator cit; size_t k;
  for (cit=non_rand_indices.begin(); cit!=non_rand_indices.end(); ++cit) {
    k   = *cit;
    L2 *= (interp_index == k) ?
      polynomialBasis[basis_index[k]][k].type2_value(x[k], key[k]) :
      polynomialBasis[basis_index[k]][k].type1_value(x[k], key[k]);
  }
  return L2;
}


inline Real InterpPolyApproximation::
type2_interpolant_gradient(const RealVector& x, size_t deriv_index,
			   size_t interp_index, const UShortArray& key,
			   const UShortArray& basis_index)
{
  // deriv_index  = desired gradient component
  // interp_index = index of gradient component used in type2 interpolation
  Real L2_grad = 1.;
  if (interp_index == deriv_index) // match in grad and type2 interp components
    for (size_t l=0; l<numVars; ++l)
      L2_grad *= (l == interp_index) ?
	polynomialBasis[basis_index[l]][l].type2_gradient(x[l], key[l]) :
	polynomialBasis[basis_index[l]][l].type1_value(x[l],    key[l]);
  else                          // mismatch in grad and type2 interp components
    for (size_t l=0; l<numVars; ++l)
      L2_grad *= (l == interp_index) ?
	polynomialBasis[basis_index[l]][l].type2_value(x[l],    key[l]) :
	polynomialBasis[basis_index[l]][l].type1_gradient(x[l], key[l]);
  return L2_grad;
}


/** All variables version. */
inline Real InterpPolyApproximation::
type2_interpolant_gradient(const RealVector& x, size_t deriv_index,
			   size_t interp_index, const UShortArray& key,
			   const UShortArray& basis_index,
			   const SizetList& non_rand_indices)
{
  // deriv_index  = desired gradient component
  // interp_index = index of gradient component used in type2 interpolation
  Real L2_grad = 1.; SizetList::const_iterator cit; size_t l;
  if (interp_index == deriv_index) // match in grad and type2 interp components
    for (cit=non_rand_indices.begin(); cit!=non_rand_indices.end(); ++cit) {
      l        = *cit;
      L2_grad *= (l == interp_index) ?
	polynomialBasis[basis_index[l]][l].type2_gradient(x[l], key[l]) :
	polynomialBasis[basis_index[l]][l].type1_value(x[l],    key[l]);
    }
  else                          // mismatch in grad and type2 interp components
    for (cit=non_rand_indices.begin(); cit!=non_rand_indices.end(); ++cit) {
      l        = *cit;
      L2_grad *= (l == interp_index) ?
	polynomialBasis[basis_index[l]][l].type2_value(x[l],    key[l]) :
	polynomialBasis[basis_index[l]][l].type1_gradient(x[l], key[l]);
    }
  return L2_grad;
}


/** All variables partial weight. */
inline Real InterpPolyApproximation::
type1_weight(const UShortArray& key, const UShortArray& basis_index, 
	     const SizetList& rand_indices)
{
  const Real3DArray& t1_wts_1d = driverRep->type1_collocation_weights_array();
  Real t1_wt_prod = 1.; SizetList::const_iterator cit; size_t j;
  for (cit=rand_indices.begin(); cit!=rand_indices.end(); ++cit)
    { j = *cit; t1_wt_prod *= t1_wts_1d[basis_index[j]][j][key[j]]; }
  return t1_wt_prod;
}


/** All variables partial weight. */
inline Real InterpPolyApproximation::
type2_weight(size_t interp_index, const UShortArray& key,
	     const UShortArray& basis_index, const SizetList& rand_indices)
{
  const Real3DArray& t1_wts_1d = driverRep->type1_collocation_weights_array();
  const Real3DArray& t2_wts_1d = driverRep->type2_collocation_weights_array();
  Real t2_wt_prod = 1.; SizetList::const_iterator cit; size_t j;
  for (cit=rand_indices.begin(); cit!=rand_indices.end(); ++cit) {
    j           = *cit;
    t2_wt_prod *= (interp_index == j) ? t2_wts_1d[basis_index[j]][j][key[j]]
                                      : t1_wts_1d[basis_index[j]][j][key[j]];
  }
  return t2_wt_prod;
}

} // namespace Pecos

#endif

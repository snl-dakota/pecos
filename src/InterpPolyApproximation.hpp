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

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  int min_coefficients() const;

  /// compute the coefficients for the expansion of multivariate Lagrange
  /// interpolation polynomials
  virtual void compute_coefficients();
  /// update the coefficients for the expansion of multivariate Lagrange
  /// interpolation polynomials
  virtual void increment_coefficients();
  /// restore the coefficients to their previous state prior to last increment
  void decrement_coefficients();
  /// restore the coefficients to a previously incremented state as
  /// identified by the current increment to the Smolyak multi index
  void restore_coefficients();
  /// finalize the coefficients by applying all previously evaluated increments
  void finalize_coefficients();
  /// store current state within storedExpType{1Coeffs,2Coeffs,1CoeffGrads},
  /// storedColloc{Key,Indices}, and storedLev{MultiIndex,Coeffs}
  void store_coefficients();
  /// augment current interpolant using
  /// storedExpType{1Coeffs,2Coeffs,1CoeffGrads}, storedColloc{Key,Indices},
  /// and storedLev{MultiIndex,Coeffs}
  void combine_coefficients(short combine_type);

  const RealVector& approximation_coefficients() const;
  void approximation_coefficients(const RealVector& approx_coeffs);

  /// size expansionType{1,2}Coeffs and expansionType1CoeffGrads
  void allocate_arrays();

  /// compute moments of response using numerical integration
  void compute_numerical_response_moments(size_t num_moments);
  /// compute moments of expansion using numerical integration
  void compute_numerical_expansion_moments(size_t num_moments);

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

  /// the type1 coefficients of the expansion for interpolating values
  RealVector expansionType1Coeffs;
  /// the type2 coefficients of the expansion for interpolating gradients
  RealMatrix expansionType2Coeffs;
  /// the gradients of the type1 expansion coefficients
  /** may be interpreted as either the gradients of the expansion
      coefficients or the coefficients of expansions for the response
      gradients.  This array is used when sensitivities of moments are
      needed with respect to variables that do not appear in the
      expansion (e.g., with respect to design variables for an
      expansion only over the random variables). */
  RealMatrix expansionType1CoeffGrads;

  /// storage of expansionType1Coeffs state for subsequent restoration
  RealVector storedExpType1Coeffs;
  /// storage of expansionType2Coeffs state for subsequent restoration
  RealMatrix storedExpType2Coeffs;
  /// storage of expansionType1CoeffGrads state for subsequent restoration
  RealMatrix storedExpType1CoeffGrads;
  /// storage of IntegrationDriver combinatorial coefficients state
  /// for subsequent restoration
  IntArray storedLevCoeffs;
  /// storage of IntegrationDriver collocation key state for
  /// subsequent restoration
  UShort3DArray storedCollocKey;
  /// storage of IntegrationDriver collocation indices state for
  /// subsequent restoration
  Sizet2DArray storedCollocIndices;

  /// the gradient of a tensor-product interpolant; a contributor to
  /// approxGradient
  RealVector tpGradient;
  /// the gradient of the mean of a tensor-product interpolant; a
  /// contributor to meanGradient
  RealVector tpMeanGrad;
  /// the gradient of the variance of a tensor-product interpolant; a
  /// contributor to varianceGradient
  RealVector tpVarianceGrad;

  /// GLOBAL_INTERPOLATION_POLYNOMIAL, PIECEWISE_NODAL_INTERPOLATION_POLYNOMIAL,
  /// or PIECEWISE_HIERARCHICAL_INTERPOLATION_POLYNOMIAL
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

  /// restore expansion{Coeffs,CoeffGrads} within increment/restore/finalize
  void restore_expansion_coefficients();

  /// compute total Sobol effects for an index within a sparse grid
  Real total_effects_integral(int set_value, const UShortArray& quad_order,
			      const UShortArray& lev_index,
			      const UShort2DArray& key,
			      const SizetArray& colloc_index);
  /// performs sorting to store constituent subsets (constituentSets)
  void get_subsets();
  /// recursively identifies constituent subsets
  void lower_sets(int plus_one_set, IntSet& top_level_set);
  /// computes partialVariance
  void partial_variance(int set_value);
  /// finds variance of sparse interpolant with respect to variables in the set
  Real partial_variance_integral(int set_value, const UShortArray& quad_order,
				 const UShortArray& lev_index,
				 const UShort2DArray& key,
				 const SizetArray& colloc_index);

  //
  //- Heading: Data
  //

  /// the constituent subsets for each superset
  std::vector<IntSet> constituentSets;
  /// the partialVariances of subset functions f_alpha
  RealVector partialVariance;
};


inline InterpPolyApproximation::
InterpPolyApproximation(short basis_type, size_t num_vars, bool use_derivs):
  PolynomialApproximation(num_vars, use_derivs), numCollocPts(0),
  basisType(basis_type)
{ }


inline InterpPolyApproximation::~InterpPolyApproximation()
{ }


inline const RealVector& InterpPolyApproximation::
approximation_coefficients() const
{
  if (basisConfigOptions.useDerivs) {
    PCerr << "Error: approximation_coefficients() not supported in "
	  << "InterpPolyApproximation for type2 coefficients." << std::endl;
    return abort_handler_t<const RealVector&>(-1);
  }
  else
    return expansionType1Coeffs;
}


inline void InterpPolyApproximation::
approximation_coefficients(const RealVector& approx_coeffs)
{
  if (basisConfigOptions.useDerivs) {
    PCerr << "Error: approximation_coefficients() not supported in "
	  << "InterpPolyApproximation for type2 coefficients." << std::endl;
    abort_handler(-1);
  }
  else
    expansionType1Coeffs = approx_coeffs;
}


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
  numericalMoments.sizeUninitialized(2);
  numericalMoments[0] = mean(x); numericalMoments[1] = variance(x);
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

} // namespace Pecos

#endif

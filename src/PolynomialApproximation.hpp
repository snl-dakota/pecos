/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        PolynomialApproximation
//- Description:  Base Class for Orthogonal/Interpolation Polynomial
//-               Approximations
//-               
//- Owner:        Mike Eldred

#ifndef POLYNOMIAL_APPROXIMATION_HPP
#define POLYNOMIAL_APPROXIMATION_HPP

#include "BasisApproximation.hpp"
#include "SurrogateDataPoint.hpp"


namespace Pecos {

class IntegrationDriver;

/// Derived approximation class for global basis polynomials.

/** The PolynomialApproximation class provides a global approximation
    based on basis polynomials.  This includes orthogonal polynomials
    used for polynomial chaos expansions and interpolation polynomials
    used for stochastic collocation. */

class PolynomialApproximation: public BasisApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  PolynomialApproximation();
  /// standard constructor
  PolynomialApproximation(size_t num_vars, short output_level);
  /// destructor
  ~PolynomialApproximation();

  //
  //- Heading: Virtual member functions
  //

  /// size Sobol arrays
  virtual void allocate_arrays();

  /// Performs global sensitivity analysis using Sobol' Indices
  virtual void compute_global_sensitivity() = 0;

  /// return the mean of the expansion, treating all variables as random
  virtual const Real& get_mean() = 0;
  /// return the mean of the expansion for a given parameter vector,
  /// treating a subset of the variables as random
  virtual const Real& get_mean(const RealVector& x) = 0;
  /// return the gradient of the expansion mean for a given parameter
  /// vector, treating all variables as random
  virtual const RealVector& get_mean_gradient() = 0;
  /// return the gradient of the expansion mean for a given parameter vector
  /// and given DVV, treating a subset of the variables as random
  virtual const RealVector& get_mean_gradient(const RealVector& x,
					      const UIntArray& dvv) = 0;

  /// return the variance of the expansion, treating all variables as random
  virtual const Real& get_variance() = 0;
  /// return the variance of the expansion for a given parameter vector,
  /// treating a subset of the variables as random
  virtual const Real& get_variance(const RealVector& x) = 0;
  /// return the gradient of the expansion variance for a given parameter
  /// vector, treating all variables as random
  virtual const RealVector& get_variance_gradient() = 0;
  /// return the gradient of the expansion variance for a given parameter
  /// vector and given DVV, treating a subset of the variables as random
  virtual const RealVector& get_variance_gradient(const RealVector& x,
						  const UIntArray& dvv) = 0;

  /// return the variance of the expansion, treating all variables as random
  virtual const Real& get_covariance(const RealVector& exp_coeffs_2) = 0;
  // return the variance of the expansion for a given parameter vector,
  // treating a subset of the variables as random
  //virtual const Real& get_covariance(const RealVector& x,
  //                                   const RealVector& exp_coeffs_2) = 0;

  //
  //- Heading: Member functions
  //

  /// set dataPoints
  void data_points(const std::list<Pecos::SurrogateDataPoint>& pts);
  /// set anchorPoint
  void anchor_point(const Pecos::SurrogateDataPoint& pt);
  /// queries the status of anchorPoint
  bool anchor() const;

  /// set expCoeffsSolnApproach
  void solution_approach(short soln_approach);
  /// get expCoeffsSolnApproach
  short solution_approach() const;

  /// set expansionCoeffFlag
  void expansion_coefficient_flag(bool coeff_flag);
  /// get expansionCoeffFlag
  bool expansion_coefficient_flag() const;

  /// set expansionCoeffGradFlag
  void expansion_coefficient_gradient_flag(bool grad_flag);
  /// get expansionCoeffGradFlag
  bool expansion_coefficient_gradient_flag() const;

  /// return sobolIndexMap 
  const IntIntMap& sobol_index_map() const;
  /// return sobolIndices
  const RealVector& sobol_indices() const;
  /// return totalSobolIndices
  const RealVector& total_sobol_indices() const;

  /// set driverRep
  void integration_driver_rep(IntegrationDriver* driver_rep);

  /// set randomVarsKey
  void random_variables_key(const BoolDeque& random_vars_key);

  /// return the number of expansion terms for a total order expansion
  /// with the provided (anisotropic) upper_bound array specification
  static size_t total_order_terms(const UShortArray& upper_bound,
				  short lower_bound_offset = -1);
  /// return the number of expansion terms for a tensor-product expansion
  /// with the provided (anisotropic) quadrature orders (default) or
  /// expansion orders (offset = true)
  static size_t tensor_product_terms(const UShortArray& order,
				     bool exp_order_offset = false);

  /// initialize expansion multi_index using a tensor-product expansion
  static void tensor_product_multi_index(const UShortArray& order,
					 UShort2DArray& multi_index,
					 bool exp_order_offset = false);
  /// initialize expansion multi_index using a total-order expansion
  /// from an upper_bound array specification
  static void total_order_multi_index(const UShortArray& upper_bound,
				      UShort2DArray& multi_index,
				      short lower_bound_offset = -1,
				      unsigned int max_terms = UINT_MAX);

  /// utility function for incrementing a set of multidimensional indices
  static void increment_indices(UShortArray& indices, const UShortArray& limits,
				bool include_limit_equality);
  /// utility function for incrementing a set of multidimensional terms
  static void increment_terms(UShortArray& terms, size_t& last_index,
			      size_t& prev_index, const size_t& term_limit,
			      bool& order_complete);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  const RealVector& approximation_coefficients() const;
  void approximation_coefficients(const RealVector& approx_coeffs);

  //
  //- Heading: Member functions
  //

  //
  //- Heading: Data
  //

  /// a special sample (often at the center of the approximation region)
  /// for which exact matching is enforced (e.g., using equality-constrained
  /// least squares regression).
  SurrogateDataPoint anchorPoint;
  /// set of samples used to build the approximation.  These sample points
  /// are fit approximately (e.g., using least squares regression); exact
  /// matching is not enforced.
  std::vector<SurrogateDataPoint> dataPoints;

  /// pointer to integration driver instance
  IntegrationDriver* driverRep;

  /// identifies the approach taken in find_coefficients():
  /// QUADRATURE, CUBATURE, SPARSE_GRID, REGRESSION, or SAMPLING
  short expCoeffsSolnApproach;
  /// SILENT_OUTPUT, QUIET_OUTPUT, NORMAL_OUTPUT, VERBOSE_OUTPUT, or
  /// DEBUG_OUTPUT
  short outputLevel;

  /// flag for calculation of expansionCoeffs from response values
  bool expansionCoeffFlag;
  /// flag for calculation of expansionCoeffGrads from response gradients
  bool expansionCoeffGradFlag;

  /// previous quadrature order;
  /// used for tracking need for expansion form updates
  UShortArray quadOrderPrev;
  /// previous Smolyak sparse grid level;
  /// used for tracking need for expansion form updates
  unsigned short ssgLevelPrev;
  /// previous Smolyak sparse grid anisotropic weighting;
  /// used for tracking need for expansion form updates
  RealVector ssgAnisoWtsPrev;

  /// array of booleans identifying the random variable subset within
  /// the active variables (used in all_variables mode)
  BoolDeque randomVarsKey;
  /// list of indices identifying the random variable subset within the active
  /// variables (used in all_variables mode; defined from randomVarsKey)
  SizetList randomIndices;
  /// list of indices identifying the non-random variable subset within the
  /// active variables (used in all_variables mode; defined from randomVarsKey)
  SizetList nonRandomIndices;

  /// expected value of the expansion
  Real expansionMean;
  /// gradient of the expected value of the expansion
  RealVector expansionMeanGrad;
  /// variance of the expansion
  Real expansionVariance;
  /// gradient of the variance of the expansion
  RealVector expansionVarianceGrad;

  /// the coefficients of the expansion
  RealVector expansionCoeffs;
  /// the gradients of the expansion coefficients
  /** may be interpreted as either the gradients of the expansion
      coefficients or the coefficients of expansions for the response
      gradients.  This array is used when sensitivities of moments are
      needed with respect to variables that do not appear in the
      expansion (e.g., with respect to design variables for an
      expansion only over the random variables). */
  RealMatrix expansionCoeffGrads;

  /// introduce mapping to unify disparate enumeration of sensitivity
  /// indices
  IntIntMap sobolIndexMap;
  /// global sensitivities as given by Sobol'
  RealVector sobolIndices;
  /// total global sensitivities as given by Sobol'
  RealVector totalSobolIndices;

private:

  //
  //- Heading: Data
  //
};


inline PolynomialApproximation::PolynomialApproximation():
  driverRep(NULL), expCoeffsSolnApproach(SAMPLING), outputLevel(NORMAL_OUTPUT),
  expansionCoeffFlag(true), expansionCoeffGradFlag(false),
  ssgLevelPrev(USHRT_MAX)
{ }


inline PolynomialApproximation::
PolynomialApproximation(size_t num_vars, short output_level):
  BasisApproximation(BaseConstructor(), num_vars), driverRep(NULL),
  expCoeffsSolnApproach(SAMPLING), outputLevel(output_level),
  expansionCoeffFlag(true), expansionCoeffGradFlag(false),
  ssgLevelPrev(USHRT_MAX)
{ }


inline PolynomialApproximation::~PolynomialApproximation()
{ }


inline void PolynomialApproximation::
data_points(const std::list<Pecos::SurrogateDataPoint>& pts)
{
  size_t i, num_pts = pts.size();
  dataPoints.resize(num_pts);
  std::list<Pecos::SurrogateDataPoint>::const_iterator cit;
  for (i=0, cit=pts.begin(); i<num_pts; ++i, ++cit)
    dataPoints[i] = *cit; 
}


inline void PolynomialApproximation::
anchor_point(const Pecos::SurrogateDataPoint& pt)
{ anchorPoint = pt; }


inline bool PolynomialApproximation::anchor() const
{ return !anchorPoint.is_null(); }


inline void PolynomialApproximation::solution_approach(short soln_approach)
{ expCoeffsSolnApproach = soln_approach; }


inline short PolynomialApproximation::solution_approach() const
{ return expCoeffsSolnApproach; }


inline void PolynomialApproximation::expansion_coefficient_flag(bool coeff_flag)
{ expansionCoeffFlag = coeff_flag; }


inline bool PolynomialApproximation::expansion_coefficient_flag() const
{ return expansionCoeffFlag; }


inline void PolynomialApproximation::
expansion_coefficient_gradient_flag(bool grad_flag)
{ expansionCoeffGradFlag = grad_flag; }


inline bool PolynomialApproximation::
expansion_coefficient_gradient_flag() const
{ return expansionCoeffGradFlag; }


inline const IntIntMap& PolynomialApproximation::sobol_index_map() const
{ return sobolIndexMap; }


inline const RealVector& PolynomialApproximation::sobol_indices() const
{ return sobolIndices; }


inline const RealVector& PolynomialApproximation::total_sobol_indices() const
{ return totalSobolIndices; }


inline void PolynomialApproximation::
integration_driver_rep(IntegrationDriver* driver_rep)
{ driverRep = driver_rep; }


inline void PolynomialApproximation::
random_variables_key(const BoolDeque& random_vars_key)
{
  randomVarsKey = random_vars_key;
  randomIndices.clear();
  nonRandomIndices.clear();
  for (size_t i=0; i<numVars; i++)
    if (random_vars_key[i])
      randomIndices.push_back(i);
    else
      nonRandomIndices.push_back(i);
}


inline void PolynomialApproximation::
increment_indices(UShortArray& indices, const UShortArray& limits,
		  bool include_limit_equality)
{
  size_t n = indices.size(), increment_index = 0;
  ++indices[increment_index];
  while ( increment_index < n &&
	  ( (  include_limit_equality && 
	       indices[increment_index] >= limits[increment_index] ) ||
	    ( !include_limit_equality && 
	       indices[increment_index] >  limits[increment_index] ) ) ) {
    indices[increment_index] = 0;
    ++increment_index;
    if (increment_index < n)
      ++indices[increment_index];
  }
}


inline void PolynomialApproximation::
increment_terms(UShortArray& terms, size_t& last_index, size_t& prev_index,
		const size_t& term_limit, bool& order_complete)
{
  bool increment_complete = false;
  while (!increment_complete) {
    terms[last_index] = 1;
    ++terms[prev_index];
    if (prev_index == 0) {
      increment_complete = true;
      if (terms[prev_index] > term_limit)
	order_complete = true;
    }
    else {
      last_index = prev_index;
      --prev_index;
      if (terms[last_index] <= terms[prev_index])
	increment_complete = true;
    }
  }
}


inline const RealVector& PolynomialApproximation::
approximation_coefficients() const
{ return expansionCoeffs; }


inline void PolynomialApproximation::
approximation_coefficients(const RealVector& approx_coeffs)
{ expansionCoeffs = approx_coeffs; }

} // namespace Pecos

#endif

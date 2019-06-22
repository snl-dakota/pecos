/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef MARGINALS_CORR_DISTRIBUTION_HPP
#define MARGINALS_CORR_DISTRIBUTION_HPP

#include "MultivariateDistribution.hpp"
#include "RandomVariable.hpp"


namespace Pecos {


/// Class for multivariate distribution based on 1D marginals + correlations

/** The MarginalsCorrDistribution models a multivariate random variable
    distribution using 1D marginals plus a correlation matrix. */

class MarginalsCorrDistribution: public MultivariateDistribution
{
public:

  //
  //- Heading: Constructors and destructor
  //

  MarginalsCorrDistribution();  ///< constructor
  ~MarginalsCorrDistribution(); ///< destructor

  //
  //- Heading: Virtual function redefinitions
  //

  /// return randomVars[i]
  const RandomVariable& random_variable(size_t i) const;
  /// return randomVars
  const std::vector<RandomVariable>& random_variables() const;
  /// return randomVars
  std::vector<RandomVariable>& random_variables();

  /// return ranVarTypes
  const ShortArray& random_variable_types() const;
  /// set ranVarTypes
  void random_variable_types(const ShortArray& rv_types);
  /// return ranVarTypes[i]
  short random_variable_type(size_t i) const;
  /// set ranVarTypes[i]
  void random_variable_type(short rv_type, size_t i);

  /// pull non-standardized distribution parameters from mv_dist
  void pull_distribution_parameters(const MultivariateDistribution& mv_dist);
  /// pull non-standardized distribution parameters from mv_dist
  void pull_distribution_parameters(MultivariateDistribution* mv_dist_rep);

  /// return activeVars
  const BitArray& active_variables() const;
  /// return activeCorr
  const BitArray& active_correlations() const;

  /// return corrMatrix
  const RealSymMatrix& correlation_matrix() const;
  /// set corrMatrix
  void correlation_matrix(const RealSymMatrix& corr);

  /// assemble means and standard deviations from RandomVariable::moments()
  RealRealPairArray moments() const;
  /// assemble means from RandomVariable::mean()
  RealVector means() const;
  /// assemble standard deviations from RandomVariable::standard_deviation()
  RealVector std_deviations() const;

  /// assemble lower and upper bounds from RandomVariable::bounds()
  RealRealPairArray bounds() const;
  /// assemble lower bounds from RandomVariable::bounds()
  RealVector lower_bounds() const;
  /// assemble upper bounds from RandomVariable::bounds()
  RealVector upper_bounds() const;

  //
  //- Heading: Member function definitions
  //

  /// set ranVarTypes and initialize randomVars
  void initialize_types(const ShortArray& rv_types,
			const BitArray& active_vars = BitArray());
  /// initializes corrMatrix and correlationFlag
  void initialize_correlations(const RealSymMatrix& corr,
			       const BitArray& active_corr = BitArray());

  /// update a scalar distribution parameter within randomVars[v]
  template <typename ValueType>
  void push_parameter(size_t v, short dist_param, ValueType value);
  /// update values for one distribution parameter across a sequence
  /// of random variables
  template <typename OrdinalType, typename ScalarType>
  void push_parameters(size_t start_v, size_t num_v, short dist_param,
    const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& values);
  /// update values for one distribution parameter across the set
  /// of random variables with matching RV type
  template <typename OrdinalType, typename ScalarType>
  void push_parameters(short rv_type, short dist_param,
    const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& values);
  /// update values for one distribution parameter across a sequence
  /// of random variables
  template <typename ValueType>
  void push_parameters(size_t start_v, size_t num_v, short dist_param,
		       const std::vector<ValueType>& values);
  /// update values for one distribution parameter across the set
  /// of random variables with matching RV type
  template <typename ValueType>
  void push_parameters(short rv_type, short dist_param,
		       const std::vector<ValueType>& values);

  /// update val from a scalar distribution parameter from randomVars[v]
  template <typename ValueType>
  void pull_parameter(size_t v, short dist_param, ValueType& val);
  /// define array of values using the identified distribution parameter
  /// across the specified range of random variables
  template <typename ValueType>
  void pull_parameters(size_t start_v, size_t num_v, short dist_param,
		       std::vector<ValueType>& values);
  /// define array of values using the identified distribution parameter
  /// across the specified range of random variables
  template <typename ValueType>
  void pull_parameters(short rv_type, short dist_param,
		       std::vector<ValueType>& values);
  /// update values for one distribution parameter across the set
  /// of random variables with matching RV type
  template <typename OrdinalType, typename ScalarType>
  void pull_parameters(short rv_type, short dist_param,
    Teuchos::SerialDenseVector<OrdinalType, ScalarType>& values);

  /// return a scalar distribution parameter from randomVars[v]
  template <typename ValueType>
  ValueType pull_parameter(size_t v, short dist_param);
  /// return array of values for one distribution parameter across the set
  /// of random variables with matching RV type
  template <typename ValueType>
  std::vector<ValueType> pull_parameters(short rv_type, short dist_param);

  /*
  /// expand corrMatrix from probabilistic variables to combined variables
  void expand_correlation_matrix(size_t num_lead_v, size_t num_prob_v,
				  size_t num_trail_v);
  /// contract corrMatrix from combined variables to probabilistic variables
  void contract_correlation_matrix(size_t num_lead_v, size_t num_prob_v,
				  size_t num_trail_v);
  */

  /// verify that randomVarsX[i].type() equals rv_type
  void check_random_variable_type(size_t i, short rv_type) const;

  /// return the univariate PDF value for a random variable
  Real pdf(Real val, size_t i) const;
  /// return the univariate log PDF value for a random variable
  Real log_pdf(Real val, size_t i) const;
  /// return the gradient of the univariate log PDF for a random variable
  Real log_pdf_gradient(Real val, size_t i) const;
  /// return the Hessian of the univariate log PDF for a random variable
  Real log_pdf_hessian(Real val, size_t i) const;

  /// draw a sample from the i-th RandomVariable
  template <typename Engine> 
  Real draw_sample(size_t i, Engine& rng) const;
  /// draw a sample from the i-th standardized RandomVariable
  template <typename Engine> 
  Real draw_standard_sample(size_t i, Engine& rng) const;

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// return the multivariate PDF value for full set of random variables
  Real pdf(const RealVector& pt) const;
  /// return the multivariate log PDF value for full set of random variables
  Real log_pdf(const RealVector& pt) const;

  /// copy marginals + correlation data between representations
  void copy_rep(MultivariateDistribution* mvd_rep);

  //
  //- Heading: Data
  //

  /// vector of types of each u-space standardized uncertain variable to
  /// which each x-space variable is transformed
  ShortArray ranVarTypes;
  /// vector of random variables encapsulating distribution parameters and
  /// statistical functions (pdf, cdf, etc.)
  std::vector<RandomVariable> randomVars;
  /// subset of randomVars that are currently active (if empty, then
  /// no subset: all variables are active)
  BitArray activeVars;

  /// matrix of random variable correlation coefficients
  RealSymMatrix corrMatrix;
  /// subset of randomVars to which the corrMatrix refers (if empty,
  /// then no subset: correlations are provided for all variables)
  BitArray activeCorr;

private:

  //
  //- Heading: Data
  //

};


inline MarginalsCorrDistribution::MarginalsCorrDistribution():
  MultivariateDistribution(BaseConstructor())
{ }


inline MarginalsCorrDistribution::~MarginalsCorrDistribution()
{ }


inline const RandomVariable& MarginalsCorrDistribution::
random_variable(size_t i) const
{ return randomVars[i]; }


inline const std::vector<RandomVariable>& MarginalsCorrDistribution::
random_variables() const
{ return randomVars; }


inline std::vector<RandomVariable>& MarginalsCorrDistribution::
random_variables()
{ return randomVars; }


inline const ShortArray& MarginalsCorrDistribution::
random_variable_types() const
{ return ranVarTypes; }


inline void MarginalsCorrDistribution::
random_variable_types(const ShortArray& rv_types)
{ ranVarTypes = rv_types; }


inline short MarginalsCorrDistribution::random_variable_type(size_t i) const
{ return ranVarTypes[i]; }


inline void MarginalsCorrDistribution::
random_variable_type(short rv_type, size_t i)
{ ranVarTypes[i] = rv_type; }


inline void MarginalsCorrDistribution::
check_random_variable_type(size_t i, short rv_type) const
{
  if (randomVars[i].type() != rv_type) {
    PCerr << "Error: inconsistent random variable type in MarginalsCorr"
	  << "Distribution::check_random_variable_type()." << std::endl;
    abort_handler(-1);
  }
}


inline const BitArray& MarginalsCorrDistribution::active_variables() const
{ return activeVars; }


inline const BitArray& MarginalsCorrDistribution::active_correlations() const
{ return activeCorr; }


inline const RealSymMatrix& MarginalsCorrDistribution::
correlation_matrix() const
{ return corrMatrix; }


inline void MarginalsCorrDistribution::
correlation_matrix(const RealSymMatrix& corr)
{ initialize_correlations(corr); } // Note: default active_corr = BitArray()


//inline const RealMatrix& MarginalsCorrDistribution::correlation_factor() const
//{ return corrCholeskyFactor; }


void MarginalsCorrDistribution::
pull_distribution_parameters(const MultivariateDistribution& mv_dist)
{ pull_distribution_parameters(mv_dist.multivar_dist_rep()); }


template <typename ValueType>
void MarginalsCorrDistribution::
push_parameter(size_t v, short dist_param, ValueType value)
{ randomVars[v].push_parameter(dist_param, value); }


template <typename OrdinalType, typename ScalarType>
void MarginalsCorrDistribution::
push_parameters(size_t start_v, size_t num_v, short dist_param,
  const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& values)
{
  // set one type of distribution parameter for a range of random variables
  // TO DO: would like to retire this version if Dakota migrates from Teuchos
  //        to std::vector for dist params

  size_t i, v, num_updates = std::min((size_t)values.length(), num_v);
  for (i=0, v=start_v; i<num_updates; ++i, ++v)
    randomVars[v].push_parameter(dist_param, values[i]);
}


template <typename ValueType>
void MarginalsCorrDistribution::
push_parameters(size_t start_v, size_t num_v, short dist_param,
		const std::vector<ValueType>& values)
{
  // set one distribution parameter type for a range of random variables

  size_t i, v, num_updates = std::min(values.size(), num_v);
  for (i=0, v=start_v; i<num_updates; ++i, ++v)
    randomVars[v].push_parameter(dist_param, values[i]);
}


template <typename OrdinalType, typename ScalarType>
void MarginalsCorrDistribution::
push_parameters(short rv_type, short dist_param,
  const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& values)
{
  // rv_type eliminates need to check for dist_param support
  // TO DO: would like to retire this version if Dakota migrates from Teuchos
  //        to std::vector for dist params

  size_t rv, num_rv = ranVarTypes.size(), cntr = 0, num_vals = values.length();
  for (rv=0; rv < num_rv && cntr < num_vals; ++rv)
    if (ranVarTypes[rv] == rv_type)
      randomVars[rv].push_parameter(dist_param, values[cntr++]);
}


template <typename ValueType>
void MarginalsCorrDistribution::
push_parameters(short rv_type, short dist_param,
		const std::vector<ValueType>& values)
{
  // rv_type eliminates need to check for dist_param support

  size_t rv, num_rv = ranVarTypes.size(), cntr = 0, num_vals = values.size();
  for (rv=0; rv < num_rv && cntr < num_vals; ++rv)
    if (ranVarTypes[rv] == rv_type)
      randomVars[rv].push_parameter(dist_param, values[cntr++]);
}


template <typename ValueType>
void MarginalsCorrDistribution::
pull_parameter(size_t v, short dist_param, ValueType& val)
{ randomVars[v].pull_parameter(dist_param, val); }


template <typename ValueType>
void MarginalsCorrDistribution::
pull_parameters(size_t start_v, size_t num_v, short dist_param,
		std::vector<ValueType>& values)
{
  // set one distribution parameter type for a range of random variables

  size_t i, v;
  values.resize(num_v);
  for (i=0, v=start_v; i<num_v; ++i, ++v)
    randomVars[v].pull_parameter(dist_param, values[i]);
}


template <typename ValueType>
void MarginalsCorrDistribution::
pull_parameters(short rv_type, short dist_param, std::vector<ValueType>& values)
{
  // rv_type eliminates need to check for dist_param support

  values.resize(std::count(ranVarTypes.begin(), ranVarTypes.end(), rv_type));

  size_t rv, num_rv = ranVarTypes.size(), cntr = 0;
  for (rv=0; rv<num_rv; ++rv)
    if (ranVarTypes[rv] == rv_type)
      randomVars[rv].pull_parameter(dist_param, values[cntr++]);
}


template <typename OrdinalType, typename ScalarType>
void MarginalsCorrDistribution::
pull_parameters(short rv_type, short dist_param,
		Teuchos::SerialDenseVector<OrdinalType, ScalarType>& values)
{
  // rv_type eliminates need to check for dist_param support

  size_t rv, num_rv = ranVarTypes.size(), cntr = 0;
  values.sizeUninitialized(
    std::count(ranVarTypes.begin(), ranVarTypes.end(), rv_type));

  for (rv=0; rv<num_rv; ++rv)
    if (ranVarTypes[rv] == rv_type)
      randomVars[rv].pull_parameter(dist_param, values[cntr++]);
}


template <typename ValueType>
ValueType MarginalsCorrDistribution::pull_parameter(size_t v, short dist_param)
{
  ValueType val;
  randomVars[v].pull_parameter(dist_param, val);
  return val;
}


template <typename ValueType>
std::vector<ValueType> MarginalsCorrDistribution::
pull_parameters(short rv_type, short dist_param)
{
  std::vector<ValueType> vals;
  pull_parameters(rv_type, dist_param, vals);
  return vals;
}


/* These APIs are not currently needed but could be useful in the future:
template <typename VectorType>
void MarginalsCorrDistribution::
push_parameters(size_t v, const ShortArray& dist_params, const VectorType& values)
{
  // set multiple distribution parameters for a single variable

  RandomVariable& random_var = randomVars[v];
  size_t i, num_params = std::min(dist_params.size(), values.length());
  for (i=0; i<num_params; ++i)
    random_var.push_parameter(dist_params[i], values[i]);
}


template <typename VectorType>
void MarginalsCorrDistribution::
push_parameters(short dist_param, const VectorType& values)
{
  // Without rv_type, query RV for dist_param support to match values to RV

  size_t rv, num_rv = randomVars.size(), cntr = 0, num_vals = values.length();
  for (rv=0; rv < num_rv && cntr < num_vals; ++rv)
    if (randomVars[rv].supports(dist_param))
      randomVars[rv].push_parameter(dist_param, values[cntr++]);
}
*/


inline RealRealPairArray MarginalsCorrDistribution::moments() const
{
  size_t i, num_v = randomVars.size();
  RealRealPairArray mom(num_v);
  for (i=0; i<num_v; ++i)
    mom[i] = randomVars[i].moments();
  return mom;
}


inline RealVector MarginalsCorrDistribution::means() const
{
  size_t i, num_v = randomVars.size();
  RealVector means(num_v, false);
  for (i=0; i<num_v; ++i)
    means[i] = randomVars[i].mean();
  return means;
}


inline RealVector MarginalsCorrDistribution::std_deviations() const
{
  size_t i, num_v = randomVars.size();
  RealVector std_devs(num_v, false);
  for (i=0; i<num_v; ++i)
    std_devs[i] = randomVars[i].standard_deviation();
  return std_devs;
}


inline RealRealPairArray MarginalsCorrDistribution::bounds() const
{
  size_t i, num_v = randomVars.size();
  RealRealPairArray bnds(num_v);
  for (i=0; i<num_v; ++i)
    bnds[i] = randomVars[i].bounds();
  return bnds;
}


inline RealVector MarginalsCorrDistribution::lower_bounds() const
{
  size_t i, num_v = randomVars.size();
  RealVector lwr_bnds(num_v, false);
  for (i=0; i<num_v; ++i)
    lwr_bnds[i] = randomVars[i].bounds().first;
  return lwr_bnds;
}


inline RealVector MarginalsCorrDistribution::upper_bounds() const
{
  size_t i, num_v = randomVars.size();
  RealVector upr_bnds(num_v, false);
  for (i=0; i<num_v; ++i)
    upr_bnds[i] = randomVars[i].bounds().second;
  return upr_bnds;
}


inline Real MarginalsCorrDistribution::pdf(Real val, size_t i) const
{ return randomVars[i].pdf(val); }


inline Real MarginalsCorrDistribution::log_pdf(Real val, size_t i) const
{ return randomVars[i].log_pdf(val); }


inline Real MarginalsCorrDistribution::
log_pdf_gradient(Real val, size_t i) const
{ return randomVars[i].log_pdf_gradient(val); }


inline Real MarginalsCorrDistribution::
log_pdf_hessian(Real val, size_t i) const
{ return randomVars[i].log_pdf_hessian(val); }


inline Real MarginalsCorrDistribution::pdf(const RealVector& pt) const
{
  // TO DO: add support for evaluation of correlated MVN density
  if (correlationFlag) {
    PCerr << "Error: MarginalsCorrDistribution::pdf() currently uses a "
	  << "product of marginal densities\n       and can only be used for "
	  << "independent random variables." << std::endl;
    abort_handler(-1);
  }
  size_t i, num_v = randomVars.size();
  Real density = 1.;
  for (i=0; i<num_v; ++i)
    density *= pdf(pt[i], i);
  return density;
}


inline Real MarginalsCorrDistribution::log_pdf(const RealVector& pt) const
{
  // TO DO: add support for evaluation of correlated MVN density
  if (correlationFlag) {
    PCerr << "Error: MarginalsCorrDistribution::log_pdf() currently uses a "
	  << "sum of log marginal densities\n       and can only be used for "
	  << "independent random variables." << std::endl;
    abort_handler(-1);
  }
  size_t i, num_v = randomVars.size();
  Real log_density = 0.;
  for (i=0; i<num_v; ++i)
    log_density += log_pdf(pt[i], i);
  return log_density;
}


template <typename Engine> 
Real MarginalsCorrDistribution::draw_sample(size_t i, Engine& rng) const
{ return randomVars[i].draw_sample(rng); }


template <typename Engine> 
Real MarginalsCorrDistribution::
draw_standard_sample(size_t i, Engine& rng) const
{
  /*
  // can only use randomVarsX[i].standard_pdf() for cases where u_type is a
  // standardized form of the x_type.  For STD_NORMAL and STD_UNIFORM, many
  // x_types can be mapped to these u_types, so use global utility fns
  // whenever there are no auxilliary parameters to manage.
  switch (ranVarTypesU[i]) {
  // these cases require static fns since U type may not correspond to X type
  case STD_NORMAL:  return  NormalRandomVariable::draw_std_sample(rng); break;
  case STD_UNIFORM: return UniformRandomVariable::draw_std_sample(rng); break;
  case STD_EXPONENTIAL:
    return ExponentialRandomVariable::draw_std_sample(rng);             break;
  // these cases can rely on correspondence between X and U types
  case STD_BETA:
    check_x_type(i, BETA);
    return randomVarsX[i].draw_standard_sample(rng); break;
  case STD_GAMMA:
    check_x_type(i, GAMMA);
    return randomVarsX[i].draw_standard_sample(rng); break;
  default: // no transformation (e.g., PCE with numerically-generated bases)
    check_x_type(i, ranVarTypesU[i]);
    return randomVarsX[i].draw_sample(rng);          break;
  }
  */

  return randomVars[i].draw_standard_sample(rng);
}

} // namespace Pecos

#endif

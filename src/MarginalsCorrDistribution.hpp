/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef MARGINALS_CORR_DISTRIBUTION_HPP
#define MARGINALS_CORR_DISTRIBUTION_HPP

#include "MultivariateDistribution.hpp"


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

  
  /// initializes randomVarsX (no transformation: u-space not needed)
  void initialize_random_variables(const ShortArray& rv_types);
  /// updates parameters within randomVarsX
  void initialize_random_variable_parameters(const RealVector& cd_l_bnds,
					     const RealVector& cd_u_bnds,
					     const AleatoryDistParams& adp,
					     const EpistemicDistParams& edp,
					     const RealVector& cs_l_bnds,
					     const RealVector& cs_u_bnds);
  /// initializes corrMatrixX and correlationFlagX
  void initialize_correlations(const RealSymMatrix& x_corr);
  /// reshape corrMatrixX for an all_variables specification
  void reshape_correlation_matrix(size_t num_lead_v, size_t num_prob_v,
				  size_t num_trail_v);

  /// return randomVars
  const std::vector<RandomVariable>& random_variables() const;

  /// return ranVarTypes
  const ShortArray& types() const;
  /// set ranVarTypes
  void types(const ShortArray& rv_types);
  /// set ranVarTypes[i]
  void type(short rv_type, size_t i);
  /// verify that randomVarsX[i].type() equals rv_type
  void check_type(size_t i, short rv_type) const;

  /// return corrMatrix
  const RealSymMatrix& correlation_matrix() const;
  /// return corrCholeskyFactor
  const RealMatrix& correlation_factor() const;

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

  /// return the univariate PDF value for a random variable
  Real pdf(Real val, size_t i) const;
  /// return the univariate log PDF value for a random variable
  Real log_pdf(Real val, size_t i) const;
  /// return the gradient of the univariate log PDF for a random variable
  Real log_pdf_gradient(Real val, size_t i) const;
  /// return the Hessian of the univariate log PDF for a random variable
  Real log_pdf_hessian(Real val, size_t i) const;

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// return the multivariate PDF value for x-space random variables
  Real pdf(const RealVector& pt) const;
  /// return the multivariate log PDF value for x-space random variables
  Real log_pdf(const RealVector& pt) const;

  //
  //- Heading: Data
  //

  /// vector of random variables encapsulating distribution parameters and
  /// statistical functions (pdf, cdf, etc.)
  std::vector<RandomVariable> randomVars;
  /// vector of types of each u-space standardized uncertain variable to
  /// which each x-space variable is transformed
  ShortArray ranVarTypes;

  /// matrix of random variable correlation coefficients
  RealSymMatrix corrMatrix;
  // cholesky factor of a modified correlation matrix (#corrMatrixX
  // is modified in transform_correlations() for use in z-space)
  //RealMatrix corrCholeskyFactor;

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


inline const std::vector<RandomVariable>& MarginalsCorrDistribution::
random_variables() const
{ return randomVars; }


inline const ShortArray& MarginalsCorrDistribution::types() const
{ return ranVarTypes; }


inline void MarginalsCorrDistribution::types(const ShortArray& rv_types)
{ ranVarTypes = rv_types; }


inline void MarginalsCorrDistribution::type(short rv_type, size_t i)
{ ranVarTypes[i] = rv_type; }


inline void MarginalsCorrDistribution::check_type(size_t i, short rv_type) const
{
  if (randomVars[i].type() != rv_type) {
    PCerr << "Error: failure in random variable type check in "
	  << "MarginalsCorrDistribution." << std::endl;
    abort_handler(-1);
  }
}


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
{
  return (mvDistRep) ? mvDistRep->randomVars[i].log_pdf_gradient(val) :
    randomVars[i].log_pdf_gradient(val);
}


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


inline const RealSymMatrix& MarginalsCorrDistribution::
correlation_matrix() const
{ return corrMatrix; }


inline const RealMatrix& MarginalsCorrDistribution::correlation_factor() const
{ return corrCholeskyFactor; }


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

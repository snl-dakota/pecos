/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "NatafTransformation.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: NatafTransformation.cpp 4768 2007-12-17 17:49:32Z mseldre $";

//#define DEBUG


namespace Pecos {


/** This procedure performs the transformation from u to x space.
    u_vars is the vector of random variables in uncorrelated standard
    normal space (u-space).  x_vars is the vector of random variables
    in the original user-defined x-space. */
void NatafTransformation::
trans_U_to_X(const RealVector& u_vars, RealVector& x_vars)
{
  if (correlationFlagX) {
    RealVector z_vars;
    trans_U_to_Z(u_vars, z_vars);
    trans_Z_to_X(z_vars, x_vars);
  }
  else // z_vars = u_vars
    trans_Z_to_X(u_vars, x_vars);

#ifdef DEBUG
  PCout << "Transforming u_vars:\n" << u_vars << "to x_vars:\n" << x_vars;
#endif
}


/** This procedure computes the transformation from u to z space.
    u_vars is the vector of random variables in uncorrelated standard
    normal space (u-space).  z_vars is the vector of random variables
    in normal space with proper correlations (z-space). */
void NatafTransformation::
trans_U_to_Z(const RealVector& u_vars, RealVector& z_vars)
{
  // corrCholeskyFactorZ: the Cholesky factor of the modified correlation matrix

  int u_len = u_vars.length();
  if (z_vars.length() != u_len)
    z_vars.size(u_len);
  z_vars.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1., corrCholeskyFactorZ,
		  u_vars, 0.);
}


/** This procedure computes the transformation from z to x space.
    z_vars is the vector of random variables in normal space with
    proper correlations (z-space).  x_vars is the vector of random
    variables in the original user-defined x-space */
void NatafTransformation::
trans_Z_to_X(const RealVector& z_vars, RealVector& x_vars)
{
  // ranVarTypesX:   the type of random variable (NORMAL/BOUNDED_NORMAL/
  //                 LOGNORMAL/BOUNDED_LOGNORMAL/UNIFORM/LOGUNIFORM/TRIANGULAR/
  //                 EXPONENTIAL/BETA/GAMMA/GUMBEL/FRECHET/WEIBULL)
  // ranVarMeansX:   means of the random variables in x-space
  // ranVarStdDevsX: standard deviations of the random variables in x-space

  // This routine performs an inverse Rackwitz-Fiessler transformation:
  //   X = F^{-1}[ Phi(Z) ]

  int z_len = z_vars.length();
  if (x_vars.length() != z_len)
    x_vars.size(z_len);
  Real dbl_inf = std::numeric_limits<Real>::infinity();

  for (int i=0; i<z_len; i++) {
    bool err_flag = false;
    switch (ranVarTypesX[i]) {
    case CONTINUOUS_DESIGN: case CONTINUOUS_STATE: case CONTINUOUS_INTERVAL:
      if (ranVarTypesU[i] == STD_UNIFORM) {
	// scale from [-1,1] to [L,U]
	const Real& lwr = ranVarLowerBndsX[i];
	x_vars[i] = lwr + (ranVarUpperBndsX[i]-lwr)*std_uniform_cdf(z_vars[i]);
      }
      else
	err_flag = true;
      break;
    case NORMAL: // unbounded normal: z = (x - mean)/stdev
      if (ranVarTypesU[i] == STD_NORMAL)
	x_vars[i] = z_vars[i]*ranVarStdDevsX[i] + ranVarMeansX[i];
      else
	err_flag = true;
      break;
    case BOUNDED_NORMAL: // bounded normal
      if (ranVarTypesU[i] == STD_NORMAL) {
	const Real& lwr = ranVarLowerBndsX[i];
	const Real& upr = ranVarUpperBndsX[i];
	const Real& z   = z_vars[i];
	if      (z ==  dbl_inf) x_vars[i] = upr;
	else if (z == -dbl_inf) x_vars[i] = lwr;
	else {
	  const Real& mean  = ranVarMeansX[i];
	  const Real& stdev = ranVarStdDevsX[i];
	  Real Phi_lms = (lwr > -dbl_inf) ? Phi((lwr-mean)/stdev) : 0.;
	  Real Phi_ums = (upr <  dbl_inf) ? Phi((upr-mean)/stdev) : 1.;
	  x_vars[i]
	    = Phi_inverse(Phi(z)*(Phi_ums - Phi_lms) + Phi_lms) * stdev + mean;
	}
      }
      else if (ranVarTypesU[i] == STD_UNIFORM) {
	const Real& lwr = ranVarLowerBndsX[i];
	const Real& upr = ranVarUpperBndsX[i];
	const Real& z   = z_vars[i];
	if      (z ==  1.) x_vars[i] = upr;
	else if (z == -1.) x_vars[i] = lwr;
	else {
	  const Real& mean  = ranVarMeansX[i];
	  const Real& stdev = ranVarStdDevsX[i];
	  Real Phi_lms = (lwr > -dbl_inf) ? Phi((lwr-mean)/stdev) : 0.;
	  Real Phi_ums = (upr <  dbl_inf) ? Phi((upr-mean)/stdev) : 1.;
	  x_vars[i] = Phi_inverse(std_uniform_cdf(z_vars[i])*
				  (Phi_ums - Phi_lms) + Phi_lms) * stdev + mean;
	}
      }
      else if (ranVarTypesU[i] == BOUNDED_NORMAL) // Golub-Welsch: no transform
	x_vars[i] = z_vars[i];
      else
	err_flag = true;
      break;
    case LOGNORMAL: // unbounded lognormal: z = (ln x - lambda)/zeta
      if (ranVarTypesU[i] == STD_NORMAL)
	x_vars[i] = std::exp(ranVarAddtlParamsX[i][0] +
			ranVarAddtlParamsX[i][1] * z_vars[i]);
      else if (ranVarTypesU[i] == LOGNORMAL) // Golub-Welsch: no transform
	x_vars[i] = z_vars[i];
      else
	err_flag = true;
      break;
    case BOUNDED_LOGNORMAL: // bounded lognormal
      if (ranVarTypesU[i] == STD_NORMAL) {
	const Real& lwr = ranVarLowerBndsX[i];
	const Real& upr = ranVarUpperBndsX[i];
	const Real& z   = z_vars[i];
	if      (z ==  dbl_inf) x_vars[i] = upr;
	else if (z == -dbl_inf) x_vars[i] = lwr;
	else {
	  const Real& lambda = ranVarAddtlParamsX[i][0];
	  const Real& zeta   = ranVarAddtlParamsX[i][1];
	  Real Phi_lms = (lwr > 0.)      ? Phi((log(lwr)-lambda)/zeta) : 0.;
	  Real Phi_ums = (upr < dbl_inf) ? Phi((log(upr)-lambda)/zeta) : 1.;
	  x_vars[i] = std::exp(Phi_inverse(Phi(z)*(Phi_ums - Phi_lms) + Phi_lms)
		    * zeta + lambda);
	}
      }
      else if (ranVarTypesU[i] == STD_UNIFORM) {
	const Real& lwr = ranVarLowerBndsX[i];
	const Real& upr = ranVarUpperBndsX[i];
	const Real& z   = z_vars[i];
	if      (z ==  1.) x_vars[i] = upr;
	else if (z == -1.) x_vars[i] = lwr;
	else {
	  const Real& lambda = ranVarAddtlParamsX[i][0];
	  const Real& zeta   = ranVarAddtlParamsX[i][1];
	  Real Phi_lms = (lwr > 0.)      ? Phi((log(lwr)-lambda)/zeta) : 0.;
	  Real Phi_ums = (upr < dbl_inf) ? Phi((log(upr)-lambda)/zeta) : 1.;
	  x_vars[i] = std::exp(Phi_inverse(std_uniform_cdf(z_vars[i])
		    * (Phi_ums-Phi_lms) + Phi_lms) * zeta + lambda);
	}
      }
      else if (ranVarTypesU[i] == BOUNDED_LOGNORMAL) // Golub-Welsch: no x-form
	x_vars[i] = z_vars[i];
      else
	err_flag = true;
      break;
    case UNIFORM: {
      const Real& lwr = ranVarLowerBndsX[i];
      Real range = ranVarUpperBndsX[i] - lwr;
      if (ranVarTypesU[i] == STD_UNIFORM) // scale from [-1,1] to [L,U]
	x_vars[i] = lwr + range*(z_vars[i]+1.)/2.;
      else if (ranVarTypesU[i] == STD_NORMAL) { // transform from std normal
	// Phi(z) = (x-L)/(U-L)
	Real normcdf = Phi(z_vars[i]);
	x_vars[i] = normcdf * range + lwr;
      }
      else
	err_flag = true;
      break;
    }
    case LOGUNIFORM: {
      const Real& lwr = ranVarLowerBndsX[i];
      Real log_lwr = log(lwr), log_range = log(ranVarUpperBndsX[i]) - log_lwr;
      if (ranVarTypesU[i] == STD_UNIFORM) // transform from uniform on [-1,1]
	x_vars[i] = lwr*std::exp((z_vars[i]+1.)*log_range/2.);
      else if (ranVarTypesU[i] == STD_NORMAL) { // transform from std normal
	// Phi(z) = (ln x - ln L)/(ln U - ln L)
	Real normcdf = Phi(z_vars[i]);
	x_vars[i] = std::exp(normcdf * log_range + log_lwr);
      }
      else if (ranVarTypesU[i] == LOGUNIFORM) // Golub-Welsch: no transform
	x_vars[i] = z_vars[i];
      else
	err_flag = true;
      break;
    }
    case TRIANGULAR:
      if (ranVarTypesU[i] == TRIANGULAR) // Golub-Welsch: no transform
	x_vars[i] = z_vars[i];
      else {
	const Real& lwr  = ranVarLowerBndsX[i];
	const Real& mode = ranVarAddtlParamsX[i][0];
	const Real& upr  = ranVarUpperBndsX[i];
	const Real& z    = z_vars[i];
	Real range = upr - lwr;
	Real zcdf;
	if (ranVarTypesU[i] == STD_UNIFORM)
	  zcdf = (z+1.)/2.;
	else if (ranVarTypesU[i] == STD_NORMAL)
	  zcdf = Phi(z);
	else
	  err_flag = true;
	// assume x < mode and then check
	Real x = lwr + std::sqrt(zcdf*range*(mode-lwr));
	Real x_pdf = 2.*(x-lwr)/range/(mode-lwr), m_pdf = 2./range;
	// check pdf value to ensure that proper equation used
	if ( x_pdf > m_pdf ) {
	  // avoid numerical probs for STD_NORMAL w/ large z > 0
	  // (zcdf indistinguishable from 1)
	  Real zcdf_comp
	    = (ranVarTypesU[i] == STD_NORMAL && z > 0.) ? Phi(-z) : 1. - zcdf;
	  x = upr - std::sqrt(zcdf_comp*range*(upr-mode));
	}
	x_vars[i] = x;
      }
      break;
    case EXPONENTIAL:
      if (ranVarTypesU[i] == STD_EXPONENTIAL) // scale: exp(-x/beta) = exp(-z)
	x_vars[i] = ranVarAddtlParamsX[i][0]*z_vars[i];
      else if (ranVarTypesU[i] == STD_NORMAL) { // transform from std normal
	const Real& z = z_vars[i];
	Real log1mnormcdf = (z > 0.) ? std::log(Phi(-z)) : bmth::log1p(-Phi(z));
	x_vars[i] = -ranVarAddtlParamsX[i][0]*log1mnormcdf;
      }
      else
	err_flag = true;
      break;
    case BETA: { // scale from std beta on [-1,1]
      const Real& lwr = ranVarLowerBndsX[i];
      const Real& upr = ranVarUpperBndsX[i];
      if (ranVarTypesU[i] == STD_BETA) // scale from [-1,1] to [L,U]
	x_vars[i] = lwr + (upr - lwr)*(z_vars[i]+1.)/2.;
      else if (ranVarTypesU[i] == STD_NORMAL) { // transform from std normal
	Real normcdf = Phi(z_vars[i]), scaled_x = std_beta_cdf_inverse(normcdf,
	  ranVarAddtlParamsX[i][0], ranVarAddtlParamsX[i][1]);
	x_vars[i] = lwr + (upr - lwr)*scaled_x;
      }
      else if (ranVarTypesU[i] == STD_UNIFORM) {
	Real scaled_x = std_beta_cdf_inverse(std_uniform_cdf(z_vars[i]),
	                ranVarAddtlParamsX[i][0], ranVarAddtlParamsX[i][1]);
	x_vars[i] = lwr + (upr - lwr)*scaled_x;
      }
      else
	err_flag = true;
      break;
    }
    case GAMMA: {
      if (ranVarTypesU[i] == STD_GAMMA)
	x_vars[i] = ranVarAddtlParamsX[i][1]*z_vars[i]; // x/beta = z
      else if (ranVarTypesU[i] == STD_NORMAL)
	x_vars[i] = gamma_cdf_inverse(Phi(z_vars[i]), ranVarAddtlParamsX[i][0],
				      ranVarAddtlParamsX[i][1]);
      else
	err_flag = true;
      break;
    }
    case GUMBEL:
      if (ranVarTypesU[i] == STD_NORMAL) {
	// Phi(z) = F(x) = e^(-e^(-alpha(x-beta)))
	const Real& z = z_vars[i];
	// avoid numerical problems for large z > 0
	// (normcdf indistinguishable from 1):
	Real lognormcdf = (z > 0.) ? bmth::log1p(-Phi(-z)) : std::log(Phi(z));
	x_vars[i] = ranVarAddtlParamsX[i][1]
	  - std::log(-lognormcdf)/ranVarAddtlParamsX[i][0];
      }
      else if (ranVarTypesU[i] == GUMBEL) // Golub-Welsch: no transform
	x_vars[i] = z_vars[i];
      else
	err_flag = true;
      break;
    case FRECHET:
      if (ranVarTypesU[i] == STD_NORMAL) {
	// Phi(z) = F(x) = e^(-(beta/x)^alpha)
	const Real& z = z_vars[i];
	// avoid numerical problems for large z > 0
	// (normcdf indistinguishable from 1):
	Real lognormcdf = (z > 0.) ? bmth::log1p(-Phi(-z)) : std::log(Phi(z));
	x_vars[i] = ranVarAddtlParamsX[i][1] *
	  std::pow(-lognormcdf, -1./ranVarAddtlParamsX[i][0]);
      }
      else if (ranVarTypesU[i] == FRECHET) // Golub-Welsch: no transform
	x_vars[i] = z_vars[i];
      else
	err_flag = true;
      break;
    case WEIBULL:
      if (ranVarTypesU[i] == STD_NORMAL) {
	// Phi(z) = F(x) = 1 - e^(-(x/beta)^alpha)
	const Real& z = z_vars[i];
	// avoid numerical problems for large z > 0
	// (normcdf indistinguishable from 1):
	Real log1mnormcdf = (z > 0.) ? std::log(Phi(-z)) : bmth::log1p(-Phi(z));
	x_vars[i] = ranVarAddtlParamsX[i][1] *
	  std::pow(-log1mnormcdf, 1./ranVarAddtlParamsX[i][0]);
      }
      else if (ranVarTypesU[i] == WEIBULL) // Golub-Welsch: no transform
	x_vars[i] = z_vars[i];
      else
	err_flag = true;
      break;
    case HISTOGRAM_BIN:
      if (ranVarTypesU[i] == STD_NORMAL)
	x_vars[i] = histogram_bin_cdf_inverse(Phi(z_vars[i]),
					      ranVarAddtlParamsX[i]);
      else if (ranVarTypesU[i] == STD_UNIFORM)
	x_vars[i] = histogram_bin_cdf_inverse(std_uniform_cdf(z_vars[i]),
					      ranVarAddtlParamsX[i]);
      else if (ranVarTypesU[i] == HISTOGRAM_BIN)
	x_vars[i] = z_vars[i];
      else
	err_flag = true;
      break;
    }
    if (err_flag) {
      PCerr << "Error: unsupported variable mapping for variable " << i
	    << " in NatafTransformation::trans_Z_to_X()" << std::endl;
      abort_handler(-1);
    }
#ifdef DEBUG
    PCout << "Z_to_X: z[" << i << "] = " << z_vars[i]
	  <<    " --> x[" << i << "] = " << x_vars[i] << std::endl;
#endif // DEBUG
  }
}


/** This procedure performs the transformation from x to u space
    u_vars is the vector of random variables in uncorrelated standard
    normal space (u-space).  x_vars is the vector of random variables
    in the original user-defined x-space. */
void NatafTransformation::
trans_X_to_U(const RealVector& x_vars, RealVector& u_vars) 
{ 
  if (correlationFlagX) {
    RealVector z_vars;
    trans_X_to_Z(x_vars, z_vars);
    trans_Z_to_U(z_vars, u_vars);
  }
  else // z_vars = u_vars
    trans_X_to_Z(x_vars, u_vars);

#ifdef DEBUG
  PCout << "Transforming x_vars:\n" << x_vars << "to u_vars:\n" << u_vars;
#endif
}


/** This procedure performs the transformation from x to z space:
    z_vars is the vector of random variables in normal space with
    proper correlations (z-space).  x_vars is the vector of random
    variables in the original user-defined x-space. */
void NatafTransformation::
trans_X_to_Z(const RealVector& x_vars, RealVector& z_vars)
{
  // ranVarTypesX:   the type of random variable (NORMAL/BOUNDED_NORMAL/
  //                 LOGNORMAL/BOUNDED_LOGNORMAL/UNIFORM/LOGUNIFORM/TRIANGULAR/
  //                 EXPONENTIAL/BETA/GAMMA/GUMBEL/FRECHET/WEIBULL)
  // ranVarMeansX:   means of the random variables in x-space
  // ranVarStdDevsX: standard deviations of the random variables in x-space

  // This routine performs a Rackwitz-Fiessler transformation:
  //   Z = Phi^{-1}[ F(X) ]

  int x_len = x_vars.length();
  if (z_vars.length() != x_len)
    z_vars.size(x_len);
  Real dbl_inf = std::numeric_limits<Real>::infinity();

  for (int i=0; i<x_len; i++) {
    bool err_flag = false;
    switch (ranVarTypesX[i]) {
    case CONTINUOUS_DESIGN: case CONTINUOUS_STATE: case CONTINUOUS_INTERVAL:
      if (ranVarTypesU[i] == STD_UNIFORM) {
	const Real& lwr = ranVarLowerBndsX[i];
	z_vars[i] = std_uniform_cdf_inverse((x_vars[i] - lwr)/
					    (ranVarUpperBndsX[i] - lwr));
      }
      else
	err_flag = true;
      break;
    case NORMAL: // unbounded normal: z = (x - mean)/stdev
      if (ranVarTypesU[i] == STD_NORMAL)
	z_vars[i] = (x_vars[i]-ranVarMeansX[i])/ranVarStdDevsX[i];
      else
	err_flag = true;
      break;
    case BOUNDED_NORMAL: // bounded normal
      if (ranVarTypesU[i] == STD_NORMAL) {
	const Real& lwr = ranVarLowerBndsX[i];
	const Real& upr = ranVarUpperBndsX[i];
	const Real& x   = x_vars[i];
	if      (x <= lwr) z_vars[i] = -dbl_inf;
	else if (x >= upr) z_vars[i] =  dbl_inf;
	else
	  z_vars[i] = Phi_inverse(bounded_normal_cdf(x, ranVarMeansX[i],
	    ranVarStdDevsX[i], lwr, upr));
      }
      else if (ranVarTypesU[i] == STD_UNIFORM) {
	const Real& lwr = ranVarLowerBndsX[i];
	const Real& upr = ranVarUpperBndsX[i];
	const Real& x   = x_vars[i];
	if      (x <= lwr) z_vars[i] = -1.;
	else if (x >= upr) z_vars[i] =  1.;
	else
	  z_vars[i] = std_uniform_cdf_inverse(bounded_normal_cdf(x,
	    ranVarMeansX[i], ranVarStdDevsX[i], lwr, upr));
      }
      else if (ranVarTypesU[i] == BOUNDED_NORMAL) // Golub-Welsch: no transform
	z_vars[i] = x_vars[i];
      else
	err_flag = true;
      break;
    case LOGNORMAL: // unbounded lognormal: z = (ln x - lambda)/zeta
      if (ranVarTypesU[i] == STD_NORMAL)
	z_vars[i] = (log(x_vars[i]) - ranVarAddtlParamsX[i][0])
	  / ranVarAddtlParamsX[i][1];
      else if (ranVarTypesU[i] == LOGNORMAL) // Golub-Welsch: no transform
	z_vars[i] = x_vars[i];
      else
	err_flag = true;
      break;
    case BOUNDED_LOGNORMAL: // bounded lognormal
      if (ranVarTypesU[i] == STD_NORMAL) {
	const Real& lwr = ranVarLowerBndsX[i];
	const Real& upr = ranVarUpperBndsX[i];
	const Real& x   = x_vars[i];
	if      (x <= lwr) z_vars[i] = -dbl_inf;
	else if (x >= upr) z_vars[i] =  dbl_inf;
	else
	  z_vars[i] = Phi_inverse(bounded_lognormal_cdf(x, ranVarMeansX[i],
	    ranVarStdDevsX[i], lwr, upr));
      }
      else if (ranVarTypesU[i] == STD_UNIFORM) {
	const Real& lwr = ranVarLowerBndsX[i];
	const Real& upr = ranVarUpperBndsX[i];
	const Real& x   = x_vars[i];
	if      (x <= lwr) z_vars[i] = -1.;
	else if (x >= upr) z_vars[i] =  1.;
	else
	  z_vars[i] = std_uniform_cdf_inverse(bounded_lognormal_cdf(x,
	    ranVarMeansX[i], ranVarStdDevsX[i], lwr, upr));
      }
      else if (ranVarTypesU[i] == BOUNDED_LOGNORMAL) // Golub-Welsch: no x-form
	z_vars[i] = x_vars[i];
      else
	err_flag = true;
      break;
    case UNIFORM: {
      // Phi(z) = (x-L)/(U-L)
      if (ranVarTypesU[i] == STD_UNIFORM) { // scale to uniform on [-1,1]
	const Real& lwr = ranVarLowerBndsX[i];
	z_vars[i] = std_uniform_cdf_inverse((x_vars[i] - lwr)/
					    (ranVarUpperBndsX[i] - lwr));
      }
      else if (ranVarTypesU[i] == STD_NORMAL) // transform to std normal
	z_vars[i] = Phi_inverse(uniform_cdf(x_vars[i], ranVarLowerBndsX[i],
					    ranVarUpperBndsX[i]));
      else
	err_flag = true;
      break;
    }
    case LOGUNIFORM: {
      if (ranVarTypesU[i] == STD_UNIFORM) { // transform to uniform on [-1,1]
	Real log_lwr = log(ranVarLowerBndsX[i]),
	     log_upr = log(ranVarUpperBndsX[i]);
	z_vars[i] = std_uniform_cdf_inverse((log(x_vars[i]) - log_lwr)/
					    (log_upr - log_lwr));
      }
      else if (ranVarTypesU[i] == STD_NORMAL) // transform to std normal
	// Phi(z) = (ln x - ln L)/(ln U - ln L)
	z_vars[i] = Phi_inverse(loguniform_cdf(x_vars[i], ranVarLowerBndsX[i],
					       ranVarUpperBndsX[i]));
      else if (ranVarTypesU[i] == LOGUNIFORM) // Golub-Welsch: no transform
	z_vars[i] = x_vars[i];
      else
	err_flag = true;
      break;
    }
    case TRIANGULAR:
      if (ranVarTypesU[i] == TRIANGULAR) // Golub-Welsch: no transform
	z_vars[i] = x_vars[i];
      else {
	Real cdf = triangular_cdf(x_vars[i], ranVarAddtlParamsX[i][0],
				  ranVarLowerBndsX[i], ranVarUpperBndsX[i]);
	if (ranVarTypesU[i] == STD_UNIFORM)
	  z_vars[i] = std_uniform_cdf_inverse(cdf);
	else if (ranVarTypesU[i] == STD_NORMAL)
	  z_vars[i] = Phi_inverse(cdf);
	else
	  err_flag = true;
      }
      break;
    case EXPONENTIAL:
      if (ranVarTypesU[i] == STD_EXPONENTIAL)// scale to std exponential exp(-z)
	z_vars[i] = x_vars[i]/ranVarAddtlParamsX[i][0];
      else if (ranVarTypesU[i] == STD_NORMAL) // transform to std normal
	z_vars[i]
	  = Phi_inverse(exponential_cdf(x_vars[i], ranVarAddtlParamsX[i][0]));
      else
	err_flag = true;
      break;
    case BETA: {
      const Real& lwr = ranVarLowerBndsX[i];
      const Real& upr = ranVarUpperBndsX[i];
      if (ranVarTypesU[i] == STD_BETA) // scale to beta on [-1,1]
	z_vars[i] = 2.*(x_vars[i] - lwr)/(upr - lwr) - 1.;
      else if (ranVarTypesU[i] == STD_NORMAL) // transform to std normal
	z_vars[i] = Phi_inverse(beta_cdf(x_vars[i], ranVarAddtlParamsX[i][0],
					 ranVarAddtlParamsX[i][1], lwr, upr));
      else if (ranVarTypesU[i] == STD_UNIFORM)
	z_vars[i] = std_uniform_cdf_inverse(beta_cdf(x_vars[i],
	  ranVarAddtlParamsX[i][0], ranVarAddtlParamsX[i][1], lwr, upr));
      else
	err_flag = true;
      break;
    }
    case GAMMA:
      if (ranVarTypesU[i] == STD_GAMMA) // scale to std gamma
	z_vars[i] = x_vars[i] / ranVarAddtlParamsX[i][1];
      else if (ranVarTypesU[i] == STD_NORMAL) // transform to std normal
	z_vars[i] = Phi_inverse(gamma_cdf(x_vars[i], ranVarAddtlParamsX[i][0],
	  ranVarAddtlParamsX[i][1]));
      else
	err_flag = true;
      break;
    case GUMBEL:
      if (ranVarTypesU[i] == STD_NORMAL)
	// Phi(z) = F(x) = e^(-e^(-alpha(x-beta)))
	z_vars[i] = Phi_inverse(gumbel_cdf(x_vars[i], ranVarAddtlParamsX[i][0],
					   ranVarAddtlParamsX[i][1]));
      else if (ranVarTypesU[i] == GUMBEL) // Golub-Welsch: no transform
	z_vars[i] = x_vars[i];
      else
	err_flag = true;
      break;
    case FRECHET:
      if (ranVarTypesU[i] == STD_NORMAL) // Phi(z) = F(x) = e^(-(beta/x)^alpha)
	z_vars[i] = Phi_inverse(frechet_cdf(x_vars[i], ranVarAddtlParamsX[i][0],
					    ranVarAddtlParamsX[i][1]));
      else if (ranVarTypesU[i] == FRECHET) // Golub-Welsch: no transform
	z_vars[i] = x_vars[i];
      else
	err_flag = true;
      break;
    case WEIBULL:
      if (ranVarTypesU[i] == STD_NORMAL)
	// Phi(z) = F(x) = 1 - e^(-(x/beta)^alpha)
	z_vars[i] = Phi_inverse(weibull_cdf(x_vars[i], ranVarAddtlParamsX[i][0],
					    ranVarAddtlParamsX[i][1]));
      else if (ranVarTypesU[i] == WEIBULL) // Golub-Welsch: no transform
	z_vars[i] = x_vars[i];
      else
	err_flag = true;
      break;
    case HISTOGRAM_BIN:
      if (ranVarTypesU[i] == STD_NORMAL)
	z_vars[i] = Phi_inverse(histogram_bin_cdf(x_vars[i],
						  ranVarAddtlParamsX[i]));
      else if (ranVarTypesU[i] == STD_UNIFORM)
	z_vars[i] = std_uniform_cdf_inverse(histogram_bin_cdf(x_vars[i],
					    ranVarAddtlParamsX[i]));
      else if (ranVarTypesU[i] == HISTOGRAM_BIN)
	z_vars[i] = x_vars[i];
      else
	err_flag = true;
      break;
    }
    if (err_flag) {
      PCerr << "Error: unsupported variable mapping for variable " << i
	    << " in NatafTransformation::trans_X_to_Z()" << std::endl;
      abort_handler(-1);
    }
#ifdef DEBUG
    PCout << "X_to_Z: x[" << i << "] = " << x_vars[i]
	  <<    " --> z[" << i << "] = " << z_vars[i] << std::endl;
#endif // DEBUG
  }
}


/** This procedure computes the transformation from z to u space.
    u_vars is the vector of random variables in uncorrelated standard
    normal space (u-space).  z_vars is the vector of random variables
    in normal space with proper correlations (z-space). */
void NatafTransformation::trans_Z_to_U(RealVector& z_vars, RealVector& u_vars)
{
  // corrCholeskyFactorZ: the Cholesky factor of the modified correlation matrix

  int z_len = z_vars.length();
  if (u_vars.length() != z_len)
    u_vars.size(z_len);

  RealSolver corr_solver;
  corr_solver.setMatrix(  Teuchos::rcp(&corrCholeskyFactorZ, false) );
  corr_solver.setVectors( Teuchos::rcp(&u_vars, false),
			  Teuchos::rcp(&z_vars, false) );
  corr_solver.solveToRefinedSolution(true);
  corr_solver.solve();
}


/** This procedure modifies the correlation matrix input by the user for use in
    the Nataf distribution model (Der Kiureghian and Liu, ASCE JEM 112:1, 1986).
    It uses empirical expressionss derived from least-squares polynomial fits
    to numerical integration data.

    \li corrMatrixX: the correlation coefficient matrix of the random variables
    provided by the user

    \li mod_corr_matrix: modified correlation matrix

    \li corrCholeskyFactorZ: Cholesky factor of the modified correlation matrix
    for use in Z_to_U and U_to_Z transformations.

    Note: The modification is exact for normal-normal, lognormal-lognormal, and
    normal-lognormal tranformations.  All other cases are approximations with
    some error as noted below. */
void NatafTransformation::transform_correlations()
{
  // ranVarTypesX:   the type of random variable.  Supported correlations are
  //                 NORMAL/LOGNORMAL/UNIFORM/EXPONENTIAL/GAMMA/GUMBEL/FRECHET/
  //                 WEIBULL.
  // ranVarMeansX:   means of the random variables in x-space
  // ranVarStdDevsX: standard deviations of the random variables in x-space
  // cf_var_i:       coefficient of variation of Xi
  // cf_var_j:       coefficient of variation of Xj

  if (!correlationFlagX)
    return;

  RealSymMatrix mod_corr_matrix(corrMatrixX); // copy

  size_t i, j,
    num_cdv = std::count(ranVarTypesX.begin(), ranVarTypesX.end(),
			 (short)CONTINUOUS_DESIGN),
    num_cdv_uv = ranVarTypesX.size() - std::count(ranVarTypesX.begin(),
      ranVarTypesX.end(), (short)CONTINUOUS_STATE);

  for (i=num_cdv; i<num_cdv_uv; i++) {
    for (j=num_cdv; j<i; j++) {

      if (ranVarTypesU[i] == STD_NORMAL && ranVarTypesU[j] == STD_NORMAL) {

	// -----------------------------
	// Der Kiureghian & Liu: Table 2
	// -----------------------------
	//
	//        Normal <--> Uniform        (Max Error = 0.0%)
	//
	if ( (ranVarTypesX[i] == UNIFORM && ranVarTypesX[j] == NORMAL) ||
	     (ranVarTypesX[i] == NORMAL  && ranVarTypesX[j] == UNIFORM) )
	  mod_corr_matrix(i,j) *= 1.023326707946488488;
	//
	//        Normal <--> Exponential    (Max Error = 0.0%)
	//
	else if ( 
	  (ranVarTypesX[i] == EXPONENTIAL && ranVarTypesX[j] == NORMAL) ||
	  (ranVarTypesX[i] == NORMAL      && ranVarTypesX[j] == EXPONENTIAL) )
	  mod_corr_matrix(i,j) *= 1.107;
	//
	//        Normal <--> Gumbel         (Max Error = 0.0%)
	//
	else if ( (ranVarTypesX[i] == GUMBEL && ranVarTypesX[j] == NORMAL) ||
		  (ranVarTypesX[i] == NORMAL && ranVarTypesX[j] == GUMBEL) )
	  mod_corr_matrix(i,j) *= 1.031;

	// -----------------------------
	// Der Kiureghian & Liu: Table 3
	// -----------------------------
	//
	//        Normal <--> Lognormal      (Exact)
	//
	else if ( (ranVarTypesX[i] == LOGNORMAL && ranVarTypesX[j] == NORMAL) ||
		  (ranVarTypesX[i] == NORMAL && ranVarTypesX[j] == LOGNORMAL)) {
	  Real cf_var = (ranVarTypesX[i] == LOGNORMAL) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  mod_corr_matrix(i,j) *= cf_var/std::sqrt(log(1. + cf_var*cf_var));
	}
	//        Normal <--> Gamma          (Max Error = 0.0%)
	//
	else if ( (ranVarTypesX[i] == GAMMA  && ranVarTypesX[j] == NORMAL) ||
		  (ranVarTypesX[i] == NORMAL && ranVarTypesX[j] == GAMMA) ) {
	  Real cf_var = (ranVarTypesX[i] == GAMMA) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  mod_corr_matrix(i,j) *= 1.001 - 0.007*cf_var + 0.118*cf_var*cf_var;
	}
	//        Normal <--> Frechet        (Max Error = 0.1%)
	//
	else if ( (ranVarTypesX[i] == FRECHET && ranVarTypesX[j] == NORMAL) ||
		  (ranVarTypesX[i] == NORMAL  && ranVarTypesX[j] == FRECHET) ) {
	  Real cf_var = (ranVarTypesX[i] == FRECHET) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  mod_corr_matrix(i,j) *= 1.03 + 0.238*cf_var + 0.364*cf_var*cf_var;
	}
	//        Normal <--> Weibull        (Max Error = 0.1%)
	//
	else if ( (ranVarTypesX[i] == WEIBULL && ranVarTypesX[j] == NORMAL) ||
		  (ranVarTypesX[i] == NORMAL  && ranVarTypesX[j] == WEIBULL) ) {
	  Real cf_var = (ranVarTypesX[i] == WEIBULL) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  mod_corr_matrix(i,j) *= 1.031 - 0.195*cf_var + 0.328*cf_var*cf_var;
	}

	// -----------------------------
	// Der Kiureghian & Liu: Table 4
	// -----------------------------
	//
	//        Uniform <--> Uniform       (Max Error = 0.0%)
	//
	else if (ranVarTypesX[i] == UNIFORM && ranVarTypesX[j] == UNIFORM) {
	  const Real& corr = corrMatrixX(i,j);
	  mod_corr_matrix(i,j) *= 1.047 - 0.047*corr*corr;
	}
	//        Uniform <--> Exponential   (Max Error = 0.0%)
	//
	else if ( 
	  (ranVarTypesX[i] == UNIFORM     && ranVarTypesX[j] == EXPONENTIAL) ||
	  (ranVarTypesX[i] == EXPONENTIAL && ranVarTypesX[j] == UNIFORM) ) {
	  const Real& corr = corrMatrixX(i,j);
	  mod_corr_matrix(i,j) *= 1.133 + 0.029*corr*corr;
	}
	//        Uniform <--> Gumbel        (Max Error = 0.0%)
	//
	else if ( (ranVarTypesX[i] == UNIFORM && ranVarTypesX[j] == GUMBEL) ||
		  (ranVarTypesX[i] == GUMBEL  && ranVarTypesX[j] == UNIFORM) ) {
	  const Real& corr = corrMatrixX(i,j);
	  mod_corr_matrix(i,j) *= 1.055 + 0.015*corr*corr;
	}
	//        Exponential <--> Exponential (Max Error = 1.5%)
	//
	else if (ranVarTypesX[i] == EXPONENTIAL &&
		 ranVarTypesX[j] == EXPONENTIAL) {
	  const Real& corr = corrMatrixX(i,j);
	  mod_corr_matrix(i,j) *= 1.229 - 0.367*corr + 0.153*corr*corr;
	}
	//        Exponential <--> Gumbel    (Max Error = 0.2%)
	//
	else if (
	  (ranVarTypesX[i] == EXPONENTIAL && ranVarTypesX[j] == GUMBEL) ||
	  (ranVarTypesX[i] == GUMBEL      && ranVarTypesX[j] == EXPONENTIAL) ) {
	  const Real& corr = corrMatrixX(i,j);
	  mod_corr_matrix(i,j) *= 1.142 - 0.154*corr + 0.031*corr*corr;
	}
	//        Gumbel <--> Gumbel         (Max Error = 0.0%)
	//
	else if (ranVarTypesX[i] == GUMBEL && ranVarTypesX[j] == GUMBEL) {
	  const Real& corr = corrMatrixX(i,j);
	  mod_corr_matrix(i,j) *= 1.064 - 0.069*corr + 0.005*corr*corr;
	}

	// -----------------------------
	// Der Kiureghian & Liu: Table 5
	// -----------------------------
	//
	//        Lognormal <--> Uniform     (Max Error = 0.7%)
	//
	else if (
	  (ranVarTypesX[i] == LOGNORMAL && ranVarTypesX[j] == UNIFORM) ||
	  (ranVarTypesX[i] == UNIFORM   && ranVarTypesX[j] == LOGNORMAL) ) {
	  const Real& corr = corrMatrixX(i,j);
	  Real cf_var = (ranVarTypesX[i] == LOGNORMAL) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  mod_corr_matrix(i,j) *= 1.019 + 0.014*cf_var + 0.249*cf_var*cf_var
	                       +  0.01*corr*corr;
	}
	//        Lognormal <--> Exponential (Max Error = 1.6%)
	//
	else if (
	  (ranVarTypesX[i] == LOGNORMAL   && ranVarTypesX[j] == EXPONENTIAL) ||
	  (ranVarTypesX[i] == EXPONENTIAL && ranVarTypesX[j] == LOGNORMAL) ) {
	  const Real& corr = corrMatrixX(i,j);
	  Real cf_var = (ranVarTypesX[i] == LOGNORMAL) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  mod_corr_matrix(i,j) *= 1.098 + 0.003*corr + 0.019*cf_var
	                       +  0.025*corr*corr + 0.303*cf_var*cf_var
	                       -  0.437*corr*cf_var;
	}
	//        Lognormal <--> Gumbel      (Max Error = 0.3%)
	//
	else if ( (ranVarTypesX[i] == GUMBEL && ranVarTypesX[j] == LOGNORMAL) ||
		  (ranVarTypesX[i] == LOGNORMAL && ranVarTypesX[j] == GUMBEL)) {
	  const Real& corr = corrMatrixX(i,j);
	  Real cf_var = (ranVarTypesX[i] == LOGNORMAL) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  mod_corr_matrix(i,j) *= 1.029 + 0.001*corr + 0.014*cf_var
                               +  0.004*corr*corr + 0.233*cf_var*cf_var
                               -  0.197*corr*cf_var;
	}
	//        Gamma <--> Uniform         (Max Error = 0.1%)
	//
	else if ( (ranVarTypesX[i] == GAMMA   && ranVarTypesX[j] == UNIFORM) ||
		  (ranVarTypesX[i] == UNIFORM && ranVarTypesX[j] == GAMMA) ) {
	  const Real& corr = corrMatrixX(i,j);
	  Real cf_var = (ranVarTypesX[i] == GAMMA) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  mod_corr_matrix(i,j) *= 1.023 - 0.007*cf_var + 0.002*corr*corr 
                               +  0.127*cf_var*cf_var;
	}
	//        Gamma <--> Exponential     (Max Error = 0.9%)
	//
	else if (
	  (ranVarTypesX[i] == GAMMA       && ranVarTypesX[j] == EXPONENTIAL) ||
	  (ranVarTypesX[i] == EXPONENTIAL && ranVarTypesX[j] == GAMMA) ) {
	  const Real& corr = corrMatrixX(i,j);
	  Real cf_var = (ranVarTypesX[i] == GAMMA) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  mod_corr_matrix(i,j) *= 1.104 + 0.003*corr - 0.008*cf_var
	                       +  0.014*corr*corr + 0.173*cf_var*cf_var
	                       -  0.296*corr*cf_var;
	}
	//        Gamma <--> Gumbel          (Max Error = 0.3%)
	//
	else if ( (ranVarTypesX[i] == GUMBEL && ranVarTypesX[j] == GAMMA) ||
		  (ranVarTypesX[i] == GAMMA  && ranVarTypesX[j] == GUMBEL) ) {
	  const Real& corr = corrMatrixX(i,j);
	  Real cf_var = (ranVarTypesX[i] == GAMMA) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  mod_corr_matrix(i,j) *= 1.031 + 0.001*corr - 0.007*cf_var
                               +  0.003*corr*corr + 0.131*cf_var*cf_var
                               -  0.132*corr*cf_var;
	}
	//        Frechet <--> Uniform       (Max Error = 2.1%)
	//
	else if ( (ranVarTypesX[i] == FRECHET && ranVarTypesX[j] == UNIFORM) ||
		  (ranVarTypesX[i] == UNIFORM && ranVarTypesX[j] == FRECHET) ) {
	  const Real& corr = corrMatrixX(i,j);
	  Real cf_var = (ranVarTypesX[i] == FRECHET) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  mod_corr_matrix(i,j) *= 1.033 + 0.305*cf_var + 0.074*corr*corr
                               +  0.405*cf_var*cf_var;
	}
	//        Frechet <--> Exponential   (Max Error = 4.5%)
	//
	else if (
	  (ranVarTypesX[i] == FRECHET     && ranVarTypesX[j] == EXPONENTIAL) ||
	  (ranVarTypesX[i] == EXPONENTIAL && ranVarTypesX[j] == FRECHET) ) {
	  const Real& corr = corrMatrixX(i,j);
	  Real cf_var = (ranVarTypesX[i] == FRECHET) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  mod_corr_matrix(i,j) *= 1.109 - 0.152*corr + 0.361*cf_var
	                       +  0.130*corr*corr + 0.455*cf_var*cf_var
	                       -  0.728*corr*cf_var;
	}
	//        Frechet <--> Gumbel        (Max Error = 1.0%)
	//
	else if ( (ranVarTypesX[i] == FRECHET && ranVarTypesX[j] == GUMBEL) ||
		  (ranVarTypesX[i] == GUMBEL  && ranVarTypesX[j] == FRECHET) ) {
	  const Real& corr = corrMatrixX(i,j);
	  Real cf_var = (ranVarTypesX[i] == FRECHET) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  mod_corr_matrix(i,j) *= 1.056 - 0.06*corr + 0.263*cf_var
                               +  0.02*corr*corr + 0.383*cf_var*cf_var
                               -  0.332*corr*cf_var;
	}
	//        Weibull <--> Uniform       (Max Error = 0.5%)
	//
	else if ( (ranVarTypesX[i] == WEIBULL && ranVarTypesX[j] == UNIFORM) ||
		  (ranVarTypesX[i] == UNIFORM && ranVarTypesX[j] == WEIBULL) ) {
	  const Real& corr = corrMatrixX(i,j);
	  Real cf_var = (ranVarTypesX[i] == WEIBULL) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  mod_corr_matrix(i,j) *= 1.061 - 0.237*cf_var - 0.005*corr*corr 
                               +  0.379*cf_var*cf_var;
	}
	//        Weibull <--> Exponential   (Max Error = 0.4%)
	//
	else if (
	  (ranVarTypesX[i] == WEIBULL     && ranVarTypesX[j] == EXPONENTIAL) ||
	  (ranVarTypesX[i] == EXPONENTIAL && ranVarTypesX[j] == WEIBULL)){
	  const Real& corr = corrMatrixX(i,j);
	  Real cf_var = (ranVarTypesX[i] == WEIBULL) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  mod_corr_matrix(i,j) *= 1.147 + 0.145*corr - 0.271*cf_var
	                       +  0.010*corr*corr + 0.459*cf_var*cf_var
	                       -  0.467*corr*cf_var;
	}
	//        Weibull <--> Gumbel        (Max Error = 0.2%)
	//
	else if ( (ranVarTypesX[i] == GUMBEL  && ranVarTypesX[j] == WEIBULL) ||
		  (ranVarTypesX[i] == WEIBULL && ranVarTypesX[j] == GUMBEL) ) {
	  const Real& corr = corrMatrixX(i,j);
	  Real cf_var = (ranVarTypesX[i] == WEIBULL) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  mod_corr_matrix(i,j) *= 1.064 + 0.065*corr - 0.21*cf_var
                               +  0.003*corr*corr + 0.356*cf_var*cf_var
                               -  0.211*corr*cf_var;
	}

	// -----------------------------
	// Der Kiureghian & Liu: Table 6
	// -----------------------------
	//
	//        Lognormal <--> Lognormal   (Exact)
	//
	else if (ranVarTypesX[i] == LOGNORMAL && ranVarTypesX[j] == LOGNORMAL) {
	  Real cf_var_i = ranVarStdDevsX[i]/ranVarMeansX[i],
	       cf_var_j = ranVarStdDevsX[j]/ranVarMeansX[j];
	  // Note: "=" instead of "*=" eliminates redundant rho multiply/divide 
	  mod_corr_matrix(i,j) = log(1.+cf_var_i*cf_var_j*corrMatrixX(i,j)) /
	    std::sqrt(log(1. + cf_var_i*cf_var_i)*log(1. + cf_var_j*cf_var_j));
	}
	//        Gamma <--> Lognormal       (Max Error = 4.0%)
	//
	else if ( (ranVarTypesX[i] == GAMMA && ranVarTypesX[j] == LOGNORMAL) ||
		  (ranVarTypesX[i] == LOGNORMAL && ranVarTypesX[j] == GAMMA) ) {
	  Real cf_var_g = (ranVarTypesX[i] == GAMMA) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  Real cf_var_l = (ranVarTypesX[i] == LOGNORMAL) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  const Real& corr = corrMatrixX(i,j);
	  mod_corr_matrix(i,j) *= 1.001 + 0.033*corr + 0.004*cf_var_l 
	                       -  0.016*cf_var_g + 0.002*corr*corr
	                       +  0.223*cf_var_l*cf_var_l
	                       +  0.13*cf_var_g*cf_var_g - 0.104*corr*cf_var_l
	                       +  0.029*cf_var_l*cf_var_g - 0.119*corr*cf_var_g;
	}
	//        Frechet <--> Lognormal     (Max Error = 4.3%)
	//
	else if (
	  (ranVarTypesX[i] == FRECHET   && ranVarTypesX[j] == LOGNORMAL) ||
	  (ranVarTypesX[i] == LOGNORMAL && ranVarTypesX[j] == FRECHET) ) {
	  const Real& corr = corrMatrixX(i,j);
	  Real cf_var_l = (ranVarTypesX[i] == LOGNORMAL) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  Real cf_var_f = (ranVarTypesX[i] == FRECHET) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  mod_corr_matrix(i,j) *= 1.026 + 0.082*corr - 0.019*cf_var_l
                               +  0.222*cf_var_f + 0.018*corr*corr
                               +  0.288*cf_var_l*cf_var_l
	                       +  0.379*cf_var_f*cf_var_f - 0.441*corr*cf_var_l
	                       +  0.126*cf_var_l*cf_var_f - 0.277*corr*cf_var_f;
	}
	//        Weibull <--> Lognormal     (Max Error = 2.4%)
	//
	else if (
	  (ranVarTypesX[i] == WEIBULL   && ranVarTypesX[j] == LOGNORMAL) ||
	  (ranVarTypesX[i] == LOGNORMAL && ranVarTypesX[j] == WEIBULL) ) {
	  Real cf_var_w = (ranVarTypesX[i] == WEIBULL) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  Real cf_var_l = (ranVarTypesX[i] == LOGNORMAL) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  const Real& corr = corrMatrixX(i,j);
	  mod_corr_matrix(i,j) *= 1.031 + 0.052*corr + 0.011*cf_var_l 
                               -  0.21*cf_var_w + 0.002*corr*corr
	                       +  0.22*cf_var_l*cf_var_l
	                       +  0.35*cf_var_w*cf_var_w + 0.005*corr*cf_var_l
	                       +  0.009*cf_var_l*cf_var_w - 0.174*corr*cf_var_w;
	}
	//        Gamma <--> Gamma           (Max Error = 4.0%)
	//
	else if (ranVarTypesX[i] == GAMMA && ranVarTypesX[j] == GAMMA) {
	  Real cf_var_i = ranVarStdDevsX[i]/ranVarMeansX[i],
	       cf_var_j = ranVarStdDevsX[j]/ranVarMeansX[j];
	  const Real& corr = corrMatrixX(i,j);
	  mod_corr_matrix(i,j) *= 1.002 + 0.022*corr
	                       -  0.012*(cf_var_i + cf_var_j) + 0.001*corr*corr
	                       +  0.125*(cf_var_i*cf_var_i + cf_var_j*cf_var_j)
	                       -  0.077*corr*(cf_var_i + cf_var_j)
	                       +  0.014*cf_var_i*cf_var_j;
	}
	//        Frechet <--> Gamma         (Max Error = 4.2%)
	//
	else if ( (ranVarTypesX[i] == FRECHET && ranVarTypesX[j] == GAMMA) ||
		  (ranVarTypesX[i] == GAMMA   && ranVarTypesX[j] == FRECHET) ) {
	  const Real& corr = corrMatrixX(i,j);
	  Real cf_var_g = (ranVarTypesX[i] == GAMMA) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  Real cf_var_f = (ranVarTypesX[i] == FRECHET) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  mod_corr_matrix(i,j) *= 1.029 + 0.056*corr - 0.03*cf_var_g
                               +  0.225*cf_var_f + 0.012*corr*corr
	                       +  0.174*cf_var_g*cf_var_g
	                       +  0.379*cf_var_f*cf_var_f - 0.313*corr*cf_var_g
	                       +  0.075*cf_var_g*cf_var_f - 0.182*corr*cf_var_f;
	}
	//        Weibull <--> Gamma         (Max Error = 4.0%)
	//
	else if ( (ranVarTypesX[i] == GAMMA   && ranVarTypesX[j] == WEIBULL) ||
		  (ranVarTypesX[i] == WEIBULL && ranVarTypesX[j] == GAMMA) ) {
	  Real cf_var_g = (ranVarTypesX[i] == GAMMA) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  Real cf_var_w = (ranVarTypesX[i] == WEIBULL) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  const Real& corr = corrMatrixX(i,j);
	  mod_corr_matrix(i,j) *= 1.032 + 0.034*corr - 0.007*cf_var_g 
                               -  0.202*cf_var_w + 0.121*cf_var_g*cf_var_g
	                       +  0.339*cf_var_w*cf_var_w - 0.006*corr*cf_var_g
                               +  0.003*cf_var_w*cf_var_g - 0.111*corr*cf_var_w;
	}
	//        Frechet <--> Frechet       (Max Error = 4.3%)
	//
	else if (ranVarTypesX[i] == FRECHET && ranVarTypesX[j] == FRECHET) {
	  const Real& corr = corrMatrixX(i,j);
	  Real corr2 = corr*corr, cf_var_i  = ranVarStdDevsX[i]/ranVarMeansX[i],
	       cf_var_i2 = cf_var_i*cf_var_i,
	       cf_var_j  = ranVarStdDevsX[j]/ranVarMeansX[j],
	       cf_var_j2 = cf_var_j*cf_var_j;
	  mod_corr_matrix(i,j) *= 1.086 + 0.054*corr
	                       +  0.104*(cf_var_i + cf_var_j) - 0.055*corr2
	                       +  0.662*(cf_var_i2 + cf_var_j2)
	                       -  0.57*corr*(cf_var_i + cf_var_j)
	                       +  0.203*cf_var_i*cf_var_j - 0.02*corr2*corr
	                       -  0.218*(cf_var_i2*cf_var_i+cf_var_j2*cf_var_j)
	                       -  0.371*corr*(cf_var_i2 + cf_var_j2)
	                       +  0.257*corr2*(cf_var_i + cf_var_j)
	                       +  0.141*cf_var_i*cf_var_j*(cf_var_i + cf_var_j);
	}
	//        Weibull <--> Frechet       (Max Error = 3.8%)
	//
	else if ( (ranVarTypesX[i] == FRECHET && ranVarTypesX[j] == WEIBULL) ||
		  (ranVarTypesX[i] == WEIBULL && ranVarTypesX[j] == FRECHET) ) {
	  const Real& corr = corrMatrixX(i,j);
	  Real cf_var_w = (ranVarTypesX[i] == WEIBULL) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  Real cf_var_f = (ranVarTypesX[i] == FRECHET) ?
	    ranVarStdDevsX[i]/ranVarMeansX[i] :
	    ranVarStdDevsX[j]/ranVarMeansX[j];
	  mod_corr_matrix(i,j) *= 1.065 + 0.146*corr + 0.241*cf_var_f
                               -  0.259*cf_var_w + 0.013*corr*corr
	                       +  0.372*cf_var_f*cf_var_f
	                       +  0.435*cf_var_w*cf_var_w + 0.005*corr*cf_var_f
	                       +  0.034*cf_var_f*cf_var_w - 0.481*corr*cf_var_w;
	}
	//        Weibull <--> Weibull       (Max Error = 2.6%)
	//
	else if (ranVarTypesX[i] == WEIBULL && ranVarTypesX[j] == WEIBULL) {
	  Real cf_var_i = ranVarStdDevsX[i]/ranVarMeansX[i],
	       cf_var_j = ranVarStdDevsX[j]/ranVarMeansX[j];
	  const Real& corr = corrMatrixX(i,j);
	  mod_corr_matrix(i,j) *= 1.063 - 0.004*corr - 0.2*(cf_var_i + cf_var_j)
	                       -  0.001*corr*corr
                               +  0.337*(cf_var_i*cf_var_i + cf_var_j*cf_var_j)
	                       +  0.007*corr*(cf_var_i + cf_var_j)
	                       -  0.007*cf_var_i*cf_var_j;
	}
      }

      //else mod_corr_matrix(i,j) is unmodified

      // fill in upper half of symmetric matrix
      mod_corr_matrix(j,i) = mod_corr_matrix(i,j);
    }
  }

  // Cholesky decomposition for modified correlation matrix
  RealSpdSolver corr_solver;
  corr_solver.setMatrix( Teuchos::rcp(&mod_corr_matrix, false) );
  corr_solver.factor(); // Cholesky factorization (LL^T) in place
  // Define corrCholeskyFactorZ to be L by assigning the lower triangle.
  size_t num_active_vars = ranVarTypesX.size();
  if (corrCholeskyFactorZ.numRows() != num_active_vars ||
      corrCholeskyFactorZ.numCols() != num_active_vars)
    corrCholeskyFactorZ.shape(num_active_vars, num_active_vars);
  for (i=0; i<num_active_vars; i++)
    for (j=0; j<=i; j++)
      corrCholeskyFactorZ(i, j) = mod_corr_matrix(i, j);
#ifdef DEBUG
  PCout << "corrCholeskyFactorZ:\n" << corrCholeskyFactorZ;
#endif

  // could pre-compute L^-1 to avoid solving L u = z repeatedly for u
  //corrCholeskyFactorZInv.shape(num_active_vars, num_active_vars);
  //corrCholeskyFactorZInv = corrCholeskyFactorZ; // copy
  //RealSolver chol_solver;
  //chol_solver.setMatrix(corrCholeskyFactorZInv); 
  //chol_solver.invert();
  //PCout << "\ncorrCholeskyFactorZInv:" << corrCholeskyFactorZInv;
}


/** This procedure tranforms a gradient vector dg/dx from the original
    user-defined x-space (where evaluations are performed) to uncorrelated
    standard normal space (u-space) through application of the Jacobian dx/du.
    x_vars is the vector of random variables in x-space. */
void NatafTransformation::
trans_grad_X_to_U(const RealVector& fn_grad_x, RealVector& fn_grad_u,
		  const RealVector& x_vars, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  RealMatrix jacobian_xu;
  jacobian_dX_dU(x_vars, jacobian_xu);
  trans_grad_X_to_U(fn_grad_x, fn_grad_u, jacobian_xu, x_dvv, cv_ids);
}


/** This procedure tranforms a gradient vector dg/dx from the original
    user-defined x-space (where evaluations are performed) to uncorrelated
    standard normal space (u-space) through application of the Jacobian dx/du.
    This overloaded form allows for the separate calculation of jacobian_xu,
    as this matrix is independent of the response function index and can be
    pulled outside response function loops. */
void NatafTransformation::
trans_grad_X_to_U(const RealVector& fn_grad_x, RealVector& fn_grad_u,
		  const RealMatrix& jacobian_xu, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  // Jacobian dimensions = length of ranVarTypesX/U = model.cv()
  int x_len = jacobian_xu.numRows();
  if (x_dvv == cv_ids) { // standard DVV
    if (fn_grad_x.length() != x_len) {
      PCerr << "Error: bad fn_grad_x dimension in NatafTransformation::"
	    << "trans_grad_X_to_U()." << std::endl;
      abort_handler(-1);
    }
    if (fn_grad_u.length() != x_len)
      fn_grad_u.size(x_len);
    fn_grad_u.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1., jacobian_xu,
		       fn_grad_x, 0.);
  }
  else { // non-standard DVV
    RealVector fn_grad_x_trans(x_len), fn_grad_u_trans(x_len, false);
    size_t i, dvv_index, num_deriv_vars = x_dvv.size();
    SizetArray dvv_index_array(x_len);
    // extract relevant DVV components from fn_grad_x
    for (i=0; i<x_len; i++) {
      dvv_index_array[i] = dvv_index = find_index(x_dvv, cv_ids[i]);
      if (dvv_index != _NPOS)
	fn_grad_x_trans[i] = fn_grad_x(dvv_index);
    }
    // perform transformation using full Jacobian
    fn_grad_u_trans.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1., jacobian_xu,
			     fn_grad_x_trans, 0.);
    // copy relevant DVV components into fn_grad_u
    if (fn_grad_u.length() != num_deriv_vars)
      fn_grad_u.size(num_deriv_vars);
    for (i=0; i<x_len; i++) {
      dvv_index = dvv_index_array[i];
      if (dvv_index != _NPOS)
	fn_grad_u(dvv_index) = fn_grad_u_trans[i];
    }
  }

#ifdef DEBUG
  PCout << "Transformed fn_grad_x:\n" << fn_grad_x
        << "to fn_grad_u:\n" << fn_grad_u;
#endif
}


/** This procedure tranforms a gradient vector dg/du from uncorrelated standard
    space (u-space) to the original user-defined x-space through application of
    the Jacobian du/dx.  x_vars is the vector of random variables in x-space. */
void NatafTransformation::
trans_grad_U_to_X(const RealVector& fn_grad_u, RealVector& fn_grad_x,
		  const RealVector& x_vars, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  RealMatrix jacobian_ux;
  jacobian_dU_dX(x_vars, jacobian_ux);
  trans_grad_U_to_X(fn_grad_u, fn_grad_x, jacobian_ux, x_dvv, cv_ids);
}


/** This procedure tranforms a gradient vector dg/du from uncorrelated standard
    space (u-space) to the original user-defined x-space through application of
    the Jacobian du/dx.  This overloaded form allows for the separate
    calculation of jacobian_ux, as this matrix is independent of the response
    function index and can be pulled outside response function loops. */
void NatafTransformation::
trans_grad_U_to_X(const RealVector& fn_grad_u, RealVector& fn_grad_x,
		  const RealMatrix& jacobian_ux, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  // Jacobian dimensions = length of ranVarTypesX/U = model.cv()
  int u_len = jacobian_ux.numRows();
  if (x_dvv == cv_ids) { // standard DVV
    if (fn_grad_u.length() != u_len) {
      PCerr << "Error: bad fn_grad_u dimension in NatafTransformation::"
	    << "trans_grad_U_to_X()." << std::endl;
      abort_handler(-1);
    }
    if (fn_grad_x.length() != u_len)
      fn_grad_x.size(u_len);
    fn_grad_x.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1., jacobian_ux,
		       fn_grad_u, 0.);
  }
  else { // non-standard DVV
    RealVector fn_grad_u_trans(u_len), fn_grad_x_trans(u_len, false);
    size_t dvv_index, num_deriv_vars = x_dvv.size();
    SizetArray dvv_index_array(u_len);
    // extract relevant DVV components from fn_grad_u
    for (int i=0; i<u_len; i++) {
      dvv_index_array[i] = dvv_index = find_index(x_dvv, cv_ids[i]);
      if (dvv_index != _NPOS)
	fn_grad_u_trans[i] = fn_grad_u(dvv_index);
    }
    // perform transformation using full Jacobian
    fn_grad_x_trans.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1., jacobian_ux,
			     fn_grad_u_trans, 0.);
    // copy relevant DVV components into fn_grad_x
    if (fn_grad_x.length() != num_deriv_vars)
      fn_grad_x.size(num_deriv_vars);
    for (int i=0; i<u_len; i++) {
      dvv_index = dvv_index_array[i];
      if (dvv_index != _NPOS)
	fn_grad_x(dvv_index) = fn_grad_x_trans[i];
    }
  }

#ifdef DEBUG
  PCout << "Transformed fn_grad_u:\n" << fn_grad_u
        << "to fn_grad_x:\n" << fn_grad_x;
#endif
}


/** This procedure multiplies a gradient vector dg/dx from the
    original user-defined x-space (where evaluations are performed)
    with the design Jacobian dx/ds of the transformation x = x(u,s) to
    form the design gradient dg/ds.  x_vars is the vector of random
    variables in x-space. */
void NatafTransformation::
trans_grad_X_to_S(const RealVector& fn_grad_x, RealVector& fn_grad_s,
		  const RealVector& x_vars, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids,
		  SizetMultiArrayConstView acv_ids,
		  const SizetArray& acv_map1_indices,
		  const ShortArray& acv_map2_targets)
{
  RealMatrix jacobian_xs;
  jacobian_dX_dS(x_vars, jacobian_xs, cv_ids, acv_ids,
                 acv_map1_indices, acv_map2_targets);
  trans_grad_X_to_S(fn_grad_x, fn_grad_s, jacobian_xs, x_dvv, cv_ids, acv_ids,
		    acv_map1_indices, acv_map2_targets);
}


/** This procedure multiplies a gradient vector dg/dx from the
    original user-defined x-space (where evaluations are performed)
    with the design Jacobian dx/ds of the transformation x = x(u,s) to
    form the design gradient dg/ds.  This overloaded form allows
    for the separate calculation of jacobian_xs, as this matrix is
    independent of the response function index and can be pulled
    outside response function loops. */
void NatafTransformation::
trans_grad_X_to_S(const RealVector& fn_grad_x, RealVector& fn_grad_s,
		  const RealMatrix& jacobian_xs, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids,
		  SizetMultiArrayConstView acv_ids,
		  const SizetArray& acv_map1_indices,
		  const ShortArray& acv_map2_targets)
{
  // Jacobian dim is x_len by s_len, where
  // > x_len = length of ranVarTypesX/U = inner model.cv()
  // > s_len = acv_map1_indices.size()  = outer model.cv()
  int x_len = jacobian_xs.numRows(), s_len = jacobian_xs.numCols();
  if (acv_map1_indices.empty() || acv_map2_targets.empty()) {
    PCerr << "Error: NatafTransformation::trans_grad_X_to_S() requires primary "
	  << "and secondary variable mappings to define S." << std::endl;
    abort_handler(-1);
  }

  bool std_dvv = (x_dvv == cv_ids);
  bool mixed_s = std::find(acv_map2_targets.begin(), acv_map2_targets.end(),
                 (short)NO_TARGET) != acv_map2_targets.end() ? true : false;
  RealVector fn_grad_x_std, fn_grad_s_std;

  // manage size of fn_grad_x input
  if (std_dvv) { // standard x-space DVV
    if (fn_grad_x.length() != x_len) {
      PCerr << "Error: bad fn_grad_x dimension in NatafTransformation::"
	    << "trans_grad_X_to_S()." << std::endl;
      abort_handler(-1);
    }
  }
  else { // non-standard x-space DVV
    fn_grad_x_std.size(x_len);
    // extract relevant DVV components from fn_grad_x
    size_t i, dvv_index;
    for (i=0; i<x_len; i++) {
      dvv_index = find_index(x_dvv, cv_ids[i]);
      if (dvv_index != _NPOS)
	fn_grad_x_std[i] = fn_grad_x(dvv_index);
    }
  }

  // manage size of fn_grad_s output
  if (mixed_s || !std_dvv)
    fn_grad_s_std.size(s_len);
  size_t i, final_s_len;
  if (std_dvv)
    final_s_len = s_len;
  else {
    final_s_len = 0;
    for (i=0; i<s_len; i++)
      if ( std::find(x_dvv.begin(), x_dvv.end(),
                     acv_ids[acv_map1_indices[i]]) != x_dvv.end() )
	final_s_len++;
  }
  if (fn_grad_s.length() != final_s_len)
    fn_grad_s.size(final_s_len);

  // perform transformation using full Jacobian
  const RealVector& fn_grad_x_trans
    = (std_dvv) ? fn_grad_x : fn_grad_x_std;
  RealVector& fn_grad_s_trans
    = (mixed_s || !std_dvv) ? fn_grad_s_std : fn_grad_s;
  fn_grad_s_trans.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1., jacobian_xs,
			   fn_grad_x_trans, 0.);

  // reassemble final fn_grad_s
  if (mixed_s || !std_dvv) {
    size_t cntr = 0;
    for (i=0; i<s_len; i++) {
      size_t acv_id = acv_ids[acv_map1_indices[i]],
	dvv_index = find_index(x_dvv, acv_id);
      if (dvv_index != _NPOS)
	fn_grad_s(cntr++) = (acv_map2_targets[i] == NO_TARGET) ?
	  fn_grad_x(dvv_index) : // no distribution parameter: if the missing
	  // fn_grad_s component is available in fn_grad_x, then use it; else it
	  // must be updated separately (as in NonDLocalReliability::dg_ds_eval)
	  fn_grad_s_trans[i];    // use the Jacobian-transformed component
      else if (std_dvv)
	fn_grad_s(cntr++) = 0.;
    }
  }

#ifdef DEBUG
  PCout << "\nfn_grad_x:\n"       << fn_grad_x
        << "\nfn_grad_x_trans:\n" << fn_grad_x_trans
        << "\njacobian_xs:\n"     << jacobian_xs
        << "\nfn_grad_s_trans:\n" << fn_grad_s_trans
        << "\nfn_grad_s:\n"       << fn_grad_s << std::endl;
#endif
}


/** This procedure tranforms a Hessian matrix from the original
    user-defined x-space (where evaluations are performed) to
    uncorrelated standard normal space (u-space).  x_vars is the
    vector of the random variables in x-space. */
void NatafTransformation::
trans_hess_X_to_U(const RealSymMatrix& fn_hess_x, RealSymMatrix& fn_hess_u,
		  const RealVector& x_vars, const RealVector& fn_grad_x,
		  const SizetArray& x_dvv, SizetMultiArrayConstView cv_ids)
{
  RealMatrix jacobian_xu;
  jacobian_dX_dU(x_vars, jacobian_xu);

  RealSymMatrixArray hessian_xu;
  bool nonlinear_vars_map = false;
  size_t i, x_len = x_vars.length();
  for (i=0; i<x_len; i++)
    if ( ( ( ranVarTypesX[i] == CONTINUOUS_DESIGN  ||
	     ranVarTypesX[i] == CONTINUOUS_STATE ||
	     ranVarTypesX[i] == UNIFORM ||
	     ranVarTypesX[i] == CONTINUOUS_INTERVAL ) &&
	   ranVarTypesU[i]   != STD_UNIFORM ) ||
	 ( ranVarTypesX[i]   == NORMAL && ranVarTypesU[i]  != STD_NORMAL ) ||
	 ( ranVarTypesX[i]   == EXPONENTIAL &&
	   ranVarTypesU[i]   != STD_EXPONENTIAL ) ||
	 ( ranVarTypesX[i]   == BETA   && ranVarTypesU[i]  != STD_BETA   ) ||
	 ( ranVarTypesX[i]   == GAMMA  && ranVarTypesU[i]  != STD_GAMMA  ) ||
	 ( ( ranVarTypesX[i] == BOUNDED_NORMAL    ||
	     ranVarTypesX[i] == LOGNORMAL         ||
	     ranVarTypesX[i] == BOUNDED_LOGNORMAL ||
	     ranVarTypesX[i] == LOGUNIFORM || ranVarTypesX[i] == TRIANGULAR ||
	     ranVarTypesX[i] == GUMBEL     || ranVarTypesX[i] == FRECHET    ||
	     ranVarTypesX[i] == WEIBULL    || ranVarTypesX[i] == HISTOGRAM_BIN )
	  && ranVarTypesX[i] != ranVarTypesU[i] ) )
      { nonlinear_vars_map = true; break; }

  if (nonlinear_vars_map) // nonlinear transformation has Hessian
    hessian_d2X_dU2(x_vars, hessian_xu);

  trans_hess_X_to_U(fn_hess_x, fn_hess_u, jacobian_xu, hessian_xu,
		    fn_grad_x, x_dvv, cv_ids);
}


/** This procedure tranforms a Hessian matrix from the original
    user-defined x-space (where evaluations are performed) to
    uncorrelated standard normal space (u-space).  This overloaded
    form allows for the separate calculation of jacobian_xu and
    hessian_xu, since these are independent of the response function
    index and can be pulled outside response function loops. */
void NatafTransformation::
trans_hess_X_to_U(const RealSymMatrix& fn_hess_x, RealSymMatrix& fn_hess_u,
		  const RealMatrix& jacobian_xu,
		  const RealSymMatrixArray& hessian_xu,
		  const RealVector& fn_grad_x, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  // Jacobian dimensions = length of ranVarTypesX/U = model.cv()
  int  x_len   = jacobian_xu.numRows();
  bool std_dvv = (x_dvv == cv_ids); // standard DVV
  bool nonlinear_vars_map = !hessian_xu.empty();

  RealSymMatrix fn_hess_x_std, fn_hess_u_std;
  RealVector fn_grad_x_std;
  SizetArray dvv_index_array;
  if (std_dvv) {
    if (fn_hess_x.numRows() != x_len) {
      PCerr << "Error: bad fn_hess_x dimension in NatafTransformation::"
	    << "trans_hess_X_to_U()." << std::endl;
      abort_handler(-1);
    }
    if ( nonlinear_vars_map &&
	 ( fn_grad_x.length() != x_len || hessian_xu.size() != x_len ) ) {
      PCerr << "Error: bad dimension in NatafTransformation::"
	    << "trans_hess_X_to_U()." << std::endl;
      abort_handler(-1);
    }
    if (fn_hess_u.numRows() != x_len)
      fn_hess_u.shape(x_len);
  }
  else { // extract relevant DVV components from fn_grad_x & fn_hess_x
    fn_hess_x_std.shape(x_len);
    fn_hess_u_std.shape(x_len);
    if (nonlinear_vars_map)
      fn_grad_x_std.size(x_len);
    size_t i, j, dvv_index_i, dvv_index_j, num_deriv_vars = x_dvv.size();
    dvv_index_array.resize(x_len);
    for (i=0; i<x_len; i++)
      dvv_index_array[i] = dvv_index_i = find_index(x_dvv, cv_ids[i]);
    if (fn_hess_u.numRows() != num_deriv_vars)
      fn_hess_u.shape(num_deriv_vars);
    // extract relevant DVV components from fn_hess_x
    for (i=0; i<x_len; i++) {
      dvv_index_i = dvv_index_array[i];
      if (dvv_index_i != _NPOS) {
	if (nonlinear_vars_map)
	  fn_grad_x_std[i] = fn_grad_x(dvv_index_i);
	for (j=0; j<x_len; j++) {
	  dvv_index_j = dvv_index_array[j];
	  if (dvv_index_j != _NPOS)
	    fn_hess_x_std(i, j) = fn_hess_x(dvv_index_i, dvv_index_j);
	}
      }
    }
  }
  const RealVector&    fn_grad_x_trans = (std_dvv) ? fn_grad_x : fn_grad_x_std;
  const RealSymMatrix& fn_hess_x_trans = (std_dvv) ? fn_hess_x : fn_hess_x_std;
  RealSymMatrix&       fn_hess_u_trans = (std_dvv) ? fn_hess_u : fn_hess_u_std;

  // transform hess_x -> hess_u
  // d^2G/dU^2 = dG/dX^T d^2X/dU^2 + dX/dU^T d^2G/dX^2 dX/dU
  // Note: G(u) may have curvature even if g(x) is linear due to first term.
  Teuchos::symMatTripleProduct(Teuchos::TRANS, 1., fn_hess_x_trans,
                               jacobian_xu, fn_hess_u_trans);
#ifdef DEBUG
  PCout << "\nfnHessU 1st term J^T H_x J:" << fn_hess_u_trans;
#endif

  if (nonlinear_vars_map) { // nonlinear transformation has Hessian
    for (int i=0; i<x_len; i++) {
      //hessian_xu[i].Scale(fn_grad_x[i]);
      //fn_hess_u += hessian_xu[i];
      const Real&          fn_grad_x_i  = fn_grad_x_trans[i];
      const RealSymMatrix& hessian_xu_i = hessian_xu[i];
      for (int j=0; j<x_len; j++)
	for (int k=0; k<=j; k++)
	  fn_hess_u_trans(j,k) += fn_grad_x_i * hessian_xu_i(j,k);
#ifdef DEBUG
      PCout << "\nhessian_xu[" << i << "]:\n" << hessian_xu_i
	    << "fn_hess_u_trans increment:\n" << fn_hess_u_trans;
#endif
    }
  }

  if (!std_dvv) { // copy relevant DVV components back into fn_hess_u
    size_t i, j, dvv_index_i, dvv_index_j;
    for (i=0; i<x_len; i++) {
      dvv_index_i = dvv_index_array[i];
      if (dvv_index_i != _NPOS) {
	for (j=0; j<x_len; j++) {
	  dvv_index_j = dvv_index_array[j];
	  if (dvv_index_j != _NPOS)
	    fn_hess_u(dvv_index_i, dvv_index_j) = fn_hess_u_trans(i, j);
	}
      }
    }
  }

#ifdef DEBUG
  PCout << "Transformed fn_hess_x:\n" << fn_hess_x
        << "to fn_hess_u:\n" << fn_hess_u;
#endif
}


/** This procedure computes the Jacobian of the transformation x(u).
    x_vars is the vector of random variables in the original
    user-defined x-space. */
void NatafTransformation::
jacobian_dX_dU(const RealVector& x_vars, RealMatrix& jacobian_xu)
{
  if (correlationFlagX) {
    // dX/dZ = diagonal
    RealMatrix jacobian_xz;
    jacobian_dX_dZ(x_vars, jacobian_xz);

    // dX/dU = dX/dZ dZ/dU = dX/dZ L = dense if variables are correlated
    int x_len = x_vars.length();
    if (jacobian_xu.numRows() != x_len || jacobian_xu.numCols() != x_len)
      jacobian_xu.shape(x_len, x_len);
    jacobian_xu.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1., jacobian_xz,
			 corrCholeskyFactorZ, 0.);
  }
  else // dX/dU = dX/dZ since dZ/dU = I
    jacobian_dX_dZ(x_vars, jacobian_xu);

#ifdef DEBUG
  PCout << "jacobian_dX_dU:\n" << jacobian_xu;
#endif
}


/** This procedure computes the Jacobian of the transformation x(z).
    x_vars is the vector of random variables in the original
    user-defined x-space. */
void NatafTransformation::
jacobian_dX_dZ(const RealVector& x_vars, RealMatrix& jacobian_xz)
{
  int x_len = x_vars.length();
  if (jacobian_xz.numRows() != x_len || jacobian_xz.numCols() != x_len)
    jacobian_xz.shape(x_len, x_len);

  // Rackwitz-Fiessler: Phi(z) = F(x)
  // d/dz -> phi(z) = f(x) dx/dz
  //         dx/dz = phi(z)/f(x)
  // dX/dZ is diagonal as defined by differentiation of trans_Z_to_X()

  RealVector z_vars;
  bool need_z = false;
  for (size_t i=0; i<x_len; i++)
    if (ranVarTypesU[i] == STD_NORMAL &&
	ranVarTypesX[i] != NORMAL && ranVarTypesX[i] != LOGNORMAL)
      { need_z = true; break; }
  if (need_z)
    trans_X_to_Z(x_vars, z_vars);

  for (int i=0; i<x_len; i++) {
    bool err_flag = false;
    switch (ranVarTypesX[i]) {
    case CONTINUOUS_DESIGN: case CONTINUOUS_STATE: case CONTINUOUS_INTERVAL:
      if (ranVarTypesU[i] == STD_UNIFORM)
	jacobian_xz(i, i) = (ranVarUpperBndsX[i] - ranVarLowerBndsX[i])/2.;
      else
	err_flag = true;
      break;
    case NORMAL: // unbounded normal: z = (x - mean)/stdev
      if (ranVarTypesU[i] == STD_NORMAL)
	jacobian_xz(i, i) = ranVarStdDevsX[i];
      else
	err_flag = true;
      break;
    case BOUNDED_NORMAL: // bounded normal
      if (ranVarTypesU[i] == STD_NORMAL)
	jacobian_xz(i, i) = phi(z_vars[i]) /
	  bounded_normal_pdf(x_vars[i], ranVarMeansX[i], ranVarStdDevsX[i],
			     ranVarLowerBndsX[i], ranVarUpperBndsX[i]);
      else if (ranVarTypesU[i] == STD_UNIFORM)
	jacobian_xz(i, i) = std_uniform_pdf() /
	  bounded_normal_pdf(x_vars[i], ranVarMeansX[i], ranVarStdDevsX[i],
			     ranVarLowerBndsX[i], ranVarUpperBndsX[i]);
      else if (ranVarTypesU[i] == BOUNDED_NORMAL) // Golub-Welsch: no transform
	jacobian_xz(i, i) = 1.;
      else
	err_flag = true;
      break;
    case LOGNORMAL: // unbounded lognormal:  z = (ln x - lambda)/zeta
      if (ranVarTypesU[i] == STD_NORMAL)
	jacobian_xz(i, i) = ranVarAddtlParamsX[i][1] * x_vars[i];
      else if (ranVarTypesU[i] == LOGNORMAL) // Golub-Welsch: no transform
	jacobian_xz(i, i) = 1.;
      else
	err_flag = true;
      break;
    case BOUNDED_LOGNORMAL: // bounded lognormal
      if (ranVarTypesU[i] == STD_NORMAL)
	jacobian_xz(i, i) = phi(z_vars[i]) /
	  bounded_lognormal_pdf(x_vars[i], ranVarMeansX[i], ranVarStdDevsX[i],
				ranVarLowerBndsX[i], ranVarUpperBndsX[i]);
      else if (ranVarTypesU[i] == STD_UNIFORM)
	jacobian_xz(i, i) = std_uniform_pdf() /
	  bounded_lognormal_pdf(x_vars[i], ranVarMeansX[i], ranVarStdDevsX[i],
				ranVarLowerBndsX[i], ranVarUpperBndsX[i]);
      else if (ranVarTypesU[i] == BOUNDED_LOGNORMAL) // Golub-Welsch: no x-form
	jacobian_xz(i, i) = 1.;
      else
	err_flag = true;
      break;
    case UNIFORM: {
      // F(x) = (x-L)/(U-L), f(x) = 1/(U-L)
      if (ranVarTypesU[i] == STD_UNIFORM) { // linear scaling
	Real scale = ranVarUpperBndsX[i] - ranVarLowerBndsX[i];
	jacobian_xz(i, i) = scale/2.;
      }
      else if (ranVarTypesU[i] == STD_NORMAL) // nonlinear transformation
	jacobian_xz(i, i) = phi(z_vars[i]) /
	  uniform_pdf(ranVarLowerBndsX[i], ranVarUpperBndsX[i]);
      else
	err_flag = true;
      break;
    }
    case LOGUNIFORM:
      if (ranVarTypesU[i] == LOGUNIFORM) // Golub-Welsch: no transform
	jacobian_xz(i, i) = 1.;
      else {
	// F(x) = (ln x - ln L)/(ln U - ln L), f(x) = 1/x/(ln U - ln L)
	if (ranVarTypesU[i] == STD_UNIFORM) {
	  Real log_range = log(ranVarUpperBndsX[i]) - log(ranVarLowerBndsX[i]);
	  jacobian_xz(i, i) = x_vars[i]*log_range/2.;
	}
	else if (ranVarTypesU[i] == STD_NORMAL)
	  jacobian_xz(i, i) = phi(z_vars[i]) /
	    loguniform_pdf(x_vars[i], ranVarLowerBndsX[i], ranVarUpperBndsX[i]);
	else
	  err_flag = true;
      }
      break;
    case TRIANGULAR:
      if (ranVarTypesU[i] == TRIANGULAR) // Golub-Welsch: no transform
	jacobian_xz(i, i) = 1.;
      else {
	Real pdf = triangular_pdf(x_vars[i], ranVarAddtlParamsX[i][0],
				  ranVarLowerBndsX[i], ranVarUpperBndsX[i]);
	if (ranVarTypesU[i] == STD_UNIFORM)
	  jacobian_xz(i, i) = std_uniform_pdf() / pdf;
	else if (ranVarTypesU[i] == STD_NORMAL)
	  jacobian_xz(i, i) = phi(z_vars[i])/pdf;
	else
	  err_flag = true;
      }
      break;
    case EXPONENTIAL: {
      const Real& beta = ranVarAddtlParamsX[i][0];
      if (ranVarTypesU[i] == STD_EXPONENTIAL) // linear scaling
	jacobian_xz(i, i) = beta;
      else if (ranVarTypesU[i] == STD_NORMAL) // nonlinear transformation
	jacobian_xz(i, i) = phi(z_vars[i]) / exponential_pdf(x_vars[i], beta);
      else
	err_flag = true;
      break;
    }
    case BETA: {
      const Real& lwr = ranVarLowerBndsX[i];
      const Real& upr = ranVarUpperBndsX[i];
      if (ranVarTypesU[i] == STD_BETA) // linear scaling
	jacobian_xz(i, i) = (upr - lwr)/2.;
      else { 
	Real pdf = beta_pdf(x_vars[i], ranVarAddtlParamsX[i][0],
			    ranVarAddtlParamsX[i][1], lwr, upr);
	if (ranVarTypesU[i] == STD_NORMAL) // nonlinear transformation
	  jacobian_xz(i, i) = phi(z_vars[i]) / pdf;
	else if (ranVarTypesU[i] == STD_UNIFORM) // nonlinear transformation
	  jacobian_xz(i, i) = std_uniform_pdf() / pdf;
	else
	  err_flag = true;
      }
      break;
    }
    case GAMMA:
      if (ranVarTypesU[i] == STD_GAMMA) // linear scaling
	jacobian_xz(i, i) = ranVarAddtlParamsX[i][1]; // beta
      else if (ranVarTypesU[i] == STD_NORMAL) // nonlinear transformation
	jacobian_xz(i, i) = phi(z_vars[i]) / gamma_pdf(x_vars[i],
	  ranVarAddtlParamsX[i][0], ranVarAddtlParamsX[i][1]); // alpha, beta
      else
	err_flag = true;
      break;
    case GUMBEL:
      if (ranVarTypesU[i] == STD_NORMAL)
	jacobian_xz(i, i) = phi(z_vars[i]) / gumbel_pdf(x_vars[i],
	  ranVarAddtlParamsX[i][0], ranVarAddtlParamsX[i][1]);
      else if (ranVarTypesU[i] == GUMBEL) // Golub-Welsch: no transform
	jacobian_xz(i, i) = 1.;
      else
	err_flag = true;
      break;
    case FRECHET:
      if (ranVarTypesU[i] == STD_NORMAL)
	jacobian_xz(i, i) = phi(z_vars[i]) / frechet_pdf(x_vars[i],
	  ranVarAddtlParamsX[i][0], ranVarAddtlParamsX[i][1]);
      else if (ranVarTypesU[i] == FRECHET) // Golub-Welsch: no transform
	jacobian_xz(i, i) = 1.;
      else
	err_flag = true;
      break;
    case WEIBULL:
      if (ranVarTypesU[i] == STD_NORMAL)
	jacobian_xz(i, i) = phi(z_vars[i]) / weibull_pdf(x_vars[i],
	  ranVarAddtlParamsX[i][0], ranVarAddtlParamsX[i][1]);
      else if (ranVarTypesU[i] == WEIBULL) // Golub-Welsch: no transform
	jacobian_xz(i, i) = 1.;
      else
	err_flag = true;
      break;
    case HISTOGRAM_BIN:
      if (ranVarTypesU[i] == STD_NORMAL)
	jacobian_xz(i, i) = phi(z_vars[i]) /
	  histogram_bin_pdf(x_vars[i], ranVarAddtlParamsX[i]);
      else if (ranVarTypesU[i] == STD_UNIFORM)
	jacobian_xz(i, i) = std_uniform_pdf() /
	  histogram_bin_pdf(x_vars[i], ranVarAddtlParamsX[i]);
      else if (ranVarTypesU[i] == HISTOGRAM_BIN)
	jacobian_xz(i, i) = 1.;
      else
	err_flag = true;
      break;
    }
    if (err_flag) {
      PCerr << "Error: unsupported variable mapping for variable " << i
	    << " in NatafTransformation::jacobian_dX_dZ()" << std::endl;
      abort_handler(-1);
    }
  }
}


/** This procedure computes the Jacobian of the transformation u(x).
    x_vars is the vector of random variables in the original
    user-defined x-space. */
void NatafTransformation::
jacobian_dU_dX(const RealVector& x_vars, RealMatrix& jacobian_ux)
{
  if (correlationFlagX) {
    // dZ/dX = diagonal
    RealMatrix jacobian_zx;
    jacobian_dZ_dX(x_vars, jacobian_zx);

    // dU/dX = dU/dZ dZ/dX = L^-1 dZ/dX = dense if variables are correlated
    // Solve as L dU/dX = dZ/dX
    RealSolver corr_solver;
    corr_solver.setMatrix( Teuchos::rcp(&corrCholeskyFactorZ, false) );
    int x_len = x_vars.length();
    if (jacobian_ux.numRows() != x_len || jacobian_ux.numCols() != x_len)
      jacobian_ux.shape(x_len, x_len);
    corr_solver.setVectors( Teuchos::rcp(&jacobian_ux, false),
                            Teuchos::rcp(&jacobian_zx, false) );
    corr_solver.solveToRefinedSolution(true);
    corr_solver.solve();
  }
  else // dU/dX = dZ/dX since dU/dZ = I
    jacobian_dZ_dX(x_vars, jacobian_ux);

#ifdef DEBUG
  PCout << "jacobian_dU_dX:\n" << jacobian_ux;
#endif
}


/** This procedure computes the Jacobian of the transformation z(x).
    x_vars is the vector of random variables in the original
    user-defined x-space. */
void NatafTransformation::
jacobian_dZ_dX(const RealVector& x_vars, RealMatrix& jacobian_zx) 
{
  int x_len = x_vars.length();
  if (jacobian_zx.numRows() != x_len || jacobian_zx.numCols() != x_len)
    jacobian_zx.shape(x_len, x_len);

  // Rackwitz-Fiessler: Phi(z) = F(x)
  // d/dx -> phi(z) dz/dx = f(x)
  //         dz/dx = f(x)/phi(z)
  // dZ/dX is diagonal as defined by differentiation of trans_X_to_Z()

  RealVector z_vars;
  bool need_z = false;
  for (size_t i=0; i<x_len; i++)
    if (ranVarTypesU[i] == STD_NORMAL &&
	ranVarTypesX[i] != NORMAL && ranVarTypesX[i] != LOGNORMAL)
      { need_z = true; break; }
  if (need_z)
    trans_X_to_Z(x_vars, z_vars);

  for (int i=0; i<x_len; i++) {
    bool err_flag = false;
    switch (ranVarTypesX[i]) {
    case CONTINUOUS_DESIGN: case CONTINUOUS_STATE: case CONTINUOUS_INTERVAL:
      if (ranVarTypesU[i] == STD_UNIFORM)
	jacobian_zx(i, i) = 2./(ranVarUpperBndsX[i] - ranVarLowerBndsX[i]);
      else
	err_flag = true;
      break;
    case NORMAL: // unbounded normal: z = (x - mean)/stdev
      if (ranVarTypesU[i] == STD_NORMAL)
	jacobian_zx(i, i) = 1./ranVarStdDevsX[i];
      else
	err_flag = true;
      break;
    case BOUNDED_NORMAL: // bounded normal
      if (ranVarTypesU[i] == STD_NORMAL)
	jacobian_zx(i, i) = bounded_normal_pdf(x_vars[i], ranVarMeansX[i],
	  ranVarStdDevsX[i], ranVarLowerBndsX[i], ranVarUpperBndsX[i]) /
	  phi(z_vars[i]);
      else if (ranVarTypesU[i] == STD_UNIFORM)
	jacobian_zx(i, i) = bounded_normal_pdf(x_vars[i], ranVarMeansX[i],
	  ranVarStdDevsX[i], ranVarLowerBndsX[i], ranVarUpperBndsX[i]) /
	  std_uniform_pdf();
      else if (ranVarTypesU[i] == BOUNDED_NORMAL) // Golub-Welsch: no transform
	jacobian_zx(i, i) = 1.;
      else
	err_flag = true;
      break;
    case LOGNORMAL: // unbounded lognormal: z = (ln x - lambda)/zeta
      if (ranVarTypesU[i] == STD_NORMAL)
	jacobian_zx(i, i) = 1./ranVarAddtlParamsX[i][1]/x_vars[i];
      else if (ranVarTypesU[i] == LOGNORMAL) // Golub-Welsch: no transform
	jacobian_zx(i, i) = 1.;
      else
	err_flag = true;
      break;
    case BOUNDED_LOGNORMAL: // bounded lognormal
      if (ranVarTypesU[i] == STD_NORMAL)
	jacobian_zx(i, i) = bounded_lognormal_pdf(x_vars[i], ranVarMeansX[i],
	  ranVarStdDevsX[i], ranVarLowerBndsX[i], ranVarUpperBndsX[i]) /
	  phi(z_vars[i]);
      else if (ranVarTypesU[i] == STD_UNIFORM)
	jacobian_zx(i, i) = bounded_lognormal_pdf(x_vars[i], ranVarMeansX[i],
	  ranVarStdDevsX[i], ranVarLowerBndsX[i], ranVarUpperBndsX[i]) /
	  std_uniform_pdf();
      else if (ranVarTypesU[i] == BOUNDED_LOGNORMAL) // Golub-Welsch: no x-form
	jacobian_zx(i, i) = 1.;
      else
	err_flag = true;
      break;
    case UNIFORM: {
      // F(x) = (x-L)/(U-L), f(x) = 1/(U-L)
      if (ranVarTypesU[i] == STD_UNIFORM) { // linear scaling
	Real scale = ranVarUpperBndsX[i] - ranVarLowerBndsX[i];
	jacobian_zx(i, i) = 2./scale;
      }
      else if (ranVarTypesU[i] == STD_NORMAL)
	jacobian_zx(i, i) = uniform_pdf(ranVarLowerBndsX[i],
					ranVarUpperBndsX[i]) / phi(z_vars[i]);
      else
	err_flag = true;
      break;
    }
    case LOGUNIFORM:
      if (ranVarTypesU[i] == LOGUNIFORM) // Golub-Welsch: no transform
	jacobian_zx(i, i) = 1.;
      else {
	// F(x) = (ln x - ln L)/(ln U - ln L), f(x) = 1/x/(ln U - ln L)
	if (ranVarTypesU[i] == STD_UNIFORM) {
	  Real log_range = log(ranVarUpperBndsX[i]) - log(ranVarLowerBndsX[i]);
	  jacobian_zx(i, i) = 2./x_vars[i]/log_range;
	}
	else if (ranVarTypesU[i] == STD_NORMAL)
	  jacobian_zx(i, i) = loguniform_pdf(x_vars[i], ranVarLowerBndsX[i],
	    ranVarUpperBndsX[i]) / phi(z_vars[i]);
	else
	  err_flag = true;
      }
      break;
    case TRIANGULAR:
      if (ranVarTypesU[i] == TRIANGULAR) // Golub-Welsch: no transform
	jacobian_zx(i, i) = 1.;
      else {
	Real pdf = triangular_pdf(x_vars[i], ranVarAddtlParamsX[i][0],
				  ranVarLowerBndsX[i], ranVarUpperBndsX[i]);
	if (ranVarTypesU[i] == STD_UNIFORM)
	  jacobian_zx(i, i) = pdf / std_uniform_pdf();
	else if (ranVarTypesU[i] == STD_NORMAL)
	  jacobian_zx(i, i) = pdf / phi(z_vars[i]);
	else
	  err_flag = true;
      }
      break;
    case EXPONENTIAL: {
      const Real& beta = ranVarAddtlParamsX[i][0];
      if (ranVarTypesU[i] == STD_EXPONENTIAL) // linear scaling
	jacobian_zx(i, i) = 1./beta;
      else if (ranVarTypesU[i] == STD_NORMAL)
	jacobian_zx(i, i) = exponential_pdf(x_vars[i], beta) / phi(z_vars[i]);
      else
	err_flag = true;
      break;
    }
    case BETA: {
      const Real& lwr = ranVarLowerBndsX[i];
      const Real& upr = ranVarUpperBndsX[i];
      if (ranVarTypesU[i] == STD_BETA) // linear scaling
	jacobian_zx(i, i) = 2./(upr - lwr);
      else {
	Real pdf = beta_pdf(x_vars[i], ranVarAddtlParamsX[i][0],
			    ranVarAddtlParamsX[i][1], lwr, upr);
	if (ranVarTypesU[i] == STD_NORMAL)
	  jacobian_zx(i, i) = pdf / phi(z_vars[i]);
	else if (ranVarTypesU[i] == STD_UNIFORM)
	  jacobian_zx(i, i) = pdf / std_uniform_pdf();
	else
	  err_flag = true;
      }
      break;
    }
    case GAMMA: {
      if (ranVarTypesU[i] == STD_GAMMA) // linear scaling
	jacobian_zx(i, i) = 1./ranVarAddtlParamsX[i][1];
      else if (ranVarTypesU[i] == STD_NORMAL)
	jacobian_zx(i, i) = gamma_pdf(x_vars[i], ranVarAddtlParamsX[i][0],
	  ranVarAddtlParamsX[i][1]) / phi(z_vars[i]);
      else
	err_flag = true;
      break;
    }
    case GUMBEL: {
      if (ranVarTypesU[i] == STD_NORMAL)
	jacobian_zx(i, i) = gumbel_pdf(x_vars[i], ranVarAddtlParamsX[i][0],
	  ranVarAddtlParamsX[i][1]) / phi(z_vars[i]);
      else if (ranVarTypesU[i] == GUMBEL) // Golub-Welsch: no transform
	jacobian_zx(i, i) = 1.;
      else
	err_flag = true;
      break;
    }
    case FRECHET: {
      if (ranVarTypesU[i] == STD_NORMAL)
	jacobian_zx(i, i) = frechet_pdf(x_vars[i], ranVarAddtlParamsX[i][0],
	  ranVarAddtlParamsX[i][1]) / phi(z_vars[i]);
      else if (ranVarTypesU[i] == FRECHET) // Golub-Welsch: no transform
	jacobian_zx(i, i) = 1.;
      else
	err_flag = true;
      break;
    }
    case WEIBULL: {
      if (ranVarTypesU[i] == STD_NORMAL)
	jacobian_zx(i, i) = weibull_pdf(x_vars[i], ranVarAddtlParamsX[i][0],
	  ranVarAddtlParamsX[i][1]) / phi(z_vars[i]);
      else if (ranVarTypesU[i] == WEIBULL) // Golub-Welsch: no transform
	jacobian_zx(i, i) = 1.;
      else
	err_flag = true;
      break;
    }
    case HISTOGRAM_BIN:
      if (ranVarTypesU[i] == STD_NORMAL)
	jacobian_zx(i, i) = histogram_bin_pdf(x_vars[i], ranVarAddtlParamsX[i])
	  / phi(z_vars[i]);
      else if (ranVarTypesU[i] == STD_UNIFORM)
	jacobian_zx(i, i) = histogram_bin_pdf(x_vars[i], ranVarAddtlParamsX[i])
	  / std_uniform_pdf();
      else if (ranVarTypesU[i] == HISTOGRAM_BIN)
	jacobian_zx(i, i) = 1.;
      else
	err_flag = true;
      break;
    }
    if (err_flag) {
      PCerr << "Error: unsupported variable mapping for variable " << i
	    << " in NatafTransformation::jacobian_dZ_dX()" << std::endl;
      abort_handler(-1);
    }
  }
}


/** This procedure computes the derivative of the original variables x
    with respect to the random variable distribution parameters s.
    This provides the design Jacobian of the transformation for use in
    computing statistical design sensitivities for OUU. */
void NatafTransformation::
jacobian_dX_dS(const RealVector& x_vars, RealMatrix& jacobian_xs,
	       SizetMultiArrayConstView cv_ids,
	       SizetMultiArrayConstView acv_ids,
	       const SizetArray& acv_map1_indices,
	       const ShortArray& acv_map2_targets)
{
  // Rectangular Jacobian = Gradient^T = num_X by num_S where num_S is the total
  // number of active continuous vars flowed down from a higher iteration level.
  // The number of distribution parameter insertions is <= num_S.
  size_t i, j, num_var_map_1c = acv_map1_indices.size();
  int x_len = x_vars.length();
  if (jacobian_xs.numRows() != x_len || jacobian_xs.numCols() != num_var_map_1c)
    jacobian_xs.shape(x_len, num_var_map_1c);

  // dX/dS is derived by differentiating trans_Z_to_X with respect to S.
  // For the uncorrelated case, u and z are constants.  For the correlated
  // case, u is a constant, but z(s) = L(s) u due to Nataf dependence on s
  // and dz/ds = dL/ds u.

  RealVector z_vars;
  trans_X_to_Z(x_vars, z_vars);

  bool need_xs = false, beta_gamma_map = false;
  // For distributions without simple closed-form CDFs (beta, gamma), dx/ds is
  // computed numerically.  If uncorrelated, then this is only needed if the
  // beta/gamma distribution parameters are design variables.  If correlated,
  // then the beta/gamma distribution parameters do not have to be design
  // variables (dx/ds for beta/gamma x will include a dz/ds contribution).
  if (correlationFlagX) {
    size_t num_cdv = std::count(ranVarTypesX.begin(), ranVarTypesX.end(),
				(short)CONTINUOUS_DESIGN),
      num_cdv_uv = ranVarTypesX.size() - std::count(ranVarTypesX.begin(),
	ranVarTypesX.end(), (short)CONTINUOUS_STATE);
    for (i=num_cdv; i<num_cdv_uv; i++) {
      if ( (ranVarTypesX[i] == BETA || ranVarTypesX[i] == GAMMA) &&
	    ranVarTypesX[i] != ranVarTypesU[i] ) {
	beta_gamma_map = true;
	for (j=num_cdv; j<num_cdv_uv; j++)
	  if (i != j && fabs(corrMatrixX(i, j)) > SMALL_NUMBER)
	    { need_xs = true; break; }
      }
    }
  }
  if ( (std::find(acv_map2_targets.begin(), acv_map2_targets.end(),
                  (short)BE_ALPHA) != acv_map2_targets.end())  ||
       (std::find(acv_map2_targets.begin(), acv_map2_targets.end(),
                  (short)BE_BETA) != acv_map2_targets.end())   ||
       (std::find(acv_map2_targets.begin(), acv_map2_targets.end(),
                  (short)GA_ALPHA) != acv_map2_targets.end()) ||
       ( beta_gamma_map && ( (std::find(acv_map2_targets.begin(),
                              acv_map2_targets.end(), (short)BE_LWR_BND)
			      != acv_map2_targets.end()) ||
			     (std::find(acv_map2_targets.begin(),
                              acv_map2_targets.end(), (short)BE_UPR_BND)
			      != acv_map2_targets.end()) ||
			     (std::find(acv_map2_targets.begin(),
                              acv_map2_targets.end(), (short)GA_BETA)
			      != acv_map2_targets.end()) ) ) )
    need_xs = true;
  // the entire numerical jacobian is computed, even though only the
  // beta/gamma rows are needed
  RealMatrix num_dx_ds, num_dz_ds;
  if (need_xs || correlationFlagX)
    numerical_design_jacobian(x_vars, need_xs, num_dx_ds,
			      correlationFlagX, num_dz_ds, cv_ids, acv_ids,
			      acv_map1_indices, acv_map2_targets);
  if (need_xs) {
    for (j=0; j<x_len; j++) {            // loop over X
      switch (ranVarTypesX[j]) {
      case BETA: case GAMMA:
	for (i=0; i<num_var_map_1c; i++) // loop over S
	  jacobian_xs(j, i) = num_dx_ds(j, i);
	break;
      }
    }
  }

  Real dbl_inf = std::numeric_limits<Real>::infinity();
  for (i=0; i<num_var_map_1c; i++) { // loop over S
    size_t cv_index = find_index(cv_ids, acv_ids[acv_map1_indices[i]]);
    // If x_dvv were passed, it would be possible to distinguish different
    // fn_grad_x components, allowing passthrough for computing fn_grad_s for
    // augmented design variables.  For now, this has to be handled spearately
    // in NonDLocalReliability::dg_ds_eval() and NonD::trans_grad_X_to_S().
    //if (cv_index == _NPOS) // augmented variable: define identity mapping
    //  jacobian_xs(dvv_index, i) = 1.;
    //else {
    if (cv_index != _NPOS) {
      short target2 = acv_map2_targets[i];
      for (j=0; j<x_len; j++) {      // loop over X
	// Jacobian row    = Z value = j
	// Jacobian column = S value = i
	switch (ranVarTypesX[j]) {
	case CONTINUOUS_DESIGN: { // x = L + (z + 1)*(U - L)/2
	  if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	    switch (target2) {
	    case CDV_LWR_BND: // Deriv of CDV w.r.t. its Lower Bound
	      jacobian_xs(j, i) = (1. - z_vars[j])/2.; break;
	    case CDV_UPR_BND: // Deriv of CDV w.r.t. its Upper Bound
	      jacobian_xs(j, i) = (z_vars[j] + 1.)/2.; break;
	    case NO_TARGET:   // can occur for all_variables Jacobians
	      jacobian_xs(j, i) = 0.;                  break;
	    default:
	      PCerr << "Error: secondary mapping failure for CONTINUOUS_DESIGN "
		    << "in NatafTransformation::jacobian_dX_dS()." << std::endl;
	      abort_handler(-1);
	      break;
	    }
	  }
	  // Deriv of CDV w.r.t. other distribution params should always be 0
	  //if (correlationFlagX)
	  //  jacobian_xs(j, i) += 
	  //    (ranVarUpperBndsX[j]-ranVarLowerBndsX[j])/2.*num_dz_ds(j, i);
	  break;
	}
	case NORMAL: { // x = z stdev + mean
	  if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	    switch (target2) {
	    case N_MEAN: // Deriv of Normal w.r.t. its Mean
	      jacobian_xs(j, i) = 1.;
	      break;
	    case N_STD_DEV: // Deriv of Normal w.r.t. its Std Deviation
	      jacobian_xs(j, i) = z_vars[j];
	      break;
	    //case N_LWR_BND: case N_UPR_BND: not supported
	    case NO_TARGET: default:
	      PCerr << "Error: secondary mapping failure for NORMAL in "
		    << "NatafTransformation::jacobian_dX_dS()." << std::endl;
	      abort_handler(-1);
	      break;
	    }
	  }
	  // Deriv of Normal w.r.t. any distribution parameter
	  if (correlationFlagX)
	    jacobian_xs(j, i) += ranVarStdDevsX[j]*num_dz_ds(j, i);
	  break;
	}
	case BOUNDED_NORMAL: { // bounded normal
	  const Real& mean  = ranVarMeansX[j];
	  const Real& stdev = ranVarStdDevsX[j];
	  const Real& lwr   = ranVarLowerBndsX[j];
	  const Real& upr   = ranVarUpperBndsX[j];
	  const Real& x     = x_vars[j]; const Real& z = z_vars[j];
	  Real lms = (lwr > -dbl_inf) ? (lwr-mean)/stdev : -dbl_inf;
	  Real ums = (upr <  dbl_inf) ? (upr-mean)/stdev :  dbl_inf;
	  Real xms = (x-mean)/stdev, phi_x = phi(xms);
	  if (ranVarTypesU[j] == STD_NORMAL) {
	    if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	      Real phi_lms = (lwr > -dbl_inf) ? phi(lms) : 0.;
	      Real phi_ums = (upr <  dbl_inf) ? phi(ums) : 0.;
	      Real normcdf_comp = (z > 0.) ? Phi(-z) : 1. - Phi(z);
	      switch (target2) {
	      case N_MEAN: // Deriv of Bounded Normal w.r.t. its Mean
		jacobian_xs(j, i)
		  = 1. - (normcdf_comp*phi_lms + Phi(z)*phi_ums)/phi_x;
		break;
	      case N_STD_DEV: // Deriv of Bounded Normal w.r.t. its Std Dev
		jacobian_xs(j, i)
		  = xms - (normcdf_comp*phi_lms*lms + Phi(z)*phi_ums*ums)/phi_x;
		break;
	      case N_LWR_BND: // Deriv of Bounded Normal w.r.t. its lower bound
		jacobian_xs(j, i) = phi_lms/phi_x*normcdf_comp;
		break;
	      case N_UPR_BND: // Deriv of Bounded Normal w.r.t. its upper bound
		jacobian_xs(j, i) = phi_ums/phi_x*Phi(z);
		break;
	      case NO_TARGET: default:
		PCerr << "Error: secondary mapping failure for BOUNDED_NORMAL "
		      << "in NatafTransformation::jacobian_dX_dS()."<<std::endl;
		abort_handler(-1);
		break;
	      }
	    }
	    // Deriv of Normal w.r.t. any distribution parameter
	    if (correlationFlagX) {
	      Real Phi_lms = (lwr > -dbl_inf) ? Phi(lms) : 0.;
	      Real Phi_ums = (upr <  dbl_inf) ? Phi(ums) : 1.;
	      jacobian_xs(j, i)
		+= stdev*phi(z)*(Phi_ums - Phi_lms)/phi_x*num_dz_ds(j, i);
	    }
	  }
	  else if (ranVarTypesU[j] == BOUNDED_NORMAL) { // *** TO DO ***
	    PCerr << "Error: mapping not yet supported in NatafTransformation::"
		  << "jacobian_dX_dS()" << std::endl;
	    abort_handler(-1);
	  }
	  break;
	}
	case LOGNORMAL: { // x = exp(lambda + z zeta)
	  const Real& x = x_vars[j];
	  if (ranVarTypesU[j] == STD_NORMAL) {
	    if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	      bool ln_err_fact = (ranVarAddtlParamsX[j].length() > 2);
	      switch (target2) {
	      case LN_MEAN: { // Deriv of Lognormal w.r.t. its Mean
		// x = exp(lambda + z zeta)
		const Real& mean = ranVarMeansX[j];
		if (ln_err_fact) // mean, error factor spec
		  jacobian_xs(j, i) = x/mean;
		else {           // mean, std deviation spec
		  const Real& stdev = ranVarStdDevsX[j];
		  const Real& zeta  = ranVarAddtlParamsX[j][1];
		  Real mean_sq = mean*mean, var = stdev*stdev;
		  jacobian_xs(j, i) = x*(zeta*mean_sq + 2.*zeta*var -
		    z_vars[j]*var)/mean/zeta/(mean_sq + var);
		}
		break;
	      }
	      case LN_STD_DEV: // Deriv of Lognormal w.r.t. its Std Deviation
		if (ln_err_fact) { // mean, error factor spec
		  PCerr << "Error: derivative with respect to LN_STD_DEV is "
			<< "unsupported for error factor specifications."
			<< std::endl;
		  abort_handler(-1);
		}
		else {             // mean, std deviation spec
		  const Real& mean  = ranVarMeansX[j];
		  const Real& stdev = ranVarStdDevsX[j];
		  const Real& zeta  = ranVarAddtlParamsX[j][1];
		  Real mean_sq = mean*mean, var = stdev*stdev;
		  jacobian_xs(j, i)
		    = x*stdev*(z_vars[j]-zeta)/zeta/(mean_sq+var);
		}
		break;
	      case LN_LAMBDA: // Deriv of Lognormal w.r.t. its Lambda
		jacobian_xs(j, i) = x;
		break;
	      case LN_ZETA:   // Deriv of Lognormal w.r.t. its Zeta
		jacobian_xs(j, i) = x*z_vars[j];
		break;
	      case LN_ERR_FACT: // Deriv of Lognormal w.r.t. its Error Factor
		if (ln_err_fact) // mean, error factor spec
		  jacobian_xs(j, i) = x*(z_vars[j] - ranVarAddtlParamsX[j][1])
		    / Phi_inverse(0.95) / ranVarAddtlParamsX[j][2];
		else {           // mean, std deviation spec
		  PCerr << "Error: derivative with respect to LN_ERR_FACT "
			<< "requires error factor specifications."
			<< std::endl;
		  abort_handler(-1);
		}
		break;
		//case LN_LWR_BND: case LN_UPR_BND: not supported
	      case NO_TARGET: default:
		PCerr << "Error: secondary mapping failure for LOGNORMAL in "
		      << "NatafTransformation::jacobian_dX_dS()." << std::endl;
		abort_handler(-1);
		break;
	      }
	    }
	    // Deriv of Lognormal w.r.t. any distribution parameter
	    if (correlationFlagX)
	      jacobian_xs(j, i) += x*ranVarAddtlParamsX[j][1]*num_dz_ds(j, i);
	  }
	  else if (ranVarTypesU[j] == LOGNORMAL) { // *** TO DO ***
	    PCerr << "Error: mapping not yet supported in NatafTransformation::"
		  << "jacobian_dX_dS()" << std::endl;
	    abort_handler(-1);
	  }
	  break;
	}
	case BOUNDED_LOGNORMAL: { // bounded lognormal
	  const Real& mean   = ranVarMeansX[j];
	  const Real& stdev  = ranVarStdDevsX[j];
	  const Real& lambda = ranVarAddtlParamsX[j][0];
	  const Real& zeta   = ranVarAddtlParamsX[j][1];
	  const Real& lwr    = ranVarLowerBndsX[j];
	  const Real& upr    = ranVarUpperBndsX[j];
	  const Real& x      = x_vars[j]; const Real& z = z_vars[j];
	  Real xms = (log(x)-lambda)/zeta, phi_xms = phi(xms);
	  Real lms = (lwr > 0.)      ? (log(lwr)-lambda)/zeta : -dbl_inf;
	  Real ums = (upr < dbl_inf) ? (log(upr)-lambda)/zeta :  dbl_inf;
	  if (ranVarTypesU[j] == STD_NORMAL) {
	    if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	      Real phi_lms = (lwr > 0.)      ? phi(lms) : 0.;
	      Real phi_ums = (upr < dbl_inf) ? phi(ums) : 0.;
	      Real dlambda_ds = 0., dzeta_ds = 0., dlwr_ds = 0., dupr_ds = 0.,
		mean_sq = mean*mean, var = stdev*stdev;
	      bool ln_err_fact = (ranVarAddtlParamsX[j].length() > 2);
	      switch (target2) {
	      case LN_MEAN: // Deriv of Bounded Lognormal w.r.t. its Mean
		if (ln_err_fact) // mean, error factor spec
		  dlambda_ds = 1./mean; //dzeta_ds = 0.;
		else {             // mean, std deviation spec
		  dlambda_ds = (1.+var/(mean_sq+var))/mean;
		  dzeta_ds   = -var/zeta/mean/(mean_sq+var);
		}
		break;
	      case LN_STD_DEV: // Deriv of Bounded LogN w.r.t. its Std Deviation
		if (ln_err_fact) {
		  PCerr << "Error: derivative with respect to LN_STD_DEV is "
			<< "unsupported for error factor specifications."
			<< std::endl;
		  abort_handler(-1);
		}
		else {
		  dlambda_ds = -stdev/(mean_sq+var);
		  dzeta_ds   = stdev/zeta/(mean_sq+var);
		}
		break;
	      case LN_LAMBDA: // Deriv of Bounded Lognormal w.r.t. its Lambda
		dlambda_ds = 1.; break;
	      case LN_ZETA:   // Deriv of Bounded Lognormal w.r.t. its Zeta
		dzeta_ds = 1.;   break;
	      case LN_ERR_FACT: // Deriv of Bounded LogN w.r.t. its Error Factor
		if (ln_err_fact) {
		  dzeta_ds   = 1./Phi_inverse(0.95)/ranVarAddtlParamsX[j][2];
		  dlambda_ds = -zeta*dzeta_ds;
		}
		else {
		  PCerr << "Error: derivative with respect to LN_ERR_FACT "
			<< "requires error factor specifications."
			<< std::endl;
		  abort_handler(-1);
		}
		break;
	      case LN_LWR_BND: // Deriv of Bounded LogN w.r.t. its Lower Bound
		dlwr_ds = 1.; break;
	      case LN_UPR_BND: // Deriv of Bounded LogN w.r.t. its Upper Bound
		dupr_ds = 1.; break;
	      case NO_TARGET: default:
		PCerr << "Error: secondary mapping failure for BOUNDED_"
		      << "LOGNORMAL in NatafTransformation::jacobian_dX_dS()."
		      << std::endl;
		abort_handler(-1);
		break;
	      }
	      Real dlms_ds = (lwr > 0.) ?
		(dlwr_ds/lwr - dlambda_ds - lms*dzeta_ds)/zeta : 0.;
	      Real dums_ds = (upr < dbl_inf) ?
		(dupr_ds/upr - dlambda_ds - ums*dzeta_ds)/zeta : 0.;
	      Real dxms_ds = Phi(z)/phi_xms*(phi_ums*dums_ds - phi_lms*dlms_ds)
		+ phi_lms/phi_xms*dlms_ds;
	      jacobian_xs(j, i) = x*(zeta*dxms_ds + dlambda_ds + xms*dzeta_ds);
	    }
	    // Deriv of Lognormal w.r.t. any distribution parameter
	    if (correlationFlagX) {
	      Real Phi_lms = (lwr > 0.)      ? Phi(lms) : 0.;
	      Real Phi_ums = (upr < dbl_inf) ? Phi(ums) : 1.;
	      jacobian_xs(j, i)
		+= (Phi_ums - Phi_lms)*phi(z)/phi_xms*num_dz_ds(j, i);
	    }
	  }
	  else if (ranVarTypesU[j] == BOUNDED_LOGNORMAL) { // *** TO DO ***
	    PCerr << "Error: mapping not yet supported in NatafTransformation::"
		  << "jacobian_dX_dS()" << std::endl;
	    abort_handler(-1);
	  }
	  break;
	}
	case UNIFORM: {
	  // to UNIFORM: x = L + (z + 1)*(U - L)/2
	  // to NORMAL:  x = L + Phi(z) (U - L)
	  const Real& z = z_vars[j];
	  if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	    switch (target2) {
	    case U_LWR_BND: // Deriv of Uniform w.r.t. its Lower Bound
	      if (ranVarTypesU[j] == STD_UNIFORM)
		jacobian_xs(j, i) = (1. - z)/2.;
	      else if (ranVarTypesU[j] == STD_NORMAL)
		jacobian_xs(j, i) = (z > 0.) ? Phi(-z) : 1. - Phi(z);
	      break;
	    case U_UPR_BND: // Deriv of Uniform w.r.t. its Upper Bound
	      if (ranVarTypesU[j] == STD_UNIFORM)
		jacobian_xs(j, i) = (z + 1.)/2.;
	      else if (ranVarTypesU[j] == STD_NORMAL)
		jacobian_xs(j, i) = Phi(z);
	      break;
	    // Uniform Mean          - TO DO
	    // Uniform Std Deviation - TO DO
	    case NO_TARGET: default:
	      PCerr << "Error: secondary mapping failure for UNIFORM in "
		    << "NatafTransformation::jacobian_dX_dS()." << std::endl;
	      abort_handler(-1);
	      break;
	    }
	  }
	  // Deriv of Uniform w.r.t. any distribution parameter
	  // Note: UNIFORM case should currently be zero, but may be
	  // nonzero in the future once correlation warping is more complete.
	  if (correlationFlagX) {
	    if (ranVarTypesU[j] == STD_UNIFORM)
	      jacobian_xs(j, i) +=
		(ranVarUpperBndsX[j]-ranVarLowerBndsX[j])/2.*num_dz_ds(j, i);
	    else if (ranVarTypesU[j] == STD_NORMAL)
	      jacobian_xs(j, i) +=
		(ranVarUpperBndsX[j]-ranVarLowerBndsX[j])*phi(z)*num_dz_ds(j,i);
	  }
	  break;
	}
	case LOGUNIFORM: {
	  // to UNIFORM: ln x = ln L + (z+1)/2 (ln U - ln L)
	  // to NORMAL:  ln x = ln L +  Phi(z) (ln U - ln L)
	  const Real& lwr = ranVarLowerBndsX[j];
	  const Real& upr = ranVarUpperBndsX[j];
	  const Real& x = x_vars[j]; const Real& z = z_vars[j];
	  if (ranVarTypesU[j] == LOGUNIFORM) { // no transformation
	    // *** TO DO ***
	    PCerr << "Error: mapping not yet supported in NatafTransformation::"
		  << "jacobian_dX_dS()" << std::endl;
	    abort_handler(-1);
	  }
	  else { // transformation to std normal or std uniform
	    if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	      switch (target2) {
	      case LU_LWR_BND: { // Deriv of Loguniform w.r.t. its Lower Bound
		if (ranVarTypesU[j] == STD_UNIFORM)
		  jacobian_xs(j, i) = x*(1.-z)/2./lwr;
		else if (ranVarTypesU[j] == STD_NORMAL) {
		  Real normcdf_comp = (z > 0.) ? Phi(-z) : 1. - Phi(z);
		  jacobian_xs(j, i) = x*normcdf_comp/lwr;
		}
		break;
	      }
	      case LU_UPR_BND: // Deriv of Loguniform w.r.t. its Upper Bound
		if (ranVarTypesU[j] == STD_UNIFORM)
		  jacobian_xs(j, i) = x*(z+1.)/2./upr;
		else if (ranVarTypesU[j] == STD_NORMAL)
		  jacobian_xs(j, i) = x*Phi(z)/upr;
		break;
	      // Loguniform Mean          - TO DO
	      // Loguniform Std Deviation - TO DO
	      case NO_TARGET: default:
		PCerr << "Error: secondary mapping failure for LOGUNIFORM in "
		      << "NatafTransformation::jacobian_dX_dS()." << std::endl;
		abort_handler(-1);
		break;
	      }
	    }
	    // Deriv of Loguniform w.r.t. any distribution parameter
	    // Note: UNIFORM case should currently be zero, but may be
	    // nonzero in the future once correlation warping is more complete.
	    if (correlationFlagX) { // not currently supported in Nataf
	      if (ranVarTypesU[j] == STD_UNIFORM)
		jacobian_xs(j, i) += x*(log(upr)-log(lwr))/2.*num_dz_ds(j, i);
	      else if (ranVarTypesU[j] == STD_NORMAL)
		jacobian_xs(j, i)
		  += x*(log(upr)-log(lwr))*phi(z)*num_dz_ds(j, i);
	    }
	  }
	  break;
	}
	case TRIANGULAR: {
	  const Real& lwr  = ranVarLowerBndsX[j];
	  const Real& mode = ranVarAddtlParamsX[j][0];
	  const Real& upr  = ranVarUpperBndsX[j];
	  const Real& x = x_vars[j]; const Real& z = z_vars[j];
	  if (ranVarTypesU[j] == TRIANGULAR) { // no transformation
	    // *** TO DO ***
	    PCerr << "Error: mapping not yet supported in NatafTransformation::"
		  << "jacobian_dX_dS()" << std::endl;
	    abort_handler(-1);
	  }
	  else { // transformation to std normal or std uniform
	    if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	      bool dist_error = false;
	      if (x < mode) {
		Real term;
		if (ranVarTypesU[j] == STD_UNIFORM)
		  term = (z+1.)/4.;
		else if (ranVarTypesU[j] == STD_NORMAL)
		  term = Phi(z)/2.;
		switch (target2) {
		case T_MODE: // Triangular Mode
		  jacobian_xs(j, i) = term*(upr-lwr)/(x-lwr);            break;
		case T_LWR_BND: // Triangular Lower Bound
		  jacobian_xs(j, i) = 1.+term*(2.*lwr-upr-mode)/(x-lwr); break;
		case T_UPR_BND: // Triangular Upper Bound
		  jacobian_xs(j, i) = term*(mode-lwr)/(x-lwr);           break;
		// Triangular Mean          - TO DO
		// Triangular Std Deviation - TO DO
		case NO_TARGET: default:
		  dist_error = true;                                     break;
		}
		// Deriv of Triangular w.r.t. any distribution parameter
		if (correlationFlagX) {// not currently supported for triangular
		  Real term_deriv;
		  if (ranVarTypesU[j] == STD_UNIFORM)
		    term_deriv = 0.25;
		  else if (ranVarTypesU[j] == STD_NORMAL)
		    term_deriv = phi(z)/2.;
		  jacobian_xs(j, i)
		    += (upr-lwr)*(mode-lwr)*term_deriv/(x-lwr)*num_dz_ds(j, i);
		}
	      }
	      else {
		Real term;
		if (ranVarTypesU[j] == STD_UNIFORM)
		  term = (1.-z)/4.;
		else if (ranVarTypesU[j] == STD_NORMAL)
		  term = (z > 0.) ? Phi(-z)/2. : (1. - Phi(z))/2.;
		switch (target2) {
		case T_MODE: // Triangular Mode
		  jacobian_xs(j, i) = term*(upr-lwr)/(upr-x);            break;
		case T_LWR_BND: // Triangular Lower Bound
		  jacobian_xs(j, i) = (upr-mode)*term/(upr-x);           break;
		case T_UPR_BND: // Triangular Upper Bound
		  jacobian_xs(j, i) = 1.-term*(2.*upr-lwr-mode)/(upr-x); break;
		// Triangular Mean          - TO DO
		// Triangular Std Deviation - TO DO
		case NO_TARGET: default:
		  dist_error = true;                                     break;
		}
		// Deriv of Triangular w.r.t. any distribution parameter
		if (correlationFlagX) {// not currently supported for triangular
		  Real term_deriv;
		  if (ranVarTypesU[j] == STD_UNIFORM)
		    term_deriv = -0.25;
		  else if (ranVarTypesU[j] == STD_NORMAL)
		    term_deriv = -phi(z)/2.;
		  jacobian_xs(j, i)
		    -= (upr-mode)*(upr-lwr)*term_deriv/(upr-x)*num_dz_ds(j, i);
		}
	      }
	      if (dist_error) {
		PCerr << "Error: secondary mapping failure for TRIANGULAR in "
		      << "NatafTransformation::jacobian_dX_dS()." << std::endl;
		abort_handler(-1);
	      }
	    }
	  }
	  break;
	}
	case EXPONENTIAL: {
	  // to EXPONENTIAL: x = beta*z
	  // to NORMAL:      Phi(z) = 1. - exp(-x/beta)
	  //                 x = -beta ln(1. - Phi(z))
	  const Real& z = z_vars[j];
	  Real normcdf_comp;
	  if (ranVarTypesU[j] == STD_NORMAL)
	    normcdf_comp = (z > 0.) ? Phi(-z) : 1. - Phi(z);
	  if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	    switch (target2) {
	    case E_BETA: // Exponential Beta
	      if (ranVarTypesU[j] == STD_EXPONENTIAL)
		jacobian_xs(j, i) = z;
	      else if (ranVarTypesU[j] == STD_NORMAL)
		jacobian_xs(j, i)
		  = (z > 0.) ? -log(normcdf_comp) : -bmth::log1p(-Phi(z));
	      break;
	    // Exponential Mean          - TO DO
	    // Exponential Std Deviation - TO DO
	    case NO_TARGET: default:
	      PCerr << "Error: secondary mapping failure for EXPONENTIAL in "
		    << "NatafTransformation::jacobian_dX_dS()." << std::endl;
	      abort_handler(-1);
	      break;
	    }
	  }
	  // Deriv of Exponential w.r.t. any distribution parameter
	  // Note: EXPONENTIAL case should currently be zero, but may be
	  // nonzero in the future once correlation warping is more complete.
	  if (correlationFlagX) {
	    const Real& beta = ranVarAddtlParamsX[j][0];
	    if (ranVarTypesU[j] == STD_EXPONENTIAL)
	      jacobian_xs(j, i) += beta*num_dz_ds(j, i);
	    else if (ranVarTypesU[j] == STD_NORMAL)
	      jacobian_xs(j, i) += beta*phi(z)/normcdf_comp*num_dz_ds(j, i);
	  }
	  break;
	}
	case BETA: {
	  if (ranVarTypesU[j] == STD_BETA && !need_xs) {
	    // x = lwr + (upr - lwr)*(z+1.)/2.;
	    if (j == cv_index) {//corresp var has deriv w.r.t. its dist param
	      switch (target2) {
	      //case BE_ALPHA: // numerically evaluated
	      //case BE_BETA:  // numerically evaluated
	      case BE_LWR_BND: // Beta Lower Bound
		jacobian_xs(j, i) = (1. - z_vars[j])/2.;
	      case BE_UPR_BND: // Beta Upper Bound
		jacobian_xs(j, i) = (z_vars[j] + 1.)/2.;
	      }
	    }
	    // Deriv of Uniform w.r.t. any distribution parameter
	    // Note: BETA case should currently be zero, but may be
	    // nonzero in the future once correlation warping is more complete.
	    if (correlationFlagX)
	      jacobian_xs(j, i)
		+= (ranVarUpperBndsX[j]-ranVarLowerBndsX[j])/2.*num_dz_ds(j, i);
	  }
	  break;
	}
	case GAMMA: {
	  if (ranVarTypesU[j] == STD_GAMMA && !need_xs) { // x = z*beta
	    if (j == cv_index) {//corresp var has deriv w.r.t. its dist param
	      switch (target2) {
	      //case GA_ALPHA: // numerically evaluated
	      case GA_BETA:    // Gamma Beta
		jacobian_xs(j, i) = z_vars[j];
	      }
	    }
	    // Deriv of Uniform w.r.t. any distribution parameter
	    // Note: GAMMA case should currently be zero, but may be
	    // nonzero in the future once correlation warping is more complete.
	    if (correlationFlagX)
	      jacobian_xs(j, i) += ranVarAddtlParamsX[j][1]*num_dz_ds(j, i);
	  }
	  break;
	}
	case GUMBEL: {
	  const Real& alpha = ranVarAddtlParamsX[j][0];
	  const Real& z = z_vars[j];
	  if (ranVarTypesU[j] == STD_NORMAL) {
	    if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	      switch (target2) {
	      case GU_ALPHA: { // Gumbel Alpha
		// x = beta - ln(-ln(Phi(z)))/alpha
		jacobian_xs(j, i) = log(-log(Phi(z)))/alpha/alpha;
		break;
	      }
	      case GU_BETA: // Gumbel Beta
		jacobian_xs(j, i) = 1.;
		break;
	      // Gumbel Mean
	        //x = x_vars[j];
		//alpha = PI/std::sqrt(6.)/ranVarStdDevsX[j];
		//num = -alpha*(x-z);
		//jacobian_xs(j, i) = -alpha*exp(num-exp(num))/phi(z);
		//break;
	      // Gumbel Standard Deviation
	        //x = x_vars[j];
		//alpha = PI/std::sqrt(6.)/ranVarStdDevsX[j];
		//num = -alpha*(x-z);
		//jacobian_xs(j, i)
		//  = num*exp(num-exp(num))/ranVarStdDevsX[j]/phi(z);
		//break;
	      case NO_TARGET: default:
		PCerr << "Error: secondary mapping failure for GUMBEL in "
		      << "NatafTransformation::jacobian_dX_dS()." << std::endl;
		abort_handler(-1);
		break;
	      }
	    }
	    // Deriv of Gumbel w.r.t. any distribution parameter
	    if (correlationFlagX) {
	      Real normcdf = Phi(z),
		lognormcdf = (z > 0.) ? bmth::log1p(-Phi(-z)) : log(normcdf);
	      jacobian_xs(j, i) -= phi(z)/alpha/normcdf/lognormcdf
		                *  num_dz_ds(j, i);
	    }
	  }
	  else if (ranVarTypesU[j] == GUMBEL) { // *** TO DO ***
	    PCerr << "Error: mapping not yet supported in NatafTransformation::"
		  << "jacobian_dX_dS()" << std::endl;
	    abort_handler(-1);
	  }
	  break;
	}
	case FRECHET: {
	  const Real& alpha = ranVarAddtlParamsX[j][0];
	  const Real& z = z_vars[j];
	  if (ranVarTypesU[j] == STD_NORMAL) {
	    if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	      switch (target2) {
	      case F_ALPHA: { // Frechet Alpha
		// x = beta (-ln(Phi(z)))^(-1/alpha)
		const Real& beta  = ranVarAddtlParamsX[j][1];
		Real num = -log(Phi(z));
		jacobian_xs(j, i) = beta/alpha/alpha*log(num)
		                  * std::pow(num,-1./alpha);
		break;
	      }
	      case F_BETA: // Frechet Beta
		jacobian_xs(j, i) = std::pow(-log(Phi(z)), -1./alpha);
		break;
	      // Frechet Mean          - TO DO
	      // Frechet Std Deviation - TO DO
	      case NO_TARGET: default:
		PCerr << "Error: secondary mapping failure for FRECHET in "
		      << "NatafTransformation::jacobian_dX_dS()." << std::endl;
		abort_handler(-1);
		break;
	      }
	    }
	    // Deriv of Frechet w.r.t. any distribution parameter
	    if (correlationFlagX) {
	      Real normcdf = Phi(z),
		lognormcdf = (z > 0.) ? bmth::log1p(-Phi(-z)) : log(normcdf);
	      jacobian_xs(j, i) -= x_vars[j]*phi(z)/alpha/normcdf/lognormcdf
		                *  num_dz_ds(j, i);
	    }
	  }
	  else if (ranVarTypesU[j] == FRECHET) { // *** TO DO ***
	    PCerr << "Error: mapping not yet supported in NatafTransformation::"
		  << "jacobian_dX_dS()" << std::endl;
	    abort_handler(-1);
	  }
	  break;
	}
	case WEIBULL: {
	  const Real& alpha = ranVarAddtlParamsX[j][0];
	  const Real& z = z_vars[j]; 
	  if (ranVarTypesU[j] == STD_NORMAL) {
	    if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	      Real log1mnormcdf
		= (z > 0.) ? log(Phi(-z)) : bmth::log1p(-Phi(z));
	      switch (target2) {
	      case W_ALPHA: { // Weibull Alpha
		// x = beta (-ln(1-Phi(z)))^(1/alpha)
		const Real& beta  = ranVarAddtlParamsX[j][1];
		jacobian_xs(j, i) = -beta/alpha/alpha*log(-log1mnormcdf)*
		  std::pow(-log1mnormcdf,1./alpha);
		break;
	      }
	      case W_BETA: // Weibull Beta
		jacobian_xs(j, i) = std::pow(-log1mnormcdf, 1./alpha);
		break;
	      // Weibull Mean          - TO DO
	      // Weibull Std Deviation - TO DO
	      case NO_TARGET: default:
		PCerr << "Error: secondary mapping failure for WEIBULL in "
		      << "NatafTransformation::jacobian_dX_dS()." << std::endl;
		abort_handler(-1);
		break;
	      }
	    }
	    // Deriv of Weibull w.r.t. any distribution parameter
	    if (correlationFlagX) {
	      Real normcdf_comp = (z > 0.) ? Phi(-z) : 1. - Phi(z);
	      Real log1mnormcdf
		= (z > 0.) ? log(normcdf_comp) : bmth::log1p(-Phi(z));
	      jacobian_xs(j, i) -= x_vars[j]*phi(z)/alpha/normcdf_comp/
		log1mnormcdf * num_dz_ds(j,i);
	    }
	  }
	  else if (ranVarTypesU[j] == WEIBULL) { // *** TO DO ***
	    PCerr << "Error: mapping not yet supported in NatafTransformation::"
		  << "jacobian_dX_dS()" << std::endl;
	    abort_handler(-1);
	  }
	  break;
	}
	case HISTOGRAM_BIN:
	  // TO DO
	  PCerr << "Error: unsupported mapping for HISTOGRAM_BIN in "
		<< "NatafTransformation::jacobian_dX_dS()." << std::endl;
	  abort_handler(-1);
	  break;
	case CONTINUOUS_INTERVAL:
	  // TO DO
	  PCerr << "Error: unsupported mapping for CONTINUOUS_INTERVAL in "
		<< "NatafTransformation::jacobian_dX_dS()." << std::endl;
	  abort_handler(-1);
	  break;
	case CONTINUOUS_STATE: { // x = L + (z + 1)*(U - L)/2
	  if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	    switch (target2) {
	    case CSV_LWR_BND: // Deriv of CSV w.r.t. its Lower Bound
	      jacobian_xs(j, i) = (1. - z_vars[j])/2.; break;
	    case CSV_UPR_BND: // Deriv of CSV w.r.t. its Upper Bound
	      jacobian_xs(j, i) = (z_vars[j] + 1.)/2.; break;
	    case NO_TARGET:   // can occur for all_variables Jacobians
	      jacobian_xs(j, i) = 0.;                  break;
	    default:
	      PCerr << "Error: secondary mapping failure for CONTINUOUS_STATE "
		    << "in NatafTransformation::jacobian_dX_dS()." << std::endl;
	      abort_handler(-1);
	      break;
	    }
	  }
	  // Deriv of CSV w.r.t. other distribution params should always be 0
	  //if (correlationFlagX)
	  //  jacobian_xs(j, i) += 
	  //    (ranVarUpperBndsX[j]-ranVarLowerBndsX[j])/2.*num_dz_ds(j, i);
	  break;
	}
	}
      }
    }
  }
}


/** This procedure computes the Hessian of the transformation x(u).
    hessian_xu is a 3D tensor modeled as an array of matrices, where
    the i_th matrix is d^2X_i/dU^2.  x_vars is the vector of random
    variables in the original user-defined x-space. */
void NatafTransformation::
hessian_d2X_dU2(const RealVector& x_vars, RealSymMatrixArray& hessian_xu)
{
  if (correlationFlagX) {
    // d^2X/dZ^2
    int x_len = x_vars.length();
    RealSymMatrixArray hessian_xz(x_len);
    hessian_d2X_dZ2(x_vars, hessian_xz);

    if (hessian_xu.size() != x_len)
      hessian_xu.resize(x_len);
    for (int i=0; i<x_len; i++) {
      // d^2X/dU^2 = dX/dZ^T d^2Z/dU^2 + dZ/dU^T d^2X/dZ^2 dZ/dU
      //           = L^T d^2X/dZ^2 L
      if (hessian_xu[i].numRows() != x_len)
	hessian_xu[i].shape(x_len);
      Teuchos::symMatTripleProduct(Teuchos::TRANS, 1., hessian_xz[i],
                                   corrCholeskyFactorZ, hessian_xu[i]);
    }
  }
  else // d^2X/dU^2 = d^2X/dZ^2 since dZ/dU = I
    hessian_d2X_dZ2(x_vars, hessian_xu);
}


/** This procedure computes the Hessian of the transformation x(z).
    hessian_xz is a 3D tensor modeled as an array of matrices, where
    the i_th matrix is d^2X_i/dZ^2.  x_vars is the vector of random
    variables in the original user-defined x-space. */
void NatafTransformation::
hessian_d2X_dZ2(const RealVector& x_vars, RealSymMatrixArray& hessian_xz)
{
  // This routine calculates the Hessian of the transformation x(z):
  //
  // For z type is STD_NORMAL:
  // d^2x/dz^2 = d/dz (dx/dz) = d/dz ( phi(z)/f(x) )
  //           = (f(x) phi'(z) - phi(z) f'(x) dx/dz)/f(x)^2
  //           = -phi(z)/f(x)^2 (z f(x) + f'(x) dx/dz) [since phi'(z)=-z phi(z)]
  //           = -dx/dz (z + f'(x)/f(x) dx/dz)
  //
  // For z type is STD_UNIFORM:
  // d^2x/dz^2 = d/dz (dx/dz) = d/dz ( u(z)/f(x) )
  //           = (f(x) u'(z) - u(z) f'(x) dx/dz)/f(x)^2
  //           = -dx/dz u(z)/f(x) f'(x)/f(x)           [since u'(z)=0]
  //           = -dx/dz dx/dz f'(x)/f(x)
  //
  // This requires the additional calculation of f'(x), the derivative of
  // the PDF.  Since GSL does not provide these, PDFs are differentiated
  // below in closed form.  For cases with f'(x) = 0 (e.g., uniform), the
  // expression can be simplified to d^2x/dz^2 = -z dx/dz

  int x_len = x_vars.length();
  if (hessian_xz.size() != x_len)
    hessian_xz.resize(x_len);

  RealVector z_vars;
  bool need_z = false;
  for (size_t i=0; i<x_len; i++)
    if (ranVarTypesU[i] == STD_NORMAL &&
	ranVarTypesX[i] != NORMAL && ranVarTypesX[i] != LOGNORMAL)
      { need_z = true; break; }
  if (need_z)
    trans_X_to_Z(x_vars, z_vars);

  for (int i=0; i<x_len; i++) {
    bool err_flag = false;
    if (hessian_xz[i].numRows() != x_len)
      hessian_xz[i].shape(x_len);
    // each d^2X_i/dZ^2 has a single entry on the diagonal as defined by
    // differentiation of jacobian_dX_dZ()
    switch (ranVarTypesX[i]) {
    case CONTINUOUS_DESIGN: case CONTINUOUS_STATE: case CONTINUOUS_INTERVAL:
      if (ranVarTypesU[i] == STD_UNIFORM)
	hessian_xz[i](i, i) = 0.;
      else
	err_flag = true;
      break;
    case NORMAL: // unbounded normal: z = (x - mean)/stdev
      if (ranVarTypesU[i] == STD_NORMAL)
	hessian_xz[i](i, i) = 0.;
      else
	err_flag = true;
      break;
    case BOUNDED_NORMAL: // bounded normal
      if (ranVarTypesU[i] == BOUNDED_NORMAL) // Golub-Welsch: no transform
	hessian_xz[i](i, i) = 0.;
      else { // nonlinear transformation
	const Real& mean  = ranVarMeansX[i];
	const Real& stdev = ranVarStdDevsX[i];
	const Real& z = z_vars[i]; const Real& x = x_vars[i];
	Real pdf = bounded_normal_pdf(x, mean, stdev, ranVarLowerBndsX[i],
				      ranVarUpperBndsX[i]),
	  pdf_deriv = pdf*(mean-x)/stdev/stdev;
	if (ranVarTypesU[i] == STD_NORMAL) {
	  Real dx_dz = phi(z) / pdf;
	  hessian_xz[i](i, i) = -dx_dz*(z + pdf_deriv*dx_dz/pdf);
	}
	else if (ranVarTypesU[i] == STD_UNIFORM) {
	  Real dx_dz = std_uniform_pdf() / pdf;
	  hessian_xz[i](i, i) = -dx_dz*dx_dz*pdf_deriv/pdf;
	}
	else
	  err_flag = true;
      }
      break;
    case LOGNORMAL: // unbounded lognormal: z = (ln x - lambda)/zeta
      if (ranVarTypesU[i] == STD_NORMAL) {
	// dx/dz = zeta x
	// d^2x/dz^2 = zeta dx/dz = zeta^2 x
	Real zeta_sq;
	lognormal_zeta_sq_from_moments(ranVarMeansX[i], ranVarStdDevsX[i],
				       zeta_sq);
	hessian_xz[i](i, i) = zeta_sq*x_vars[i];
	// causes DIFFs due to loss of precision from squaring of sqrt:
	//const Real& zeta    = ranVarAddtlParamsX[i][1];
	//hessian_xz[i](i, i) = zeta*zeta*x_vars[i];
      }
      else if (ranVarTypesU[i] == LOGNORMAL) // Golub-Welsch: no transform
	hessian_xz[i](i, i) = 0.;
      else
	err_flag = true;
      break;
    case BOUNDED_LOGNORMAL: // bounded lognormal
      if (ranVarTypesU[i] == BOUNDED_LOGNORMAL) // Golub-Welsch: no x-form
	hessian_xz[i](i, i) = 0.;
      else { // nonlinear transformation
	const Real& lambda = ranVarAddtlParamsX[i][0];
	const Real& zeta   = ranVarAddtlParamsX[i][1];
	const Real& mean   = ranVarMeansX[i];
	const Real& stdev  = ranVarStdDevsX[i];
	const Real& x = x_vars[i]; const Real& z = z_vars[i];
	Real xms = (log(x)-lambda)/zeta,
	     pdf = bounded_lognormal_pdf(x, mean, stdev, ranVarLowerBndsX[i],
					 ranVarUpperBndsX[i]),
	     pdf_deriv = -pdf*(zeta+xms)/x/zeta;
	if (ranVarTypesU[i] == STD_NORMAL) {
	  Real dx_dz = phi(z) / pdf;
	  hessian_xz[i](i, i) = -dx_dz*(z + pdf_deriv*dx_dz/pdf);
	}
	else if (ranVarTypesU[i] == STD_UNIFORM) {
	  Real dx_dz = std_uniform_pdf() / pdf;
	  hessian_xz[i](i, i) = -dx_dz*dx_dz*pdf_deriv/pdf;
	}
        else
	  err_flag = true;
      }
      break;
    case UNIFORM:
      //  F(x) = (x-L)/(U-L)
      //  f(x) = 1/(U-L)
      // f'(x) = 0.
      if (ranVarTypesU[i] == STD_UNIFORM) // linear scaling
	hessian_xz[i](i, i) = 0.;
      else if (ranVarTypesU[i] == STD_NORMAL) { // nonlinear transformation
	const Real& z = z_vars[i];
	hessian_xz[i](i, i) = -z*phi(z) /
	  uniform_pdf(ranVarLowerBndsX[i], ranVarUpperBndsX[i]);
      }
      else
	err_flag = true;
      break;
    case LOGUNIFORM:
      if (ranVarTypesU[i] == LOGUNIFORM) // Golub-Welsch: no transform
	hessian_xz[i](i, i) = 0.;
      else {
	//  F(x) = (ln x - ln L)/(ln U - ln L)
	//  f(x) =  1/(ln U - ln L)/x
	// f'(x) = -1/(ln U - ln L)/x^2
	if (ranVarTypesU[i] == STD_UNIFORM) {
	  Real log_range = log(ranVarUpperBndsX[i]) - log(ranVarLowerBndsX[i]);
	  hessian_xz[i](i, i) = x_vars[i]*log_range*log_range/4.;
	}
	else if (ranVarTypesU[i] == STD_NORMAL) {
	  const Real& z = z_vars[i]; const Real& x = x_vars[i];
	  Real pdf =
	    loguniform_pdf(x, ranVarLowerBndsX[i], ranVarUpperBndsX[i]),
	    pdf_deriv = -pdf/x, dx_dz = phi(z)/pdf;
	  hessian_xz[i](i, i) = -dx_dz*(z + pdf_deriv*dx_dz/pdf);
	}
	else
	  err_flag = true;
      }
      break;
    case TRIANGULAR:
      if (ranVarTypesU[i] == TRIANGULAR) // Golub-Welsch: no transform
	hessian_xz[i](i, i) = 0.;
      else {
	//             x < M                           x > M
	//  F(x): (x-L)^2/(U-L)/(M-L)    (M-L)/(U-L) - (x+M-2U)(x-M)/(U-L)/(U-M)
	//  f(x): 2(x-L)/(U-L)/(M-L)     2(U-x)/(U-L)/(U-M)
	// f'(x): 2/(U-L)/(M-L)          -2/(U-L)/(U-M)
	// Note: at x=M, F(x) and f(x) are continuous but f'(x) is not
	const Real& lwr  = ranVarLowerBndsX[i];
	const Real& mode = ranVarAddtlParamsX[i][0];
	const Real& upr  = ranVarUpperBndsX[i];
	const Real& x    = x_vars[i];
	Real pdf, pdf_deriv, range = upr - lwr;
	if (x < mode) {
	  pdf_deriv = 2./range/(mode-lwr);
	  pdf = (x-lwr)*pdf_deriv;
	}
	else if (x > mode) {
	  pdf_deriv = -2./range/(upr-mode);
	  pdf = (x-upr)*pdf_deriv;
	}
	else { // x == mode
	  pdf_deriv = 0.; // f'(x) is undefined: use 0.
	  pdf = 2./range;
	}
	if (ranVarTypesU[i] == STD_UNIFORM)
	  hessian_xz[i](i, i) = -pdf_deriv/4./std::pow(pdf, 3);
	else if (ranVarTypesU[i] == STD_NORMAL) {
	  const Real& z = z_vars[i];
	  Real dx_dz = phi(z)/pdf;
	  hessian_xz[i](i, i) = -dx_dz*(z + pdf_deriv*dx_dz/pdf);
	}
	else
	  err_flag = true;
      }
      break;
    case EXPONENTIAL:
      //  F(x) = 1. - e^(-x/beta)
      //  f(x) = e^(-x/beta) / beta
      // f'(x) = - e^(-x/beta) / beta^2
      if (ranVarTypesU[i] == STD_EXPONENTIAL) // linear scaling
	hessian_xz[i](i, i) = 0.;
      else if (ranVarTypesU[i] == STD_NORMAL) { // nonlinear transformation
	const Real& beta = ranVarAddtlParamsX[i][0]; const Real& z = z_vars[i];
	Real pdf = exponential_pdf(x_vars[i], beta), pdf_deriv = -pdf/beta,
	  dx_dz = phi(z)/pdf;
	hessian_xz[i](i, i) = -dx_dz*(z + pdf_deriv*dx_dz/pdf);
      }
      else
	err_flag = true;
      break;
    case BETA:
      //  F(x) = gsl
      //  f(x) = gsl
      // f'(x) = f(x) ((alpha-1)/(x-lwr) - (beta-1)/(upr-x))
      if (ranVarTypesU[i] == STD_BETA) // linear scaling
	hessian_xz[i](i, i) = 0.;
      else { // nonlinear transformation
	const Real& alpha = ranVarAddtlParamsX[i][0];
	const Real& beta  = ranVarAddtlParamsX[i][1];
	const Real& lwr   = ranVarLowerBndsX[i];
	const Real& upr   = ranVarUpperBndsX[i];
	const Real& z = z_vars[i]; const Real& x = x_vars[i];
	Real pdf = beta_pdf(x, alpha, beta, lwr, upr),
	  pdf_deriv = pdf*((alpha-1.)/(x-lwr) - (beta-1.)/(upr-x));
	if (ranVarTypesU[i] == STD_NORMAL) {
	  Real dx_dz = phi(z) / pdf;
	  hessian_xz[i](i, i) = -dx_dz*(z + pdf_deriv*dx_dz/pdf);
	}
	else if (ranVarTypesU[i] == STD_UNIFORM) {
	  Real dx_dz = std_uniform_pdf() / pdf;
	  hessian_xz[i](i, i) = -dx_dz*dx_dz*pdf_deriv/pdf;
	}
	else
	  err_flag = true;
      }
      break;
    case GAMMA:
      //  F(x) = gsl
      //  f(x) = beta^(-alpha) x^(alpha-1) e^(-x/beta) / GammaFn(alpha)
      // f'(x) = beta^(-alpha)/GammaFn(alpha) (e^(-x/beta) (alpha-1) x^(alpha-2)
      //                                       - x^(alpha-1) e^(-x/beta)/beta)
      if (ranVarTypesU[i] == STD_GAMMA) // linear scaling
	hessian_xz[i](i, i) = 0.;
      else if (ranVarTypesU[i] == STD_NORMAL) { // nonlinear transformation
	const Real& alpha = ranVarAddtlParamsX[i][0];
	const Real& beta  = ranVarAddtlParamsX[i][1];
	const Real& z = z_vars[i]; const Real& x = x_vars[i];
	Real pdf = gamma_pdf(x, alpha, beta),
	  pdf_deriv = std::pow(beta,-alpha) / gamma_function(alpha) *
	    (std::exp(-x/beta)*(alpha-1.) *
	     std::pow(x,alpha-2.) - std::pow(x,alpha-1.) *
	     std::exp(-x/beta)/beta), //gamma_pdf_deriv(x, alpha, beta),
	  dx_dz = phi(z) / pdf;
	hessian_xz[i](i, i) = -dx_dz*(z + pdf_deriv*dx_dz/pdf);
      }
      else
	err_flag = true;
      break;
    case GUMBEL:
      if (ranVarTypesU[i] == STD_NORMAL) {
	//  F(x) = e^(-e^(-alpha*(x-u)))
	//  f(x) = alpha e^(-alpha*(x-u)) F(x)
	// f'(x) = alpha (e^(-alpha*(x-u)) f(x) - alpha F(x) e^(-alpha*(x-u)))
	const Real& alpha = ranVarAddtlParamsX[i][0];
	const Real& beta  = ranVarAddtlParamsX[i][1];
	const Real& z     = z_vars[i];
	Real num = std::exp(-alpha*(x_vars[i]-beta)), cdf = std::exp(-num),
	  pdf = alpha*num*cdf, pdf_deriv = alpha*(num*pdf - alpha*cdf*num),
	  dx_dz = phi(z)/pdf;
	hessian_xz[i](i, i) = -dx_dz*(z + pdf_deriv*dx_dz/pdf);
      }
      else if (ranVarTypesU[i] == GUMBEL) // Golub-Welsch: no transform
	hessian_xz[i](i, i) = 0.;
      else
	err_flag = true;
      break;
    case FRECHET:
      if (ranVarTypesU[i] == STD_NORMAL) {
	//  F(x) = e^(-(beta/x)^alpha)
	//  f(x) = F(x) alpha (beta/x)^(alpha-1) beta/x^2
	//       = F(x) alpha/beta (beta/x)^(alpha+1)
	// f'(x) = alpha/beta ((beta/x)^(alpha+1) f(x) -
	//                     F(x) (alpha+1)/beta (beta/x)^(alpha+2))
	const Real& alpha = ranVarAddtlParamsX[i][0];
	const Real& beta  = ranVarAddtlParamsX[i][1];
	const Real& z     = z_vars[i];
	Real num = beta/x_vars[i], cdf = std::exp(-std::pow(num,alpha)),
	  pdf = alpha/beta*std::pow(num,alpha+1.)*cdf, dx_dz = phi(z)/pdf,
	  pdf_deriv = alpha/beta*(std::pow(num,alpha+1.)*pdf - 
				  cdf*(alpha+1.)/beta * std::pow(num,alpha+2.));
	hessian_xz[i](i, i) = -dx_dz*(z + pdf_deriv*dx_dz/pdf);
      }
      else if (ranVarTypesU[i] == FRECHET) // Golub-Welsch: no transform
	hessian_xz[i](i, i) = 0.;
      else
	err_flag = true;
      break;
    case WEIBULL:
      if (ranVarTypesU[i] == STD_NORMAL) {
	//  F(x) = 1.-e^(-(x/beta)^alpha)
	//  f(x) = alpha/beta e^(-(x/beta)^alpha) (x/beta)^(alpha-1)
	// f'(x) = alpha/beta (e^(-(x/beta)^alpha) (alpha-1)/beta
	//                     (x/beta)^(alpha-2) - (x/beta)^(alpha-1) f(x))
	const Real& alpha = ranVarAddtlParamsX[i][0];
	const Real& beta  = ranVarAddtlParamsX[i][1];
	const Real& z     = z_vars[i];
	Real num = x_vars[i]/beta, num2 = std::exp(-std::pow(num,alpha)),
	  pdf = alpha/beta*num2*std::pow(num,alpha-1.), dx_dz = phi(z)/pdf,
	  pdf_deriv = alpha/beta*(num2*(alpha-1.)/beta*std::pow(num,alpha-2.) -
				  std::pow(num,alpha-1.)*pdf);
	//Real cdf = 1.-num2;
	hessian_xz[i](i, i) = -dx_dz*(z + pdf_deriv*dx_dz/pdf);
      }
      else if (ranVarTypesU[i] == WEIBULL) // Golub-Welsch: no transform
	hessian_xz[i](i, i) = 0.;
      else
	err_flag = true;
      break;
    case HISTOGRAM_BIN:
      if (ranVarTypesU[i] == STD_NORMAL) {
	const Real& z = z_vars[i];
	hessian_xz[i](i, i) = -phi(z) * z /
	  histogram_bin_pdf(x_vars[i], ranVarAddtlParamsX[i]); // pdf_deriv = 0.
      }
      else if (ranVarTypesU[i] == STD_UNIFORM) // transformation is piece-wise
	hessian_xz[i](i, i) = 0.;              // linear with discont. Jacobian
      else if (ranVarTypesU[i] == HISTOGRAM_BIN)
	hessian_xz[i](i, i) = 0.;
      else
	err_flag = true;
      break;
    }
    if (err_flag) {
      PCerr << "Error: unsupported variable mapping for variable " << i
	    << " in NatafTransformation::hessian_d2X_dZ2()" << std::endl;
      abort_handler(-1);
    }
  }
}

} // namespace Pecos

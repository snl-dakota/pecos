/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:	 NatafTransformation
//- Description: Base class for nonlinear distribution transformations
//- Owner:       Mike Eldred
//- Checked by:
//- Version:

//#include "Teuchos_SerialSpdDenseSolver.h"
#ifdef PECOS_GSL
#include "gsl/gsl_sf_gamma.h"
#endif

static const char rcsId[]="@(#) $Id: NatafTransformation.C 4768 2007-12-17 17:49:32Z mseldre $";

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
}


/** This procedure computes the transformation from u to z space.
    u_vars is the vector of random variables in uncorrelated standard
    normal space (u-space).  z_vars is the vector of random variables
    in normal space with proper correlations (z-space). */
void NatafTransformation::
trans_U_to_Z(const RealVector& u_vars, RealVector& z_vars)
{
  // corrCholeskyFactorZ: the Cholesky factor of the modified correlation matrix

  int u_len = u_vars.Length();
  if (z_vars.Length() != u_len)
    z_vars.Size(u_len);
  z_vars.Multiply('N', 'N', 1., corrCholeskyFactorZ, u_vars, 0.);
}


/** This procedure computes the transformation from z to x space.
    z_vars is the vector of random variables in normal space with
    proper correlations (z-space).  x_vars is the vector of random
    variables in the original user-defined x-space */
void NatafTransformation::
trans_Z_to_X(const RealVector& z_vars, RealVector& x_vars)
{
  // ranVarTypesX:   the type of random variable (NORMAL/LOGNORMAL/UNIFORM/
  //                 LOGUNIFORM/TRIANGULAR/EXPONENTIAL/BETA/GAMMA/GUMBEL/
  //                 FRECHET/WEIBULL)
  // ranVarMeansX:   means of the random variables in x-space
  // ranVarStdDevsX: standard deviations of the random variables in x-space

  // This routine performs an inverse Rackwitz-Fiessler transformation:
  //   X = F^{-1}[ Phi(Z) ]

  int z_len = z_vars.Length();
  if (x_vars.Length() != z_len)
    x_vars.Size(z_len);

  for (int i=0; i<z_len; i++) {
    switch (ranVarTypesX[i]) {
    case DESIGN: case STATE: {
      // scale from [-1,1] to [L,U]
      const Real& lwr = ranVarLowerBndsX(i);
      x_vars(i) = lwr + (ranVarUpperBndsX(i) - lwr)*(z_vars(i)+1.)/2.;
      break;
    }
    case NORMAL: {
      // z = (x - mu)/sigma
      x_vars(i) = z_vars(i)*ranVarStdDevsX(i) + ranVarMeansX(i);
      break;
    }
    case LOGNORMAL: {
      // z = (ln x - lamba)/zeta
      const Real& mean = ranVarMeansX(i);
      Real cf_var = ranVarStdDevsX(i)/mean, zeta_sq = log(1 + cf_var*cf_var),
	lambda = log(mean) - zeta_sq/2.;
      x_vars(i) = exp(lambda + sqrt(zeta_sq)*z_vars(i));
      break;
    }
    case UNIFORM: {
      const Real& lwr = ranVarLowerBndsX(i);
      Real range = ranVarUpperBndsX(i) - lwr;
      if (extendedUSpace) // scale from [-1,1] to [L,U]
	x_vars(i) = lwr + range*(z_vars(i)+1.)/2.;
      else { // transform from std normal
	// Phi(z) = (x-L)/(U-L)
	Real normcdf = Phi(z_vars(i));
	x_vars(i) = normcdf * range + lwr;
      }
      break;
    }
    case LOGUNIFORM: {
      const Real& lwr = ranVarLowerBndsX(i);
      Real log_lwr = log(lwr), log_range = log(ranVarUpperBndsX(i)) - log_lwr;
      if (extendedUSpace) // transform from uniform on [-1,1]
	x_vars(i) = lwr*exp((z_vars(i)+1.)*log_range/2.);
      else { // transform from std normal
	// Phi(z) = (ln x - ln L)/(ln U - ln L)
	Real normcdf = Phi(z_vars(i));
	x_vars(i) = exp(normcdf * log_range + log_lwr);
      }
      break;
    }
    case TRIANGULAR: {
      const Real& lwr  = ranVarLowerBndsX(i);
      const Real& mode = ranVarAddtlParamsX[i](0);
      const Real& upr  = ranVarUpperBndsX(i);
      Real normcdf = Phi(z_vars(i)), range = upr - lwr;
      // assume x < mode and then check
      Real x = lwr + sqrt(normcdf*range*(mode-lwr));
      Real x_pdf = 2.*(x-lwr)/range/(mode-lwr), m_pdf = 2./range;
      // check pdf value to ensure that proper equation used
      if ( x_pdf > m_pdf )
	x = upr - sqrt((1.-normcdf)*range*(upr-mode));
      x_vars(i) = x;
      break;
    }
    case EXPONENTIAL: {
      const Real& beta = ranVarAddtlParamsX[i](0);
      if (extendedUSpace) // scale from std exponential exp(-z)
	x_vars(i) = beta*z_vars(i); // exp(-x/beta) = exp(-z)
      else // transform from std normal
	x_vars(i) = -beta*log(1.-Phi(z_vars(i)));
      break;
    }
    case BETA: { // scale from std beta on [-1,1]
      const Real& lwr = ranVarLowerBndsX(i);
      const Real& upr = ranVarUpperBndsX(i);
      if (extendedUSpace) // scale from [-1,1] to [L,U]
	x_vars(i) = lwr + (upr - lwr)*(z_vars(i)+1.)/2.;
#ifdef PECOS_GSL
      else { // transform from std normal
	const Real& alpha = ranVarAddtlParamsX[i](0);
	const Real& beta  = ranVarAddtlParamsX[i](1);
	Real normcdf = Phi(z_vars(i));
	// GSL does not support beta CDF inverse, therefore a special Newton
	// solve has been implemented for the inversion (which uses the GSL
	// CDF/PDF fns)
	Real scaled_x = cdf_beta_Pinv(normcdf, alpha, beta);
	x_vars(i) = lwr + (upr - lwr)*scaled_x;
      }
#endif // PECOS_GSL
      break;
    }
    case GAMMA: {
      const Real& beta = ranVarAddtlParamsX[i](1);
      if (extendedUSpace)
	x_vars(i) = beta*z_vars(i); // x/beta = z
#ifdef PECOS_GSL
      else {
	const Real& alpha = ranVarAddtlParamsX[i](0);
	Real normcdf = Phi(z_vars(i));
	// GSL gamma passes alpha followed by beta
	x_vars(i) = gsl_cdf_gamma_Pinv(normcdf, alpha, beta);
      }
#endif // PECOS_GSL
      break;
    }
    case GUMBEL: {
      const Real& alpha = ranVarAddtlParamsX[i](0);
      const Real& beta  = ranVarAddtlParamsX[i](1);
      Real normcdf = Phi(z_vars(i));
      x_vars(i) = beta - log(-log(normcdf))/alpha;
      break;
    }
    case FRECHET: {
      const Real& alpha = ranVarAddtlParamsX[i](0);
      const Real& beta  = ranVarAddtlParamsX[i](1);
      Real normcdf = Phi(z_vars(i));
      x_vars(i) = beta*pow(-log(normcdf),-1./alpha);
      break;
    }
#ifdef PECOS_GSL
    case WEIBULL: {
      const Real& alpha = ranVarAddtlParamsX[i](0);
      const Real& beta  = ranVarAddtlParamsX[i](1);
      Real normcdf = Phi(z_vars(i));
      // GSL weibull passes beta followed by alpha
      x_vars(i)	= gsl_cdf_weibull_Pinv(normcdf, beta, alpha);
      break;
    }
#endif // PECOS_GSL
    }
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
}


/** This procedure performs the transformation from x to z space:
    z_vars is the vector of random variables in normal space with
    proper correlations (z-space).  x_vars is the vector of random
    variables in the original user-defined x-space. */
void NatafTransformation::
trans_X_to_Z(const RealVector& x_vars, RealVector& z_vars)
{
  // ranVarTypesX:   the type of random variable (NORMAL/LOGNORMAL/UNIFORM/
  //                 LOGUNIFORM/TRIANGULAR/EXPONENTIAL/BETA/GAMMA/GUMBEL/
  //                 FRECHET/WEIBULL)
  // ranVarMeansX:   means of the random variables in x-space
  // ranVarStdDevsX: standard deviations of the random variables in x-space

  // This routine performs a Rackwitz-Fiessler transformation:
  //   Z = Phi^{-1}[ F(X) ]

  int x_len = x_vars.Length();
  if (z_vars.Length() != x_len)
    z_vars.Size(x_len);

  for (int i=0; i<x_len; i++) {
    switch (ranVarTypesX[i]) {
    case DESIGN: case STATE: {
      const Real& lwr = ranVarLowerBndsX(i);
      z_vars(i) = 2.*(x_vars(i) - lwr)/(ranVarUpperBndsX(i) - lwr) - 1.;
      break;
    }
    case NORMAL: {
      // z = (x - mu)/sigma
      z_vars(i)	= (x_vars(i)-ranVarMeansX(i))/ranVarStdDevsX(i);
      break;
    }
    case LOGNORMAL: {
      // z = (ln x - lamba)/zeta
      const Real& mean = ranVarMeansX(i);
      Real cf_var = ranVarStdDevsX(i)/mean, zeta_sq = log(1. + cf_var*cf_var),
	lambda = log(mean) - zeta_sq/2.;
      z_vars(i) = (log(x_vars(i)) - lambda)/sqrt(zeta_sq);
      break;
    }
    case UNIFORM: {
      // Phi(z) = (x-L)/(U-L)
      const Real& lwr = ranVarLowerBndsX(i);
      if (extendedUSpace) // scale to uniform on [-1,1]
	z_vars(i) = 2.*(x_vars(i) - lwr)/(ranVarUpperBndsX(i) - lwr) - 1.;
      else { // transform to std normal
	Real cdf  = (x_vars(i) - lwr)/(ranVarUpperBndsX(i) - lwr);
	z_vars(i) = Phi_inverse(cdf);
      }
      break;
    }
    case LOGUNIFORM: {
      Real log_lwr = log(ranVarLowerBndsX(i)),
	   log_upr = log(ranVarUpperBndsX(i));
      if (extendedUSpace) // transform to uniform on [-1,1]
	z_vars(i) = 2.*(log(x_vars(i)) - log_lwr)/(log_upr - log_lwr) - 1.;
      else { // transform to std normal
	// Phi(z) = (ln x - ln L)/(ln U - ln L)
	Real cdf = (log(x_vars(i)) - log_lwr)/(log_upr - log_lwr);
	z_vars(i) = Phi_inverse(cdf);
      }
      break;
    }
    case TRIANGULAR: {
      const Real& lwr  = ranVarLowerBndsX(i);
      const Real& mode = ranVarAddtlParamsX[i](0);
      const Real& upr  = ranVarUpperBndsX(i);
      const Real& x    = x_vars(i);
      Real range       = upr - lwr;
      Real cdf = (x < mode) ? pow(x-lwr,2.)/range/(mode-lwr) :
	((mode-lwr) - (x+mode-2*upr)*(x-mode)/(upr-mode))/range;
      z_vars(i) = Phi_inverse(cdf);
      break;
    }
    case EXPONENTIAL: {
      const Real& beta = ranVarAddtlParamsX[i](0);
      if (extendedUSpace) // scale to std exponential exp(-z)
	z_vars(i) = x_vars(i)/beta;
      else { // transform to std normal
	Real cdf  = 1. - exp(-x_vars(i)/beta);
	z_vars(i) = Phi_inverse(cdf);
      }
      break;
    }
    case BETA: {
      const Real& lwr = ranVarLowerBndsX(i);
      const Real& upr = ranVarUpperBndsX(i);
      if (extendedUSpace) // scale to beta on [-1,1]
	z_vars(i) = 2.*(x_vars(i) - lwr)/(upr - lwr) - 1.;
#ifdef PECOS_GSL
      else { // transform to std normal
	const Real& alpha = ranVarAddtlParamsX[i](0);
	const Real& beta  = ranVarAddtlParamsX[i](1);
	Real scaled_x = (x_vars(i)-lwr)/(upr - lwr);
	// GSL beta passes alpha followed by beta
	Real cdf = gsl_cdf_beta_P(scaled_x, alpha, beta);
	z_vars(i) = Phi_inverse(cdf);
      }
#endif // PECOS_GSL
      break;
    }
    case GAMMA: {
      const Real& beta = ranVarAddtlParamsX[i](1);
      if (extendedUSpace) // scale to std gamma
	z_vars(i) = x_vars(i)/beta;
#ifdef PECOS_GSL
      else { // transform to std normal
	const Real& alpha = ranVarAddtlParamsX[i](0);
	// GSL gamma passes alpha followed by beta
	Real cdf  = gsl_cdf_gamma_P(x_vars(i), alpha, beta);
	z_vars(i) = Phi_inverse(cdf);
      }
#endif // PECOS_GSL
      break;
    }
    case GUMBEL: {
      const Real& alpha = ranVarAddtlParamsX[i](0);
      const Real& beta  = ranVarAddtlParamsX[i](1);
      Real cdf  = exp(-exp(-alpha*(x_vars(i) - beta)));
      z_vars(i) = Phi_inverse(cdf);
      break;
    }
    case FRECHET: {
      const Real& alpha = ranVarAddtlParamsX[i](0);
      const Real& beta  = ranVarAddtlParamsX[i](1);
      Real cdf  = exp(-pow(beta/x_vars(i), alpha));
      z_vars(i) = Phi_inverse(cdf);
      break;
    }
#ifdef PECOS_GSL
    case WEIBULL: {
      const Real& alpha = ranVarAddtlParamsX[i](0);
      const Real& beta  = ranVarAddtlParamsX[i](1);
      // GSL weibull passes beta followed by alpha
      Real cdf  = gsl_cdf_weibull_P(x_vars(i), beta, alpha);
      z_vars(i) = Phi_inverse(cdf);
      break;
    }
#endif // PECOS_GSL
    }
  }
}


/** This procedure computes the transformation from z to u space.
    u_vars is the vector of random variables in uncorrelated standard
    normal space (u-space).  z_vars is the vector of random variables
    in normal space with proper correlations (z-space). */
void NatafTransformation::trans_Z_to_U(RealVector& z_vars, RealVector& u_vars)
{
  // corrCholeskyFactorZ: the Cholesky factor of the modified correlation matrix

  int z_len = z_vars.Length();
  if (u_vars.Length() != z_len)
    u_vars.Size(z_len);

  Epetra_SerialDenseSolver corr_solver;
  corr_solver.SetMatrix(corrCholeskyFactorZ);
  corr_solver.SetVectors(u_vars, z_vars);
  corr_solver.SolveToRefinedSolution(true);
  corr_solver.Solve();
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
void NatafTransformation::trans_correlations()
{
  // ranVarTypesX:   the type of random variable.  Supported correlations are
  //                 NORMAL/LOGNORMAL/UNIFORM/GAMMA/GUMBEL/FRECHET/WEIBULL.
  // ranVarMeansX:   means of the random variables in x-space
  // ranVarStdDevsX: standard deviations of the random variables in x-space
  // cf_var_i:       coefficient of variation of Xi
  // cf_var_j:       coefficient of variation of Xj

  RealSymMatrix mod_corr_matrix(corrMatrixX); // copy

  size_t i, j;
  for (i=numDesignVars; i<numDesignVars+numUncertainVars; i++) {
    for (j=numDesignVars; j<i; j++) {

      if (ranVarTypesU[i] == NORMAL && ranVarTypesU[j] == NORMAL) {

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
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
	  mod_corr_matrix(i,j) *= cf_var/sqrt(log(1. + cf_var*cf_var));
	}
	//        Normal <--> Gamma          (Max Error = 0.0%)
	//
	else if ( (ranVarTypesX[i] == GAMMA  && ranVarTypesX[j] == NORMAL) ||
		  (ranVarTypesX[i] == NORMAL && ranVarTypesX[j] == GAMMA) ) {
	  Real cf_var = (ranVarTypesX[i] == GAMMA) ?
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
	  mod_corr_matrix(i,j) *= 1.001 - 0.007*cf_var + 0.118*cf_var*cf_var;
	}
	//        Normal <--> Frechet        (Max Error = 0.1%)
	//
	else if ( (ranVarTypesX[i] == FRECHET && ranVarTypesX[j] == NORMAL) ||
		  (ranVarTypesX[i] == NORMAL  && ranVarTypesX[j] == FRECHET) ) {
	  Real cf_var = (ranVarTypesX[i] == FRECHET) ?
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
	  mod_corr_matrix(i,j) *= 1.03 + 0.238*cf_var + 0.364*cf_var*cf_var;
	}
	//        Normal <--> Weibull        (Max Error = 0.1%)
	//
	else if ( (ranVarTypesX[i] == WEIBULL && ranVarTypesX[j] == NORMAL) ||
		  (ranVarTypesX[i] == NORMAL  && ranVarTypesX[j] == WEIBULL) ) {
	  Real cf_var = (ranVarTypesX[i] == WEIBULL) ?
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
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
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
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
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
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
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
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
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
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
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
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
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
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
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
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
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
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
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
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
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
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
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
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
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
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
	  Real cf_var_i = ranVarStdDevsX(i)/ranVarMeansX(i),
	       cf_var_j = ranVarStdDevsX(j)/ranVarMeansX(j);
	  // Note: "=" instead of "*=" eliminates redundant rho multiply/divide 
	  mod_corr_matrix(i,j) = log(1.+cf_var_i*cf_var_j*corrMatrixX(i,j)) /
	    sqrt(log(1. + cf_var_i*cf_var_i)*log(1. + cf_var_j*cf_var_j));
	}
	//        Gamma <--> Lognormal       (Max Error = 4.0%)
	//
	else if ( (ranVarTypesX[i] == GAMMA && ranVarTypesX[j] == LOGNORMAL) ||
		  (ranVarTypesX[i] == LOGNORMAL && ranVarTypesX[j] == GAMMA) ) {
	  Real cf_var_g = (ranVarTypesX[i] == GAMMA) ?
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
	  Real cf_var_l = (ranVarTypesX[i] == LOGNORMAL) ?
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
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
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
	  Real cf_var_f = (ranVarTypesX[i] == FRECHET) ?
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
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
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
	  Real cf_var_l = (ranVarTypesX[i] == LOGNORMAL) ?
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
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
	  Real cf_var_i = ranVarStdDevsX(i)/ranVarMeansX(i),
	       cf_var_j = ranVarStdDevsX(j)/ranVarMeansX(j);
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
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
	  Real cf_var_f = (ranVarTypesX[i] == FRECHET) ?
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
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
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
	  Real cf_var_w = (ranVarTypesX[i] == WEIBULL) ?
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
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
	  Real corr2 = corr*corr, cf_var_i  = ranVarStdDevsX(i)/ranVarMeansX(i),
	       cf_var_i2 = cf_var_i*cf_var_i,
	       cf_var_j  = ranVarStdDevsX(j)/ranVarMeansX(j),
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
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
	  Real cf_var_f = (ranVarTypesX[i] == FRECHET) ?
	    ranVarStdDevsX(i)/ranVarMeansX(i) :
	    ranVarStdDevsX(j)/ranVarMeansX(j);
	  mod_corr_matrix(i,j) *= 1.065 + 0.146*corr + 0.241*cf_var_f
                               -  0.259*cf_var_w + 0.013*corr*corr
	                       +  0.372*cf_var_f*cf_var_f
	                       +  0.435*cf_var_w*cf_var_w + 0.005*corr*cf_var_f
	                       +  0.034*cf_var_f*cf_var_w - 0.481*corr*cf_var_w;
	}
	//        Weibull <--> Weibull       (Max Error = 2.6%)
	//
	else if (ranVarTypesX[i] == WEIBULL && ranVarTypesX[j] == WEIBULL) {
	  Real cf_var_i = ranVarStdDevsX(i)/ranVarMeansX(i),
	       cf_var_j = ranVarStdDevsX(j)/ranVarMeansX(j);
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
  Epetra_SerialSpdDenseSolver corr_solver;
  corr_solver.SetMatrix(mod_corr_matrix);
  corr_solver.Factor(); // Cholesky factorization (LL^T) in place
  // Define corrCholeskyFactorZ to be L by assigning the lower triangle.
  size_t num_active_vars = numDesignVars + numUncertainVars + numStateVars;
  if (corrCholeskyFactorZ.M() != num_active_vars ||
      corrCholeskyFactorZ.N() != num_active_vars)
    corrCholeskyFactorZ.Shape(num_active_vars, num_active_vars);
  for (i=0; i<num_active_vars; i++)
    for (j=0; j<=i; j++)
      corrCholeskyFactorZ(i, j) = mod_corr_matrix(i, j);
  //Cout << "\ncorrCholeskyFactorZ:" << corrCholeskyFactorZ;

  // could pre-compute L^-1 to avoid solving L u = z repeatedly for u
  //corrCholeskyFactorZInv.Shape(num_active_vars, num_active_vars);
  //corrCholeskyFactorZInv = corrCholeskyFactorZ; // copy
  //Epetra_SerialDenseSolver chol_solver;
  //chol_solver.SetMatrix(corrCholeskyFactorZInv); 
  //chol_solver.Invert();
  //Cout << "\ncorrCholeskyFactorZInv:" << corrCholeskyFactorZInv;
}


/** This procedure tranforms a gradient vector dg/dx from the original
    user-defined x-space (where evaluations are performed) to uncorrelated
    standard normal space (u-space) through application of the Jacobian dx/du.
    x_vars is the vector of random variables in x-space. */
void NatafTransformation::
trans_grad_X_to_U(const RealVector& fn_grad_x, RealVector& fn_grad_u,
		  const RealVector& x_vars,    const IntArray& x_dvv)
{
  RealMatrix jacobian_xu;
  jacobian_dX_dU(x_vars, jacobian_xu);
  trans_grad_X_to_U(fn_grad_x, fn_grad_u, jacobian_xu, x_dvv);
}


/** This procedure tranforms a gradient vector dg/dx from the original
    user-defined x-space (where evaluations are performed) to uncorrelated
    standard normal space (u-space) through application of the Jacobian dx/du.
    This overloaded form allows for the separate calculation of jacobian_xu,
    as this matrix is independent of the response function index and can be
    pulled outside response function loops. */
void NatafTransformation::
trans_grad_X_to_U(const RealVector& fn_grad_x,   RealVector& fn_grad_u,
		  const RealMatrix& jacobian_xu, const IntArray& x_dvv)
{
  // Jacobian dimensions = length of ranVarTypesX/U = iteratedModel.cv()
  int x_len = jacobian_xu.M();
  const IntArray& cv_ids = iteratedModel.continuous_variable_ids();
  if (x_dvv == cv_ids) { // standard DVV
    if (fn_grad_x.Length() != x_len) {
      Cerr << "Error: bad fn_grad_x dimension in NatafTransformation::"
	   << "trans_grad_X_to_U()." << endl;
      abort_handler(-1);
    }
    if (fn_grad_u.Length() != x_len)
      fn_grad_u.Size(x_len);
    fn_grad_u.Multiply('T', 'N', 1., jacobian_xu, fn_grad_x, 0.);
  }
  else { // non-standard DVV
    RealVector fn_grad_x_trans(x_len), fn_grad_u_trans(x_len);
    size_t i, dvv_index, num_deriv_vars = x_dvv.length();
    SizetArray dvv_index_array(x_len);
    // extract relevant DVV components from fn_grad_x
    for (i=0; i<x_len; i++) {
      dvv_index_array[i] = dvv_index = x_dvv.index(cv_ids[i]);
      if (dvv_index != _NPOS)
	fn_grad_x_trans(i) = fn_grad_x(dvv_index);
    }
    // perform transformation using full Jacobian
    fn_grad_u_trans.Multiply('T', 'N', 1., jacobian_xu, fn_grad_x_trans, 0.);
    // copy relevant DVV components into fn_grad_u
    if (fn_grad_u.Length() != num_deriv_vars)
      fn_grad_u.Size(num_deriv_vars);
    for (i=0; i<x_len; i++) {
      dvv_index = dvv_index_array[i];
      if (dvv_index != _NPOS)
	fn_grad_u(dvv_index) = fn_grad_u_trans(i);
    }
  }
}


/** This procedure tranforms a gradient vector dg/du from uncorrelated standard
    space (u-space) to the original user-defined x-space through application of
    the Jacobian du/dx.  x_vars is the vector of random variables in x-space. */
void NatafTransformation::
trans_grad_U_to_X(const RealVector& fn_grad_u, RealVector& fn_grad_x,
		  const RealVector& x_vars,    const IntArray& x_dvv)
{
  RealMatrix jacobian_ux;
  jacobian_dU_dX(x_vars, jacobian_ux);
  trans_grad_U_to_X(fn_grad_u, fn_grad_x, jacobian_ux, x_dvv);
}


/** This procedure tranforms a gradient vector dg/du from uncorrelated standard
    space (u-space) to the original user-defined x-space through application of
    the Jacobian du/dx.  This overloaded form allows for the separate
    calculation of jacobian_ux, as this matrix is independent of the response
    function index and can be pulled outside response function loops. */
void NatafTransformation::
trans_grad_U_to_X(const RealVector& fn_grad_u,   RealVector& fn_grad_x,
		  const RealMatrix& jacobian_ux, const IntArray& x_dvv)
{
  // Jacobian dimensions = length of ranVarTypesX/U = iteratedModel.cv()
  int u_len = jacobian_ux.M();
  const IntArray& cv_ids = iteratedModel.continuous_variable_ids();
  if (x_dvv == cv_ids) { // standard DVV
    if (fn_grad_u.Length() != u_len) {
      Cerr << "Error: bad fn_grad_u dimension in NatafTransformation::"
	   << "trans_grad_U_to_X()." << endl;
      abort_handler(-1);
    }
    if (fn_grad_x.Length() != u_len)
      fn_grad_x.Size(u_len);
    fn_grad_x.Multiply('T', 'N', 1., jacobian_ux, fn_grad_u, 0.);
  }
  else { // non-standard DVV
    RealVector fn_grad_u_trans(u_len), fn_grad_x_trans(u_len);
    size_t dvv_index, num_deriv_vars = x_dvv.length();
    SizetArray dvv_index_array(u_len);
    // extract relevant DVV components from fn_grad_u
    for (int i=0; i<u_len; i++) {
      dvv_index_array[i] = dvv_index = x_dvv.index(cv_ids[i]);
      if (dvv_index != _NPOS)
	fn_grad_u_trans(i) = fn_grad_u(dvv_index);
    }
    // perform transformation using full Jacobian
    fn_grad_x_trans.Multiply('T', 'N', 1., jacobian_ux, fn_grad_u_trans, 0.);
    // copy relevant DVV components into fn_grad_x
    if (fn_grad_x.Length() != num_deriv_vars)
      fn_grad_x.Size(num_deriv_vars);
    for (int i=0; i<u_len; i++) {
      dvv_index = dvv_index_array[i];
      if (dvv_index != _NPOS)
	fn_grad_x(dvv_index) = fn_grad_x_trans(i);
    }
  }
}


/** This procedure multiplies a gradient vector dg/dx from the
    original user-defined x-space (where evaluations are performed)
    with the design Jacobian dx/ds of the transformation x = x(u,s) to
    form the design gradient dg/ds.  x_vars is the vector of random
    variables in x-space. */
void NatafTransformation::
trans_grad_X_to_S(const RealVector& fn_grad_x,
		  RealVector& fn_grad_s,
		  const RealVector& x_vars, const IntArray& x_dvv)
{
  RealMatrix jacobian_xs;
  jacobian_dX_dS(x_vars, jacobian_xs);
  trans_grad_X_to_S(fn_grad_x, fn_grad_s, jacobian_xs, x_dvv);
}


/** This procedure multiplies a gradient vector dg/dx from the
    original user-defined x-space (where evaluations are performed)
    with the design Jacobian dx/ds of the transformation x = x(u,s) to
    form the design gradient dg/ds.  This overloaded form allows
    for the separate calculation of jacobian_xs, as this matrix is
    independent of the response function index and can be pulled
    outside response function loops. */
void NatafTransformation::
trans_grad_X_to_S(const RealVector& fn_grad_x,   RealVector& fn_grad_s,
		  const RealMatrix& jacobian_xs, const IntArray& x_dvv)
{
  // Jacobian dim is x_len by s_len, where
  // > x_len = length of ranVarTypesX/U        = inner model.cv()
  // > s_len = primaryACVarMapIndices.length() = outer model.cv()
  int x_len = jacobian_xs.M(), s_len = jacobian_xs.N();
  if (primaryACVarMapIndices.empty() || secondaryACVarMapTargets.empty()) {
    Cerr << "Error: NatafTransformation::trans_grad_X_to_S() requires primary "
	 << "and secondary variable mappings to define S." << endl;
    abort_handler(-1);
  }

  const IntArray&  cv_ids = iteratedModel.continuous_variable_ids();
  const IntArray& acv_ids = iteratedModel.all_continuous_variable_ids();
  bool std_dvv = (x_dvv == cv_ids);
  bool mixed_s = secondaryACVarMapTargets.contains(NO_TARGET);
  RealVector fn_grad_x_std, fn_grad_s_std;

  // manage size of fn_grad_x input
  if (std_dvv) { // standard x-space DVV
    if (fn_grad_x.Length() != x_len) {
      Cerr << "Error: bad fn_grad_x dimension in NatafTransformation::"
	   << "trans_grad_X_to_S()." << endl;
      abort_handler(-1);
    }
  }
  else { // non-standard x-space DVV
    fn_grad_x_std.Size(x_len);
    // extract relevant DVV components from fn_grad_x
    size_t i, dvv_index;
    for (i=0; i<x_len; i++) {
      dvv_index = x_dvv.index(cv_ids[i]);
      if (dvv_index != _NPOS)
	fn_grad_x_std(i) = fn_grad_x(dvv_index);
    }
  }

  // manage size of fn_grad_s output
  if (mixed_s || !std_dvv)
    fn_grad_s_std.Size(s_len);
  size_t i, final_s_len;
  if (std_dvv)
    final_s_len = s_len;
  else {
    final_s_len = 0;
    for (i=0; i<s_len; i++)
      if (x_dvv.contains(acv_ids[primaryACVarMapIndices[i]]))
	final_s_len++;
  }
  if (fn_grad_s.Length() != final_s_len)
    fn_grad_s.Size(final_s_len);

  // perform transformation using full Jacobian
  const RealVector& fn_grad_x_trans
    = (std_dvv) ? fn_grad_x : fn_grad_x_std;
  RealVector& fn_grad_s_trans
    = (mixed_s || !std_dvv) ? fn_grad_s_std : fn_grad_s;
  fn_grad_s_trans.Multiply('T', 'N', 1., jacobian_xs, fn_grad_x_trans, 0.);

  // reassemble final fn_grad_s
  if (mixed_s || !std_dvv) {
    size_t cntr = 0;
    for (i=0; i<s_len; i++) {
      int acv_id = acv_ids[primaryACVarMapIndices[i]];
      size_t dvv_index = x_dvv.index(acv_id);
      if (dvv_index != _NPOS)
	fn_grad_s(cntr++) = (secondaryACVarMapTargets[i] == NO_TARGET) ?
	  fn_grad_x(dvv_index) : // no distribution parameter: if the missing
	  // fn_grad_s component is available in fn_grad_x, then use it; else it
	  // must be updated separately (as in NonDLocalReliability::dg_ds_eval)
	  fn_grad_s_trans(i);    // use the Jacobian-transformed component
      else if (std_dvv)
	fn_grad_s(cntr++) = 0.;
    }
  }

#ifdef DEBUG
  Cout << "\nfn_grad_x:\n"       << fn_grad_x
       << "\nfn_grad_x_trans:\n" << fn_grad_x_trans
       << "\njacobian_xs:\n"     << jacobian_xs
       << "\nfn_grad_s_trans:\n" << fn_grad_s_trans
       << "\nfn_grad_s:\n"       << fn_grad_s << endl;
#endif
}


/** This procedure tranforms a Hessian matrix from the original
    user-defined x-space (where evaluations are performed) to
    uncorrelated standard normal space (u-space).  x_vars is the
    vector of the random variables in x-space. */
void NatafTransformation::
trans_hess_X_to_U(const RealSymMatrix& fn_hess_x, RealSymMatrix& fn_hess_u,
		  const RealVector& x_vars,	  const RealVector& fn_grad_x,
		  const IntArray& x_dvv)
{
  RealMatrix jacobian_xu;
  jacobian_dX_dU(x_vars, jacobian_xu);

  RealSymMatrixArray hessian_xu;
  if ( (!extendedUSpace && numUncertainVars != ranVarTypesX.count(NORMAL)) ||
       ( extendedUSpace && numUncertainVars != ranVarTypesX.count(NORMAL) +
	 ranVarTypesX.count(UNIFORM) + ranVarTypesX.count(EXPONENTIAL) +
	 ranVarTypesX.count(BETA)    + ranVarTypesX.count(GAMMA) ) )
    hessian_d2X_dU2(x_vars, hessian_xu); // nonlinear transformation has Hessian

  trans_hess_X_to_U(fn_hess_x, fn_hess_u, jacobian_xu, hessian_xu,
		    fn_grad_x, x_dvv);
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
		  const RealVector& fn_grad_x,    const IntArray& x_dvv)
{
  // Jacobian dimensions = length of ranVarTypesX/U = iteratedModel.cv()
  int x_len = jacobian_xu.M();
  const IntArray& cv_ids = iteratedModel.continuous_variable_ids();
  bool std_dvv = (x_dvv == cv_ids); // standard DVV
  bool nonlinear_vars_map = !hessian_xu.empty();

  RealSymMatrix fn_hess_x_std, fn_hess_u_std;
  RealVector fn_grad_x_std;
  SizetArray dvv_index_array;
  if (std_dvv) {
    if (fn_hess_x.M() != x_len) {
      Cerr << "Error: bad fn_hess_x dimension in NatafTransformation::"
	   << "trans_hess_X_to_U()." << endl;
      abort_handler(-1);
    }
    if ( nonlinear_vars_map &&
	 ( fn_grad_x.Length() != x_len || hessian_xu.length() != x_len ) ) {
      Cerr << "Error: bad dimension in NatafTransformation::"
	   << "trans_hess_X_to_U()." << endl;
      abort_handler(-1);
    }
    if (fn_hess_u.M() != x_len)
      fn_hess_u.Shape(x_len);
  }
  else { // extract relevant DVV components from fn_grad_x & fn_hess_x
    fn_hess_x_std.Shape(x_len);
    fn_hess_u_std.Shape(x_len);
    if (nonlinear_vars_map)
      fn_grad_x_std.Size(x_len);
    size_t i, j, dvv_index_i, dvv_index_j, num_deriv_vars = x_dvv.length();
    dvv_index_array.reshape(x_len);
    for (i=0; i<x_len; i++)
      dvv_index_array[i] = dvv_index_i = x_dvv.index(cv_ids[i]);
    if (fn_hess_u.M() != num_deriv_vars)
      fn_hess_u.Shape(num_deriv_vars);
    // extract relevant DVV components from fn_hess_x
    for (i=0; i<x_len; i++) {
      dvv_index_i = dvv_index_array[i];
      if (dvv_index_i != _NPOS) {
	if (nonlinear_vars_map)
	  fn_grad_x_std(i) = fn_grad_x(dvv_index_i);
	for (j=0; j<x_len; j++) {
	  dvv_index_j = dvv_index_array[j];
	  if (dvv_index_j != _NPOS)
	    fn_hess_x_std(i, j) = fn_hess_x(dvv_index_i, dvv_index_j);
	}
      }
    }
  }
  const RealVector& fn_grad_x_trans
    = (std_dvv) ? fn_grad_x : fn_grad_x_std;
  const RealSymMatrix& fn_hess_x_trans
    = (std_dvv) ? fn_hess_x : fn_hess_x_std;
  RealSymMatrix& fn_hess_u_trans
    = (std_dvv) ? fn_hess_u : fn_hess_u_std;

  // transform hess_x -> hess_u
  // d^2G/dU^2 = dG/dX^T d^2X/dU^2 + dX/dU^T d^2G/dX^2 dX/dU
  // Note: G(u) may have curvature even if g(x) is linear due to first term.
  RealMatrix fn_hess_x_by_xu(x_len, x_len);
  fn_hess_x_by_xu.Multiply('L',       1., fn_hess_x_trans, jacobian_xu, 0.);
  fn_hess_u_trans.Multiply('T', 'N',  1., jacobian_xu, fn_hess_x_by_xu, 0.);
  //Cout << "\nfnHessU 1st term J^T H_x J:" << fn_hess_u;
  //Cout << "\nfnGradX:" << fn_grad_x;

  if (nonlinear_vars_map) { // nonlinear transformation has Hessian
    for (int i=0; i<x_len; i++) {
      //Cout << "\nhessian_xu[" << i << "]:" << hessian_xu[i];
      //hessian_xu[i].Scale(fn_grad_x(i));
      //fn_hess_u += hessian_xu[i];
      const double&                      fn_grad_x_i  = fn_grad_x_trans(i);
      const RealSymMatrix& hessian_xu_i = hessian_xu[i];
      for (int j=0; j<x_len; j++)
	for (int k=0; k<x_len; k++)
	  fn_hess_u_trans(j,k) += fn_grad_x_i * hessian_xu_i(j,k);
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
    int x_len = x_vars.Length();
    if (jacobian_xu.M() != x_len || jacobian_xu.N() != x_len)
      jacobian_xu.Shape(x_len, x_len);
    jacobian_xu.Multiply('N', 'N', 1., jacobian_xz, corrCholeskyFactorZ, 0.);
  }
  else // dX/dU = dX/dZ since dZ/dU = I
    jacobian_dX_dZ(x_vars, jacobian_xu);
}


/** This procedure computes the Jacobian of the transformation x(z).
    x_vars is the vector of random variables in the original
    user-defined x-space. */
void NatafTransformation::
jacobian_dX_dZ(const RealVector& x_vars, RealMatrix& jacobian_xz) 
{
  int x_len = x_vars.Length();
  if (jacobian_xz.M() != x_len || jacobian_xz.N() != x_len)
    jacobian_xz.Shape(x_len, x_len);

  // Rackwitz-Fiessler: Phi(z) = F(x)
  // d/dz -> phi(z) = f(x) dx/dz
  //         dx/dz = phi(z)/f(x)
  // dX/dZ is diagonal as defined by differentiation of trans_Z_to_X()

  RealVector z_vars;
  if ( ranVarTypesX.contains(TRIANGULAR) || ranVarTypesX.contains(GUMBEL)    ||
       ranVarTypesX.contains(FRECHET)    || ranVarTypesX.contains(WEIBULL)   ||
       ( !extendedUSpace && ( ranVarTypesX.contains(UNIFORM) ||
	  ranVarTypesX.contains(LOGUNIFORM)  ||
	  ranVarTypesX.contains(EXPONENTIAL) || ranVarTypesX.contains(BETA)  ||
	  ranVarTypesX.contains(GAMMA) ) ) )
    trans_X_to_Z(x_vars, z_vars);

  for (int i=0; i<x_len; i++) {
    switch (ranVarTypesX[i]) {
    case DESIGN: case STATE:
      jacobian_xz(i, i) = (ranVarUpperBndsX(i) - ranVarLowerBndsX(i))/2.;
      break;
    case NORMAL: {
      // z = (x - mu)/sigma
      jacobian_xz(i, i) = ranVarStdDevsX(i);
      break;
    }
    case LOGNORMAL: {
      // z = (ln x - lamba)/zeta
      Real cf_var = ranVarStdDevsX(i)/ranVarMeansX(i),
	zeta = sqrt(log(1 + cf_var*cf_var));
      jacobian_xz(i, i) = zeta*x_vars(i);
      break;
    }
    case UNIFORM: {
      // F(x) = (x-L)/(U-L), f(x) = 1/(U-L)
      Real scale = ranVarUpperBndsX(i) - ranVarLowerBndsX(i);
      if (extendedUSpace) // linear scaling
	jacobian_xz(i, i) = scale/2.;
      else { // nonlinear transformation
	Real pdf = 1./scale;
	jacobian_xz(i, i) = phi(z_vars(i))/pdf;
      }
      break;
    }
    case LOGUNIFORM: {
      // F(x) = (ln x - ln L)/(ln U - ln L), f(x) = 1/x/(ln U - ln L)
      Real log_range = log(ranVarUpperBndsX(i)) - log(ranVarLowerBndsX(i));
      if (extendedUSpace)
	jacobian_xz(i, i) = x_vars(i)*log_range/2.;
      else {
	Real pdf = 1./x_vars(i)/log_range;
	jacobian_xz(i, i) = phi(z_vars(i))/pdf;
      }
      break;
    }
    case TRIANGULAR: {
      const Real& lwr  = ranVarLowerBndsX(i);
      const Real& mode = ranVarAddtlParamsX[i](0);
      const Real& upr  = ranVarUpperBndsX(i);
      const Real& x    = x_vars(i);
      Real range       = upr - lwr;
      Real pdf = (x < mode) ? 2.*(x-lwr)/range/(mode-lwr)
                            : 2.*(upr-x)/range/(upr-mode);
      jacobian_xz(i, i) = phi(z_vars(i))/pdf;
      break;
    }
    case EXPONENTIAL: {
      const Real& beta = ranVarAddtlParamsX[i](0);
      if (extendedUSpace) // linear scaling
	jacobian_xz(i, i) = beta;
      else { // nonlinear transformation
	Real pdf = exp(-x_vars(i)/beta)/beta;
	jacobian_xz(i, i) = phi(z_vars(i))/pdf;
      }
      break;
    }
    case BETA: {
      const Real& lwr = ranVarLowerBndsX(i);
      const Real& upr = ranVarUpperBndsX(i);
      Real scale = upr - lwr;
      if (extendedUSpace) // linear scaling
	jacobian_xz(i, i) = scale/2.;
#ifdef PECOS_GSL
      else { // nonlinear transformation
	const Real& alpha = ranVarAddtlParamsX[i](0);
	const Real& beta  = ranVarAddtlParamsX[i](1);
	Real scaled_x = (x_vars(i)-lwr)/scale;
	// GSL beta passes alpha followed by beta
	Real pdf = gsl_ran_beta_pdf(scaled_x, alpha, beta)/scale;
	jacobian_xz(i, i) = phi(z_vars(i))/pdf;
      }
#endif // PECOS_GSL
      break;
    }
    case GAMMA: {
      const Real& beta = ranVarAddtlParamsX[i](1);
      if (extendedUSpace) // linear scaling
	jacobian_xz(i, i) = beta;
#ifdef PECOS_GSL
      else { // nonlinear transformation
	const Real& alpha = ranVarAddtlParamsX[i](0);
	// GSL gamma passes alpha followed by beta
	Real pdf = gsl_ran_gamma_pdf(x_vars(i), alpha, beta);
	jacobian_xz(i, i) = phi(z_vars(i))/pdf;
      }
#endif // PECOS_GSL
      break;
    }
    case GUMBEL: {
      const Real& alpha = ranVarAddtlParamsX[i](0);
      const Real& beta  = ranVarAddtlParamsX[i](1);
      Real num = exp(-alpha*(x_vars(i)-beta)), pdf = alpha*num*exp(-num);
      jacobian_xz(i, i) = phi(z_vars(i))/pdf;
      break;
    }
    case FRECHET: {
      const Real& alpha = ranVarAddtlParamsX[i](0);
      const Real& beta  = ranVarAddtlParamsX[i](1);
      Real num = beta/x_vars(i),
	   pdf = alpha/beta*pow(num,alpha+1.)*exp(-pow(num,alpha));
      jacobian_xz(i, i) = phi(z_vars(i))/pdf;
      break;
    }
#ifdef PECOS_GSL
    case WEIBULL: {
      const Real& alpha = ranVarAddtlParamsX[i](0);
      const Real& beta  = ranVarAddtlParamsX[i](1);
      // GSL weibull passes beta followed by alpha
      Real pdf = gsl_ran_weibull_pdf(x_vars(i), beta, alpha);
      jacobian_xz(i, i) = phi(z_vars(i))/pdf;
      break;
    }
#endif // PECOS_GSL
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
    Epetra_SerialDenseSolver corr_solver;
    corr_solver.SetMatrix(corrCholeskyFactorZ);
    int x_len = x_vars.Length();
    if (jacobian_ux.M() != x_len || jacobian_ux.N() != x_len)
      jacobian_ux.Shape(x_len, x_len);
    corr_solver.SetVectors(jacobian_ux, jacobian_zx);
    corr_solver.SolveToRefinedSolution(true);
    corr_solver.Solve();
  }
  else // dU/dX = dZ/dX since dU/dZ = I
    jacobian_dZ_dX(x_vars, jacobian_ux);
}


/** This procedure computes the Jacobian of the transformation z(x).
    x_vars is the vector of random variables in the original
    user-defined x-space. */
void NatafTransformation::
jacobian_dZ_dX(const RealVector& x_vars, RealMatrix& jacobian_zx) 
{
  int x_len = x_vars.Length();
  if (jacobian_zx.M() != x_len || jacobian_zx.N() != x_len)
    jacobian_zx.Shape(x_len, x_len);

  // Rackwitz-Fiessler: Phi(z) = F(x)
  // d/dx -> phi(z) dz/dx = f(x)
  //         dz/dx = f(x)/phi(z)
  // dZ/dX is diagonal as defined by differentiation of trans_X_to_Z()

  RealVector z_vars;
  if ( ranVarTypesX.contains(TRIANGULAR) || ranVarTypesX.contains(GUMBEL)    ||
       ranVarTypesX.contains(FRECHET)    || ranVarTypesX.contains(WEIBULL)   ||
       ( !extendedUSpace && ( ranVarTypesX.contains(UNIFORM) ||
          ranVarTypesX.contains(LOGUNIFORM)  ||
	  ranVarTypesX.contains(EXPONENTIAL) || ranVarTypesX.contains(BETA)  ||
	  ranVarTypesX.contains(GAMMA) ) ) )
    trans_X_to_Z(x_vars, z_vars);

  for (int i=0; i<x_len; i++) {
    switch (ranVarTypesX[i]) {
    case DESIGN: case STATE:
      jacobian_zx(i, i) = 2./(ranVarUpperBndsX(i) - ranVarLowerBndsX(i));
      break;
    case NORMAL:
      // z = (x - mu)/sigma
      jacobian_zx(i, i) = 1./ranVarStdDevsX(i);
      break;
    case LOGNORMAL: {
      // z = (ln x - lamba)/zeta
      Real cf_var = ranVarStdDevsX(i)/ranVarMeansX(i),
	zeta = sqrt(log(1 + cf_var*cf_var));
      jacobian_zx(i, i) = 1./zeta/x_vars(i);
      break;
    }
    case UNIFORM: {
      // F(x) = (x-L)/(U-L), f(x) = 1/(U-L)
      Real scale = ranVarUpperBndsX(i) - ranVarLowerBndsX(i);
      if (extendedUSpace) // linear scaling
	jacobian_zx(i, i) = 2./scale;
      else {
	Real pdf = 1./scale;
	jacobian_zx(i, i) = pdf/phi(z_vars(i));
      }
      break;
    }
    case LOGUNIFORM: {
      // F(x) = (ln x - ln L)/(ln U - ln L), f(x) = 1/x/(ln U - ln L)
      Real log_range = log(ranVarUpperBndsX(i)) - log(ranVarLowerBndsX(i));
      if (extendedUSpace)
	jacobian_zx(i, i) = 2./x_vars(i)/log_range;
      else {
	Real pdf = 1./x_vars(i)/log_range;
	jacobian_zx(i, i) = pdf/phi(z_vars(i));
      }
      break;
    }
    case TRIANGULAR: {
      const Real& lwr  = ranVarLowerBndsX(i);
      const Real& mode = ranVarAddtlParamsX[i](0);
      const Real& upr  = ranVarUpperBndsX(i);
      const Real& x    = x_vars(i);
      Real range       = upr - lwr;
      Real pdf = (x < mode) ? 2.*(x-lwr)/range/(mode-lwr)
	                    : 2.*(upr-x)/range/(upr-mode);
      jacobian_zx(i, i) = pdf/phi(z_vars(i));
      break;
    }
    case EXPONENTIAL: {
      const Real& beta = ranVarAddtlParamsX[i](0);
      if (extendedUSpace) // linear scaling
	jacobian_zx(i, i) = 1./beta;
      else {
	Real pdf = exp(-x_vars(i)/beta)/beta;
	jacobian_zx(i, i) = pdf/phi(z_vars(i));
      }
      break;
    }
#ifdef PECOS_GSL
    case BETA: {
      const Real& lwr = ranVarLowerBndsX(i);
      const Real& upr = ranVarUpperBndsX(i);
      Real scale = upr - lwr;
      if (extendedUSpace) // linear scaling
	jacobian_zx(i, i) = 2./scale;
      else {
	const Real& alpha = ranVarAddtlParamsX[i](0);
	const Real& beta  = ranVarAddtlParamsX[i](1);
	Real scaled_x = (x_vars(i)-lwr)/scale;
	// GSL beta passes alpha followed by beta
	Real pdf = gsl_ran_beta_pdf(scaled_x, alpha, beta)/scale;
	jacobian_zx(i, i) = pdf/phi(z_vars(i));
      }
      break;
    }
    case GAMMA: {
      const Real& beta = ranVarAddtlParamsX[i](1);
      if (extendedUSpace) // linear scaling
	jacobian_zx(i, i) = 1./beta;
      else {
	const Real& alpha = ranVarAddtlParamsX[i](0);
	// GSL gamma passes alpha followed by beta
	Real pdf = gsl_ran_gamma_pdf(x_vars(i), alpha, beta);
	jacobian_zx(i, i) = pdf/phi(z_vars(i));
      }
      break;
    }
#endif // PECOS_GSL
    case GUMBEL: {
      const Real& alpha = ranVarAddtlParamsX[i](0);
      const Real& beta  = ranVarAddtlParamsX[i](1);
      Real num = exp(-alpha*(x_vars(i)-beta)), pdf = alpha*num*exp(-num);
      jacobian_zx(i, i) = pdf/phi(z_vars(i));
      break;
    }
    case FRECHET: {
      const Real& alpha = ranVarAddtlParamsX[i](0);
      const Real& beta  = ranVarAddtlParamsX[i](1);
      Real num = beta/x_vars(i),
	   pdf = alpha/beta*pow(num,alpha+1.)*exp(-pow(num,alpha));
      jacobian_zx(i, i) = pdf/phi(z_vars(i));
      break;
    }
#ifdef PECOS_GSL
    case WEIBULL: {
      const Real& alpha = ranVarAddtlParamsX[i](0);
      const Real& beta  = ranVarAddtlParamsX[i](1);
      // GSL weibull passes beta followed by alpha
      Real pdf = gsl_ran_weibull_pdf(x_vars(i), beta, alpha);
      jacobian_zx(i, i) = pdf/phi(z_vars(i));
      break;
    }
#endif // PECOS_GSL
    }
  }
}


/** This procedure computes the derivative of the original variables x
    with respect to the random variable distribution parameters s.
    This provides the design Jacobian of the transformation for use in
    computing statistical design sensitivities for OUU. */
void NatafTransformation::
jacobian_dX_dS(const RealVector& x_vars, RealMatrix& jacobian_xs)
{
  // Rectangular Jacobian = Gradient^T = num_X by num_S where num_S is the total
  // number of active continuous vars flowed down from a higher iteration level.
  // The number of distribution parameter insertions is <= num_S.
  size_t i, j, num_var_map_1c = primaryACVarMapIndices.length();
  int x_len = x_vars.Length();
  if (jacobian_xs.M() != x_len || jacobian_xs.N() != num_var_map_1c)
    jacobian_xs.Shape(x_len, num_var_map_1c);

  // dX/dS is derived by differentiating trans_Z_to_X with respect to S.
  // For the uncorrelated case, u and z are constants.  For the correlated
  // case, u is a constant, but z(s) = L(s) u due to Nataf dependence on s
  // and dz/ds = dL/ds u.

  RealVector z_vars;
  trans_X_to_Z(x_vars, z_vars);

  bool need_xs = false;
  // For distributions without simple closed-form CDFs (beta, gamma), dx/ds is
  // computed numerically.  If uncorrelated, then this is only needed if the
  // beta/gamma distribution parameters are design variables.  If correlated,
  // then the beta/gamma distribution parameters do not have to be design
  // variables (dx/ds for beta/gamma x will include a dz/ds contribution).
  if (correlationFlagX && !extendedUSpace)
    for (i=numDesignVars; i<numDesignVars+numUncertainVars; i++)
      if (ranVarTypesX[i] == BETA ||  ranVarTypesX[i] == GAMMA)
	for (j=numDesignVars; j<numDesignVars+numUncertainVars; j++)
	  if (i != j && fabs(corrMatrixX(i, j)) > 1.e-25)
	    { need_xs = true; break; }
  if ( secondaryACVarMapTargets.contains(B_ALPHA)  ||
       secondaryACVarMapTargets.contains(B_BETA)   ||
       secondaryACVarMapTargets.contains(GA_ALPHA) ||
       ( !extendedUSpace && ( secondaryACVarMapTargets.contains(B_LWR_BND) ||
			      secondaryACVarMapTargets.contains(B_UPR_BND) ||
			      secondaryACVarMapTargets.contains(GA_BETA) ) ) )
    need_xs = true;
  // the entire numerical jacobian is computed, even though only the
  // beta/gamma rows are needed
  RealMatrix num_dx_ds, num_dz_ds;
  if (need_xs || correlationFlagX)
    numerical_design_jacobian(x_vars, need_xs, num_dx_ds,
			      correlationFlagX, num_dz_ds);
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

  const IntArray&  cv_ids = iteratedModel.continuous_variable_ids();
  const IntArray& acv_ids = iteratedModel.all_continuous_variable_ids();
  for (i=0; i<num_var_map_1c; i++) { // loop over S
    size_t cv_index = cv_ids.index(acv_ids[primaryACVarMapIndices[i]]);
    // If x_dvv were passed, it would be possible to distinguish different
    // fn_grad_x components, allowing passthrough for computing fn_grad_s for
    // augmented design variables.  For now, this has to be handled separately
    // in NonDLocalReliability::dg_ds_eval() and
    // NatafTransformation::trans_grad_X_to_S().
    //if (cv_index == _NPOS) // augmented variable: define identity mapping
    //  jacobian_xs(dvv_index, i) = 1.;
    //else {
    if (cv_index != _NPOS) {
      short target2 = secondaryACVarMapTargets[i];
      for (j=0; j<x_len; j++) {      // loop over X
	// Jacobian row    = Z value = j
	// Jacobian column = S value = i
	switch (ranVarTypesX[j]) {
	case DESIGN: { // x = L + (z + 1)*(U - L)/2
	  if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	    switch (target2) {
	    case CDV_LWR_BND: // Deriv of CDV w.r.t. its Lower Bound
	      jacobian_xs(j, i) = (1. - z_vars(j))/2.; break;
	    case CDV_UPR_BND: // Deriv of CDV w.r.t. its Upper Bound
	      jacobian_xs(j, i) = (z_vars(j) + 1.)/2.; break;
	    case NO_TARGET:   // can occur for all_variables Jacobians
	      jacobian_xs(j, i) = 0.;                  break;
	    default:
	      Cerr << "Error: secondary mapping failure in NatafTransformation"
		   << "::jacobian_dX_dS()." << endl;
	      abort_handler(-1);
	      break;
	    }
	  }
	  // Deriv of CDV w.r.t. other distribution params should always be 0
	  //if (correlationFlagX)
	  //  jacobian_xs(j, i) += 
	  //    (ranVarUpperBndsX(j)-ranVarLowerBndsX(j))/2.*num_dz_ds(j, i);
	  break;
	}
	case NORMAL: { // x = z sigma + mu
	  if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	    switch (target2) {
	    case N_MEAN: // Deriv of Normal w.r.t. its Mean
	      jacobian_xs(j, i) = 1.;
	      break;
	    case N_STD_DEV: // Deriv of Normal w.r.t. its Std Deviation
	      jacobian_xs(j, i) = z_vars(j);
	      break;
	    //case N_LWR_BND: case N_UPR_BND:
	      // Normal Lower Bound - NOT SUPPORTED
	      // Normal Upper Bound - NOT SUPPORTED
	      //break;
	    case NO_TARGET: default:
	      Cerr << "Error: secondary mapping failure in NatafTransformation"
		   << "::jacobian_dX_dS()." << endl;
	      abort_handler(-1);
	      break;
	    }
	  }
	  // Deriv of Normal w.r.t. any distribution parameter
	  if (correlationFlagX)
	    jacobian_xs(j, i) += ranVarStdDevsX(j)*num_dz_ds(j, i);
	  break;
	}
	case LOGNORMAL: { // x = exp(lamba + z zeta)
	  const Real& x = x_vars(j);
	  if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	    switch (target2) {
	    case LN_MEAN: { // Deriv of Lognormal w.r.t. its Mean
	      // x = exp(lamba + z zeta)
	      const Real& mean = ranVarMeansX(j);
	      if (ranVarAddtlParamsX[j].Length()) // mean, error factor spec
		jacobian_xs(j, i) = x/mean;
	      else {
		// mean, std dev specification
		const Real& stdev = ranVarStdDevsX(j);
		Real mean_sq = mean*mean, var = stdev*stdev,
		  zeta = sqrt(log(1. + var/mean_sq));
		jacobian_xs(j, i) = x*(zeta*mean_sq + 2.*zeta*var -
		  z_vars(j)*var)/mean/zeta/(mean_sq + var);
	      }
	      break;
	    }
	    case LN_STD_DEV: { // Deriv of Lognormal w.r.t. its Std Deviation
	      const Real& mean  = ranVarMeansX(j);
	      const Real& stdev = ranVarStdDevsX(j);
	      Real mean_sq = mean*mean, var = stdev*stdev,
		zeta = sqrt(log(1. + var/mean_sq));
	      jacobian_xs(j, i) = x*stdev*(z_vars(j) - zeta)/zeta/(mean_sq+var);
	      break;
	    }
	    case LN_ERR_FACT: { // Deriv of Lognormal w.r.t. its Error Factor
	      const Real& err = ranVarAddtlParamsX[j](0);
	      jacobian_xs(j, i) = x/1.645/err*(z_vars(j) - log(err)/1.645);
	      break;
	    }
	    //case LN_LWR_BND: case LN_UPR_BND:
	      // Lognormal Lower Bound - NOT SUPPORTED
	      // Lognormal Upper Bound - NOT SUPPORTED
	      //break;
	    case NO_TARGET: default:
	      Cerr << "Error: secondary mapping failure in NatafTransformation"
		   << "::jacobian_dX_dS()." << endl;
	      abort_handler(-1);
	      break;
	    }
	  }
	  // Deriv of Lognormal w.r.t. any distribution parameter
	  if (correlationFlagX) {
	    Real cf_var = ranVarStdDevsX(j)/ranVarMeansX(j),
	      zeta = sqrt(log(1. + cf_var*cf_var));
	    jacobian_xs(j, i) += x*zeta*num_dz_ds(j, i);
	  }
	  break;
	}
	case UNIFORM: {
	  //  extendedUSpace: x = L + (z + 1)*(U - L)/2
	  // !extendedUSpace: x = L + Phi(z) (U - L)
	  const Real& z = z_vars(j);
	  if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	    switch (target2) {
	    case U_LWR_BND: // Deriv of Uniform w.r.t. its Lower Bound
	      jacobian_xs(j, i)	= (extendedUSpace) ? (1. - z)/2. : 1. - Phi(z);
	      break;
	    case U_UPR_BND: // Deriv of Uniform w.r.t. its Upper Bound
	      jacobian_xs(j, i) = (extendedUSpace) ? (z + 1.)/2. : Phi(z);
	      break;
	    // Uniform Mean          - TO DO
	    // Uniform Std Deviation - TO DO
	    case NO_TARGET: default:
	      Cerr << "Error: secondary mapping failure in NatafTransformation"
		   << "::jacobian_dX_dS()." << endl;
	      abort_handler(-1);
	      break;
	    }
	  }
	  // Deriv of Uniform w.r.t. any distribution parameter
	  // Note: extendedUSpace case should currently be zero, but may be
	  // nonzero in the future once correlation warping is more complete.
	  if (correlationFlagX)
	    jacobian_xs(j, i) += (extendedUSpace) ?
	      (ranVarUpperBndsX(j)-ranVarLowerBndsX(j))/2.*num_dz_ds(j, i) :
	      (ranVarUpperBndsX(j)-ranVarLowerBndsX(j))*phi(z)*num_dz_ds(j, i);
	  break;
	}
	case LOGUNIFORM: {
	  //  extendedUSpace: ln x = ln L + (z+1)/2 (ln U - ln L)
	  // !extendedUSpace: ln x = ln L +  Phi(z) (ln U - ln L)
	  const Real& lwr = ranVarLowerBndsX(j);
	  const Real& upr = ranVarUpperBndsX(j);
	  const Real& x = x_vars(j); const Real& z = z_vars(j);
	  if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	    switch (target2) {
	    case LU_LWR_BND: // Deriv of Loguniform w.r.t. its Lower Bound
	      jacobian_xs(j, i)
		= (extendedUSpace) ? x*(1.-z)/2./lwr : x*(1.-Phi(z))/lwr;
	      break;
	    case LU_UPR_BND: // Deriv of Loguniform w.r.t. its Upper Bound
	      jacobian_xs(j, i)
		= (extendedUSpace) ? x*(z+1.)/2./upr : x*Phi(z)/upr;
	      break;
	    // Loguniform Mean          - TO DO
	    // Loguniform Std Deviation - TO DO
	    case NO_TARGET: default:
	      Cerr << "Error: secondary mapping failure in NatafTransformation"
		   << "::jacobian_dX_dS()." << endl;
	      abort_handler(-1);
	      break;
	    }
	  }
	  // Deriv of Loguniform w.r.t. any distribution parameter
	  // Note: extendedUSpace case should currently be zero, but may be
	  // nonzero in the future once correlation warping is more complete.
	  if (correlationFlagX) // not currently supported in Nataf
	    jacobian_xs(j, i) += (extendedUSpace) ?
	      x*(log(upr)-log(lwr))/2.*num_dz_ds(j, i): 
	      x*(log(upr)-log(lwr))*phi(z)*num_dz_ds(j, i);
	  break;
	}
	case TRIANGULAR: {
	  const Real& lwr  = ranVarLowerBndsX(j);
	  const Real& mode = ranVarAddtlParamsX[j](0);
	  const Real& upr  = ranVarUpperBndsX(j);
	  const Real& x = x_vars(j); const Real& z = z_vars(j);
	  if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	    switch (target2) {
	    case T_MODE: // Triangular Mode
	      jacobian_xs(j, i) = (x < mode) ? Phi(z)*(upr-lwr)/2./(x-lwr) :
		(Phi(z)-1.)*(upr-lwr)/2./(x-upr);
	      break;
	    case T_LWR_BND: // Triangular Lower Bound
	      jacobian_xs(j, i) = (x < mode) ? 
		1. + Phi(z)*(2.*lwr-upr-mode)/2./(x-lwr) :
		(upr-mode)*(Phi(z)-1.)/2./(x-upr);
	      break;
	    case T_UPR_BND: // Triangular Upper Bound
	      jacobian_xs(j, i) = (x < mode) ? Phi(z)*(mode-lwr)/2./(x-lwr) :
		(2.*x-lwr-mode-Phi(z)*(2.*upr-lwr-mode))/2./(x-upr);
	      break;
	    // Triangular Mean          - TO DO
	    // Triangular Std Deviation - TO DO
	    case NO_TARGET: default:
	      Cerr << "Error: secondary mapping failure in NatafTransformation"
		   << "::jacobian_dX_dS()." << endl;
	      abort_handler(-1);
	      break;
	    }
	  }
	  // Deriv of Triangular w.r.t. any distribution parameter
	  if (correlationFlagX) // not currently supported for triangular
	    jacobian_xs(j, i) += (x < mode) ?
	      (upr-lwr)*(mode-lwr)*phi(z)/2./(x-lwr)*num_dz_ds(j, i) :
	      (upr-mode)*(upr-lwr)*phi(z)/2./(upr-x)*num_dz_ds(j, i);
	  break;
	}
	case EXPONENTIAL: {
	  //  extendedUSpace: x = beta*z
	  // !extendedUSpace: Phi(z) = 1. - exp(-x/beta)
	  //                  x = -beta ln(1. - Phi(z))
	  const Real& z = z_vars(j);
	  Real normcdf_comp;
	  if (!extendedUSpace)
	    normcdf_comp = 1. - Phi(z);
	  if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	    switch (target2) {
	    case E_BETA: // Exponential Beta
	      jacobian_xs(j, i) = (extendedUSpace) ? z : -log(normcdf_comp);
	      break;
	    // Exponential Mean          - TO DO
	    // Exponential Std Deviation - TO DO
	    case NO_TARGET: default:
	      Cerr << "Error: secondary mapping failure in NatafTransformation"
		   << "::jacobian_dX_dS()." << endl;
	      abort_handler(-1);
	      break;
	    }
	  }
	  // Deriv of Exponential w.r.t. any distribution parameter
	  // Note: extendedUSpace case should currently be zero, but may be
	  // nonzero in the future once correlation warping is more complete.
	  if (correlationFlagX) {
	    const Real& beta = ranVarAddtlParamsX[j](0);
	    jacobian_xs(j, i) += (extendedUSpace) ? beta*num_dz_ds(j, i) :
	      beta*phi(z)/normcdf_comp*num_dz_ds(j, i);
	  }
	  break;
	}
	case BETA: {
	  if (extendedUSpace && !need_xs) { // x = lwr + (upr - lwr)*(z+1.)/2.;
	    if (j == cv_index) {//corresp var has deriv w.r.t. its dist param
	      switch (target2) {
	      //case B_ALPHA: // numerically evaluated
	      //case B_BETA:  // numerically evaluated
	      case B_LWR_BND: // Beta Lower Bound
		jacobian_xs(j, i) = (1. - z_vars(j))/2.;
	      case B_UPR_BND: // Beta Upper Bound
		jacobian_xs(j, i) = (z_vars(j) + 1.)/2.;
	      }
	    }
	    // Deriv of Uniform w.r.t. any distribution parameter
	    // Note: extendedUSpace case should currently be zero, but may be
	    // nonzero in the future once correlation warping is more complete.
	    if (correlationFlagX)
	      jacobian_xs(j, i)
		+= (ranVarUpperBndsX(j)-ranVarLowerBndsX(j))/2.*num_dz_ds(j, i);
	  }
	  break;
	}
	case GAMMA: {
	  if (extendedUSpace && !need_xs) { // x = z*beta
	    if (j == cv_index) {//corresp var has deriv w.r.t. its dist param
	      switch (target2) {
	      //case GA_ALPHA: // numerically evaluated
	      case GA_BETA:    // Gamma Beta
		jacobian_xs(j, i) = z_vars(j);
	      }
	    }
	    // Deriv of Uniform w.r.t. any distribution parameter
	    // Note: extendedUSpace case should currently be zero, but may be
	    // nonzero in the future once correlation warping is more complete.
	    if (correlationFlagX)
	      jacobian_xs(j, i) += ranVarAddtlParamsX[j](1)*num_dz_ds(j, i);
	  }
	  break;
	}
	case GUMBEL: {
	  const Real& alpha = ranVarAddtlParamsX[j](0);
	  const Real& z = z_vars(j);
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
	      //x = x_vars(j);
	      //alpha = Pi/sqrt(6.)/ranVarStdDevsX(j);
	      //num = -alpha*(x-z);
	      //jacobian_xs(j, i) = -alpha*exp(num-exp(num))/phi(z);
	      //break;
	    // Gumbel Standard Deviation
	      //x = x_vars(j);
	      //alpha = Pi/sqrt(6.)/ranVarStdDevsX(j);
	      //num = -alpha*(x-z);
	      //jacobian_xs(j, i)
	      //  = num*exp(num-exp(num))/ranVarStdDevsX(j)/phi(z);
	      //break;
	    case NO_TARGET: default:
	      Cerr << "Error: secondary mapping failure in NatafTransformation"
		   << "::jacobian_dX_dS()." << endl;
	      abort_handler(-1);
	      break;
	    }
	  }
	  // Deriv of Gumbel w.r.t. any distribution parameter
	  if (correlationFlagX) {
	    Real normcdf = Phi(z);
	    jacobian_xs(j, i) -= phi(z)/alpha/normcdf/log(normcdf)
	                       * num_dz_ds(j, i);
	  }
	  break;
	}
	case FRECHET: {
	  const Real& alpha = ranVarAddtlParamsX[j](0);
	  const Real& z = z_vars(j);
	  if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	    switch (target2) {
	    case F_ALPHA: { // Frechet Alpha
	      // x = beta (-ln(Phi(z)))^(-1/alpha)
	      const Real& beta  = ranVarAddtlParamsX[j](1);
	      Real num = -log(Phi(z));
	      jacobian_xs(j, i) = beta/alpha/alpha*log(num)*pow(num,-1./alpha);
	      break;
	    }
	    case F_BETA: // Frechet Beta
	      jacobian_xs(j, i) = pow(-log(Phi(z)), -1./alpha);
	      break;
	    // Frechet Mean          - TO DO
	    // Frechet Std Deviation - TO DO
	    case NO_TARGET: default:
	      Cerr << "Error: secondary mapping failure in NatafTransformation"
		   << "::jacobian_dX_dS()." << endl;
	      abort_handler(-1);
	      break;
	    }
	  }
	  // Deriv of Frechet w.r.t. any distribution parameter
	  if (correlationFlagX) {
	    Real normcdf = Phi(z);
	    jacobian_xs(j, i) -= x_vars(j)*phi(z)/alpha/normcdf/log(normcdf)
	                       * num_dz_ds(j, i);
	  }
	  break;
	}
	case WEIBULL: {
	  const Real& alpha = ranVarAddtlParamsX[j](0);
	  const Real& z = z_vars(j); 
	  if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	    switch (target2) {
	    case W_ALPHA: { // Weibull Alpha
	      // x = beta (-ln(1-Phi(z)))^(1/alpha)
	      const Real& beta  = ranVarAddtlParamsX[j](1);
	      Real num = -log(1.-Phi(z));
	      jacobian_xs(j, i) = -beta/alpha/alpha*log(num)*pow(num,1./alpha);
	      break;
	    }
	    case W_BETA: // Weibull Beta
	      jacobian_xs(j, i) = pow(-log(1.-Phi(z)), 1./alpha);
	      break;
	    // Weibull Mean          - TO DO
	    // Weibull Std Deviation - TO DO
	    case NO_TARGET: default:
	      Cerr << "Error: secondary mapping failure in NatafTransformation"
		   << "::jacobian_dX_dS()." << endl;
	      abort_handler(-1);
	      break;
	    }
	  }
	  // Deriv of Weibull w.r.t. any distribution parameter
	  if (correlationFlagX) {
	    Real num = 1.-Phi(z);
	    jacobian_xs(j, i) -= x_vars(j)*phi(z)/alpha/num/log(num)
	                       * num_dz_ds(j,i);
	  }
	  break;
	}
	case STATE: { // x = L + (z + 1)*(U - L)/2
	  if (j == cv_index) { // corresp var has deriv w.r.t. its dist param
	    switch (target2) {
	    case CSV_LWR_BND: // Deriv of CSV w.r.t. its Lower Bound
	      jacobian_xs(j, i) = (1. - z_vars(j))/2.; break;
	    case CSV_UPR_BND: // Deriv of CSV w.r.t. its Upper Bound
	      jacobian_xs(j, i) = (z_vars(j) + 1.)/2.; break;
	    case NO_TARGET:   // can occur for all_variables Jacobians
	      jacobian_xs(j, i) = 0.;                  break;
	    default:
	      Cerr << "Error: secondary mapping failure in NatafTransformation"
		   << "::jacobian_dX_dS()." << endl;
	      abort_handler(-1);
	      break;
	    }
	  }
	  // Deriv of CSV w.r.t. other distribution params should always be 0
	  //if (correlationFlagX)
	  //  jacobian_xs(j, i) += 
	  //    (ranVarUpperBndsX(j)-ranVarLowerBndsX(j))/2.*num_dz_ds(j, i);
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
    int x_len = x_vars.Length();
    RealSymMatrixArray hessian_xz(x_len);
    hessian_d2X_dZ2(x_vars, hessian_xz);

    if (hessian_xu.length() != x_len)
      hessian_xu.reshape(x_len);
    for (int i=0; i<x_len; i++) {
      // d^2X/dU^2 = dX/dZ^T d^2Z/dU^2 + dZ/dU^T d^2X/dZ^2 dZ/dU
      //           = L^T d^2X/dZ^2 L
      RealMatrix hess_xzi_L(x_len, x_len);
      hess_xzi_L.Multiply('L', 1., hessian_xz[i], corrCholeskyFactorZ, 0.);
      if (hessian_xu[i].M() != x_len)
	hessian_xu[i].Shape(x_len);
      hessian_xu[i].Multiply('T', 'N', 1., corrCholeskyFactorZ, hess_xzi_L, 0.);
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
  // d^2x/dz^2 = d/dz (dx/dz) = d/dz ( phi(z)/f(x) )
  //           = (f(x) phi'(z) - phi(z) f'(x) dx/dz)/f(x)^2
  //           = -phi(z)/f(x)^2 (z f(x) + f'(x) dx/dz)
  //           = -dx/dz (z + f'(x)/f(x) dx/dz)
  //
  // This requires the additional calculation of f'(x), the derivative of
  // the PDF.  Since GSL does not provide these, PDFs are differentiated
  // below in closed form.  For cases with f'(x) = 0 (e.g., uniform), the
  // expression can be simplified to d^2x/dz^2 = -z dx/dz

  int x_len = x_vars.Length();
  if (hessian_xz.length() != x_len)
    hessian_xz.reshape(x_len);

  RealVector z_vars;
  if ( ranVarTypesX.contains(TRIANGULAR) || ranVarTypesX.contains(GUMBEL)   ||
       ranVarTypesX.contains(FRECHET)    || ranVarTypesX.contains(WEIBULL)  ||
       ( !extendedUSpace && ( ranVarTypesX.contains(UNIFORM) ||
          ranVarTypesX.contains(LOGUNIFORM)  ||
	  ranVarTypesX.contains(EXPONENTIAL) || ranVarTypesX.contains(BETA) ||
	  ranVarTypesX.contains(GAMMA) ) ) )
    trans_X_to_Z(x_vars, z_vars);

  for (int i=0; i<x_len; i++) {
    if (hessian_xz[i].M() != x_len)
      hessian_xz[i].Shape(x_len);
    // each d^2X_i/dZ^2 has a single entry on the diagonal as defined by
    // differentiation of jacobian_dX_dZ()
    switch (ranVarTypesX[i]) {
    case DESIGN: case STATE:
      hessian_xz[i](i, i) = 0.;
      break;
    case NORMAL:
      // z = (x - mu)/sigma
      hessian_xz[i](i, i) = 0.;
      break;
    case LOGNORMAL: {
      // z = (ln x - lamba)/zeta
      // dx/dz = zeta x
      // d^2x/dz^2 = zeta dx/dz = zeta^2 x
      Real cf_var = ranVarStdDevsX(i)/ranVarMeansX(i),
	zeta_sq = log(1 + cf_var*cf_var);
      hessian_xz[i](i, i) = zeta_sq*x_vars(i);
      break;
    }
    case UNIFORM: {
      //  F(x) = (x-L)/(U-L)
      //  f(x) = 1/(U-L)
      // f'(x) = 0.
      if (extendedUSpace) // linear scaling
	hessian_xz[i](i, i) = 0.;
      else { // nonlinear transformation
	const Real& z = z_vars(i);
	Real pdf = 1./(ranVarUpperBndsX(i)-ranVarLowerBndsX(i));
	hessian_xz[i](i, i) = -z*phi(z)/pdf;
      }
      break;
    }
    case LOGUNIFORM: {
      //  F(x) = (ln x - ln L)/(ln U - ln L)
      //  f(x) =  1/(ln U - ln L)/x
      // f'(x) = -1/(ln U - ln L)/x^2
      const Real& x = x_vars(i);
      Real log_range = log(ranVarUpperBndsX(i)) - log(ranVarLowerBndsX(i));
      if (extendedUSpace)
	hessian_xz[i](i, i) = x*log_range*log_range/4.;
      else {
	const Real& z = z_vars(i);
	Real pdf = 1./x/log_range, pdf_deriv = -pdf/x, dx_dz = phi(z)/pdf;
	hessian_xz[i](i, i) = -dx_dz*(z + pdf_deriv*dx_dz/pdf);
      }
      break;
    }
    case TRIANGULAR: {
      //             x < M                           x > M
      //  F(x): (x-L)^2/(U-L)/(M-L)    (M-L)/(U-L) - (x+M-2U)(x-M)/(U-L)/(U-M)
      //  f(x): 2(x-L)/(U-L)/(M-L)     2(U-x)/(U-L)/(U-M)
      // f'(x): 2/(U-L)/(M-L)          -2/(U-L)/(U-M)
      // Note: at x=M, F(x) and f(x) are continuous but f'(x) is not
      const Real& lwr  = ranVarLowerBndsX(i);
      const Real& mode = ranVarAddtlParamsX[i](0);
      const Real& upr  = ranVarUpperBndsX(i);
      const Real& z    = z_vars(i);
      const Real& x    = x_vars(i);
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
      Real dx_dz = phi(z)/pdf;
      hessian_xz[i](i, i) = -dx_dz*(z + pdf_deriv*dx_dz/pdf);
      break;
    }
    case EXPONENTIAL: {
      //  F(x) = 1. - e^(-x/beta)
      //  f(x) = e^(-x/beta) / beta
      // f'(x) = - e^(-x/beta) / beta^2
      if (extendedUSpace) // linear scaling
	hessian_xz[i](i, i) = 0.0;
      else { // nonlinear transformation
	const Real& beta = ranVarAddtlParamsX[i](0);
	const Real& z = z_vars(i); const Real& x = x_vars(i);
	Real pdf = exp(-x/beta)/beta, pdf_deriv = -pdf/beta, dx_dz = phi(z)/pdf;
	hessian_xz[i](i, i) = -dx_dz*(z + pdf_deriv*dx_dz/pdf);
      }
      break;
    }
    case BETA: {
      //  F(x) = gsl
      //  f(x) = gsl
      // f'(x) = f(x) ((alpha-1)/(x-lwr) - (beta-1)/(upr-x))
      if (extendedUSpace) // linear scaling
	hessian_xz[i](i, i) = 0.0;
#ifdef PECOS_GSL
      else { // nonlinear transformation
	const Real& alpha = ranVarAddtlParamsX[i](0);
	const Real& beta  = ranVarAddtlParamsX[i](1);
	const Real& lwr   = ranVarLowerBndsX(i);
	const Real& upr   = ranVarUpperBndsX(i);
	const Real& z = z_vars(i); const Real& x = x_vars(i);
	Real scale = upr - lwr, scaled_x = (x-lwr)/scale;
	// GSL beta passes alpha followed by beta
	Real pdf = gsl_ran_beta_pdf(scaled_x, alpha, beta)/scale;
	Real pdf_deriv = pdf*((alpha-1.)/(x-lwr) - (beta-1.)/(upr-x));
	Real dx_dz = phi(z)/pdf;
	hessian_xz[i](i, i) = -dx_dz*(z + pdf_deriv*dx_dz/pdf);
      }
#endif // PECOS_GSL
      break;
    }
    case GAMMA: {
      //  F(x) = gsl
      //  f(x) = beta^(-alpha) x^(alpha-1) e^(-x/beta) / GammaFn(alpha)
      // f'(x) = beta^(-alpha)/GammaFn(alpha) (e^(-x/beta) (alpha-1) x^(alpha-2)
      //                                       - x^(alpha-1) e^(-x/beta)/beta)
      if (extendedUSpace) // linear scaling
	hessian_xz[i](i, i) = 0.0;
#ifdef PECOS_GSL
      else { // nonlinear transformation
	const Real& alpha = ranVarAddtlParamsX[i](0);
	const Real& beta  = ranVarAddtlParamsX[i](1);
	const Real& z = z_vars(i); const Real& x = x_vars(i);
	// GSL gamma passes alpha followed by beta
	Real pdf = gsl_ran_gamma_pdf(x, alpha, beta);
	Real pdf_deriv = pow(beta,-alpha)/gsl_sf_gamma(alpha)
	               * (exp(-x/beta)*(alpha-1.)*pow(x,alpha-2.)
	               - pow(x,alpha-1.)*exp(-x/beta)/beta);
	Real dx_dz = phi(z)/pdf;
	hessian_xz[i](i, i) = -dx_dz*(z + pdf_deriv*dx_dz/pdf);
      }
#endif // PECOS_GSL
      break;
    }
    case GUMBEL: {
      //  F(x) = e^(-e^(-alpha*(x-u)))
      //  f(x) = alpha e^(-alpha*(x-u)) F(x)
      // f'(x) = alpha (e^(-alpha*(x-u)) f(x) - alpha F(x) e^(-alpha*(x-u)))
      const Real& alpha = ranVarAddtlParamsX[i](0);
      const Real& beta  = ranVarAddtlParamsX[i](1);
      const Real& z     = z_vars(i);
      Real num = exp(-alpha*(x_vars(i)-beta)), cdf = exp(-num),
	pdf = alpha*num*cdf, pdf_deriv = alpha*(num*pdf - alpha*cdf*num),
	dx_dz = phi(z)/pdf;
      hessian_xz[i](i, i) = -dx_dz*(z + pdf_deriv*dx_dz/pdf);
      break;
    }
    case FRECHET: {
      //  F(x) = e^(-(beta/x)^alpha)
      //  f(x) = F(x) alpha (beta/x)^(alpha-1) beta/x^2
      //       = F(x) alpha/beta (beta/x)^(alpha+1)
      // f'(x) = alpha/beta ((beta/x)^(alpha+1) f(x) -
      //                     F(x) (alpha+1)/beta (beta/x)^(alpha+2))
      const Real& alpha = ranVarAddtlParamsX[i](0);
      const Real& beta  = ranVarAddtlParamsX[i](1);
      const Real& z     = z_vars(i);
      Real num = beta/x_vars(i), cdf = exp(-pow(num,alpha)),
	pdf = alpha/beta*pow(num,alpha+1.)*cdf, dx_dz = phi(z)/pdf,
	pdf_deriv = alpha/beta*(pow(num,alpha+1.)*pdf -	cdf*(alpha+1.)/beta*
				pow(num,alpha+2.));
      hessian_xz[i](i, i) = -dx_dz*(z + pdf_deriv*dx_dz/pdf);
      break;
    }
    case WEIBULL: {
      //  F(x) = 1.-e^(-(x/beta)^alpha)
      //  f(x) = alpha/beta e^(-(x/beta)^alpha) (x/beta)^(alpha-1)
      // f'(x) = alpha/beta (e^(-(x/beta)^alpha) (alpha-1)/beta
      //                     (x/beta)^(alpha-2) - (x/beta)^(alpha-1) f(x))
      const Real& alpha = ranVarAddtlParamsX[i](0);
      const Real& beta  = ranVarAddtlParamsX[i](1);
      const Real& z     = z_vars(i);
      Real num = x_vars(i)/beta, num2 = exp(-pow(num,alpha)),
	pdf = alpha/beta*num2*pow(num,alpha-1.), dx_dz = phi(z)/pdf,
	pdf_deriv = alpha/beta*(num2*(alpha-1.)/beta*pow(num,alpha-2.) -
				pow(num,alpha-1.)*pdf);
      //Real cdf = 1.-num2;
      hessian_xz[i](i, i) = -dx_dz*(z + pdf_deriv*dx_dz/pdf);
      break;
    }
    }
  }
}

} // namespace Pecos

/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "FourierInverseTransformation.hpp"

#ifdef HAVE_DFFTPACK
#define ZFFTI_F77 F77_FUNC(zffti,ZFFTI)
extern "C" void ZFFTI_F77(int& n, double* wsave);

#define ZFFTB_F77 F77_FUNC(zfftb,ZFFTB)
// WJB ToDo: ask MSE about second arg -- comments suggest non-const
extern "C" void ZFFTB_F77(int& n, const void* scary_deref, double* wsave);
#endif

static const char rcsId[]="@(#) $Id: FourierInverseTransformation.cpp 4768 2007-12-17 17:49:32Z mseldre $";

//#define DEBUG


namespace Pecos {

void FourierInverseTransformation::
compute_samples(size_t num_samples, size_t seed)
{
  // sigma_j^2 = psd(omega_j)*deltaOmega_j
  size_t i, num_terms = omegaSequence.length();
  sigmaSequence.sizeUninitialized(num_terms);
  for (i=0; i<num_terms; i++)
    sigmaSequence[i] = sqrt(psdSequence[i]*deltaOmega);

  inverseSamples.shapeUninitialized(num_samples, num_terms);

  if (fourierMethod == "shinozuka_deodatis")
    compute_samples_shinozuka_deodatis(num_samples, seed);
  else if (fourierMethod == "grigoriu")
    compute_samples_grigoriu(num_samples, seed);
}


void FourierInverseTransformation::
compute_samples_shinozuka_deodatis(size_t num_samples, size_t seed)
{
  // Function to generate num_samples independent samples of a zero-mean,
  // stationary, real-valued Gaussian process using the FFT
  // algorithm. The model is (Shinozuka and Deodatis):
  //
  //          m
  // Xm(t) = SUM sqrt(2) * s_k * cos(v_k*t + Psi_k)
  //         k=1
  //
  // where  Psi_k ~ iid U(0,2 pi) are random variables
  //        v_k   = (k-1)*delv is the frequency discretization
  //        s_k^2 = g(v_k)*delv is a discretization of one-sided PSD g(v)

  /*
  // sample Psi ~ U(0,2*pi)
  rand('seed',seed);
  Psi=2*pi*rand(m,num_samples);

  // form num_samples samples of B vector
  IMAG=sqrt(-1);
  for i=1:num_samples,
    B(:,i)=sqrt(2)*s.*exp(IMAG*Psi(:,i));
  end

  // use ifft to get samples of X
  X=m*real(ifft(B));
  */

  // Generate num_terms*num_samples LHS samples for Psi ~ iid U(0, 2.*Pi).
  // This should be more efficient than generating num_terms points each time
  // within the num_samples loop.
  int num_terms = psdSequence.length(),
    num_terms_samples = num_terms*num_samples;
  RealVector uuv_l_bnds(1, false), uuv_u_bnds(1, false);
  uuv_l_bnds[0] = 0.; uuv_u_bnds[0] = 2.*Pi;
  RealMatrix Psi_samples(num_terms_samples, 1);
#ifdef HAVE_LHS
  lhsSampler.generate_uniform_samples(uuv_l_bnds, uuv_u_bnds, num_terms_samples,
				      seed, Psi_samples);
#else
  PCerr << "Error: LHS required for uniform sample generation." << std::endl;
#endif

  size_t i, j;
  ComplexArray B(num_terms); // ComplexVector fails to instantiate
  for (i=0; i<num_samples; i++) {
    for (j=0; j<num_terms; j++)
      //Real A = sigmaSequence[j]*sqrt(2.);
      //B[j] = std::complex<Real>(A*cos(Psi_ij), A*sin(Psi_ij)); // Euler
      B[j] = std::polar(sigmaSequence[j]*sqrt(2.),
			Psi_samples(j + i * num_terms, 0));
    compute_ifft_sample_set(B, i);
  }
}


void FourierInverseTransformation::
compute_samples_grigoriu(size_t num_samples, size_t seed)
{
  // Function to generate num_samples independent samples of a zero-mean,
  // stationary, real-valued Gaussian process using the FFT
  // algorithm. The model is (Grigoriu):
  //
  //          m
  // Xm(t) = SUM s_k * [V_k * cos(v_k*t) + W_k * sin(v_k*t) ]
  //         k=1
  //
  // where  V_k,W_k ~ iid N(0,1) are random variables
  //        v_k     = (k-1)*delv is the frequency discretization
  //        s_k^2   = g(v_k)*delv is a discretization of one-sided PSD g(v)

  /*
  // form num_samples samples of B vector
  randn('seed',seed);
  IMAG=sqrt(-1);
  for i=1:num_samples,
    V=randn(m,1);W=randn(m,1);  // V, W ~ iid N(0,1)
    A=s.*sqrt(V.^2 + W.^2);     // A ~ Rayleigh
    Psi=-atan2(W,V);            // Psi ~ U(-pi,pi)
    B(:,i)=A.*exp(IMAG*Psi);
  end
  // use ifft to get samples of X
  X=m*real(ifft(B));
  */

  // Generate num_terms*num_samples LHS samples for V, W ~ iid N(0,1).
  // This should be more efficient than generating num_terms points each time
  // within the num_samples loop.
  int num_terms = psdSequence.length(),
    num_terms_samples = num_terms*num_samples;
  RealVector zero_means(2, false), unit_std_devs(2, false), empty_ra;
  zero_means[0]    = zero_means[1]    = 0.;
  unit_std_devs[0] = unit_std_devs[1] = 1.;
  RealMatrix VW_samples(num_terms_samples, 2);
#ifdef HAVE_LHS
  lhsSampler.generate_normal_samples(zero_means, unit_std_devs, empty_ra,
    empty_ra, num_terms_samples, seed, VW_samples);
#else
  PCerr << "Error: LHS required for standard normal sample generation."
	<< std::endl;
#endif

  size_t i, j;
  ComplexArray B(num_terms); // ComplexVector fails to instantiate
  for (i=0; i<num_samples; i++) {
    for (j=0; j<num_terms; j++) {
      size_t ij = j + i * num_terms;
      const Real& v_ij = VW_samples(ij, 0);
      const Real& w_ij = VW_samples(ij, 1);
      //Real A = sigmaSequence[j]*sqrt(v_ij*v_ij + w_ij*w_ij); // A ~ Rayleigh
      //Real Psi = -atan2(w_ij, v_ij);                       // Psi ~ U(-pi,pi)
      //B[j] = std::complex<Real>(A*cos(Psi), A*sin(Psi));   // Euler's formula
      B[j] = std::polar(sigmaSequence[j]*sqrt(v_ij*v_ij + w_ij*w_ij),
			-atan2(w_ij, v_ij));
    }
    compute_ifft_sample_set(B, i);
  }
}


void FourierInverseTransformation::
compute_ifft_sample_set(const ComplexArray& B, size_t i)
{
  int num_terms = psdSequence.length();
#ifdef HAVE_DFFTPACK
  double* wsave = new double [4*num_terms+15];
  ZFFTI_F77(num_terms, wsave);
  ZFFTB_F77(num_terms, &B[0], wsave); // transforms in place
  delete [] wsave;
#else
  PCerr << "Error: DFFTPACK required for inverse FFT." << std::endl;
#endif
  for (size_t j=0; j<num_terms; j++)
    inverseSamples(i,j) = num_terms*B[j].real();
}

} // namespace Pecos

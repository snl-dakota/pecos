/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/** \file pecos_ifft_sd.cpp
    \brief A driver program for PECOS */

#include "DataTransformation.hpp"


/// A driver program for PECOS.

/** Generates PSD samples for a MATLAB example problem using the
    Shinozuka and Deodatis iFFT algorithm. */

int main(int argc, char* argv[])
{
  // Instantiate/initialize the data transformation instance which manages
  // the ProbabilityTransformation and BasisFunction instances.
  Pecos::DataTransformation ifft_transform("inverse_fourier");

  // Constants for this problem
  Real   vbar  = 5000.; // cut-off frequency (rad/s)
  Real   T     = 10.;   // stop time (sec)
  size_t nseed = 314;   // random number seed
  size_t ns    = 100;   // number of samples

  ifft_transform.initialize(vbar, T);
  ifft_transform.power_spectral_density("band_limited_white_noise", 2000.);
  ifft_transform.compute_samples(ns, nseed); // compute samples

  // generate samples (analytical PSD)
  //[t,G1]=gsspfft01(p,psdfun,vbar,T,nseed,ns);

  return 0;
}

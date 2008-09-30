/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/** \file pecos_ifft_g.cpp
    \brief A driver program for PECOS */

#include "DataTransformation.hpp"


/// A driver program for PECOS.

/** Generates PSD samples for a MATLAB example problem using the
    Grigoriu iFFT algorithm. */

int main(int argc, char* argv[])
{
  // Instantiate/initialize the data transformation instance which manages
  // the ProbabilityTransformation and BasisFunction instances.
  Pecos::DataTransformation ifft_transform("inverse_fourier");

  // Constants for this problem
  Real   vbar  = 10000.; // cut-off frequency (rad/s)
  Real   T     = 5.;     // stop time (sec)
  size_t nseed = 314;    // random number seed
  size_t ns    = 50;     // number of samples

  ifft_transform.initialize(vbar, T);
  ifft_transform.power_spectral_density("markov", 1000.); // user-defined psd
  ifft_transform.compute_samples(ns, nseed);

  // generate samples (analytical PSD)
  //[t,G1]=gsspfft01(p,psdfun,vbar,T,nseed,ns);

  return 0;
}

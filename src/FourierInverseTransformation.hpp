/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef FOURIER_INVERSE_TRANSFORMATION_HPP
#define FOURIER_INVERSE_TRANSFORMATION_HPP

#include "InverseTransformation.hpp"
#include "LHSDriver.hpp"


namespace Pecos {


/// Class for iFFT data transformation.

/** The FourierInverseTransformation employs an inverse fast Fourier
    transform (iFFT) to map from the frequency domain to the time domain. */

class FourierInverseTransformation: public InverseTransformation
{
public:

  //
  //- Heading: Constructors and destructor
  //

  FourierInverseTransformation(const String& data_trans_type); ///< constructor
  ~FourierInverseTransformation();                             ///< destructor

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void compute_samples(size_t num_samples, size_t seed);

private:

  //
  //- Heading: Utility routines
  //

  /// perform iFFT using Shinozuka-Deodatis algorithm
  void compute_samples_shinozuka_deodatis(size_t num_samples, size_t seed);
  /// perform iFFT using Grigoriu algorithm
  void compute_samples_grigoriu(size_t num_samples, size_t seed);
  /// use DFFTPACK or FFTW to map B vector into the i-th inverseSamples vector
  void compute_ifft_sample_set(const ComplexArray& B, size_t i);

  //
  //- Heading: Data
  //

  /// sequence of standard deviations computed from psdSequence
  RealVector sigmaSequence;

  /// iFFT approach: "shinozuka_deodatis" or "grigoriu"
  String fourierMethod;

  /// LHS wrapper for generating normal or uniform sample sets
  LHSDriver lhsSampler;
};


inline FourierInverseTransformation::
FourierInverseTransformation(const String& data_trans_type):
  fourierMethod(data_trans_type), lhsSampler("lhs", IGNORE_RANKS, false)
{ }


inline FourierInverseTransformation::~FourierInverseTransformation()
{ }

} // namespace Pecos

#endif

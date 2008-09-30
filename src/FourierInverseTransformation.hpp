/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef FOURIER_INVERSE_TRANSFORMATION_HPP
#define FOURIER_INVERSE_TRANSFORMATION_HPP

#include "InverseTransformation.hpp"


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

  FourierInverseTransformation();  ///< constructor
  ~FourierInverseTransformation(); ///< destructor

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void compute_samples(size_t num_samples, size_t seed);

private:

  //
  //- Heading: Utility routines
  //

  void compute_samples_shinozuka_deodatis(size_t num_samples, size_t seed);
  void compute_samples_grigoriu(size_t num_samples, size_t seed);

  //
  //- Heading: Data
  //

  /// sequence of standard deviations computed from psdSequence
  RealVector sigmaSequence;
};


inline FourierInverseTransformation::FourierInverseTransformation()
{ }


inline FourierInverseTransformation::~FourierInverseTransformation()
{ }

} // namespace Pecos

#endif

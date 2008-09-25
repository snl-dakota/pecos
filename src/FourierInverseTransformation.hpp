/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef FOURIER_INVERSE_TRANSFORMATION_HPP
#define FOURIER_INVERSE_TRANSFORMATION_HPP

#include "DataTransformation.hpp"


namespace Pecos {


/// Class for iFFT data transformation.

/** The FourierInverseTransformation employs an inverse fast Fourier
    transform (iFFT) to map from the frequency domain to the time domain. */

class FourierInverseTransformation: public DataTransformation
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


private:

  //
  //- Heading: Utility routines
  //

};


inline FourierInverseTransformation::FourierInverseTransformation():
  DataTransformation(BaseConstructor())
{ }


inline FourierInverseTransformation::~FourierInverseTransformation()
{ }

} // namespace Pecos

#endif

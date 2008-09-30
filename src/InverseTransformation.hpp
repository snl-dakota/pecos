/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef INVERSE_TRANSFORMATION_HPP
#define INVERSE_TRANSFORMATION_HPP

#include "DataTransformation.hpp"


namespace Pecos {


/// Class for inverse data transformation.

/** The InverseTransformation employs an inverse transform to map from
    the frequency domain to the time domain. */

class InverseTransformation: public DataTransformation
{
public:

  //
  //- Heading: Constructors and destructor
  //

  InverseTransformation();  ///< constructor
  ~InverseTransformation(); ///< destructor

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void initialize(const Real& total_t, const Real& w_bar);

  void power_spectral_density(const String& psd_name, Real param = 0.);
  //void power_spectral_density(fn_ptr);
  void power_spectral_density(const RealPairArray& psd);

  //void correlation_function(const String& fn_name, Real param = 0.);
  //void correlation_function(fn_ptr);
  //void correlation_function(const RealPairArray& corr_fn);

  //
  //- Heading: Data
  //

  /// total time window
  Real totalTime;
  /// time increment
  Real deltaTime;
  /// discretized time sequence
  RealVector timeSequence;

  /// cut-off frequency (rad/s)
  Real omegaBar;
  /// frequency increment (rad/s)
  Real deltaOmega;
  /// discretized frequency sequence
  RealVector omegaSequence;

  /// PSD sequence (frequency domain)
  RealVector psdSequence;
  /// final set of inverse samples (time domain)
  RealMatrix inverseSamples;

private:

  //
  //- Heading: Utility routines
  //

};


inline InverseTransformation::InverseTransformation():
  DataTransformation(BaseConstructor())
{ }


inline InverseTransformation::~InverseTransformation()
{ }

} // namespace Pecos

#endif

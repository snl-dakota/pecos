/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef DATA_TRANSFORMATION_HPP
#define DATA_TRANSFORMATION_HPP

#include "pecos_global_defs.hpp"
#include "ProbabilityTransformation.hpp"
#include "BasisFunction.hpp"


namespace Pecos {


/// Base class for forward/inverse transformations between time and
/// frequency domain data

/** The base class for data transformations based on forward/inverse
    mappings between the time and frequency domain based on
    spectral/FFT, Karhunen-Loeve, and sampling-based approaches. */

class DataTransformation
{
public:

  /// default constructor
  DataTransformation();
  /// standard constructor for envelope
  DataTransformation(const String& data_trans_type);
  /// copy constructor
  DataTransformation(const DataTransformation& data_trans);

  /// destructor
  virtual ~DataTransformation();

  /// assignment operator
  DataTransformation operator=(const DataTransformation& data_trans);

  //
  //- Heading: Virtual functions
  //

  /// set scalar data
  virtual void initialize(const Real& total_t, const Real& w_bar);

  /// set PSD to standard embedded function
  virtual void power_spectral_density(const String& psd_name, Real param = 0.);
  // define PSD from a user-defined function
  //virtual void power_spectral_density(fn_ptr);
  /// pass a discretized PSD directly
  virtual void power_spectral_density(const RealVector& psd);

  /// compute InverseTransformation::inverseSamples
  virtual void compute_samples(size_t num_samples, size_t seed);

  //
  //- Heading: Member functions
  //

protected:

  //
  //- Heading: Constructors
  //

  /// constructor initializes the base class part of letter classes
  /// (BaseConstructor overloading avoids infinite recursion in the
  /// derived class constructors - Coplien, p. 139)
  DataTransformation(BaseConstructor);

  //
  //- Heading: Member functions
  //


  //
  //- Heading: Data members
  //

  /// nonlinear variable transformation
  ProbabilityTransformation probTransform;
  /// set of Fourier, eigen, or polynomial basis functions
  BasisFunctionArray basisFns;

private:

  //
  //- Heading: Member functions
  //

  /// Used only by the standard envelope constructor to initialize
  /// transRep to the appropriate derived type.
  DataTransformation* get_data_trans(const String& data_trans_type);

  //
  //- Heading: Data members
  //

  /// pointer to the letter (initialized only for the envelope)
  DataTransformation* dataTransRep;
  /// number of objects sharing dataTransRep
  int referenceCount;
};

} // namespace Pecos

#endif

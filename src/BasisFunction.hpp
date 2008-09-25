/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 BasisFunction
//- Description: Base class for basis functions
//- Owner:	 Mike Eldred
//- Checked by:
//- Version:

#ifndef BASIS_FUNCTION_HPP
#define BASIS_FUNCTION_HPP

#include "pecos_global_defs.hpp"


namespace Pecos {


/// Base class for basis functions used for projection of random
/// variables through time or space

/** The base class for basis functions defined from Fourier functions,
    eigenfunctions, or polynomial functions. */

class BasisFunction
{
public:

  /// default constructor
  BasisFunction();
  /// standard constructor for envelope
  BasisFunction(const String& basis_fn_type);
  /// copy constructor
  BasisFunction(const BasisFunction& basis_fn);

  /// destructor
  virtual ~BasisFunction();

  /// assignment operator
  BasisFunction operator=(const BasisFunction& basis_fn);

  //
  //- Heading: Virtual functions
  //


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
  BasisFunction(BaseConstructor);

  //
  //- Heading: Member functions
  //


  //
  //- Heading: Data members
  //


private:

  //
  //- Heading: Member functions
  //

  /// Used only by the standard envelope constructor to initialize
  /// basisFnRep to the appropriate derived type.
  BasisFunction* get_basis_fn(const String& basis_fn_type);

  //
  //- Heading: Data members
  //

  /// pointer to the letter (initialized only for the envelope)
  BasisFunction* basisFnRep;
  /// number of objects sharing basisFnRep
  int referenceCount;
};

} // namespace Pecos

#endif

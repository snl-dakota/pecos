/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "BasisFunction.hpp"

static const char rcsId[]="@(#) $Id: BasisFunction.cpp 4768 2007-12-17 17:49:32Z mseldre $";

namespace Pecos {


/** This constructor is the one which must build the base class data
    for all derived classes.  get_basis_fn() instantiates a derived
    class letter and the derived constructor selects this base class
    constructor in its initialization list (to avoid recursion in the
    base class constructor calling get_basis_fn() again).  Since the
    letter IS the representation, its rep pointer is set to NULL (an
    uninitialized pointer causes problems in ~BasisFunction). */
BasisFunction::BasisFunction(BaseConstructor):
  basisFnRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "BasisFunction::BasisFunction(BaseConstructor) called to "
        << "build base class for letter." << std::endl;
#endif
}


/** The default constructor: basisFnRep is NULL in this case.  This
    makes it necessary to check for NULL in the copy constructor,
    assignment operator, and destructor. */
BasisFunction::BasisFunction(): basisFnRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "BasisFunction::BasisFunction() called to build empty "
        << "envelope." << std::endl;
#endif
}


/** Envelope constructor only needs to extract enough data to properly
    execute get_basis_fn, since BasisFunction(BaseConstructor)
    builds the actual base class data for the derived basis functions. */
BasisFunction::BasisFunction(const String& basis_fn_type):
  referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "BasisFunction::BasisFunction(string&) called to "
        << "instantiate envelope." << std::endl;
#endif

  // Set the rep pointer to the appropriate derived type
  basisFnRep = get_basis_fn(basis_fn_type);
  if ( !basisFnRep ) // bad type or insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize basisFnRep to the 
    appropriate derived type. */
BasisFunction* BasisFunction::get_basis_fn(const String& basis_fn_type)
{
#ifdef REFCOUNT_DEBUG
  PCout << "Envelope instantiating letter in get_basis_fn(string&)."
        << std::endl;
#endif

  if (basis_fn_type == "fourier")
    return NULL;//new FourierBasisFunction();
  else if (basis_fn_type == "eigen")
    return NULL;//new SVDLeftEigenBasisFunction();
  //else if (basis_fn_type == "orthogonal_polynomial")
  //  return new OrthogPolyBasisFunction();
  else {
    PCerr << "Error: BasisFunction type " << basis_fn_type << " not available."
	  << std::endl;
    return NULL;
  }
}


/** Copy constructor manages sharing of basisFnRep and incrementing
    of referenceCount. */
BasisFunction::BasisFunction(const BasisFunction& basis_fn)
{
  // Increment new (no old to decrement)
  basisFnRep = basis_fn.basisFnRep;
  if (basisFnRep) // Check for an assignment of NULL
    basisFnRep->referenceCount++;

#ifdef REFCOUNT_DEBUG
  PCout << "BasisFunction::BasisFunction(BasisFunction&)" << std::endl;
  if (basisFnRep)
    PCout << "basisFnRep referenceCount = " << basisFnRep->referenceCount
	  << std::endl;
#endif
}


/** Assignment operator decrements referenceCount for old basisFnRep,
    assigns new basisFnRep, and increments referenceCount for new
    basisFnRep. */
BasisFunction BasisFunction::operator=(const BasisFunction& basis_fn)
{
  if (basisFnRep != basis_fn.basisFnRep) { // normal case: old != new
    // Decrement old
    if (basisFnRep) // Check for null pointer
      if (--basisFnRep->referenceCount == 0) 
	delete basisFnRep;
    // Assign and increment new
    basisFnRep = basis_fn.basisFnRep;
    if (basisFnRep) // Check for an assignment of NULL
      basisFnRep->referenceCount++;
  }
  // else if assigning same rep, then do nothing since referenceCount
  // should already be correct

#ifdef REFCOUNT_DEBUG
  PCout << "BasisFunction::operator=(BasisFunction&)" << std::endl;
  if (basisFnRep)
    PCout << "basisFnRep referenceCount = " << basisFnRep->referenceCount
	  << std::endl;
#endif

  return *this; // calls copy constructor since returned by value
}


/** Destructor decrements referenceCount and only deletes basisFnRep
    when referenceCount reaches zero. */
BasisFunction::~BasisFunction()
{ 
  // Check for NULL pointer 
  if (basisFnRep) {
    --basisFnRep->referenceCount;
#ifdef REFCOUNT_DEBUG
    PCout << "basisFnRep referenceCount decremented to "
	  << basisFnRep->referenceCount << std::endl;
#endif
    if (basisFnRep->referenceCount == 0) {
#ifdef REFCOUNT_DEBUG
      PCout << "deleting basisFnRep" << std::endl;
#endif
      delete basisFnRep;
    }
  }
}

} // namespace Pecos

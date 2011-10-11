/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "BasisApproximation.hpp"
#include "NodalInterpPolyApproximation.hpp"
#include "HierarchInterpPolyApproximation.hpp"
#include "OrthogPolyApproximation.hpp"

static const char rcsId[]="@(#) $Id: BasisApproximation.cpp 4768 2007-12-17 17:49:32Z mseldre $";

namespace Pecos {


/** This constructor is the one which must build the base class data
    for all derived classes.  get_basis_approx() instantiates a derived
    class letter and the derived constructor selects this base class
    constructor in its initialization list (to avoid recursion in the
    base class constructor calling get_basis_approx() again).  Since the
    letter IS the representation, its rep pointer is set to NULL (an
    uninitialized pointer causes problems in ~BasisApproximation). */
BasisApproximation::BasisApproximation(BaseConstructor, size_t num_vars):
  numVars(num_vars), basisApproxRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "BasisApproximation::BasisApproximation(BaseConstructor) called to "
        << "build base class for letter." << std::endl;
#endif
}


/** The default constructor: basisApproxRep is NULL in this case.  This
    makes it necessary to check for NULL in the copy constructor,
    assignment operator, and destructor. */
BasisApproximation::BasisApproximation():
  basisApproxRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "BasisApproximation::BasisApproximation() called to build empty "
        << "envelope." << std::endl;
#endif
}


/** Envelope constructor only needs to extract enough data to properly
    execute get_basis_approx, since BasisApproximation(BaseConstructor)
    builds the actual base class data for the derived basis functions. */
BasisApproximation::
BasisApproximation(short basis_type, const UShortArray& approx_order,
		   size_t num_vars, bool use_derivs):
  referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "BasisApproximation::BasisApproximation(string&) called to "
        << "instantiate envelope." << std::endl;
#endif

  // Set the rep pointer to the appropriate derived type
  basisApproxRep
    = get_basis_approx(basis_type, approx_order, num_vars, use_derivs);
  if ( !basisApproxRep ) // bad type or insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize basisApproxRep to the 
    appropriate derived type. */
BasisApproximation* BasisApproximation::
get_basis_approx(short basis_type, const UShortArray& approx_order,
		 size_t num_vars, bool use_derivs)
{
#ifdef REFCOUNT_DEBUG
  PCout << "Envelope instantiating letter in get_basis_approx(string&)."
        << std::endl;
#endif

  switch (basis_type) {
  case GLOBAL_INTERPOLATION_POLYNOMIAL:
  case PIECEWISE_NODAL_INTERPOLATION_POLYNOMIAL:
    return new NodalInterpPolyApproximation(basis_type, num_vars, use_derivs);
    break;
  case PIECEWISE_HIERARCHICAL_INTERPOLATION_POLYNOMIAL:
    return new HierarchInterpPolyApproximation(basis_type,num_vars,use_derivs);
    break;
  case GLOBAL_ORTHOGONAL_POLYNOMIAL: case PIECEWISE_ORTHOGONAL_POLYNOMIAL:
    return new OrthogPolyApproximation(approx_order, num_vars, use_derivs);
    break;
  //case FOURIER_BASIS:
  //  return new FourierBasisApproximation(num_vars);             break;
  //case EIGEN_BASIS:
  //  return new SVDLeftEigenBasisApproximation(num_vars);        break;
  default:
    PCerr << "Error: BasisApproximation type " << basis_type
	  << " not available." << std::endl;
    return NULL; break;
  }
}


/** Copy constructor manages sharing of basisApproxRep and incrementing
    of referenceCount. */
BasisApproximation::BasisApproximation(const BasisApproximation& basis_approx)
{
  // Increment new (no old to decrement)
  basisApproxRep = basis_approx.basisApproxRep;
  if (basisApproxRep) // Check for an assignment of NULL
    basisApproxRep->referenceCount++;

#ifdef REFCOUNT_DEBUG
  PCout << "BasisApproximation::BasisApproximation(BasisApproximation&)"
	<< std::endl;
  if (basisApproxRep)
    PCout << "basisApproxRep referenceCount = "
	  << basisApproxRep->referenceCount << std::endl;
#endif
}


/** Assignment operator decrements referenceCount for old basisApproxRep,
    assigns new basisApproxRep, and increments referenceCount for new
    basisApproxRep. */
BasisApproximation BasisApproximation::
operator=(const BasisApproximation& basis_approx)
{
  if (basisApproxRep != basis_approx.basisApproxRep) { // std case: old != new
    // Decrement old
    if (basisApproxRep) // Check for null pointer
      if (--basisApproxRep->referenceCount == 0) 
	delete basisApproxRep;
    // Assign and increment new
    basisApproxRep = basis_approx.basisApproxRep;
    if (basisApproxRep) // Check for an assignment of NULL
      basisApproxRep->referenceCount++;
  }
  // else if assigning same rep, then do nothing since referenceCount
  // should already be correct

#ifdef REFCOUNT_DEBUG
  PCout << "BasisApproximation::operator=(BasisApproximation&)" << std::endl;
  if (basisApproxRep)
    PCout << "basisApproxRep referenceCount = "
	  << basisApproxRep->referenceCount << std::endl;
#endif

  return *this; // calls copy constructor since returned by value
}


/** Destructor decrements referenceCount and only deletes basisApproxRep
    when referenceCount reaches zero. */
BasisApproximation::~BasisApproximation()
{ 
  // Check for NULL pointer 
  if (basisApproxRep) {
    --basisApproxRep->referenceCount;
#ifdef REFCOUNT_DEBUG
    PCout << "basisApproxRep referenceCount decremented to "
	  << basisApproxRep->referenceCount << std::endl;
#endif
    if (basisApproxRep->referenceCount == 0) {
#ifdef REFCOUNT_DEBUG
      PCout << "deleting basisApproxRep" << std::endl;
#endif
      delete basisApproxRep;
    }
  }
}


Real BasisApproximation::value(const RealVector& x)
{
  if (!basisApproxRep) {
    PCerr << "Error: value() not available for this basis approximation "
	  << "type." << std::endl;
    abort_handler(-1);
  }

  return basisApproxRep->value(x);
}


const RealVector& BasisApproximation::gradient(const RealVector& x)
{
  if (!basisApproxRep) {
    PCerr << "Error: gradient() not available for this basis approximation "
	  << "type." << std::endl;
    abort_handler(-1);
  }

  return basisApproxRep->gradient(x);
}


const RealSymMatrix& BasisApproximation::hessian(const RealVector& x)
{
  if (!basisApproxRep) {
    PCerr << "Error: hessian() not available for this basis approximation "
	  << "type." << std::endl;
    abort_handler(-1);
  }
    
  return basisApproxRep->hessian(x);
}


int BasisApproximation::min_coefficients() const
{
  if (!basisApproxRep) { // no default implementation
    PCerr << "Error: min_coefficients() not defined for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }

  return basisApproxRep->min_coefficients(); // fwd to letter
}


void BasisApproximation::compute_coefficients()
{
  if (basisApproxRep)
    basisApproxRep->compute_coefficients(); 
  else {
    PCerr << "Error: compute_coefficients() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}


void BasisApproximation::increment_coefficients()
{
  if (basisApproxRep)
    basisApproxRep->increment_coefficients(); 
  else {
    PCerr << "Error: increment_coefficients() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}


void BasisApproximation::decrement_coefficients()
{
  if (basisApproxRep)
    basisApproxRep->decrement_coefficients(); 
  else {
    PCerr << "Error: decrement_coefficients() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}


void BasisApproximation::restore_coefficients()
{
  if (basisApproxRep)
    basisApproxRep->restore_coefficients(); 
  else {
    PCerr << "Error: restore_coefficients() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}


void BasisApproximation::finalize_coefficients()
{
  if (basisApproxRep)
    basisApproxRep->finalize_coefficients(); 
  else {
    PCerr << "Error: finalize_coefficients() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}


void BasisApproximation::store_coefficients()
{
  if (basisApproxRep)
    basisApproxRep->store_coefficients(); 
  else {
    PCerr << "Error: store_coefficients() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}


void BasisApproximation::combine_coefficients(short combine_type)
{
  if (basisApproxRep)
    basisApproxRep->combine_coefficients(combine_type);
  else {
    PCerr << "Error: combine_coefficients() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}


bool BasisApproximation::restore_available()
{ return (basisApproxRep) ? basisApproxRep->restore_available() : false; }


size_t BasisApproximation::restoration_index()
{
  if (!basisApproxRep) {
    PCerr << "Error: restoration_index() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
  return basisApproxRep->restoration_index(); 
}


size_t BasisApproximation::finalization_index(size_t i)
{
  if (!basisApproxRep) {
    PCerr << "Error: finalization_index(size_t) not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
  return basisApproxRep->finalization_index(i); 
}


void BasisApproximation::print_coefficients(std::ostream& s) const
{
  if (basisApproxRep)
    basisApproxRep->print_coefficients(s);
  else {
    PCerr << "Error: print_coefficients() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}


const RealVector& BasisApproximation::approximation_coefficients() const
{
  if (!basisApproxRep) {
    PCerr << "Error: approximation_coefficients() not available for this "
	  << "basis approximation type." << std::endl;
    abort_handler(-1);
  }
   
  return basisApproxRep->approximation_coefficients(); // fwd to letter
}


void BasisApproximation::
approximation_coefficients(const RealVector& approx_coeffs)
{
  if (basisApproxRep)
    basisApproxRep->approximation_coefficients(approx_coeffs); // fwd to letter
  else {
    PCerr << "Error: approximation_coefficients() not available for this "
	  << "basis approximation type." << std::endl;
    abort_handler(-1);
  }
}

} // namespace Pecos

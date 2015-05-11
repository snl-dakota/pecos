/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 RandomVariable
//- Description: Implementation code for RandomVariable class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "RandomVariable.hpp"
#include "NormalRandomVariable.hpp"
#include "BoundedNormalRandomVariable.hpp"
#include "LognormalRandomVariable.hpp"
#include "BoundedLognormalRandomVariable.hpp"
#include "UniformRandomVariable.hpp"
#include "LoguniformRandomVariable.hpp"
#include "TriangularRandomVariable.hpp"
#include "BetaRandomVariable.hpp"
#include "HistogramBinRandomVariable.hpp"
#include "ExponentialRandomVariable.hpp"
#include "GammaRandomVariable.hpp"
#include "GumbelRandomVariable.hpp"
#include "FrechetRandomVariable.hpp"
#include "WeibullRandomVariable.hpp"

static const char rcsId[]="@(#) $Id: RandomVariable.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


/** This constructor is the one which must build the base class data
    for all derived classes.  get_random_variable() instantiates a derived
    class letter and the derived constructor selects this base class
    constructor in its initialization list (to avoid recursion in the
    base class constructor calling get_random_variable() again).  Since the
    letter IS the representation, its rep pointer is set to NULL (an
    uninitialized pointer causes problems in ~RandomVariable). */
RandomVariable::RandomVariable(BaseConstructor):
  ranVarType(NO_TYPE), ranVarRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "RandomVariable::RandomVariable(BaseConstructor) called to build "
        << "base class for letter." << std::endl;
#endif
}


/** The default constructor: ranVarRep is NULL in this case.  This
    makes it necessary to check for NULL in the copy constructor,
    assignment operator, and destructor. */
RandomVariable::RandomVariable(): ranVarRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "RandomVariable::RandomVariable() called to build empty envelope."
	<< std::endl;
#endif
}


/** Envelope constructor only needs to extract enough data to properly
    execute get_random_variable, since RandomVariable(BaseConstructor)
    builds the actual base class data for the derived basis functions. */
RandomVariable::RandomVariable(short ran_var_type):
  referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "RandomVariable::RandomVariable(short) called to instantiate "
	<< "envelope." << std::endl;
#endif

  // Set the rep pointer to the appropriate derived type
  ranVarRep = get_random_variable(ran_var_type);
  if ( !ranVarRep ) // bad type or insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize ranVarRep to the 
    appropriate derived type. */
RandomVariable* RandomVariable::get_random_variable(short ran_var_type)
{
#ifdef REFCOUNT_DEBUG
  PCout << "Envelope instantiating letter in get_random_variable(short)."
	<< std::endl;
#endif

  RandomVariable* ran_var_rep;
  switch (ran_var_type) {
  case STD_NORMAL: case NORMAL: ran_var_rep = new NormalRandomVariable(); break;
  case BOUNDED_NORMAL:   ran_var_rep = new BoundedNormalRandomVariable(); break;
  case LOGNORMAL:        ran_var_rep = new LognormalRandomVariable();     break;
  case BOUNDED_LOGNORMAL:
    ran_var_rep = new BoundedLognormalRandomVariable();                   break;
  case STD_UNIFORM: case UNIFORM:
    ran_var_rep = new UniformRandomVariable();                            break;
  case LOGUNIFORM:          ran_var_rep = new LoguniformRandomVariable(); break;
  case TRIANGULAR:          ran_var_rep = new TriangularRandomVariable(); break;
  case STD_EXPONENTIAL: case EXPONENTIAL:
    ran_var_rep = new ExponentialRandomVariable();                        break;
  case STD_BETA:  case BETA:  ran_var_rep = new BetaRandomVariable();     break;
  case STD_GAMMA: case GAMMA: ran_var_rep = new GammaRandomVariable();    break;
  case GUMBEL:                ran_var_rep = new GumbelRandomVariable();   break;
  case FRECHET:               ran_var_rep = new FrechetRandomVariable();  break;
  case WEIBULL:               ran_var_rep = new WeibullRandomVariable();  break;
  case HISTOGRAM_BIN:
    ran_var_rep = new HistogramBinRandomVariable();                       break;
  default:
    PCerr << "Error: RandomVariable type " << ran_var_type << " not available."
	  << std::endl;
    ran_var_rep = NULL;                                                   break;
  }

  // derived classes cover multiple ranVarTypes, so best to update here:
  if (ran_var_rep)
    ran_var_rep->ranVarType = ran_var_type;

  return ran_var_rep;
}


/** Copy constructor manages sharing of ranVarRep and incrementing
    of referenceCount. */
RandomVariable::RandomVariable(const RandomVariable& ran_var)
{
  // Increment new (no old to decrement)
  ranVarRep = ran_var.ranVarRep;
  if (ranVarRep) // Check for an assignment of NULL
    ++ranVarRep->referenceCount;

#ifdef REFCOUNT_DEBUG
  PCout << "RandomVariable::RandomVariable(RandomVariable&)"
	<< std::endl;
  if (ranVarRep)
    PCout << "ranVarRep referenceCount = " << ranVarRep->referenceCount
	  << std::endl;
#endif
}


/** Assignment operator decrements referenceCount for old ranVarRep,
    assigns new ranVarRep, and increments referenceCount for new
    ranVarRep. */
RandomVariable RandomVariable::operator=(const RandomVariable& ran_var)
{
  if (ranVarRep != ran_var.ranVarRep) { // std case: old != new
    // Decrement old
    if (ranVarRep) // Check for null pointer
      if (--ranVarRep->referenceCount == 0) 
	delete ranVarRep;
    // Assign and increment new
    ranVarRep = ran_var.ranVarRep;
    if (ranVarRep) // Check for an assignment of NULL
      ++ranVarRep->referenceCount;
  }
  // else if assigning same rep, then do nothing since referenceCount
  // should already be correct

#ifdef REFCOUNT_DEBUG
  PCout << "RandomVariable::operator=(RandomVariable&)" << std::endl;
  if (ranVarRep)
    PCout << "ranVarRep referenceCount = " << ranVarRep->referenceCount
	  << std::endl;
#endif

  return *this; // calls copy constructor since returned by value
}


/** Destructor decrements referenceCount and only deletes ranVarRep
    when referenceCount reaches zero. */
RandomVariable::~RandomVariable()
{ 
  // Check for NULL pointer 
  if (ranVarRep) {
    --ranVarRep->referenceCount;
#ifdef REFCOUNT_DEBUG
    PCout << "ranVarRep referenceCount decremented to "
	  << ranVarRep->referenceCount << std::endl;
#endif
    if (ranVarRep->referenceCount == 0) {
#ifdef REFCOUNT_DEBUG
      PCout << "deleting ranVarRep" << std::endl;
#endif
      delete ranVarRep;
    }
  }
}


Real RandomVariable::cdf(Real x) const
{
  if (!ranVarRep) {
    PCerr << "Error: cdf() not supported for this random variable type."
	  << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->cdf(x); // forward to letter
}


Real RandomVariable::cdf_inverse(Real p) const
{
  if (!ranVarRep) {
    PCerr << "Error: cdf_inverse() not supported for this random variable type."
	  << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->cdf(p); // forward to letter
}


Real RandomVariable::pdf(Real x) const
{
  if (!ranVarRep) {
    PCerr << "Error: pdf() not supported for this random variable type."
	  << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->pdf(x); // forward to letter
}


Real RandomVariable::pdf_gradient(Real x) const
{
  if (!ranVarRep) {
    PCerr << "Error: pdf_gradient() not supported for this random variable "
	  << "type." << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->pdf_gradient(x); // forward to letter
}


Real RandomVariable::pdf_hessian(Real x) const
{
  if (!ranVarRep) {
    PCerr << "Error: pdf_hessian() not supported for this random variable type."
	  << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->pdf_hessian(x); // forward to letter
}

} // namespace Pecos

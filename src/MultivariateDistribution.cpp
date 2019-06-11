/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "pecos_stat_util.hpp"
#include "MarginalsCorrDistribution.hpp"
#include "MultivariateNormalDistribution.hpp"

static const char rcsId[]="@(#) $Id: MultivariateDistribution.cpp 4768 2007-12-17 17:49:32Z mseldre $";

namespace Pecos {


/** This constructor is the one which must build the base class data
    for all derived classes.  get_mv_dist() instantiates a derived
    class letter and the derived constructor selects this base class
    constructor in its initialization list (to avoid recursion in the
    base class constructor calling get_mv_dist() again).  Since the
    letter IS the representation, its rep pointer is set to NULL (an
    uninitialized pointer causes problems in ~MultivariateDistribution). */
MultivariateDistribution::MultivariateDistribution(BaseConstructor):
  correlationFlag(false), mvDistRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "MultivariateDistribution::MultivariateDistribution(Base"
        << "Constructor) called to build base class for letter." << std::endl;
#endif
}


/** The default constructor: mvDistRep is NULL in this case.  This
    makes it necessary to check for NULL in the copy constructor,
    assignment operator, and destructor. */
MultivariateDistribution::MultivariateDistribution():
  mvDistRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "MultivariateDistribution::MultivariateDistribution() called to "
        << "build empty envelope." << std::endl;
#endif
}


/** Envelope constructor only needs to extract enough data to properly
    execute get_mv_dist, since MultivariateDistribution(BaseConstructor)
    builds the actual base class data for the derived transformations. */
MultivariateDistribution::
MultivariateDistribution(short mv_dist_type): referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "MultivariateDistribution::MultivariateDistribution(string&) "
        << "called to instantiate envelope." << std::endl;
#endif

  // Set the rep pointer to the appropriate derived type
  mvDistRep = get_distribution(mv_dist_type);
  if ( !mvDistRep ) // bad type or insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize mvDistRep to the 
    appropriate derived type. */
MultivariateDistribution* MultivariateDistribution::
get_distribution(short mv_dist_type)
{
#ifdef REFCOUNT_DEBUG
  PCout << "Envelope instantiating letter in get_distribution(string&)."
	<< std::endl;
#endif

  MultivariateDistribution* mvd_rep;
  switch (mv_dist_type) {
  case MARGINALS_CORRELATIONS:
    mvd_rep = new MarginalsCorrDistribution();      break;
  case MULTIVARIATE_NORMAL:
    mvd_rep = new MultivariateNormalDistribution(); break;
  //case JOINT_KDE:
  //  mvd_rep = new JointKDEDistribution();         break;
  //case GAUSSIAN_COPULA:
  //  mvd_rep = new CopulaDistribution<Gaussian>(); break; // if templated...
  //etc.
  default:
    PCerr << "Error: MultivariateDistribution type " << mv_dist_type
	  << " not available." << std::endl;
    mvd_rep = NULL;
  }

  // some derived classes (especially template classes) cover multiple
  // ranVarTypes, so override ctor assignments for those cases:
  if (mvd_rep)
    mvd_rep->mvDistType = mv_dist_type;

  return mvd_rep;
}


/** Copy constructor manages sharing of mvDistRep and incrementing
    of referenceCount. */
MultivariateDistribution::
MultivariateDistribution(const MultivariateDistribution& mv_dist)
{
  // Increment new (no old to decrement)
  mvDistRep = mv_dist.mvDistRep;
  if (mvDistRep) // Check for an assignment of NULL
    ++mvDistRep->referenceCount;

#ifdef REFCOUNT_DEBUG
  PCout << "MultivariateDistribution::MultivariateDistribution("
        << "MultivariateDistribution&)" << std::endl;
  if (mvDistRep)
    PCout << "mvDistRep referenceCount = " << mvDistRep->referenceCount
	  << std::endl;
#endif
}


/** Assignment operator decrements referenceCount for old mvDistRep, assigns
    new mvDistRep, and increments referenceCount for new mvDistRep. */
MultivariateDistribution MultivariateDistribution::
operator=(const MultivariateDistribution& mv_dist)
{
  if (mvDistRep != mv_dist.mvDistRep) { // normal case: old != new
    // Decrement old
    if (mvDistRep) // Check for null pointer
      if (--mvDistRep->referenceCount == 0) 
	delete mvDistRep;
    // Assign and increment new
    mvDistRep = mv_dist.mvDistRep;
    if (mvDistRep) // Check for an assignment of NULL
      ++mvDistRep->referenceCount;
  }
  // else if assigning same rep, then do nothing since referenceCount
  // should already be correct

#ifdef REFCOUNT_DEBUG
  PCout << "MultivariateDistribution::operator=(MultivariateDistribution&)"
        << std::endl;
  if (mvDistRep)
    PCout << "mvDistRep referenceCount = " << mvDistRep->referenceCount
	  << std::endl;
#endif

  return *this; // calls copy constructor since returned by value
}


/** Destructor decrements referenceCount and only deletes mvDistRep
    when referenceCount reaches zero. */
MultivariateDistribution::~MultivariateDistribution()
{ 
  // Check for NULL pointer 
  if (mvDistRep) {
    --mvDistRep->referenceCount;
#ifdef REFCOUNT_DEBUG
    PCout << "mvDistRep referenceCount decremented to " 
	  << mvDistRep->referenceCount << std::endl;
#endif
    if (mvDistRep->referenceCount == 0) {
#ifdef REFCOUNT_DEBUG
      PCout << "deleting mvDistRep" << std::endl;
#endif
      delete mvDistRep;
    }
  }
}


const RandomVariable& MultivariateDistribution::random_variable(size_t i) const
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: random_variable(size_t) not supported for this "
	  << "multivariate distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->random_variable(i);
}


const std::vector<RandomVariable>& MultivariateDistribution::
random_variables() const
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: random_variables() not supported for this multivariate "
	  << "distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->random_variables();
}


std::vector<RandomVariable>& MultivariateDistribution::
random_variables()
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: random_variables() not supported for this multivariate "
	  << "distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->random_variables();
}


const ShortArray& MultivariateDistribution::types() const
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: types() not supported for this multivariate distribution "
	  << "type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->types();
}


short MultivariateDistribution::type(size_t i) const
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: type(size_t) not supported for this multivariate "
	  << "distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->type(i);
}


const BitArray& MultivariateDistribution::active_variables() const
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: active_variables() not supported for this multivariate "
	  << "distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->active_variables();
}


const RealSymMatrix& MultivariateDistribution::correlation_matrix() const
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: correlation_matrix() not supported for this multivariate "
	  << "distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->correlation_matrix();
}


const BitArray& MultivariateDistribution::active_correlations() const
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: active_correlations() not supported for this multivariate "
	  << "distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->active_correlations();
}


Real MultivariateDistribution::pdf(const RealVector& pt) const
{
  if (mvDistRep)
    return mvDistRep->pdf(pt);
  else { // forward to letter
    PCerr << "Error: pdf(RealVector) not supported for this multivariate "
	  << "distribution type." << std::endl;
    abort_handler(-1);
    return 0.;
  }
}


Real MultivariateDistribution::log_pdf(const RealVector& pt) const
{
  if (mvDistRep)
    return mvDistRep->pdf(pt);
  else // default implementation (exponential-based distribs will override)
    return std::log(pdf(pt));
}



/** This function provides a deep copy (the copy constructor supports
    shallow copies with shared reps) and is commonly used to publish
    tranformation data when the Model variables are in a transformed
    space (e.g., u-space) and x-space data may not be generated
    directly.  This allows for the use of inverse transformations to
    return the transformed space variables to their original states.
void MultivariateDistribution::
copy(const MultivariateDistribution& mv_dist)
{
  if (mvDistRep) // target is envelope
    mvDistRep->copy(mv_dist);
  else { // target is letter
    if (mv_dist.mvDistRep) { // source is envelope
      correlationFlag = mv_dist.mvDistRep->correlationFlagX;
      //...
    }
    else { // source is letter
      correlationFlag = mv_dist.correlationFlagX;
      //...
    }
  }
}
*/

} // namespace Pecos

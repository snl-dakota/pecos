/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 IntervalRandomVariable
//- Description: A random variable described by discrete values + probabilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef INTERVAL_RANDOM_VARIABLE_HPP
#define INTERVAL_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for interval random variables

/** Manages basic probability assignments (BPAs: pairings of
    intervals with probabilities) for types int and real. */

template <typename T>
class IntervalRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  IntervalRandomVariable();
  /// alternate constructor
  IntervalRandomVariable(const std::map<std::pair<T, T>, Real>& bpa);
  /// destructor
  ~IntervalRandomVariable();

  //
  //- Heading: Virtual function redefinitions
  //

  //Real mean() const;
  //Real median() const;
  //Real mode() const;
  //Real standard_deviation() const;
  //Real variance() const;
  
  //RealRealPair moments() const;
  RealRealPair bounds() const;

  //Real coefficient_of_variation() const;

  void pull_parameter(short dist_param,
		      std::map<std::pair<T, T>, Real>& bpa) const;
  void push_parameter(short dist_param,
		      const std::map<std::pair<T, T>, Real>& bpa);

  //
  //- Heading: Member functions
  //

  void update(const std::map<std::pair<T, T>, Real>& bpa);

  //
  //- Heading: Static member functions (global utilities)
  //

  //static void moments_from_params(const std::map<std::pair<T, T>, Real>& bpa,
  //				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Data
  //

  /// value-prob pairs for int values within a set
  std::map<std::pair<T, T>, Real> intervalBPA;
};


//// GENERIC ////


template <typename T>
IntervalRandomVariable<T>::IntervalRandomVariable():
  RandomVariable(BaseConstructor())
{ }


template <typename T>
IntervalRandomVariable<T>::
IntervalRandomVariable(const std::map<std::pair<T, T>, Real>& bpa):
  RandomVariable(BaseConstructor()), intervalBPA(bpa)
{ }


template <typename T>
IntervalRandomVariable<T>::~IntervalRandomVariable()
{ }


template <typename T>
void IntervalRandomVariable<T>::
update(const std::map<std::pair<T, T>, Real>& bpa)
{ intervalBPA = bpa; }
// specializations could be used for assigning ranVarType, or could employ
// std::is_same for type identification.  Simplest: ranVarType assigned at
// bottom of RandomVariable::get_random_variable().


template <typename T>
void IntervalRandomVariable<T>::
pull_parameter(short dist_param, std::map<std::pair<T, T>, Real>& bpa) const
{
  // could specialize template, but case aggregation seems adequate

  switch (dist_param) {
  case CIU_BPA: case DIU_BPA:  bpa = intervalBPA;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in IntervalRandomVariable::pull_parameter(T)." << std::endl;
    abort_handler(-1); break;
  }
}


template <typename T>
void IntervalRandomVariable<T>::
push_parameter(short dist_param, const std::map<std::pair<T, T>, Real>& bpa)
{
  // could specialize template, but case aggregation seems adequate

  switch (dist_param) {
  case CIU_BPA: case DIU_BPA:  intervalBPA = bpa;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in IntervalRandomVariable::push_parameter(T)." << std::endl;
    abort_handler(-1); break;
  }
}

/*
template <typename T>
RealRealPair IntervalRandomVariable<T>::moments() const
{
  RealRealPair moms;
  moments_from_params(intervalBPA, moms.first, moms.second);
  return moms;
}


template <typename T>
Real IntervalRandomVariable<T>::mean() const
{ return moments().first; }


//template <typename T>
//Real IntervalRandomVariable<T>::median() const
//{ return inverse_cdf(.5); } // default


template <typename T>
Real IntervalRandomVariable<T>::standard_deviation() const
{ return moments().second; }


template <typename T>
Real IntervalRandomVariable<T>::variance() const
{ Real std_dev = moments().second; return std_dev * std_dev; }


template <typename T>
Real IntervalRandomVariable<T>::coefficient_of_variation() const
{ RealRealPair mom = moments(); return mom.second / mom.first; }


template <typename T>
Real IntervalRandomVariable<T>::mode() const
{
  Real mode, mode_prob;
  typename std::map<std::pair<T, T>, Real>::const_iterator
    cit = intervalBPA.begin();
  mode = (Real)cit->first;  mode_prob = cit->second;  ++cit;
  for (; cit != intervalBPA.end(); ++cit)
    if (cit->second > mode_prob)
      { mode = (Real)cit->first;  mode_prob = cit->second; }
  return mode;
}
*/


template <typename T>
RealRealPair IntervalRandomVariable<T>::bounds() const
{
  // Note: intervals are ordered by map, which will order first by lower bound
  // and then by upper bound.  Would be sufficient to scan upper bounds, but
  // go ahead and scan both.
  typename std::map<std::pair<T, T>, Real>::const_iterator
    cit = intervalBPA.begin();
  typename std::pair<T, T> global_bnds = cit->first; // copy
  ++cit;
  for (; cit!=intervalBPA.end(); ++cit) {
    const std::pair<T, T>& bnds = cit->first;
    if (bnds.first  < global_bnds.first)  global_bnds.first  = bnds.first;
    if (bnds.second > global_bnds.second) global_bnds.second = bnds.second;
  }
  return RealRealPair((Real)global_bnds.first, (Real)global_bnds.second);
}


/*
/// for T-valued histogram, return a real-valued mean and std dev
template <typename T>
void IntervalRandomVariable<T>::
moments_from_params(const std::map<std::pair<T, T>, Real>& bpa,
		    Real& mean, Real& std_dev)
{
  // In point histogram case, (x,y) and (x,c) are equivalent since bins
  // have zero-width.  Assume normalization (prob values sum to 1.).
  mean = 0.;
  Real val, prod, raw2 = 0.;
  typename std::map<std::pair<T, T>, Real>::const_iterator cit;
  for (cit = bpa.begin(); cit != bpa.end(); ++cit) {
    val   = (Real)cit->first;
    prod  = cit->second * val; // prob * val
    mean += prod;
    raw2 += prod * val;        // prob * val^2
  }
  std_dev = std::sqrt(raw2 - mean * mean);
}
*/

//// SPECIALIZATIONS ////


} // namespace Pecos

#endif

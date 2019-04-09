/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 RangeVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef HISTOGRAM_PT_RANDOM_VARIABLE_HPP
#define HISTOGRAM_PT_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived RandomVariable class for range variables.

/** This is distinct from UniformRandomVariable in that it is used for
    non-random types without associated probability densities.  As
    such, some statistical operations are suppressed. */

template <typename T>
class RangeVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  RangeVariable();
  /// alternate constructor
  RangeVariable(T lwr, T upr);
  /// destructor
  ~RangeVariable();

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

  void pull_parameter(short dist_param, T& val) const;
  void push_parameter(short dist_param, T  val);

  //
  //- Heading: Member functions
  //

  void update(T lwr, T upr);

  //
  //- Heading: Static member functions (global utilities)
  //

  //static void moments_from_params(T lwr, T upr, Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Data
  //

  /// lower bound of range variable
  T lowerBnd;
  /// upper bound of range variable
  T upperBnd;
};


//// GENERIC ////


template <typename T>
RangeVariable<T>::RangeVariable():
  RandomVariable(BaseConstructor())
{ }


template <typename T>
RangeVariable<T>::RangeVariable(T lwr, T upr):
  RandomVariable(BaseConstructor())
{ update(lwr, upr); }


template <typename T>
RangeVariable<T>::~RangeVariable()
{ }


template <typename T>
void RangeVariable<T>::update(T lwr, T upr)
{ lowerBnd = lwr; upperBnd = upr; }
// specializations used for assigning ranVarType, but could also employ
// std::is_same for type identification


template <typename T>
void RangeVariable<T>::pull_parameter(short dist_param, T& val) const
{
  // could specialize template, but case aggregation seems adequate

  switch (dist_param) {
  case CR_LWR_BND: case DR_LWR_BND:  val = lowerBound;  break;
  case CR_UPR_BND: case DR_UPR_BND:  val = upperBound;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in RangeVariable::pull_parameter(T)." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


template <typename T>
void RangeVariable<T>::push_parameter(short dist_param, T val)
{
  // could specialize template, but case aggregation seems adequate

  switch (dist_param) {
  case CR_LWR_BND: case DR_LWR_BND:  lowerBound = val;  break;
  case CR_UPR_BND: case DR_UPR_BND:  upperBound = val;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in RangeVariable::push_parameter(T)." << std::endl;
    abort_handler(-1); break;
  }
}


/*
template <typename T>
RealRealPair RangeVariable<T>::moments() const
{
  RealRealPair moms;
  moments_from_params(lowerBound, upperBound, moms.first, moms.second);
  return moms;
}


template <typename T>
Real RangeVariable<T>::mean() const
{ return moments().first; }


//template <typename T>
//Real RangeVariable<T>::median() const
//{ return inverse_cdf(.5); } // default


template <typename T>
Real RangeVariable<T>::standard_deviation() const
{ return moments().second; }


template <typename T>
Real RangeVariable<T>::variance() const
{ Real std_dev = moments().second; return std_dev * std_dev; }


template <typename T>
Real RangeVariable<T>::coefficient_of_variation() const
{ RealRealPair mom = moments(); return mom.second / mom.first; }


template <typename T>
Real RangeVariable<T>::mode() const
{
  Real mode, mode_cnt;
  typename std::map<T, Real>::const_iterator cit = valueCountPairs.begin();
  mode = (Real)cit->first;  mode_cnt = cit->second;  ++cit;
  for (; cit != valueCountPairs.end(); ++cit)
    if (cit->second > mode_cnt)
      { mode = (Real)cit->first;  mode_cnt = cit->second; }
  return mode;
}


/// for T-valued histogram, return a real-valued mean and std dev
template <typename T>
void RangeVariable<T>::
moments_from_params(const std::map<T, Real>& vals_cnts,
		    Real& mean, Real& std_dev)
{
  // in point case, (x,y) and (x,c) are equivalent since bins have zero-width.
  // assume normalization (counts sum to 1.).
  mean = 0.;
  Real val, prod, raw2 = 0.;
  typename std::map<T, Real>::const_iterator cit;
  for (cit = vals_cnts.begin(); cit != vals_cnts.end(); ++cit) {
    val   = (Real)cit->first;
    prod  = cit->second * val; // normalized count * val
    mean += prod;
    raw2 += prod * val;        // normalized count * val^2
  }
  std_dev = std::sqrt(raw2 - mean * mean);
}
*/


template <typename T>
RealRealPair RangeVariable<T>::bounds() const
{ return RealRealPair((Real)lowerBound, (Real)upperBound); }

} // namespace Pecos

#endif

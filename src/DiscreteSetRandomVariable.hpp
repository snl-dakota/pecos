/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 DiscreteSetRandomVariable
//- Description: A random variable described by discrete values + probabilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef DISCRETE_SET_RANDOM_VARIABLE_HPP
#define DISCRETE_SET_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for discrete set random variables.

/** Manages value-probability pairings for types int, string, and real.
    String values are managed by index rather than value, requiring
    template specializations. */

template <typename T>
class DiscreteSetRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  DiscreteSetRandomVariable();
  /// alternate constructor
  DiscreteSetRandomVariable(const std::map<T, Real>& vals_probs);
  /// destructor
  ~DiscreteSetRandomVariable();

  //
  //- Heading: Virtual function redefinitions
  //

  Real mean() const;
  //Real median() const;
  Real mode() const;
  Real standard_deviation() const;
  Real variance() const;
  
  RealRealPair moments() const;
  RealRealPair bounds() const;

  Real coefficient_of_variation() const;

  void pull_parameter(short dist_param, std::map<T, Real>& val) const;
  void push_parameter(short dist_param, const std::map<T, Real>& val);

  void copy_parameters(const RandomVariable& rv);

  //
  //- Heading: Member functions
  //

  void update(const std::map<T, Real>& vals_probs);

  //
  //- Heading: Static member functions (global utilities)
  //

  static void moments_from_params(const std::map<T, Real>& vals_probs,
				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Data
  //

  /// value-prob pairs for int values within a set
  std::map<T, Real> valueProbPairs;
};


//// GENERIC ////


template <typename T>
DiscreteSetRandomVariable<T>::DiscreteSetRandomVariable():
  RandomVariable(BaseConstructor())
{ }


template <typename T>
DiscreteSetRandomVariable<T>::
DiscreteSetRandomVariable(const std::map<T, Real>& vals_probs):
  RandomVariable(BaseConstructor()), valueProbPairs(vals_probs)
{ }


template <typename T>
DiscreteSetRandomVariable<T>::~DiscreteSetRandomVariable()
{ }


template <typename T>
void DiscreteSetRandomVariable<T>::update(const std::map<T, Real>& vals_probs)
{ valueProbPairs = vals_probs; }
// specializations could be used for assigning ranVarType, or could employ
// std::is_same for type identification.  Simplest: ranVarType assigned at
// bottom of RandomVariable::get_random_variable().


template <typename T>
void DiscreteSetRandomVariable<T>::
pull_parameter(short dist_param, std::map<T, Real>& val) const
{
  // could specialize template, but case aggregation seems adequate

  switch (dist_param) {
  case H_PT_INT_PAIRS:    case H_PT_STR_PAIRS:    case H_PT_REAL_PAIRS:
  case DUSI_VALUES_PROBS: case DUSS_VALUES_PROBS: case DUSR_VALUES_PROBS:
    val = valueProbPairs; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in DiscreteSetRandomVariable::pull_parameter(T)." << std::endl;
    abort_handler(-1); break;
  }
}


template <typename T>
void DiscreteSetRandomVariable<T>::
push_parameter(short dist_param, const std::map<T, Real>& val)
{
  // could specialize template, but case aggregation seems adequate

  switch (dist_param) {
  case H_PT_INT_PAIRS:    case H_PT_STR_PAIRS:    case H_PT_REAL_PAIRS:
  case DUSI_VALUES_PROBS: case DUSS_VALUES_PROBS: case DUSR_VALUES_PROBS:
    valueProbPairs = val; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in DiscreteSetRandomVariable::push_parameter(T)." << std::endl;
    abort_handler(-1); break;
  }
}


template <typename T>
void DiscreteSetRandomVariable<T>::copy_parameters(const RandomVariable& rv)
{
  switch (ranVarType) {
  case HISTOGRAM_PT_INT:
    rv.pull_parameter(H_PT_INT_PAIRS,    valueProbPairs); break;
  case DISCRETE_UNCERTAIN_SET_INT:
    rv.pull_parameter(DUSI_VALUES_PROBS, valueProbPairs); break;
  case HISTOGRAM_PT_STRING:
    rv.pull_parameter(H_PT_STR_PAIRS,    valueProbPairs); break;
  case DISCRETE_UNCERTAIN_SET_STRING:
    rv.pull_parameter(DUSS_VALUES_PROBS, valueProbPairs); break;
  case HISTOGRAM_PT_REAL:
    rv.pull_parameter(H_PT_REAL_PAIRS,   valueProbPairs); break;
  case DISCRETE_UNCERTAIN_SET_REAL:
    rv.pull_parameter(DUSR_VALUES_PROBS, valueProbPairs); break;
  default:
    PCerr << "Error: update failure for RandomVariable type " << rv.type()
	  << " in DiscreteSetRandomVariable::copy_parameters(T)." << std::endl;
    abort_handler(-1); break;
  }
}


template <typename T>
RealRealPair DiscreteSetRandomVariable<T>::moments() const
{
  RealRealPair moms;
  moments_from_params(valueProbPairs, moms.first, moms.second);
  return moms;
}


template <typename T>
Real DiscreteSetRandomVariable<T>::mean() const
{ return moments().first; }


//template <typename T>
//Real DiscreteSetRandomVariable<T>::median() const
//{ return inverse_cdf(.5); } // default


template <typename T>
Real DiscreteSetRandomVariable<T>::standard_deviation() const
{ return moments().second; }


template <typename T>
Real DiscreteSetRandomVariable<T>::variance() const
{ Real std_dev = moments().second; return std_dev * std_dev; }


template <typename T>
Real DiscreteSetRandomVariable<T>::coefficient_of_variation() const
{ RealRealPair mom = moments(); return mom.second / mom.first; }


template <typename T>
Real DiscreteSetRandomVariable<T>::mode() const
{
  Real mode, mode_prob;
  typename std::map<T, Real>::const_iterator cit = valueProbPairs.begin();
  mode = (Real)cit->first;  mode_prob = cit->second;  ++cit;
  for (; cit != valueProbPairs.end(); ++cit)
    if (cit->second > mode_prob)
      { mode = (Real)cit->first;  mode_prob = cit->second; }
  return mode;
}


template <typename T>
RealRealPair DiscreteSetRandomVariable<T>::bounds() const
{
  RealRealPair bnds;
  bnds.first  = (Real)valueProbPairs.begin()->first;   // lower bound
  bnds.second = (Real)(--valueProbPairs.end())->first; // upper bound
  return bnds;
}


/// for T-valued histogram, return a real-valued mean and std dev
template <typename T>
void DiscreteSetRandomVariable<T>::
moments_from_params(const std::map<T, Real>& vals_probs,
		    Real& mean, Real& std_dev)
{
  // In point histogram case, (x,y) and (x,c) are equivalent since bins
  // have zero-width.  Assume normalization (prob values sum to 1.).
  mean = 0.;
  Real val, prod, raw2 = 0.;
  typename std::map<T, Real>::const_iterator cit;
  for (cit = vals_probs.begin(); cit != vals_probs.end(); ++cit) {
    val   = (Real)cit->first;
    prod  = cit->second * val; // prob * val
    mean += prod;
    raw2 += prod * val;        // prob * val^2
  }
  std_dev = std::sqrt(raw2 - mean * mean);
}


//// SPECIALIZATIONS ////


// for string vars, moments/bounds are based on weighting of set indices


template <>
inline Real DiscreteSetRandomVariable<String>::mode() const
{
  Real mode, mode_prob;
  SRMCIter cit = valueProbPairs.begin();
  mode = 0.;  mode_prob = cit->second;  ++cit;
  for (size_t index=1; cit!=valueProbPairs.end(); ++cit, ++index)
    if (cit->second > mode_prob)
      { mode = (Real)index;  mode_prob = cit->second; }
  return mode;
}


template <>
inline RealRealPair DiscreteSetRandomVariable<String>::bounds() const
{
  RealRealPair bnds;
  bnds.first  = 0.;                                 // index lower bound
  bnds.second = (Real)(valueProbPairs.size() - 1); // index upper bound
  return bnds;
}


template <>
inline void DiscreteSetRandomVariable<String>::
moments_from_params(const StringRealMap& s_prs, Real& mean, Real& std_dev)
{
  // in point case, (x,y) and (x,c) are equivalent since bins have zero-width.
  // assume normalization (probs sum to 1.).
  mean = 0.;
  Real val, prod, raw2 = 0.;  size_t index = 0;  SRMCIter cit;
  for (cit = s_prs.begin(); cit != s_prs.end(); ++cit, ++index) {
    val   = (Real)index;
    prod  = cit->second * val; // normalized prob * val
    mean += prod;
    raw2 += prod * val;        // normalized prob * val^2
  }
  std_dev = std::sqrt(raw2 - mean * mean);
}

} // namespace Pecos

#endif

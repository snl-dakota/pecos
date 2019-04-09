/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 HistogramPtRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef HISTOGRAM_PT_RANDOM_VARIABLE_HPP
#define HISTOGRAM_PT_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for discrete histogram random variables.

/** Manages value-count pairings for types int, string, and real.
    String values are managed by index rather than value, requiring
    template specializations. */

template <typename T>
class HistogramPtRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  HistogramPtRandomVariable();
  /// alternate constructor
  HistogramPtRandomVariable(const std::map<T, Real>& vals_cnts);
  /// destructor
  ~HistogramPtRandomVariable();

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

  //
  //- Heading: Member functions
  //

  template <typename T>
  void update(const std::map<T, Real>& vals_cnts);

  //
  //- Heading: Static member functions (global utilities)
  //

  static void moments_from_params(const std::map<T, Real>& vals_cnts,
				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Data
  //

  /// value-count pairs for int values within a set
  std::map<T, Real> valueCountPairs;
};


//// GENERIC ////


template <typename T>
HistogramPtRandomVariable<T>::HistogramPtRandomVariable():
  RandomVariable(BaseConstructor())
{ }


template <typename T>
HistogramPtRandomVariable<T>::
HistogramPtRandomVariable(const std::map<T, Real>& vals_cnts):
  RandomVariable(BaseConstructor())
{ update(vals_cnts); }


template <typename T>
HistogramPtRandomVariable<T>::~HistogramPtRandomVariable()
{ }


template <typename T>
void HistogramPtRandomVariable<T>::update(const std::map<T, Real>& vals_cnts)
{ valueCountPairs = vals_cnts; }
// specializations used for assigning ranVarType, but could also employ
// std::is_same for type identification


template <typename T>
void HistogramPtRandomVariable<T>::
pull_parameter(short dist_param, std::map<T, Real>& val) const
{
  // could specialize template, but seems adequate

  switch (dist_param) {
  case H_PT_INT_PAIRS: case H_PT_STR_PAIRS: case H_PT_REAL_PAIRS:
    val = valueCountPairs; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in HistogramPtRandomVariable::pull_parameter(T)." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


template <typename T>
void HistogramPtRandomVariable<T>::
push_parameter(short dist_param, const std::map<T, Real>& val) // *** TO DO: proliferate
{
  // could specialize template, but seems adequate

  switch (dist_param) {
  case H_PT_INT_PAIRS: case H_PT_STR_PAIRS: case H_PT_REAL_PAIRS:
    valueCountPairs = val; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in HistogramPtRandomVariable::push_parameter(T)." << std::endl;
    abort_handler(-1); break;
  }
}


template <typename T>
RealRealPair HistogramPtRandomVariable<T>::moments() const
{
  RealRealPair moms;
  moments_from_params(valueCountPairs, moms.first, moms.second);
  return moms;
}


template <typename T>
Real HistogramPtRandomVariable<T>::mean() const
{ return moments().first; }


//template <typename T>
//Real HistogramPtRandomVariable<T>::median() const
//{ return inverse_cdf(.5); } // default


template <typename T>
Real HistogramPtRandomVariable<T>::standard_deviation() const
{ return moments().second; }


template <typename T>
Real HistogramPtRandomVariable<T>::variance() const
{ Real std_dev = moments().second; return std_dev * std_dev; }


template <typename T>
Real HistogramPtRandomVariable<T>::coefficient_of_variation() const
{ RealRealPair mom = moments(); return mom.second / mom.first; }


template <typename T>
Real HistogramPtRandomVariable<T>::mode() const
{
  Real mode, mode_cnt;
  typename std::map<T, Real>::const_iterator cit = valueCountPairs.begin();
  mode = (Real)cit->first;  mode_cnt = cit->second;  ++cit;
  for (; cit != valueCountPairs.end(); ++cit)
    if (cit->second > mode_cnt)
      { mode = (Real)cit->first;  mode_cnt = cit->second; }
  return mode;
}


template <typename T>
RealRealPair HistogramPtRandomVariable<T>::bounds() const
{
  RealRealPair bnds;
  bnds.first  = (Real)valueCountPairs.begin()->first;   // lower bound
  bnds.second = (Real)(--valueCountPairs.end())->first; // upper bound
  return bnds;
}


/// for T-valued histogram, return a real-valued mean and std dev
template <typename T>
void HistogramPtRandomVariable<T>::
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


//// SPECIALIZATIONS ////


template <>
void HistogramPtRandomVariable<int>::update(const IntRealMap& vals_cnts)
{ valueCountPairs = vals_cnts;  ranVarType = HISTOGRAM_PT_INT; }


template <>
void HistogramPtRandomVariable<String>::update(const StringRealMap& vals_cnts)
{ valueCountPairs = vals_cnts;  ranVarType = HISTOGRAM_PT_STRING; }


template <>
void HistogramPtRandomVariable<Real>::update(const RealRealMap& vals_cnts)
{ valueCountPairs = vals_cnts;  ranVarType = HISTOGRAM_PT_REAL; }


// for string vars, moments/bounds are based on weighting of set indices


template <>
Real HistogramPtRandomVariable<String>::mode() const
{
  Real mode, mode_cnt;
  SRMCIter cit = valueCountPairs.begin();
  mode = 0.;  mode_cnt = cit->second;  ++cit;
  for (size_t index=1; cit!=valueCountPairs.end(); ++cit, ++index)
    if (cit->second > mode_cnt)
      { mode = (Real)index;  mode_cnt = cit->second; }
  return mode;
}


template <>
RealRealPair HistogramPtRandomVariable<String>::bounds() const
{
  RealRealPair bnds;
  bnds.first  = 0.;                                 // index lower bound
  bnds.second = (Real)(valueCountPairs.size() - 1); // index upper bound
  return bnds;
}


template <>
void HistogramPtRandomVariable<String>::
moments_from_params(const StringRealMap& s_prs, Real& mean, Real& std_dev)
{
  // in point case, (x,y) and (x,c) are equivalent since bins have zero-width.
  // assume normalization (counts sum to 1.).
  mean = 0.;
  Real val, prod, raw2 = 0.;  size_t index = 0;  SRMCIter cit;
  for (cit = s_prs.begin(); cit != s_prs.end(); ++cit, ++index) {
    val   = (Real)index;
    prod  = cit->second * val; // normalized count * val
    mean += prod;
    raw2 += prod * val;        // normalized count * val^2
  }
  std_dev = std::sqrt(raw2 - mean * mean);
}

} // namespace Pecos

#endif

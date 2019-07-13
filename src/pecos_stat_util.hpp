/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_STAT_UTIL_HPP
#define PECOS_STAT_UTIL_HPP

#include "pecos_data_types.hpp"
#include "BoundedNormalRandomVariable.hpp"
#include "BoundedLognormalRandomVariable.hpp"
#include "LoguniformRandomVariable.hpp"
#include "TriangularRandomVariable.hpp"
#include "BetaRandomVariable.hpp"
#include "GammaRandomVariable.hpp"
#include "GumbelRandomVariable.hpp"
#include "FrechetRandomVariable.hpp"
#include "WeibullRandomVariable.hpp"
#include "HistogramBinRandomVariable.hpp"
#include "InvGammaRandomVariable.hpp"
#include "PoissonRandomVariable.hpp"
#include "BinomialRandomVariable.hpp"
#include "NegBinomialRandomVariable.hpp"
#include "GeometricRandomVariable.hpp"
#include "HypergeometricRandomVariable.hpp"
#include "RangeVariable.hpp"
#include "SetVariable.hpp"
#include "IntervalRandomVariable.hpp"
#include "DiscreteSetRandomVariable.hpp"


namespace Pecos {

void bins_to_xy_pairs(const RealRealMap& h_bin_prs,
		      RealArray& x_val, RealArray& y_val);

void intervals_to_xy_pairs(const RealRealPairRealMap& ci_bpa,
			   RealArray& x_val, RealArray& y_val);

void intervals_to_xy_pairs(const IntIntPairRealMap& di_bpa,
			   RealArray& x_val, RealArray& y_val);


inline Real gamma_function(Real x)
{ return bmth::tgamma(x); }


inline void int_range_to_xy_pairs(int l_bnd,        int u_bnd,
				  RealArray& x_val, RealArray& y_val)
{
  // supports either discrete integer range or range of set indices

  int i, num_params = u_bnd - l_bnd + 1;
  x_val.resize(num_params);  y_val.assign(num_params, 1.);
  for (i=0; i<num_params; ++i)
    x_val[i] = (Real)(l_bnd + i);
}


template <typename T>
void set_to_xy_pairs(const std::set<T>& values,
		     RealArray& x_val, RealArray& y_val)
{
  int i, num_params = values.size();
  x_val.resize(num_params);  y_val.assign(num_params, 1.);
  typename std::set<T>::const_iterator cit;
  for (cit=values.begin(), i=0; cit!=values.end(); ++cit, ++i)
    x_val[i] = (Real)(*cit); // value
}


template <typename T>
void map_to_xy_pairs(const std::map<T, Real>& vals_probs,
		     RealArray& x_val, RealArray& y_val)
{
  int i, num_params = vals_probs.size();
  x_val.resize(num_params);  y_val.resize(num_params);
  typename std::map<T, Real>::const_iterator cit;
  for (cit=vals_probs.begin(), i=0; cit!=vals_probs.end(); ++cit, ++i) {
    x_val[i] = (Real)cit->first;  // value
    y_val[i] =       cit->second; // probability
  }
}


template <typename T>
void map_indices_to_xy_pairs(const std::map<T, Real>& vals_probs,
			     RealArray& x_val, RealArray& y_val)
{
  int i, num_params = vals_probs.size();
  x_val.resize(num_params);  y_val.resize(num_params);
  typename std::map<T, Real>::const_iterator cit;
  for (cit=vals_probs.begin(), i=0; cit!=vals_probs.end(); ++cit, ++i) {
    x_val[i] = (Real)i;     // index rather than value
    y_val[i] = cit->second; // probability
  }
}


inline void
design_state_subset(const std::vector<RandomVariable>& random_vars,
		    BitArray& subset, size_t start_set, size_t num_set)
{
  size_t i, num_rv = random_vars.size(), end_set = start_set + num_set;
  subset.resize(num_rv, false); // init bits to false
  for (i=start_set; i<end_set; ++i)
    // activate design + state vars
    switch (random_vars[i].type()) {
    case CONTINUOUS_RANGE:    case DISCRETE_RANGE: case DISCRETE_SET_INT:
    case DISCRETE_SET_STRING: case DISCRETE_SET_REAL:
      subset.set(i); break;
    }
}


inline void
uncertain_subset(const std::vector<RandomVariable>& random_vars,
		 BitArray& subset)
{
  size_t i, num_rv = random_vars.size();
  subset.resize(num_rv, true); // init bits to true
  for (i=0; i<num_rv; ++i)
    // deactivate complement of uncertain vars
    switch (random_vars[i].type()) {
    case CONTINUOUS_RANGE:    case DISCRETE_RANGE: case DISCRETE_SET_INT:
    case DISCRETE_SET_STRING: case DISCRETE_SET_REAL:
      subset.reset(i); break;
    }
}


inline void
aleatory_uncertain_subset(const std::vector<RandomVariable>& random_vars,
			  BitArray& subset)
{
  size_t i, num_rv = random_vars.size();
  subset.resize(num_rv, true); // init bits to true
  for (i=0; i<num_rv; ++i)
    // deactivate complement of aleatory uncertain vars
    switch (random_vars[i].type()) {
    case CONTINUOUS_RANGE:    case DISCRETE_RANGE: case DISCRETE_SET_INT:
    case DISCRETE_SET_STRING: case DISCRETE_SET_REAL:
    case CONTINUOUS_INTERVAL_UNCERTAIN: case DISCRETE_INTERVAL_UNCERTAIN:
    case DISCRETE_UNCERTAIN_SET_INT:    case DISCRETE_UNCERTAIN_SET_STRING:
    case DISCRETE_UNCERTAIN_SET_REAL:
      subset.reset(i); break;
    }
}


inline void
epistemic_uncertain_subset(const std::vector<RandomVariable>& random_vars,
			   BitArray& subset)
{
  size_t i, num_rv = random_vars.size();
  subset.resize(num_rv, false); // init bits to false
  for (i=0; i<num_rv; ++i)
    // activate epistemic uncertain vars
    switch (random_vars[i].type()) {
    case CONTINUOUS_INTERVAL_UNCERTAIN: case DISCRETE_INTERVAL_UNCERTAIN:
    case DISCRETE_UNCERTAIN_SET_INT:    case DISCRETE_UNCERTAIN_SET_STRING:
    case DISCRETE_UNCERTAIN_SET_REAL:
      subset.set(i); break;
    }
}

} // namespace Pecos

#endif

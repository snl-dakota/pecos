/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_STAT_UTIL_HPP
#define PECOS_STAT_UTIL_HPP

#include "pecos_data_types.hpp"
#include <boost/math/distributions.hpp>
#include <boost/math/special_functions/sqrt1pm1.hpp> // includes expm1,log1p
#include <boost/random/uniform_real.hpp>

namespace bmth = boost::math;
namespace bmp  = bmth::policies;


namespace Pecos {

// -----------------------------------
// Non-default boost math/policy types
// -----------------------------------

// continuous random variable types:
typedef bmth::
  normal_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  normal_dist;
typedef bmth::
  lognormal_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  lognormal_dist;
typedef bmth::
  triangular_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  triangular_dist;
typedef bmth::
  exponential_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  exponential_dist;
typedef bmth::
  beta_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  beta_dist;
typedef bmth::
  gamma_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  gamma_dist;
typedef bmth::
  inverse_gamma_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  inv_gamma_dist;
typedef bmth::
  extreme_value_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  extreme_value_dist;
typedef bmth::
  weibull_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  weibull_dist;
// discrete random variable types:
typedef bmth::
  poisson_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  poisson_dist;
typedef bmth::
  binomial_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  binomial_dist;
typedef bmth::
  negative_binomial_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  negative_binomial_dist;
typedef bmth::
  geometric_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  geometric_dist;
typedef bmth::
  hypergeometric_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  hypergeometric_dist;
// distributions used in statistical utilities (e.g., confidence intervals):
typedef bmth::
  chi_squared_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  chi_squared_dist;
typedef bmth::
  students_t_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  students_t_dist;
typedef bmth::
  fisher_f_distribution< Real,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  fisher_f_dist;


// -------------------------
// Variable subset utilities
// -------------------------

class RandomVariable; // fwd declaration


/// define BitArray corresponding to design/state subset of random variables
void design_state_subset(const std::vector<RandomVariable>& random_vars,
			 BitArray& subset, size_t start_set, size_t num_set);
/// define BitArray corresponding to uncertain subset of random variables
void uncertain_subset(const std::vector<RandomVariable>& random_vars,
		      BitArray& subset);
/// define BitArray corresponding to aleatory uncertain subset of random vars
void aleatory_uncertain_subset(const std::vector<RandomVariable>& random_vars,
			       BitArray& subset);
/// define BitArray corresponding to epistemic uncertain subset of random vars
void epistemic_uncertain_subset(const std::vector<RandomVariable>& random_vars,
				BitArray& subset);


// ----------------------------------
// Distribution parameter conversions
// ----------------------------------

inline void bins_to_xy_pairs(const RealRealMap& h_bin_prs,
			     RealArray& x_val, RealArray& y_val)
{
  // histogram bins: pairs are defined from an abscissa in the first field
  // and a count (not a density) in the second field.  This distinction is
  // important for unequal bin widths.

  size_t i, num_params = h_bin_prs.size(), end = num_params - 1;
  RRMCIter cit = h_bin_prs.begin();

  // Assume already normalized with sum = 1
  //Real sum = 0.;
  //RRMCIter end = --h_bin_prs.end(); // last y from DAKOTA must be zero
  //for (; cit!=end; ++cit)
  //  sum += cit->second;

  // LHS requires accumulation of CDF with first y at 0 and last y at 1
  x_val.resize(num_params);  y_val.resize(num_params);
  y_val[0] = 0.;
  for (i=0; i<end; ++i, ++cit) {
    x_val[i]   = cit->first;
    y_val[i+1] = y_val[i] + cit->second/* /sum */;
  }
  x_val[end] = cit->first; // last prob value (cit->second) must be zero
}


/// generic case (for integer types as x values)
template <typename T>
void intervals_to_xy_pairs(const std::map<std::pair<T,T>, Real>& di_bpa,
			   std::vector<T>& x_val, RealArray& y_val)
{
  // x_sort_unique contains ALL of the unique integer values for this
  // x_sort_unique contains ALL of the unique integer values for this
  // discrete interval variable in increasing order.  For example, if
  // there are 3 intervals for a variable and the bounds are (1,4),
  // (3,6), and (9,10), x_sorted will be (1,2,3,4,5,6,9,10).
  typename std::map<std::pair<T,T>, Real>::const_iterator cit;
  typename std::set<T> x_sort_unique;
  T l_bnd, u_bnd, val;
  for (cit=di_bpa.begin(); cit!=di_bpa.end(); ++cit) {
    const std::pair<T,T>& bounds = cit->first;
    u_bnd = bounds.second;
    for (val=bounds.first; val<=u_bnd; ++val)
      x_sort_unique.insert(val);
  }
  // copy sorted set<T> to x_val
  size_t j, num_params = x_sort_unique.size(), index;
  x_val.resize(num_params);
  typename std::set<T>::iterator it = x_sort_unique.begin();
  for (j=0; j<num_params; ++j, ++it)
    x_val[j] = *it;

  // Calculate probability densities and account for overlapping intervals.
  // Loop over the original intervals and see where they fall relative to
  // the new, sorted intervals for the density calculation.
  y_val.assign(num_params, 0.);
  for (cit=di_bpa.begin(); cit!=di_bpa.end(); ++cit) {
    const std::pair<T,T>& bounds = cit->first;
    l_bnd = bounds.first;  u_bnd = bounds.second;
    Real di_density = cit->second / (u_bnd - l_bnd + 1); // prob/#integers
    it = x_sort_unique.find(l_bnd);
    if (it == x_sort_unique.end()) {
      PCerr << "Error: lower bound not found in sorted set within LHSDriver "
	    << "mapping of discrete interval uncertain variable."<< std::endl;
      abort_handler(-1);
    }
    index = std::distance(x_sort_unique.begin(), it);
    for (val=l_bnd; val<=u_bnd; ++val, ++index)
      y_val[index] += di_density;
  }
#ifdef DEBUG
  for (j=0; j<num_params; ++j)
    PCout << "diuv: x_val[" << j << "] is " << x_val[j]
	  << " y_val[" << j << "] is " << y_val[j] << '\n';
#endif // DEBUG
}


/// specialization for Real x values
template <>
inline void intervals_to_xy_pairs(const RealRealPairRealMap& ci_bpa,
				  RealArray& x_val, RealArray& y_val)
{
  // x_sort_unique is a set with ALL of the interval bounds for this variable
  // in increasing order and unique.  For example, if there are 2 intervals
  // for a variable, and the bounds are (1,4) and (3,6), x_sorted will be
  // (1, 3, 4, 6).  If the intervals are contiguous, e.g. one interval is
  // (1,3) and the next is (3,5), x_sort_unique is (1,3,5).
  RRPRMCIter cit;  RealSet x_sort_unique;
  for (cit=ci_bpa.begin(); cit!=ci_bpa.end(); ++cit) {
    const RealRealPair& bounds = cit->first;
    x_sort_unique.insert(bounds.first);
    x_sort_unique.insert(bounds.second);
  }
  // convert sorted RealSet to x_val
  size_t j, num_params = x_sort_unique.size();
  x_val.resize(num_params);  y_val.resize(num_params);
  RSIter it = x_sort_unique.begin();
  for (j=0; j<num_params; ++j, ++it)
    x_val[j] = *it;

  // Calculate the probability densities, and account for the cases where
  // there are intervals that are overlapping.  This section of code goes
  // through the original intervals and see where they fall relative to the
  // new, sorted intervals for the density calculation.
  RealVector prob_dens(num_params); // initialize to 0.
  for (cit=ci_bpa.begin(); cit!=ci_bpa.end(); ++cit) {
    const RealRealPair& bounds = cit->first;
    Real l_bnd = bounds.first, u_bnd = bounds.second;
    Real ci_density = cit->second / (u_bnd - l_bnd);
    int cum_int_index = 0;
    while (l_bnd > x_val[cum_int_index])
      ++cum_int_index;
    ++cum_int_index;
    while (cum_int_index < num_params && x_val[cum_int_index] <= u_bnd)
      { prob_dens[cum_int_index] += ci_density; ++cum_int_index; }
  }

  // put the densities in a cumulative format necessary for LHS histograms.
  // Note that x_val and y_val are defined as Real* for input to f77.
  y_val[0] = 0.;
  for (j=1; j<num_params; ++j)
    y_val[j] = (prob_dens[j] > 0.0) ?
      y_val[j-1] + prob_dens[j] * (x_val[j] - x_val[j-1]) :
      y_val[j-1] + 0.0001; // handle case where there is a gap
  // normalize if necessary
  if (y_val[num_params-1] != 1.) {
    Real y_total = y_val[num_params-1];
    for (j=1; j<num_params; ++j)
      y_val[j] /= y_total;
  }
#ifdef DEBUG
  for (j=0; j<num_params; ++j)
    PCout << "ciuv: x_val[" << j << "] is " << x_val[j]
	  << " y_val[" << j << "] is " << y_val[j] << '\n';
#endif // DEBUG
}


template <typename T>
void intervals_to_xy_map(const std::map<std::pair<T,T>, Real>& bpa,
			 std::map<T, Real>& xy_vals)
{
  // existing conversion for continuous / discrete BPA
  std::vector<T> x_vals;  RealArray y_vals;
  intervals_to_xy_pairs(bpa, x_vals, y_vals);

  // flatten arrays into a single map
  size_t i, len = x_vals.size();
  for (i=0; i<len; ++i)
    xy_vals[x_vals[i]] = y_vals[i];
}


template <typename T>
void intervals_to_xy_vector(const std::map<std::pair<T,T>, Real>& bpa,
			    RealVector& xy_vals)
{
  // existing conversion for continuous / discrete BPA
  std::vector<T> x_vals;  RealArray y_vals;
  intervals_to_xy_pairs(bpa, x_vals, y_vals);

  // flatten arrays into a single vector
  size_t i, len = x_vals.size(), cntr = 0;
  xy_vals.sizeUninitialized(2*len);
  for (i=0; i<len; ++i) {
    xy_vals[cntr] = (Real)x_vals[i];  ++cntr;
    xy_vals[cntr] =       y_vals[i];  ++cntr;
  }
}


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

} // namespace Pecos

#endif

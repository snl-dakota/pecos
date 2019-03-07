/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_MATH_UTIL_HPP
#define PECOS_MATH_UTIL_HPP

#include "pecos_data_types.hpp"


namespace Pecos {

/** Append to combined_mi based on append_mi. */
inline void append_multi_index(const UShort2DArray& append_mi,
			       UShort2DArray& combined_mi)
{
  if (combined_mi.empty())
    combined_mi = append_mi;
  else {
    size_t i, num_app_mi = append_mi.size();
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = append_mi[i];
      if (std::find(combined_mi.begin(), combined_mi.end(), search_mi) ==
	  combined_mi.end())
	combined_mi.push_back(search_mi);
    }
  }
}


/** Append to combined_mi based on append_mi. */
inline void append_multi_index(const UShortArraySet& append_mi,
			       UShort2DArray& combined_mi)
{
  UShortArraySet::const_iterator cit;
  for (cit=append_mi.begin(); cit!=append_mi.end(); ++cit) {
    const UShortArray& search_mi = *cit;
    if (std::find(combined_mi.begin(), combined_mi.end(), search_mi) ==
	combined_mi.end())
      combined_mi.push_back(search_mi);
  }
}


/** Append append_mi to combined_mi, and update append_mi_map
    (SizetArray) and append_mi_map_ref to facilitate related
    aggregations without repeated searching. */
inline void append_multi_index(const UShort2DArray& append_mi,
			       UShort2DArray& combined_mi,
			       SizetArray& append_mi_map,
			       size_t& append_mi_map_ref)
{
  size_t i, num_app_mi = append_mi.size();
  append_mi_map.resize(num_app_mi);
  if (combined_mi.empty()) {
    combined_mi = append_mi;
    append_mi_map_ref = 0;
    for (i=0; i<num_app_mi; ++i)
      append_mi_map[i] = i;
  }
  else {
    append_mi_map_ref = combined_mi.size();
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = append_mi[i];
      size_t index = find_index(combined_mi, search_mi);
      if (index == _NPOS) { // search_mi does not yet exist in multi_index
	append_mi_map[i] = combined_mi.size();
	combined_mi.push_back(search_mi);
      }
      else
	append_mi_map[i] = index;
    }
  }
}


/** Append append_mi to combined_mi, and update append_mi_map
    (SizetSet) and append_mi_map_ref to facilitate related
    aggregations without repeated searching. */
inline void append_multi_index(const UShort2DArray& append_mi,
			       UShort2DArray& combined_mi,
			       SizetSet& append_mi_map,
			       size_t& append_mi_map_ref)
{
  size_t i, num_app_mi = append_mi.size();
  append_mi_map.clear();
  if (combined_mi.empty()) {
    combined_mi = append_mi;
    append_mi_map_ref = 0;
    for (i=0; i<num_app_mi; ++i)
      append_mi_map.insert(i);
  }
  else {
    append_mi_map_ref = combined_mi.size();
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = append_mi[i];
      size_t index = find_index(combined_mi, search_mi);
      if (index == _NPOS) { // search_mi does not yet exist in multi_index
	append_mi_map.insert(combined_mi.size());
	combined_mi.push_back(search_mi);
      }
      else
	append_mi_map.insert(index);
    }
  }
}


/** Append to combined_mi based on append_mi and previously defined
    append_mi_map and append_mi_map_ref.  If necessary, update
    append_mi_map and append_mi_map_ref. */
inline void append_multi_index(const UShort2DArray& append_mi,
			       SizetArray& append_mi_map,
			       size_t& append_mi_map_ref,
			       UShort2DArray& combined_mi)
{
  if (combined_mi.empty())
    combined_mi = append_mi; // assume append_mi_map{,_ref} are up to date
  else {
    size_t i, num_app_mi = append_mi.size(), num_mi = combined_mi.size();
    if (num_mi == append_mi_map_ref) { // current mi corresponds to ref
      for (i=0; i<num_app_mi; ++i)
	if (append_mi_map[i] >= append_mi_map_ref)
	  combined_mi.push_back(append_mi[i]);
    }
    else if (num_mi > append_mi_map_ref) { // mi has grown since ref taken
      for (i=0; i<num_app_mi; ++i)
	if (append_mi_map[i] >= append_mi_map_ref) { // previously appended
	  const UShortArray& search_mi = append_mi[i];
	  // search from reference pt forward
	  UShort2DArray::iterator it, it_start = combined_mi.begin();
	  std::advance(it_start, append_mi_map_ref);
	  it = std::find(it_start, combined_mi.end(), search_mi);
	  if (it == combined_mi.end()) { // still an append: update map, append
	    append_mi_map[i] = combined_mi.size();
	    combined_mi.push_back(append_mi[i]);
	  }
	  else // no longer an append: only update map
	    append_mi_map[i] = append_mi_map_ref + std::distance(it_start, it);
	}
      append_mi_map_ref = num_mi; // reference point now updated
    }
    else { // combined_mi is not allowed to shrink since ref taken
      PCerr << "Error: combined_mi inconsistent with reference size in "
	    << "append_multi_index()." << std::endl;
      abort_handler(-1);
    }
  }
}


// The following variants maintain a separation between ref_mi and combined_mi,
// rather than updating in place.  In current use cases, append_mi_map has
// provided sufficient bookkeeping to allow in-place updating.

/*  Create combined_mi by appending append_mi to ref_mi.
inline void append_multi_index(const UShort2DArray& ref_mi,
                               const UShort2DArray& append_mi,
		               UShort2DArray& combined_mi)
{
  if (ref_mi.empty())
    combined_mi = append_mi;
  else {
    combined_mi = ref_mi;
    size_t i, num_app_mi = append_mi.size();
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = append_mi[i];
      if (std::find(combined_mi.begin(), combined_mi.end(),
		    search_mi) == combined_mi.end())
	combined_mi.push_back(search_mi);
    }
  }
}
*/

/*  Append append_mi to ref_mi to create combined_mi, and update
    append_mi_map and append_mi_map_ref to facilitate related
    aggregations without repeated searching.
inline void append_multi_index(const UShort2DArray& ref_mi,
                               const UShort2DArray& append_mi,
		               UShort2DArray& combined_mi,
                               SizetArray& append_mi_map,
			       size_t& append_mi_map_ref)
{
  size_t i, num_app_mi = append_mi.size();
  append_mi_map.resize(num_app_mi);
  if (ref_mi.empty()) {
    combined_mi = append_mi;
    append_mi_map_ref = 0;
    for (i=0; i<num_app_mi; ++i)
      append_mi_map[i] = i;
  }
  else {
    combined_mi = ref_mi;
    append_mi_map_ref = combined_mi.size();
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = app_mi[i];
      size_t index = find_index(combined_mi, search_mi);
      if (index == _NPOS) { // search_mi does not yet exist in multi_index
	append_mi_map[i] = combined_mi.size();
	combined_mi.push_back(search_mi);
      }
      else
	append_mi_map[i] = index;
    }
  }
}
*/

/*  Append append_mi to ref_mi to create combined_mi using previously
    defined append_mi_map and append_mi_map_ref.  If necessary, update
    append_mi_map and append_mi_map_ref.
inline void append_multi_index(const UShort2DArray& ref_mi,
                               const UShort2DArray& append_mi,
			       SizetArray& append_mi_map,
			       size_t& append_mi_map_ref,
			       UShort2DArray& combined_mi)
{
  if (ref_mi.empty())
    combined_mi = append_mi; // assume append_mi_map{,_ref} are up to date
  else {
    combined_mi = ref_mi;
    size_t i, num_app_mi = append_mi.size(), num_mi = combined_mi.size();
    if (num_mi == append_mi_map_ref) { // current mi corresponds to ref
      for (i=0; i<num_app_mi; ++i)
	if (append_mi_map[i] >= append_mi_map_ref)
	  combined_mi.push_back(append_mi[i]);
    }
    else if (num_mi > append_mi_map_ref) { // mi has grown since ref taken
      for (i=0; i<num_app_mi; ++i)
	if (append_mi_map[i] >= append_mi_map_ref) { // previously appended
	  const UShortArray& search_mi = append_mi[i];
	  // search from reference pt forward
	  UShort2DArray::iterator it, it_start = combined_mi.begin();
	  std::advance(it_start, append_mi_map_ref);
	  it = std::find(it_start, combined_mi.end(), search_mi);
	  if (it == combined_mi.end()) { // still an append: update map, append
	    append_mi_map[i] = combined_mi.size();
	    combined_mi.push_back(append_mi[i]);
	  }
	  else // no longer an append: only update map
	    append_mi_map[i] = append_mi_map_ref + std::distance(it_start, it);
	}
      append_mi_map_ref = num_mi; // reference point now updated
    }
    else { // combined_mi is not allowed to shrink since ref taken
      PCerr << "Error: ref_mi inconsistent with reference size in "
	    << "append_multi_index()." << std::endl;
      abort_handler(-1);
    }
  }
}


// assess whether pareto array is dominated by combined_pareto
bool SharedOrthogPolyApproxData::
assess_dominance(const UShort2DArray& pareto,
		 const UShort2DArray& combined_pareto)
{
  bool new_dominated = true, i_dominated, j_dominated;
  size_t i, j, num_p = pareto.size(),
    num_combined_p = combined_pareto.size();
  for (i=0; i<num_p; ++i) {
    const UShortArray& pareto_i = pareto[i];
    i_dominated = false;
    for (j=0; j<num_combined_p; ++j) {
      assess_dominance(pareto_i, combined_pareto[j], i_dominated, j_dominated);
      if (i_dominated) break;
    }
    if (!i_dominated) {
      new_dominated = false;
#ifdef DEBUG
      PCout << "Nondominated new pareto member =\n" << pareto_i;
#else
      break;
#endif // DEBUG
    }
  }
  return new_dominated;
}
*/


/** Pareto dominance (typical definition): multi_index "a" weakly dominates
    multi_index "b" iff a_i >= b_i for all i and a_i > b_i for at least 
    one dimension.  Here we add the additional distinction of a challenger
    versus an incumbent: tie goes to the incumbent (the challenger is not
    added redundantly to the Pareto set). */
inline void assess_dominance(const UShortArray& new_order,
			     const UShortArray& ref_order,
			     bool& new_dominated, bool& ref_dominated)
{
  // can't use std::vector::operator< (used for component-wise sorting)
  size_t i, n = new_order.size();
  bool equal = true, ref_dominated_temp = true;
  new_dominated = true;
  for (i=0; i<n; ++i)
    if (new_order[i] > ref_order[i])
      { equal = false; new_dominated = false; }
    else if (ref_order[i] > new_order[i])
      { equal = false; ref_dominated_temp = false; }
  // asymmetric logic since incumbent wins a tie
  ref_dominated = (!equal && ref_dominated_temp);
}


/** "Strong" Pareto dominance: multi_index "a" strongly dominates multi_index
    "b" iff a_i > b_i for all i.  This case needs no notion of challenger
    versus incumbent. */
inline void assess_strong_dominance(const UShortArray& order_a,
				    const UShortArray& order_b,
				    bool& a_dominated, bool& b_dominated)
{
  // can't use std::vector::operator< (used for component-wise sorting)
  size_t i, n = order_a.size();
  a_dominated = b_dominated = true;
  for (i=0; i<n; ++i)
    if (order_a[i] == order_b[i])
      { a_dominated = b_dominated = false; break; }
    else if (order_a[i] > order_b[i])
      a_dominated = false;
    else // order_b[i] > order_a[i]
      b_dominated = false;
}


/** Pareto dominance (typical definition): return true if new_order dominates
    ref_order (using typical dominance definition). */
inline bool new_dominates_reference(const UShortArray& new_order,
				    const UShortArray& ref_order)
{
  size_t i, n = new_order.size();
  bool equal = true, ref_dominated = true;
  for (i=0; i<n; ++i)
    if (new_order[i] > ref_order[i])
      equal = false;
    else if (ref_order[i] > new_order[i])
      { equal = false; ref_dominated = false; }
  // asymmetric logic since incumbent wins a tie
  return (!equal && ref_dominated);
}


/** Pareto dominance (typical definition): return true if new_order dominates
    all entries in ref_orders (using typical dominance definition). */
inline bool new_dominates_reference(const UShortArray&   new_order,
				    const UShort2DArray& ref_orders)
{
  size_t i, num_orders = ref_orders.size();
  for (i=0; i<num_orders; ++i)
    if (!new_dominates_reference(new_order, ref_orders[i]))
      { return false; break; }
  return true;
}

} // namespace Pecos

#endif

/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:       DiscrepancyCalculator
//- Description: Utility for computing discrepancies between a truth model
//-              and an approximation.
//- Owner:       Mike Eldred
//- Checked by:
//- Version: $Id: DiscrepancyCalculator.hpp 7024 2010-10-16 01:24:42Z mseldre $

#ifndef DISCREPANCY_CALCULATOR_H
#define DISCREPANCY_CALCULATOR_H

#include "pecos_data_types.hpp"


namespace Pecos {

class SurrogateData; // forward declaration


/// Class for discrepancy calculations

/** The DiscrepancyCalculator class provides common functions for
    computing additive and multiplicative discrepancies. */

class DiscrepancyCalculator
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  DiscrepancyCalculator();
  /// destructor
  ~DiscrepancyCalculator();

  //
  //- Heading: Member functions
  //

  /// check for numerical issues with multiplicative discrepancy calculations
  static bool check_multiplicative(Real truth_fn, Real approx_fn,
				   short corr_order);
  /// check for numerical issues with multiplicative discrepancy calculations
  static bool check_multiplicative(const RealVector& truth_fns,
				   const RealVector& approx_fns,
				   short corr_order);

  /// compute additive 0th order correction between truth and approximate values
  static void compute_additive(Real truth_fn, Real approx_fn, Real& discrep_fn);
  /// compute additive 1st order correction between truth and
  /// approximate gradients
  static void compute_additive(const RealVector& truth_grad,
			       const RealVector& approx_grad,
			       RealVector& discrep_grad);
  /// compute additive 2nd order correction between truth and
  /// approximate Hessians
  static void compute_additive(const RealSymMatrix& truth_hess,
			       const RealSymMatrix& approx_hess,
			       RealSymMatrix& discrep_hess);
  /// compute additive corrections between truth and approximate responses
  /// for all requested orders
  static void compute_additive(Real truth_fn, const RealVector& truth_grad,
			       const RealSymMatrix& truth_hess,
			       Real approx_fn, const RealVector& approx_grad,
			       const RealSymMatrix& approx_hess,
			       Real& discrep_fn, RealVector& discrep_grad,
			       RealSymMatrix& discrep_hess,
			       short data_bits = 7);

  /// compute multiplicative 0th order correction between truth and
  /// approximate values
  static void compute_multiplicative(Real truth_fn, Real approx_fn,
				     Real& discrep_fn);
  /// compute multiplicative 1st order correction between truth and
  /// approximate gradients
  static void compute_multiplicative(Real truth_fn,
				     const RealVector& truth_grad,
				     Real approx_fn,
				     const RealVector& approx_grad,
				     RealVector& discrep_grad);
  /// compute multiplicative 2nd order correction between truth and
  /// approximate Hessians
  static void compute_multiplicative(Real truth_fn,
				     const RealVector& truth_grad,
				     const RealSymMatrix& truth_hess,
				     Real approx_fn,
				     const RealVector& approx_grad,
				     const RealSymMatrix& approx_hess,
				     RealSymMatrix& discrep_hess);
  /// compute multiplicative corrections between truth and approximate
  /// responses for all required orders
  static void compute_multiplicative(Real truth_fn,
				     const RealVector& truth_grad,
				     const RealSymMatrix& truth_hess,
				     Real approx_fn,
				     const RealVector& approx_grad,
				     const RealSymMatrix& approx_hess,
				     Real& discrep_fn, RealVector& discrep_grad,
				     RealSymMatrix& discrep_hess,
				     short data_bits = 7);

  /// compute discrepancies for arrays of SurrogateDataResp
  static void compute(const SDRArray& hf_sdr_array,
		      const SDRArray& lf_sdr_array,
		      SDRArray& delta_sdr_array, short combine_type);
  /// compute discrepancies between two model keys in surr_data and store
  /// results in the discrepancy key
  static void compute(SurrogateData& surr_data, const UShortArray& delta_key,
		      short combine_type);

  /// define a model key including data group, model form, and resolution
  /// level indices
  static void form_key(unsigned short group, unsigned short form,
		       unsigned short lev, UShortArray& key);
  /// decrement an incoming model key to correspond to the next lower
  /// resolution or fidelity within a model sequence
  static bool decrement_key(UShortArray& key);
  /// test whether key is an aggregated (e.g., discrepancy) key
  static bool aggregated_key(const UShortArray& key);
  /// aggregate two model keys to indicate a data combination
  /// (e.g., a discrepancy)
  static void aggregate_keys(const UShortArray& key1, const UShortArray& key2,
			     UShortArray& aggregate_key);
  /// extract the constituent keys from an aggregated key
  static void extract_keys(const UShortArray& aggregate_key, UShortArray& key1,
			   UShortArray& key2);

  /*
  /// function for applying additive correction to an approximate response
  void apply_additive(const Variables& vars, Response& approx_response);
  /// function for applying multiplicative correction to an approximate response
  void apply_multiplicative(const Variables& vars,
			    Response& approx_response);

  /// update correctionType
  void correction_type(short corr_type);
  /// return correctionType
  short correction_type() const;
  /// update correctionOrder
  void correction_order(short order);
  /// return correctionOrder
  short correction_order() const;
  /// update dataOrder
  void data_order(short order);
  /// return dataOrder
  short data_order() const;
  */

private:

  //
  //- Heading: Convenience functions
  //

  //
  //- Heading: Data
  //

  /*
  /// approximation correction approach to be used: NO_CORRECTION,
  /// ADDITIVE_CORRECTION, MULTIPLICATIVE_CORRECTION, or COMBINED_CORRECTION.
  short correctionType;
  /// approximation correction order to be used: 0, 1, or 2
  short correctionOrder;
  /// order of correction data in 3-bit format: overlay of 1 (value),
  /// 2 (gradient), and 4 (Hessian)
  short dataOrder;
  */
};


inline DiscrepancyCalculator::DiscrepancyCalculator()
  //correctionType(NO_CORRECTION), correctionOrder(0), dataOrder(1)
{ }


inline DiscrepancyCalculator::~DiscrepancyCalculator()
{ }


inline void DiscrepancyCalculator::
compute_additive(Real truth_fn, const RealVector& truth_grad,
		 const RealSymMatrix& truth_hess,
		 Real approx_fn, const RealVector& approx_grad,
		 const RealSymMatrix& approx_hess,
		 Real& discrep_fn, RealVector& discrep_grad,
		 RealSymMatrix& discrep_hess, short data_bits)
{
  if (data_bits & 1)
    compute_additive(truth_fn, approx_fn, discrep_fn);
  if (data_bits & 2)
    compute_additive(truth_grad, approx_grad, discrep_grad);
  if (data_bits & 4)
    compute_additive(truth_hess, approx_hess, discrep_hess);
}


inline void DiscrepancyCalculator::
compute_multiplicative(Real truth_fn, const RealVector& truth_grad,
		       const RealSymMatrix& truth_hess,
		       Real approx_fn, const RealVector& approx_grad,
		       const RealSymMatrix& approx_hess,
		       Real& discrep_fn, RealVector& discrep_grad,
		       RealSymMatrix& discrep_hess, short data_bits)
{
  if (data_bits & 1)
    compute_multiplicative(truth_fn, approx_fn, discrep_fn);
  if (data_bits & 2)
    compute_multiplicative(truth_fn, truth_grad, approx_fn, approx_grad,
			   discrep_grad);
  if (data_bits & 4)
    compute_multiplicative(truth_fn, truth_grad, truth_hess, approx_fn,
			   approx_grad, approx_hess, discrep_hess);
}


inline void DiscrepancyCalculator::
form_key(unsigned short group, unsigned short form, unsigned short lev,
	 UShortArray& key)
{ key.resize(3); key[0] = group; key[1] = form; key[2] = lev; }


/*
inline void DiscrepancyCalculator::
form_key(unsigned short lev, unsigned short form, bool multilevel,
	 UShortArray& model_key)
{
  if (multilevel) // model form is fixed @ HF; lev enumerates the levels
    { model_key.resize(2); model_key[0] = form; model_key[1] = lev; }
  else            // lev enumerates the model forms; levels are ignored
    { model_key.resize(1); model_key[0] = lev; } // mi_key[1] = _NPOS;
}


inline void DiscrepancyCalculator::
modified_lf_key(const UShortArray& hf_key, UShortArray& lf_key)
{
  if (hf_key.back() > 0) {
    // decrement trailing index
    lf_key = hf_key; --lf_key.back();
    // append the HF key in order to tag a particular (discrepancy) pairing
    lf_key.insert(lf_key.end(), hf_key.begin(), hf_key.end());
  }
  else
    lf_key.clear();
}
*/


inline bool DiscrepancyCalculator::decrement_key(UShortArray& key)
{
  // decrement the active index, if present, to create a key within the same
  // group id but with the next lower resolution in the sequence

  if (key.size() != 3) { // don't allow aggregated keys
    PCerr << "Error: wrong size for {group,form,lev} format in Discrepancy"
	  << "Calculator::decrement_key()" << std::endl;
    abort_handler(-1);    
  }

  // Logic is fragile in that it fails if a fixed model index (index that is
  // not part of the sequence) is assigned a value other than 0 or USHRT_MAX
  // > precedence given to lev for this reason, as form is more likely to have
  //   a non-zero/inf fixed value (see NonDExpansion::configure_sequence())
  // > more robust approach would be to pass in a multilev boolean
  unsigned short &form = key[1], &lev = key[2];
  if      (lev  && lev  != USHRT_MAX)
    { --lev;  return true; }
  else if (form && form != USHRT_MAX)
    { --form; return true; }
  //else no op (already at coarsest resolution / lowest fidelity)
  return false;

  // Old logic for {form} | {form,lev} format was simply --key.back();
}


inline bool DiscrepancyCalculator::aggregated_key(const UShortArray& key)
{
  size_t len = key.size();
  switch (len) {
  case 0: case 3: // no key or single model (not aggregated)
    return false; break;
  case 1: case 2:
    PCerr << "Error: invalid key size for {group,{form,lev}} format in "
	  << "DiscrepancyCalculator::aggregated_key()" << std::endl;
    abort_handler(-1);
    return false; break;
  default: {
    if (len-1 % 2) { // expect {form,lev} pairs following group
      PCerr << "Error: invalid key size for {group,{form,lev}} format in "
	    << "DiscrepancyCalculator::aggregated_key()" << std::endl;
      abort_handler(-1);
    }
    return true;  break;
  }
  }
}


inline void DiscrepancyCalculator::
aggregate_keys(const UShortArray& key1, const UShortArray& key2,
	       UShortArray& aggregate_key)
{
  // extract and verify consistency in group number
  unsigned short group = USHRT_MAX;
  bool empty1 = key1.empty(), empty2 = key2.empty();
  if (!empty1 && !empty2) {
    group = key1.front();
    if (group != key2.front()) {
      PCerr << "Error: mismatch in group ids in DiscrepancyCalculator::"
	    << "aggregate_keys()" << std::endl;
      abort_handler(-1);
    }
  }
  else if (!empty1) group = key1.front();
  else if (!empty2) group = key2.front();
  else {
    PCerr << "Error: neither key defined in DiscrepancyCalculator::"
	  << "aggregate_keys(key1, key2)" << std::endl;
    abort_handler(-1);    
  }
  
  // form aggregate of group + HF form/lev + LF form/lev
  aggregate_key.resize(1);  aggregate_key[0] = group;
  if (!empty1)
    aggregate_key.insert(aggregate_key.end(), key1.begin()+1, key1.end());
  if (!empty2)
    aggregate_key.insert(aggregate_key.end(), key2.begin()+1, key2.end());
}


inline void DiscrepancyCalculator::
extract_keys(const UShortArray& aggregate_key, UShortArray& key1,
	     UShortArray& key2)
{
  if (aggregate_key.empty())
    { key1.clear(); key2.clear(); return; }
  
  // extract one or two aggregated keys
  unsigned short group = aggregate_key.front(),
    len = aggregate_key.size() - 1, num_keys = len / 2;
  switch (num_keys) {
  case 1:
    key1 = aggregate_key;   key2.clear();  break;
  case 2: { // normal case
    UShortArray::const_iterator start1 = aggregate_key.begin() + 1,
      end1 = start1 + 2, end2 = end1 + 2;
    key1.assign(1, group);  key1.insert(key1.end(), start1, end1);
    key2.assign(1, group);  key2.insert(key2.end(),   end1, end2);
    break;
  }
  default:
    PCerr << "Error: bad aggregate key size in DiscrepancyCalculator::"
	  << "extract_keys()" << std::endl;
    abort_handler(-1);
    break;
  }
}

} // namespace Pecos

#endif

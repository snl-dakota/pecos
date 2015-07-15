/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 HistogramBinRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef HISTOGRAM_BIN_RANDOM_VARIABLE_HPP
#define HISTOGRAM_BIN_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for gumbel random variables.

/** Manages the binPairs mapping. */

class HistogramBinRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  HistogramBinRandomVariable();
  /// alternate constructor
  HistogramBinRandomVariable(const RealRealMap& bin_prs);
  /// destructor
  ~HistogramBinRandomVariable();

  //
  //- Heading: Virtual function redefinitions
  //

  Real cdf(Real x) const;
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  Real pdf_hessian(Real x) const;

  RealRealPair moments() const;
  RealRealPair bounds() const;

  //
  //- Heading: Member functions
  //

  void update(const RealRealMap& bin_prs);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, const RealRealMap& bin_prs);
  static Real cdf(Real x, const RealRealMap& bin_prs);
  static Real inverse_cdf(Real cdf, const RealRealMap& bin_prs);

  static Real pdf(Real x, const RealVector& bin_prs);
  static Real cdf(Real x, const RealVector& bin_prs);
  static Real inverse_cdf(Real cdf, const RealVector& bin_prs);

  static void moments_from_params(const RealRealMap& bin_prs,
				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Data
  //

  /// pairings of abscissas to bin counts
  RealRealMap binPairs;
};


inline HistogramBinRandomVariable::HistogramBinRandomVariable():
  RandomVariable(BaseConstructor())
{ ranVarType = HISTOGRAM_BIN; }


inline HistogramBinRandomVariable::
HistogramBinRandomVariable(const RealRealMap& bin_prs):
  RandomVariable(BaseConstructor()), binPairs(bin_prs)
{ ranVarType = HISTOGRAM_BIN; }


inline HistogramBinRandomVariable::~HistogramBinRandomVariable()
{ }


inline Real HistogramBinRandomVariable::cdf(Real x, const RealRealMap& bin_prs)
{
  size_t i, num_bins = bin_prs.size() - 1;
  RRMCIter cit = bin_prs.begin();
  if (x <= cit->first)
    return 0.;
  else if (x >= (--bin_prs.end())->first)
    return 1.;
  else {
    Real p_cdf = 0., count, upr, lwr;
    for (i=0; i<num_bins; ++i) {
      lwr = cit->first; count = cit->second; ++cit;
      upr = cit->first;
      if (x < upr) {
	p_cdf += count * (x - lwr) / (upr - lwr);
	break;
      }
      else
	p_cdf += count;
    }
    return p_cdf;
  }
}


inline Real HistogramBinRandomVariable::cdf(Real x) const
{ return cdf(x, binPairs); }


inline Real HistogramBinRandomVariable::ccdf(Real x) const
{
  size_t i, num_bins = binPairs.size() - 1;
  RRMCIter cit = binPairs.begin();
  if (x <= cit->first)
    return 1.;
  else if (x >= (--binPairs.end())->first)
    return 0.;
  else {
    Real p_ccdf = 1., count, upr, lwr;
    for (i=0; i<num_bins; ++i) {
      lwr = cit->first; count = cit->second; ++cit;
      upr = cit->first;
      if (x < upr) {
	p_ccdf -= count * (x - lwr) / (upr - lwr);
	break;
      }
      else
	p_ccdf -= count;
    }
    return p_ccdf;
  }
}


inline Real HistogramBinRandomVariable::
inverse_cdf(Real p_cdf, const RealRealMap& bin_prs)
{
  size_t i, num_bins = bin_prs.size() - 1;
  RRMCIter cit = bin_prs.begin();
  if (p_cdf <= 0.)
    return cit->first;               // lower bound abscissa
  else if (p_cdf >= 1.)
    return (--bin_prs.end())->first; // upper bound abscissa
  else {
    Real upr_cdf = 0., count, upr, lwr;
    for (i=0; i<num_bins; ++i) {
      lwr = cit->first; count = cit->second; ++cit;
      upr = cit->first;
      upr_cdf += count;
      if (p_cdf < upr_cdf)
	return upr - (upr_cdf - p_cdf) / count * (upr - lwr);
    }
    // If still no return due to numerical roundoff, return upper bound
    return (--bin_prs.end())->first; // upper bound abscissa
  }
}


inline Real HistogramBinRandomVariable::inverse_cdf(Real p_cdf) const
{ return inverse_cdf(p_cdf, binPairs); }


inline Real HistogramBinRandomVariable::inverse_ccdf(Real p_ccdf) const
{
  size_t i, num_bins = binPairs.size() - 1;
  RRMCIter cit = binPairs.begin();
  if (p_ccdf >= 1.)
    return cit->first;                // lower bound abscissa
  else if (p_ccdf <= 0.)
    return (--binPairs.end())->first; // upper bound abscissa
  else {
    Real upr_ccdf = 1., count, upr, lwr;
    for (i=0; i<num_bins; ++i) {
      lwr = cit->first; count = cit->second; ++cit;
      upr = cit->first;
      upr_ccdf -= count;
      if (p_ccdf > upr_ccdf)
	return upr - (p_ccdf - upr_ccdf) / count * (upr - lwr);
    }
    // If still no return due to numerical roundoff, return upper bound
    return (--binPairs.end())->first; // upper bound abscissa
  }
}


inline Real HistogramBinRandomVariable::pdf(Real x, const RealRealMap& bin_prs)
{
  // The PDF is discontinuous at the bin bounds.  To resolve this, this
  // function employs an (arbitrary) convention: closed/inclusive lower bound
  // and open/exclusive upper bound.
  size_t i, num_bins = bin_prs.size() - 1;
  RRMCIter cit = bin_prs.begin();
  if (x < cit->first || x >= (--bin_prs.end())->first) // outside support
    return 0.;
  else {
    Real lwr, upr, count;
    for (i=0; i<num_bins; ++i) {
      lwr = cit->first; count = cit->second; ++cit;
      upr = cit->first;
      if (x < upr) // return density
	return count / (upr - lwr);
    }
    return 0.;
  }
}


inline Real HistogramBinRandomVariable::pdf(Real x) const
{ return pdf(x, binPairs); }


inline Real HistogramBinRandomVariable::pdf_gradient(Real x) const
{ return 0.; }


inline Real HistogramBinRandomVariable::pdf_hessian(Real x) const
{ return 0.; }


inline RealRealPair HistogramBinRandomVariable::moments() const
{
  Real mean, std_dev;
  moments_from_params(binPairs, mean, std_dev);
  return RealRealPair(mean, std_dev);
}


inline RealRealPair HistogramBinRandomVariable::bounds() const
{ return RealRealPair(binPairs.begin()->first, (--binPairs.end())->first); }


inline void HistogramBinRandomVariable::update(const RealRealMap& bin_prs)
{ binPairs = bin_prs; }


inline Real HistogramBinRandomVariable::pdf(Real x, const RealVector& bin_prs)
{
  // Need this case to be fast for usage in NumericGenOrthogPolynomial...
  /* Cleaner, but induces unnecessary copy overhead:
  RealRealMap bin_prs_rrm;
  copy_data(bin_prs_rv, bin_prs_rrm);
  return HistogramBinRandomVariable::pdf(x, bin_prs_rrm);
  */

  size_t num_bins = bin_prs.length() / 2 - 1;
  if (x < bin_prs[0] || x >= bin_prs[2*num_bins])
    return 0.;
  else {
    Real upr;
    for (int i=0; i<num_bins; ++i) {
      upr = bin_prs[2*i+2];
      if (x < upr) // return density = count / (upr - lwr);
	return bin_prs[2*i+1] / (upr - bin_prs[2*i]);
    }
    return 0.;
  }
}


inline Real HistogramBinRandomVariable::cdf(Real x, const RealVector& bin_prs)
{
  /* Cleaner, but induces unnecessary copy overhead:
  RealRealMap bin_prs_rrm;
  copy_data(bin_prs_rv, bin_prs_rrm);
  return HistogramBinRandomVariable::cdf(x, bin_prs_rrm);
  */

  size_t num_bins = bin_prs.length() / 2 - 1;
  if (x <= bin_prs[0])
    return 0.;
  else if (x >= bin_prs[2*num_bins])
    return 1.;
  else {
    Real p_cdf = 0., count, upr, lwr;
    for (int i=0; i<num_bins; ++i) {
      count = bin_prs[2*i+1]; upr = bin_prs[2*i+2];
      if (x < upr) {
        lwr = bin_prs[2*i];
	p_cdf += count * (x - lwr) / (upr - lwr);
	break;
      }
      else
	p_cdf += count;
    }
    return p_cdf;
  }
}


inline Real HistogramBinRandomVariable::
inverse_cdf(Real p_cdf, const RealVector& bin_prs)
{
  /* Cleaner, but induces unnecessary copy overhead:
  RealRealMap bin_prs_rrm;
  copy_data(bin_prs_rv, bin_prs_rrm);
  return HistogramBinRandomVariable::inverse_cdf(p_cdf, bin_prs_rrm);
  */

  size_t num_bins = bin_prs.length() / 2 - 1;
  if (p_cdf <= 0.)
    return bin_prs[0];          // lower bound abscissa
  else if (p_cdf >= 1.)
    return bin_prs[2*num_bins]; // upper bound abscissa
  else {
    Real upr_cdf = 0., count, upr, lwr;
    for (int i=0; i<num_bins; ++i) {
      count = bin_prs[2*i+1];
      upr_cdf += count;
      if (p_cdf < upr_cdf) {
        upr = bin_prs[2*i+2]; lwr = bin_prs[2*i];
	return upr - (upr_cdf - p_cdf) / count * (upr - lwr);
      }
    }
    // If still no return due to numerical roundoff, return upper bound
    return bin_prs[2*num_bins]; // upper bound abscissa
  }
}


inline void HistogramBinRandomVariable::
moments_from_params(const RealRealMap& bin_prs, Real& mean, Real& std_dev)
{
  // in bin case, (x,y) and (x,c) are not equivalent since bins have non-zero
  // width -> assume (x,c) count/frequency-based (already converted from (x,y)
  // skyline/density-based) with normalization (counts sum to 1.)
  mean = std_dev = 0.;
  size_t i, num_bins = bin_prs.size() - 1;
  //RRMCIter it_lwr = bin_prs.begin();
  //RRMCIter it_upr = ++bin_prs.begin();
  //RRMCIter it_end = --(bin_prs.end());
  //for (; it_lwr != it_end; ++it_lwr, ++it_upr) {
  //  lwr = it_lwr->first;
  //  count = it_lwr->second;
  //  upr = it_upr->first;
  //  ...
  //}
  Real lwr, count, upr;//, clu;
  RRMCIter cit = bin_prs.begin();
  for (i=0; i<num_bins; ++i) {
    lwr = cit->first; count = cit->second; ++cit;
    upr = cit->first;
    //clu = count * (lwr + upr);
    mean    += count * (lwr + upr); // clu
    std_dev += count * (upr*upr + upr*lwr + lwr*lwr); // upr*clu + count*lwr*lwr
  }
  mean   /= 2.;
  std_dev = std::sqrt(std_dev/3. - mean*mean);
}

} // namespace Pecos

#endif

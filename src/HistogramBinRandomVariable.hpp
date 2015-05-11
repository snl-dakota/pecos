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

/** Manages alpha and beta parameters. */

class HistogramBinRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  HistogramBinRandomVariable();
  /// alternate constructor
  HistogramBinRandomVariable(const RealRealMap& hist_bin_prs);
  /// destructor
  ~HistogramBinRandomVariable();

  //
  //- Heading: Member functions
  //

  Real cdf(Real x) const;
  Real cdf_inverse(Real p) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  Real pdf_hessian(Real x) const;

protected:

  //
  //- Heading: Data
  //

  /// alpha parameter of gumbel random variable
  RealRealMap binPairs;
};


inline HistogramBinRandomVariable::HistogramBinRandomVariable():
  RandomVariable(BaseConstructor())
{ }


inline HistogramBinRandomVariable::
HistogramBinRandomVariable(const RealRealMap& hist_bin_prs):
  RandomVariable(BaseConstructor()), binPairs(hist_bin_prs)
{ }


inline HistogramBinRandomVariable::~HistogramBinRandomVariable()
{ }


inline Real HistogramBinRandomVariable::cdf(Real x) const
{ return histogram_bin_cdf(x, binPairs); }


inline Real HistogramBinRandomVariable::cdf_inverse(Real p) const
{ return histogram_bin_cdf_inverse(p, binPairs); }


inline Real HistogramBinRandomVariable::pdf(Real x) const
{ return histogram_bin_pdf(x, binPairs); }


inline Real HistogramBinRandomVariable::pdf_gradient(Real x) const
{ return 0.; }


inline Real HistogramBinRandomVariable::pdf_hessian(Real x) const
{ return 0.; }

} // namespace Pecos

#endif

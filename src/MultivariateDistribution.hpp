/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef MULTIVARIATE_DISTRIBUTION_HPP
#define MULTIVARIATE_DISTRIBUTION_HPP

#include "pecos_data_types.hpp"

namespace Pecos {


/// Base class for all multivariate distributions

/** The base class for multivariate distributions, including joint
    distributions and distributions based on marginals plus correlation. */

class MultivariateDistribution
{
public:

  /// default constructor
  MultivariateDistribution();
  /// standard constructor for envelope
  MultivariateDistribution(short mv_dist_type);
  /// copy constructor
  MultivariateDistribution(const MultivariateDistribution& mv_dist);

  /// destructor
  virtual ~MultivariateDistribution();

  /// assignment operator
  MultivariateDistribution operator=(const MultivariateDistribution& mv_dist);

  //
  //- Heading: Virtual functions
  //

  /// return the multivariate PDF value for the random variables
  virtual Real pdf(const RealVector& pt) const = 0;
  /// return the multivariate log PDF value for the random variables
  virtual Real log_pdf(const RealVector& pt) const = 0;

  //
  //- Heading: Member functions
  //

  // perform a deep copy of incoming mv_dist
  //void copy(const MultivariateDistribution& mv_dist);

  /// return correlationFlag
  bool correlation() const;

  /// function to check modelRep (does this envelope contain a letter)
  bool is_null() const;

protected:

  //
  //- Heading: Constructors
  //

  /// constructor initializes the base class part of letter classes
  /// (BaseConstructor overloading avoids infinite recursion in the
  /// derived class constructors - Coplien, p. 139)
  MultivariateDistribution(BaseConstructor);
  
  //
  //- Heading: Data members
  //

  /// flag for indicating if correlation exists among the random variables
  bool correlationFlag;

private:

  //
  //- Heading: Member functions
  //

  /// Used only by the standard envelope constructor to initialize
  /// mvDistRep to the appropriate derived type.
  MultivariateDistribution* get_distribution(short mv_dist_type);

  //
  //- Heading: Data members
  //

  /// pointer to the letter (initialized only for the envelope)
  MultivariateDistribution* mvDistRep;
  /// number of objects sharing mvDistRep
  int referenceCount;
};


inline bool MultivariateDistribution::correlation() const
{ return (mvDistRep) ? mvDistRep->correlationFlag : correlationFlag; }


inline bool MultivariateDistribution::is_null() const
{ return (mvDistRep) ? false : true; }

} // namespace Pecos

#endif

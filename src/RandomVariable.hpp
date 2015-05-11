/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 RandomVariable
//- Description: 
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef RANDOM_VARIABLE_HPP
#define RANDOM_VARIABLE_HPP

#include "pecos_data_types.hpp"
#include "pecos_stat_util.hpp"

namespace Pecos {


/// base class for random variable classes

/** This class enables cdf(), cdf_inverse(), pdf(), pdf_gradient(),
    and pdf_hessian() from contained distribution parameters. */

// ExtremeValueBase: Gumbel, Frechet, Weibull
//   > alpha, beta, no bounds
// EmpiricalBase: histogram bin, KDE, ...

class RandomVariable
{
public:

  //
  //- Heading: Constructors, destructor, and operator=
  //

  /// default constructor
  RandomVariable();
  /// standard constructor for envelope
  RandomVariable(short ran_var_type);
  /// copy constructor
  RandomVariable(const RandomVariable& ran_var);

  /// destructor
  virtual ~RandomVariable();

  /// assignment operator
  RandomVariable operator=(const RandomVariable& ran_var);

  //
  //- Heading: Virtual functions
  //

  /// return the cumulative distribution function value of the random
  /// variable at x
  virtual Real cdf(Real x) const;
  /// return the x value corresponding to a cumulative probability
  virtual Real cdf_inverse(Real p) const;

  /// return the probability density function value of the random
  /// variable at x
  virtual Real pdf(Real x) const;
  /// return the gradient of the probability density function value of
  /// the random variable at x
  virtual Real pdf_gradient(Real x) const;
  /// return the hessian of the probability density function value of
  /// the random variable at x
  virtual Real pdf_hessian(Real x) const;

  //
  //- Heading: Member functions
  //

  /// return the complementary cumulative distribution function value
  /// of the random variable at x
  Real ccdf(Real x) const;
  /// return the x value corresponding to a complementary cumulative probability
  Real ccdf_inverse(Real p) const;

  /// set ranVarType
  void type(short ran_var_type);
  /// get ranVarType
  short type() const;

  /// returns ranVarRep for access to derived class member functions
  /// that are not mapped to the base level
  RandomVariable* random_variable_rep() const;

protected:

  //
  //- Heading: Constructors
  //

  /// constructor initializes the base class part of letter classes
  /// (BaseConstructor overloading avoids infinite recursion in the
  /// derived class constructors - Coplien, p. 139)
  RandomVariable(BaseConstructor);

  //
  //- Heading: Member functions
  //

  //
  //- Heading: Data
  //

  /// enumeration value indicating type of random variable
  short ranVarType;

private:

  //
  //- Heading: Member functions
  //

  /// Used only by the standard envelope constructor to initialize
  /// ranVarRep to the appropriate derived type.
  RandomVariable* get_random_variable(short ran_var_type);

  //
  //- Heading: Data members
  //

  /// pointer to the letter (initialized only for the envelope)
  RandomVariable* ranVarRep;
  /// number of objects sharing ranVarRep
  int referenceCount;
};


inline Real RandomVariable::ccdf(Real x) const
{ return (ranVarRep) ? 1. - ranVarRep->cdf(x) : 1. - cdf(x); }


inline Real RandomVariable::ccdf_inverse(Real p) const
{ return (ranVarRep) ? ranVarRep->cdf_inverse(1. - p) : cdf_inverse(1. - p); }


inline void RandomVariable::type(short ran_var_type)
{
  if (ranVarRep) ranVarRep->ranVarType = ran_var_type;
  else           ranVarType = ran_var_type;
}


inline short RandomVariable::type() const
{ return (ranVarRep) ? ranVarRep->ranVarType : ranVarType; }


inline RandomVariable* RandomVariable::random_variable_rep() const
{ return ranVarRep; }

} // namespace Pecos

#endif

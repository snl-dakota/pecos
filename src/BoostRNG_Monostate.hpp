/*  ______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */
 
//- Class:       BoostRNG_Monostate
//- Description: Wrapper for various implementations of Random Number Generators
//- Owner:       Laura Swiler, Dave Gay, and Bill Bohnhoff 
//- Checked by:
//- Version: $Id$

#ifndef BOOST_RNG_MONOSTATE_HPP
#define BOOST_RNG_MONOSTATE_HPP

#include "pecos_data_types.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>


namespace Pecos {

typedef Real (*Rfunc)();

class BoostRNG_Monostate
{
private:

  //
  //- Heading: Data
  //

  static unsigned int rngSeed;
  static boost::mt19937 rnumGenerator;
  static boost::uniform_real<> uniDist;
  static boost::variate_generator<boost::mt19937&,boost::uniform_real<> > uniMT;

public:

  //
  //- Heading: Member functions
  //

  /// set randomSeed
  static void seed(unsigned int rng_seed);
  /// return randomSeed
  static unsigned int seed();

  static Real random_num1();

  //
  //- Heading: Data
  //

  /// Cached function pointers
  static Rfunc randomNum;
  static Rfunc randomNum2;
};


inline void BoostRNG_Monostate::seed(unsigned int rng_seed)
{ rngSeed = rng_seed; rnumGenerator.seed(rng_seed); }

inline unsigned int BoostRNG_Monostate::seed()
{ return rngSeed; }

inline Real BoostRNG_Monostate::random_num1()
{ return uniMT(); }

unsigned int BoostRNG_Monostate::rngSeed(41u); // 41 used in the Boost examples

boost::mt19937 BoostRNG_Monostate::rnumGenerator( BoostRNG_Monostate::seed() );

boost::uniform_real<> BoostRNG_Monostate::uniDist(0, 1);

boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
  BoostRNG_Monostate::uniMT(BoostRNG_Monostate::rnumGenerator,
                            BoostRNG_Monostate::uniDist);
 
Real (*BoostRNG_Monostate::randomNum)()  = BoostRNG_Monostate::random_num1;
Real (*BoostRNG_Monostate::randomNum2)() = BoostRNG_Monostate::random_num1;
} // namespace Pecos


#ifdef HAVE_LHS
#define rnum1 FC_FUNC(rnumlhs1,RNUMLHS1)
#define rnum2 FC_FUNC(rnumlhs2,RNUMLHS2)
#endif  // HAVE_LHS

extern "C" Pecos::Real rnum1(void), rnum2(void);

// inline leads to DIFFs in the Regression suite
/* inline */ Pecos::Real rnum1(void)
{
  //PCout << "running Boost MT" << "\n";
  return Pecos::BoostRNG_Monostate::randomNum();
}

/* inline */ Pecos::Real rnum2(void)
{
  // clone of rnum1
  //PCout << "running Boost MT" << "\n";
  return Pecos::BoostRNG_Monostate::randomNum2();
}

#endif  // BOOST_RNG_MONOSTATE_HPP


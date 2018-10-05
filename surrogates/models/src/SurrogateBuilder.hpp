/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef SURROGATE_BUILDER_HPP
#define SURROGATE_BUILDER_HPP

#include "Approximation.hpp"
#include "OptionsList.hpp"

namespace Surrogates {

/**
   \class SurrogateBuilder
   \brief Abstract base class for methods used to construct 
   * surrogates
*/
class SurrogateBuilder {
public:

  /// Function that is being emulated
  Function* targetFunction_;
  
  /// Constructor
  SurrogateBuilder(){};
  
  /// Destructor
  virtual ~SurrogateBuilder(){};

  /** \brief Build a surrogate to a set of specifications.
   *
   * \param[in] opts Spefications used to build the surrogate
   *
   */
  virtual void build(OptionsList &opts, Approximation &approx,
                     OptionsList &result) = 0;

  void set_target_function(Function &target_function){
    targetFunction_ = &target_function;};
};

}; // namespace surrogates

#endif //SURROGATE_BUILDER_HPP

/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include <VariableTransformation.hpp>

namespace Surrogates {

VariableTransformation::VariableTransformation(){}

VariableTransformation::~VariableTransformation(){}

void VariableTransformation::
set_variables(const std::shared_ptr<Surrogates::Variables> &vars){
  vars_ = vars;
}

int  VariableTransformation::num_vars(){
  return vars_->num_vars();
}

} // namespace Surrogates

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

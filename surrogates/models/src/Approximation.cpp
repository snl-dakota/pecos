#include "Approximation.hpp"

namespace Surrogates {

void Approximation::
set_variable_transformation(const std::shared_ptr<VariableTransformation> &var_transform){
  varTransform_ = var_transform;
}

std::shared_ptr<VariableTransformation> Approximation::
get_variable_transformation(){
  return varTransform_;
}

int Approximation::num_vars(){
  if (!varTransform_)
    throw(std::runtime_error("Variable transform has not been set"));
  return varTransform_->num_vars();
}

} // namespace Surrogates

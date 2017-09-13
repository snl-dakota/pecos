#include <Variables.hpp>

namespace Surrogates {

Variables::Variables() : numVars_(0){}

Variables::~Variables(){}

void Variables::set_options(const OptionsList &opts){
}

void Variables::get_options(OptionsList &opts){
}

void Variables::set_num_vars(int num_vars){
  numVars_ = num_vars;
}

int Variables::num_vars(){return numVars_;}

} // namespace Surrogates

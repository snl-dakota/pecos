#include <Function.hpp>

namespace Surrogates {

Function::Function() : numQOI_(0) {}

Function::~Function(){}

void Function::set_options(const OptionsList &opts) {
  opts_=opts;
  numQOI_ = opts_.get<int>("num_qoi");
}

void Function::get_options(OptionsList &opts) {
  opts=opts_;
}

void Function::gradient(const RealMatrix &samples, int qoi, RealMatrix &gradients) {
  throw(std::string("This Function type does not provide gradients"));
}

void Function::jacobian(const RealVector &sample, RealMatrix &jacobian) {
  throw(std::string("This Function type does not provide a jacobian"));
}

void Function::hessian(const RealMatrix &samples, int qoi, RealMatrixList &hessians) {
  throw(std::string("This Function type does not provide gradients"));
}

} // namespace Surrogates

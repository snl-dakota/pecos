#include <CppFunction.hpp>

namespace Surrogates {

CppFunction::CppFunction(){}

CppFunction::~CppFunction(){}

void CppFunction::
set_function(void (*function)(const Real* sample, Real* func_vals,
			      const OptionsList &opts)){
  targetFunction_ = function;
}

void CppFunction::value(const RealVector &sample, RealVector &values){
  resize_if_needed(values,numQOI_);
  targetFunction_(sample.values(),values.values(),opts_);
}

void CppFunction::value(const RealMatrix &samples, RealMatrix &values){
  resize_if_needed(values,samples.numCols(),numQOI_);
  RealVector val(numQOI_,false); //zeroOut = false
  for (int i=0; i<samples.numCols(); ++i){
    const RealVector sample = Teuchos::getCol(Teuchos::View,
					      const_cast<RealMatrix&>(samples),i);
    value(sample,val);
    for(int j = 0; j < numQOI_; j++)
      values(i,j) = val(j);
  }
}

} // namespace Surrogates

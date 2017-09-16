#include <vector>
#include "std_vector_example_type_defs.hpp"

namespace Pecos{

class MyClass {
private:
  RealArray result;
public:
  MyClass(){};
  virtual ~MyClass(){};
  virtual const RealArray& get(unsigned short n) {
    result.resize(n);
    return result;};
  void set(const Other::ShortArray &vec){
    std::cout << "[" ;
    for (size_t i=0; i< vec.size(); ++i)
      std::cout << vec[i] << " " ;
    std::cout << "]\n" ;
  }
};
}

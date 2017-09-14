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
};
}

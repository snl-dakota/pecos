#include <vector>
#include <iostream>
#include "std_vector_example_type_defs.hpp"

namespace Animal{

  void print_animal(Pecos::animal_enum n) {
    switch(n){
    case Pecos::Monkey :
      std::cout << "Monkey" << std::endl;
      break;
    case Pecos::Gorilla :
      std::cout << "Gorilla" << std::endl;
      break;
    default:
      std::cout << "Enum not found" << std::endl;
    }
  }

}


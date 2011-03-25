#include "HierarchicalLinearBasisTestClass.hpp"
#include <cppunit/ui/text/TestRunner.h>

int main(int argc, char** argv)
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(HierarchicalLinearBasisTestClass::suite());
  bool success = runner.run();

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
  
}

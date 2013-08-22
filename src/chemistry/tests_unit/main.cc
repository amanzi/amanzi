/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <UnitTest++.h>
#include <TestReporterStdout.h>

#include "chemistry_output.hh"
#include "chemistry_verbosity.hh"
#include "chemistry_containers.hh"

// create a global ChemistryOutput* pointer in the Amanzi::chemisry
// namespace that can be used by an other chemistry object
namespace Amanzi {
namespace AmanziChemistry {
extern ChemistryOutput* chem_out;
}  // end namespace AmanziChemistry
}  // end namespace Amanzi


int main(int argc, char* argv[]) {
  namespace ac = Amanzi::AmanziChemistry;

  ac::SetupDefaultChemistryOutput();

  int status = UnitTest::RunAllTests();

  delete ac::chem_out;

  return status;
}

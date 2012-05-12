/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <UnitTest++.h>
#include <TestReporterStdout.h>

#include "chemistry_output.hh"
#include "chemistry_verbosity.hh"
#include "chemistry_containers.hh"

// create a global ChemistryOutput* pointer in the amanzi::chemisry
// namespace that can be used by an other chemistry object
namespace amanzi {
namespace chemistry {
ChemistryOutput* chem_out = NULL;
}  // end namespace chemistry
}  // end namespace amanzi


int main(int argc, char* argv[]) {
  namespace ac = amanzi::chemistry;

  ac::SetupDefaultChemistryOutput();

  int status = UnitTest::RunAllTests();

  delete ac::chem_out;

  return status;
}

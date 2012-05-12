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
ChemistryOutput* chem_out;
}  // end namespace chemistry
}  // end namespace amanzi


int main(int argc, char* argv[]) {
  namespace ac = amanzi::chemistry;

  ac::OutputOptions output_options;
  output_options.use_stdout = true;
  output_options.file_name = "chemistry-unit-test-results.txt";
  output_options.verbosity_levels.push_back(ac::strings::kVerbosityVerbose);

  ac::chem_out = new ac::ChemistryOutput();
  ac::chem_out->Initialize(output_options);

  int status = UnitTest::RunAllTests();
  delete ac::chem_out;
  return status;
}

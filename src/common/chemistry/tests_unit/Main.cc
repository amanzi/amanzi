#include <UnitTest++.h>
#include <TestReporterStdout.h>

#include "VerboseObject_objs.hh"
#include "VerboseObject.hh"

#include "chemistry_verbosity.hh"

// create a global ChemistryOutput* pointer in the Amanzi::chemisry
// namespace that can be used by an other chemistry object
namespace Amanzi {
namespace AmanziChemistry {
extern VerboseObject* chem_out;
}  // end namespace AmanziChemistry
}  // end namespace Amanzi


int main(int argc, char* argv[]) {

  Teuchos::ParameterList plist;
  Amanzi::AmanziChemistry::chem_out = new Amanzi::VerboseObject("ChemistryPK", plist); 

  int status = UnitTest::RunAllTests();

  delete Amanzi::AmanziChemistry::chem_out;

  return status;
}

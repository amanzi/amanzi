#include <mpi.h>
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "global_verbosity.hh"

#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"

#include "MeshFactory.hh"
#include "Mesh_STK.hh"

#include "field.hh"
#include "field_composite_vector.hh"
#include "state.hh"
#include "constant_temperature.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;


struct TestState {
  Epetra_MpiComm *comm;
  Teuchos::RCP<Mesh> mesh;
  Teuchos::RCP<State> state;

  TestState() {
    comm = new Epetra_MpiComm(MPI_COMM_WORLD);
    MeshFactory mesh_fact(comm);
    mesh = mesh_fact(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
    state = Teuchos::rcp(new State());
    state->RegisterDomainMesh(mesh);
  }
  ~TestState() { delete comm; }
};


Teuchos::EVerbosityLevel ATS::VerbosityLevel::level_ = Teuchos::VERB_MEDIUM;

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  TestState test;

  Teuchos::ParameterList plist;
  Teuchos::updateParametersFromXmlFile(std::string("const_temp.xml"), &plist);

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector("soln"));
  plist.set("PK name", "constant temperature");
  Amanzi::Energy::ConstantTemperature ct(plist,soln);

  ct.setup(test.state.ptr());

  test.state->Setup();
  ct.initialize(test.state.ptr());
  test.state->Initialize();
  ct.commit_state(0., test.state);
}

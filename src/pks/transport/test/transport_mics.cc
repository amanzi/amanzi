/*
The transport component of the Amanzi code, serial unit tests.
License: BSD
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

// TPLs
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "State.hh"

#include "Transport_PK.hh"


/* **************************************************************** 
 * Test Init() procedure in the constructor.
 * ************************************************************* */
TEST(CONSTRUCTOR) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: read transport XML file" << std::endl;
#ifdef HAVE_MPI
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  std::string xmlFileName = "test/transport_mics.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);
 
  /* create an MSTK mesh framework */
  ParameterList region_list = plist->get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(Simple);

  MeshFactory factory(comm);
  factory.preference(pref);
  RCP<const Mesh> mesh = factory(0.0,0.0,0.0, 1.0,1.0,1.0, 1, 2, 1, gm); 
 
  //MeshAudit audit(mesh);
  //audit.Verify();

  /* create a simple state and populate it */
  std::vector<std::string> component_names;
  component_names.push_back("Component 0");
  component_names.push_back("Component 1");

  Teuchos::ParameterList state_list = plist->sublist("State");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
  S->set_time(0.0);
  S->set_intermediate_time(0.0);

  Transport_PK TPK(plist, S, "Transport", component_names);
  TPK.Setup(S.ptr());
  TPK.CreateDefaultState(mesh, 2);

  /* initialize a transport process kernel from a transport state */
  TPK.Initialize(S.ptr());

  double cfl = TPK.cfl();
  CHECK(0 < cfl && cfl <= 1.0);
 
  delete comm;
}




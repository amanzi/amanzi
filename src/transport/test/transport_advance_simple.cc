/*
  The transport component of the Amanzi code, serial unit tests.
  License: BSD
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "Mesh.hh"

#include "State.hh"
#include "Transport_PK.hh"


TEST(ADVANCE_WITH_SIMPLE) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: advance using simple mesh" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  string xmlFileName = "test/transport_advance_simple.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  /* create an SIMPLE mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, (Epetra_MpiComm *)comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(Simple);

  MeshFactory meshfactory(comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory(0.0,0.0,0.0, 1.0,1.0,1.0, 20, 1, 1, gm); 

  /* create a simple state and populate it */
  Amanzi::VerboseObject::hide_line_prefix = false;

  std::vector<std::string> component_names;
  component_names.push_back("Component 0");

  RCP<State> S = rcp(new State());
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Transport_PK TPK(plist, S, component_names);
  TPK.CreateDefaultState(mesh, 2);

  /* modify the default state for the problem at hand */
  std::string passwd("state"); 
  Teuchos::RCP<Epetra_MultiVector> 
      flux = S->GetFieldData("darcy_flux", passwd)->ViewComponent("face", false);

  AmanziGeometry::Point velocity(1.0, 0.0, 0.0);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    (*flux)[0][f] = velocity * normal;
  }

  /* initialize a transport process kernel */
  TPK.InitPK();
  TPK.PrintStatistics();

  /* advance the state */
  double T = 0.0;
  Teuchos::RCP<Epetra_MultiVector> 
      tcc = S->GetFieldData("total_component_concentration", passwd)->ViewComponent("cell", false);

  int iter = 0;
  while (T < 1.0) {
    double dT = TPK.CalculateTransportDt();
    TPK.Advance(dT);
    TPK.CommitState(S);
    T += dT;
    iter++;

    if (iter < 10) {
      printf("T=%6.2f  C_0(x):", T);
      for (int k = 0; k < 15; k++) printf("%7.4f", (*tcc)[0][k]); cout << endl;
    }

    for (int k = 0; k < 19; k++) {
      CHECK(((*tcc)[0][k] - (*tcc)[0][k+1]) > -1e-15);
    }
  }

  for (int k = 0; k < 20; k++) {  // check that the final state is constant
    CHECK_CLOSE((*tcc)[0][k], 1.0, 1e-6);
  }

  delete comm;
}





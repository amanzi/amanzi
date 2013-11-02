/*
  The transport component of the Amanzi code, serial unit tests.
  License: BSD
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "MeshFactory.hh"
#include "GMVMesh.hh"

#include "State.hh"
#include "Transport_PK.hh"


/* **************************************************************** */
TEST(ADVANCE_WITH_SUBCYCLING) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;
  using namespace Amanzi::AmanziGeometry;

cout << "Test: Subcycling on a 2D square mesh" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  string xmlFileName = "test/transport_subcycling.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  /* create an MSTK mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, (Epetra_MpiComm *)comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory("test/rect2D_10x10_ss.exo", gm);
  
  /* create a simple state and populate it */
  std::vector<std::string> component_names;
  component_names.push_back("Component 0");
  component_names.push_back("Component 1");

  RCP<State> S = rcp(new State());
  S->RegisterDomainMesh(mesh);

  Transport_PK TPK(plist, S, component_names);
  TPK.CreateDefaultState(mesh, 2);

  /* modify the default state for the problem at hand */
  std::string passwd("state"); 
  Teuchos::RCP<Epetra_MultiVector> 
      flux = S->GetFieldData("darcy_flux", passwd)->ViewComponent("face", false);

  AmanziGeometry::Point velocity(1.0, 1.0);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    (*flux)[0][f] = velocity * normal;
  }

  /* initialize a transport process kernel */
  TPK.InitPK();
  TPK.PrintStatistics();

  /* advance the state */
  Teuchos::RCP<Epetra_MultiVector>
      tcc = S->GetFieldData("total_component_concentration", passwd)->ViewComponent("cell", false);

  double T = 0.0;
  int iter = 0;
  while (T < 1.0) {
    // imitation of a small time step relative to flow time step
    double dT = TPK.CalculateTransportDt();  
    double dT_MPC = dT * 7.7;

    TPK.Advance(dT_MPC);
    TPK.CommitState(S);
    T += dT_MPC;
    iter++;

    if (iter < 5) {
      printf("T=%8.4f  C_0(x):", T);
      for (int k = 0; k < 9; k++) {
        int k1 = 9 - k;  // reflects cell numbering in the exodus file
        printf("%7.4f", (*tcc)[0][k1]); 
      }
      printf("\n");
    }

    // for (int k = 0; k < 8; k++)
    //   CHECK( ((*tcc)[0][k+1] - (*tcc)[0][k]) > -1e-15 );
    if (iter == 5) {
      GMV::open_data_file(*mesh, (std::string)"transport.gmv");
      GMV::start_data();
      GMV::write_cell_data(*tcc, 0, "component0");
      GMV::write_cell_data(*tcc, 1, "component1");
      GMV::close_data_file();
    }
  }

  /* check that the final state is constant */
  for (int k = 0; k < 10; k++) 
    CHECK_CLOSE(1.0, (*tcc)[0][k], 1e-6);

  delete comm;
}






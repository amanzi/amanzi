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
// DEPRECATED #include "Teuchos_XMLParameterListHelpers.hpp"

#include "Mesh_simple.hh"
#include "MeshAudit.hh"

#include "State.hpp"
#include "Transport_PK.hpp"


/* **************************************************************** 
 * Test Init() procedure in the constructor.
 * ************************************************************* */
TEST(CONSTRUCTOR) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: read transport XML file" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  ParameterList parameter_list;
  string xmlFileName = "test/transport_mics.xml";
  // DEPRECATED  updateParametersFromXmlFile(xmlFileName, &parameter_list);

  ParameterXMLFileReader xmlreader(xmlFileName);
  parameter_list = xmlreader.getParameters();  
 
  /* create an MSTK mesh framework */
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, comm);
  RCP<Mesh> mesh = rcp(new Mesh_simple(0.0,0.0,0.0, 1.0,1.0,1.0, 1, 2, 1, comm, gm)); 
 
  //MeshAudit audit(mesh);
  //audit.Verify();

  /* create a Transport state with two components */
  int num_components = 2;
  State mpc_state(num_components, 0, mesh);
  RCP<Transport_State> TS = rcp(new Transport_State(mpc_state));

  ParameterList transport_list =  parameter_list.get<Teuchos::ParameterList>("Transport");
  Transport_PK TPK(transport_list, TS);
  TPK.InitPK();
  TPK.PrintStatistics();

  double cfl = TPK.cfl();
  CHECK(0 < cfl && cfl <= 1.0);
 
  delete comm;
}




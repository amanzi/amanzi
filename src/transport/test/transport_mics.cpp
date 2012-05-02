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
#include "Teuchos_XMLParameterListHelpers.hpp"

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
  updateParametersFromXmlFile(xmlFileName, &parameter_list);
 
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
  TPK.printStatistics();

  double cfl = TPK.cfl();
  CHECK(0 < cfl && cfl <= 1.0);
 
  delete comm;
}


/* **************************************************************** 
 * Test the geometric objects of an orthogonal mesh
 * ************************************************************* */
TEST(FACES_VOLUMES) {
  using namespace Teuchos;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: faces and volumes" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  // read parameter list
  ParameterList parameter_list;
  string xmlFileName = "test/transport_mics.xml";
  updateParametersFromXmlFile(xmlFileName, &parameter_list);
 
  // create an simple mesh framework
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, (Epetra_MpiComm *)comm);
  RCP<Mesh> mesh = rcp(new Mesh_simple(0.0,0.0,0.0, 1.0,1.0,1.0, 1, 2, 1, comm, gm)); 
 
  // create a transport state with two components
  int num_components = 2;
  State mpc_state(num_components, 0, mesh);
  RCP<Transport_State> TS = rcp(new Transport_State(mpc_state));
  ParameterList transport_list =  parameter_list.get<Teuchos::ParameterList>("Transport");
  Transport_PK TPK(transport_list, TS);

  // printing face areas
  int f;
  double area;
  const Epetra_Map& face_map = mesh->face_map(true);

  cout << "Face areas: ";
  for (f=face_map.MinLID(); f<=face_map.MaxLID(); f++) { 
    area = mesh->face_area(f);
    std::cout << area << " ";
    CHECK_CLOSE(area, 0.75, 0.25);
  }
  std::cout << endl;

  /* printing cell volumes */
  int  c;
  double volume;
  const Epetra_Map& cell_map = mesh->cell_map(true);

  std::cout << "Cell volumes: ";
  for (c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++) { 
    volume = mesh->cell_volume(c);
    std::cout << volume << " ";
    CHECK_EQUAL(volume, 0.5);
  }
  std::cout << endl;

  delete comm;
}




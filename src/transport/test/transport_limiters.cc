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

#include "State.hh"
#include "Transport_PK.hh"


/* **************************************************************** 
 * Test LimiterBarthJespersen()routine.
 * ************************************************************* */
TEST(LIMITER_BARTH_JESPERSEN) {
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

  // read parameter list
  ParameterList parameter_list;
  string xmlFileName = "test/transport_limiters.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  parameter_list = xmlreader.getParameters();  
 
  // create an MSTK mesh framework
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory factory(comm);
  factory.preference(pref);
  RCP<Mesh> mesh = factory(0.0,0.0,0.0, 1.0,1.0,1.0, 3, 4, 7, gm); 

  // create a Transport state with two components
  int num_components = 1;
  State mpc_state(num_components, 0, mesh);
  RCP<Transport_State> TS = rcp(new Transport_State(mpc_state));

  Point u(1.0, 0.0, 0.0);
  TS->AnalyticDarcyFlux(u);
  TS->AnalyticPorosity();
  TS->AnalyticWaterSaturation();
  TS->AnalyticWaterDensity();

  // create transport PK  
  ParameterList transport_list = parameter_list.get<Teuchos::ParameterList>("Transport");
  Transport_PK TPK(transport_list, TS);
  TPK.InitPK();
  TPK.PrintStatistics();
  double dT = TPK.CalculateTransportDt();  // We call it to identify upwind cells.

  // create a linear field
  const Epetra_Map& cmap = mesh->cell_map(false);
  RCP<Epetra_Vector> scalar_field = rcp(new Epetra_Vector(cmap));
  RCP<Epetra_MultiVector> gradient =Teuchos::rcp(new  Epetra_MultiVector(cmap, 3));

  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    (*scalar_field)[c] = 5.0 - xc[0] - 0.5 * xc[1] - 0.2 * xc[2];
    (*gradient)[0][c] = -1.0;
    (*gradient)[1][c] = -0.5;
    (*gradient)[2][c] = -0.2;
  }

  // calculate and verify limiters, 
  RCP<Epetra_Vector> limiter = rcp(new Epetra_Vector(cmap));
  TPK.LimiterBarthJespersen(0, scalar_field, gradient, limiter);

  for (int c = 0; c < ncells - 1; c++) {  // the corner cell gives limiter=0
    CHECK_CLOSE(1.0, (*limiter)[c], 1e-6);
  }
 
  delete comm;
}




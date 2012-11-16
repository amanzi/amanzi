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
// DEPRECATED #include "Teuchos_XMLParameterListHelpers.hpp"

#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "gmv_mesh.hh"

#include "State.hpp"
#include "Transport_PK.hpp"


Amanzi::AmanziGeometry::Point f_velocity(const Amanzi::AmanziGeometry::Point& x, double t ) { 
  return Amanzi::AmanziGeometry::Point(1.0, 0.5);
}


/* **************************************************************** */
TEST(ADVANCE_WITH_2D_MESH) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;
  using namespace Amanzi::AmanziGeometry;

cout << "Test: 2D transport on a square mesh for long time" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm  *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm  *comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  ParameterList parameter_list;
  string xmlFileName = "test/transport_2D_long.xml";
  // DEPRECATED updateParametersFromXmlFile(xmlFileName, &parameter_list);

  ParameterXMLFileReader xmlreader(xmlFileName);
  parameter_list = xmlreader.getParameters();

  /* create an MSTK mesh framework */
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, (Epetra_MpiComm *)comm);
  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(comm);
  meshfactory.preference(pref);
  RCP<Mesh> mesh = meshfactory("test/rect2D_50x50_ss.exo", gm);
  
  /* create a MPC state with one component */
  int num_components = 1;
  State mpc_state(num_components, 0, mesh);
 
  /* create a transport state from the MPC state and populate it */
  RCP<Transport_State> TS = rcp(new Transport_State(mpc_state));

  TS->AnalyticDarcyFlux(f_velocity);
  TS->AnalyticPorosity();
  TS->AnalyticWaterSaturation();
  TS->AnalyticWaterDensity();

  /* initialize a transport process kernel from the transport state */
  ParameterList transport_list =  parameter_list.get<Teuchos::ParameterList>("Transport");
  Transport_PK TPK(transport_list, TS);
  TPK.InitPK();
  TPK.PrintStatistics();
 
  /* advance the transport state */
  int iter, k;
  double T = 0.0;
  RCP<Transport_State> TS_next = TPK.transport_state_next();

  RCP<Epetra_MultiVector> tcc = TS->total_component_concentration();
  RCP<Epetra_MultiVector> tcc_next = TS_next->total_component_concentration();

  iter = 0;
  bool flag = true;
  while (T < 0.3) {
    double dT = TPK.CalculateTransportDt();
    TPK.Advance(dT);
    T += dT;
    iter++;

    if (T>0.2 && flag) {
      flag = false;
      if (TPK.MyPID == 0) {
        GMV::open_data_file(*mesh, (std::string)"transport.gmv");
        GMV::start_data();
        GMV::write_cell_data(*tcc_next, 0, "component0");
        //GMV::write_cell_data(*tcc_next, 1, "component1");
        //GMV::write_cell_data(*tcc_next, 2, "component2");
        GMV::close_data_file();
      }
      break;
    }

    *tcc = *tcc_next;
  }
  TPK.CheckTracerBounds(*tcc_next, 0, 0.0, 1.0, AmanziTransport::TRANSPORT_LIMITER_TOLERANCE);
 
  delete comm;
}






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
#include "gmv_mesh.hh"

#include "State.hh"
#include "Transport_PK.hh"


Amanzi::AmanziGeometry::Point f_velocity(const Amanzi::AmanziGeometry::Point& x, double t) { 
  return Amanzi::AmanziGeometry::Point(1.0, 1.0);
}


/* **************************************************************** */
TEST(ADVANCE_WITH_SUBCYCLING) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;
  using namespace Amanzi::AmanziGeometry;

cout << "Test: Subcycling on a 2D square mesh" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm  *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm  *comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  ParameterList parameter_list;
  string xmlFileName = "test/transport_subcycling.xml";
  // DEPRECATED  updateParametersFromXmlFile(xmlFileName, &parameter_list);

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
  RCP<Mesh> mesh = meshfactory("test/rect2D_10x10_ss.exo", gm);
  
  /* create a MPC state with two component */
  int num_components = 2;
  State mpc_state(num_components, 0, mesh);
 
  /* create a transport state from the MPC state and populate it */
  RCP<Transport_State> TS = rcp(new Transport_State(mpc_state));

  TS->AnalyticDarcyFlux(f_velocity);
  TS->AnalyticPorosity();
  TS->AnalyticWaterSaturation();
  TS->AnalyticWaterDensity();

  /* initialize a transport process kernel from a transport state */
  ParameterList transport_list = parameter_list.get<Teuchos::ParameterList>("Transport");
  Transport_PK TPK(transport_list, TS);
  TPK.InitPK();
  TPK.PrintStatistics();

  /* advance the state */
  int iter, k;
  double T = 0.0;
  RCP<Transport_State> TS_next = TPK.transport_state_next();
  RCP<Epetra_MultiVector> tcc = TS->total_component_concentration();
  RCP<Epetra_MultiVector> tcc_next = TS_next->total_component_concentration();

  iter = 0;
  while (T < 1.0) {
    // imitation of a small time step relative to flow time step
    double dT = TPK.CalculateTransportDt();  
    double dT_MPC = dT * 7.7;

    TPK.Advance(dT_MPC);
    T += dT_MPC;
    iter++;

    if (iter < 5) {
      printf("T=%8.4f  C_0(x):", T);
      for( int k=0; k<9; k++ ) {
        int k1 = 9 - k;  // reflects cell numbering in the exodus file
        printf("%7.4f", (*tcc_next)[0][k1]); 
      }
      printf("\n");
    }

    //for( int k=0; k<8; k++ )
    //  CHECK( ((*tcc_next)[0][k+1] - (*tcc_next)[0][k]) > -1e-15 );
    if (iter == 15) {
      GMV::open_data_file(*mesh, (std::string)"transport.gmv");
      GMV::start_data();
      GMV::write_cell_data(*tcc_next, 0, "component0");
      GMV::write_cell_data(*tcc_next, 1, "component1");
      GMV::close_data_file();
    }

    *tcc = *tcc_next;
  }

  /* check that the final state is constant */
  for (int k=0; k<10; k++) 
    CHECK_CLOSE(1.0, (*tcc_next)[0][k], 1e-6);

  delete comm;
}






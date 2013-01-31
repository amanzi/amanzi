#include <iostream>
#include <cstdlib>
#include <cmath>
#include "UnitTest++.h"
#include <vector>

#include "MeshFactory.hh"

#include "State.hpp"
#include "Transport_PK.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"


TEST(ADVANCE_WITH_STK) {
  using namespace Teuchos;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: advance with STK" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  // read parameter list
  ParameterList parameter_list;
  string xmlFileName = "test/transport_advance_stk.xml";
  // DEPRECATED updateParametersFromXmlFile(xmlFileName, &parameter_list);

  ParameterXMLFileReader xmlreader(xmlFileName);
  parameter_list = xmlreader.getParameters();

    // create an MSTK mesh framework 
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, (Epetra_MpiComm *)comm);
  FrameworkPreference pref;
  pref.clear();
  pref.push_back(STKMESH);

  MeshFactory meshfactory(comm);
  meshfactory.preference(pref);
  RCP<Mesh> mesh = meshfactory("test/hex_3x3x3_ss.exo", gm);
  
  // create a transport state with two component 
  int num_components = 2;
  State mpc_state(num_components, 0, mesh);
  RCP<Transport_State> TS = rcp(new Transport_State(mpc_state));

  Point u(1.0, 0.0, 0.0);
  TS->AnalyticDarcyFlux(u);
  TS->AnalyticPorosity();
  TS->AnalyticWaterSaturation();
  TS->AnalyticWaterDensity();

  ParameterList transport_list = parameter_list.get<Teuchos::ParameterList>("Transport");
  Transport_PK TPK(transport_list, TS);
  TPK.InitPK();

  // advance the state
  double dT = TPK.CalculateTransportDt();
  TPK.Advance(dT);

  // printing cell concentration  
  int i, k;
  double T = 0.0;
  RCP<Transport_State> TS_next = TPK.transport_state_next();
  RCP<Epetra_MultiVector> tcc = TS->total_component_concentration();
  RCP<Epetra_MultiVector> tcc_next = TS_next->total_component_concentration();

  for (i = 0; i < 50; i++) {
    dT = TPK.CalculateTransportDt();
    TPK.Advance(dT);
    T += dT;

    if (i < 10) {
      printf("T=%6.2f  C_0(x):", T);
      for (int k=0; k<9; k++) printf("%7.4f", (*tcc_next)[0][k]); cout << endl;
    }
     *tcc = *tcc_next;
  }

  // check that the final state is constant  
  for (int k=0; k<4; k++) 
    CHECK_CLOSE((*tcc_next)[0][k], 1.0, 1e-6);

  delete gm;
  delete comm;
}
 



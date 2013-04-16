#include <cstdlib>
#include <iostream>

#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Mesh_MSTK.hh"
#include "state.hh"

#include "overland.hh"

#include "flow_test_class.hh"

#include "../my_macro.hh"

/* **************************************************************** */
TEST(ADVANCE_WITH_SIMPLE) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: advance using simple mesh" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm  *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm  *comm = new Epetra_SerialComm();
#endif

  // read parameter list
  Teuchos::ParameterList parameter_list;
  //string xmlFileName = "test/overland_advance_simple.xml";
  string xmlFileName = "test/overland_advance_simple_2Dtestcase.xml";
  updateParametersFromXmlFile(xmlFileName, &parameter_list);

  // create an SIMPLE mesh framework
  Teuchos::ParameterList region_list =
    parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, (Epetra_MpiComm *)comm);
  //Teuchos::RCP<Mesh> mesh = Teuchos::rcp(new Mesh_MSTK(0.0,0.0, 182.87999414784019,1.0, 100, 1, comm, gm));
  Teuchos::RCP<Mesh> mesh = Teuchos::rcp(new Mesh_MSTK(0.0,0.0,1620,1000.0,81,50, comm, gm));

  // create and initialize the test class
  FlowTestOne test(parameter_list, mesh, 1);
  test.initialize();

  // advance the state
  int iter, k;
  double p = 0.0;

  Teuchos::RCP<const CompositeVector> pressure =
    test.S1->GetFieldData("overland_pressure");
  double time = test.S1->time();

  iter = 0;
  if (iter < 10) {
    printf( "time=%6.2f  p(x):", time );
    for( int k=0; k<15; k++ ) printf("%7.4f", (*pressure)("cell",k));
    cout << endl;
  }

  int n = 0 ;

  // test 1D
  //double t_final = 3600. ;

  // test 2D
  double t_final = 10400. ;

  double L1, L2;
  while (time < t_final ) {
    double dtime = test.FPK->get_dt();
    test.FPK->advance(dtime);
    time += dtime;
    iter++;
    //if ( iter==3 ) { exit(0); }
    //if ( time>1800 ) { exit(0); }
    test.commit_step();
   }
  delete comm;
}

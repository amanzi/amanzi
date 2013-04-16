#include <cstdlib>
#include <iostream>

#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Mesh_simple.hh"
#include "state.hh"

#include "richards.hh"

#include "flow_test_class.hh"

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
  string xmlFileName = "test/richards_advance_simple.xml";
  updateParametersFromXmlFile(xmlFileName, &parameter_list);

  // create an SIMPLE mesh framework
  Teuchos::ParameterList region_list =
    parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, (Epetra_MpiComm *)comm);
  Teuchos::RCP<Mesh> mesh = Teuchos::rcp(new Mesh_simple(0.0,0.0,0.0, 1.0,1.0,1.0, 20, 20, 2, comm, gm));

  // create and initialize the test class
  FlowTestOne test(parameter_list, mesh, 1);
  test.initialize();

  // advance the state
  int iter, k;
  double p = 0.0;

  Teuchos::RCP<const CompositeVector> pressure =
    test.S1->GetFieldData("pressure");
  double time = test.S1->time();

  iter = 0;
  if (iter < 10) {
    printf( "time=%6.2f  p(x):", time );
    for( int k=0; k<15; k++ ) printf("%7.4f", (*pressure)("cell",k));
    cout << endl;
  }

  double L1, L2;
  while (time < 1.0) {
    double dtime = test.FPK->get_dt();
    test.FPK->advance(dtime);
    time += dtime;
    iter++;

    if (iter < 10) {
      printf( "time=%6.2f  p(x):", time );
      for( int k=0; k<15; k++ ) printf("%7.4f", (*pressure)("cell",k));
      cout << endl;
    }

    test.evaluate_error_pressure(time, &L1, &L2);
    CHECK_CLOSE(0., L2, 1e-6);
    CHECK_CLOSE(0., L1, 1e-6);
    test.commit_step();
  }
  delete comm;
}

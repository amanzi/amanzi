#include <cstdlib>
#include <iostream>

#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Mesh_simple.hh"
#include "State.hh"

#include "MeshFactory.hh"
#include "Mesh.hh"
#include "Domain.hh"
#include "GeometricModel.hh"

#include "advection_diffusion.hh"

#include "energy_test_class.hh"

/* **************************************************************** */
void RunTest(std::string filename, std::string testname) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: advance using simple mesh" << std::endl;
#ifdef HAVE_MPI
  Epetra_MpiComm  *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm  *comm = new Epetra_SerialComm();
#endif

  // read parameter list
  Teuchos::RCP<Teuchos::ParameterList> parameter_list = Teuchos::getParametersFromXmlFile(filename);

  // create an SIMPLE mesh framework
  Teuchos::ParameterList region_list =
    parameter_list->get<Teuchos::ParameterList>("Regions");
  //GeometricModelPtr gm = new GeometricModel(3, region_list, (Epetra_MpiComm *)comm);
  GeometricModelPtr gm  = new GeometricModel(3, region_list, (Epetra_MpiComm *)comm);
  Teuchos::RCP<Mesh> mesh = Teuchos::rcp(new Mesh_simple(0.0, 0.0, 0.0, 10.0, 1.0, 1.0, 20, 1, 1, comm, gm));

  // create and initialize the test class
  Teuchos::RCP<EnergyTest> test;
  if (testname == "one") {
    test = Teuchos::rcp(new EnergyTestOne(parameter_list, mesh, 1));
  } else if (testname == "step") {
    test = Teuchos::rcp(new EnergyTestStep(parameter_list, mesh, 1));
  } else if (testname == "diffused step") {
    test = Teuchos::rcp(new EnergyTestDiffusedStep(parameter_list, mesh, 1));
  } else if (testname == "advected diffused step") {
    test = Teuchos::rcp(new EnergyTestAdvDiffusedStep(parameter_list, mesh, 1));
  }
  test->initialize();

  // advance the state
  int iter, k;
  double T = 0.0;


  const Epetra_MultiVector& temp = *test->S1->GetFieldData("temperature")->ViewComponent("cell",false);

  iter = 0;
  if (iter < 10) {
    printf( "T=%6.2f", T);
    if (test->S1->GetFieldData("temperature")->HasComponent("face")) {
      printf(" F(x=0)=%6.2f", (*test->S1->GetFieldData("temperature")->ViewComponent("face",false))[0][80]);
    }
    printf(" C(x):");
    for( int k=0; k<20; k++ ) printf("%7.4f", temp[0][k]);
    std::cout << std::endl;
  }

  double L1, L2;
  while (T < 1.0) {
    double dT = test->EPK->get_dt();
    test->S1->advance_cycle();
    test->S1->advance_time(dT);
    test->EPK->advance(dT);
    T += dT;
    iter++;

    if (iter < 10) {
      printf( "T=%6.2f", T);
      if (test->S1->GetFieldData("temperature")->HasComponent("face")) {
        printf(" F(x=0)=%6.2f  ", (*test->S1->GetFieldData("temperature")->ViewComponent("face",false))[0][80]);
      }
      printf("C(x):");
      for( int k=0; k<20; k++ ) printf("%7.4f", temp[0][k]);
      std::cout << std::endl;
    }

    test->evaluate_error_temp(T, &L1, &L2);
    // CHECK_CLOSE(0., L2, 1e-6);
    // CHECK_CLOSE(0., L1, 1e-6);
    test->commit_step(dT);
  }
  delete comm;
}

// solution: u = 1
TEST(ADV_DIFF_ONE_FV) {
  std::cout << "Advance ADV-DIFF problem with solution = 1" << std::endl;
  RunTest("test/advection_diffusion_advance_simple_fv.xml", "one");
}

// solution (forward time upwind space): slightly diffused interface
TEST(ADV_DIFF_ADVECTED_STEP_FV) {
  std::cout << "Advance ADV-DIFF problem with an advected step function" << std::endl;
  RunTest("test/advection_diffusion_advance_simple_fv.xml", "step");
}

// solution (backward time upwind space): implicit advection is MORE diffusive then explicit
TEST(ADV_DIFF_ADV_DIFFUSION_FV_IMPLICIT) {
  std::cout << "Advance ADV problem with advected and  diffused flux across the domain" << std::endl;
  RunTest("test/advection_diffusion_advance_simple_fv_implicit.xml", "step");
}

// solution: diffusion
TEST(ADV_DIFF_DIFFUSION_FV) {
  std::cout << "Advance ADV-DIFF problem with diffused flux across the domain" << std::endl;
  RunTest("test/advection_diffusion_advance_simple_fv.xml", "diffused step");
}

// solution: advection and diffusion
TEST(ADV_DIFF_ADV_DIFFUSION_FV) {
  std::cout << "Advance ADV-DIFF problem with advected and  diffused flux across the domain" << std::endl;
  RunTest("test/advection_diffusion_advance_simple_fv.xml", "advected diffused step");
}

// solution: advection and diffusion, but with constant flux.  not the same
// solution as the previous, as in the above the diffusive flux changes as a function
// of time
TEST(ADV_DIFF_ADVECTED_STEP_FV_NEUMANN) {
  std::cout << "Advance ADV-DIFF problem with advected flux across the domain, using a Neumann BC" << std::endl;
  RunTest("test/advection_diffusion_advance_simple_fv_neumann.xml", "advected diffused step");
}


// repeat all with MFD
TEST(ADV_DIFF_ONE_MFD) {
  std::cout << "Advance ADV-DIFF problem with solution = 1" << std::endl;
  RunTest("test/advection_diffusion_advance_simple_mfd.xml", "one");
}

TEST(ADV_DIFF_ADVECTED_STEP_MFD) {
  std::cout << "Advance ADV-DIFF problem with an advected step function" << std::endl;
  RunTest("test/advection_diffusion_advance_simple_mfd.xml", "step");
}

TEST(ADV_DIFF_ADVECTED_STEP_MFD_IMPLICIT) {
  std::cout << "Advance ADV-DIFF problem with an advected step function" << std::endl;
  RunTest("test/advection_diffusion_advance_simple_mfd_implicit.xml", "step");
}

TEST(ADV_DIFF_DIFFUSION_MFD) {
  std::cout << "Advance ADV-DIFF problem with diffused flux across the domain" << std::endl;
  RunTest("test/advection_diffusion_advance_simple_mfd.xml", "diffused step");
}


TEST(ADV_DIFF_ADV_DIFFUSION_MFD) {
  std::cout << "Advance ADV-DIFF problem with advected and  diffused flux across the domain" << std::endl;
  RunTest("test/advection_diffusion_advance_simple_mfd.xml", "advected diffused step");
}

TEST(ADV_DIFF_ADVECTED_STEP_MFD_NEUMANN) {
  std::cout << "Advance ADV-DIFF problem with diffused flux across the domain, using a Neumann BC" << std::endl;
  RunTest("test/advection_diffusion_advance_simple_mfd_neumann.xml", "advected diffused step");
}

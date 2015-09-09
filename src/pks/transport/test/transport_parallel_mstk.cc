/*
  The transport component of the Amanzi code, serial unit tests.
  License: BSD
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "State.hh"

// Transport
#include "Transport_PK.hh"


double f_step(const Amanzi::AmanziGeometry::Point& x, double t) { 
  if (x[0] <= t) return 1;
  return 0;
}


TEST(ADVANCE_WITH_MSTK_PARALLEL) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: advance with MSTK in parallel" << std::endl;
#ifdef HAVE_MPI
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  std::string xmlFileName = "test/transport_parallel_mstk.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  /* create an MSTK mesh framework */
  ParameterList region_list = plist->get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, (Epetra_MpiComm *)comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory("test/hex_3x3x3_ss.exo", gm);
  
  /* create a simple state and populate it */
  std::vector<std::string> component_names;
  component_names.push_back("Component 0");
  component_names.push_back("Component 1");

  Teuchos::ParameterList state_list = plist->sublist("State");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
  S->set_time(0.0);
  S->set_intermediate_time(0.0);

  Transport_PK TPK(plist, S, "Transport", component_names);
  TPK.Setup();
  TPK.CreateDefaultState(mesh, 2);
  S->InitializeFields();
  S->InitializeEvaluators();

  /* modify the default state for the problem at hand */
  std::string passwd("state"); 
  Teuchos::RCP<Epetra_MultiVector> 
      flux = S->GetFieldData("darcy_flux", passwd)->ViewComponent("face", false);

  AmanziGeometry::Point velocity(1.0, 0.0, 0.0);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    (*flux)[0][f] = velocity * normal;
  }

  Teuchos::RCP<Epetra_MultiVector>
      tcc = S->GetFieldData("total_component_concentration", passwd)->ViewComponent("cell", false);

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    (*tcc)[0][c] = f_step(xc, 0.0);
  }

  /* initialize a transport process kernel */
  TPK.Initialize();

  /* advance the state */
  double t_old(0.0), t_new, dt;
  dt = TPK.CalculateTransportDt();
  t_new = t_old + dt;
  TPK.AdvanceStep(t_old, t_new);

  /* print cell concentrations */
  int iter = 0;
  while(t_new < 1.0) {
    dt = TPK.CalculateTransportDt();
    t_new = t_old + dt;

    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new);

    t_old = t_new;
    iter++;

    if (iter < 10 && TPK.MyPID == 3) {
      printf("T=%7.2f  C_0(x):", t_new);
      for (int k = 0; k < 2; k++) printf("%7.4f", (*tcc)[0][k]); std::cout << std::endl;
    }
  }

  //for (int k=0; k<12; k++) 
  //  CHECK_CLOSE((*tcc_next)[0][k], 1.0, 1e-6);
}
 
 



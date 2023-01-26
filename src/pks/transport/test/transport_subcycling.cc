/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport

*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

// TPLs
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "GMVMesh.hh"
#include "State.hh"

// Transport
#include "TransportExplicit_PK.hh"

/* **************************************************************** */
TEST(ADVANCE_WITH_SUBCYCLING)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: Subcycling on a 2D square mesh" << std::endl;
#ifdef HAVE_MPI
  Comm_ptr_type comm = Amanzi::getDefaultComm();
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  std::string xmlFileName = "test/transport_subcycling.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  /* create a mesh framework */
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  RCP<const Mesh> mesh = meshfactory.create("test/rect2D_10x10_ss.exo");

  /* create a simple state and populate it */
  std::vector<std::string> component_names;
  component_names.push_back("Component 0");
  component_names.push_back("Component 1");

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  TransportExplicit_PK TPK(plist, S, "transport", component_names);
  TPK.Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  S->set_time(0.0);
  S->set_intermediate_time(0.0);

  /* modify the default state for the problem at hand */
  std::string passwd("state");
  auto& flux = *S->GetW<CompositeVector>("volumetric_flow_rate", passwd).ViewComponent("face");

  AmanziGeometry::Point velocity(1.0, 1.0);
  int nfaces_owned =
    mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->getFaceNormal(f);
    flux[0][f] = velocity * normal;
  }

  /* initialize a transport process kernel */
  TPK.Initialize();

  /* advance the state */
  auto tcc =
    S->GetW<CompositeVector>("total_component_concentration", passwd).ViewComponent("cell");

  double t_old(0.0), t_new(0.0), dt;
  int iter = 0;
  while (t_new < 1.0) {
    // imitation of a small time step relative to flow time step
    dt = TPK.StableTimeStep(-1);
    double dt_MPC = dt * 7.7;
    t_new = t_old + dt_MPC;

    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;
    iter++;

    if (iter < 5) {
      printf("T=%8.4f  C_0(x):", t_new);
      for (int k = 0; k < 9; k++) {
        int k1 = 9 - k; // reflects cell numbering in the exodus file
        printf("%7.4f", (*tcc)[0][k1]);
      }
      printf("\n");
    }

    // for (int k = 0; k < 8; k++)
    //   CHECK( ((*tcc)[0][k+1] - (*tcc)[0][k]) > -1e-15 );
    if (iter == 5) {
      GMV::open_data_file(*mesh, (std::string) "transport.gmv");
      GMV::start_data();
      GMV::write_cell_data(*tcc, 0, "component0");
      GMV::write_cell_data(*tcc, 1, "component1");
      GMV::close_data_file();
    }
  }

  /* check that the final state is constant */
  for (int k = 0; k < 10; k++) CHECK_CLOSE(1.0, (*tcc)[0][k], 1e-6);
}

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

// TPLs
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "IO.hh"
#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "Point.hh"
#include "State.hh"

// Transport
#include "TransportExplicit_PK.hh"

/* **************************************************************** */
TEST(DIFFUSION_GAS_SMILES) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "\nTest: diffusion gas (Smiles)" << std::endl;
  Comm_ptr_type comm = Amanzi::getDefaultComm();

  std::string xmlFileName = "test/transport_diffusion_smiles.xml";
  auto plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create an MSTK mesh framework
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<const Mesh> mesh = meshfactory.create(0.0,0.0, 10.0,1.0, 51, 1); 

  // create a simple state and populate it
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  std::vector<std::string> component_names;
  component_names.push_back("Tritium(l)");
  component_names.push_back("Tritium(g)");

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  TransportExplicit_PK TPK(plist, S, "transport", component_names);
  TPK.Setup();
  TPK.CreateDefaultState(mesh, 1);
  S->InitializeFields();
  S->InitializeEvaluators();
  S->set_time(0.0);
  S->set_intermediate_time(0.0);
  S->set_initial_time(0.0);
  S->set_final_time(0.0);

  // modify the default state for the problem at hand
  std::string passwd("state"); 
  double sl(0.01);
  S->GetW<CompositeVector>("prev_saturation_liquid", passwd).PutScalar(sl);
  S->GetW<CompositeVector>("saturation_liquid", passwd).PutScalar(sl);
  auto& flux = *S->GetW<CompositeVector>("volumetric_flow_rate", passwd).ViewComponent("face");
  flux.PutScalar(0.0);

  auto& tcc = *S->GetW<CompositeVector>("total_component_concentration", passwd).ViewComponent("cell");

  // initialize a transport process kernel
  Amanzi::VerboseObject::global_hide_line_prefix = true;
  TPK.Initialize();

  // advance the state
  int iter = 0;
  double dt0(2 * 3.1558e+6), t_old(0.0), t_new(0.0), dt, T1(3.1558e+8);
  while (t_new < T1) {
    dt = std::min(TPK.StableTimeStep(), T1 - t_old);
    dt = std::min(dt, dt0);
    t_new = t_old + dt;

    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;
    iter++;
    printf("AAA: %g %16.10g %16.10g ... %16.10g   %16.10g\n", t_new, tcc[0][0], tcc[0][1], tcc[0][25], tcc[1][25]);
  }

  auto vo = Teuchos::rcp(new Amanzi::VerboseObject("Smiles", *plist));
  WriteStateStatistics(*S, *vo);

  // check for values
  double tcc_eff = tcc[0][25] + tcc[1][25] * (1 - sl) / sl;
  printf("Normalized effective liquid concentration: %g \n", tcc_eff);
  CHECK_CLOSE(0.14, tcc_eff, 0.002);

  // check for bounds
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int c = 0; c < ncells_owned; ++c) {
    CHECK(tcc[0][c] >= 0.0 && tcc[0][c] <= 1.0);
    CHECK(tcc[1][c] >= 0.0 && tcc[1][c] <= 1.0);
  }
}




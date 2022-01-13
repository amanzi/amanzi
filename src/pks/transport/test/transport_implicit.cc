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
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "State.hh"

// Transport
#include "TransportImplicit_PK.hh"


TEST(ADVANCE_WITH_MESH_FRAMEWORK) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::string framework_name = "MSTK";
  Framework framework = Framework::MSTK;
  
  std::cout << "Test: implicit advance "<< std::endl;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  // read parameter list
  std::string xmlFileName("test/transport_implicit.xml");
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(framework);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  RCP<const Mesh> mesh;

  mesh = meshfactory.create("test/hex_3x3x3_ss.exo");
  // mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 10);
  
  // create a simple state and populate it
  Amanzi::VerboseObject::global_hide_line_prefix = false;

  std::vector<std::string> component_names;
  component_names.push_back("Component 0");

  Teuchos::ParameterList state_list = plist->sublist("state");  
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Teuchos::ParameterList pk_tree = plist->sublist("cycle driver").sublist("pk_tree")
                                        .sublist("transport implicit");

  // create the global solution vector
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  TransportImplicit_PK TPK(pk_tree, plist, S, soln);

  TPK.Setup();
  TPK.CreateDefaultState(mesh, 2);
  S->InitializeFields();
  S->InitializeEvaluators();
  S->set_time(0.0);
  S->set_intermediate_time(0.0);

  // modify the default state for the problem at hand
  std::string passwd("state"); 
  auto& flux = *S->GetW<CompositeVector>("darcy_flux", passwd).ViewComponent("face");

  AmanziGeometry::Point velocity(1.0, 1.0, 0.0);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    flux[0][f] = velocity * normal;
  }

  // initialize a transport process kernel
  TPK.Initialize();

  // advance the state
  double t_old(0.0), t_new, dt;
  dt = 0.01;
  t_new = t_old + dt;
  TPK.AdvanceStep(t_old, t_new);
  TPK.CommitStep(t_old, t_new, Tags::DEFAULT);
  
  //printing cell concentration
  auto tcc = S->GetW<CompositeVector>("total_component_concentration", passwd).ViewComponent("cell");

  while(t_new < 1.2) {
    t_new = t_old + dt;
    
    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new, Tags::DEFAULT);
    
    t_old = t_new;
 
    if (t_new < 0.4) {
      printf("T=%6.2f  C_0(x):", t_new);
      for (int k = 0; k < 9; k++) printf("%7.4f", (*tcc)[0][k]); std::cout << std::endl;
    }
  }

  // check that the final state is constant
  for (int k = 0; k < 4; k++) 
    CHECK_CLOSE((*tcc)[0][k], 1.0, 1e-5);
}
 



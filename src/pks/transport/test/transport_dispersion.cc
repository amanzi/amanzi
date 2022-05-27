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
#include "GMVMesh.hh"
#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "Point.hh"
#include "State.hh"

// Transport
#include "TransportExplicit_PK.hh"


double f_step(const Amanzi::AmanziGeometry::Point& x, double t) { 
  if (x[0] <= 1 + t) return 1;
  return 0;
}

double f_smooth(const Amanzi::AmanziGeometry::Point& x, double t) { 
  return 0.5 - atan(50*(x[0]-5-t)) / M_PI;
}

double f_cubic(const Amanzi::AmanziGeometry::Point& x, double t) {
  if (x[0] < 1 + t) return 1;
  if (x[0] > 3 + t) return 0;
  double z = (x[0]-1-t) / 2;
  return 2*z*z*z - 3*z*z + 1;
}


/* **************************************************************** */
TEST(DISPERSION) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: dispersion" << std::endl;
#ifdef HAVE_MPI
  Comm_ptr_type comm = Amanzi::getDefaultComm();
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  std::string xmlFileName = "test/transport_dispersion.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  /* create an MSTK mesh framework */
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(pref);
  int nx = 20;
  RCP<const Mesh> mesh = meshfactory.create(0.0,0.0,0.0, 5.0,1.0,1.0, nx, 10, 1); 

  /* create a simple state and populate it */
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  std::vector<std::string> component_names;
  component_names.push_back("Component 0");

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

  /* modify the default state for the problem at hand */
  std::string passwd("state"); 
  auto& flux = *S->GetW<CompositeVector>("volumetric_flow_rate", passwd).ViewComponent("face");

  AmanziGeometry::Point velocity(1.0, 0.0, 0.0);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    flux[0][f] = velocity * normal;
  }

  auto tcc = S->GetW<CompositeVector>("total_component_concentration", passwd).ViewComponent("cell");

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    (*tcc)[0][c] = f_step(xc, 0.0);
  }

  S->GetW<double>("const_fluid_density", passwd) = 1.0;

  /* initialize a transport process kernel */
  Amanzi::VerboseObject::global_hide_line_prefix = true;
  TPK.Initialize();

  /* advance the state */
  double dt0;
  dt0 = TPK.StableTimeStep();

  int iter = 0;
  double t_old(0.0), t_new(0.0), dt, T1(1.0);
  while (t_new < T1) {
    dt = std::min(TPK.StableTimeStep(), T1 - t_old);
    dt = std::min(dt, dt0);
    t_new = t_old + dt;
    // for (int k = 0; k < nx; k++) printf("%10.8f\n", (*tcc)[0][k]); 
    // printf("\n");

    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;
    iter++;

    TPK.VV_CheckTracerBounds(*tcc, 0, 0.0, 1.0, 1e-12);
  }
}


/* **************************************************************** */
TEST(DIFFUSION) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "\nTest: diffusion" << std::endl;
#ifdef HAVE_MPI
  Comm_ptr_type comm = Amanzi::getDefaultComm();
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  std::string xmlFileName = "test/transport_diffusion.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  /* create an MSTK mesh framework */
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(pref);
  RCP<const Mesh> mesh = meshfactory.create(0.0,0.0, 1.0,1.0, 20, 20); 

  /* create a simple state and populate it */
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  std::vector<std::string> component_names;
  component_names.push_back("Component 0");

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

  /* modify the default state for the problem at hand */
  std::string passwd("state"); 
  auto& flux = *S->GetW<CompositeVector>("volumetric_flow_rate", passwd).ViewComponent("face");

  AmanziGeometry::Point velocity(0.5, 0.0);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    flux[0][f] = velocity * normal;
  }

  auto tcc = S->GetW<CompositeVector>("total_component_concentration", passwd).ViewComponent("cell");

  S->GetW<double>("const_fluid_density", passwd) = 1.0;

  /* initialize a transport process kernel */
  Amanzi::VerboseObject::global_hide_line_prefix = true;
  TPK.Initialize();

  /* advance the state */
  double dt0;
  dt0 = TPK.StableTimeStep();

  int iter = 0;
  double t_old(0.0), t_new(0.0), dt, T1(1.0);
  while (t_new < T1) {
    dt = std::min(TPK.StableTimeStep(), T1 - t_old);
    dt = std::min(dt, dt0);
    t_new = t_old + dt;
    // for (int k = 0; k < nx; k++) printf("%10.8f\n", (*tcc)[0][k]); 
    // printf("\n");

    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;
    iter++;

    TPK.VV_CheckTracerBounds(*tcc, 0, 0.0, 1.0, 1e-12);
  }
 
  GMV::open_data_file(*mesh, (std::string)"transport.gmv");
  GMV::start_data();
  GMV::write_cell_data(*tcc, 0, "Component_0");
  GMV::close_data_file();
}


/* **************************************************************** */
TEST(GAS_DIFFUSION) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "\nTest: gas diffusion" << std::endl;
#ifdef HAVE_MPI
  Comm_ptr_type comm = Amanzi::getDefaultComm();
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  std::string xmlFileName = "test/transport_diffusion_gas.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  /* create an MSTK mesh framework */
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(pref);
  RCP<const Mesh> mesh = meshfactory.create(0.0,0.0, 1.0,1.0, 21, 21); 

  /* create a simple state and populate it */
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  std::vector<std::string> component_names;
  component_names.push_back("Component 0");
  component_names.push_back("Component 1");

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

  /* modify the default state for the problem at hand */
  std::string passwd("state"); 
  S->GetW<CompositeVector>("prev_saturation_liquid", passwd).PutScalar(0.4);
  S->GetW<CompositeVector>("saturation_liquid", passwd).PutScalar(0.4);
  auto& flux = *S->GetW<CompositeVector>("volumetric_flow_rate", passwd).ViewComponent("face");

  AmanziGeometry::Point velocity(0.1, 0.0);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    flux[0][f] = velocity * normal;
  }

  auto& tcc = *S->GetW<CompositeVector>("total_component_concentration", passwd).ViewComponent("cell");

  S->GetW<double>("const_fluid_density", passwd) = 1.0;

  /* initialize a transport process kernel */
  Amanzi::VerboseObject::global_hide_line_prefix = true;
  TPK.Initialize();

  /* advance the state */
  double dt0;
  dt0 = TPK.StableTimeStep();

  int iter = 0;
  double t_old(0.0), t_new(0.0), dt, T1(0.99);
  while (t_new < T1) {
    dt = std::min(TPK.StableTimeStep(), T1 - t_old);
    dt = std::min(dt, dt0);
    t_new = t_old + dt;

    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;
    iter++;

    TPK.VV_CheckTracerBounds(tcc, 0, 0.0, 1.0, 1e-12);
    TPK.VV_CheckTracerBounds(tcc, 1, 0.0, 1.0, 1e-12);
  }

  // check symmetry of diffusion
  CHECK_CLOSE(tcc[1][199], tcc[1][241], 1e-12);
  CHECK_CLOSE(tcc[1][218], tcc[1][222], 1e-12);

  // check for bounds
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int c = 0; c < ncells_owned; ++c) {
    CHECK(tcc[0][c] >= 0.0 && tcc[0][c] <= 1.0);
    CHECK(tcc[1][c] >= 0.0 && tcc[1][c] <= 1.0);
  }
 
  GMV::open_data_file(*mesh, (std::string)"transport.gmv");
  GMV::start_data();
  GMV::write_cell_data(tcc, 0, "Component_0");
  GMV::write_cell_data(tcc, 1, "Component_1");
  GMV::close_data_file();
}




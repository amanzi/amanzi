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
#include "Transport_PK.hh"


double f_step(const Amanzi::AmanziGeometry::Point& x, double t) { 
  if (x[0] <= 1 + t) return 1;
  return 0;
}

double f_smooth(const Amanzi::AmanziGeometry::Point& x, double t) { 
  return 0.5 - atan(50*(x[0]-5-t)) / M_PI;
}

double f_cubic(const Amanzi::AmanziGeometry::Point& x, double t) {
  if( x[0] < 1 + t ) return 1;
  if( x[0] > 3 + t ) return 0;
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
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  std::string xmlFileName = "test/transport_dispersion.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  /* create an MSTK mesh framework */
  ParameterList region_list = plist->get<Teuchos::ParameterList>("Regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, comm));

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory factory(comm);
  factory.preference(pref);
  int nx = 20;
  RCP<const Mesh> mesh = factory(0.0,0.0,0.0, 5.0,1.0,1.0, nx, 10, 1, gm); 

  /* create a simple state and populate it */
  Amanzi::VerboseObject::hide_line_prefix = true;

  std::vector<std::string> component_names;
  component_names.push_back("Component 0");

  Teuchos::ParameterList state_list = plist->sublist("State");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
  S->set_time(0.0);
  S->set_intermediate_time(0.0);
  S->set_initial_time(0.0);
  S->set_final_time(0.0);

  Transport_PK TPK(plist, S, "Transport", component_names);
  TPK.Setup();
  TPK.CreateDefaultState(mesh, 1);
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

  *(S->GetScalarData("fluid_density", passwd)) = 1.0;

  /* initialize a transport process kernel */
  Amanzi::VerboseObject::hide_line_prefix = true;
  TPK.Initialize();

  /* advance the state */
  double dt0;
  dt0 = TPK.CalculateTransportDt();

  int i, k, iter = 0;
  double t_old(0.0), t_new(0.0), dt, T1(1.0);
  while (t_new < T1) {
    dt = std::min(TPK.CalculateTransportDt(), T1 - t_old);
    dt = std::min(dt, dt0);
    t_new = t_old + dt;
    // for (int k = 0; k < nx; k++) printf("%10.8f\n", (*tcc)[0][k]); 
    // printf("\n");

    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new);

    t_old = t_new;
    iter++;

    TPK.VV_CheckTracerBounds(*tcc, 0, 0.0, 1.0, 1e-12);
  }
 
  delete comm;
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
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  std::string xmlFileName = "test/transport_diffusion.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  /* create an MSTK mesh framework */
  ParameterList region_list = plist->get<Teuchos::ParameterList>("Regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, comm));

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory factory(comm);
  factory.preference(pref);
  RCP<const Mesh> mesh = factory(0.0,0.0, 1.0,1.0, 20, 20, gm); 

  /* create a simple state and populate it */
  Amanzi::VerboseObject::hide_line_prefix = true;

  std::vector<std::string> component_names;
  component_names.push_back("Component 0");

  Teuchos::ParameterList state_list = plist->sublist("State");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
  S->set_time(0.0);
  S->set_intermediate_time(0.0);
  S->set_initial_time(0.0);
  S->set_final_time(0.0);

  Transport_PK TPK(plist, S, "Transport", component_names);
  TPK.Setup();
  TPK.CreateDefaultState(mesh, 1);
  S->InitializeFields();
  S->InitializeEvaluators();

  /* modify the default state for the problem at hand */
  std::string passwd("state"); 
  Teuchos::RCP<Epetra_MultiVector> 
      flux = S->GetFieldData("darcy_flux", passwd)->ViewComponent("face", false);

  AmanziGeometry::Point velocity(0.5, 0.0);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    (*flux)[0][f] = velocity * normal;
  }

  Teuchos::RCP<Epetra_MultiVector> 
      tcc = S->GetFieldData("total_component_concentration", passwd)->ViewComponent("cell", false);

  *(S->GetScalarData("fluid_density", passwd)) = 1.0;

  /* initialize a transport process kernel */
  Amanzi::VerboseObject::hide_line_prefix = true;
  TPK.Initialize();

  /* advance the state */
  double dt0;
  dt0 = TPK.CalculateTransportDt();

  int i, k, iter = 0;
  double t_old(0.0), t_new(0.0), dt, T1(1.0);
  while (t_new < T1) {
    dt = std::min(TPK.CalculateTransportDt(), T1 - t_old);
    dt = std::min(dt, dt0);
    t_new = t_old + dt;
    // for (int k = 0; k < nx; k++) printf("%10.8f\n", (*tcc)[0][k]); 
    // printf("\n");

    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new);

    t_old = t_new;
    iter++;

    TPK.VV_CheckTracerBounds(*tcc, 0, 0.0, 1.0, 1e-12);
  }
 
  GMV::open_data_file(*mesh, (std::string)"transport.gmv");
  GMV::start_data();
  GMV::write_cell_data(*tcc, 0, "Component_0");
  GMV::close_data_file();

  delete comm;
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
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  std::string xmlFileName = "test/transport_diffusion_gas.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  /* create an MSTK mesh framework */
  ParameterList region_list = plist->get<Teuchos::ParameterList>("Regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, comm));

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory factory(comm);
  factory.preference(pref);
  RCP<const Mesh> mesh = factory(0.0,0.0, 1.0,1.0, 21, 21, gm); 

  /* create a simple state and populate it */
  Amanzi::VerboseObject::hide_line_prefix = true;

  std::vector<std::string> component_names;
  component_names.push_back("Component 0");
  component_names.push_back("Component 1");

  Teuchos::ParameterList state_list = plist->sublist("State");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
  S->set_time(0.0);
  S->set_intermediate_time(0.0);
  S->set_initial_time(0.0);
  S->set_final_time(0.0);

  Transport_PK TPK(plist, S, "Transport", component_names);
  TPK.Setup();
  TPK.CreateDefaultState(mesh, 1);
  S->InitializeFields();
  S->InitializeEvaluators();

  /* modify the default state for the problem at hand */
  std::string passwd("state"); 
  S->GetFieldData("prev_saturation_liquid", passwd)->PutScalar(0.4);
  S->GetFieldData("saturation_liquid", passwd)->PutScalar(0.4);
  Teuchos::RCP<Epetra_MultiVector> flux = S->GetFieldData("darcy_flux", passwd)->ViewComponent("face", false);

  AmanziGeometry::Point velocity(0.1, 0.0);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    (*flux)[0][f] = velocity * normal;
  }

  Epetra_MultiVector& tcc = *S->GetFieldData("total_component_concentration", passwd)->ViewComponent("cell", false);

  *(S->GetScalarData("fluid_density", passwd)) = 1.0;

  /* initialize a transport process kernel */
  Amanzi::VerboseObject::hide_line_prefix = true;
  TPK.Initialize();

  /* advance the state */
  double dt0;
  dt0 = TPK.CalculateTransportDt();

  int i, k, iter = 0;
  double t_old(0.0), t_new(0.0), dt, T1(1.0);
  while (t_new < T1) {
    dt = std::min(TPK.CalculateTransportDt(), T1 - t_old);
    dt = std::min(dt, dt0);
    t_new = t_old + dt;

    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new);

    t_old = t_new;
    iter++;

    TPK.VV_CheckTracerBounds(tcc, 0, 0.0, 1.0, 1e-12);
    TPK.VV_CheckTracerBounds(tcc, 1, 0.0, 1.0, 1e-12);
  }

  // check symmetry of diffusion
  CHECK_CLOSE(tcc[1][199], tcc[1][241], 1e-12);
  CHECK_CLOSE(tcc[1][218], tcc[1][222], 1e-12);
 
  GMV::open_data_file(*mesh, (std::string)"transport.gmv");
  GMV::start_data();
  GMV::write_cell_data(tcc, 0, "Component_0");
  GMV::write_cell_data(tcc, 1, "Component_1");
  GMV::close_data_file();

  delete comm;
}




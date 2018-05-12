/*
  The transport component of the Amanzi code, serial unit tests.
  License: BSD
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

// TPLs
#include "Epetra_MpiComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "LeastSquare.hh"
#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "Point.hh"
#include "State.hh"

// Transport
#include "Transport_PK.hh"

double f_cubic(const Amanzi::AmanziGeometry::Point& x, double t) {
  if (x[0] < 1 + t) return 1.0;
  if (x[0] > 3 + t) return 0.0;
  double z = (x[0] - 1 - t) / 2;
  return 2 * z * z * z - 3 * z * z + 1;
}

double f_cubic_unit(const Amanzi::AmanziGeometry::Point& x, double t) {
  if (x[0] < 0.2 + t) return 1.0;
  if (x[0] > 0.5 + t) return 0.0;
  double z = (x[0] - 0.2 - t) / 0.3;
  return 2 * z * z * z - 3 * z * z + 1;
}


TEST(CONVERGENCE_ANALYSIS_DONOR) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "TEST: convergence analysis, donor scheme, orthogonal meshes" << std::endl;

  /* read parameter list */
  std::string xmlFileName = "test/transport_convergence.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  /* convergence estimate */
  std::vector<double> h;
  std::vector<double> L1error, L2error;

  for (int nx = 20; nx < 321; nx *= 2) {
    /* create a SIMPLE mesh framework */
    ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
        Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, &comm));

    MeshFactory meshfactory(&comm);
    meshfactory.preference(FrameworkPreference({Framework::MSTK, Framework::Simple}));
    RCP<const Mesh> mesh = meshfactory(0.0,0.0,0.0, 5.0,1.0,1.0, nx,2,2, gm);

    /* create a simple state and populate it */
    Amanzi::VerboseObject::hide_line_prefix = false;
    Amanzi::VerboseObject::global_default_level = Teuchos::VERB_NONE;

    std::vector<std::string> component_names;
    component_names.push_back("Component 0");

    Teuchos::ParameterList state_list = plist->sublist("state");
    RCP<State> S = rcp(new State(state_list));
    S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
    S->set_time(0.0);
    S->set_intermediate_time(0.0);

    Transport_PK TPK(plist, S, "transport", component_names);
    TPK.Setup(S.ptr());
    TPK.CreateDefaultState(mesh, 1);
    S->InitializeFields();
    S->InitializeEvaluators();

    /* modify the default state for the problem at hand */
    std::string passwd("state"); 
    Teuchos::RCP<Epetra_MultiVector> 
        flux = S->GetFieldData("darcy_flux", passwd)->ViewComponent("face", true);

    AmanziGeometry::Point velocity(1.0, 0.0, 0.0);
    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
    for (int f = 0; f < nfaces; f++) {
      const AmanziGeometry::Point& normal = mesh->face_normal(f);
      (*flux)[0][f] = velocity * normal;
    }

    Teuchos::RCP<Epetra_MultiVector> 
        tcc = S->GetFieldData("total_component_concentration", passwd)->ViewComponent("cell", true);

    int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
    for (int c = 0; c < ncells; c++) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
      (*tcc)[0][c] = f_cubic(xc, 0.0);
    }

    *(S->GetScalarData("fluid_density", passwd)) = 1.0;

    /* initialize a transport process kernel */
    TPK.Initialize(S.ptr());
    TPK.spatial_disc_order = TPK.temporal_disc_order = 1;
 
    /* advance the state */
    int iter = 0;
    double t_old(0.0), t_new, dt, T1(1.0);
    while (t_old < T1) {
      dt = std::min(TPK.StableTimeStep(), T1 - t_old);
      t_new = t_old + dt;

      TPK.AdvanceStep(t_old, t_new);
      TPK.CommitStep(t_old, t_new, S);

      t_old = t_new;
      iter++;
    }

    /* calculate L1 and L2 errors */
    double L1, L2;
    TPK.CalculateLpErrors(f_cubic, t_new, (*tcc)(0), &L1, &L2);
    if (MyPID == 0)
      printf("nx=%3d  L1 error=%7.5f  L2 error=%7.5f  dT=%7.4f\n", nx, L1, L2, T1 / iter);

    h.push_back(5.0 / nx);
    L1error.push_back(L1);
    L2error.push_back(L2);
  }

  double L1rate = Amanzi::Utils::bestLSfit(h, L1error);
  double L2rate = Amanzi::Utils::bestLSfit(h, L2error);
  if (MyPID == 0)
    printf("convergence rates: %8.2f %20.2f\n", L1rate, L2rate);

  CHECK_CLOSE(L1rate, 1.0, 0.1);
  CHECK_CLOSE(L2rate, 1.0, 0.1);
}


/* **************************************************************** */
TEST(CONVERGENCE_ANALYSIS_DONOR_SUBCYCLING) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0)
    std::cout << "\nTEST: convergence analysis, donor scheme, orthogonal meshes with subcycling" << std::endl;

  /* read parameter list */
  std::string xmlFileName = "test/transport_convergence.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  /* convergence estimate */
  std::vector<double> h;
  std::vector<double> L1error, L2error;

  for (int nx = 20; nx < 321; nx *= 2) {
    /* create a SIMPLE mesh framework */
    ParameterList region_list = plist->sublist("regions");
    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
        Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, &comm));

    MeshFactory meshfactory(&comm);
    meshfactory.preference(FrameworkPreference({Framework::MSTK, Framework::Simple}));
    RCP<const Mesh> mesh = meshfactory(0.0,0.0,0.0, 5.0,1.0,1.0, nx,2,2, gm);

    /* create a simple state and populate it */
    Amanzi::VerboseObject::hide_line_prefix = false;
    Amanzi::VerboseObject::global_default_level = Teuchos::VERB_NONE;

    std::vector<std::string> component_names;
    component_names.push_back("Component 0");

    Teuchos::ParameterList state_list = plist->sublist("state");
    RCP<State> S = rcp(new State(state_list));
    S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
    S->set_time(0.0);
    S->set_intermediate_time(0.0);

    Transport_PK TPK(plist, S, "transport", component_names);
    TPK.Setup(S.ptr());
    TPK.CreateDefaultState(mesh, 1);
    S->InitializeFields();
    S->InitializeEvaluators();

    /* modify the default state for the problem at hand */
    std::string passwd("state"); 
    Teuchos::RCP<Epetra_MultiVector> 
        flux = S->GetFieldData("darcy_flux", passwd)->ViewComponent("face", true);

    AmanziGeometry::Point velocity(1.0, 0.0, 0.0);
    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
    for (int f = 0; f < nfaces; f++) {
      const AmanziGeometry::Point& normal = mesh->face_normal(f);
      (*flux)[0][f] = velocity * normal;
    }

    Teuchos::RCP<Epetra_MultiVector> 
        tcc = S->GetFieldData("total_component_concentration", passwd)->ViewComponent("cell", true);

    int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
    for (int c = 0; c < ncells; c++) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
      (*tcc)[0][c] = f_cubic(xc, 0.0);
    }

    *(S->GetScalarData("fluid_density", passwd)) = 1.0;

    /* initialize a transport process kernel */
    TPK.Initialize(S.ptr());
    TPK.spatial_disc_order = TPK.temporal_disc_order = 1;
 
    /* advance the state */
    int ncycles = 0, iter = 0;
    double t_old(0.0), t_new, dt, T1(1.0);
    while (t_old < T1) {
      dt = std::min(TPK.StableTimeStep(), T1 - t_old);
      dt = dt * 7.7;
      t_new = t_old + dt;

      TPK.AdvanceStep(t_old, t_new);
      TPK.CommitStep(t_old ,t_new, S);

      t_old = t_new;
      ncycles += TPK.nsubcycles;
      iter++;
    }

    /* calculate L1 and L2 errors */
    double L1, L2;
    ncycles /= iter;
    TPK.CalculateLpErrors(f_cubic, t_new, (*tcc)(0), &L1, &L2);
    if (MyPID == 0)
      printf("nx=%3d  L1 error=%7.5f  L2 error=%7.5f  dT=%7.4f  (average # subcycles %d)\n", 
             nx, L1, L2, T1 / iter, ncycles);

    h.push_back(5.0 / nx);
    L1error.push_back(L1);
    L2error.push_back(L2);
  }

  double L1rate = Amanzi::Utils::bestLSfit(h, L1error);
  double L2rate = Amanzi::Utils::bestLSfit(h, L2error);
  if (MyPID == 0)
    printf("convergence rates: %8.2f %20.2f\n", L1rate, L2rate);

  CHECK_CLOSE(L1rate, 1.0, 0.1);
  CHECK_CLOSE(L2rate, 1.0, 0.1);
}


/* **************************************************************** */
TEST(CONVERGENCE_ANALYSIS_2ND) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) 
    std::cout << "\nTest: Convergence analysis, 2nd order scheme" << std::endl;

  /* read parameter list */
  std::string xmlFileName = "test/transport_convergence.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  /* create a SIMPLE mesh framework */
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, &comm));
 
  /* convergence estimate */
  double dt0;
  std::vector<double> h;
  std::vector<double> L1error, L2error;

  for (int nx = 20; nx < 161; nx *= 2) {
    MeshFactory meshfactory(&comm);
    meshfactory.preference(FrameworkPreference({Framework::MSTK, Framework::Simple}));
    RCP<const Mesh> mesh = meshfactory(0.0,0.0,0.0, 5.0,1.0,1.0, nx, 2, 1, gm); 

    /* create a simple state and populate it */
    Amanzi::VerboseObject::hide_line_prefix = false;
    Amanzi::VerboseObject::global_default_level = Teuchos::VERB_NONE;

    std::vector<std::string> component_names;
    component_names.push_back("Component 0");

    Teuchos::ParameterList state_list = plist->sublist("state");
    RCP<State> S = rcp(new State(state_list));
    S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

    Transport_PK TPK(plist, S, "transport", component_names);
    TPK.Setup(S.ptr());
    TPK.CreateDefaultState(mesh, 1);
    S->InitializeFields();
    S->InitializeEvaluators();

    /* modify the default state for the problem at hand */
    std::string passwd("state"); 
    Teuchos::RCP<Epetra_MultiVector> 
        flux = S->GetFieldData("darcy_flux", passwd)->ViewComponent("face", true);

    AmanziGeometry::Point velocity(1.0, 0.0, 0.0);
    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
    for (int f = 0; f < nfaces; f++) {
      const AmanziGeometry::Point& normal = mesh->face_normal(f);
      (*flux)[0][f] = velocity * normal;
    }

    Teuchos::RCP<Epetra_MultiVector> 
        tcc = S->GetFieldData("total_component_concentration", passwd)->ViewComponent("cell", true);

    int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
    for (int c = 0; c < ncells; c++) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
      (*tcc)[0][c] = f_cubic(xc, 0.0);
    }

    *(S->GetScalarData("fluid_density", passwd)) = 1.0;

    /* initialize a transport process kernel */
    TPK.Initialize(S.ptr());
    TPK.spatial_disc_order = TPK.temporal_disc_order = 2;
 
    /* advance the state */
    if (nx == 20) dt0 = TPK.StableTimeStep();
    else dt0 /= 2;

    int iter = 0;
    double t_old(0.0), t_new(0.0), dt, T1(2.0);
    while (t_new < T1) {
      dt = std::min(TPK.StableTimeStep(), T1 - t_old);
      dt = std::min(dt, dt0);
      t_new = t_old + dt;

      TPK.AdvanceStep(t_old, t_new);
      TPK.CommitStep(t_old, t_new, S);

      t_old = t_new;

      if (TPK.internal_tests_) {
        TPK.VV_CheckTracerBounds(*tcc, 0, 0.0, 1.0, 1e-12);
      }
      iter++;
    }
    //for (int k=0; k<nx; k++) std::cout << (*tcc)[0][k] << std::endl;

    double L1, L2;  // L1 and L2 errors
    TPK.CalculateLpErrors(f_cubic, t_new, (*tcc)(0), &L1, &L2);
    if (MyPID == 0)
      printf("nx=%3d  L1 error=%10.8f  L2 error=%10.8f  dT=%7.4f\n", nx, L1, L2, T1 / iter);

    h.push_back(5.0 / nx);
    L1error.push_back(L1);
    L2error.push_back(L2);
  }

  double L1rate = Amanzi::Utils::bestLSfit(h, L1error);
  double L2rate = Amanzi::Utils::bestLSfit(h, L2error);
  if (MyPID == 0)
    printf("convergence rates: %8.2f %20.2f\n", L1rate, L2rate);

  CHECK_CLOSE(2.0, L1rate, 0.4);
  CHECK_CLOSE(2.0, L2rate, 0.5);
}


/* **************************************************************** */
TEST(CONVERGENCE_ANALYSIS_DONOR_POLY) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0)
    std::cout << "\nTEST: convergence analysis, donor scheme, polygonal meshes" << std::endl;

  /* read parameter list */
  std::string xmlFileName = "test/transport_convergence_poly.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  /* convergence estimate */
  std::vector<double> h;
  std::vector<double> L1error, L2error;

  for (int loop = 0; loop < 3; loop++) {
    /* create a mesh framework */
    ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
        Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, &comm));

    MeshFactory meshfactory(&comm);
    meshfactory.preference(FrameworkPreference({Framework::MSTK, Framework::STKMESH}));
    RCP<const Mesh> mesh;
    if (loop == 0) {
      mesh = meshfactory("test/median15x16.exo", gm);
    } else if (loop == 1) {
      mesh = meshfactory("test/median32x33.exo", gm);
    } else if (loop == 2) {
      mesh = meshfactory("test/median63x64.exo", gm);
    }

    /* create a simple state and populate it */
    Amanzi::VerboseObject::hide_line_prefix = false;
    Amanzi::VerboseObject::global_default_level = Teuchos::VERB_NONE;

    std::vector<std::string> component_names;
    component_names.push_back("Component 0");

    Teuchos::ParameterList state_list = plist->sublist("state");
    RCP<State> S = rcp(new State(state_list));
    S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
    S->set_time(0.0);
    S->set_intermediate_time(0.0);

    Transport_PK TPK(plist, S, "transport", component_names);
    TPK.Setup(S.ptr());
    TPK.CreateDefaultState(mesh, 1);
    S->InitializeFields();
    S->InitializeEvaluators();

    /* modify the default state for the problem at hand */
    std::string passwd("state"); 
    Teuchos::RCP<Epetra_MultiVector> 
        flux = S->GetFieldData("darcy_flux", passwd)->ViewComponent("face", true);

    AmanziGeometry::Point velocity(1.0, 0.0);
    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
    for (int f = 0; f < nfaces; f++) {
      const AmanziGeometry::Point& normal = mesh->face_normal(f);
      (*flux)[0][f] = velocity * normal;
    }

    Teuchos::RCP<Epetra_MultiVector> 
        tcc = S->GetFieldData("total_component_concentration", passwd)->ViewComponent("cell", true);

    int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
    for (int c = 0; c < ncells; c++) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
      (*tcc)[0][c] = f_cubic_unit(xc, 0.0);
    }

    *(S->GetScalarData("fluid_density", passwd)) = 1.0;

    /* initialize a transport process kernel */
    TPK.Initialize(S.ptr());
    TPK.spatial_disc_order = TPK.temporal_disc_order = 1;
 
    /* advance the state */
    int iter = 0;
    double t_old(0.0), t_new, dt, T1(0.2);
    while (t_old < T1) {
      dt = std::min(TPK.StableTimeStep(), T1 - t_old);
      t_new = t_old + dt;

      TPK.AdvanceStep(t_old, t_new);
      TPK.CommitStep(t_old, t_new, S);

      t_old = t_new;
      iter++;
    }

    /* calculate L1 and L2 errors */
    double L1, L2;
    TPK.CalculateLpErrors(f_cubic_unit, t_new, (*tcc)(0), &L1, &L2);
    int nx = 16 * (loop + 1);
    if (MyPID == 0)
      printf("nx=%3d  L1 error=%7.5f  L2 error=%7.5f  dT=%7.4f\n", nx, L1, L2, T1 / iter);

    h.push_back(1.0 / nx);
    L1error.push_back(L1);
    L2error.push_back(L2);
  }

  double L1rate = Amanzi::Utils::bestLSfit(h, L1error);
  double L2rate = Amanzi::Utils::bestLSfit(h, L2error);
  if (MyPID == 0)
    printf("convergence rates: %5.2f %17.2f\n", L1rate, L2rate);

  CHECK_CLOSE(L1rate, 1.0, 0.1);
  CHECK_CLOSE(L2rate, 1.0, 0.1);
}


/* **************************************************************** */
TEST(CONVERGENCE_ANALYSIS_2ND_POLY) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0)
    std::cout << "\nTEST: convergence analysis, 2nd order scheme, polygonal meshes" << std::endl;

  /* read parameter list */
  std::string xmlFileName = "test/transport_convergence_poly.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  /* convergence estimate */
  std::vector<double> h;
  std::vector<double> L1error, L2error;

  for (int loop = 0; loop < 3; loop++) {
    /* create a mesh framework */
    ParameterList region_list = plist->sublist("regions");
    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
        Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, &comm));

    MeshFactory meshfactory(&comm);
    meshfactory.preference(FrameworkPreference({Framework::MSTK, Framework::STKMESH}));
    RCP<const Mesh> mesh;
    if (loop == 0) {
      mesh = meshfactory("test/median15x16.exo", gm);
    } else if (loop == 1) {
      mesh = meshfactory("test/median32x33.exo", gm);
    } else if (loop == 2) {
      mesh = meshfactory("test/median63x64.exo", gm);
    }

    /* create a simple state and populate it */
    Amanzi::VerboseObject::hide_line_prefix = false;
    Amanzi::VerboseObject::global_default_level = Teuchos::VERB_NONE;

    std::vector<std::string> component_names;
    component_names.push_back("Component 0");

    Teuchos::ParameterList state_list = plist->sublist("state");
    RCP<State> S = rcp(new State(state_list));
    S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
    S->set_time(0.0);
    S->set_intermediate_time(0.0);

    Transport_PK TPK(plist, S, "transport", component_names);
    TPK.Setup(S.ptr());
    TPK.CreateDefaultState(mesh, 1);
    S->InitializeFields();
    S->InitializeEvaluators();

    /* modify the default state for the problem at hand */
    std::string passwd("state"); 
    Teuchos::RCP<Epetra_MultiVector> 
        flux = S->GetFieldData("darcy_flux", passwd)->ViewComponent("face", true);

    AmanziGeometry::Point velocity(1.0, 0.0);
    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
    for (int f = 0; f < nfaces; f++) {
      const AmanziGeometry::Point& normal = mesh->face_normal(f);
      (*flux)[0][f] = velocity * normal;
    }

    Teuchos::RCP<Epetra_MultiVector> 
        tcc = S->GetFieldData("total_component_concentration", passwd)->ViewComponent("cell", true);

    int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
    for (int c = 0; c < ncells; c++) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
      (*tcc)[0][c] = f_cubic_unit(xc, 0.0);
    }

    *(S->GetScalarData("fluid_density", passwd)) = 1.0;

    /* initialize a transport process kernel */
    TPK.Initialize(S.ptr());
    TPK.spatial_disc_order = TPK.temporal_disc_order = 2;
 
    /* advance the state */
    int iter = 0;
    double t_old(0.0), t_new, dt, T1(0.2);
    while (t_old < T1) {
      dt = std::min(TPK.StableTimeStep(), T1 - t_old);
      t_new = t_old + dt;

      TPK.AdvanceStep(t_old, t_new);
      TPK.CommitStep(t_old, t_new, S);

      t_old = t_new;
      iter++;
    }

    /* calculate L1 and L2 errors */
    double L1, L2;
    TPK.CalculateLpErrors(f_cubic_unit, t_new, (*tcc)(0), &L1, &L2);
    int nx = 16 * (loop + 1);
    if (MyPID == 0)
      printf("nx=%3d  L1 error=%7.5f  L2 error=%7.5f  dT=%7.4f\n", nx, L1, L2, T1 / iter);

    h.push_back(1.0 / nx);
    L1error.push_back(L1);
    L2error.push_back(L2);
  }

  double L1rate = Amanzi::Utils::bestLSfit(h, L1error);
  double L2rate = Amanzi::Utils::bestLSfit(h, L2error);
  if (MyPID == 0)
    printf("convergence rates: %5.2f %17.2f\n", L1rate, L2rate);

  CHECK_CLOSE(2.0, L1rate, 0.4);
  CHECK_CLOSE(2.0, L2rate, 0.5);
}


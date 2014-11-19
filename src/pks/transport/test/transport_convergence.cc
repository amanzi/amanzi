/*
  The transport component of the Amanzi code, serial unit tests.
  License: BSD
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Epetra_SerialComm.h"

#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "Point.hh"

#include "State.hh"
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

  std::cout << "TEST: convergence analysis, donor scheme, orthogonal meshes" << std::endl;
  Epetra_SerialComm* comm = new Epetra_SerialComm();

  /* read parameter list */
  std::string xmlFileName = "test/transport_convergence.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  /* convergence estimate */
  std::vector<double> h;
  std::vector<double> L1error, L2error;

  for (int nx = 20; nx < 321; nx *= 2) {
    /* create a SIMPLE mesh framework */
    ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
    GeometricModelPtr gm = new GeometricModel(3, region_list, (Epetra_MpiComm *)comm);

    FrameworkPreference pref;
    pref.clear();
    pref.push_back(Simple);
    
    MeshFactory meshfactory((Epetra_MpiComm *)comm);
    meshfactory.preference(pref);
    RCP<const Mesh> mesh = meshfactory(0.0,0.0,0.0, 5.0,1.0,1.0, nx,2,2, gm);

    /* create a simple state and populate it */
    Amanzi::VerboseObject::hide_line_prefix = false;
    Amanzi::VerboseObject::global_default_level = Teuchos::VERB_NONE;

    std::vector<std::string> component_names;
    component_names.push_back("Component 0");

    RCP<State> S = rcp(new State());
    S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
    S->set_time(0.0);
    S->set_intermediate_time(0.0);

    Transport_PK TPK(plist, S, component_names);
    TPK.CreateDefaultState(mesh, 1);

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
      (*tcc)[0][c] = f_cubic(xc, 0.0);
    }

    S->GetFieldData("porosity", passwd)->PutScalar(1.0);
    *(S->GetScalarData("fluid_density", passwd)) = 1.0;

    /* initialize a transport process kernel */
    TPK.Initialize(S.ptr());
    TPK.spatial_disc_order = TPK.temporal_disc_order = 1;
    if (nx == 20) TPK.PrintStatistics();
 
    /* advance the state */
    int iter = 0;
    double dummy_dT, T = 0.0, T1 = 1.0;
    while (T < T1) {
      double dT = std::min(TPK.CalculateTransportDt(), T1 - T);
      TPK.Advance(dT, dummy_dT);
      TPK.CommitState(dT, S.ptr());
      T += dT;
      iter++;
    }

    /* calculate L1 and L2 errors */
    double L1, L2;
    TPK.CalculateLpErrors(f_cubic, T, (*tcc)(0), &L1, &L2);
    printf("nx=%3d  L1 error=%7.5f  L2 error=%7.5f  dT=%7.4f\n", nx, L1, L2, T1 / iter);

    h.push_back(5.0 / nx);
    L1error.push_back(L1);
    L2error.push_back(L2);

    delete gm;
  }

  double L1rate = Amanzi::Transport::bestLSfit(h, L1error);
  double L2rate = Amanzi::Transport::bestLSfit(h, L2error);
  printf("convergence rates: %8.2f %20.2f\n", L1rate, L2rate);

  CHECK_CLOSE(L1rate, 1.0, 0.1);
  CHECK_CLOSE(L2rate, 1.0, 0.1);

  delete comm;
}


/* **************************************************************** */
TEST(CONVERGENCE_ANALYSIS_DONOR_SUBCYCLING) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "\nTEST: convergence analysis, donor scheme, orthogonal meshes with subcycling" << std::endl;
  Epetra_SerialComm* comm = new Epetra_SerialComm();

  /* read parameter list */
  std::string xmlFileName = "test/transport_convergence.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  /* convergence estimate */
  std::vector<double> h;
  std::vector<double> L1error, L2error;

  for (int nx = 20; nx < 321; nx *= 2) {
    /* create a SIMPLE mesh framework */
    ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
    GeometricModelPtr gm = new GeometricModel(3, region_list, (Epetra_MpiComm *)comm);

    FrameworkPreference pref;
    pref.clear();
    pref.push_back(Simple);
    
    MeshFactory meshfactory((Epetra_MpiComm *)comm);
    meshfactory.preference(pref);
    RCP<const Mesh> mesh = meshfactory(0.0,0.0,0.0, 5.0,1.0,1.0, nx,2,2, gm);

    /* create a simple state and populate it */
    Amanzi::VerboseObject::hide_line_prefix = false;
    Amanzi::VerboseObject::global_default_level = Teuchos::VERB_NONE;

    std::vector<std::string> component_names;
    component_names.push_back("Component 0");

    RCP<State> S = rcp(new State());
    S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
    S->set_time(0.0);
    S->set_intermediate_time(0.0);

    Transport_PK TPK(plist, S, component_names);
    TPK.CreateDefaultState(mesh, 1);

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
      (*tcc)[0][c] = f_cubic(xc, 0.0);
    }

    S->GetFieldData("porosity", passwd)->PutScalar(1.0);
    *(S->GetScalarData("fluid_density", passwd)) = 1.0;

    /* initialize a transport process kernel */
    TPK.Initialize(S.ptr());
    TPK.spatial_disc_order = TPK.temporal_disc_order = 1;
    if (nx == 20) TPK.PrintStatistics();
 
    /* advance the state */
    int ncycles = 0, iter = 0;
    double dummy_dT, T = 0.0, T1 = 1.0;
    while (T < T1) {
      double dT = std::min(TPK.CalculateTransportDt(), T1 - T);
      dT = dT * 7.7;

      TPK.Advance(dT, dummy_dT);
      TPK.CommitState(dT, S.ptr());
      T += dT;

      ncycles += TPK.nsubcycles;
      iter++;
    }

    /* calculate L1 and L2 errors */
    double L1, L2;
    ncycles /= iter;
    TPK.CalculateLpErrors(f_cubic, T, (*tcc)(0), &L1, &L2);
    printf("nx=%3d  L1 error=%7.5f  L2 error=%7.5f  dT=%7.4f  (average # subcycles %d)\n", 
           nx, L1, L2, T1 / iter, ncycles);

    h.push_back(5.0 / nx);
    L1error.push_back(L1);
    L2error.push_back(L2);

    delete gm;
  }

  double L1rate = Amanzi::Transport::bestLSfit(h, L1error);
  double L2rate = Amanzi::Transport::bestLSfit(h, L2error);
  printf("convergence rates: %8.2f %20.2f\n", L1rate, L2rate);

  CHECK_CLOSE(L1rate, 1.0, 0.1);
  CHECK_CLOSE(L2rate, 1.0, 0.1);

  delete comm;
}


/* **************************************************************** */
TEST(CONVERGENCE_ANALYSIS_2ND) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "\nTest: Convergence analysis, 2nd order scheme" << std::endl;
  Epetra_SerialComm  *comm = new Epetra_SerialComm();

  /* read parameter list */
  std::string xmlFileName = "test/transport_convergence.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  /* create a SIMPLE mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, (Epetra_MpiComm *)comm);
 
  /* convergence estimate */
  std::vector<double> h;
  std::vector<double> L1error, L2error;

  for (int nx = 20; nx < 161; nx *= 2) {
    FrameworkPreference pref;
    pref.clear();
    pref.push_back(Simple);
    
    MeshFactory meshfactory((Epetra_MpiComm *)comm);
    meshfactory.preference(pref);
    RCP<const Mesh> mesh = meshfactory(0.0,0.0,0.0, 5.0,1.0,1.0, nx, 2, 1, gm); 

    /* create a simple state and populate it */
    Amanzi::VerboseObject::hide_line_prefix = false;
    Amanzi::VerboseObject::global_default_level = Teuchos::VERB_NONE;

    std::vector<std::string> component_names;
    component_names.push_back("Component 0");

    RCP<State> S = rcp(new State());
    S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

    Transport_PK TPK(plist, S, component_names);
    TPK.CreateDefaultState(mesh, 1);

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
      (*tcc)[0][c] = f_cubic(xc, 0.0);
    }

    S->GetFieldData("porosity", passwd)->PutScalar(1.0);
    *(S->GetScalarData("fluid_density", passwd)) = 1.0;

    /* initialize a transport process kernel */
    TPK.Initialize(S.ptr());
    TPK.spatial_disc_order = TPK.temporal_disc_order = 2;
    if (nx == 20) TPK.PrintStatistics();
 
    /* advance the state */
    double dummy_dT, dT, dT0;
    if (nx == 20) dT0 = TPK.CalculateTransportDt();
    else dT0 /= 2;

    int iter = 0;
    double T = 0.0, T1 = 2.0;
    while (T < T1) {
      dT = std::min(TPK.CalculateTransportDt(), T1 - T);
      dT = std::min(dT, dT0);

      TPK.Advance(dT, dummy_dT);
      TPK.CommitState(dT, S.ptr());
      T += dT;
      if (TPK.internal_tests) {
        TPK.VV_CheckTracerBounds(*tcc, 0, 0.0, 1.0, 1e-12);
      }
      iter++;
    }
    //for (int k=0; k<nx; k++) std::cout << (*tcc)[0][k] << std::endl;

    double L1, L2;  // L1 and L2 errors
    TPK.CalculateLpErrors(f_cubic, T, (*tcc)(0), &L1, &L2);
    printf("nx=%3d  L1 error=%10.8f  L2 error=%10.8f  dT=%7.4f\n", nx, L1, L2, T1 / iter);

    h.push_back(5.0 / nx);
    L1error.push_back(L1);
    L2error.push_back(L2);
  }

  double L1rate = Amanzi::Transport::bestLSfit(h, L1error);
  double L2rate = Amanzi::Transport::bestLSfit(h, L2error);
  printf("convergence rates: %8.2f %20.2f\n", L1rate, L2rate);

  CHECK_CLOSE(2.0, L1rate, 0.3);
  CHECK_CLOSE(2.0, L2rate, 0.5);

  delete gm;
  delete comm;
}


/* **************************************************************** */
TEST(CONVERGENCE_ANALYSIS_DONOR_POLY) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "\nTEST: convergence analysis, donor scheme, polygonal meshes" << std::endl;
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  /* read parameter list */
  std::string xmlFileName = "test/transport_convergence_poly.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  /* convergence estimate */
  std::vector<double> h;
  std::vector<double> L1error, L2error;

  for (int loop = 0; loop < 3; loop++) {
    /* create a mesh framework */
    ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
    GeometricModelPtr gm = new GeometricModel(2, region_list, comm);

    FrameworkPreference pref;
    pref.clear();
    pref.push_back(MSTK);
    pref.push_back(STKMESH);
    
    MeshFactory meshfactory(comm);
    meshfactory.preference(pref);
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

    RCP<State> S = rcp(new State());
    S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
    S->set_time(0.0);
    S->set_intermediate_time(0.0);

    Transport_PK TPK(plist, S, component_names);
    TPK.CreateDefaultState(mesh, 1);

    /* modify the default state for the problem at hand */
    std::string passwd("state"); 
    Teuchos::RCP<Epetra_MultiVector> 
        flux = S->GetFieldData("darcy_flux", passwd)->ViewComponent("face", false);

    AmanziGeometry::Point velocity(1.0, 0.0);
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
      (*tcc)[0][c] = f_cubic_unit(xc, 0.0);
    }

    S->GetFieldData("porosity", passwd)->PutScalar(1.0);
    *(S->GetScalarData("fluid_density", passwd)) = 1.0;

    /* initialize a transport process kernel */
    TPK.Initialize(S.ptr());
    TPK.spatial_disc_order = TPK.temporal_disc_order = 1;
    if (loop == 0) TPK.PrintStatistics();
 
    /* advance the state */
    int iter = 0;
    double dummy_dT, T = 0.0, T1 = 0.2;
    while (T < T1) {
      double dT = std::min(TPK.CalculateTransportDt(), T1 - T);
      TPK.Advance(dT, dummy_dT);
      TPK.CommitState(dT, S.ptr());
      T += dT;
      iter++;
    }

    /* calculate L1 and L2 errors */
    double L1, L2;
    TPK.CalculateLpErrors(f_cubic_unit, T, (*tcc)(0), &L1, &L2);
    int nx = 16 * (loop + 1);
    printf("nx=%3d  L1 error=%7.5f  L2 error=%7.5f  dT=%7.4f\n", nx, L1, L2, T1 / iter);

    h.push_back(1.0 / nx);
    L1error.push_back(L1);
    L2error.push_back(L2);

    delete gm;
  }

  double L1rate = Amanzi::Transport::bestLSfit(h, L1error);
  double L2rate = Amanzi::Transport::bestLSfit(h, L2error);
  printf("convergence rates: %5.2f %17.2f\n", L1rate, L2rate);

  CHECK_CLOSE(L1rate, 1.0, 0.1);
  CHECK_CLOSE(L2rate, 1.0, 0.1);

  delete comm;
}


/* **************************************************************** */
TEST(CONVERGENCE_ANALYSIS_2ND_POLY) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "\nTEST: convergence analysis, 2nd order scheme, polygonal meshes" << std::endl;
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  /* read parameter list */
  std::string xmlFileName = "test/transport_convergence_poly.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  /* convergence estimate */
  std::vector<double> h;
  std::vector<double> L1error, L2error;

  for (int loop = 0; loop < 3; loop++) {
    /* create a mesh framework */
    ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
    GeometricModelPtr gm = new GeometricModel(2, region_list, comm);

    FrameworkPreference pref;
    pref.clear();
    pref.push_back(MSTK);
    pref.push_back(STKMESH);
    
    MeshFactory meshfactory(comm);
    meshfactory.preference(pref);
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

    RCP<State> S = rcp(new State());
    S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
    S->set_time(0.0);
    S->set_intermediate_time(0.0);

    Transport_PK TPK(plist, S, component_names);
    TPK.CreateDefaultState(mesh, 1);

    /* modify the default state for the problem at hand */
    std::string passwd("state"); 
    Teuchos::RCP<Epetra_MultiVector> 
        flux = S->GetFieldData("darcy_flux", passwd)->ViewComponent("face", false);

    AmanziGeometry::Point velocity(1.0, 0.0);
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
      (*tcc)[0][c] = f_cubic_unit(xc, 0.0);
    }

    S->GetFieldData("porosity", passwd)->PutScalar(1.0);
    *(S->GetScalarData("fluid_density", passwd)) = 1.0;

    /* initialize a transport process kernel */
    TPK.Initialize(S.ptr());
    TPK.spatial_disc_order = TPK.temporal_disc_order = 2;
    if (loop == 0) TPK.PrintStatistics();
 
    /* advance the state */
    int iter = 0;
    double dummy_dT, T = 0.0, T1 = 0.2;
    while (T < T1) {
      double dT = std::min(TPK.CalculateTransportDt(), T1 - T);
      TPK.Advance(dT, dummy_dT);
      TPK.CommitState(dT, S.ptr());
      T += dT;
      iter++;
    }

    /* calculate L1 and L2 errors */
    double L1, L2;
    TPK.CalculateLpErrors(f_cubic_unit, T, (*tcc)(0), &L1, &L2);
    int nx = 16 * (loop + 1);
    printf("nx=%3d  L1 error=%7.5f  L2 error=%7.5f  dT=%7.4f\n", nx, L1, L2, T1 / iter);

    h.push_back(1.0 / nx);
    L1error.push_back(L1);
    L2error.push_back(L2);

    delete gm;
  }

  double L1rate = Amanzi::Transport::bestLSfit(h, L1error);
  double L2rate = Amanzi::Transport::bestLSfit(h, L2error);
  printf("convergence rates: %5.2f %17.2f\n", L1rate, L2rate);

  CHECK_CLOSE(2.0, L1rate, 0.3);
  CHECK_CLOSE(2.0, L2rate, 0.5);

  delete comm;
}


/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <iostream>
#include <string>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_MpiComm.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Darcy_PK.hh"

#include "exceptions.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::AmanziFlow;

class DarcyProblem {
 public:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh;

  Teuchos::RCP<State> S;
  AmanziFlow::Darcy_PK* DPK;
  Teuchos::ParameterList dp_list;
  std::string passwd;

  Epetra_MpiComm* comm;
  int MyPID;

  DarcyProblem() {
    passwd = "state"; 
    comm = new Epetra_MpiComm(MPI_COMM_WORLD);
    MyPID = comm->MyPID();    
  };

  ~DarcyProblem() {
    delete DPK;
    delete comm;
  }

  void Init(const string xmlFileName, const char* meshExodus) {
    Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
    Teuchos::ParameterList plist = xmlreader.getParameters();

    /* create a MSTK mesh framework */
    Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
    GeometricModelPtr gm = new GeometricModel(3, region_list, comm);

    FrameworkPreference pref;
    pref.clear();
    pref.push_back(MSTK);

    MeshFactory meshfactory(comm);
    meshfactory.preference(pref);
    mesh = meshfactory(meshExodus, gm);

    /* create Darcy process kernel */
    Teuchos::ParameterList state_list = plist.sublist("State");
    S = Teuchos::rcp(new State(state_list));
    S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));
    S->set_time(0.0);

    DPK = new Darcy_PK(plist, S);
    S->Setup();
    S->InitializeFields();
    DPK->CreateDefaultState(mesh);

    S->CheckAllFieldsInitialized();

    Teuchos::ParameterList& flow_list = plist.get<Teuchos::ParameterList>("Flow");
    dp_list = flow_list.get<Teuchos::ParameterList>("Darcy Problem");
 }

  void createBClist(const char* type, const char* bc_x, 
                    Teuchos::Array<std::string>& regions, double value) {
    std::string func_list_name;
    if (!strcmp(type, "pressure")) {
      func_list_name = "boundary pressure";
    } else if (!strcmp(type, "static head")) {
      func_list_name = "water table elevation";
    } else if (!strcmp(type, "mass flux")) {
      func_list_name = "outward mass flux";
    }
    Teuchos::ParameterList& bc_list = dp_list.get<Teuchos::ParameterList>("boundary conditions");
    Teuchos::ParameterList& type_list = bc_list.get<Teuchos::ParameterList>(type);

    Teuchos::ParameterList& bc_sublist = type_list.sublist(bc_x);
    bc_sublist.set("regions", regions);

    Teuchos::ParameterList& bc_sublist_named = bc_sublist.sublist(func_list_name);
    Teuchos::ParameterList& function_list = bc_sublist_named.sublist("function-constant");
    function_list.set("value", value);
  }

  double calculatePressureCellError(double p0, AmanziGeometry::Point& pressure_gradient) {
    Epetra_MultiVector& pressure = *S->GetFieldData("pressure", passwd)->ViewComponent("cell", false);

    double error_L2 = 0.0;
    int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (int c = 0; c < ncells; c++) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
      double pressure_exact = p0 + pressure_gradient * xc;
      //      if (MyPID==0) cout << c << " " << pressure[0][c] << " exact=" <<  pressure_exact << endl;
      error_L2 += std::pow(pressure[0][c] - pressure_exact, 2.0);
    }
    return sqrt(error_L2);
  }

  double calculatePressureFaceError(double p0, AmanziGeometry::Point& pressure_gradient) {
    Epetra_MultiVector& lambda = *S->GetFieldData("pressure", passwd)->ViewComponent("face", false);

    double error_L2 = 0.0;
    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
    for (int f = 0; f < nfaces; f++) {
      const AmanziGeometry::Point& xf = mesh->face_centroid(f);
      double pressure_exact = p0 + pressure_gradient * xf;
      //      cout << f << " " << lambda[0][f] << " exact=" << pressure_exact << endl;
      error_L2 += std::pow(lambda[0][f] - pressure_exact, 2.0);
    }
    return sqrt(error_L2);
  }

  double calculateDarcyFluxError(AmanziGeometry::Point& velocity_exact) {
    Epetra_MultiVector& flux = *S->GetFieldData("darcy_flux", passwd)->ViewComponent("face", false);

    double error_L2 = 0.0;
    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
    for (int f = 0; f < nfaces; f++) {
      const AmanziGeometry::Point& normal = mesh->face_normal(f);
      // cout << f << " " << flux[0][f] << " exact=" << velocity_exact * normal << endl;
      error_L2 += std::pow(flux[0][f] - velocity_exact * normal, 2.0);
    }
    return sqrt(error_L2);
  }

  double calculateDarcyDivergenceError() {
    Teuchos::RCP<CompositeVector> cv = S->GetFieldData("darcy_flux", passwd);
    cv->ScatterMasterToGhosted("face");
    Epetra_MultiVector& flux = *cv->ViewComponent("face", true);

    double error_L2 = 0.0;
    int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

    for (int c = 0; c < ncells_owned; c++) {
      AmanziMesh::Entity_ID_List faces;
      std::vector<int> dirs;

      mesh->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();

      double div = 0.0;
      for (int i = 0; i < nfaces; i++) {
        int f = faces[i];
        div += flux[0][f] * dirs[i];
      }
      error_L2 += div*div / mesh->cell_volume(c);
    }
    return sqrt(error_L2);
  }
};


SUITE(Darcy_PK) {
/* ******************************************************************
* Testing the mesh of hexahedra.
****************************************************************** */
  TEST_FIXTURE(DarcyProblem, DirichletDirichlet) {
    Init("test/flow_darcy_misc.xml", "test/hexes.exo");

    if (MyPID == 0) std::cout <<"Darcy PK on hexes: Dirichlet-Dirichlet" << std::endl;

    double p0 = 1.0;
    AmanziGeometry::Point pressure_gradient(0.0, 0.0, -1.0);
    AmanziGeometry::Point velocity(3);
    velocity = -(pressure_gradient - DPK->rho() * DPK->gravity()) / DPK->mu();

    Teuchos::Array<std::string> regions(1);  // modify boundary conditions
    regions[0] = string("Top side");
    createBClist("pressure", "BC 1", regions, 0.0);

    regions[0] = string("Bottom side");
    createBClist("pressure", "BC 2", regions, 1.0);

    DPK->ResetParameterList(dp_list);
    DPK->InitPK();  // setup the problem
    DPK->InitSteadyState(0.0, 1.0);
    DPK->AdvanceToSteadyState(0.0, 1.0);
    DPK->CommitState(S);

    double errorP = calculatePressureCellError(p0, pressure_gradient);  // error checks
    CHECK(errorP < 1.0e-8);
    double errorL = calculatePressureFaceError(p0, pressure_gradient);
    CHECK(errorL < 1.0e-8);
    double errorU = calculateDarcyFluxError(velocity);
    CHECK(errorU < 1.0e-8);
    double errorDiv = calculateDarcyDivergenceError();
    if (MyPID == 0)
        std::cout << "Error: " << errorP << " " << errorL << " " << errorU << " " << errorDiv << std::endl;
  }


  TEST_FIXTURE(DarcyProblem, DirichletNeumann) {
    Init("test/flow_darcy_misc.xml", "test/hexes.exo");
    if (MyPID == 0) std::cout <<"Darcy PK on hexes: Dirichlet-Neumann" << std::endl;

    double p0 = 1.0;
    AmanziGeometry::Point pressure_gradient(0.0, 0.0, -1.0);
    AmanziGeometry::Point velocity(3);
    velocity = -(pressure_gradient - DPK->rho() * DPK->gravity()) / DPK->mu();
    double u0 = DPK->rho() * velocity * AmanziGeometry::Point(0.0, 0.0, 1.0);

    Teuchos::Array<std::string> regions(1);  // modify boundary conditions
    regions[0] = string("Top side");
    createBClist("mass flux", "BC 1", regions, u0);

    regions[0] = string("Bottom side");
    createBClist("pressure", "BC 2", regions, 1.0);
    DPK->ResetParameterList(dp_list);

    DPK->InitPK();  // setup the problem
    DPK->InitSteadyState(0.0, 1.0);
    DPK->AdvanceToSteadyState(0.0, 1.0);
    DPK->CommitState(S);

    double error = calculatePressureCellError(p0, pressure_gradient);
    CHECK(error < 1.0e-8);
    error = calculatePressureFaceError(p0, pressure_gradient);
    CHECK(error < 1.0e-8);
    error = calculateDarcyFluxError(velocity);
    CHECK(error < 1.0e-8);
  }

  TEST_FIXTURE(DarcyProblem, StaticHeadDirichlet) {
    Init("test/flow_darcy_misc.xml", "test/hexes.exo");
    if (MyPID == 0) std::cout <<"Darcy PK on hexes: StaticHead-Dirichlet" << std::endl;

    double p0 = 2.0;
    AmanziGeometry::Point pressure_gradient(0.0, 0.0, -1.0);
    AmanziGeometry::Point velocity(3);
    velocity = -(pressure_gradient - DPK->rho() * DPK->gravity()) / DPK->mu();

    Teuchos::Array<std::string> regions(1);  // modify boundary conditions
    regions[0] = string("Top side");
    createBClist("pressure", "BC 1", regions, 1.0);

    regions[0] = string("Bottom side");
    createBClist("static head", "BC 2", regions, 0.25);
    DPK->ResetParameterList(dp_list);

    DPK->InitPK();  // setup the problem
    DPK->InitSteadyState(0.0, 1.0);
    DPK->AdvanceToSteadyState(0.0, 1.0);
    DPK->CommitState(S);

    double error = calculatePressureCellError(p0, pressure_gradient);  // error checks
    CHECK(error < 1.0e-8);
    error = calculatePressureFaceError(p0, pressure_gradient);
    CHECK(error < 1.0e-8);
    error = calculateDarcyFluxError(velocity);
    CHECK(error < 1.0e-8);
  }

/* ******************************************************************
* Testing the mesh of prisms.
****************************************************************** */
  TEST_FIXTURE(DarcyProblem, DDprisms) {
    Init("test/flow_darcy_misc.xml", "test/prisms.exo");
    if (MyPID == 0) std::cout << "\nDarcy PK on prisms: Dirichlet-Dirichlet" << std::endl;

    double p0 = 1.0;
    AmanziGeometry::Point pressure_gradient(0.0, 0.0, -1.0);
    AmanziGeometry::Point velocity(3);
    velocity = -(pressure_gradient - DPK->rho() * DPK->gravity()) / DPK->mu();

    Teuchos::Array<std::string> regions(1);  // modify boundary conditions
    regions[0] = string("Top side");
    createBClist("pressure", "BC 1", regions, 0.0);

    regions[0] = string("Bottom side");
    createBClist("pressure", "BC 2", regions, 1.0);
    DPK->ResetParameterList(dp_list);

    DPK->InitPK();  // setup the problem
    DPK->InitSteadyState(0.0, 1.0);
    DPK->AdvanceToSteadyState(0.0, 1.0);
    DPK->CommitState(S);

    double errorP = calculatePressureCellError(p0, pressure_gradient);  // error checks
    CHECK(errorP < 1.0e-8);
    double errorL = calculatePressureFaceError(p0, pressure_gradient);
    CHECK(errorL < 1.0e-8);
    double errorU = calculateDarcyFluxError(velocity);
    CHECK(errorU < 1.0e-8);
    double errorDiv = calculateDarcyDivergenceError();
    if (MyPID == 0)
        std::cout << "Error: " << errorP << " " << errorL << " " << errorU << " " << errorDiv << std::endl;
  }

/* ******************************************************************
* Testing the mesh of tetrahedra.
****************************************************************** */
  TEST_FIXTURE(DarcyProblem, DNtetrahedra) {
    Init("test/flow_darcy_misc.xml", "test/tetrahedra.exo");
    if (MyPID == 0) std::cout << "\nDarcy PK on tets: Dirichlet-Neumann" << std::endl;

    double p0 = 1.0;
    AmanziGeometry::Point pressure_gradient(0.0, 0.0, -1.0);
    AmanziGeometry::Point velocity(3);
    velocity = -(pressure_gradient - DPK->rho() * DPK->gravity()) / DPK->mu();
    double u0 = DPK->rho() * velocity * AmanziGeometry::Point(0.0, 0.0, 1.0);

    Teuchos::Array<std::string> regions(1);  // modify boundary conditions
    regions[0] = string("Top side");
    createBClist("mass flux", "BC 1", regions, u0);

    regions[0] = string("Bottom side");
    createBClist("pressure", "BC 2", regions, 1.0);
    DPK->ResetParameterList(dp_list);

    DPK->InitPK();  // setup the problem
    DPK->InitSteadyState(0.0, 1.0);
    DPK->AdvanceToSteadyState(0.0, 1.0);
    DPK->CommitState(S);

    double errorP = calculatePressureCellError(p0, pressure_gradient);
    CHECK(errorP < 1.0e-8);
    double errorL = calculatePressureFaceError(p0, pressure_gradient);
    CHECK(errorL < 1.0e-8);
    double errorU = calculateDarcyFluxError(velocity);
    CHECK(errorU < 1.0e-8);
    if (MyPID == 0)
        std::cout << "Error: " << errorP << " " << errorL << " " << errorU << std::endl;
  }

/* ******************************************************************
* Testing the mesh of prisms and hexahedra.
****************************************************************** */
  TEST_FIXTURE(DarcyProblem, DDmixed) {
    Init("test/flow_darcy_misc.xml", "test/mixed.exo");
    if (MyPID == 0) std::cout << "\nDarcy PK on mixed mesh: Dirichlet-Dirichlet" << std::endl;

    double p0 = 1.0;
    AmanziGeometry::Point pressure_gradient(0.0, 0.0, -1.0);
    AmanziGeometry::Point velocity(3);
    velocity = -(pressure_gradient - DPK->rho() * DPK->gravity()) / DPK->mu();

    Teuchos::Array<std::string> regions(1);  // modify boundary conditions
    regions[0] = string("Top side");
    createBClist("pressure", "BC 1", regions, 0.0);

    regions[0] = string("Bottom side");
    createBClist("pressure", "BC 2", regions, 1.0);
    DPK->ResetParameterList(dp_list);

    DPK->InitPK();  // setup the problem
    DPK->InitSteadyState(0.0, 1.0);
    DPK->AdvanceToSteadyState(0.0, 1.0);
    DPK->CommitState(S);

    double errorP = calculatePressureCellError(p0, pressure_gradient);  // error checks
    CHECK(errorP < 1.0e-8);
    double errorL = calculatePressureFaceError(p0, pressure_gradient);
    CHECK(errorL < 1.0e-8);
    double errorU = calculateDarcyFluxError(velocity);
    CHECK(errorU < 1.0e-8);
    double errorDiv = calculateDarcyDivergenceError();
    if (MyPID == 0)
        std::cout << "Error: " << errorP << " " << errorL << " " << errorU << " " << errorDiv << std::endl;
  }

}




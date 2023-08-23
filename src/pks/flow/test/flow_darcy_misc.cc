/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Epetra_SerialComm.h"
#include "AmanziComm.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "exceptions.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "MeshException.hh"

// Flow
#include "Darcy_PK.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Flow;

class DarcyProblem {
 public:
  Teuchos::RCP<Teuchos::ParameterList> plist;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh;

  Teuchos::RCP<State> S;
  Flow::Darcy_PK* DPK;
  std::string passwd;
  double mu, rho;
  AmanziGeometry::Point gravity;

  Comm_ptr_type comm;
  int MyPID;

  Framework frameworks[2];
  std::vector<std::string> framework_name;

  DarcyProblem()
  {
    passwd = "";
    comm = Amanzi::getDefaultComm();
    MyPID = comm->MyPID();

    frameworks[0] = Framework::MSTK;
    frameworks[1] = Framework::STK;
    framework_name.push_back("MSTK");
    framework_name.push_back("stk:mesh");
  };

  ~DarcyProblem() { delete DPK; }

  int Init(const std::string xmlFileName, const char* meshExodus, const Framework& framework)
  {
    if (framework == Framework::STK && comm->NumProc() > 1) return 1;

    // create a MSTK mesh framework
    plist = Teuchos::getParametersFromXmlFile(xmlFileName);
    Teuchos::ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, regions_list, *comm));

    Preference pref({ framework });

    MeshFactory meshfactory(comm, gm);
    try {
      meshfactory.set_preference(pref);
    } catch (Message& e) {
      return 1;
    }
    mesh = meshfactory.create(meshExodus);

    /* create Darcy process kernel */
    Teuchos::ParameterList state_list = plist->sublist("state");
    S = Teuchos::rcp(new State(state_list));
    S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

    Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
    DPK = new Darcy_PK(plist, "flow", S, soln);
    DPK->Setup();
    S->Setup();
    S->set_time(0.0);

    S->InitializeFields();
    S->InitializeEvaluators();

    // create other parameters
    mu = S->Get<double>("const_fluid_viscosity");
    rho = S->Get<double>("const_fluid_density");
    gravity = S->Get<AmanziGeometry::Point>("gravity");

    return 0;
  }

  void createBClist(const char* type,
                    const char* bc_x,
                    Teuchos::Array<std::string>& regions,
                    double value)
  {
    std::string func_list_name;
    if (!strcmp(type, "pressure")) {
      func_list_name = "boundary pressure";
    } else if (!strcmp(type, "static head")) {
      func_list_name = "water table elevation";
    } else if (!strcmp(type, "mass flux")) {
      func_list_name = "outward mass flux";
    }
    Teuchos::ParameterList& dp_list = plist->sublist("PKs").get<Teuchos::ParameterList>("flow");

    Teuchos::ParameterList& bc_list = dp_list.get<Teuchos::ParameterList>("boundary conditions");
    Teuchos::ParameterList& type_list = bc_list.get<Teuchos::ParameterList>(type);

    Teuchos::ParameterList& bc_sublist = type_list.sublist(bc_x);
    bc_sublist.set("regions", regions).set("spatial distribution method", "none");

    if (!strcmp(type, "static head")) {
      Teuchos::ParameterList& bc_sublist_named = bc_sublist.sublist("static head")
                                                   .sublist("function-static-head")
                                                   .set("p0", 1.0)
                                                   .set("density", 2.0)
                                                   .set("gravity", 2.0)
                                                   .set("space dimension", 3)
                                                   .sublist(func_list_name);
      Teuchos::ParameterList& function_list = bc_sublist_named.sublist("function-constant");
      function_list.set("value", value);
    } else {
      Teuchos::ParameterList& bc_sublist_named = bc_sublist.sublist(func_list_name);
      Teuchos::ParameterList& function_list = bc_sublist_named.sublist("function-constant");
      function_list.set("value", value);
    }
  }

  double calculatePressureCellError(double p0, AmanziGeometry::Point& pressure_gradient)
  {
    const auto& pressure = *S->Get<CompositeVector>("pressure").ViewComponent("cell");

    double error_L2 = 0.0;
    int ncells =
      mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
    for (int c = 0; c < ncells; c++) {
      const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
      double pressure_exact = p0 + pressure_gradient * xc;
      // if (MyPID==0) std::cout << c << " " << pressure[0][c] << " exact=" <<  pressure_exact << std::endl;
      error_L2 += std::pow(pressure[0][c] - pressure_exact, 2.0);
    }
    return sqrt(error_L2);
  }

  double calculatePressureFaceError(double p0, AmanziGeometry::Point& pressure_gradient)
  {
    const auto& lambda = *S->Get<CompositeVector>("pressure").ViewComponent("face");

    double error_L2 = 0.0;
    int nfaces =
      mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
    for (int f = 0; f < nfaces; f++) {
      const AmanziGeometry::Point& xf = mesh->getFaceCentroid(f);
      double pressure_exact = p0 + pressure_gradient * xf;
      // std::cout << f << " " << lambda[0][f] << " exact=" << pressure_exact << std::endl;
      error_L2 += std::pow(lambda[0][f] - pressure_exact, 2.0);
    }
    return sqrt(error_L2);
  }

  double calculateDarcyFluxError(AmanziGeometry::Point& velocity_exact)
  {
    auto& flowrate = *S->Get<CompositeVector>("volumetric_flow_rate").ViewComponent("face");

    double error_L2 = 0.0;
    int nfaces =
      mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
    for (int f = 0; f < nfaces; f++) {
      const AmanziGeometry::Point& normal = mesh->getFaceNormal(f);
      error_L2 += std::pow(flowrate[0][f] - velocity_exact * normal, 2.0);
      // if (MyPID == 0) std::cout << f << " " << flowrate[0][f] << " exact=" << velocity_exact * normal << std::endl;
    }
    return sqrt(error_L2);
  }

  double calculateDarcyDivergenceError()
  {
    auto& cv = S->GetW<CompositeVector>("volumetric_flow_rate", Tags::DEFAULT, passwd);
    cv.ScatterMasterToGhosted("face");
    Epetra_MultiVector& flux = *cv.ViewComponent("face", true);

    double error_L2 = 0.0;
    int ncells_owned =
      mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);

    for (int c = 0; c < ncells_owned; c++) {
      AmanziMesh::Entity_ID_List faces;
      std::vector<int> dirs;

      mesh->getCellFacesAndDirections(c, &faces, &dirs);
      int nfaces = faces.size();

      double div = 0.0;
      for (int i = 0; i < nfaces; i++) {
        int f = faces[i];
        div += flux[0][f] * dirs[i];
      }
      error_L2 += div * div / mesh->getCellVolume(c);
    }
    return sqrt(error_L2);
  }
};


SUITE(Darcy_PK)
{
  /* ******************************************************************
* Testing the mesh of hexahedra.
****************************************************************** */
  TEST_FIXTURE(DarcyProblem, DirichletDirichlet)
  {
    for (int i = 0; i < 2; i++) {
      int ierr = Init("test/flow_darcy_misc.xml", "test/hexes.exo", frameworks[i]);
      if (MyPID == 0)
        std::cout << "\nDarcy PK on tets: Dirichlet-Dirichlet"
                  << ", mesh framework: " << framework_name[i] << std::endl;
      if (ierr == 1) continue;

      Teuchos::Array<std::string> regions(1); // modify boundary conditions
      regions[0] = std::string("Top");
      createBClist("pressure", "BC 1", regions, 0.0);

      regions[0] = std::string("Bottom");
      createBClist("pressure", "BC 2", regions, 1.0);
      // DPK->ResetParameterList(dp_list);

      DPK->Initialize();
      S->CheckAllFieldsInitialized();

      DPK->SolveFullySaturatedProblem(S->GetW<CompositeVector>("pressure", passwd), false);
      DPK->CommitStep(0.0, 1.0, Tags::DEFAULT);

      // calculate errors
      double p0 = 1.0;
      AmanziGeometry::Point pressure_gradient(0.0, 0.0, -1.0);
      AmanziGeometry::Point velocity(3);
      velocity = -(pressure_gradient - rho * gravity) / mu;

      double errorP = calculatePressureCellError(p0, pressure_gradient);
      CHECK(errorP < 1.0e-8);
      double errorL = calculatePressureFaceError(p0, pressure_gradient);
      CHECK(errorL < 1.0e-8);
      double errorU = calculateDarcyFluxError(velocity);
      CHECK(errorU < 1.0e-8);
      double errorDiv = calculateDarcyDivergenceError();
      if (MyPID == 0)
        std::cout << "Error: " << errorP << " " << errorL << " " << errorU << " " << errorDiv
                  << std::endl;
    }
  }

  TEST_FIXTURE(DarcyProblem, DirichletNeumann)
  {
    for (int i = 0; i < 2; i++) {
      int ierr = Init("test/flow_darcy_misc.xml", "test/hexes.exo", frameworks[i]);
      if (MyPID == 0)
        std::cout << "\nDarcy PK on tets: Dirichlet-Neumann"
                  << ", mesh framework: " << framework_name[i] << std::endl;
      if (ierr == 1) continue;

      Teuchos::Array<std::string> regions(1); // modify boundary conditions
      regions[0] = std::string("Top");
      createBClist("mass flux", "BC 1", regions, -20.0);

      regions[0] = std::string("Bottom");
      createBClist("pressure", "BC 2", regions, 1.0);

      DPK->Initialize();
      S->CheckAllFieldsInitialized();

      DPK->SolveFullySaturatedProblem(S->GetW<CompositeVector>("pressure", passwd), false);
      DPK->CommitStep(0.0, 1.0, Tags::DEFAULT);

      // calculate errors
      double p0 = 1.0;
      AmanziGeometry::Point pressure_gradient(0.0, 0.0, -1.0);
      AmanziGeometry::Point velocity(3);
      velocity = -(pressure_gradient - rho * gravity) / mu;

      double errorP = calculatePressureCellError(p0, pressure_gradient);
      CHECK(errorP < 1.0e-8);
      double errorL = calculatePressureFaceError(p0, pressure_gradient);
      CHECK(errorL < 1.0e-8);
      double errorU = calculateDarcyFluxError(velocity);
      CHECK(errorU < 1.0e-8);
      if (MyPID == 0)
        std::cout << "Error: " << errorP << " " << errorL << " " << errorU << std::endl;
    }
  }

  TEST_FIXTURE(DarcyProblem, StaticHeadDirichlet)
  {
    for (int i = 0; i < 2; i++) {
      int ierr = Init("test/flow_darcy_misc.xml", "test/hexes.exo", frameworks[i]);
      if (MyPID == 0)
        std::cout << "\nDarcy PK on tets: StaticHead-Dirichlet"
                  << ", mesh framework: " << framework_name[i] << std::endl;
      if (ierr == 1) continue;

      Teuchos::Array<std::string> regions(1); // modify boundary conditions
      regions[0] = std::string("Top");
      createBClist("pressure", "BC 1", regions, 1.0);

      regions[0] = std::string("Bottom");
      createBClist("static head", "BC 2", regions, 0.25);
      // DPK->ResetParameterList(dp_list);

      DPK->Initialize();
      S->CheckAllFieldsInitialized();

      DPK->SolveFullySaturatedProblem(S->GetW<CompositeVector>("pressure", passwd), false);
      DPK->CommitStep(0.0, 1.0, Tags::DEFAULT);

      // calculate errors
      double p0 = 2.0;
      AmanziGeometry::Point pressure_gradient(0.0, 0.0, -1.0);
      AmanziGeometry::Point velocity(3);
      velocity = -(pressure_gradient - rho * gravity) / mu;

      double errorP = calculatePressureCellError(p0, pressure_gradient);
      CHECK(errorP < 1.0e-8);
      double errorL = calculatePressureFaceError(p0, pressure_gradient);
      CHECK(errorL < 1.0e-8);
      double errorU = calculateDarcyFluxError(velocity);
      CHECK(errorU < 1.0e-8);
      if (MyPID == 0)
        std::cout << "Error: " << errorP << " " << errorL << " " << errorU << std::endl;
    }
  }

  /* ******************************************************************
* Testing the mesh of prisms.
****************************************************************** */
  TEST_FIXTURE(DarcyProblem, DDprisms)
  {
    for (int i = 0; i < 2; i++) {
      int ierr = Init("test/flow_darcy_misc.xml", "test/prisms.exo", frameworks[i]);
      if (MyPID == 0)
        std::cout << "\nDarcy PK on tets: Dirichlet-Dirichlet"
                  << ", mesh framework: " << framework_name[i] << std::endl;
      if (ierr == 1) continue;

      Teuchos::Array<std::string> regions(1); // modify boundary conditions
      regions[0] = std::string("Top");
      createBClist("pressure", "BC 1", regions, 0.0);

      regions[0] = std::string("Bottom");
      createBClist("pressure", "BC 2", regions, 1.0);
      // DPK->ResetParameterList(dp_list);

      DPK->Initialize();
      S->CheckAllFieldsInitialized();

      DPK->SolveFullySaturatedProblem(S->GetW<CompositeVector>("pressure", passwd), false);
      DPK->CommitStep(0.0, 1.0, Tags::DEFAULT);

      // calculate errors
      double p0 = 1.0;
      AmanziGeometry::Point pressure_gradient(0.0, 0.0, -1.0);
      AmanziGeometry::Point velocity(3);
      velocity = -(pressure_gradient - rho * gravity) / mu;

      double errorP = calculatePressureCellError(p0, pressure_gradient);
      CHECK(errorP < 1.0e-8);
      double errorL = calculatePressureFaceError(p0, pressure_gradient);
      CHECK(errorL < 1.0e-8);
      double errorU = calculateDarcyFluxError(velocity);
      CHECK(errorU < 1.0e-8);
      double errorDiv = calculateDarcyDivergenceError();
      if (MyPID == 0)
        std::cout << "Error: " << errorP << " " << errorL << " " << errorU << " " << errorDiv
                  << std::endl;
    }
  }

  /* ******************************************************************
* Testing the mesh of tetrahedra.
****************************************************************** */
  TEST_FIXTURE(DarcyProblem, DNtetrahedra)
  {
    for (int i = 0; i < 2; i++) {
      int ierr = Init("test/flow_darcy_misc.xml", "test/tetrahedra.exo", frameworks[i]);
      if (MyPID == 0)
        std::cout << "\nDarcy PK on tets: Dirichlet-Neumann"
                  << ", mesh framework: " << framework_name[i] << std::endl;
      if (ierr == 1) continue;

      Teuchos::Array<std::string> regions(1); // modify boundary conditions
      regions[0] = std::string("Top");
      createBClist("mass flux", "BC 1", regions, -20.0);

      regions[0] = std::string("Bottom");
      createBClist("pressure", "BC 2", regions, 1.0);
      // DPK->ResetParameterList(dp_list);

      DPK->Initialize();
      S->CheckAllFieldsInitialized();

      DPK->SolveFullySaturatedProblem(S->GetW<CompositeVector>("pressure", passwd), false);
      DPK->CommitStep(0.0, 1.0, Tags::DEFAULT);

      // calculate errors
      double p0 = 1.0;
      AmanziGeometry::Point pressure_gradient(0.0, 0.0, -1.0);
      AmanziGeometry::Point velocity(3);
      velocity = -(pressure_gradient - rho * gravity) / mu;

      double errorP = calculatePressureCellError(p0, pressure_gradient);
      CHECK(errorP < 1.0e-8);
      double errorL = calculatePressureFaceError(p0, pressure_gradient);
      CHECK(errorL < 1.0e-8);
      double errorU = calculateDarcyFluxError(velocity);
      CHECK(errorU < 1.0e-8);
      if (MyPID == 0)
        std::cout << "Error: " << errorP << " " << errorL << " " << errorU << std::endl;
    }
  }

  /* ******************************************************************
* Testing the mesh of prisms and hexahedra.
****************************************************************** */
  TEST_FIXTURE(DarcyProblem, DDmixed)
  {
    for (int i = 0; i < 2; i++) {
      int ierr = Init("test/flow_darcy_misc.xml", "test/mixed.exo", frameworks[i]);
      if (MyPID == 0)
        std::cout << "\nDarcy PK on mixed mesh: Dirichlet-Dirichlet"
                  << ", mesh framework: " << framework_name[i] << std::endl;
      if (ierr == 1) continue;

      Teuchos::Array<std::string> regions(1); // modify boundary conditions
      regions[0] = std::string("Top");
      createBClist("pressure", "BC 1", regions, 0.0);

      regions[0] = std::string("Bottom");
      createBClist("pressure", "BC 2", regions, 1.0);
      // DPK->ResetParameterList(dp_list);

      DPK->Initialize();
      S->CheckAllFieldsInitialized();

      DPK->SolveFullySaturatedProblem(S->GetW<CompositeVector>("pressure", passwd), false);
      DPK->CommitStep(0.0, 1.0, Tags::DEFAULT);

      // calculate errors
      double p0 = 1.0;
      AmanziGeometry::Point pressure_gradient(0.0, 0.0, -1.0);
      AmanziGeometry::Point velocity(3);
      velocity = -(pressure_gradient - rho * gravity) / mu;

      double errorP = calculatePressureCellError(p0, pressure_gradient);
      CHECK(errorP < 1.0e-8);
      double errorL = calculatePressureFaceError(p0, pressure_gradient);
      CHECK(errorL < 1.0e-8);
      double errorU = calculateDarcyFluxError(velocity);
      CHECK(errorU < 1.0e-8);
      double errorDiv = calculateDarcyDivergenceError();
      if (MyPID == 0)
        std::cout << "Error: " << errorP << " " << errorL << " " << errorU << " " << errorDiv
                  << std::endl;
    }
  }
} // SUITE

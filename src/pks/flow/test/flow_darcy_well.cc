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

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "State.hh"

// Flow
#include "Darcy_PK.hh"


/* *********************************************************************
* Two tests with different time step controllers.
********************************************************************* */
void
RunTestDarcyWell(std::string controller, bool fit)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: 2D well model: " << controller << std::endl;

  // read parameter list
  std::string xmlFileName = controller;
  Teuchos::RCP<ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create an MSTK mesh framework
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(2, regions_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);
  //pref.push_back(Framework::MOAB);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  RCP<const Mesh> mesh = meshfactory.create(-10.0, -5.0, 10.0, 0.0, 101, 50);

  /*
  Teuchos::ParameterList mesh_info_list;
  std::string filename = controller.replace(controller.size()-4, 4, "_mesh");
  mesh_info_list.set<std::string>("filename", filename);
  Teuchos::RCP<Amanzi::MeshInfo> mesh_info = Teuchos::rcp(new Amanzi::MeshInfo(mesh_info_list, &comm));
  mesh_info->WriteMeshCentroids(*mesh);
  */

  // create a simple state and populate it
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  auto soln = Teuchos::rcp(new TreeVector());
  auto DPK = Teuchos::rcp(new Darcy_PK(plist, "flow", S, soln));
  DPK->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // modify the default state for the problem at hand
  // -- permeability
  std::string passwd("");
  auto& K = *S->GetW<CompositeVector>("permeability", "permeability").ViewComponent("cell");
  double diff_in_perm = 0.0;

  if (!S->GetRecord("permeability").initialized()) {
    for (int c = 0; c < K.MyLength(); c++) {
      const AmanziGeometry::Point xc = mesh->getCellCentroid(c);
      K[0][c] = 0.1 + std::sin(xc[0]) * 0.02;
      K[1][c] = 2.0 + std::cos(xc[1]) * 0.4;
    }
    S->GetRecordW("permeability", "permeability").set_initialized();
  } else {
    for (int c = 0; c < K.MyLength(); c++) {
      const AmanziGeometry::Point xc = mesh->getCellCentroid(c);
      diff_in_perm += abs(K[0][c] - (0.1 + std::sin(xc[0]) * 0.02)) +
                      abs(K[1][c] - (2.0 + std::cos(xc[1]) * 0.4));
    }
    std::cout << "diff_in_perm " << diff_in_perm << "\n";
    CHECK(diff_in_perm < 1.0e-12);
  }


  // -- fluid density and viscosity
  S->GetW<double>("const_fluid_density", "state") = 1.0;
  S->GetRecordW("const_fluid_density", "state").set_initialized();

  S->GetW<double>("const_fluid_viscosity", "state") = 1.0;
  S->GetRecordW("const_fluid_viscosity", "state").set_initialized();

  // initialize the Darcy process kernel
  DPK->Initialize();

  std::string filename = controller.replace(controller.size() - 4, 4, "_flow2D.gmv");

  // transient solution
  double t_old(0.0), t_new, dt(0.5);
  for (int n = 0; n < 10; n++) {
    t_new = t_old + dt;

    DPK->AdvanceStep(t_old, t_new);
    DPK->CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;

    const auto& p = *S->Get<CompositeVector>("pressure").ViewComponent("cell");
    if (MyPID == 0) {
      // filename = controller.replace(controller.size()-5, 5, "_flow2D.gmv");
      // GMV::open_data_file(*mesh, filename);
      GMV::open_data_file(*mesh, (std::string) "flow2D.gmv");
      GMV::start_data();
      GMV::write_cell_data(p, 0, "pressure");
      GMV::close_data_file();
    }

    dt = DPK->get_dt();

    for (int c = 0; c < K.MyLength(); c++) {
      if (fit) {
        const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
        if (fabs(xc[0]) < 0.05 && fabs(xc[1] + 2.475) < 0.05) {
          // use quadratic approximation in time. This may capture bugs in the well model.
          double quad =
            -0.00068556 * (n - 4) * (n - 9) + 0.021082 * n * (n - 9) - 0.032341 * n * (n - 4);
          CHECK_CLOSE(p[0][c], quad, fabs(quad) * 0.1);
        }
      }
    }
  }
}


TEST(FLOW_2D_DARCY_WELL_STANDARD)
{
  RunTestDarcyWell("test/flow_darcy_well.xml", true);
}

TEST(FLOW_2D_DARCY_WELL_FROMFILE)
{
  RunTestDarcyWell("test/flow_darcy_well_fromfile.xml", false);
}

TEST(FLOW_2D_DARCY_WELL_ADAPRIVE)
{
  RunTestDarcyWell("test/flow_darcy_well_adaptive.xml", false);
}


/* **************************************************************** */
void
Run_3D_DarcyWell(std::string controller)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: 3D Darcy flow, two wells" << std::endl;

  // read parameter list
  std::string xmlFileName = controller;
  Teuchos::RCP<ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create an MSTK mesh framework
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, regions_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);
  //pref.push_back(Framework::MOAB);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  RCP<const Mesh> mesh = meshfactory.create(-10.0, -1.0, -5.0, 10.0, 1.0, 0.0, 101, 10, 25);

  // create a simple state and populate it
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<Darcy_PK> DPK = Teuchos::rcp(new Darcy_PK(plist, "flow", S, soln));
  DPK->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // modify the default state for the problem at hand
  // -- permeability
  std::string passwd("");
  auto& K = *S->GetW<CompositeVector>("permeability", "permeability").ViewComponent("cell");

  for (int c = 0; c < K.MyLength(); c++) {
    K[0][c] = 0.1;
    K[1][c] = 2.0;
    K[2][c] = 2.0;
  }
  S->GetRecordW("permeability", "permeability").set_initialized();

  // -- fluid density and viscosity
  S->GetW<double>("const_fluid_density", "state") = 1.0;
  S->GetRecordW("const_fluid_density", "state").set_initialized();

  S->GetW<double>("const_fluid_viscosity", "state") = 1.0;
  S->GetRecordW("const_fluid_viscosity", "state").set_initialized();

  // initialize the Darcy process kernel
  DPK->Initialize();

  std::string filename = controller.replace(controller.size() - 4, 4, "_flow3D.gmv");

  // transient solution
  double t_old(0.0), t_new, dt(0.5);
  for (int n = 0; n < 2; n++) {
    t_new = t_old + dt;

    DPK->AdvanceStep(t_old, t_new);
    DPK->CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;

    if (MyPID == 0) {
      const auto& p = *S->Get<CompositeVector>("pressure").ViewComponent("cell");
      GMV::open_data_file(*mesh, filename);
      GMV::start_data();
      GMV::write_cell_data(p, 0, "pressure");
      GMV::close_data_file();
    }
  }
}


//TEST(FLOW_3D_DARCY_WELL)
//{
//  Run_3D_DarcyWell("test/flow_darcy_well_3D.xml");
//}


/* **************************************************************** */
TEST(FLOW_3D_DARCY_PEACEMAN_WELL)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: 3D Darcy flow, one well" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/flow_darcy_1well_peaceman_3D.xml";
  Teuchos::RCP<ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create an MSTK mesh framework
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, regions_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);
  //pref.push_back(Framework::MOAB);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  RCP<const Mesh> mesh = meshfactory.create(-55.0, -55.0, -2., 55.0, 55.0, 0.0, 23, 23, 4);

  // create a simple state and populate it
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  std::string passwd("");
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<Darcy_PK> DPK = Teuchos::rcp(new Darcy_PK(plist, "flow", S, soln));
  DPK->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // modify the default state for the problem at hand
  // -- permeability
  auto& K = *S->GetW<CompositeVector>("permeability", "permeability").ViewComponent("cell");

  for (int c = 0; c < K.MyLength(); c++) {
    K[0][c] = 10.;
    K[1][c] = 10.;
    K[2][c] = 1.;
  }
  S->GetRecordW("permeability", "permeability").set_initialized();

  // -- fluid density and viscosity
  S->GetW<double>("const_fluid_density", "state") = 1.0;
  S->GetRecordW("const_fluid_density", "state").set_initialized();

  S->GetW<double>("const_fluid_viscosity", "state") = 1.0;
  S->GetRecordW("const_fluid_viscosity", "state").set_initialized();

  // -- porosity
  S->GetW<CompositeVector>("porosity", "porosity").PutScalar(1.0);
  S->GetRecordW("porosity", "porosity").set_initialized();

  // initialize the Darcy process kernel
  DPK->Initialize();

  //std::string filename = controller.replace(controller.size()-4, 4, "_flow3D.gmv");
  std::string filename = "flow_darcy_well_peaceman_3D.gmv";

  // steady_state solution
  DPK->SolveFullySaturatedProblem(S->GetW<CompositeVector>("pressure", passwd), true);

  const auto& p = *S->Get<CompositeVector>("pressure").ViewComponent("cell");
  Epetra_MultiVector err_p(p), p_exact(p);

  auto gravity = S->Get<AmanziGeometry::Point>("gravity");

  double Q = -100.0;
  double rw = 0.1;
  double k = 10.0;
  double h = 5.0;
  double pw = 10.0;
  double depth = 2.5;

  int ncells =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  double err = 0.;
  double sol = 0.;
  for (int c = 0; c < ncells; c++) {
    const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
    double r = sqrt(xc[0] * xc[0] + xc[1] * xc[1]);
    double p_ex;
    p_ex = pw + gravity[2] * (xc[2] + depth);
    if (r > 1e-3) {
      p_ex = p_ex + Q / (2 * M_PI * k * h) * (log(r) - log(rw));
    } else {
      p_ex = p[0][c];
    }

    p_exact[0][c] = p_ex;

    double vol = mesh->getCellVolume(c);
    err += (p_ex - p[0][c]) * (p_ex - p[0][c]) * vol;

    err_p[0][c] = abs(p_ex - p[0][c]);

    sol += p[0][c] * p[0][c] * vol;
  }

  // for (double r=55; r < 75; r+=0.25) {
  //   double p_ex = pw + Q/(2*M_PI*k*h)*(log(r) - log(rw));
  //   std::cout<<r<<" "<<p_ex<<" "<<Q/(2*M_PI*k*h)<<" "<<log(rw)<<"\n";
  // }

  err = sqrt(err);
  sol = sqrt(sol);
  err /= sol;
  std::cout << "Error: " << err << "\n";

  CHECK(err < 0.02);

  if (MyPID == 0) {
    GMV::open_data_file(*mesh, filename);
    GMV::start_data();
    GMV::write_cell_data(p, 0, "pressure");
    GMV::write_cell_data(p_exact, 0, "exact");
    GMV::write_cell_data(err_p, 0, "error");
    if (S->HasRecord("well_index")) {
      const auto& wi = *S->Get<CompositeVector>("well_index").ViewComponent("cell");
      GMV::write_cell_data(wi, 0, "well_index");
    }
    GMV::close_data_file();
  }
}

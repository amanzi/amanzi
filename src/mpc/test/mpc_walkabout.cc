#include <iostream>
#include "stdlib.h"
#include "math.h"


#include <Epetra_MpiComm.h>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

#include "CycleDriver.hh"
#include "eos_registration.hh"
#include "MeshFactory.hh"
#include "Mesh.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_flow_registration.hh"
#include "State.hh"
#include "wrm_flow_registration.hh"


TEST(MPC_WALKABOUT_2D) {
using namespace Amanzi;

  auto comm = Amanzi::getDefaultComm();
  
  // read the main parameter list
  std::string xmlFileName = "test/mpc_walkabout_2D.xml";
  Teuchos::RCP<Teuchos::ParameterList> glist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = glist->sublist("regions");
  Teuchos::RCP<AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new AmanziGeometry::GeometricModel(2, region_list, *comm));

  // create mesh
  auto mesh_list = Teuchos::sublist(glist, "mesh");
  AmanziMesh::MeshFactory meshfactory(comm, gm, mesh_list);

  meshfactory.set_preference(AmanziMesh::Preference({AmanziMesh::Framework::MSTK}));
  auto mesh = meshfactory.create("test/mpc_walkabout_2D.exo");

  Teuchos::ParameterList state_plist = glist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);
  
  // use cycle driver to create and initialize state
  ObservationData obs_data;    
  CycleDriver cycle_driver(glist, S, comm, obs_data);
  S = cycle_driver.Go();

  // verify no-flow at selected points using existing S
  std::cout << "Start test of 2D Walkabout\n";
  AmanziGeometry::Point xv(2);
  std::vector<AmanziGeometry::Point> xyz, velocity;
  cycle_driver.walkabout()->CalculateDarcyVelocity(S, xyz, velocity);

  if (comm->NumProc() == 1) {
    std::vector<int> list = {1, 2, 3, 16, 17, 18};
    for (int v : list) { 
      mesh->node_get_coordinates(v, &xv);
      CHECK_CLOSE(0.0, xv * velocity[v], 1e-14);
    }
  }

  // verify velocity at all points
  int nnodes = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  auto& flow = *S->GetFieldData("darcy_flux", "flow")->ViewComponent("face");
  auto& pres = *S->GetFieldData("pressure", "flow")->ViewComponent("cell");
  
  // -- overwite with constant velocity
  AmanziGeometry::Point vel(1.0, 2.0);
  for (int f = 0; f < nfaces; ++f) {
    flow[0][f] = vel * mesh->face_normal(f);
  }

  // -- overwite with linear pressure
  for (int c = 0; c < ncells; ++c) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    pres[0][c] = 1.0 + xc[0] + 2 * xc[1];
  }

  // -- check recovered velocity
  Teuchos::ParameterList& wlist = glist->sublist("walkabout data");
  auto walkabout = Teuchos::rcp(new Amanzi::WalkaboutCheckpoint(wlist, *S));

  walkabout->CalculateDarcyVelocity(S, xyz, velocity);

  for (int v = 0; v < nnodes; ++v) {
    CHECK(norm(vel - velocity[v]) < 1e-10);
  }

  // -- check interpolated pressure and saturation
  std::vector<int> material_ids;
  std::vector<double> porosity, saturation, pressure, isotherm_kd;

  walkabout->CalculateData(S, xyz, velocity,
                           porosity, saturation, pressure, isotherm_kd, material_ids);

  for (int v = 0; v < nnodes; ++v) {
    mesh->node_get_coordinates(v, &xv);
    CHECK_CLOSE(1.0 + xv[0] + 2 * xv[1], pressure[v], 0.15);  // some tets have 1 neighboor
    CHECK_CLOSE(1.0, saturation[v], 1e-10);
  }

  // create walkabout file without pk
  walkabout->disable(false);
  walkabout->WriteDataFile(S, Teuchos::null);
}


TEST(MPC_WALKABOUT_3D) {
using namespace Amanzi;

  auto comm = Amanzi::getDefaultComm();
  
  // read the main parameter list
  std::string xmlFileName = "test/mpc_walkabout_3D.xml";
  Teuchos::RCP<Teuchos::ParameterList> glist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = glist->sublist("regions");
  Teuchos::RCP<AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new AmanziGeometry::GeometricModel(3, region_list, *comm));

  // create mesh
  auto mesh_list = Teuchos::sublist(glist, "mesh");
  AmanziMesh::MeshFactory meshfactory(comm, gm, mesh_list);

  meshfactory.set_preference(AmanziMesh::Preference({AmanziMesh::Framework::MSTK}));
  auto mesh = meshfactory.create("test/mpc_walkabout_tet5.exo");
  
  Teuchos::ParameterList state_plist = glist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);
  
  // use cycle driver to create and initialize state
  ObservationData obs_data;    
  CycleDriver cycle_driver(glist, S, comm, obs_data);
  cycle_driver.Go();

  // verify velocity at all points
  // -- overwrite flow & pressure
  std::cout << "Start test of 3D Walkabout\n";
  AmanziGeometry::Point vel(1.0, 2.0, 3.0);
  int nnodes = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  auto& flow = *S->GetFieldData("darcy_flux", "flow")->ViewComponent("face");
  auto& pres = *S->GetFieldData("pressure", "flow")->ViewComponent("cell");

  for (int f = 0; f < nfaces; ++f) {
    flow[0][f] = vel * mesh->face_normal(f);
  }

  pres.PutScalar(1.0);

  // -- check recovered velocity
  Teuchos::ParameterList& wlist = glist->sublist("walkabout data");
  auto walkabout = Teuchos::rcp(new Amanzi::WalkaboutCheckpoint(wlist, *S));

  std::vector<AmanziGeometry::Point> xyz, velocity;
  walkabout->CalculateDarcyVelocity(S, xyz, velocity);

  for (int v = 0; v < nnodes; ++v) {
    CHECK(norm(vel - velocity[v]) < 1e-10);
  }

  // -- check interpolated pressure
  std::vector<int> material_ids;
  std::vector<double> porosity, saturation, pressure, isotherm_kd;

  walkabout->CalculateData(S, xyz, velocity,
                           porosity, saturation, pressure, isotherm_kd, material_ids);

  for (int v = 0; v < nnodes; ++v) {
    CHECK_CLOSE(1.0, pressure[v], 1e-10);
    CHECK_CLOSE(1.0, saturation[v], 1e-10);
  }

  // verify other quantities at selected point on main diagonal
  AmanziGeometry::Point x0(0.0, 0.0, 0.0), xv(3);
  AmanziGeometry::Point x1(1.0, 1.0, 1.0);
  AmanziGeometry::Point x2(2.0, 2.0, 2.0);
  AmanziGeometry::Point x3(3.0, 3.0, 3.0);

  for (int v = 0; v < nnodes; ++v) {
    mesh->node_get_coordinates(v, &xv);
    if (norm(xv - x0) < 1e-10) {
      CHECK_CLOSE(0.2, porosity[v], 1e-10);
      CHECK_EQUAL(1000, material_ids[v]);
    } else if (norm(xv - x1) < 1e-10) {
      CHECK_CLOSE(0.3, porosity[v], 1e-10);
      CHECK_EQUAL(2000, material_ids[v]);
    } else if (norm(xv - x2) < 1e-10) {
      CHECK_CLOSE(0.5, porosity[v], 1e-10);
      CHECK_EQUAL(3000, material_ids[v]);
    } else if (norm(xv - x3) < 1e-10) {
      CHECK_CLOSE(0.6, porosity[v], 1e-10);
      CHECK_EQUAL(3000, material_ids[v]);
    }
  }

  // create walkabout file without flow pk
  walkabout->disable(false);
  walkabout->WriteDataFile(S, Teuchos::null);
}



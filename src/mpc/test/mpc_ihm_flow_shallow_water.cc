/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>
#include "stdlib.h"
#include "math.h"

// TPLs
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "CycleDriver.hh"
#include "eos_reg.hh"
#include "Mesh.hh"
#include "MeshExtractedManifold.hh"
#include "MeshFactory.hh"
#include "models_flow_reg.hh"
#include "models_shallow_water_reg.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_flow_reg.hh"
#include "pks_mpc_reg.hh"
#include "pks_shallow_water_reg.hh"
#include "State.hh"

// General
#define _USE_MATH_DEFINES
#include "math.h"


TEST(MPC_DRIVER_IHM_FLOW_SHALLOW_WATER_DAM_BREAK)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  std::string xmlInFileName = "test/mpc_ihm_flow_shallow_water.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  // create mesh
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  mesh_list->set<bool>("request edges", true);
  mesh_list->set<bool>("request faces", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({ Framework::MSTK }));
  auto mesh = factory.create(0.0, 0.0, 0.0, 10.0, 1.0, 1.0, 80, 1, 40);

  // deform mesh (if needed)
  int nnodes = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::NODE,
                                    Amanzi::AmanziMesh::Parallel_kind::OWNED);
  AmanziMesh::Entity_ID_View nodeids("nodeids", nnodes);
  AmanziMesh::Point_View new_positions("new_positions", nnodes);
  for (int n = 0; n < nnodes; ++n) {
    nodeids[n] = n;

    AmanziGeometry::Point node_crd;
    node_crd = mesh->getNodeCoordinate(n);
    new_positions[n] = node_crd;
  }
  AmanziMesh::MeshAlgorithms::deform(*mesh, nodeids, new_positions);

  // create dummy observation data object
  Amanzi::ObservationData obs_data;

  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  // create additional mesh for SW
  std::vector<std::string> names;
  names.push_back("surface");

  auto mesh_surface = factory.create(mesh, { "TopSurface" }, AmanziMesh::Entity_kind::FACE, true);

  S->RegisterMesh("surface", mesh_surface);

  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  cycle_driver.Go();

  // calculate the fluid pressure at the top of the subsurface at final time
  const auto& p = *S->Get<CompositeVector>("pressure").ViewComponent("face");
  int nfaces =
    mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  double p_top_avg = 0.0, top_surface_area = 0.0;

  for (int f = 0; f < nfaces; ++f) {
    const AmanziGeometry::Point& xf = mesh->getFaceCentroid(f);
    if (std::abs(xf[2] - 1.0) < 1.e-12) {
      p_top_avg += p[0][f] * mesh->getFaceArea(f);
      top_surface_area += mesh->getFaceArea(f);
    }
  }

  p_top_avg = p_top_avg / top_surface_area;

  std::cout << "average subsurface pressure at surface: " << p_top_avg << std::endl;

  // calculate the average ponded depth at final time
  double h_avg;
  S->Get<CompositeVector>("surface-ponded_depth").MeanValue(&h_avg);
  const double rho = S->Get<double>("const_fluid_density");
  const double patm = S->Get<double>("atmospheric_pressure");
  double g = norm(S->Get<AmanziGeometry::Point>("gravity"));

  std::cout << "average surface ponded height: " << h_avg << std::endl;
  std::cout << "average hydrostatic pressure at surface: " << patm + rho * g * h_avg << std::endl;

  // compare the values
  CHECK(std::abs(p_top_avg - (patm + rho * g * h_avg)) / (patm + rho * g * h_avg) < 1.0e-2);
}

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
#include "MeshAudit.hh"
#include "energy_tcm_registration.hh"
#include "energy_iem_registration.hh"
#include "eos_registration.hh"
#include "Mesh.hh"
#include "MeshExtractedManifold.hh"
#include "MeshFactory.hh"
#include "mpc_pks_registration.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_energy_registration.hh"
#include "pks_flow_registration.hh"
#include "pks_transport_registration.hh"
#include "State.hh"
#include "wrm_flow_registration.hh"


TEST(MPC_DRIVER_FLOW_MATRIX_FRACTURE)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  // setup a piecewice linear solution with a jump
  std::string xmlInFileName = "test/mpc_coupled_flow_source.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  // create mesh
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  mesh_list->set<bool>("request edges",true); 
  mesh_list->set<bool>("request faces",true); 
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({ Framework::MSTK }));
  auto mesh = factory.create(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 11, 11, 10);

  // create dummy observation data object
  Amanzi::ObservationData obs_data;

  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  // create mesh for fracture
  std::vector<std::string> names;
  names.push_back("fracture");
  auto mesh_fracture_framework = Teuchos::rcp(new MeshExtractedManifold(
    mesh, "fracture", AmanziMesh::Entity_kind::FACE, comm, gm, mesh_list));
  auto mesh_fracture = Teuchos::rcp(new Mesh(mesh_fracture_framework, Teuchos::rcp(new Amanzi::AmanziMesh::MeshFrameworkAlgorithms()), mesh_list)); 

  S->RegisterMesh("fracture", mesh_fracture);

  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  cycle_driver.Go();

  // verify solution symmetry
  int ncells =
    mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL, Amanzi::AmanziMesh::Parallel_kind::OWNED);
  const auto& pc = *S->Get<CompositeVector>("pressure").ViewComponent("cell");

  Amanzi::AmanziGeometry::Point x0(5.0, 5.0, 5.0);
  std::vector<double> dist({ 1.03752, 1.88568, 2.77273, 3.67058 });
  for (double d : dist) {
    double pmin(1e+10), pmax(-1e+10);
    for (int c = 0; c < ncells; ++c) {
      const auto& xc = mesh->getCellCentroid(c);
      if (std::fabs(norm(xc - x0) - d) < 1e-4 && std::fabs(xc[2] - 4.5) < 1e-6) {
        pmin = std::min(pmin, pc[0][c]);
        pmax = std::max(pmax, pc[0][c]);
      }
    }
    CHECK_CLOSE(pmin, pmax, 1e-5);
  }
}

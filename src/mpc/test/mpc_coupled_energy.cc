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
#include "eos_reg.hh"
#include "evaluators_flow_reg.hh"
#include "evaluators_mpc_reg.hh"
#include "Mesh.hh"
#include "MeshExtractedManifold.hh"
#include "MeshFactory.hh"
#include "models_energy_reg.hh"
#include "models_flow_reg.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_energy_reg.hh"
#include "pks_mpc_reg.hh"
#include "pks_flow_reg.hh"
#include "pks_transport_reg.hh"
#include "State.hh"


TEST(MPC_DRIVER_ENERGY_MATRIX_FRACTURE)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  // setup a piecewice linear solution with a jump
  std::string xmlInFileName = "test/mpc_coupled_energy.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  // create mesh
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  ;
  mesh_list->set<bool>("request edges", true);
  mesh_list->set<bool>("request faces", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({ Framework::MSTK }));
  auto mesh = factory.create(400.0, 0.0, 0.0, 600.0, 100.0, 1000.0, 20, 5, 100);

  // create dummy observation data object
  Amanzi::ObservationData obs_data;

  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  //create additional mesh for fracture
  std::vector<std::string> names;
  names.push_back("fracture");
  // auto mesh_fracture = factory.create(mesh, names, AmanziMesh::FACE);
  auto mesh_fracture_mf = Teuchos::rcp(
    new MeshExtractedManifold(mesh, "fracture", AmanziMesh::FACE, comm, gm, mesh_list));
  auto mesh_fracture = Teuchos::rcp(
    new Mesh(mesh_fracture_mf, Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms()), mesh_list));

  S->RegisterMesh("fracture", mesh_fracture);

  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  S = cycle_driver.Go();

  // verification
  // -- baounds preservation
  int ncells =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  const auto& T_m = *S->Get<CompositeVector>("temperature").ViewComponent("cell");
  for (int c = 0; c < ncells; ++c) {
    double T = T_m[0][c];
    CHECK(T < 573.15 && T > 473.15);
  }

  ncells =
    mesh_fracture->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  const auto& T_f = *S->Get<CompositeVector>("fracture-temperature").ViewComponent("cell");
  for (int c = 0; c < ncells; ++c) {
    double T = T_f[0][c];
    CHECK(T < 573.15 && T > 473.15);
  }

  // T increase in fracture
  std::cout << "T integral increase in fracture: " << std::endl;
  std::string label = obs_data.observationLabels()[0];
  double tmp(0.0);
  for (auto& quad : obs_data[label]) {
    CHECK(quad.value > tmp);
    tmp = quad.value;
    quad.print(std::cout);
  }

  // T increase in matrix
  std::cout << "\nT integral increase in matrix: " << std::endl;
  label = obs_data.observationLabels()[1];
  tmp = 0.0;
  for (auto& quad : obs_data[label]) {
    CHECK(quad.value > tmp);
    tmp = quad.value;
    quad.print(std::cout);
  }
}

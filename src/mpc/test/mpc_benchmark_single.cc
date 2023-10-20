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
#include "Mesh.hh"
#include "MeshExtractedManifold.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "models_flow_reg.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_flow_reg.hh"
#include "pks_mpc_reg.hh"
#include "pks_transport_reg.hh"
#include "State.hh"


void
MPC_CoupledFlowTransport(const std::string& xmlfile, const std::string& exofile)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  // setup a piecewice linear solution with a jump
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlfile);

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  // create mesh
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  MeshFactory factory(comm, gm, mesh_list);
  mesh_list->set<bool>("request edges", true);
  mesh_list->set<bool>("request faces", true);
  factory.set_preference(Preference({ Framework::MSTK }));
  auto mesh = factory.create(exofile);

  // create dummy observation data object
  Amanzi::ObservationData obs_data;

  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));

  S->RegisterMesh("domain", mesh);

  /*
  Amanzi::MeshAudit mesh_auditor(mesh);
  int status = mesh_auditor.Verify();
  if (status != 0) {
    Errors::Message msg("Mesh Audit could not verify correctness of mesh.");
    Exceptions::amanzi_throw(msg);
  }
  */

  //create additional mesh for fracture
  std::vector<std::string> names;
  names.push_back("fracture");

  // auto mesh_fracture = factory.create(mesh, names, AmanziMesh::Entity_kind::FACE);
  auto mesh_fracture_framework = Teuchos::rcp(new MeshExtractedManifold(
    mesh, "fracture", AmanziMesh::Entity_kind::FACE, comm, gm, mesh_list));
  auto mesh_fracture =
    Teuchos::rcp(new Mesh(mesh_fracture_framework,
                          Teuchos::rcp(new Amanzi::AmanziMesh::MeshFrameworkAlgorithms()),
                          mesh_list));

  S->RegisterMesh("fracture", mesh_fracture);

  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  cycle_driver.Go();
}


TEST(MPC_COUPLED_FLOW_TRANSPORT)
{
  MPC_CoupledFlowTransport("test/mpc_benchmark_single.xml", "test/single_fracture_tet.exo");
  // MPC_CoupledFlowTransport("test/mpc_benchmark_regular_1.xml", "test/regular_fracture_ref0.exo");
}

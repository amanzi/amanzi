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
#include "bilinear_form_reg.hh"
#include "CycleDriver.hh"
#include "MeshAudit.hh"
#include "eos_reg.hh"
#include "evaluators_reg.hh"
#include "evaluators_flow_reg.hh"
#include "Mesh.hh"
#include "MeshExtractedManifold.hh"
#include "MeshFactory.hh"
#include "models_flow_reg.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_flow_reg.hh"
#include "pks_mechanics_reg.hh"
#include "pks_mpc_reg.hh"
#include "State.hh"

// MPC
#include "EvaluatorAperture.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

void
RunTest(const std::string xmlInFileName)
{
  Comm_ptr_type comm = Amanzi::getDefaultComm();

  // setup a piecewice linear solution with a jump
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, region_list, *comm));

  // create mesh
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  mesh_list->set<bool>("request edges", true);
  mesh_list->set<bool>("request faces", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({ Framework::MSTK }));
  auto mesh = factory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 6, 6, 6);

  // create dummy observation data object
  Amanzi::ObservationData obs_data;

  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_plist));
  S->RegisterMesh("domain", mesh);

  // create additional mesh for fracture
  std::vector<std::string> names;
  names.push_back("fracture");

  auto mesh_fracture_fw = Teuchos::rcp(new MeshExtractedManifold(
    mesh, "fracture", AmanziMesh::Entity_kind::FACE, comm, gm, mesh_list));
  auto mesh_fracture = Teuchos::rcp(
    new Mesh(mesh_fracture_fw, Teuchos::rcp(new AmanziMesh::MeshAlgorithms()), mesh_list));

  S->RegisterMesh("fracture", mesh_fracture);

  CycleDriver cycle_driver(plist, S, comm, obs_data);
  S = cycle_driver.Go();
}


TEST(MPC_DRIVER_COUPLED_MECHANICS_FLOW)
{
  RunTest("test/mpc_coupled_mechanics_flow.xml");
}

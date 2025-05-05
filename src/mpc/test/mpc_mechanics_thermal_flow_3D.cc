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
#include "models_energy_reg.hh"
#include "models_flow_reg.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_energy_reg.hh"
#include "pks_flow_reg.hh"
#include "pks_mechanics_reg.hh"
#include "pks_mpc_reg.hh"
#include "State.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

void
RunTest(const std::string xmlInFileName)
{
  Comm_ptr_type comm = Amanzi::getDefaultComm();

  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  // create mesh
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({ Framework::MSTK }));
  auto mesh_parent = factory.create(-4.0, -4.0, -4.0, 4.0, 4.0, 4.0, 8, 8, 8);

  std::vector<std::string> names;
  names.push_back("DoughnutLogical");
  auto mesh = factory.create(mesh_parent, names, AmanziMesh::Entity_kind::CELL);

  // create state
  Teuchos::ParameterList state_plist = plist->sublist("state");
  auto S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  Amanzi::ObservationData obs_data;
  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  S = cycle_driver.Go();

  // verify pressure dynamics
  std::string label = obs_data.observationLabels()[0];
  int nobs = obs_data[label].size();

  for (int i = 0; i < nobs; ++i) CHECK(obs_data[label][i].value > 9.0e+5);
  CHECK(obs_data[label][0].value < obs_data[label][nobs / 2].value);
  CHECK(obs_data[label][nobs - 1].value < obs_data[label][nobs / 2].value);
}


TEST(MPC_DRIVER_THERMAL_ELASTICITY_3D)
{
  RunTest("test/mpc_mechanics_thermal_flow_3D.xml");
}

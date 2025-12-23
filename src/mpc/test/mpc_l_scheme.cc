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

#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "CycleDriver.hh"
#include "eos_reg.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "models_energy_reg.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_energy_reg.hh"
#include "pks_flow_reg.hh"
#include "pks_mpc_reg.hh"
#include "pks_transport_reg.hh"
#include "State.hh"
#include "models_flow_reg.hh"


TEST(MPC_L_SCHEME)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;

  auto comm = Amanzi::getDefaultComm();

  // read the main parameter list
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile("test/mpc_l_scheme.xml");

  std::string model("liquid water FEHM"), formula("600.0");
  auto& ev = plist->sublist("state").sublist("evaluators");
  ev.sublist("molar_density_liquid").sublist("EOS parameters").set<std::string>("eos type", model);
  ev.sublist("viscosity_liquid").sublist("EOS parameters").set<std::string>("eos type", model);

  auto& ic = plist->sublist("state").sublist("initial conditions");
  ic.sublist("temperature").sublist("function").sublist("EntireDomain").sublist("function")
    .sublist("function-exprtk").set<std::string>("formula", formula);

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, region_list, *comm));

  // create mesh
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 100.0, 20.0, 20.0, 50, 4, 4);
  AMANZI_ASSERT(!mesh.is_null());

  // create dummy observation data object
  Amanzi::ObservationData obs_data;

  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  cycle_driver.Go();

  const auto& T_c = *S->Get<CompositeVector>("temperature").ViewComponent("cell");
  for (int c = 0; c < T_c.MyLength(); ++c) {
    CHECK(300.0 < T_c[0][c] && T_c[0][c] < 600.0);
  }
}


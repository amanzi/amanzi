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


#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

#include "IO.hh"
#include "CycleDriver.hh"
#include "eos_reg.hh"
#include "GeometricModel.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_transport_reg.hh"
#include "State.hh"


TEST(MPC_DRIVER_TRANSPORT)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = Amanzi::getDefaultComm();

  // For now create one geometric model from all the regions in the spec
  auto glist1 = Teuchos::getParametersFromXmlFile("test/mpc_transport.xml");
  Teuchos::ParameterList region_list = glist1->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(2, region_list, *comm));

  // create mesh
  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);
  //pref.push_back(Framework::MOAB);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory.create("test/mpc_transport_mesh_10x10.exo");
  AMANZI_ASSERT(!mesh.is_null());

  // create dummy observation data object
  Amanzi::ObservationData obs_data;

  Teuchos::ParameterList state_plist1 = glist1->sublist("state");
  auto S1 = Teuchos::rcp(new Amanzi::State(state_plist1));
  S1->RegisterMesh("domain", mesh);

  Amanzi::CycleDriver cd1(glist1, S1, comm, obs_data);
  cd1.Go();

  int ncells = S1->GetMesh()->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,
                                             Amanzi::AmanziMesh::Parallel_kind::OWNED);
  auto tcc1_c = *S1->Get<CompositeVector>("total_component_concentration").ViewComponent("cell");
  for (int c = 0; c < ncells; ++c) {
    CHECK(tcc1_c[0][c] >= 0.0 && tcc1_c[0][c] <= 1.0);
    CHECK(tcc1_c[1][c] >= 0.0 && tcc1_c[1][c] <= 1.0);
  }

  // WriteStateStatistics(*S1);

  // -----------------
  // re-initialize tcc
  // -----------------
  auto glist2 = Teuchos::getParametersFromXmlFile("test/mpc_transport_initialize.xml");
  Teuchos::ParameterList state_plist2 = glist2->sublist("state");
  auto S2 = Teuchos::rcp(new Amanzi::State(state_plist2));
  S2->RegisterMesh("domain", mesh);

  Amanzi::CycleDriver cd2(glist2, S2, comm, obs_data);
  cd2.Go();

  auto tcc2_c = *S2->Get<CompositeVector>("total_component_concentration").ViewComponent("cell");
  for (int c = 0; c < ncells; ++c) {
    CHECK(tcc2_c[0][c] >= 0.0 && tcc2_c[0][c] <= 1.0);
    CHECK(tcc2_c[1][c] >= 0.0 && tcc2_c[1][c] <= 1.0);
  }
}

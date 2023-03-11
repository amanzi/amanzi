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
#include "eos_registration.hh"
#include "Mesh.hh"
#include "MeshExtractedManifold.hh"
#include "MeshFactory.hh"
#include "models_shallow_water_reg.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_mpc_reg.hh"
#include "pks_shallow_water_reg.hh"
#include "State.hh"


TEST(MPC_DRIVER_SHALLOW_WATER)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  // read the main parameter list
  std::string xmlInFileName = "test/mpc_shallow_water.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  // create mesh
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({ Framework::MSTK }));
  auto mesh = factory.create(0.0, 0.0, 10.0, 10.0, 20, 20, true, true);

  // create dummy observation data object
  double vmin;
  Amanzi::ObservationData obs_data;

  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("surface", mesh);

  {
    Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
    try {
      cycle_driver.Go();
      S->Get<CompositeVector>("surface-ponded_depth").MinValue(&vmin);
    } catch (...) {
      CHECK(false);
    }
  }
  S = Teuchos::null;

  CHECK(vmin > 0.0);

  // checking that we created only one pk
  CHECK(PKFactory::num_pks == 1);
  std::cout << PKFactory::list_pks << std::endl;
}

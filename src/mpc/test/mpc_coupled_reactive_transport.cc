/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>
#include <cstdlib>
#include <cmath>

// TPLs
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "CycleDriver.hh"
#include "eos_reg.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_chemistry_reg.hh"
#include "pks_mpc_reg.hh"
#include "pks_transport_reg.hh"
#include "State.hh"


TEST(MPC_DRIVER_COUPLED_REACTIVE_TRANSPORT)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  // setup a piecewice linear solution with a jump
  std::string xmlInFileName = "test/mpc_coupled_reactive_transport.xml";
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
  auto mesh = factory.create("test/single_fracture_tet.exo");

  // create dummy observation data object
  Amanzi::ObservationData obs_data;

  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  //create additional mesh for fracture
  std::vector<std::string> names;
  names.push_back("fracture");
  auto mesh_list2 = Teuchos::sublist(plist, "mesh", true);
  mesh_list2->set<bool>("request edges", false);
  mesh_list2->set<bool>("request faces", true);
  MeshFactory factory2(comm, gm, mesh_list2);
  auto mesh_fracture = factory2.create(mesh, names, AmanziMesh::Entity_kind::FACE);

  S->RegisterMesh("fracture", mesh_fracture);

  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  cycle_driver.Go();
}

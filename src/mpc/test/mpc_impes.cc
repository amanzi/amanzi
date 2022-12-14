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
#include "UnitTest++.h"


#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "State.hh"
#include "CycleDriver.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_pressuresaturation_registration.hh"
#include "pks_pressure_registration.hh"
#include "pks_saturation_registration.hh"

#include "MeshFactory.hh"
#include "Mesh.hh"
#include "Domain.hh"
#include "GeometricModel.hh"


TEST(NEW_DRIVER_COUPLED_MULTIPHASE)
{
  using namespace std;

#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  std::string xmlInFileName = "test/mpc_impes.xml";

  // read the main parameter list
  Teuchos::ParameterList driver_parameter_list;
  Teuchos::ParameterXMLFileReader xmlreader(xmlInFileName);
  driver_parameter_list = xmlreader.getParameters();

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  // create mesh
  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 8, 8);
  AMANZI_ASSERT(!mesh.is_null());

  // create dummy observation data object
  Amanzi::ObservationData obs_data;
  Teuchos::RCP<Teuchos::ParameterList> glist_rcp = Teuchos::rcp(new Teuchos::ParameterList(plist));

  Teuchos::ParameterList state_plist = glist_rcp->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  Amanzi::CycleDriver cycle_driver(glist_rcp, S, comm, obs_data);
  cycle_driver.Go();
}

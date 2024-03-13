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

#include "CycleDriver.hh"
#include "eos_reg.hh"
#include "MeshFactory.hh"
#include "Mesh.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_dummy_reg.hh"
#include "State.hh"


TEST(NEW_DRIVER_DUMMY_PK)
{
  using namespace std;

  auto comm = Amanzi::getDefaultComm();

  std::string xmlInFileName = "test/mpc_dummy.xml";

  // read the main parameter list
  Teuchos::ParameterList plist;
  Teuchos::ParameterXMLFileReader xmlreader(xmlInFileName);
  plist = xmlreader.getParameters();

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList reg_params = plist.sublist("regions");

  int spdim = 2;
  auto geom_model_ptr =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(spdim, reg_params, *comm));

  // Amanzi::AmanziGeometry::Domain *simdomain_ptr = new Amanzi::AmanziGeometry::Domain(spdim);

  // simdomain_ptr->Add_Geometric_Model(geom_model_ptr);

  // ---------------- MESH -----------------------------------------------
  int rank(0), ierr, aerr;

  // get the Mesh sublist
  Teuchos::ParameterList mesh_parameter_list = plist.sublist("mesh");

  // Create a mesh factory for this geometric model
  Amanzi::AmanziMesh::MeshFactory meshfactory(comm, geom_model_ptr, Teuchos::null);

  // get the Mesh sublist
  ierr = 0;
  Teuchos::ParameterList mesh_params = plist.sublist("mesh");

  Teuchos::ParameterList unstr_mesh_params = mesh_params.sublist("unstructured");

  // Decide on which mesh framework to use
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

  if (unstr_mesh_params.isSublist("generate mesh")) { // If Read parameters are specified
    Teuchos::ParameterList gen_params = unstr_mesh_params.sublist("generate mesh");
    ierr = 0;

    try {
      // create the mesh by internal generation
      mesh = meshfactory.create(gen_params);

    } catch (const std::exception& e) {
      std::cerr << rank << ": error: " << e.what() << std::endl;
      ierr++;
    }

    comm->SumAll(&ierr, &aerr, 1);
    if (aerr > 0) { exit(-aerr); }

  } else { // If Generate parameters are specified
    std::cerr << rank << ": error: "
              << "Neither Read nor Generate options specified for mesh" << std::endl;
    throw std::exception();
  }

  AMANZI_ASSERT(!mesh.is_null());


  // create dummy observation data object
  Amanzi::ObservationData obs_data;
  Teuchos::RCP<Teuchos::ParameterList> glist_rcp = Teuchos::rcp(new Teuchos::ParameterList(plist));

  Teuchos::ParameterList state_plist = glist_rcp->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  Amanzi::CycleDriver cycle_driver(glist_rcp, S, comm, obs_data);

  S = cycle_driver.Go();

  double dt_last;
  dt_last = cycle_driver.get_dt();
  int cycle = S->get_cycle();


  //std::cout<< cycle << " "<<dt_last<<"\n";

  CHECK(cycle == 200);
  CHECK(fabs(dt_last - 1980.32) < 1);
}

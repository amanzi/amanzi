#include <iostream>
#include "stdlib.h"
#include "math.h"

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

#include "CycleDriver.hh"
#include "eos_registration.hh"
#include "MeshFactory.hh"
#include "Mesh.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_flow_registration.hh"
#include "State.hh"
#include "wrm_flow_registration.hh"


TEST(MPC_DRIVER_FLOW) {

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  
  // read the main parameter list
  std::string xmlFileName = "test/mpc_driver_flow.xml";
  Teuchos::RCP<Teuchos::ParameterList> glist = Teuchos::getParametersFromXmlFile(xmlFileName);
  glist->sublist("cycle driver").sublist("time periods")
        .sublist("TP 0").set<double>("end period time", 1.0); 

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = glist->sublist("regions");
  Teuchos::RCP<AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new AmanziGeometry::GeometricModel(2, region_list, &comm));

  // create mesh
  MeshFactory meshfactory(&comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<AmanziMesh::Mesh> mesh = meshfactory(0.0, 0.0, 216.0, 120.0, 27, 15, gm);

  // use cycle driver to create and initialize state
  ObservationData obs_data;    
  CycleDriver cycle_driver(glist, mesh, &comm, obs_data);
  auto S = cycle_driver.Go();

  // overwrite flow & pressure
  AmanziGeometry::Point vel(1.0, 2.0);
  int nnodes = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  auto& flow = *S->GetFieldData("darcy_flux", "flow")->ViewComponent("face");
  auto& pres = *S->GetFieldData("pressure", "flow")->ViewComponent("cell");

  for (int f = 0; f < nfaces; ++f) {
    flow[0][f] = vel * mesh->face_normal(f);
  }

  pres.PutScalar(1.0);

  // verify recovered velocity
  glist->set<std::string>("file name base", "walkabout")
        .set<int>("file name digits", 5);
  auto walkabout = Teuchos::rcp(new Amanzi::WalkaboutCheckpoint(*glist, &comm));

  std::vector<AmanziGeometry::Point> xyz, velocity;
  walkabout->CalculateDarcyVelocity(S, xyz, velocity);

  for (int v = 0; v < nnodes; ++v) {
    CHECK(norm(vel - velocity[v]) < 1e-10);
  }

  // verify recovered pressure
  std::vector<int> material_ids;
  std::vector<double> porosity, saturation, pressure, isotherm_kd;

  walkabout->CalculateData(S, xyz, velocity,
                           porosity, saturation, pressure, isotherm_kd, material_ids);

  for (int v = 0; v < nnodes; ++v) {
    CHECK_CLOSE(1.0, pressure[v], 1e-10);
  }

  // create walkabout file
  walkabout->disable(false);
  walkabout->WriteWalkabout(S);
}



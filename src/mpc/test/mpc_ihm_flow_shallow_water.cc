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
#include "eos_registration.hh"
#include "Mesh.hh"
#include "MeshExtractedManifold.hh"
#include "MeshFactory.hh"
#include "mpc_pks_registration.hh"
#include "numerical_flux_registration.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_flow_registration.hh"
#include "pks_shallow_water_registration.hh"
#include "State.hh"
#include "wrm_flow_registration.hh"

// General
#define _USE_MATH_DEFINES
#include "math.h"


TEST(MPC_DRIVER_IHM_FLOW_SHALLOW_WATER_DAM_BREAK) {

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  
  std::string xmlInFileName = "test/mpc_ihm_flow_shallow_water.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);
  
  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));
  
  // create mesh
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({Framework::MSTK}));
  auto mesh = factory.create(0.0, 0.0, 0.0, 10.0, 1.0, 1.0, 80, 1, 40, true, true);
  
  // deform mesh (if needed)
  AmanziMesh::Entity_ID_List nodeids;
  AmanziGeometry::Point_List new_positions, final_positions;
  int nnodes = mesh->num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::OWNED);
  
  for (int n = 0; n < nnodes; ++n) {
    nodeids.push_back(n);
    
    AmanziGeometry::Point node_crd;
    mesh->node_get_coordinates(n, &node_crd);
    double x = node_crd[0], y = node_crd[1], z = node_crd[2];
    
    // if (std::abs(z - 10) < 1.e-12) {
    //   node_crd[2] += 0.0;
    // }
    new_positions.push_back(node_crd);
  }
  mesh->deform(nodeids, new_positions, false, &final_positions);
  
  // create dummy observation data object
  Amanzi::ObservationData obs_data;
  
  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);
  
  // create additional mesh for SW
  std::vector<std::string> names;
  names.push_back("surface");
  
  auto mesh_surface = factory.create(mesh, { "TopSurface" }, AmanziMesh::FACE, true, true, true);
  
  S->RegisterMesh("surface", mesh_surface);
  
  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  cycle_driver.Go();
  
  // calculate the fluid pressure at the top of the subsurface at final time
  const Epetra_MultiVector &p = *S->GetFieldData("pressure")->ViewComponent("face");
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  double p_top_avg = 0.0, top_surface_area = 0.0;
  
  for (int f = 0; f < nfaces; ++f) {
    const AmanziGeometry::Point &xf = mesh->face_centroid(f);
    if (std::abs(xf[2] - 1.0) < 1.e-12) {
      p_top_avg += p[0][f] * mesh->face_area(f);
      top_surface_area += mesh->face_area(f);
    }
  }
  
  p_top_avg = p_top_avg/top_surface_area;
  
  std::cout<<"average subsurface pressure at surface: "<<p_top_avg<<std::endl;
  
  // calculate the average ponded depth at final time
  double h_avg;
  S->GetFieldData("surface-ponded_depth")->MeanValue(&h_avg);
  const double rho = *S->GetScalarData("const_fluid_density");
  const double patm = *S->GetScalarData("atmospheric_pressure");
  double tmp[1];
  S->GetConstantVectorData("gravity")->Norm2(tmp);
  double g = tmp[0];
  
  std::cout<<"average surface ponded height: "<<h_avg<<std::endl;
  std::cout<<"average hydrostatic pressure at surface: "<<patm + rho*g*h_avg<<std::endl;
  
  // compare the values
  CHECK(std::abs(p_top_avg - (patm + rho*g*h_avg))/(patm + rho*g*h_avg) < 1.0e-2);
}




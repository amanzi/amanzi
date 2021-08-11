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
#include "mpc_pks_registration.hh"
#include "numerical_flux_registration.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_flow_registration.hh"
#include "pks_shallow_water_registration.hh"
#include "State.hh"
#include "wrm_flow_registration.hh"
#include "InputAnalysis.hh"
#include "LeastSquare.hh"
#include "exceptions.hh"
#include "VerboseObject.hh"

// General
#define _USE_MATH_DEFINES
#include "math.h"

double pi = 3.141592653589793;

TEST(MPC_DRIVER_IHM_FLOW_SHALLOW_WATER_TEST) {

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  
  // manufactured solution for the surface-subsurface coupling
  std::string xmlInFileName = "test/mpc_ihm_flow_shallow_water_test_2.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);
  
  std::vector<double> L_inf_h, L1_h, L_inf_u, L1_u, L_inf_p, L2_p, h_max;
  
  for (int NN = 40; NN <= 80; NN = NN + 20) {
  
    // For now create one geometric model from all the regions in the spec
    Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
    auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));
    
    // create mesh
    auto mesh_list = Teuchos::sublist(plist, "mesh", true);
    MeshFactory factory(comm, gm, mesh_list);
    factory.set_preference(Preference({Framework::MSTK}));
    auto mesh = factory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, NN, 1, NN, true, true);
    
    // deform mesh
    AmanziMesh::Entity_ID_List nodeids;
    AmanziGeometry::Point_List new_positions, final_positions;
    int nnodes = mesh->num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::OWNED);
    
    for (int n = 0; n < nnodes; ++n) {
      nodeids.push_back(n);
      
      AmanziGeometry::Point node_crd;
      mesh->node_get_coordinates(n, &node_crd);
      double x = node_crd[0], y = node_crd[1], z = node_crd[2];
      
      if (std::abs(z - 10) < 1.e-12) {
        //      node_crd[2] += (1.0/200)*std::max(0.0, 100 - ((x - 50) * (x - 50) + (y - 50) * (y - 50)));
        node_crd[2] += 0.0;
        
      }
      new_positions.push_back(node_crd);
    }
    mesh->deform(nodeids, new_positions, false, &final_positions);
    
    // create dummy observation data object
    Amanzi::ObservationData obs_data;
    
    Teuchos::ParameterList state_plist = plist->sublist("state");
    Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
    S->RegisterMesh("domain", mesh);
    
    //create additional mesh for SW
    std::vector<std::string> names;
    names.push_back("surface");
    
    //   auto mesh_surface = Teuchos::rcp(new MeshExtractedManifold(
    //       mesh, "TopSurface", AmanziMesh::FACE, comm, gm, mesh_list, true, true, true));
    auto mesh_surface = factory.create(mesh, { "TopSurface" }, AmanziMesh::FACE, true, true, true);
    
    S->RegisterMesh("surface", mesh_surface);
    
    Teuchos::RCP<TreeVector> soln = Teuchos::rcp (new TreeVector());
    
    Teuchos::ParameterList sw_list = plist->sublist("PKs").sublist("shallow water");
    Teuchos::ParameterList flow_list = plist->sublist("PKs").sublist("flow");
    
    
    Amanzi::InputAnalysis analysis(mesh, "domain");
    analysis.Init(*plist);
    analysis.RegionAnalysis();
    analysis.OutputBCs();
    
    Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
    cycle_driver.Go();
    
    const Epetra_MultiVector &h = *S->GetFieldData("surface-ponded_depth")->ViewComponent("cell");
    const Epetra_MultiVector &vel = *S->GetFieldData("surface-velocity")->ViewComponent("cell");
    const Epetra_MultiVector &p = *S->GetFieldData("pressure")->ViewComponent("cell");
    const Epetra_MultiVector &flux = *S->GetFieldData("darcy_flux")->ViewComponent("face", true);
    const double rho = *S->GetScalarData("const_fluid_density");
    
    double time = S->time();
    
    std::cout<<"end time: "<<time<<std::endl;

    // calculate ponded_depth and velocity error
    double err_max_h = 0.0, err_L1_h = 0.0, err_max_u = 0.0, err_L1_u = 0.0;
    int ncells_owned = mesh_surface->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
    std::cout<<ncells_owned<<std::endl;
    for (int c = 0; c < ncells_owned; ++c) {
      const Amanzi::AmanziGeometry::Point &xc = mesh_surface -> cell_centroid(c);
      double x = xc[0], y = xc[1];
//      err_max_h = std::max(err_max_h, std::abs( h[0][c] - ( 1.0 + std::sin(pi*xc[0])*std::exp(-time) ) ));
//      err_L1_h += (std::abs( h[0][c] - ( 1.0 + std::sin(pi*xc[0])*std::exp(-time) ) )) * mesh_surface->cell_volume(c);
//
//      err_max_u = std::max(err_max_u, std::abs(vel[0][c] - (0.0) ) );
//      err_L1_u += std::abs(vel[0][c] - ( 0.0 )) * mesh_surface->cell_volume(c);
      
      err_max_h = std::max(err_max_h, std::abs( h[0][c] - ( 1.0 + x ) ));
      err_L1_h += (std::abs( h[0][c] - ( 1.0 + x ) ) ) * mesh_surface->cell_volume(c);

      err_max_u = std::max(err_max_u, std::abs(vel[0][c] - (0.0) ) );
      err_L1_u += std::abs(vel[0][c] - ( 0.0 )) * mesh_surface->cell_volume(c);
    }
    L_inf_h.push_back(err_max_h);
    L1_h.push_back(err_L1_h);
    L_inf_u.push_back(err_max_u);
    L1_u.push_back(err_L1_u);
    std::cout<<"Computed error h (L_1): "<<err_L1_h<<std::endl;
    std::cout<<"Computed error h (L_inf): "<<err_max_h<<std::endl;
    std::cout<<"Computed error u (L_1): "<<err_L1_u<<std::endl;
    std::cout<<"Computed error u (L_inf): "<<err_max_u<<std::endl;
    
    // calculate pressure error
    double err_max_p = 0.0, err_L2_p = 0.0;
    ncells_owned = mesh->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
    std::cout<<ncells_owned<<std::endl;
    for (int c = 0; c < ncells_owned; ++c) {
      const Amanzi::AmanziGeometry::Point &xc = mesh -> cell_centroid(c);
      double x = xc[0], y = xc[1], z = xc[2];
//      err_max_p = std::max(err_max_p, std::abs(p[0][c] - (1.0 + std::exp(-1.0*time)*std::sin(pi*z) + std::exp(-1.0*time)*std::sin(pi*x) ) ) );
//      err_L2_p += std::pow( std::abs(p[0][c] - ( 1.0 + std::exp(-1.0*time)*std::sin(pi*z) + std::exp(-1.0*time)*std::sin(pi*x) ) ), 2) * mesh->cell_volume(c);
      
      err_max_p = std::max(err_max_p, std::abs(p[0][c] - ( 3.0 + x - 2.0*z ) ) );
      err_L2_p += std::pow( std::abs(p[0][c] - ( 3.0 + x - 2.0*z ) ), 2) * mesh->cell_volume(c);
    }
    err_L2_p = std::sqrt(err_L2_p);
    L2_p.push_back(err_L2_p);
    L_inf_p.push_back(err_max_p);
    h_max.push_back(1.0/NN);
    std::cout<<"Computed error p (L_2): "<<err_L2_p<<std::endl;
    std::cout<<"Computed error p (L_inf): "<<err_max_p<<std::endl;
    
  }
  
  std::cout<<"h error L_inf: "<<L_inf_h[0]<<" : "<<L_inf_h[1]<<" : "<<L_inf_h[2]<<std::endl;
  std::cout<<"h error L_1: "<<L1_h[0]<<" : "<<L1_h[1]<<" : "<<L1_h[2]<<std::endl;
  std::cout<<"u error L_inf: "<<L_inf_u[0]<<" : "<<L_inf_u[1]<<" : "<<L_inf_u[2]<<std::endl;
  std::cout<<"u error L_1: "<<L1_u[0]<<" : "<<L1_u[1]<<" : "<<L1_u[2]<<std::endl;
  std::cout<<"p error L_inf: "<<L_inf_p[0]<<" : "<<L_inf_p[1]<<" : "<<L_inf_p[2]<<std::endl;
  std::cout<<"p error L_2: "<<L2_p[0]<<" : "<<L2_p[1]<<" : "<<L2_p[2]<<std::endl;
    
  double p_rate = Amanzi::Utils::bestLSfit(h_max, L2_p);
  double p_rate_L_inf = Amanzi::Utils::bestLSfit(h_max, L_inf_p);
  double h_rate_L_inf = Amanzi::Utils::bestLSfit(h_max, L_inf_h);
  double h_rate_L1 = Amanzi::Utils::bestLSfit(h_max, L1_h);
  double u_rate_L_inf = Amanzi::Utils::bestLSfit(h_max, L_inf_u);
  double u_rate_L1 = Amanzi::Utils::bestLSfit(h_max, L1_u);
  
  std::cout<<"order of convergence p (L_inf): "<<p_rate_L_inf<<std::endl;
  std::cout<<"order of convergence p (L_2): "<<p_rate<<std::endl;
  std::cout<<"order of convergence h (L_inf): "<<h_rate_L_inf<<std::endl;
  std::cout<<"order of convergence h (L_1): "<<h_rate_L1<<std::endl;
  std::cout<<"order of convergence u (L_inf): "<<u_rate_L_inf<<std::endl;
  std::cout<<"order of convergence u (L_1): "<<u_rate_L1<<std::endl;
  
}



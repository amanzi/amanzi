/*
  Shallow water PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
*/

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "IO.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "OutputXDMF.hh"

// Amanzi::ShallowWater
#include "ShallowWater_PK.hh"
#include "ShallowWater_Helper.hh"

// General
#define _USE_MATH_DEFINES
#include "math.h"

using namespace Amanzi;

// ---- ---- ---- ---- ---- ---- ---- ----
// Inital Conditions
// ---- ---- ---- ---- ---- ---- ---- ----
void
dry_bed_setIC(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
              Teuchos::RCP<Amanzi::State>& S, int icase)
{
  double pi = M_PI;
  const double rho = S->Get<double>("const_fluid_density");
  const double patm = S->Get<double>("atmospheric_pressure");

  double g = norm(S->Get<AmanziGeometry::Point>("gravity"));

  int ncells_wghost = mesh->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::ALL);
  int nnodes_wghost = mesh->num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::ALL);
  std::string passwd = "state";

  auto& B_c = *S->GetW<CompositeVector>("surface-bathymetry", Tags::DEFAULT, passwd).ViewComponent("cell");
  auto& B_n = *S->GetW<CompositeVector>("surface-bathymetry", Tags::DEFAULT, passwd).ViewComponent("node");
  auto& h_c = *S->GetW<CompositeVector>("surface-ponded_depth", Tags::DEFAULT, passwd).ViewComponent("cell");
  auto& ht_c = *S->GetW<CompositeVector>("surface-total_depth", Tags::DEFAULT, passwd).ViewComponent("cell");
  auto& vel_c = *S->GetW<CompositeVector>("surface-velocity", Tags::DEFAULT, passwd).ViewComponent("cell");
  auto& q_c = *S->GetW<CompositeVector>("surface-discharge", Tags::DEFAULT, "surface-discharge").ViewComponent("cell");
  auto& p_c = *S->GetW<CompositeVector>("surface-ponded_pressure", Tags::DEFAULT,"surface-ponded_pressure").ViewComponent("cell");

  // Define bathymetry at the cell vertices (Bn)
  for (int n = 0; n < nnodes_wghost; ++n) {
    Amanzi::AmanziGeometry::Point node_crd;
    
    mesh->node_get_coordinates(n, &node_crd); // Coordinate of current node
    
    double x = node_crd[0], y = node_crd[1];
    
    B_n[0][n] = 0.0;
    if ((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5) < 0.35*0.35 + 1.e-12) {
      B_n[0][n] = 0.8;
    }
  }
  
  S->Get<CompositeVector>("surface-bathymetry").ScatterMasterToGhosted("node");
  
  for (int c = 0; c < ncells_wghost; ++c) {
    const Amanzi::AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    
    Amanzi::AmanziMesh::Entity_ID_List cfaces, cnodes, cedges;
    mesh->cell_get_faces(c, &cfaces);
    mesh->cell_get_nodes(c, &cnodes);
    mesh->cell_get_edges(c, &cedges);
    
    int nfaces_cell = cfaces.size();
    
    B_c[0][c] = 0.0;
    
    // Compute cell averaged bathymetry (Bc)
    for (int f = 0; f < nfaces_cell; ++f) {
      Amanzi::AmanziGeometry::Point x0, x1;
      int edge = cfaces[f];
      
      Amanzi::AmanziMesh::Entity_ID_List face_nodes;
      mesh->face_get_nodes(edge, &face_nodes);
      int n0 = face_nodes[0], n1 = face_nodes[1];
      
      mesh->node_get_coordinates(n0, &x0);
      mesh->node_get_coordinates(n1, &x1);
      
      double area = norm((xc - x0) ^ (xc - x1)) / 2.0;
      
      B_c[0][c] += (area / mesh->cell_volume(c)) * (B_n[0][n0] + B_n[0][n1]) / 2.0;
    }
    
    // start with dry bed
    ht_c[0][c] = std::max(0.0, B_c[0][c]);
    
    h_c[0][c] = ht_c[0][c] - B_c[0][c];
    vel_c[0][c] = 0.0;
    vel_c[1][c] = 0.0;
    q_c[0][c] = h_c[0][c] * vel_c[0][c];
    q_c[1][c] = h_c[0][c] * vel_c[1][c];
    p_c[0][c] = patm + rho * g * h_c[0][c];
  }
  
  S->Get<CompositeVector>("surface-bathymetry").ScatterMasterToGhosted("cell");
  S->Get<CompositeVector>("surface-total_depth").ScatterMasterToGhosted("cell");
  S->Get<CompositeVector>("surface-ponded_depth").ScatterMasterToGhosted("cell");
  S->Get<CompositeVector>("surface-velocity").ScatterMasterToGhosted("cell");
  S->Get<CompositeVector>("surface-discharge").ScatterMasterToGhosted("cell");
}


/* **************************************************************** */
void
RunTest(int icase)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::ShallowWater;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  int MyPID = comm->MyPID();

  if (MyPID == 0) {
    std::cout
      << "Test: Shallow water: Flow over dry bed with B > 0"
      << std::endl;
  }

  // Read parameter list
  std::string xmlFilename;
  xmlFilename = "test/shallow_water_dry_bed.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFilename);

  // Create a mesh framework
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));

  // Creat a mesh
  bool request_faces = true, request_edges = false;
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));

  RCP<Mesh> mesh;
  if (icase == 1) {
    mesh = meshfactory.create ("test/triangular16.exo");
  } else if (icase == 2) {
    mesh = meshfactory.create ("test/median15x16.exo");
  }

  // Create a state
  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterMesh("surface", mesh);

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());

  Teuchos::ParameterList sw_list = plist->sublist("PKs").sublist("shallow water");

  // Create a shallow water PK
  ShallowWater_PK SWPK(sw_list, plist, S, soln);
  SWPK.Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  S->set_time(0.0);
  SWPK.Initialize();

  dry_bed_setIC(mesh, S, icase);
  S->CheckAllFieldsInitialized();

  const auto& B = *S->Get<CompositeVector>("surface-bathymetry").ViewComponent("cell");
  const auto& Bn = *S->Get<CompositeVector>("surface-bathymetry").ViewComponent("node");
  const auto& hh = *S->Get<CompositeVector>("surface-ponded_depth").ViewComponent("cell");
  const auto& ht = *S->Get<CompositeVector>("surface-total_depth").ViewComponent("cell");
  const auto& vel = *S->Get<CompositeVector>("surface-velocity").ViewComponent("cell");
  const auto& p = *S->Get<CompositeVector>("surface-ponded_pressure").ViewComponent("cell");

  // Create screen io
  auto vo = Teuchos::rcp(new Amanzi::VerboseObject("ShallowWater", *plist));
  WriteStateStatistics(*S, *vo);

  // Advance in time
  double t_old(0.0), t_new(0.0), dt;

  // Initialize io
  Teuchos::ParameterList iolist;
  std::string fname;
  fname = "SW_sol";
  iolist.get<std::string>("file name base", fname);
  OutputXDMF io(iolist, mesh, true, false);

  std::string passwd("state");

  int iter = 0;
  std::vector<double> dx, Linferror, L1error, L2error;

  double Tend = 0.5;

  int ncells_wghost = mesh->num_entities(
    Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::ALL);
  
  while ((t_new < Tend) && (iter >= 0)) {
    double t_out = t_new;

    if (iter % 100 == 0) {
      IO_Fields(t_out, iter, MyPID, io, *S, nullptr, nullptr);
    }

    dt = SWPK.get_dt();

    t_new = t_old + dt;

    SWPK.AdvanceStep(t_old, t_new);
    SWPK.CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;
    iter += 1;
    
    // check depth positivity
    for (int c = 0; c < ncells_wghost; ++c) {
      CHECK(hh[0][c] >= 0.0);
    }
  } // time loop

  if (MyPID == 0) { std::cout << "Time-stepping finished. " << std::endl; }
  std::cout<<"current time: "<<t_new<<", dt = "<<dt<<std::endl;
  
  double t_out = t_new;
  IO_Fields(t_out, iter, MyPID, io, *S, nullptr, nullptr);
}

TEST(SHALLOW_WATER_DRY_BED)
{
  RunTest(1); // triangular mesh
  RunTest(2); // hexagonal mesh
}

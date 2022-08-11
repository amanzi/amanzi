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
#include "ShallowWater_PK.hh"
#include "OutputXDMF.hh"

// General
#define _USE_MATH_DEFINES
#include "math.h"

using namespace Amanzi;

double mass_eq = 0.0;

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

  int ncells_wghost = mesh->num_entities(
    Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::ALL);
  int nnodes_wghost = mesh->num_entities(
    Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::ALL);
  std::string passwd = "state";

  auto& B_c =
    *S->GetW<CompositeVector>("surface-bathymetry", Tags::DEFAULT, passwd)
       .ViewComponent("cell");
  auto& B_n =
    *S->GetW<CompositeVector>("surface-bathymetry", Tags::DEFAULT, passwd)
       .ViewComponent("node");
  auto& h_c =
    *S->GetW<CompositeVector>("surface-ponded_depth", Tags::DEFAULT, passwd)
       .ViewComponent("cell");
  auto& ht_c =
    *S->GetW<CompositeVector>("surface-total_depth", Tags::DEFAULT, passwd)
       .ViewComponent("cell");
  //need to implement auto& ht_n = *S->GetW<CompositeVector>("surface-total_depth", Tags::DEFAULT, passwd).ViewComponent("cell");

  auto& vel_c =
    *S->GetW<CompositeVector>("surface-velocity", Tags::DEFAULT, passwd)
       .ViewComponent("cell");
  auto& q_c = *S->GetW<CompositeVector>(
                  "surface-discharge", Tags::DEFAULT, "surface-discharge")
                 .ViewComponent("cell");
  auto& p_c = *S->GetW<CompositeVector>("surface-ponded_pressure",
                                        Tags::DEFAULT,
                                        "surface-ponded_pressure")
                 .ViewComponent("cell");

  // Define bathymetry at the cell vertices (Bn)
  for (int n = 0; n < nnodes_wghost; ++n) {
    Amanzi::AmanziGeometry::Point node_crd;

    mesh->node_get_coordinates(n, &node_crd); // Coordinate of current node

    double x = node_crd[0], y = node_crd[1];
      
    B_n[0][n] = 0.0;
    if (std::abs(x - 0.5) <= 0.3 && std::abs(y - 0.5) <= 0.3) {
      B_n[0][n] = 0.5 * std::sin(pi * x) * std::sin(pi * y);
    }
    
    //B_n[0][n] = std::sin(pi * x) * std::sin(pi * y);

    //B_n[0][n] = x;
    if (x >= 0.3) { 
     // B_n[0][n] = 0.0 * std::sin(0.5 * pi * x);
     // B_n[0][n] = x;
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
    h_c[0][c] = 0.0;
    ht_c[0][c] = 0.0;

    // Compute cell averaged bathymetrt (Bc)
    for (int f = 0; f < nfaces_cell; ++f) {
      Amanzi::AmanziGeometry::Point x0, x1;
      int edge = cfaces[f];

      Amanzi::AmanziMesh::Entity_ID_List face_nodes;
      mesh->face_get_nodes(edge, &face_nodes);

      mesh->node_get_coordinates(face_nodes[0], &x0);
      mesh->node_get_coordinates(face_nodes[1], &x1);

      Amanzi::AmanziGeometry::Point tria_edge0, tria_edge1;

      tria_edge0 = xc - x0;
      tria_edge1 = xc - x1;

      Amanzi::AmanziGeometry::Point area_cross_product =
        (0.5) * tria_edge0 ^ tria_edge1;

      double area = norm(area_cross_product);

      B_c[0][c] += (area / mesh->cell_volume(c)) *
                   (B_n[0][face_nodes[0]] + B_n[0][face_nodes[1]]) / 2.0;

      ht_c[0][c] += (area / mesh->cell_volume(c)) * (std::max(0.5, B_n[0][face_nodes[0]])  + std::max(0.5, B_n[0][face_nodes[1]]) ) / 2.0;
    }

    ht_c[0][c] = std::max(0.5, B_c[0][c]);
    
    //ht_c[0][c] = h_c[0][c] + B_c[0][c];
    
    if ((xc[0] - 0.2)*(xc[0] - 0.2) + (xc[1] - 0.5)*(xc[1] - 0.5) < 0.1*0.1) {
      ht_c[0][c] += 0.1;
    }
    h_c[0][c] = ht_c[0][c] - B_c[0][c];
    
    vel_c[0][c] = 0.0;
    vel_c[1][c] = 0.0;
    q_c[0][c] = h_c[0][c] * vel_c[0][c];
    q_c[1][c] = h_c[0][c] * vel_c[1][c];
    p_c[0][c] = patm + rho * g * h_c[0][c];
  }
  

  // triangular mesh only; well-balanced IC [Liu et al.' 18]
  for (int c = 0; c < ncells_wghost; ++c) {
    AmanziMesh::Entity_ID_List cnodes;
    mesh->cell_get_nodes(c, &cnodes);
   
    double B13, B12, B23; 
   
    B13 = std::max( std::max(B_n[0][cnodes[0]], B_n[0][cnodes[1]]), B_n[0][cnodes[2]] );
    B23 = std::min( std::min(B_n[0][cnodes[0]], B_n[0][cnodes[1]]), B_n[0][cnodes[2]] );
    B12 = (B_n[0][cnodes[0]] + B_n[0][cnodes[1]] + B_n[0][cnodes[2]] ) - (B13 + B23); 
    
    Amanzi::AmanziGeometry::Point xv13, xv12, xv23;
    int i13, i12, i23;
    for (int i = 0; i < cnodes.size(); ++i) {
       if(std::abs(B13 - B_n[0][cnodes[i]]) < 1.e-15 ) {
         mesh->node_get_coordinates(cnodes[i], &xv13);
         i13 = cnodes[i];
       }
     }
     
     for (int i = 0; i < cnodes.size(); ++i) {
       if (std::abs(B12 - B_n[0][cnodes[i]]) < 1.e-15 && i13 != cnodes[i]) {
         mesh->node_get_coordinates(cnodes[i], &xv12);
         i12 = cnodes[i];
       }
     }
 
     for (int i = 0; i < cnodes.size(); ++i) {
       if (cnodes[i] != i13 && cnodes[i] != i12) {
         mesh->node_get_coordinates(cnodes[i], &xv23);
         i23 = cnodes[i];
       }
     }
 
 
     // store sorted nodes and coordinates
     std::vector<double> Bi(3);
     Bi[0] = B13;
     Bi[1] = B12;
     Bi[2] = B23;
 
     std::vector<AmanziGeometry::Point> xvi(3);
     xvi[0] = xv13;
     xvi[1] = xv12;
     xvi[2] = xv23;

     double w = 0.8; // free water surface
    
     // fully flooded cell
     if (w >= B13 && w >= B12 && w >= B23) {
       h_c[0][c] = (w - B_c[0][c]);
     } else if (w < B23) {
       h_c[0][c] = 0.0;
     } else {
       // type 1 cell
       if (w >= B23 && w < B12 ) {
         h_c[0][c] = std::pow((w - B23), 3.0) / (3.0 * (B13 - B23) * (B12 - B23)); 
       } else if (w >= B12 && w < B13) {
         h_c[0][c] = ( -w*w*w + 3.0*B13*w*w - 3.0*(B23*B13 + B12*B13 - B12*B23)*w + B23*B23*B13) / (B13 - B23) + B12*(B12 + B23);
         h_c[0][c] /= 3.0*(B13 - B12);
       }
     }

     ht_c[0][c] = h_c[0][c] + B_c[0][c];

     ht_c[0][c] = std::max(0.0, B_c[0][c]);
     
     //ht_c[0][c]  = std::max(0.8, B_c[0][c]);
     const Amanzi::AmanziGeometry::Point &xc = mesh->cell_centroid(c);
     if ((xc[0] - 0.1)*(xc[0] - 0.1) + (xc[1] - 0.1)*(xc[1] - 0.1) < 0.05*0.05) {
       ht_c[0][c] += 0.0;
     }
     if (c == 179) {
      ht_c[0][c] += 0.00;
      std::cout<<"cell = "<<c<<"; center = "<<xc[0]<<", "<<xc[1]<<std::endl;
     }
     h_c[0][c] = ht_c[0][c] - B_c[0][c];
      
      double weq = 0.80, hceq; // equilibrium free water surface
      // fully flooded cell
      if (weq >= B13 && weq >= B12 && weq >= B23) {
        hceq = (weq - B_c[0][c]);
      } else if (weq < B23) {
        hceq = 0.0;
      } else {
        // type 1 cell
        if (weq >= B23 && weq < B12 ) {
          hceq = std::pow((weq - B23), 3.0) / (3.0 * (B13 - B23) * (B12 - B23));
        } else if (weq >= B12 && weq < B13) {
          hceq = ( -weq*weq*weq + 3.0*B13*weq*weq - 3.0*(B23*B13 + B12*B13 - B12*B23)*weq + B23*B23*B13) / (B13 - B23) + B12*(B12 + B23);
          hceq /= 3.0*(B13 - B12);
       }
      }

      mass_eq += hceq * mesh->cell_volume(c);

 
  }
  
  S->Get<CompositeVector>("surface-bathymetry").ScatterMasterToGhosted("cell");
  S->Get<CompositeVector>("surface-total_depth").ScatterMasterToGhosted("cell");
  S->Get<CompositeVector>("surface-ponded_depth")
    .ScatterMasterToGhosted("cell");
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
      << "Test: Shallow water: Dry bed well balancing"
      << std::endl;
  }

  // Read parameter list
  std::string xmlFilename;
  
  xmlFilename = "test/shallow_water_dry_bed_wb.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist =
    Teuchos::getParametersFromXmlFile(xmlFilename);

  // Create a mesh framework
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(
    new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));

  // Creat a mesh
  bool request_faces = true, request_edges = false;
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));

  RCP<Mesh> mesh;
  //mesh = meshfactory.create ("test/triangular16.exo");
  mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 50, 50, request_faces, request_edges);
  //mesh = meshfactory.create ("test/median32x33.exo");

  // Other polygonal meshes
  //  RCP<Mesh> mesh = meshfactory.create ("test/median15x16.exo");
  //  RCP<Mesh> mesh = meshfactory.create ("test/random40.exo");

  // Create a state

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterMesh("surface", mesh);

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());

  Teuchos::ParameterList sw_list =
    plist->sublist("PKs").sublist("shallow water");

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

  int ncells_wghost = mesh->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::ALL);\

  const auto& B =
    *S->Get<CompositeVector>("surface-bathymetry").ViewComponent("cell");
  const auto& Bn =
    *S->Get<CompositeVector>("surface-bathymetry").ViewComponent("node");
  const auto& hh =
    *S->Get<CompositeVector>("surface-ponded_depth").ViewComponent("cell");
  const auto& ht =
    *S->Get<CompositeVector>("surface-total_depth").ViewComponent("cell");
  const auto& vel =
    *S->Get<CompositeVector>("surface-velocity").ViewComponent("cell");
  const auto& q =
    *S->Get<CompositeVector>("surface-discharge").ViewComponent("cell");
  const auto& p =
    *S->Get<CompositeVector>("surface-ponded_pressure").ViewComponent("cell");

  // Create a pid vector
  Epetra_MultiVector pid(B);

  for (int c = 0; c < pid.MyLength(); ++c) { pid[0][c] = MyPID; }

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

  double Tend;
  //Tend = 50000.0;
  Tend = 20.0;
  
  while ((t_new < Tend) && (iter >= 0)) {
    double t_out = t_new;

    if (iter % 1000 == 0) {
      io.InitializeCycle(t_out, iter, "");

      io.WriteVector(*hh(0), "depth", AmanziMesh::CELL);
      io.WriteVector(*ht(0), "total_depth", AmanziMesh::CELL);
      io.WriteVector(*vel(0), "vx", AmanziMesh::CELL);
      io.WriteVector(*vel(1), "vy", AmanziMesh::CELL);
      io.WriteVector(*q(0), "qx", AmanziMesh::CELL);
      io.WriteVector(*q(1), "qy", AmanziMesh::CELL);
      io.WriteVector(*B(0), "B", AmanziMesh::CELL);
      io.WriteVector(*p(0), "hyd pressure", AmanziMesh::CELL);
      io.WriteVector(*Bn(0), "B_n", AmanziMesh::NODE);
      io.WriteVector(*pid(0), "pid", AmanziMesh::CELL);

      io.FinalizeCycle();
    }

    dt = SWPK.get_dt();

    t_new = t_old + dt;

    SWPK.AdvanceStep(t_old, t_new);
    SWPK.CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;
    iter += 1;
    
    if (iter % 1000 == 0) {
    	std::cout<<"current time: "<<t_new<<", dt = "<<dt<<std::endl;
      double mass_end = 0.0;
      for (int c = 0; c < ncells_wghost; ++c) {
        mass_end += hh[0][c] * mesh->cell_volume(c);
      }
      std::cout<<"Equilibrium mass = "<<mass_eq<<"; Actual calculated = "<<mass_end<<"; difference = "<<std::abs(mass_end - mass_eq)<<std::endl;

      std::cout<<" ------------------------------------------------------------------------ "<<std::endl;
    }
  } // time loop

  if (MyPID == 0) { std::cout << "Time-stepping finished. " << std::endl; }
	std::cout<<"current time: "<<t_new<<", dt = "<<dt<<std::endl;
	
  double t_out = t_new;

  io.InitializeCycle(t_out, iter, "");

  io.WriteVector(*hh(0), "depth", AmanziMesh::CELL);
  io.WriteVector(*ht(0), "total_depth", AmanziMesh::CELL);
  io.WriteVector(*vel(0), "vx", AmanziMesh::CELL);
  io.WriteVector(*vel(1), "vy", AmanziMesh::CELL);
  io.WriteVector(*q(0), "qx", AmanziMesh::CELL);
  io.WriteVector(*q(1), "qy", AmanziMesh::CELL);
  io.WriteVector(*B(0), "B", AmanziMesh::CELL);
  io.WriteVector(*Bn(0), "B_n", AmanziMesh::NODE);
  io.WriteVector(*p(0), "hyd pressure", AmanziMesh::CELL);
  io.WriteVector(*pid(0), "pid", AmanziMesh::CELL);

  io.FinalizeCycle();
}

TEST(SHALLOW_WATER_DRY_BED_WB)
{
  RunTest(1);
}

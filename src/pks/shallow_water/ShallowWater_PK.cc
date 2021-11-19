/*
 Shallow water PK
 
 Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
 Amanzi is released under the three-clause BSD License.
 The terms of use and "as is" disclaimer for this license are
 provided in the top-level COPYRIGHT file.
 
 Author: Svetlana Tokareva (tokareva@lanl.gov)
*/

#include <algorithm>
#include <cmath>
#include <vector>

#include "PK_DomainFunctionFactory.hh"

#include "CompositeVector.hh"
#include "Geometry.hh"

// Amanzi::ShallowWater
#include "DischargeEvaluator.hh"
#include "HydrostaticPressureEvaluator.hh"
#include "NumericalFluxFactory.hh"
#include "ShallowWater_PK.hh"

namespace Amanzi {
namespace ShallowWater {
//--------------------------------------------------------------
// Standard constructor
//--------------------------------------------------------------
ShallowWater_PK::ShallowWater_PK(Teuchos::ParameterList& pk_tree,
                                 const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                 const Teuchos::RCP<State>& S,
                                 const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, glist, S, soln),
    S_(S),
    soln_(soln),
    glist_(glist),
    passwd_("state"),
    iters_(0)
{
  std::string pk_name = pk_tree.name();
  auto found = pk_name.rfind("->");
  if (found != std::string::npos) pk_name.erase(0, found + 2);

  // Create miscellaneous lists.
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  sw_list_ = Teuchos::sublist(pk_list, pk_name, true);

  // domain name
  domain_ = sw_list_->template get<std::string>("domain name", "surface");

  cfl_ = sw_list_->get<double>("cfl", 0.1);
  max_iters_ = sw_list_->get<int>("number of reduced cfl cycles", 10);

  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = sw_list_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("ShallowWater", vlist));
}


//--------------------------------------------------------------
// Register fields and field evaluators with the state
// Conservative variables: (h, hu, hv)
//--------------------------------------------------------------
void ShallowWater_PK::Setup(const Teuchos::Ptr<State>& S)
{
  mesh_ = S->GetMesh(domain_);
  dim_ = mesh_->space_dimension();

  // domain name
  velocity_key_ = Keys::getKey(domain_, "velocity");
  discharge_key_ = Keys::getKey(domain_, "discharge");
  ponded_depth_key_ = Keys::getKey(domain_, "ponded_depth");
  total_depth_key_ = Keys::getKey(domain_, "total_depth");
  bathymetry_key_ = Keys::getKey(domain_, "bathymetry");
  hydrostatic_pressure_key_ = Keys::getKey(domain_, "ponded_pressure");

  //-------------------------------
  // constant fields
  //-------------------------------
  if (!S->HasField("gravity")) {
    S->RequireConstantVector("gravity", passwd_, 2);
  }
  
  // required for calculating hydrostatic pressure
  if (!S->HasField("const_fluid_density")) {
    S->RequireScalar("const_fluid_density", passwd_);
  }
  
  if (!S->HasField("atmospheric_pressure")) {
    S->RequireScalar("atmospheric_pressure", passwd_);
  }

  //-------------------------------
  // primary fields
  //-------------------------------

  // ponded_depth_key_
  if (!S->HasField(ponded_depth_key_)) {
    std::vector<std::string> names({"cell", "node"});
    std::vector<int> ndofs(2, 1);
    std::vector<AmanziMesh::Entity_kind> locations({AmanziMesh::CELL, AmanziMesh::NODE});
    
    S->RequireField(ponded_depth_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);
    
//    S->RequireField(ponded_depth_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
//      ->SetComponent("cell", AmanziMesh::CELL, 1);
    AddDefaultPrimaryEvaluator_(ponded_depth_key_);
  }

  // total_depth_key_
  if (!S->HasField(total_depth_key_)) {
    std::vector<std::string> names({"cell", "node"});
    std::vector<int> ndofs(2, 1);
    std::vector<AmanziMesh::Entity_kind> locations({AmanziMesh::CELL, AmanziMesh::NODE});
    
    S->RequireField(total_depth_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);
//    S->RequireField(total_depth_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
//      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // velocity
  if (!S->HasField(velocity_key_)) {
    std::vector<std::string> names({"cell", "node"});
    std::vector<int> ndofs(2, 2);
    std::vector<AmanziMesh::Entity_kind> locations({AmanziMesh::CELL, AmanziMesh::NODE});
    
    S->RequireField(velocity_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);
    
//    S->RequireField(velocity_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
//      ->SetComponent("cell", AmanziMesh::CELL, 2);
    AddDefaultPrimaryEvaluator_(velocity_key_);
  }

  // discharge
  if (!S->HasField(discharge_key_)) {
    std::vector<std::string> names({"cell", "node"});
    std::vector<int> ndofs(2, 2);
    std::vector<AmanziMesh::Entity_kind> locations({AmanziMesh::CELL, AmanziMesh::NODE});
    
    S->RequireField(discharge_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);
//    S->RequireField(discharge_key_, discharge_key_)->SetMesh(mesh_)->SetGhosted(true)
//      ->SetComponent("cell", AmanziMesh::CELL, 2);

    Teuchos::ParameterList elist;
    auto eval = Teuchos::rcp(new DischargeEvaluator(elist));
    S->SetFieldEvaluator(discharge_key_, eval);
  }

  // bathymetry
  if (!S->HasField(bathymetry_key_)) {
    std::vector<std::string> names({"cell", "node"});
    std::vector<int> ndofs(2, 1);
    std::vector<AmanziMesh::Entity_kind> locations({AmanziMesh::CELL, AmanziMesh::NODE});
    
    S->RequireField(bathymetry_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);
  }

  //-------------------------------
  // secondary fields
  //-------------------------------

  // hydrostatic pressure
  if (!S->HasField(hydrostatic_pressure_key_)) {
    S->RequireField(hydrostatic_pressure_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    
    Teuchos::ParameterList elist;
    auto eval = Teuchos::rcp(new HydrostaticPressureEvaluator(elist));
    S->SetFieldEvaluator(hydrostatic_pressure_key_, eval);
  }
}


//--------------------------------------------------------------
// Initialize internal data
//--------------------------------------------------------------
void ShallowWater_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  // Create BC objects
  Teuchos::RCP<ShallowWaterBoundaryFunction> bc;
  Teuchos::RCP<Teuchos::ParameterList>
      bc_list = Teuchos::rcp(new Teuchos::ParameterList(sw_list_->sublist("boundary conditions", false)));

  bcs_.clear();

  // -- velocity
  if (bc_list->isSublist("velocity")) {
    PK_DomainFunctionFactory<ShallowWaterBoundaryFunction > bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("velocity");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);

        bc = bc_factory.Create(spec, "velocity", AmanziMesh::FACE, Teuchos::null);
        bc->set_bc_name("velocity");
        bc->set_type(WhetStone::DOF_Type::VECTOR);
        bcs_.push_back(bc);
      }
    }
  }
  
  // source term
  if (sw_list_->isSublist("source terms")) {
    PK_DomainFunctionFactory<PK_DomainFunction> factory(mesh_, S_);
    auto src_list = sw_list_->sublist("source terms");
    for (auto it = src_list.begin(); it != src_list.end(); ++it) {
      std::string name = it->first;
      if (src_list.isSublist(name)) {
        Teuchos::ParameterList& spec = src_list.sublist(name);
        
        srcs_.push_back(factory.Create(spec, "source", AmanziMesh::CELL, Teuchos::null));
      }
    }
  }

  // gravity
  double tmp[1];
  S_->GetConstantVectorData("gravity", "state")->Norm2(tmp);
  g_ = tmp[0];

  // numerical flux
  Teuchos::ParameterList model_list;
  model_list.set<std::string>("numerical flux", sw_list_->get<std::string>("numerical flux", "central upwind"))
            .set<double>("gravity", g_);
  NumericalFluxFactory nf_factory;
  numerical_flux_ = nf_factory.Create(model_list);

  // default
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  if (!S_->GetField(bathymetry_key_, passwd_)->initialized()) {
    InitializeField_(S_.ptr(), passwd_, bathymetry_key_, 0.0);
  }
  
  const auto& B_n = *S_->GetFieldData(bathymetry_key_)->ViewComponent("node");
  const auto& B_c = *S_->GetFieldData(bathymetry_key_)->ViewComponent("cell");
  
  // compute B_c from B_n for well balanced scheme (Beljadid et. al. 2016)
  S_->GetFieldData(bathymetry_key_)->ScatterMasterToGhosted("node");

  for (int c = 0; c < ncells_owned; ++c) {
    const Amanzi::AmanziGeometry::Point &xc = mesh_->cell_centroid(c);
        
    Amanzi::AmanziMesh::Entity_ID_List cfaces;
    mesh_->cell_get_faces(c, &cfaces);
    int nfaces_cell = cfaces.size();
        
    B_c[0][c] = 0.0;
        
    // compute cell averaged bathymery (Bc)
    for (int f = 0; f < nfaces_cell; ++f) {
    
      Amanzi::AmanziGeometry::Point x0, x1;
      int edge = cfaces[f];
            
      AmanziMesh::Entity_ID_List face_nodes;
      mesh_->face_get_nodes(edge, &face_nodes);
                        
      mesh_->node_get_coordinates(face_nodes[0], &x0);
      mesh_->node_get_coordinates(face_nodes[1], &x1);

      AmanziGeometry::Point area_cross_product = (xc - x0) ^ (xc - x1);
      double area = norm(area_cross_product) / 2;
            
      B_c[0][c] += (area / mesh_->cell_volume(c)) * (B_n[0][face_nodes[0]] + B_n[0][face_nodes[1]]) / 2;
    }
  }
  // redistribute the result
  S_->GetFieldData(bathymetry_key_)->ScatterMasterToGhosted("cell");
  // initialize h from ht or ht from h
  if (!S_->GetField(ponded_depth_key_, passwd_)->initialized()) {
    const auto& h_c = *S_->GetFieldData(ponded_depth_key_)->ViewComponent("cell");
    auto& ht_c = *S_->GetFieldData(total_depth_key_, passwd_)->ViewComponent("cell");

    for (int c = 0; c < ncells_owned; c++) {
      h_c[0][c] = ht_c[0][c] - B_c[0][c];
    }

    S_->GetField(ponded_depth_key_, passwd_)->set_initialized();
  }
  
  if (!S_->GetField(total_depth_key_, passwd_)->initialized()) {
    const auto& h_c = *S_->GetFieldData(ponded_depth_key_)->ViewComponent("cell");
    auto& ht_c = *S_->GetFieldData(total_depth_key_, passwd_)->ViewComponent("cell");

    for (int c = 0; c < ncells_owned; c++) {
      ht_c[0][c] = h_c[0][c] + B_c[0][c];
    }

    S_->GetField(total_depth_key_, passwd_)->set_initialized();
  }

  InitializeField_(S_.ptr(), passwd_, velocity_key_, 0.0);
  InitializeField_(S_.ptr(), passwd_, discharge_key_, 0.0);

  // secondary fields
  S_->GetFieldEvaluator(hydrostatic_pressure_key_)->HasFieldChanged(S.ptr(), passwd_);
  
  // summary of initialization
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Shallow water PK was initialized." << std::endl;
  }
  //--------------------------------------------------------------
  // Compute P1 element basis values
  //--------------------------------------------------------------

  // calculate volume quadrature points, values etc.
  int nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
  phi_.resize(ncells_owned);
  phi_x_.resize(ncells_owned);
  phi_y_.resize(ncells_owned);
  weights_vol_.resize(ncells_owned);
  
  for (int c = 0; c < ncells_owned; ++c) {
    AmanziMesh::Entity_ID_List cnodes;
    mesh_->cell_get_nodes(c, &cnodes);
    
    // 1. volume quadrature information.
    // 1a. calculate weights
//    int order = 5;
//    int n_points, start_position;
//    n_points = WhetStone::q2d_order[order][0];
//    start_position = WhetStone::q2d_order[order][1];
//    weights_vol_[c].resize(n_points);
//    // find weights of quadrature points
//    for (int qp = 0; qp < n_points; ++qp) {
//      weights_vol_[c][qp] = WhetStone::q2d_weights[qp+start_position] * mesh_->cell_volume(c);
//    }
    
    phi_[c].resize(nnodes_owned);
    phi_x_[c].resize(nnodes_owned);
    phi_y_[c].resize(nnodes_owned);

    // Guassian quadrature points on rectangle (4 points)
    int n_points = 4;
    double quadrature_points[4][2] = {
      {0.211324865405187, 0.211324865405187},
      {0.788675134594813, 0.211324865405187},
      {0.788675134594813, 0.788675134594813},
      {0.211324865405187, 0.788675134594813},
    };
    
    weights_vol_[c].resize(n_points);
    
    for (int j = 0; j < cnodes.size(); ++j) {
      // 1b. construct volume quadrature points
      std::vector<AmanziGeometry::Point> quad_nodes_vol(n_points); // to store physical coordinates of quadrature points
      std::vector<AmanziGeometry::Point> coords(4); // to store physical coordinates of triangular (size 3) /quadrilateral (size 4) cell
      for (int i = 0; i < cnodes.size(); ++i) { // find coordinates of triangle/quadrilateral vertices in mesh
        mesh_->node_get_coordinates(cnodes[i], &coords[i]);
      }
      // find physical coordinates of quadrature points
      for (int qp = 0; qp < n_points; ++qp) {
//        quad_nodes_vol[qp] = (1.0 - WhetStone::q2d_points[qp+start_position][0] - WhetStone::q2d_points[qp+start_position][1] )*coords[0] + (WhetStone::q2d_points[qp+start_position][0])*coords[1] + (WhetStone::q2d_points[qp+start_position][1])*coords[2];
        
        double x = quadrature_points[qp][0], y = quadrature_points[qp][1];
        
        quad_nodes_vol[qp] = (1.0 - x )*(1.0 - y)*coords[0] + x*(1.0 - y)*coords[1] + x*y*coords[2] + (1.0 - x)*y*coords[3];
      }
      
      weights_vol_[c].resize(n_points);
      double Jacobian[2][2];
      for (int  qp = 0; qp < n_points; ++qp) {
        double x = quadrature_points[qp][0], y = quadrature_points[qp][1];
        
        Jacobian[0][0] = -(1.0 - y)*coords[0][0] + (1.0 - y)*coords[1][0] + y*coords[2][0] - y*coords[3][0];
        Jacobian[0][1] = -(1.0 - x )*coords[0][0] - x*coords[1][0] + x*coords[2][0] + (1.0 - x)*coords[3][0];
        Jacobian[1][0] = -(1.0 - y)*coords[0][1] + (1.0 - y)*coords[1][1] + y*coords[2][1] - y*coords[3][1];
        Jacobian[1][1] = -(1.0 - x )*coords[0][1] - x*coords[1][1] + x*coords[2][1] + (1.0 - x)*coords[3][1];
        
        weights_vol_[c][qp] = (0.5)*(0.5)*std::abs((Jacobian[0][0]*Jacobian[1][1] - Jacobian[1][0]*Jacobian[0][1]));
      }
      
      phi_[c][cnodes[j]].resize(n_points);
      phi_x_[c][cnodes[j]].resize(n_points);
      phi_y_[c][cnodes[j]].resize(n_points);

      // 1c. calculate basis values, grad values at volume quadrature points
      for (int qp = 0; qp < quad_nodes_vol.size(); ++qp) {
//        phi_[c][cnodes[j]][qp] = basis_value(cnodes[j], c, quad_nodes_vol[qp]);
//        std::vector<double> grad = basis_grad(cnodes[j], c, quad_nodes_vol[qp]);
//        phi_x_[c][cnodes[j]][qp] = grad[0];
//        phi_y_[c][cnodes[j]][qp] = grad[1];

        phi_[c][cnodes[j]][qp] = basis_value_quad(cnodes[j], c, quad_nodes_vol[qp]);
        std::vector<double> grad = basis_grad_quad(cnodes[j], c, quad_nodes_vol[qp]);
        phi_x_[c][cnodes[j]][qp] = grad[0];
        phi_y_[c][cnodes[j]][qp] = grad[1];
      }
    }
  }

  // calculate face quadrature points, values etc.
  /*
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  phi_face_.resize(nfaces_wghost); // phi_face_[f][j][i] stores the value at the face quadrature point quad_nodes_face[i]. basis is 1 on vertex j and face f.
  weights_face_.resize(nfaces_wghost);
  
  // 2. face quadrature information
  // 2a. calculate face quadrature weights
  int n_points_face = 3;
  for (int f = 0; f < nfaces_wghost; ++f) {
    weights_face_[f].resize(n_points_face); // 3 points order 2*3 - 1
  
    for (int qpf = 0; qpf < n_points_face; ++qpf) {
      weights_face_[f][qpf] = WhetStone::q1d_weights[n_points_face-1][qpf] * mesh_->face_area(f);
    }
  }
  
  // 2b. construct face quadrature points
  for (int f = 0; f < nfaces_wghost; ++f) {
    phi_face_[f].resize(nnodes_owned);

    AmanziMesh::Entity_ID_List fnodes;
    mesh_->face_get_nodes(f, &fnodes);
    std::vector<AmanziGeometry::Point> edge_coords(2); // coordinates of face vertices
    for (int i = 0; i < fnodes.size(); ++i) {
      mesh_->node_get_coordinates(fnodes[i], &edge_coords[i]);
    }

    std::vector<AmanziGeometry::Point> quad_nodes_face(n_points_face);
    for (int qpf = 0; qpf < n_points_face; ++qpf) {
      quad_nodes_face[qpf] = (1.0 - WhetStone::q1d_points[n_points_face-1][qpf] )*edge_coords[0] + (WhetStone::q1d_points[n_points_face-1][qpf])*edge_coords[1];
    }

    // 2c. calculate basis values at face quadrature points
    for (int i = 0; i < fnodes.size(); ++i) {
      phi_face_[f][fnodes[i]].resize(n_points_face);

      int i_2 = (i + 1) % 2;

      for (int qpf = 0; qpf < n_points_face; ++qpf) {
        if (std::abs(edge_coords[i_2][0] - edge_coords[i][0]) < 1.e-12) {
          phi_face_[f][fnodes[i]][qpf] = 1.0 - (quad_nodes_face[qpf][1] - edge_coords[i][1])/(edge_coords[i_2][1] - edge_coords[i][1]);
        }
        else {
          phi_face_[f][fnodes[i]][qpf] = 1.0 - (quad_nodes_face[qpf][0] - edge_coords[i][0])/(edge_coords[i_2][0] - edge_coords[i][0]);
        }
      }
    }
  }
  */
  
} // Initialize()


//--------------------------------------------------------------
// Advance conservative variables: (h, hu, hv)
//--------------------------------------------------------------
bool ShallowWater_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  double dt = t_new - t_old;
  iters_++;

  bool failed = false;
  double eps2 = 1e-12;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
  
  // distribute data to ghost cells
  S_->GetFieldData(total_depth_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(ponded_depth_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(velocity_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(discharge_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(total_depth_key_)->ScatterMasterToGhosted("node");
  S_->GetFieldData(ponded_depth_key_)->ScatterMasterToGhosted("node");
  S_->GetFieldData(velocity_key_)->ScatterMasterToGhosted("node");
  S_->GetFieldData(discharge_key_)->ScatterMasterToGhosted("node");

  // save a copy of primary and conservative fields
  Epetra_MultiVector& B_c = *S_->GetFieldData(bathymetry_key_, passwd_)->ViewComponent("cell", true);
  Epetra_MultiVector& B_n = *S_->GetFieldData(bathymetry_key_, passwd_)->ViewComponent("node", true);
  Epetra_MultiVector& h_c = *S_->GetFieldData(ponded_depth_key_, passwd_)->ViewComponent("cell", true);
  Epetra_MultiVector& h_n = *S_->GetFieldData(ponded_depth_key_, passwd_)->ViewComponent("node", true);
  Epetra_MultiVector& ht_c = *S_->GetFieldData(total_depth_key_, passwd_)->ViewComponent("cell", true);
  Epetra_MultiVector& ht_n = *S_->GetFieldData(total_depth_key_, passwd_)->ViewComponent("node", true);
  Epetra_MultiVector& vel_c = *S_->GetFieldData(velocity_key_, passwd_)->ViewComponent("cell", true);
  Epetra_MultiVector& vel_n = *S_->GetFieldData(velocity_key_, passwd_)->ViewComponent("node", true);

  S_->GetFieldEvaluator(discharge_key_)->HasFieldChanged(S_.ptr(), passwd_);
  Epetra_MultiVector& q_c = *S_->GetFieldData(discharge_key_, discharge_key_)->ViewComponent("cell", true);
  Epetra_MultiVector& q_n = *S_->GetFieldData(discharge_key_, discharge_key_)->ViewComponent("node", true);

  S_->GetFieldEvaluator(hydrostatic_pressure_key_)->HasFieldChanged(S_.ptr(), passwd_);
  
  // create copies of primary fields
  Epetra_MultiVector h_c_tmp(h_c);
  Epetra_MultiVector vel_c_tmp(vel_c);

  // update boundary conditions
  if (bcs_.size() > 0)
      bcs_[0]->Compute(t_old, t_new);
  
  // update source (external) terms
  for (int i = 0; i < srcs_.size(); ++i) {
    srcs_[i]->Compute(t_old, t_new);
  }
  
  // compute source (external) values
  // coupling submodel="rate" returns volumetric flux [m^3/s] integrated over
  // the time step in the last (the second) component of local data vector
  total_source_ = 0.0;
  std::vector<double> ext_S_cell(ncells_owned, 0.0);
  for (int  i = 0; i < srcs_.size(); ++i) {
    for (auto it = srcs_[i]->begin(); it != srcs_[i]->end(); ++it) {
      int c = it->first;
      ext_S_cell[c] = it->second[0];  // data unit is [m/s]
      total_source_ += it->second[0] * mesh_->cell_volume(c) * dt; // data unit is [m^3]
    }
  }
  
  // Shallow water equations have the form
  // U_t + F_x(U) + G_y(U) = S(U)

  std::vector<std::vector<double>> U(3);
  for (int m = 0; m < 3; ++m) {
    U[m].resize(nnodes_owned);
  }
  for (int i = 0; i < nnodes_owned; ++i) {
    U[0][i] = ht_n[0][i] - B_n[0][i];
    U[1][i] = q_n[0][i];
    U[2][i] = q_n[1][i];
  }

  std::vector<std::vector<double>> U_pr(3);
  U_pr = U;
  
  std::vector<std::vector<double>> phi(3);                                      // phi[m][j]
  std::vector<double> phi_tmp(3), Phi_total(3), beta(3, 0.0), sum_max(3, 0.0);  // phi[m][j] = phi_tmp[m], sum_max[m] eq (7) denominator
  std::vector<double> phi_beta_cell(3, 0.0);                                    // sum used in eq (4)
  double dual_cell_vol = 0.0;
  
  // 1. predictor step
  for (int i = 0; i < nnodes_owned; ++i){
    AmanziGeometry::Point node_coordinates;
    mesh_->node_get_coordinates(i, &node_coordinates);
  
    std::fill(phi_beta_cell.begin(), phi_beta_cell.end(), 0.0);
    dual_cell_vol = 0.0;
    
    
    // loop over cells joined to the vertex i
    AmanziMesh::Entity_ID_List ncells, cnodes;
    mesh_->node_get_cells(i, Amanzi::AmanziMesh::Parallel_type::ALL, &ncells);
    
    for (int K = 0; K < ncells.size(); ++K) {
      mesh_->cell_get_nodes(ncells[K], &cnodes);
      
      for (int m = 0; m < 3; ++m) {
        phi[m].resize(cnodes.size());
      }
      
      std::fill(Phi_total.begin(), Phi_total.end(), 0.0);
      for (int j = 0; j < cnodes.size(); ++j) {
        phi_tmp = ResidualsLF(ncells[K], cnodes[j], U); // eq (10)
        
        for (int m = 0; m < 3; ++m) {
          Phi_total[m] += phi_tmp[m];     // eq(8)
          phi[m][j] = phi_tmp[m];         // eq (10)
        }
      }
      
      std::fill(beta.begin(), beta.end(), 0.0);
      std::fill(sum_max.begin(), sum_max.end(), 0.0);
      for (int m = 0; m < 3; ++m) {
        
        if (std::abs(Phi_total[m]) > 0.0) {
          
          for (int j = 0; j < cnodes.size(); ++j) {
            
            if (cnodes[j] == i) {
              beta[m] = std::max(0.0, phi[m][j]/Phi_total[m]);
            }
            sum_max[m] +=  std::max(0.0, phi[m][j]/Phi_total[m]);
          }
        }
        
        if (std::abs(sum_max[m]) > 0.0) {
          beta[m] = beta[m] / sum_max[m];
        }
        phi_beta_cell[m] += beta[m] * Phi_total[m]; // eq(6)
        
      }
      dual_cell_vol += (1.0/cnodes.size())*mesh_->cell_volume(ncells[K]);
      
      //        std::vector<double> phi_beta_cell_tmp = ResidualsLF(ncells[K], i, U);
      //        for (int m = 0; m < 3; ++m) {
      //          phi_beta_cell[m] += phi_beta_cell_tmp[m];
      //        }
      
    } // K (cell) loop
    
    
    for (int m = 0; m < 3; ++m) {
      U_pr[m][i] = -(dt/dual_cell_vol)*phi_beta_cell[m] + U[m][i];
      // boundary condition
      if (std::abs(node_coordinates[0] - 0.0) < 1.e-12 || std::abs(node_coordinates[0] - 1.0) < 1.e-12 || std::abs(node_coordinates[1] - 0.0) < 1.e-12 || std::abs(node_coordinates[1] - 1.0) < 1.e-12) {
        //U_pr[0][i] = 0.5;
//        U_pr[1][i] = 0.0;
//        U_pr[2][i] = 0.0;
      }
    }
    
  } // i (total DOF) loop
  // predictor step ends
  
  for (int i = 0; i < nnodes_owned; ++i) {
    h_n[0][i] = U_pr[0][i];
    q_n[0][i] = U_pr[1][i];
    q_n[1][i] = U_pr[2][i];
    ht_n[0][i] = h_n[0][i] + B_n[0][i];
    double h = h_n[0][i];
    double factor = (2.0*h)/(h*h + std::max(h*h, 1.e-14));
    vel_n[0][i] = factor * q_n[0][i];
    vel_n[1][i] = factor * q_n[1][i];
  }
  
/*
  // 2. corrector step begins
  std::vector<std::vector<double>> U_new = U_pr;
  
  for (int i = 0; i < nnodes_owned; ++i){
    
    AmanziGeometry::Point node_coordinates;
    mesh_->node_get_coordinates(i, &node_coordinates);
    
    // compute on interior nodes only

    {
      
      std::fill(phi_beta_cell.begin(), phi_beta_cell.end(), 0.0);
      dual_cell_vol = 0.0;
      
      // loop over cells joined to the vertex i
      AmanziMesh::Entity_ID_List ncells, cnodes;
      mesh_->node_get_cells(i, Amanzi::AmanziMesh::Parallel_type::ALL, &ncells);
      
      for (int K = 0; K < ncells.size(); ++K) {
        mesh_->cell_get_nodes(ncells[K], &cnodes);
        
        for (int m = 0; m < 3; ++m) {
          phi[m].resize(cnodes.size());
        }
        
        std::fill(Phi_total.begin(), Phi_total.end(), 0.0);
        for (int j = 0; j < cnodes.size(); ++j) {
          phi_tmp = ResidualsTimeSpace(ncells[K], cnodes[j], U, U_pr, dt); // eq (9)
          
          for (int m = 0; m < 3; ++m) {
            Phi_total[m] += phi_tmp[m];     // eq(8)
            phi[m][j] = phi_tmp[m];         // eq (9)
          }
        }
        
        std::fill(beta.begin(), beta.end(), 0.0);
        std::fill(sum_max.begin(), sum_max.end(), 0.0);
        for (int m = 0; m < 3; ++m) {
          
          if (std::abs(Phi_total[m]) > 0.0) {
            
            for (int j = 0; j < cnodes.size(); ++j) {
              
              if (cnodes[j] == i) {
                beta[m] = std::max(0.0, phi[m][j]/Phi_total[m]);
              }
              sum_max[m] +=  std::max(0.0, phi[m][j]/Phi_total[m]);
            }
          }
          
          if (std::abs(sum_max[m]) > 0.0) {
            beta[m] = beta[m] / sum_max[m];
          }
          
          phi_beta_cell[m] += beta[m] * Phi_total[m]; // eq(6)
        }
        dual_cell_vol += (1.0/3)*mesh_->cell_volume(ncells[K]);
      } // K (cell) loop
      
      for (int m = 0; m < 3; ++m) {
        U_new[m][i] = -(dt/dual_cell_vol)*phi_beta_cell[m] + U_pr[m][i];
        
        if (std::abs(node_coordinates[0] - 0.0) < 1.e-12 || std::abs(node_coordinates[0] - 1.0) < 1.e-12 || std::abs(node_coordinates[1] - 0.0) < 1.e-12 || std::abs(node_coordinates[1] - 1.0) < 1.e-12) {
          U_new[1][i] = 0.0;
          U_new[2][i] = 0.0;
        }
        
      }
    }
  }// i (total DOF) loop
  // corrector step ends
*/

//  for (int i = 0; i < nnodes_owned; ++i) {
//    h_n[0][i] = U_new[0][i];
//    q_n[0][i] = U_new[1][i];
//    q_n[1][i] = U_new[2][i];
//    ht_n[0][i] = ht_n[0][i] + B_n[0][i];
//    double h = h_n[0][i];
//    double factor = (2.0*h)/(h*h + std::max(h*h, 1.e-14));
//    vel_n[0][i] = factor * q_n[0][i];
//    vel_n[1][i] = factor * q_n[1][i];
//  }
  
  // compute cell averaged quantities
  for (int c = 0; c < ncells_owned; ++c) {
    const Amanzi::AmanziGeometry::Point &xc = mesh_->cell_centroid(c);
        
    Amanzi::AmanziMesh::Entity_ID_List cfaces;
    mesh_->cell_get_faces(c, &cfaces);
    int nfaces_cell = cfaces.size();
        
    ht_c[0][c] = 0.0;
    h_c[0][c] = 0.0;
    vel_c[0][c] = 0.0;
    vel_c[1][c] = 0.0;
    q_c[0][c] = 0.0;
    q_c[1][c] = 0.0;
    
    for (int f = 0; f < nfaces_cell; ++f) {
    
      Amanzi::AmanziGeometry::Point x0, x1;
      int edge = cfaces[f];
            
      AmanziMesh::Entity_ID_List face_nodes;
      mesh_->face_get_nodes(edge, &face_nodes);
                        
      mesh_->node_get_coordinates(face_nodes[0], &x0);
      mesh_->node_get_coordinates(face_nodes[1], &x1);

      AmanziGeometry::Point area_cross_product = (xc - x0) ^ (xc - x1);
      double area = norm(area_cross_product) / 2;

      h_c[0][c] += (area / mesh_->cell_volume(c)) * (h_n[0][face_nodes[0]] + h_n[0][face_nodes[1]]) / 2;
      q_c[0][c] += (area / mesh_->cell_volume(c)) * (q_n[0][face_nodes[0]] + q_n[0][face_nodes[1]]) / 2;
      q_c[1][c] += (area / mesh_->cell_volume(c)) * (q_n[1][face_nodes[0]] + q_n[1][face_nodes[1]]) / 2;
      vel_c[0][c] += (area / mesh_->cell_volume(c)) * (vel_n[0][face_nodes[0]] + vel_n[0][face_nodes[1]]) / 2;
      vel_c[1][c] += (area / mesh_->cell_volume(c)) * (vel_n[1][face_nodes[0]] + vel_n[1][face_nodes[1]]) / 2;
    }
    ht_c[0][c] = h_c[0][c] + B_c[0][c];
  }
  return failed;
}
  

//--------------------------------------------------------------
// Advance conservative variables: (h, hu, hv)
//--------------------------------------------------------------
void ShallowWater_PK::CommitStep(
    double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(S->GetFieldEvaluator(velocity_key_))->SetFieldAsChanged(S.ptr());
  Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(S->GetFieldEvaluator(ponded_depth_key_))->SetFieldAsChanged(S.ptr());
}


//--------------------------------------------------------------
// Lax-Friedrichs residual
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::ResidualsLF(int K, int j, std::vector<std::vector<double> > U) // input argument must contain coefficients of the basis expansion. cell c, node i, U[m][] contains the DOFs for mth component
{
  // calculate -\int_K F \cdot \nabla \phi_j + \int_{\partial K} (F \cdot n) \phi_j - \int_{K} S \phi_j + \alpha (U_m - Ubar) [eq(10)]
  std::vector<double> integral(3, 0.0); // value to retrun
  
  AmanziMesh::Entity_ID_List cnodes; //DOFs of the cell (nodes)
  mesh_->cell_get_nodes(K, &cnodes);

  Epetra_MultiVector& ht_n = *S_->GetFieldData(total_depth_key_, passwd_)->ViewComponent("node", true);
  
  // 1. calculate volume integral
  std::vector<double> Uqp(3), Sqp(3), Uxqp(3), Uyqp(3);
  std::vector<std::vector<double>> flux(2);
  std::vector<double> phi_j_grad(2);
  double h, qx, qy; // values conservative variables
  double h_x, h_y, qx_x, qx_y, qy_x, qy_y; // gradients
  double phi_j;
  for (int qp = 0; qp < weights_vol_[K].size(); ++qp ) {
    Uqp = EvalSol_vol(U, qp, K);
    Sqp = EvalPhySource_vol(U, qp, K);

    phi_j = phi_[K][j][qp];
    
    Uxqp = EvalSol_x_vol(U, qp, K);
    Uyqp = EvalSol_y_vol(U, qp, K);
    
    h = Uqp[0], qx = Uqp[1], qy = Uqp[2];
    
    h_x = Uxqp[0], h_y = Uyqp[0];
    qx_x = Uxqp[1], qx_y = Uyqp[1];
    qy_x = Uxqp[2], qy_y = Uyqp[2];
    
    integral[0] += (qx_x + qy_y)*phi_j*weights_vol_[K][qp];
    
    integral[1] += (2*qx*qx_x/h - qx*qx*h_x/(h*h) + g_*h*h_x + qy*qx_y/h + qx*qy_y/h - qx*qy*h_y/(h*h))*phi_j*weights_vol_[K][qp];

    integral[2] += (qy*qx_x/h + qx*qy_x/h - qx*qy*h_x/(h*h) + 2*qy*qy_y/h - qy*qy*h_y/(h*h) + g_*h*h_y)*phi_j*weights_vol_[K][qp];
    
    for (int m = 0; m < 3; ++m) {
      integral[m] += Sqp[m] * phi_j * weights_vol_[K][qp]; // bathymetry term
    }
    
//    double ht_x = 0.0, ht_y = 0.0;
//    for (int i = 0; i < cnodes.size(); ++i) {
//      ht_x += ht_n[0][cnodes[i]] * phi_x_[K][cnodes[i]][qp];
//      ht_y += ht_n[0][cnodes[i]] * phi_y_[K][cnodes[i]][qp];
//    }
//
//    integral[1] += (2*qx*qx_x/h - qx*qx*h_x/(h*h) + g_*h*ht_x + qy*qx_y/h + qx*qy_y/h - qx*qy*h_y/(h*h))*phi_j*weights_vol_[K][qp];
//
//    integral[2] += (qy*qx_x/h + qx*qy_x/h - qx*qy*h_x/(h*h) + 2*qy*qy_y/h - qy*qy*h_y/(h*h) + g_*h*ht_y)*phi_j*weights_vol_[K][qp];
  }

  // 3. compute viscoisty term
  // 3a. compute U_bar i.e. average
  std::vector<double> U_avg(3, 0.0);
  for (int m = 0; m < 3; ++m) {
    for (int i = 0; i < cnodes.size(); ++i) {
      U_avg[m] += U[m][cnodes[i]];
    }
    U_avg[m] = U_avg[m]/cnodes.size();
  }
  // 3b. compute artificial viscosity \alpha
//  double h, qx, qy, vx, vy,
  double vx, vy, S, Smax = 0.0;
  for (int i = 0; i < cnodes.size(); ++i) {
    h = U[0][cnodes[i]];
    qx = U[1][cnodes[i]];
    qy = U[2][cnodes[i]];
    vx = 2.0 * h * qx / (h*h + std::max(h*h, 1.e-14));
    vy = 2.0 * h * qy / (h*h + std::max(h*h, 1.e-14));
    
    S = std::max(std::abs(vx) + std::sqrt(g_*h), std::abs(vy) + std::sqrt(g_*h));
    Smax = std::max(Smax, S);
  }
  
  AmanziMesh::Entity_ID_List cfaces;
  mesh_->cell_get_faces(K, &cfaces);
  double lmax = 0.0;
  for (int f = 0; f < cfaces.size(); ++f) {
    double farea = mesh_->face_area(cfaces[f]);
    lmax = std::max(lmax, farea);
  }
  
  double alpha = 1.0* lmax * Smax;
  
  for (int m = 0; m < 3; ++m) {
    integral[m] += alpha * (U[m][j] - U_avg[m]);
  }
  
  return integral;
}


//--------------------------------------------------------------
// Time-space residuals for time-stepping scheme
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::ResidualsTimeSpace(int K, int j, std::vector<std::vector<double> > U, std::vector<std::vector<double> > U_pr, double dt)
{
  // compute time step integral \int_{K} ((u - u_star)/\Delta t)  * \phi_i dx
  std::vector<double> integral(3, 0.0);

  AmanziMesh::Entity_ID_List cnodes; //DOFs of the cell (nodes)
  mesh_->cell_get_nodes(K, &cnodes);

  std::vector<double> Uqp(3);
  std::vector<double> U_prqp(3);
  double phi_j;
  for (int qp = 0; qp < weights_vol_[K].size(); ++qp ) {
    Uqp = EvalSol_vol(U, qp, K);
    U_prqp = EvalSol_vol(U_pr, qp, K);

    phi_j = phi_[K][j][qp];

    for (int m = 0; m < 3; ++m) {
    integral[m] += (1.0/dt) * (U_prqp[m] - Uqp[m]) * phi_j * weights_vol_[K][qp];
    }
  }
  
  std::vector<double> phi_K_j, phi_K_j_star;

  phi_K_j = ResidualsLF(K, j, U);
  phi_K_j_star = ResidualsLF(K, j, U_pr);
  
  for (int m = 0; m < 3; ++m) {
    integral[m] += 0.5 * (phi_K_j[m] + phi_K_j_star[m]);
  }
  
  return integral;
}


//--------------------------------------------------------------
// Evaluate solution at quadrature point
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::EvalSol_vol(std::vector<std::vector<double>> U, int qp, int c) // point x_qp lies in cell c
{
  std::vector<double> eval_sol(3, 0.0);
  
  // find nodes of the cell c and loop over them to find \sum_{i \in K} U_i(t) \phi_i(x,y)
  AmanziMesh::Entity_ID_List cnodes;
  mesh_->cell_get_nodes(c, &cnodes);
  for (int m = 0; m < 3; ++m) {
    for (int i = 0; i < cnodes.size(); ++i) {
      eval_sol[m] += U[m][cnodes[i]] * phi_[c][cnodes[i]][qp];
    }
  }
  return eval_sol;
}


//--------------------------------------------------------------
// Evaluate solution gradient (x) at quadrature point
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::EvalSol_x_vol(std::vector<std::vector<double>> U, int qp, int c) // point x_qp lies in cell c
{
  std::vector<double> eval_sol(3, 0.0);
  
  // find nodes of the cell c and loop over them to find \sum_{i \in K} U_i(t) \phi_i(x,y)
  AmanziMesh::Entity_ID_List cnodes;
  mesh_->cell_get_nodes(c, &cnodes);
  for (int m = 0; m < 3; ++m) {
    for (int i = 0; i < cnodes.size(); ++i) {
      eval_sol[m] += U[m][cnodes[i]] * phi_x_[c][cnodes[i]][qp];
    }
  }
  return eval_sol;
}


//--------------------------------------------------------------
// Evaluate solution gradient (y) at quadrature point
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::EvalSol_y_vol(std::vector<std::vector<double>> U, int qp, int c) // point x_qp lies in cell c
{
  std::vector<double> eval_sol(3, 0.0);
  
  // find nodes of the cell c and loop over them to find \sum_{i \in K} U_i(t) \phi_i(x,y)
  AmanziMesh::Entity_ID_List cnodes;
  mesh_->cell_get_nodes(c, &cnodes);
  for (int m = 0; m < 3; ++m) {
    for (int i = 0; i < cnodes.size(); ++i) {
      eval_sol[m] += U[m][cnodes[i]] * phi_y_[c][cnodes[i]][qp];
    }
  }
  return eval_sol;
}

//--------------------------------------------------------------
// Evaluate solution at face quadrature point
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::EvalSol_face(std::vector<std::vector<double>> U, int qpf, int f) // point x_qp lies on face f
{
  std::vector<double> eval_sol(3, 0.0);
  
  // find nodes of the face and loop over them to calcylate \sum_{i \in f \subset \partial K} U_i(t) \phi_i(x,y)
  AmanziMesh::Entity_ID_List fnodes;
  mesh_->face_get_nodes(f, &fnodes);
  for (int m = 0; m < 3; ++m) {
    for (int i = 0; i < fnodes.size(); ++i) {
      eval_sol[m] += U[m][fnodes[i]] * phi_face_[f][fnodes[i]][qpf];
    }
  }
  return eval_sol;
}


//--------------------------------------------------------------
// Evaluate physical source at quadrature point
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::EvalPhySource_vol(std::vector<std::vector<double>> U, int qp, int c) // point x_qp lies in cell c
{
  std::vector<double> eval_sol(3, 0.0);
  
  Epetra_MultiVector& B_n = *S_->GetFieldData(bathymetry_key_, passwd_)->ViewComponent("node", true);

  AmanziMesh::Entity_ID_List cnodes;
  mesh_->cell_get_nodes(c, &cnodes);
  std::vector<double> phi_i_grad(2);
  double h = 0.0;
  for (int i = 0; i < cnodes.size(); ++i) {
    h += U[0][cnodes[i]] * phi_[c][cnodes[i]][qp];
  }
  
  double B_x = 0.0, B_y = 0.0;
  for (int i = 0; i < cnodes.size(); ++i) {
    B_x += B_n[0][cnodes[i]] * phi_x_[c][cnodes[i]][qp];
    B_y += B_n[0][cnodes[i]] * phi_y_[c][cnodes[i]][qp];
  }
  
  eval_sol[0] = 0.0;
  eval_sol[1] = g_ * h * B_x;
  eval_sol[2] = g_ * h * B_y;
  
  return eval_sol;
}

//--------------------------------------------------------------
// Basis value (change to barycentric implementation)
//--------------------------------------------------------------
double ShallowWater_PK::basis_value(int i, int c, AmanziGeometry::Point x) // DOF (vertex), cell, evaluation point
{
  std::vector<AmanziGeometry::Point> x_vertex(3), x_vertex_tria(3);
  
  AmanziMesh::Entity_ID_List cnodes;
  mesh_->cell_get_nodes(c, &cnodes);
  
  for (int j = 0; j < cnodes.size(); ++j) { // construct vertices of the plane over triangle
    mesh_->node_get_coordinates(cnodes[j], &x_vertex_tria[j]);
    AmanziGeometry::Point x_tmp = x_vertex_tria[j];
    double x0 = x_tmp[0], x1 = x_tmp[1], x2 = 0.0;

    if (cnodes[j] == i) {
      x2 = 1.0;
    }
    AmanziGeometry::Point x_tmp_2(x0, x1, x2);
    x_vertex[j] = x_tmp_2;
  }
  AmanziGeometry::Point x0, x1, x2;
  
  x0 = x_vertex[0]; // vertices of plane over triangle element
  x1 = x_vertex[1];
  x2 = x_vertex[2];
  
  AmanziGeometry::Point edge_0 = x0 - x1, edge_1 = x2 - x1;
  AmanziGeometry::Point n = edge_0^edge_1; // normal to plane
  
  return -( (x[0] - x1[0])*n[0] + (x[1] - x1[1])*n[1] )/ (n[2]) + x1[2]; // (x1 - x) \cdot n = 0
}


//--------------------------------------------------------------
// Basis gradient value (change to barycentric implementation)
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::basis_grad(int i, int c, AmanziGeometry::Point x)
{
  std::vector<double> grad(2);
  std::vector<AmanziGeometry::Point> x_vertex(3), x_vertex_tria(3);
  
  AmanziMesh::Entity_ID_List cnodes;
  mesh_->cell_get_nodes(c, &cnodes);
  
  for (int j = 0; j < cnodes.size(); ++j) { // construct vertices of the plane over triangle
    mesh_->node_get_coordinates(cnodes[j], &x_vertex_tria[j]);
    AmanziGeometry::Point x_tmp = x_vertex_tria[j];
    double x0 = x_tmp[0], x1 = x_tmp[1], x2 = 0.0;

    if (cnodes[j] == i) {
      x2 = 1.0;
    }
    AmanziGeometry::Point x_tmp_2(x0, x1, x2);
    x_vertex[j] = x_tmp_2;
  }
  
  AmanziGeometry::Point x0, x1, x2;
  
  x0 = x_vertex[0]; // vertices of plane on triangle element
  x1 = x_vertex[1];
  x2 = x_vertex[2];
  
  AmanziGeometry::Point edge_0 = x0 - x1, edge_1 = x2 - x1;
  AmanziGeometry::Point n = edge_0^edge_1;
  
  grad[0] = -n[0]/n[2];
  grad[1] = -n[1]/n[2];
  
  return grad;
}


//--------------------------------------------------------------
// Basis value for quadrilaterals
//--------------------------------------------------------------
double ShallowWater_PK::basis_value_quad(int i, int c, AmanziGeometry::Point x) // DOF (vertex), cell, evaluation point
{
  std::vector<AmanziGeometry::Point> x_vertex(4), x_vertex_quad(4);
  WhetStone::DenseVector F(4);
  
  AmanziMesh::Entity_ID_List cnodes;
  mesh_->cell_get_nodes(c, &cnodes);
  
  for (int j = 0; j < cnodes.size(); ++j) {
    mesh_->node_get_coordinates(cnodes[j], &x_vertex_quad[j]);
    AmanziGeometry::Point x_tmp = x_vertex_quad[j];
    double x0 = x_tmp[0], x1 = x_tmp[1], x2 = 0.0; // (x0, x1, x2) = (x0, x1, \phi(x0,x1))

    if (cnodes[j] == i) {
      x2 = 1.0;
    }
    AmanziGeometry::Point x_tmp_2(x0, x1, x2);
    x_vertex[j] = x_tmp_2;
  }
  
  WhetStone::DenseMatrix M(4,4);
  M.PutScalar(0.0);
  
  for (int i = 0; i < 4; ++i) {
    AmanziGeometry::Point x_tmp = x_vertex[i];
    M(i,0) = x_tmp[0] * x_tmp[1];
    M(i,1) = x_tmp[0];
    M(i,2) = x_tmp[1];
    M(i,3) = 1.0;
    F(i) = x_tmp[2];
  }
  
  M.Inverse();
  
  WhetStone::DenseVector Coeffs(4);
  M.Multiply(F, Coeffs, false);

  return Coeffs(0)*x[0]*x[1] + Coeffs(1)*x[0] + Coeffs(2)*x[1] + Coeffs(3);
}


//--------------------------------------------------------------
// Basis gradient value for quadrilaterals
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::basis_grad_quad(int i, int c, AmanziGeometry::Point x) // DOF (vertex), cell, evaluation point
{
  std::vector<double> grad(2);

  std::vector<AmanziGeometry::Point> x_vertex(4), x_vertex_quad(4);
  WhetStone::DenseVector F(4);
  
  AmanziMesh::Entity_ID_List cnodes;
  mesh_->cell_get_nodes(c, &cnodes);
  
  for (int j = 0; j < cnodes.size(); ++j) { // construct vertices of the plane over triangle
    mesh_->node_get_coordinates(cnodes[j], &x_vertex_quad[j]);
    AmanziGeometry::Point x_tmp = x_vertex_quad[j];
    double x0 = x_tmp[0], x1 = x_tmp[1], x2 = 0.0; // (x0, x1, x2) = (x0, x1, \phi(x0,x1))

    if (cnodes[j] == i) {
      x2 = 1.0;
    }
    AmanziGeometry::Point x_tmp_2(x0, x1, x2);
    x_vertex[j] = x_tmp_2;
  }
  
  WhetStone::DenseMatrix M(4,4);
  M.PutScalar(0.0);
  
  for (int i = 0; i < 4; ++i) {
    AmanziGeometry::Point x_tmp = x_vertex[i];
    M(i,0) = x_tmp[0] * x_tmp[1];
    M(i,1) = x_tmp[0];
    M(i,2) = x_tmp[1];
    M(i,3) = 1.0;
    F(i) = x_tmp[2];
  }
  
  M.Inverse();
  
  WhetStone::DenseVector Coeffs(4);
  M.Multiply(F, Coeffs, false);
  
  grad[0] = Coeffs(0)*x[1] + Coeffs(1);
  grad[1] = Coeffs(0)*x[0] + Coeffs(2);

  return grad;
}


//--------------------------------------------------------------
// Physical source term S(U) = (0, -ghB_x, -ghB_y)
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::PhysicalSource(const std::vector<double>& U)
{
  double h, u, v, qx, qy;
  double eps2 = 1e-12;

  // SW conservative variables: (h, hu, hv)
  h  = U[0];
  qx = U[1];
  qy = U[2];
  u  = 2.0 * h * qx / (h*h + std::fmax(h*h, eps2));
  v  = 2.0 * h * qy / (h*h + std::fmax(h*h, eps2));

  double dBathx = 0.0, dBathy = 0.0;

  std::vector<double> S(3);
  S[0] = 0.0;
  S[1] = -g_ * h * dBathx;
  S[2] = -g_ * h * dBathy;

  return S;
}


//--------------------------------------------------------------
// Calculation of time step limited by the CFL condition
//--------------------------------------------------------------
double ShallowWater_PK::get_dt()
{
  double d, vn, dt = 1.e10;

  Epetra_MultiVector& h_c = *S_->GetFieldData(ponded_depth_key_, passwd_)->ViewComponent("cell", true);
  Epetra_MultiVector& vel_c = *S_->GetFieldData(velocity_key_, passwd_)->ViewComponent("cell", true);

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  AmanziMesh::Entity_ID_List cfaces;
  
  for (int c = 0; c < ncells_owned; c++) {
    const Amanzi::AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

    mesh_->cell_get_faces(c, &cfaces);

    for (int n = 0; n < cfaces.size(); ++n) {
      int f = cfaces[n];
      double farea = mesh_->face_area(f);
      const auto& xf = mesh_->face_centroid(f);
      const auto& normal = mesh_->face_normal(f);
  
      double h = h_c[0][c];
      double vx = vel_c[0][c];
      double vy = vel_c[1][c];

      // computing local (cell, face) time step using Kurganov's estimate d / (2a)
      vn = (vx * normal[0] + vy * normal[1]) / farea;
      d = norm(xc - xf);
      dt = std::min(d / (2 * (std::abs(vn) + std::sqrt(g_ * h))), dt);
    }
  }
  
  double dt_min;
  mesh_->get_comm()->MinAll(&dt, &dt_min, 1);

  if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "stable dt=" << dt_min << ", cfl=" << cfl_ << std::endl;
  }

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH && iters_ == max_iters_) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "switching from reduced to regular cfl=" << cfl_ << std::endl;
  }

  if (iters_ < max_iters_)
    return 0.1 * cfl_ * dt_min;
  else
    return cfl_ * dt_min;
}


//--------------------------------------------------------------
// Bathymetry (Linear construction for a rectangular cell)
//--------------------------------------------------------------
/*
double ShallowWater_PK::BathymetryRectangularCellValue(
    int c, const AmanziGeometry::Point& xp, const Epetra_MultiVector& B_n)
{
  double x = xp[0], y = xp[1];
    
  AmanziMesh::Entity_ID_List nodes, faces;
    
  mesh_->cell_get_faces(c, &faces);
  mesh_->cell_get_nodes(c, &nodes);
    
  double dx = mesh_->face_area(faces[0]), dy = mesh_->face_area(faces[1]);
    
  Amanzi::AmanziGeometry::Point xl;
  mesh_->node_get_coordinates(nodes[0], &xl); // Lower left corner of the cell (for rectangular cell)
 
  // Values of B at the corners of the cell
  double B1 = B_n[0][nodes[0]];
  double B2 = B_n[0][nodes[1]];
  double B3 = B_n[0][nodes[3]];
  double B4 = B_n[0][nodes[2]];
    
  double xln = xl[0], yln = xl[1];
  double B_rec = B1 + (B2 - B1)*(xp[0] - xln)/dx + (B3 - B1)*(xp[1] - yln)/dy + (B4 - B2 - B3 + B1)*(xp[0] - xln)*(xp[1] - yln)/(dx*dy) ;
  return B_rec;
}
*/


//--------------------------------------------------------------
// Bathymetry (Evaluate value at edge midpoint for a polygonal cell)
//--------------------------------------------------------------
double ShallowWater_PK::BathymetryEdgeValue(int e, const Epetra_MultiVector& B_n)
{
  AmanziMesh::Entity_ID_List nodes;
  mesh_->face_get_nodes(e, &nodes);

  return (B_n[0][nodes[0]] + B_n[0][nodes[1]]) / 2;
}


//--------------------------------------------------------------
// Physical flux in x-direction
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::PhysFlux_x(std::vector<double> U)
{
  std::vector<double> F;

  F.resize(3);

  double h, u, v, qx, qy, g = g_;
  double eps = 1.e-6;

  // SW conservative variables: (h, hu, hv)

  h  = U[0];
  qx = U[1];
  qy = U[2];
  u  = 2.*h*qx/(h*h + std::max(h*h,eps*eps));
  v  = 2.*h*qy/(h*h + std::max(h*h,eps*eps));

  // Form vector of x-fluxes F(U) = (hu, hu^2 + 1/2 gh^2, huv)

  F[0] = h*u;
  F[1] = h*u*u+0.5*g*h*h;
  F[2] = h*u*v;

  return F;
}


//--------------------------------------------------------------
// Physical flux in y-direction
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::PhysFlux_y(std::vector<double> U)
{
  std::vector<double> G;

  G.resize(3);

  double h, u, v, qx, qy, g = g_;
  double eps = 1.e-6;

  // SW conservative variables: (h, hu, hv)

  h  = U[0];
  qx = U[1];
  qy = U[2];
  u  = 2.*h*qx/(h*h + std::max(h*h,eps*eps));
  v  = 2.*h*qy/(h*h + std::max(h*h,eps*eps));

  // Form vector of y-fluxes G(U) = (hv, huv, hv^2 + 1/2 gh^2)

  G[0] = h*v;
  G[1] = h*u*v;
  G[2] = h*v*v+0.5*g*h*h;

  return G;
}


//--------------------------------------------------------------
// Barycentric coordinates (work in progress)
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::get_barycentric(std::vector<AmanziGeometry::Point> vertices, AmanziGeometry::Point x) // triangle clockwise vertices 0, 1, 2
{
  std::vector<double> bary_coords(3, 0.0); // return value
  AmanziGeometry::Point edge_0, edge_1;
  edge_0 = vertices[2] - vertices[0];
  edge_1 = vertices[1] - vertices[0];
  double area = 0.5 * norm(edge_0^edge_1);
  
  for (int m = 0; m < 3; ++m) {
    edge_0 = x - vertices[m];
    int n = (m +1) % 3;
    edge_0 = vertices[n] - vertices[m];
    edge_1 = x - vertices[m];
    double sub_area = 0.5 * norm(edge_0 ^ edge_1);
    bary_coords[m] = sub_area/area;
  }
  
}


//--------------------------------------------------------------
// Error diagnostics
//--------------------------------------------------------------
bool ShallowWater_PK::ErrorDiagnostics_(int c, double h, double B, double ht)
{
  if (h < 0.0) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "negative height in cell " << c 
               << ", total=" << ht 
               << ", bathymetry=" << B
               << ", height=" << h << std::endl;
    return true;
  }
  return false;
}

}  // namespace ShallowWater
}  // namespace Amanzi

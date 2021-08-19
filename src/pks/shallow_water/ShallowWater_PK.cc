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
    S->RequireField(ponded_depth_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    AddDefaultPrimaryEvaluator_(ponded_depth_key_);
  }

  // total_depth_key_
  if (!S->HasField(total_depth_key_)) {
    S->RequireField(total_depth_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // velocity
  if (!S->HasField(velocity_key_)) {
    S->RequireField(velocity_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 2);
    AddDefaultPrimaryEvaluator_(velocity_key_);
  }

  // discharge
  if (!S->HasField(discharge_key_)) {
    S->RequireField(discharge_key_, discharge_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 2);

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

  // reconstruction
  Teuchos::ParameterList plist = sw_list_->sublist("reconstruction");
  
  total_depth_grad_ = Teuchos::rcp(new Operators::ReconstructionCell(mesh_));
  total_depth_grad_->Init(plist);

  discharge_x_grad_ = Teuchos::rcp(new Operators::ReconstructionCell(mesh_));
  discharge_x_grad_->Init(plist);

  discharge_y_grad_ = Teuchos::rcp(new Operators::ReconstructionCell(mesh_));
  discharge_y_grad_->Init(plist);

  limiter_ = Teuchos::rcp(new Operators::LimiterCell(mesh_));
  limiter_->Init(plist);

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
}


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

  // save a copy of primary and conservative fields
  Epetra_MultiVector& B_c = *S_->GetFieldData(bathymetry_key_, passwd_)->ViewComponent("cell", true);
  Epetra_MultiVector& B_n = *S_->GetFieldData(bathymetry_key_, passwd_)->ViewComponent("node", true);
  Epetra_MultiVector& h_c = *S_->GetFieldData(ponded_depth_key_, passwd_)->ViewComponent("cell", true);
  Epetra_MultiVector& ht_c = *S_->GetFieldData(total_depth_key_, passwd_)->ViewComponent("cell", true);
  Epetra_MultiVector& vel_c = *S_->GetFieldData(velocity_key_, passwd_)->ViewComponent("cell", true);

  S_->GetFieldEvaluator(discharge_key_)->HasFieldChanged(S_.ptr(), passwd_);
  Epetra_MultiVector& q_c = *S_->GetFieldData(discharge_key_, discharge_key_)->ViewComponent("cell", true);

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
    U[0][i] = 0.5; // place holder
    U[1][i] = 0.0; // place holder
    U[2][i] = 0.0; // place holder
  }

  std::vector<std::vector<double>> U_pr(3);
  U_pr = U;
  
  // 1. predictor step
  for (int i = 0; i < nnodes_owned; ++i){
    
    std::vector<double> phi_beta_cell(3, 0.0); // sum used in eq (4)
    double dual_cell_vol = 0.0;
    
    AmanziMesh::Entity_ID_List ncells;
    mesh_->node_get_cells(i, AmanziMesh::Parallel_type::ALL, &ncells);
    
    // loop over cells joined to the vertex i
    for (int K = 0; K < ncells.size(); ++K) {
      
      // boundary conditions
      if (ncells[K] == -1) {
        U[0][i] = 0.5;
        U[1][i] = 0.0;
        U[2][i] = 0.0;
        phi_beta_cell[0] = 0.0;
        phi_beta_cell[1] = 0.0;
        phi_beta_cell[2] = 0.0;
        break;
      }
      else {
        
        AmanziMesh::Entity_ID_List cnodes;
        mesh_->cell_get_nodes(K, &cnodes);
        
        std::vector<std::vector<double>> phi(3); // phi[m][j]
        for (int m = 0; m < 3; ++m) {
          phi[m].resize(cnodes.size());
        }
        
        std::vector<double> Phi_total(3, 0.0); // Phi_total[m]
        std::vector<double> phi_tmp(3, 0.0); // phi[m][j] = phi_tmp[m]
        for (int j = 0; j < cnodes.size(); ++j) {
          phi_tmp = ResidualsLF(ncells[K], cnodes[j], U); // eq (10)
          
          for (int m = 0; m < 3; ++m) {
            Phi_total[m] += phi_tmp[m];     // eq(8)
            phi[m][j] = phi_tmp[m]; // eq (10)
          }
        }
        
        // calculate distribution coefficients beta
        std::vector<double> beta(3, 0.0); // beta[m]
        std::vector<double> sum_max(3, 0.0); // sum_max[m]
        for (int j = 0; j < cnodes.size(); ++j) {
          for (int m = 0; m < 3; ++m) {
            
            if (std::abs(Phi_total[m]) > 0.0) {
              
              if (cnodes[j] == i) {
                beta[m] = std::max(0.0, phi[m][j]/Phi_total[m]);
              }
              sum_max[m] += std::max(0.0, std::max(0.0, phi[m][j]/Phi_total[m]) );
            }
//            std::cout<<"sum_max: "<<m<<": "<<sum_max[m]<<std::endl;
//            std::cout<<"Phi_total: "<<m<<": "<<Phi_total[m]<<std::endl;
          }
        }
        
        for (int m = 0; m < 3; ++m) {
          if (std::abs(sum_max[m]) > 0.0) {
            beta[m] = beta[m] / sum_max[m];
          }
        }
        
        
        // calculate contribution to node i
        for (int m = 0; m < 3; ++m) {
          phi_beta_cell[m] += beta[m] * Phi_total[m]; // eq(6)
        }
        dual_cell_vol += mesh_->cell_volume(ncells[K]);
      }
    } // K (cell) loop
    
    for (int m = 0; m < 3; ++m) {
      U_pr[m][i] = -(dt/dual_cell_vol)*phi_beta_cell[m] + U[m][i];
    }
    
    std::cout<<"i: "<<i<<" out of "<<nnodes_owned<<std::endl;
    
  } // i (total DOF) loop
  
  for (int c = 0; c < ncells_owned; ++c) {
  
    AmanziMesh::Entity_ID_List cnodes;
    mesh_->cell_get_nodes(c, &cnodes);
    
    ht_c[0][c] = 0.0;
    h_c[0][c] = 0.0;
    vel_c[0][c] = 0.0;
    vel_c[1][c] = 0.0;
    q_c[0][c] = 0.0;
    q_c[1][c] = 0.0;
    
    for (int i = 0; i < cnodes.size(); ++i) {
      h_c[0][c] += U_pr[0][cnodes[i]];
      q_c[0][c] += U_pr[1][cnodes[i]];
      q_c[1][c] += U_pr[1][cnodes[i]];
    }
    
    h_c[0][c] = h_c[0][c]/cnodes.size();
    q_c[0][c] = q_c[0][c]/cnodes.size();
    q_c[1][c] = q_c[1][c]/cnodes.size();
    
    ht_c[0][c] = h_c[0][c] + B_c[0][c];
    double h = h_c[0][c];
    double factor = (2.0*h)/(h*h + std::max(h*h, 1.e-14));
    vel_c[0][c] = factor * q_c[0][c];
    vel_c[1][c] = factor * q_c[1][c];
    
    std::cout<<"h: "<<h_c[0][c]<<", v:"<<vel_c[0][c]<<", "<<vel_c[1][c]<<std::endl;
    
  }
  
  return failed;
  
}
  
  
  //--------------------------------------------------------------//--------------------------------------------------------------
  //--------------------------------------------------------------//--------------------------------------------------------------
  
//
//  int orientation;
//  AmanziMesh::Entity_ID_List cfaces, cnodes;
//
//  std::vector<double> FL, FR, FNum(3), FNum_rot, FS(3);  // fluxes
//  std::vector<double> S;                                 // source term
//  std::vector<double> UL(3), UR(3), U, U_new(3);         // numerical fluxes
//
//  // Simplest first-order form
//  // U_i^{n+1} = U_i^n - dt/vol * (F_{i+1/2}^n - F_{i-1/2}^n) + dt * S_i
//
//  for (int c = 0; c < ncells_owned; c++) {
//
//    const Amanzi::AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
//
//    double vol = mesh_->cell_volume(c);
//
//    mesh_->cell_get_faces(c, &cfaces);
//    mesh_->cell_get_nodes(c, &cnodes);
//
//    for (int i = 0; i < 3; i++) FS[i] = 0.0;
//
//    for (int n = 0; n < cfaces.size(); ++n) {
//      int f = cfaces[n];
//      double farea = mesh_->face_area(f);
//      const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);
//      AmanziGeometry::Point normal = mesh_->face_normal(f, false, c, &orientation);
//      normal /= farea;
//
//      double ht_rec = total_depth_grad_->getValue(c, xcf);
//
//      int edge = cfaces[n];
//      double B_rec = BathymetryEdgeValue(edge, B_n);
//
//      if (ht_rec < B_rec) {
//        ht_rec = ht_c[0][c];
//        B_rec = B_c[0][c];
//      }
//      double h_rec = ht_rec - B_rec;
//      failed = ErrorDiagnostics_(c, h_rec, B_rec, ht_rec);
//      if (failed) return failed;
//
//      double qx_rec = discharge_x_grad_->getValue(c, xcf);
//      double qy_rec = discharge_y_grad_->getValue(c, xcf);
//
//      double h2 = h_rec * h_rec;
//      double factor = 2.0 * h_rec / (h2 + std::fmax(h2, eps2));
//      double vx_rec = factor * qx_rec;
//      double vy_rec = factor * qy_rec;
//
//      // rotating velocity to the face-based coordinate system
//      double vn, vt;
//      vn =  vx_rec * normal[0] + vy_rec * normal[1];
//      vt = -vx_rec * normal[1] + vy_rec * normal[0];
//
//      UL[0] = h_rec;
//      UL[1] = h_rec * vn;
//      UL[2] = h_rec * vt;
//
//      int cn = WhetStone::cell_get_face_adj_cell(*mesh_, c, f);
//
//      if (cn == -1) {
//        if (bcs_.size() > 0 && bcs_[0]->bc_find(f)) {
//          for (int i = 0; i < 3; ++i) UR[i] = bcs_[0]->bc_value(f)[i];
//            vn =  UR[1] * normal[0] + UR[2] * normal[1];
//            vt = -UR[1] * normal[1] + UR[2] * normal[0];
//            UR[1] = UR[0] * vn;
//            UR[2] = UR[0] * vt;
//        }
//        else {
//          UR = UL;
//        }
//
//      } else {
//        const Amanzi::AmanziGeometry::Point& xcn = mesh_->cell_centroid(cn);
//
//        ht_rec = total_depth_grad_->getValue(cn, xcf);
//
//        B_rec = BathymetryEdgeValue(edge, B_n);
//
//        if (ht_rec < B_rec) {
//          ht_rec = ht_c[0][cn];
//          B_rec = B_c[0][cn];
//        }
//        h_rec = ht_rec - B_rec;
//        failed = ErrorDiagnostics_(cn, h_rec, B_rec, ht_rec);
//        if (failed) return failed;
//
//        qx_rec = discharge_x_grad_->getValue(cn, xcf);
//        qy_rec = discharge_y_grad_->getValue(cn, xcf);
//
//        h2 = h_rec * h_rec;
//        factor = 2.0 * h_rec / (h2 + std::fmax(h2, eps2));
//        vx_rec = factor * qx_rec;
//        vy_rec = factor * qy_rec;
//
//        vn =  vx_rec * normal[0] + vy_rec * normal[1];
//        vt = -vx_rec * normal[1] + vy_rec * normal[0];
//
//        UR[0] = h_rec;
//        UR[1] = h_rec * vn;
//        UR[2] = h_rec * vt;
//      }
//
//      FNum_rot = numerical_flux_->Compute(UL, UR);
//      // FNum_rot = NumericalFlux_x(UL, UR);
//
//      FNum[0] = FNum_rot[0];
//      FNum[1] = FNum_rot[1] * normal[0] - FNum_rot[2] * normal[1];
//      FNum[2] = FNum_rot[1] * normal[1] + FNum_rot[2] * normal[0];
//
//      // update accumulated cell-based flux
//      for (int i = 0; i < 3; i++) {
//        FS[i] += FNum[i] * farea;
//      }
//    }
//
//    double h, u, v, qx, qy;
//    U.resize(3);
//
//    h  = h_c[0][c];
//
//    qx = q_c[0][c];
//    qy = q_c[1][c];
//    double factor = 2.0 * h / (h * h + std::fmax(h * h, eps2));
//
//    u = factor * qx;
//    v = factor * qy;
//
//    U[0] = h;
//    U[1] = h * u;
//    U[2] = h * v;
//
//    S = NumericalSource(U, c);
//
//    for (int i = 0; i < 3; i++) {
//      U_new[i] = U[i] - dt / vol * FS[i] + dt * S[i];
//    }
//
//    U_new[0] += dt * ext_S_cell[c];
//
//    // transform to conservative variables
//    h  = U_new[0];
//    qx = U_new[1];
//    qy = U_new[2];
//
//    h_c_tmp[0][c] = h;
//    factor = 2.0 * h / (h * h + std::fmax(h * h, eps2));
//
//    vel_c_tmp[0][c] = factor * qx;
//    vel_c_tmp[1][c] = factor * qy;
//  } // c
//
//  h_c = h_c_tmp;
//  vel_c = vel_c_tmp;
//
//  for (int c = 0; c < ncells_owned; c++) {
//    ht_c[0][c] = h_c_tmp[0][c] + B_c[0][c];
//  }
//
//  return failed;
//}


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
std::vector<double> ShallowWater_PK::ResidualsLF(int K, int j, std::vector<std::vector<double> >& U) // input argument must contain coefficients of the basis expansion. cell c, node i, U[m][] contains the DOFs for mth component
{
  // calculate -\int_K F \cdot \nabla \phi_j + \int_{\partial K} (F \cdot n) \phi_j - \int_{K} S \phi_j + \alpha (U_m - Ubar) [eq(10)]
  
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  
  std::vector<double> integral(3, 0.0); // value to retrun
    
  AmanziMesh::Entity_ID_List cnodes; //DOFs of the cell (nodes)
  mesh_->cell_get_nodes(K, &cnodes);
  
  // 1a. construct volume quadrature points
  std::vector<AmanziGeometry::Point> quad_nodes_vol(3);
  std::vector<double> weights_vol(3);
  std::vector<AmanziGeometry::Point> coords(3);
  for (int i = 0; i < cnodes.size(); ++i) { // find coordinates of triangle vertices
    mesh_->node_get_coordinates(cnodes[i], &coords[i]);
  }
  // find physical coordinates of quadrature points (3 quadrature points for now; P1 triangle elements)
  for (int i = 0; i < cnodes.size(); ++i) {
    quad_nodes_vol[i] = (1 - WhetStone::q2d_points[i+1][0] - WhetStone::q2d_points[i+1][1] )*coords[0] + (WhetStone::q2d_points[i+1][0])*coords[1] + (WhetStone::q2d_points[i+1][1])*coords[2];
    weights_vol[i] = WhetStone::q2d_weights[i+1] * mesh_->cell_volume(K);
//    if (j == 1) {
//      std::cout<<"vertices: "<<coords[i][0]<<", "<<coords[i][1]<<std::endl;
//      std::cout<<"quad vol: "<<quad_nodes_vol[i][0]<<", "<<quad_nodes_vol[i][1]<<std::endl;
//      std::cout<<"weights: "<<weights_vol[i]<<std::endl;
//    }
  }
  
  // 1b. evaluate volume integral using quadrature points
  for (int qp = 0; qp < quad_nodes_vol.size(); ++qp ) {
    
    std::vector<double> Uqp(3);
    std::vector<std::vector<double>> flux(2);
    std::vector<double> phi_j_grad(2);
    double phi_j;
    
    AmanziGeometry::Point x_qp = quad_nodes_vol[qp];
    
    Uqp = EvalSol(U, x_qp);
    
    flux[0] = PhysFlux_x(Uqp);
    flux[1] = PhysFlux_y(Uqp);

    phi_j_grad = basis_grad(j, K, x_qp);
    phi_j = basis_value(j, K, x_qp);
//    std::cout<<"j, K, x_qp: "<<j<<", "<<K<<"; "<<x_qp[0]<<", "<<x_qp[1]<<std::endl;
//    std::cout<<"phi_j: "<<phi_j<<std::endl;

    for (int m = 0; m < 3; ++m) {
    integral[m] += -(flux[0][m]*phi_j_grad[0] + flux[1][m]*phi_j_grad[1]) * weights_vol[qp];
    }
//    integral += -PhysicalSource(Uqp) * phi_j * weights_vol[qp]; // bathymetry term
  }
    
  // 2a. construct face quadrature points
  std::vector<AmanziGeometry::Point> quad_nodes_face(2);
  std::vector<double> weights_face(2);
  AmanziMesh::Entity_ID_List cfaces;
  mesh_->cell_get_faces(K, &cfaces);
  
  for (int f = 0; f < cfaces.size(); ++f) {
    
    AmanziMesh::Entity_ID_List fnodes;
    mesh_->face_get_nodes(cfaces[f], &fnodes);
    
    std::vector<AmanziGeometry::Point> edge_coords(2);
    for (int i = 0; i < fnodes.size(); ++i) { // find coordinates of triangle edge nodes
      mesh_->node_get_coordinates(fnodes[i], &edge_coords[i]);
    }
    // find physical coordinates of face quadrature points (2 face quadrature points for now; P1 triangle elements)
    for (int i = 0; i < fnodes.size(); ++i) {
      quad_nodes_face[i] = (1 - WhetStone::q1d_points[1][i] )*edge_coords[0] + (WhetStone::q2d_points[1][i])*edge_coords[1];
      weights_face[i] = WhetStone::q1d_weights[1][i] * mesh_->face_area(cfaces[f]);
    }
    
    // 2b. evaluate face integral using quadrature points
    int orientation;
    AmanziGeometry::Point n = mesh_->face_normal(cfaces[f],false, K, &orientation);
    double farea = mesh_->face_area(cfaces[f]);
    n /= farea;
    
    // loop over quadrature points in dK
    for (int qp = 0; qp < quad_nodes_face.size(); ++qp) {
      
      std::vector<double> Uqp;
      std::vector<std::vector<double>> flux(2);
      double phi_i;

      AmanziGeometry::Point x_qp = quad_nodes_face[qp];

      Uqp = EvalSol(U,x_qp);
      
      flux[0] = PhysFlux_x(Uqp);
      flux[1] = PhysFlux_y(Uqp);

      phi_i = basis_value(j, K, x_qp);

      for (int m = 0; m < 3; ++m) {
        integral[m] += (flux[0][m]*n[0] + flux[1][m]*n[1]) * phi_i * weights_face[qp];
      }
    }
    
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
  double h, qx, qy, vx, vy, S, Smax = 0.0;
  for (int i = 0; i < cnodes.size(); ++i) {
    h = U[0][cnodes[i]];
    qx = U[1][cnodes[i]];
    qy = U[2][cnodes[i]];
    vx = 2.0 * h * qx / (h*h + std::max(h*h, 1.e-14));
    vy = 2.0 * h * qy / (h*h + std::max(h*h, 1.e-14));
    
    S = std::max(std::abs(vx) + std::sqrt(g_*h), std::abs(vy) + std::sqrt(g_*h));
    Smax = std::max(Smax, S);
  }
  
  double lmax = 0.0;
  for (int f = 0; f < cfaces.size(); ++f) {
    double farea = mesh_->face_area(cfaces[f]);
    lmax = std::max(lmax, farea);
  }
  
  double alpha = lmax * Smax;
  
  for (int m = 0; m < 3; ++m) {
    integral[m] += alpha * (U[m][j] - U_avg[m]);
  }
  
  return integral;
}


//--------------------------------------------------------------
// Time-space residuals for time-stepping scheme
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::ResidualsTimeSpace(int K, int j, std::vector<std::vector<double> >& U, std::vector<std::vector<double> >& U_pr, double dt)
{
  
  // compute time step integral \int_{K} ((u - u_star)/\Delta t)  * \phi_i dx
  std::vector<double> integral(3, 0.0);
  
  AmanziMesh::Entity_ID_List cnodes; //DOFs of the cell (nodes)
  mesh_->cell_get_nodes(K, &cnodes);
  
  // 1a. construct volume quadrature points
  std::vector<AmanziGeometry::Point> quad_nodes_vol(3);
  std::vector<double> weights_vol(3);
  std::vector<AmanziGeometry::Point> coords(3);
  for (int i = 0; i < cnodes.size(); ++i) { // find coordinates of triangle vertices
    mesh_->node_get_coordinates(cnodes[i], &coords[i]);
  }
  // find physical coordinates of quadrature points (3 quadrature points for now; P1 triangle elements)
  for (int i = 0; i < cnodes.size(); ++i) {
    quad_nodes_vol[i] = (1 - WhetStone::q2d_points[i+1][0] - WhetStone::q2d_points[i+1][1] )*coords[0] + (WhetStone::q2d_points[i+1][0])*coords[1] + (WhetStone::q2d_points[i+1][1])*coords[2];
    weights_vol[i] = WhetStone::q2d_weights[i+1];
  }
  
  // 1b. evaluate volume integral using quadrature points
  for (int qp = 0; qp < quad_nodes_vol.size(); ++qp ) {
    
    std::vector<double> Uqp(3);
    std::vector<double> U_prqp(3);
    double phi_j;
    
    AmanziGeometry::Point x_qp = quad_nodes_vol[qp];
    
    Uqp = EvalSol(U, x_qp);
    U_prqp = EvalSol(U, x_qp);

    phi_j = basis_value(j, K, x_qp);

    for (int m = 0; m < 3; ++m) {
    integral[m] += (1.0/dt) * (Uqp[m] - U_prqp[m]) * phi_j * weights_vol[qp];
    }
//    integral += -PhysicalSource(Uqp) * phi_j * weights_vol[qp]; // bathymetry term
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
std::vector<double> ShallowWater_PK::EvalSol(std::vector<std::vector<double>> U, AmanziGeometry::Point x_qp)
{
  std::vector<double> eval_sol(3, 0.0);
  
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int c;
  
  for (int K = 0; K < ncells_owned; ++K) { // find the cell where the quadrature point is
    AmanziMesh::Entity_ID_List cnodes;
    mesh_->cell_get_nodes (K, &cnodes);
    
    std::vector<AmanziGeometry::Point> coords(3);
    for (int i = 0; i < 3; ++i) {
      mesh_->node_get_coordinates(cnodes[i], &coords[i]);
    }
    
    if (AmanziGeometry::point_in_polygon(x_qp, coords) == true) {
      c = K;
      break;
    }
  }
  
  // find nodes of the cell c and loop over them to find \sum_{i \in K} U_i(t) \phi_i(x,y)
  AmanziMesh::Entity_ID_List cnodes;
  mesh_->cell_get_nodes(c, &cnodes);
  for (int m = 0; m < 3; ++m) {
    for (int i = 0; i < 3; ++i) {
      eval_sol[m] += U[m][cnodes[i]] * basis_value(cnodes[i], c, x_qp);
    }
  }
  
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
  
  for (int j = 0; j < 3; ++j) { // construct vertices of the plane over triangle
    mesh_->node_get_coordinates(cnodes[j], &x_vertex_tria[j]);
    AmanziGeometry::Point x_tmp = x_vertex_tria[j];
    double x0 = x_tmp[0], x1 = x_tmp[1], x2 = 0.0;

    if (cnodes[j] == i) {
      x2 = 1.0;
    }
    else {
      x2 = 0.0;
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
  
  return -( (x[0] - x1[0])*n[0] + (x[1] - x1[1])*n[1] )/ (n[2]); // (x1 - x) \cdot n = 0
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
  
  for (int j = 0; j < 3; ++j) { // construct vertices of the plane over triangle
    mesh_->node_get_coordinates(cnodes[j], &x_vertex_tria[j]);
    AmanziGeometry::Point x_tmp = x_vertex_tria[j];
    double x0 = x_tmp[0], x1 = x_tmp[1], x2 = 0.0;

    if (cnodes[j] == i) {
      x2 = 1.0;
    }
    else {
      x2 = 0.0;
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


//--------------------------------------------------------------------
// Discretization of the source term (well-balanced for lake at rest)
//--------------------------------------------------------------------
std::vector<double> ShallowWater_PK::NumericalSource(
    const std::vector<double>& U, int c)
{
  Epetra_MultiVector& B_c = *S_->GetFieldData(bathymetry_key_, passwd_)->ViewComponent("cell", true);
  Epetra_MultiVector& B_n = *S_->GetFieldData(bathymetry_key_, passwd_)->ViewComponent("node", true);
  Epetra_MultiVector& ht_c = *S_->GetFieldData(total_depth_key_, passwd_)->ViewComponent("cell", true);
  Epetra_MultiVector& h_c = *S_->GetFieldData(ponded_depth_key_, passwd_)->ViewComponent("cell", true);

  const Amanzi::AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  AmanziMesh::Entity_ID_List cfaces;
  mesh_->cell_get_faces(c, &cfaces);

  int orientation;
  double S1(0.0), S2(0.0);
  double vol = mesh_->cell_volume(c);

  for (int n = 0; n < cfaces.size(); ++n) {
    int f = cfaces[n];
    const auto& normal = mesh_->face_normal(f, false, c, &orientation);
    const auto& xcf = mesh_->face_centroid(f);
    double farea = mesh_->face_area(f);

    double ht_rec = total_depth_grad_->getValue(c, xcf);
      
    int e = cfaces[n];
    double B_rec = BathymetryEdgeValue(e, B_n);

    if (ht_rec < B_rec) {
      ht_rec = ht_c[0][c];
      B_rec = B_c[0][c];
    }
    
    // Kurganov Acta Numerica 2018
    S1 += (-B_rec * h_c[0][c] ) * normal[0];
    S2 += (-B_rec * h_c[0][c] ) * normal[1];
  }
    
  std::vector<double> S(3);
  S[0] = 0.0;
  S[1] = g_ / vol * S1;
  S[2] = g_ / vol * S2;

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
  u  = 2.*h*qx/(h*h + fmax(h*h,eps*eps));
  v  = 2.*h*qy/(h*h + fmax(h*h,eps*eps));

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
  u  = 2.*h*qx/(h*h + fmax(h*h,eps*eps));
  v  = 2.*h*qy/(h*h + fmax(h*h,eps*eps));

  // Form vector of y-fluxes G(U) = (hv, huv, hv^2 + 1/2 gh^2)

  G[0] = h*v;
  G[1] = h*u*v;
  G[2] = h*v*v+0.5*g*h*h;

  return G;
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

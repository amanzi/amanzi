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

// Amanzi::ShallowWater
#include "ShallowWater_PK.hh"
#include "DischargeEvaluator.hh"

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
    passwd_("state")
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

  //-------------------------------
  // constant fields
  //-------------------------------

  if (!S->HasField("gravity")) {
    S->RequireConstantVector("gravity", passwd_, 2);
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
}


//--------------------------------------------------------------
// Initilize internal data
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
    PK_DomainFunctionFactory<ShallowWaterBoundaryFunction > bc_factory(mesh_);

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

  // gravity
  Epetra_Vector& gvec = *S_->GetConstantVectorData("gravity", "state");
  g_ = std::fabs(gvec[1]);

  // reconstruction
  Teuchos::ParameterList plist = sw_list_->sublist("reconstruction");

  total_depth_grad_ = Teuchos::rcp(new Operators::ReconstructionCell(mesh_));
  total_depth_grad_->Init(plist);

  velocity_x_grad_ = Teuchos::rcp(new Operators::ReconstructionCell(mesh_));
  velocity_x_grad_->Init(plist);

  velocity_y_grad_ = Teuchos::rcp(new Operators::ReconstructionCell(mesh_));
  velocity_y_grad_->Init(plist);

  discharge_x_grad_ = Teuchos::rcp(new Operators::ReconstructionCell(mesh_));
  discharge_x_grad_->Init(plist);

  discharge_y_grad_ = Teuchos::rcp(new Operators::ReconstructionCell(mesh_));
  discharge_y_grad_->Init(plist);

  limiter_ = Teuchos::rcp(new Operators::LimiterCell(mesh_));
  limiter_->Init(plist);

  // default
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  InitializeField_(S_.ptr(), passwd_, bathymetry_key_, 0.0);
  InitializeField_(S_.ptr(), passwd_, ponded_depth_key_, 1.0);

  if (!S_->GetField(total_depth_key_, passwd_)->initialized()) {
    const auto& h_c = *S_->GetFieldData(ponded_depth_key_)->ViewComponent("cell");
    const auto& B_c = *S_->GetFieldData(bathymetry_key_)->ViewComponent("cell");
    auto& ht_c = *S_->GetFieldData(total_depth_key_, passwd_)->ViewComponent("cell");

    for (int c = 0; c < ncells_owned; c++) {
      ht_c[0][c] = h_c[0][c] + B_c[0][c];
    }

    S_->GetField(total_depth_key_, passwd_)->set_initialized();
  }

  InitializeField_(S_.ptr(), passwd_, velocity_key_, 0.0);
  InitializeField_(S_.ptr(), passwd_, discharge_key_, 0.0);

  // static fields
  S_->GetFieldData(bathymetry_key_)->ScatterMasterToGhosted("cell");

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

  bool failed = false;
  double eps2 = 1e-12;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  // save a copy of primary and conservative fields
  Epetra_MultiVector& B_c = *S_->GetFieldData(bathymetry_key_, passwd_)->ViewComponent("cell", true);
  Epetra_MultiVector& B_n = *S_->GetFieldData(bathymetry_key_, passwd_)->ViewComponent("node", true);
  Epetra_MultiVector& h_c = *S_->GetFieldData(ponded_depth_key_, passwd_)->ViewComponent("cell", true);
  Epetra_MultiVector& ht_c = *S_->GetFieldData(total_depth_key_, passwd_)->ViewComponent("cell", true);
  Epetra_MultiVector& vel_c = *S_->GetFieldData(velocity_key_, passwd_)->ViewComponent("cell", true);
    
  S_->GetFieldEvaluator(discharge_key_)->HasFieldChanged(S_.ptr(), passwd_);
  Epetra_MultiVector& q_c = *S_->GetFieldData(discharge_key_, discharge_key_)->ViewComponent("cell", true);

  // create copies of primary fields
  Epetra_MultiVector h_c_tmp(h_c);
  Epetra_MultiVector vel_c_tmp(vel_c);

  // distribute data to ghost cells
  S_->GetFieldData(total_depth_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(ponded_depth_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(velocity_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(discharge_key_)->ScatterMasterToGhosted("cell");

  // limited reconstructions
  auto tmp1 = S_->GetFieldData(total_depth_key_, passwd_)->ViewComponent("cell", true);
  total_depth_grad_->ComputeGradient(tmp1);
  limiter_->ApplyLimiter(tmp1, 0, total_depth_grad_->gradient());
  limiter_->gradient()->ScatterMasterToGhosted("cell");

  auto tmp3 = S_->GetFieldData(velocity_key_, passwd_)->ViewComponent("cell", true);
  velocity_x_grad_->ComputeGradient(tmp3, 0);
  limiter_->ApplyLimiter(tmp3, 0, velocity_x_grad_->gradient());
  limiter_->gradient()->ScatterMasterToGhosted("cell");

  auto tmp4 = S_->GetFieldData(velocity_key_, passwd_)->ViewComponent("cell", true);
  velocity_y_grad_->ComputeGradient(tmp4, 1);
  limiter_->ApplyLimiter(tmp4, 1, velocity_y_grad_->gradient());
  limiter_->gradient()->ScatterMasterToGhosted("cell");

  auto tmp5 = S_->GetFieldData(discharge_key_, discharge_key_)->ViewComponent("cell", true);
  discharge_x_grad_->ComputeGradient(tmp5, 0);
  limiter_->ApplyLimiter(tmp5, 0, discharge_x_grad_->gradient());
  limiter_->gradient()->ScatterMasterToGhosted("cell");

  auto tmp6 = S_->GetFieldData(discharge_key_, discharge_key_)->ViewComponent("cell", true);
  discharge_y_grad_->ComputeGradient(tmp6, 1);
  limiter_->ApplyLimiter(tmp6, 1, discharge_y_grad_->gradient());
  limiter_->gradient()->ScatterMasterToGhosted("cell");

  // update boundary conditions
  if (bcs_.size() > 0)
      bcs_[0]->Compute(t_old, t_new);

  // Shallow water equations have the form
  // U_t + F_x(U) + G_y(U) = S(U)

  int orientation;
  AmanziMesh::Entity_ID_List cfaces, cnodes;

  std::vector<double> FL, FR, FNum(3), FNum_rot, FS(3);  // fluxes
  std::vector<double> S;                                 // source term
  std::vector<double> UL(3), UR(3), U, U_new(3);         // numerical fluxes

  // Simplest first-order form
  // U_i^{n+1} = U_i^n - dt/vol * (F_{i+1/2}^n - F_{i-1/2}^n) + dt * S_i

  for (int c = 0; c < ncells_owned; c++) {
      
    const Amanzi::AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

    double vol = mesh_->cell_volume(c);

    mesh_->cell_get_faces(c, &cfaces);
    mesh_->cell_get_nodes(c, &cnodes);
    
    for (int i = 0; i < 3; i++) FS[i] = 0.0;

    for (int n = 0; n < cfaces.size(); ++n) {
      int f = cfaces[n];
      double farea = mesh_->face_area(f);
      const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);
      AmanziGeometry::Point normal = mesh_->face_normal(f, false, c, &orientation);
      normal /= farea;

      double ht_rec = total_depth_grad_->getValue(c, xcf);
  
      int edge = cfaces[n];
      double B_rec = Bathymetry_edge_value(c, edge, xcf, B_n);
        

      if (ht_rec < B_rec) {
        ht_rec = ht_c[0][c];
        B_rec = B_c[0][c];
      }
      double h_rec = ht_rec - B_rec;
      ErrorDiagnostics_(c, h_rec, ht_rec, B_rec);

      double vx_rec = velocity_x_grad_->getValue(c, xcf);
      double vy_rec = velocity_y_grad_->getValue(c, xcf);
      double qx_rec = discharge_x_grad_->getValue(c, xcf);
      double qy_rec = discharge_y_grad_->getValue(c, xcf);

      double h2 = h_rec * h_rec;
      double factor = 2.0 * h_rec / (h2 + std::fmax(h2, eps2));
      vx_rec = factor * qx_rec;
      vy_rec = factor * qy_rec;

      // rotating velovity to the face-based coordinate system
      double vn, vt;
      vn =  vx_rec * normal[0] + vy_rec * normal[1];
      vt = -vx_rec * normal[1] + vy_rec * normal[0];

      UL[0] = h_rec;
      UL[1] = h_rec * vn;
      UL[2] = h_rec * vt;

      int cn = WhetStone::cell_get_face_adj_cell(*mesh_, c, f);

      if (cn == -1) {
          if (bcs_.size() > 0 && bcs_[0]->bc_find(f)){
              for (int i = 0; i < 3; ++i){ UR[i] = bcs_[0]->bc_value(f)[i]; } }

        else
          UR = UL;

      }
      else {
          
        const Amanzi::AmanziGeometry::Point& xcn = mesh_->cell_centroid(cn);
          
        ht_rec = total_depth_grad_->getValue(cn, xcf);
          
        B_rec = Bathymetry_edge_value(cn, edge, xcf, B_n);
          
        if (ht_rec < B_rec) {
            std::cout<<"ht_rec < B_rec"<<std::endl;
          ht_rec = ht_c[0][cn];
          B_rec = B_c[0][cn];
        }
        h_rec = ht_rec - B_rec;
        ErrorDiagnostics_(cn, h_rec, ht_rec, B_rec);

        vx_rec = velocity_x_grad_->getValue(cn, xcf);
        vy_rec = velocity_y_grad_->getValue(cn, xcf);
        qx_rec = discharge_x_grad_->getValue(cn, xcf);
        qy_rec = discharge_y_grad_->getValue(cn, xcf);

        h2 = h_rec * h_rec;
        factor = 2.0 * h_rec / (h2 + std::fmax(h2, eps2));
        vx_rec = factor * qx_rec;
        vy_rec = factor * qy_rec;

        vn =  vx_rec * normal[0] + vy_rec * normal[1];
        vt = -vx_rec * normal[1] + vy_rec * normal[0];

        UR[0] = h_rec;
        UR[1] = h_rec * vn;
        UR[2] = h_rec * vt;
      }

      FNum_rot = NumericalFlux_x(UL, UR);

      FNum[0] = FNum_rot[0];
      FNum[1] = FNum_rot[1] * normal[0] - FNum_rot[2] * normal[1];
      FNum[2] = FNum_rot[1] * normal[1] + FNum_rot[2] * normal[0];
        

      // update accumulated cell-based flux
      for (int i = 0; i < 3; i++) {
        FS[i] += FNum[i] * farea;
      }
    }

    double h, u, v, qx, qy;
    U.resize(3);

    h  = h_c[0][c];
      
    qx = q_c[0][c];
    qy = q_c[1][c];
    double factor = 2.0 * h / (h * h + std::fmax(h * h, eps2));
    u = factor * qx;
    v = factor * qy;

    U[0] = h;
    U[1] = h * u;
    U[2] = h * v;

    S = NumericalSource(U, c);

    for (int i = 0; i < 3; i++) {
      U_new[i] = U[i] - dt / vol * FS[i] + dt * S[i];
    }
      
    // transform to conservative variables
    h  = U_new[0];
    qx = U_new[1];
    qy = U_new[2];

    h_c_tmp[0][c] = h;
    factor = 2.0 * h / (h * h + std::fmax(h * h, eps2));
    vel_c_tmp[0][c] = factor * qx;
    vel_c_tmp[1][c] = factor * qy;
      
  } // c

  h_c = h_c_tmp;
  vel_c = vel_c_tmp;

  for (int c = 0; c < ncells_owned; c++) {
    ht_c[0][c] = h_c_tmp[0][c] + B_c[0][c];
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


//==========================================================================
//
// Discretization: numerical fluxes, source terms, etc
//
//==========================================================================

//--------------------------------------------------------------
// minmod function
//--------------------------------------------------------------
double minmod(double a, double b)
{
  double m;

  if (a*b > 0) {
    if (std::fabs(a) < std::fabs(b)) {
      m = a;
    } else {
      m = b;
    }
  } else {
    m = 0.;
  }

  return m;
}


//--------------------------------------------------------------
// Physical flux F(U) in the x-direction:
// F(U) = (hu, hu^2 + 1/2 gh^2, huv)
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::PhysicalFlux_x(const std::vector<double>& U)
{
  std::vector<double> F(3);

  double h, h2, u, v;
  double eps2 = 1e-12;

  // transform from conservative (h, hu, hv) to primary (h, u, v) variables
  h = U[0];
  h2 = h * h;
  u = 2.0 * h * U[1] / (h2 + std::fmax(h2, eps2));
  v = 2.0 * h * U[2] / (h2 + std::fmax(h2, eps2));

  F[0] = h * u;
  F[1] = h * u * u + 0.5 * g_ * h2;
  F[2] = h * u * v;

  return F;
}


//--------------------------------------------------------------
// Physical flux G(U) in the y-direction. Since the SW system has
// rotational invariance, this flux is not used.
// G(U) = (hv, huv, hv^2 + 1/2 gh^2)
//--------------------------------------------------------------
/*
std::vector<double> ShallowWater_PK::PhysicalFlux_y(const std::vector<double>& U)
{
  std::vector<double> G(3);
  double h, u, v, qx, qy;
  double eps = 1.e-6;
  // transform from conservative (h, hu, hv) to primary (h, u, v) variables
  h  = U[0];
  qx = U[1];
  qy = U[2];
  u  = 2.*h*qx/(h*h + std::fmax(h*h,eps*eps));
  v  = 2.*h*qy/(h*h + std::fmax(h*h,eps*eps));
  G[0] = h*v;
  G[1] = h*u*v;
  G[2] = h*v*v+0.5*g_*h*h;
  return G;
}
*/


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
// Numerical flux in the x-direction
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::NumericalFlux_x(
    std::vector<double>& UL, std::vector<double>& UR)
{
  return NumericalFlux_x_Rusanov(UL, UR);
//  return NumericalFlux_x_CentralUpwind(UL, UR);
}


//--------------------------------------------------------------
// Rusanov numerical flux -- very simple but very diffusive
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::NumericalFlux_x_Rusanov(
    const std::vector<double>& UL, const std::vector<double>& UR)
{
  std::vector<double> FL, FR, F(3);

  double hL, uL, vL, hR, uR, vR, factor;
  double eps = 1.e-6;

  // SW conservative variables: (h, hu, hv)

  hL = UL[0];
  factor = 2.0 * hL / (hL * hL + std::fmax(hL * hL, eps * eps));
  uL = factor * UL[1];
  vL = factor * UL[2];

  hR = UR[0];
  factor = 2.0 * hR / (hR * hR + std::fmax(hR * hR, eps * eps));
  uR = factor * UR[1];
  vR = factor * UR[2];

  FL = PhysicalFlux_x(UL);
  FR = PhysicalFlux_x(UR);

  double SL, SR, Smax;

  SL = std::max(std::fabs(uL) + std::sqrt(g_*hL),std::fabs(vL) + std::sqrt(g_*hL));
  SR = std::max(std::fabs(uR) + std::sqrt(g_*hR),std::fabs(vR) + std::sqrt(g_*hR));

  Smax = std::max(SL,SR);

  for (int i = 0; i < 3; i++) {
    F[i] = 0.5*(FL[i]+FR[i]) - 0.5*Smax*(UR[i]-UL[i]);
  }

  return F;
}


//--------------------------------------------------------------
//Central-upwind numerical flux (Kurganov, Acta Numerica 2018)
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::NumericalFlux_x_CentralUpwind(
    const std::vector<double>& UL, const std::vector<double>& UR)
{
  std::vector<double> FL, FR, F(3), U_star(3), dU(3);

  double hL, uL, vL, hR, uR, vR, qxL, qyL, qxR, qyR;
  double apx, amx, apy, amy;
  double ap, am;
  double eps = 1.e-6;

  // SW conservative variables: (h, hu, hv)

  hL  = UL[0];
  qxL = UL[1];
  qyL = UL[2];
  uL  = 2.*hL*qxL/(hL*hL + std::fmax(hL*hL,eps*eps));
  vL  = 2.*hL*qyL/(hL*hL + std::fmax(hL*hL,eps*eps));

  hR  = UR[0];
  qxR = UR[1];
  qyR = UR[2];
  uR  = 2.*hR*qxR/(hR*hR + std::fmax(hR*hR,eps*eps));
  vR  = 2.*hR*qyR/(hR*hR + std::fmax(hR*hR,eps*eps));

  apx = std::max(std::max(std::fabs(uL)+std::sqrt(g_*hL),std::fabs(uR)+std::sqrt(g_*hR)),0.);
  apy = std::max(std::max(std::fabs(vL)+std::sqrt(g_*hL),std::fabs(vR)+std::sqrt(g_*hR)),0.);
  ap  = std::max(apx,apy);

  amx = std::min(std::min(std::fabs(uL)-std::sqrt(g_*hL),std::fabs(uR)-std::sqrt(g_*hR)),0.);
  amy = std::min(std::min(std::fabs(vL)-std::sqrt(g_*hL),std::fabs(vR)-std::sqrt(g_*hR)),0.);
  am  = std::min(amx,amy);

  FL = PhysicalFlux_x(UL);
  FR = PhysicalFlux_x(UR);

  for (int i = 0; i < 3; i++) {
    U_star[i] = (ap*UR[i] - am*UL[i] - (FR[i] - FL[i])) / (ap - am + eps);
  }

  for (int i = 0; i < 3; i++) {
    dU[i] = minmod(UR[i]-U_star[i],U_star[i]-UL[i]);
  }

  for (int i = 0; i < 3; i++) {
    F[i] = (ap*FL[i] - am*FR[i] + ap*am*(UR[i] - UL[i] - dU[i])) / (ap - am + eps);
  }

  return F;
}

 
//--------------------------------------------------------------------
// discretization of the source term (well-balanced for lake at rest)
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
    double B_rec = Bathymetry_edge_value(c, e, xcf, B_n);

    if (ht_rec < B_rec) {
      ht_rec = ht_c[0][c];
      B_rec = B_c[0][c];
    }

    S1 += (-B_rec * ht_rec + B_rec * B_rec / 2) * normal[0];
    S2 += (-B_rec * ht_rec + B_rec * B_rec / 2) * normal[1];
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
  double h, u, v, S, vol, dx, dt;
  double eps = 1.e-6;

  auto& h_c = *S_->GetFieldData(ponded_depth_key_)->ViewComponent("cell");
  auto& vel_c = *S_->GetFieldData(velocity_key_)->ViewComponent("cell");

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  dt = 1.e10;
  for (int c = 0; c < ncells_owned; c++) {
    h = h_c[0][c];
    u = vel_c[0][c];
    v = vel_c[1][c];
    S = std::max(std::fabs(u), std::fabs(v)) + std::sqrt(g_*h) + eps;
    vol = mesh_->cell_volume(c);
    dx = std::sqrt(vol);
    dt = std::min(dt, dx/S);
  }

  double dt_min;
  mesh_->get_comm()->MinAll(&dt, &dt_min, 1);

  if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "dt = " << dt << ", dt_min = " << dt_min << std::endl;
  }

  return cfl_ * dt_min;
}


//--------------------------------------------------------------
// Bathymetry (Linear construction for a rectangular cell)
//--------------------------------------------------------------
/*
double ShallowWater_PK::Bathymetry_rectangular_cell_value(int c, const AmanziGeometry::Point& xp, Epetra_MultiVector& B_n)
{
  double x = xp[0], y = xp[1];
    
  AmanziMesh::Entity_ID_List nodes, faces;
    
  mesh_->cell_get_faces(c, &faces);
  mesh_ -> cell_get_nodes(c, &nodes);
    
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

double ShallowWater_PK::Bathymetry_edge_value(int c, int e, const AmanziGeometry::Point& xp, Epetra_MultiVector& B_n)
{
  double x = xp[0], y = xp[1];
  AmanziMesh::Entity_ID_List face_nodes;
  mesh_ ->face_get_nodes(e, &face_nodes);

  double B_rec = ( B_n[0][face_nodes[0]] + B_n[0][face_nodes[1]] ) /2.0;
    
  return B_rec;
}


//--------------------------------------------------------------
// Error diagnostics
//--------------------------------------------------------------
void ShallowWater_PK::ErrorDiagnostics_(int c, double h, double B, double ht)
{
  if (h < 0.0) {
    Errors::Message msg;
    msg << "Shallow water PK: negative height: "
        << "\n  c = " << c
        << "\n  ht_rec = " << ht
        << "\n  B_rec  = " << B
        << "\n  h_rec  = " << h << "\n";
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace ShallowWater
}  // namespace Amanzi

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
#include "DischargeEvaluator.hh"
#include "HydrostaticPressureEvaluator.hh"
#include "NumericalFluxFactory.hh"
#include "ShallowWater_PK.hh"

namespace Amanzi {
namespace ShallowWater {

double inverse_with_tolerance(double h);

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
  total_depth_key_ = Keys::getKey(domain_, "total_depth");
  ponded_depth_key_ = Keys::getKey(domain_, "ponded_depth");
  prev_ponded_depth_key_ = Keys::getKey(domain_, "prev_ponded_depth");
  bathymetry_key_ = Keys::getKey(domain_, "bathymetry");
  hydrostatic_pressure_key_ = Keys::getKey(domain_, "ponded_pressure");

  riemann_flux_key_ = Keys::getKey(domain_, "riemann_flux");

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
  // -- ponded depth 
  if (!S->HasField(ponded_depth_key_)) {
    S->RequireField(ponded_depth_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    AddDefaultPrimaryEvaluator_(ponded_depth_key_);
  }

  // -- total depth
  if (!S->HasField(total_depth_key_)) {
    S->RequireField(total_depth_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // -- velocity
  if (!S->HasField(velocity_key_)) {
    S->RequireField(velocity_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 2);
    AddDefaultPrimaryEvaluator_(velocity_key_);
  }

  // -- discharge
  if (!S->HasField(discharge_key_)) {
    S->RequireField(discharge_key_, discharge_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 2);

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", discharge_key_);
    auto eval = Teuchos::rcp(new DischargeEvaluator(elist));
    S->SetFieldEvaluator(discharge_key_, eval);
  }

  // -- bathymetry
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
  // -- hydrostatic pressure
  if (!S->HasField(hydrostatic_pressure_key_)) {
    S->RequireField(hydrostatic_pressure_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    
    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", hydrostatic_pressure_key_);
    auto eval = Teuchos::rcp(new HydrostaticPressureEvaluator(elist));
    S->SetFieldEvaluator(hydrostatic_pressure_key_, eval);
  }

  // -- riemann flux
  if (!S->HasField(riemann_flux_key_)) {
    S->RequireField(riemann_flux_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);
  }

  // -- previous state of ponded depth (for coupling)
  if (!S->HasField(prev_ponded_depth_key_)) {
    S->RequireField(prev_ponded_depth_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->GetField(prev_ponded_depth_key_, passwd_)->set_io_vis(false);
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
        
    // compute cell averaged bathymery (Bc)
    Amanzi::AmanziGeometry::Point x0, x1;
    AmanziMesh::Entity_ID_List face_nodes;

    double tmp(0.0);
    for (int f = 0; f < nfaces_cell; ++f) {
      int edge = cfaces[f];
            
      mesh_->face_get_nodes(edge, &face_nodes);
                        
      mesh_->node_get_coordinates(face_nodes[0], &x0);
      mesh_->node_get_coordinates(face_nodes[1], &x1);

      AmanziGeometry::Point area = (xc - x0) ^ (xc - x1);
      tmp += norm(area) * (B_n[0][face_nodes[0]] + B_n[0][face_nodes[1]]) / 4;
    }
    B_c[0][c] = tmp / mesh_->cell_volume(c);
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

  InitializeField_(S_.ptr(), passwd_, riemann_flux_key_, 0.0);
  InitializeFieldFromField_(prev_ponded_depth_key_, ponded_depth_key_, false);
  
  // ---- ----
  Teuchos::RCP<TreeVector> tmp_h = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> tmp_q = Teuchos::rcp(new TreeVector());

  soln_->PushBack(tmp_h);
  soln_->PushBack(tmp_q);

  soln_h_ = S->GetFieldData(ponded_depth_key_, passwd_);
  soln_q_ = S->GetFieldData(discharge_key_, discharge_key_);

  soln_->SubVector(0)->SetData(soln_h_);
  soln_->SubVector(1)->SetData(soln_q_);
    
  // summary of initialization
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    double bmin, bmax; 
    B_c.MinValue(&bmin);
    B_c.MinValue(&bmax);
    *vo_->os() << "bathymetry min=" << bmin << ", max=" << bmax 
               << "\nShallow water PK was initialized." << std::endl;
  }
}


//--------------------------------------------------------------
// Auxiliary initialization technique.
//--------------------------------------------------------------
void ShallowWater_PK::InitializeFieldFromField_(
    const std::string& field0, const std::string& field1, bool call_evaluator)
{
  if (S_->HasField(field0)) {
    if (S_->GetField(field0)->owner() == passwd_) {
      if (!S_->GetField(field0, passwd_)->initialized()) {
        if (call_evaluator)
            S_->GetFieldEvaluator(field1)->HasFieldChanged(S_.ptr(), passwd_);

        const CompositeVector& f1 = *S_->GetFieldData(field1);
        CompositeVector& f0 = *S_->GetFieldData(field0, passwd_);
        f0 = f1;

        S_->GetField(field0, passwd_)->set_initialized();

        if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
            *vo_->os() << "initialized " << field0 << " to " << field1 << std::endl;
      }
    }
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

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  S_->GetFieldEvaluator(discharge_key_)->HasFieldChanged(S_.ptr(), passwd_);

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
  Epetra_MultiVector& riemann_f = *S_->GetFieldData(riemann_flux_key_, passwd_)->ViewComponent("face", true);
    
  S_->GetFieldEvaluator(discharge_key_)->HasFieldChanged(S_.ptr(), passwd_);
  Epetra_MultiVector& q_c = *S_->GetFieldData(discharge_key_, discharge_key_)->ViewComponent("cell", true);

  // create copies of primary fields
  *S_->GetFieldData(prev_ponded_depth_key_, passwd_)->ViewComponent("cell", true) = h_c;
  
  Epetra_MultiVector& h_old = *soln_->SubVector(0)->Data()->ViewComponent("cell");
  Epetra_MultiVector& q_old = *soln_->SubVector(1)->Data()->ViewComponent("cell");
  
  for (int c = 0; c < ncells_owned; c++) {
    double factor = inverse_with_tolerance(h_old[0][c]);
    h_c[0][c] = h_old[0][c];
    q_c[0][c] = q_old[0][c];
    q_c[1][c] = q_old[1][c];
    vel_c[0][c] = factor * q_old[0][c];
    vel_c[1][c] = factor * q_old[1][c];
    ht_c[0][c] = h_c[0][c] + B_c[0][c];
  }

  // Shallow water equations have the form
  // U_t + F_x(U) + G_y(U) = S(U)
    
  // initialize time integrator
  auto ti_method = Explicit_TI::forward_euler;
//  auto ti_method = Explicit_TI::midpoint;
//  auto ti_method = Explicit_TI::tvd_3rd_order;
//  auto ti_method = Explicit_TI::runge_kutta_4th_order;
  
  auto soln_new = Teuchos::rcp(new TreeVector(*soln_));
  *soln_new = *soln_;
  Explicit_TI::RK<TreeVector> rk1(*this, ti_method, *soln_);
  rk1.TimeStep(t_old, dt, *soln_, *soln_new);
  *soln_ = *soln_new;
  
  Epetra_MultiVector& h_temp = *soln_->SubVector(0)->Data()->ViewComponent("cell");
  Epetra_MultiVector& q_temp = *soln_->SubVector(1)->Data()->ViewComponent("cell");
  
  // update solution
  for (int c = 0; c < ncells_owned; c++) {
    double factor = inverse_with_tolerance(h_temp[0][c]);
    h_c[0][c] = h_temp[0][c];
    q_c[0][c] = q_temp[0][c];
    q_c[1][c] = q_temp[1][c];
    vel_c[0][c] = factor * q_temp[0][c];
    vel_c[1][c] = factor * q_temp[1][c];
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
  S_->GetFieldEvaluator(hydrostatic_pressure_key_)->HasFieldChanged(S_.ptr(), passwd_);
  
  Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(S->GetFieldEvaluator(velocity_key_))->SetFieldAsChanged(S.ptr());
  Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(S->GetFieldEvaluator(ponded_depth_key_))->SetFieldAsChanged(S.ptr());
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

    double ht_rec = total_depth_grad_->getValue(c, xcf);
    double B_rec = BathymetryEdgeValue(f, B_n);

    if (ht_rec < B_rec) {
      ht_rec = ht_c[0][c];
      B_rec = B_c[0][c];
    }
    
    // Polygonal meshes Beljadid et. al. 2016
    S1 += (0.5) * (ht_rec - B_rec) * (ht_rec - B_rec) * normal[0];
    S2 += (0.5) * (ht_rec - B_rec) * (ht_rec - B_rec) * normal[1];
  }
  
  Epetra_MultiVector& ht_grad = *total_depth_grad_->gradient()->ViewComponent("cell", true);
  
  S1 /= vol;
  S2 /= vol;
  S1 -= ht_grad[0][c] * U[0];
  S2 -= ht_grad[1][c] * U[0];
  S1 *= g_;
  S2 *= g_;

  std::vector<double> S(3);
  
  S[0] = 0.0;
  S[1] = S1;
  S[2] = S2;
  
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

//  double vol = std::sqrt(mesh_->cell_volume(0));
//  return 1.e-4*vol;
  
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
// Inversion operation protected for small values
//--------------------------------------------------------------
double inverse_with_tolerance(double h)
{
  double eps(1e-6), eps2(1e-12);  // hard-coded tolerances

  if (h > eps) return 1.0 / h;

  double h2 = h * h;
  return 2 * h / (h2 + std::fmax(h2, eps2));
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

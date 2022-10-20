/*
 Shallow water PK

 Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
 Amanzi is released under the three-clause BSD License.
 The terms of use and "as is" disclaimer for this license are
 provided in the top-level COPYRIGHT file.

 Author: Svetlana Tokareva (tokareva@lanl.gov)
*/

#include <algorithm>
#include <cfloat>
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

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

//--------------------------------------------------------------
// Standard constructor
//--------------------------------------------------------------
ShallowWater_PK::ShallowWater_PK(Teuchos::ParameterList& pk_tree,
                                 const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                 const Teuchos::RCP<State>& S,
                                 const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, glist, S, soln), 
    glist_(glist),
    soln_(soln),
    S_(S),
    passwd_(""),
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
void
ShallowWater_PK::Setup()
{
  mesh_ = S_->GetMesh(domain_);
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
  if (!S_->HasRecord("gravity")) {
    S_->Require<AmanziGeometry::Point>("gravity", Tags::DEFAULT, "state");
  }

  // required for calculating hydrostatic pressure
  if (!S_->HasRecord("const_fluid_density")) {
    S_->Require<double>("const_fluid_density", Tags::DEFAULT, "state");
  }

  if (!S_->HasRecord("atmospheric_pressure")) {
    S_->Require<double>("atmospheric_pressure", Tags::DEFAULT, "state");
  }

  //-------------------------------
  // primary fields
  //-------------------------------
  // -- ponded depth
  if (!S_->HasRecord(ponded_depth_key_)) {
    S_->Require<CV_t, CVS_t>(ponded_depth_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    AddDefaultPrimaryEvaluator_(ponded_depth_key_);
  }

  // -- total depth
  if (!S_->HasRecord(total_depth_key_)) {
    S_->Require<CV_t, CVS_t>(total_depth_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // -- velocity
  if (!S_->HasRecord(velocity_key_)) {
    S_->Require<CV_t, CVS_t>(velocity_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 2);
    AddDefaultPrimaryEvaluator_(velocity_key_);
  }

  // -- discharge
  if (!S_->HasRecord(discharge_key_)) {
    S_->Require<CV_t, CVS_t>(discharge_key_, Tags::DEFAULT, discharge_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 2);

    Teuchos::ParameterList elist(discharge_key_);
    elist.set<std::string>("my key", discharge_key_).set<std::string>("tag", Tags::DEFAULT.get());
    auto eval = Teuchos::rcp(new DischargeEvaluator(elist));
    S_->SetEvaluator(discharge_key_, Tags::DEFAULT, eval);
  }

  // -- bathymetry
  if (!S_->HasRecord(bathymetry_key_)) {
    std::vector<std::string> names({ "cell", "node" });
    std::vector<int> ndofs(2, 1);
    std::vector<AmanziMesh::Entity_kind> locations({ AmanziMesh::CELL, AmanziMesh::NODE });

    S_->Require<CV_t, CVS_t>(bathymetry_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);
  }

  //-------------------------------
  // secondary fields
  //-------------------------------
  // -- hydrostatic pressure
  if (!S_->HasRecord(hydrostatic_pressure_key_)) {
    S_->Require<CV_t, CVS_t>(hydrostatic_pressure_key_, Tags::DEFAULT, hydrostatic_pressure_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist(hydrostatic_pressure_key_);
    elist.set<std::string>("my key", hydrostatic_pressure_key_).set<std::string>("tag", "");
    auto eval = Teuchos::rcp(new HydrostaticPressureEvaluator(elist));
    S_->SetEvaluator(hydrostatic_pressure_key_, Tags::DEFAULT, eval);
  }

  // -- riemann flux
  if (!S_->HasRecord(riemann_flux_key_)) {
    S_->Require<CV_t, CVS_t>(riemann_flux_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);
  }

  // -- previous state of ponded depth (for coupling)
  if (!S_->HasRecord(prev_ponded_depth_key_)) {
    S_->Require<CV_t, CVS_t>(prev_ponded_depth_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->GetRecordW(prev_ponded_depth_key_, passwd_).set_io_vis(false);
  }
}


//--------------------------------------------------------------
// Initialize internal data
//--------------------------------------------------------------
void
ShallowWater_PK::Initialize()
{
  // Create BC objects
  Teuchos::RCP<ShallowWaterBoundaryFunction> bc;
  Teuchos::RCP<Teuchos::ParameterList> bc_list =
    Teuchos::rcp(new Teuchos::ParameterList(sw_list_->sublist("boundary conditions", false)));

  bcs_.clear();

  // velocity BC is required on the faces while ponded depth is required on
  // nodes
  // -- velocity BC
  if (bc_list->isSublist("velocity")) {
    PK_DomainFunctionFactory<ShallowWaterBoundaryFunction> bc_factory(mesh_, S_);

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

  // -- ponded depth BC
  if (bc_list->isSublist("ponded depth")) {
    PK_DomainFunctionFactory<ShallowWaterBoundaryFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("ponded depth");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);

        bc = bc_factory.Create(spec, "ponded depth", AmanziMesh::NODE, Teuchos::null);
        bc->set_bc_name("ponded depth");
        bc->set_type(WhetStone::DOF_Type::SCALAR);
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
  g_ = norm(S_->Get<AmanziGeometry::Point>("gravity"));

  // numerical flux
  Teuchos::ParameterList model_list;
  model_list.set<std::string>("numerical flux", sw_list_->get<std::string>("numerical flux", "central upwind"))
    .set<double>("gravity", g_);
  NumericalFluxFactory nf_factory;
  numerical_flux_ = nf_factory.Create(model_list);

  // reconstruction
  Teuchos::ParameterList plist = sw_list_->sublist("reconstruction");

  total_depth_grad_ = Teuchos::rcp(new Operators::ReconstructionCellLinear(mesh_));
  total_depth_grad_->Init(plist);

  discharge_x_grad_ = Teuchos::rcp(new Operators::ReconstructionCellLinear(mesh_));
  discharge_x_grad_->Init(plist);

  discharge_y_grad_ = Teuchos::rcp(new Operators::ReconstructionCellLinear(mesh_));
  discharge_y_grad_->Init(plist);

  use_limiter_ = sw_list_->get<bool>("use limiter", true);
  if (use_limiter_ == true) {
    limiter_ = Teuchos::rcp(new Operators::LimiterCell(mesh_));
    limiter_->Init(plist);
  }

  // default
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  if (!S_->GetRecord(bathymetry_key_).initialized()) {
    InitializeCVField(S_, *vo_, bathymetry_key_, Tags::DEFAULT, passwd_, 0.0);
  }

  const auto& B_n = *S_->Get<CV_t>(bathymetry_key_).ViewComponent("node");
  auto& B_c = *S_->GetW<CV_t>(bathymetry_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");

  // compute B_c from B_n for well balanced scheme (Beljadid et. al. 2016)
  S_->Get<CV_t>(bathymetry_key_).ScatterMasterToGhosted("node");

  double cell_area_max = 0.0;
  for (int c = 0; c < ncells_owned; ++c) {
    const Amanzi::AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

    cell_area_max = std::max(cell_area_max, mesh_->cell_volume(c));
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
  S_->Get<CV_t>(bathymetry_key_).ScatterMasterToGhosted("cell");

  // calculate cell area square (used as a tolerance)  
  cell_area2_max_ = cell_area_max * cell_area_max;

  // initialize h from ht or ht from h
  if (!S_->GetRecord(ponded_depth_key_, Tags::DEFAULT).initialized()) {
    auto& h_c = *S_->GetW<CV_t>(ponded_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");
    auto& ht_c = *S_->GetW<CV_t>(total_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");

    for (int c = 0; c < ncells_owned; c++) { h_c[0][c] = ht_c[0][c] - B_c[0][c]; }

    S_->GetRecordW(ponded_depth_key_, Tags::DEFAULT, passwd_).set_initialized();
  }

  if (!S_->GetRecord(total_depth_key_).initialized()) {
    const auto& h_c = *S_->Get<CV_t>(ponded_depth_key_).ViewComponent("cell");
    auto& ht_c = *S_->GetW<CV_t>(total_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");

    for (int c = 0; c < ncells_owned; c++) { ht_c[0][c] = h_c[0][c] + B_c[0][c]; }

    S_->GetRecordW(total_depth_key_, Tags::DEFAULT, passwd_).set_initialized();
  }

  InitializeCVField(S_, *vo_, velocity_key_, Tags::DEFAULT, passwd_, 0.0);
  InitializeCVField(S_, *vo_, discharge_key_, Tags::DEFAULT, discharge_key_, 0.0);

  // secondary fields
  S_->GetEvaluator(hydrostatic_pressure_key_).Update(*S_, passwd_);

  InitializeCVField(S_, *vo_, riemann_flux_key_, Tags::DEFAULT, passwd_, 0.0);
  InitializeFieldFromField_(prev_ponded_depth_key_, ponded_depth_key_, false);

  // soln_ is the TreeVector of conservative variables [h hu hv]
  Teuchos::RCP<TreeVector> tmp_h = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> tmp_q = Teuchos::rcp(new TreeVector());

  soln_->PushBack(tmp_h);
  soln_->PushBack(tmp_q);

  auto soln_h = S_->GetPtrW<CV_t>(ponded_depth_key_, Tags::DEFAULT, passwd_);
  auto soln_q = S_->GetPtrW<CV_t>(discharge_key_, Tags::DEFAULT, discharge_key_);

  soln_->SubVector(0)->SetData(soln_h);
  soln_->SubVector(1)->SetData(soln_q);

  // temporal discretization order
  temporal_disc_order_ = sw_list_->get<int>("temporal discretization order", 1);

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
void
ShallowWater_PK::InitializeFieldFromField_(const std::string& field0,
                                           const std::string& field1,
                                           bool call_evaluator)
{
  if (S_->HasRecord(field0)) {
    if (S_->GetRecord(field0).owner() == passwd_) {
      if (!S_->GetRecord(field0).initialized()) {
        if (call_evaluator) S_->GetEvaluator(field1).Update(*S_, passwd_);

        const auto& f1 = S_->Get<CV_t>(field1);
        auto& f0 = S_->GetW<CV_t>(field0, Tags::DEFAULT, passwd_);
        f0 = f1;

        S_->GetRecordW(field0, passwd_).set_initialized();

        Teuchos::OSTab tab = vo_->getOSTab();
        if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
          *vo_->os() << "initialized " << field0 << " to " << field1 << std::endl;
      }
    }
  }
}


//--------------------------------------------------------------
// Advance conservative variables: (h, hu, hv)
//--------------------------------------------------------------
bool
ShallowWater_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  double dt = t_new - t_old;
  iters_++;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  S_->GetEvaluator(discharge_key_).Update(*S_, passwd_);

  // distribute data to ghost cells
  S_->Get<CV_t>(total_depth_key_).ScatterMasterToGhosted("cell");
  S_->Get<CV_t>(ponded_depth_key_).ScatterMasterToGhosted("cell");
  S_->Get<CV_t>(velocity_key_).ScatterMasterToGhosted("cell");
  S_->Get<CV_t>(discharge_key_).ScatterMasterToGhosted("cell");

  // save a copy of primary and conservative fields
  auto& B_c = *S_->GetW<CV_t>(bathymetry_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
  auto& h_c = *S_->GetW<CV_t>(ponded_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
  auto& ht_c = *S_->GetW<CV_t>(total_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
  auto& vel_c = *S_->GetW<CV_t>(velocity_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);

  S_->GetEvaluator(discharge_key_).Update(*S_, passwd_);
  auto& q_c = *S_->GetW<CV_t>(discharge_key_, Tags::DEFAULT, discharge_key_).ViewComponent("cell", true);

  // create copies of primary fields
  StateArchive archive(S_, vo_);
  archive.Add({ ponded_depth_key_ }, { discharge_key_ },
              { ponded_depth_key_, velocity_key_ }, 
              Tags::DEFAULT, "shallow_water");

  Epetra_MultiVector& h_old = *soln_->SubVector(0)->Data()->ViewComponent("cell");
  Epetra_MultiVector& q_old = *soln_->SubVector(1)->Data()->ViewComponent("cell");

  for (int c = 0; c < ncells_owned; ++c) {
    double factor = inverse_with_tolerance(h_old[0][c], cell_area2_max_);
    h_c[0][c] = h_old[0][c];
    q_c[0][c] = q_old[0][c];
    q_c[1][c] = q_old[1][c];
    vel_c[0][c] = factor * q_old[0][c];
    vel_c[1][c] = factor * q_old[1][c];
    ht_c[0][c] = h_c[0][c] + B_c[0][c];
  }

  // update source (external) terms
  for (int i = 0; i < srcs_.size(); ++i) { srcs_[i]->Compute(t_old, t_new); }

  // compute total source value for each time step
  total_source_ = 0.0;
  for (int i = 0; i < srcs_.size(); ++i) {
    for (auto it = srcs_[i]->begin(); it != srcs_[i]->end(); ++it) {
      int c = it->first;
      total_source_ += it->second[0] * mesh_->cell_volume(c) * dt; // data unit is [m^3]
    }
  }

  // Shallow water equations have the form
  // U_t + F_x(U) + G_y(U) = S(U)

  // initialize time integrator
  auto ti_method = Explicit_TI::forward_euler;
  if (temporal_disc_order_ == 2) {
    ti_method = Explicit_TI::midpoint;
  } else if (temporal_disc_order_ == 3) {
    ti_method = Explicit_TI::tvd_3rd_order;
  } else if (temporal_disc_order_ == 4) {
    ti_method = Explicit_TI::runge_kutta_4th_order;
  }

  // To use evaluators, we need to overwrite state variables. Since,
  // soln_ is a shallow copy of state variables, we cannot use it.
  // Instead, we need its deep copy. This is by design of the DAG.
  int ierr(0);
  try {
    auto soln_old = Teuchos::rcp(new TreeVector(*soln_));
    auto soln_new = Teuchos::rcp(new TreeVector(*soln_));

    Explicit_TI::RK<TreeVector> rk1(*this, ti_method, *soln_old);
    rk1.TimeStep(t_old, dt, *soln_old, *soln_new);

    *soln_ = *soln_new;
    VerifySolution_(*soln_); 
  } catch (...) {
    ierr = -1;
  }

  // recover state in case of error
  if (ierr < 0) {
    archive.Restore(passwd_);
    return true;
  }

  Epetra_MultiVector& h_temp = *soln_->SubVector(0)->Data()->ViewComponent("cell");
  Epetra_MultiVector& q_temp = *soln_->SubVector(1)->Data()->ViewComponent("cell");

  // update solution
  for (int c = 0; c < ncells_owned; ++c) {
    double factor = inverse_with_tolerance(h_temp[0][c], cell_area2_max_);
    h_c[0][c] = h_temp[0][c];
    q_c[0][c] = q_temp[0][c];
    q_c[1][c] = q_temp[1][c];
    vel_c[0][c] = factor * q_temp[0][c];
    vel_c[1][c] = factor * q_temp[1][c];
    ht_c[0][c] = h_c[0][c] + B_c[0][c];
  }

  // For consistency with other flow models, we need to track previous h
  // which was placed earlier in the archive.
  S_->GetW<CV_t>(prev_ponded_depth_key_, Tags::DEFAULT, passwd_) = archive.get(ponded_depth_key_);

  Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CV_t, CVS_t>>(
    S_->GetEvaluatorPtr(ponded_depth_key_, Tags::DEFAULT))->SetChanged();

  Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CV_t, CVS_t>>(
    S_->GetEvaluatorPtr(velocity_key_, Tags::DEFAULT))->SetChanged();

  // min-max values
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    double hmin(DBL_MAX), hmax(DBL_MIN);
    for (int c = 0; c < ncells_owned; ++c) {
      hmin = std::min(hmin, h_c[0][c]);
      hmax = std::max(hmax, h_c[0][c]);
    }

    double qmin(DBL_MAX), qmax(DBL_MIN);
    auto& riemann_f = *S_->GetW<CompositeVector>(riemann_flux_key_, passwd_).ViewComponent("face");
    for (int f = 0; f < nfaces_owned; ++f) {
      double flux = riemann_f[0][f] / mesh_->face_area(f);
      qmin = std::min(qmin, flux);
      qmax = std::max(qmax, flux);
    }

    double outmin[2], inmin[2] = { hmin, qmin };
    double outmax[2], inmax[2] = { hmax, qmax };
    mesh_->get_comm()->MinAll(inmin, outmin, 2);
    mesh_->get_comm()->MaxAll(inmax, outmax, 2);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "min/max(h): " << outmin[0] << "/" << outmax[0] 
               << ", min/max(flux): " << outmin[1] << "/" << outmax[1] << std::endl;
  }

  return false;
}


//--------------------------------------------------------------
// Advance conservative variables: (h, hu, hv)
//--------------------------------------------------------------
void
ShallowWater_PK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  S_->GetEvaluator(hydrostatic_pressure_key_).Update(*S_, passwd_);
}


//--------------------------------------------------------------
// Total Depth ht = h + B (Evaluate value at edge midpoint for a polygonal cell)
// Reconstruct if necessary for positivity
//--------------------------------------------------------------
double
ShallowWater_PK::TotalDepthEdgeValue(
    int c, int e, double htc, double Bc, const Epetra_MultiVector& B_n)
{
  double ht_edge(0.0); // value to return

  auto& ht_grad = *total_depth_grad_->data()->ViewComponent("cell", true);
  // set ht_grad = 0 for first order reconstruction (for debugging purposes)
  // ht_grad[0][c] = 0.0;
  // ht_grad[1][c] = 0.0;

  const auto& xc = mesh_->cell_centroid(c);
  const auto& xf = mesh_->face_centroid(e);
  Amanzi::AmanziMesh::Entity_ID_List cnodes, cfaces;
  mesh_->cell_get_nodes(c, &cnodes);

  bool cell_is_dry, cell_is_fully_flooded, cell_is_partially_wet;
  cell_is_partially_wet = false;
  cell_is_dry = false;
  cell_is_fully_flooded = false;

  // calculate maximum bathymetry value on cell nodes
  double Bmax = 0.0;
  for (int i = 0; i < cnodes.size(); ++i) { 
    Bmax = std::max(B_n[0][cnodes[i]], Bmax); 
  }

  // characterize cell based on [Beljadid et al.' 16]
  if ((htc > Bmax) && (htc - Bc > 0.0)) {
    cell_is_fully_flooded = true;
  } else if (std::abs(htc - Bc) < 1.e-15) {
    cell_is_dry = true;
  } else {
    cell_is_partially_wet = true;
  }

  // depth poisitivity based on [Beljadid et al.' 2016, Computers and Fluids]
  if (cell_is_fully_flooded) {
    ht_edge = total_depth_grad_->getValue(c, xf);
    if (ht_edge - BathymetryEdgeValue(e, B_n) < 0.0) { 
      ht_edge = htc; 
    }
  } else if (cell_is_dry) {
    ht_edge = BathymetryEdgeValue(e, B_n);
  } else if (cell_is_partially_wet) {
    Amanzi::AmanziMesh::Entity_ID_List cfaces;
    mesh_->cell_get_faces(c, &cfaces);

    double mu_eps_sum = 0.0;

    for (int f = 0; f < cfaces.size(); ++f) {
      Amanzi::AmanziGeometry::Point x0, x1;
      int edge = cfaces[f];

      Amanzi::AmanziMesh::Entity_ID_List face_nodes;
      mesh_->face_get_nodes(edge, &face_nodes);
      int n0 = face_nodes[0], n1 = face_nodes[1];

      mesh_->node_get_coordinates(n0, &x0);
      mesh_->node_get_coordinates(n1, &x1);

      double area = norm((xc - x0) ^ (xc - x1)) / 2.0;

      double epsilon;

      if ((htc < B_n[0][n0]) && (htc < B_n[0][n1])) {
        epsilon = 0.0;
      } else if ((htc >= B_n[0][n0]) && (htc >= B_n[0][n1])) {
        epsilon = 1.0;
      } else {
        epsilon = 0.5;
      }

      mu_eps_sum += (area / mesh_->cell_volume(c)) * (epsilon);
    }

    Amanzi::AmanziMesh::Entity_ID_List face_nodes;
    mesh_->face_get_nodes(e, &face_nodes);

    ht_edge = 0.0;
    for (int i = 0; i < face_nodes.size(); ++i) {
      if (htc < B_n[0][face_nodes[i]]) {
        ht_edge += B_n[0][face_nodes[i]];
      } else {
        ht_edge += B_n[0][face_nodes[i]] + ((htc - Bc) / mu_eps_sum);
      }
    }
    ht_edge /= 2.0;
    ht_grad[0][c] = 0.0;
    ht_grad[1][c] = 0.0;
  }

  return ht_edge;
}


//--------------------------------------------------------------------
// Discretization of the source term (well-balanced for lake at rest)
//--------------------------------------------------------------------
std::vector<double>
ShallowWater_PK::NumericalSource(int c, double htc, double Bc, const Epetra_MultiVector& B_n)
{
  AmanziMesh::Entity_ID_List cfaces, cnodes;
  mesh_->cell_get_faces(c, &cfaces);
  mesh_->cell_get_nodes(c, &cnodes);

  int orientation;
  double S1(0.0), S2(0.0);
  double vol = mesh_->cell_volume(c);

  for (int n = 0; n < cfaces.size(); ++n) {
    int f = cfaces[n];
    const auto& normal = mesh_->face_normal(f, false, c, &orientation);
		
    double ht_rec = TotalDepthEdgeValue(c, f, htc, Bc, B_n);
    double B_rec = BathymetryEdgeValue(f, B_n);

    // Polygonal meshes [Beljadid et al.' 2016, Computers and Fluids]
    S1 += (0.5) * (ht_rec - B_rec) * (ht_rec - B_rec) * normal[0];
    S2 += (0.5) * (ht_rec - B_rec) * (ht_rec - B_rec) * normal[1];
  }

  auto& ht_grad = *total_depth_grad_->data()->ViewComponent("cell", true);

  S1 /= vol;
  S2 /= vol;
  S1 -= ht_grad[0][c] * (htc - Bc);
  S2 -= ht_grad[1][c] * (htc - Bc);
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
double
ShallowWater_PK::get_dt()
{
  double d, d_min = 1.e10, vn, dt = 1.e10, dt_dry = 1.e-1;

  const auto& h_c = *S_->Get<CV_t>(ponded_depth_key_).ViewComponent("cell", true);
  const auto& vel_c = *S_->Get<CV_t>(velocity_key_).ViewComponent("cell", true);

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
      d_min = std::min(d_min, d);

      dt = std::min(d / std::max((2.0 * (std::abs(vn) + std::sqrt(g_ * h))), 1.e-12), dt);
    }
  }

  // reduce dt_min when dt is too large for completely dry conditions (h = 0, qx = 0, qy = 0)
  if (dt >= d_min * 1.e8) { dt = d_min * dt_dry; }

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

  if (iters_ < max_iters_) {
    return 0.1 * cfl_ * dt_min;
  } else {
    return cfl_ * dt_min;
  }
}


//--------------------------------------------------------------
// Throws an error if solution u is not valid.
//--------------------------------------------------------------
void ShallowWater_PK::VerifySolution_(TreeVector& u) 
{
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  const auto& h_c = *u.SubVector(0)->Data()->ViewComponent("cell");

  int ierr(0);
  for (int c = 0; c < ncells_owned; ++c) {
    if (h_c[0][c] < 0.0) {
      ierr = -1;
      break;
    }
  }

  int ierr_tmp(ierr);
  mesh_->get_comm()->MinAll(&ierr_tmp, &ierr, 1);
  if (ierr < 0) {
    Errors::Message msg;
    msg << "Negative ponded depth.\n";
    Exceptions::amanzi_throw(msg);
  }
}


//--------------------------------------------------------------
// Bathymetry (Evaluate value at edge midpoint for a polygonal cell)
//--------------------------------------------------------------
double
ShallowWater_PK::BathymetryEdgeValue(int e, const Epetra_MultiVector& B_n)
{
  AmanziMesh::Entity_ID_List nodes;
  mesh_->face_get_nodes(e, &nodes);

  return (B_n[0][nodes[0]] + B_n[0][nodes[1]]) / 2.0;
}


//--------------------------------------------------------------
// Inversion operation protected for small values
//--------------------------------------------------------------
double
inverse_with_tolerance(double h, double tol)
{
  double eps(1.e-4); // hard-coded tolerances

  if (h > eps) return 1.0 / h;

  double h2 = h * h;
  return 2 * h / (h2 + std::max(h2, tol));
}

} // namespace ShallowWater
} // namespace Amanzi

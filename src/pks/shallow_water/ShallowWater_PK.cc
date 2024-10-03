/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Svetlana Tokareva (tokareva@lanl.gov)
*/

/*
 Shallow water PK

 Authors: Svetlana Tokareva (tokareva@lanl.gov)
          Giacomo Capodaglio (gcapodaglio@lanl.gov)
          Naren Vohra (vohra@lanl.gov)
*/

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>

#include "CompositeVector.hh"
#include "PK_DomainFunctionFactory.hh"
#include "StateArchive.hh"

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
  : PK(pk_tree, glist, S, soln), glist_(glist), soln_(soln), S_(S), passwd_(""), iters_(0)
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
  cfl_positivity_ = sw_list_->get<double>("depth positivity cfl", 0.95);
  velocity_desingularization_eps_ = sw_list_->get<double>("desingularization epsilon", 1.e-04);
  source_key_ = sw_list_->get<std::string>("source key", "");

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
  dim_ = mesh_->getSpaceDimension();

  // domain name
  velocity_key_ = Keys::getKey(domain_, "velocity");
  discharge_key_ = Keys::getKey(domain_, "discharge");
  total_depth_key_ = Keys::getKey(domain_, "total_depth");
  bathymetry_key_ = Keys::getKey(domain_, "bathymetry");
  hydrostatic_pressure_key_ = Keys::getKey(domain_, "ponded_pressure");
  riemann_flux_key_ = Keys::getKey(domain_, "riemann_flux");
  SetupPrimaryVariableKeys();
  SetupExtraEvaluatorsKeys();

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
  // -- ponded depth for shallow water, wetted area for pipe flow
  if (!S_->HasRecord(primary_variable_key_)) {
    S_->Require<CV_t, CVS_t>(primary_variable_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    AddDefaultPrimaryEvaluator(S_, primary_variable_key_);
  }

  // -- velocity
  if (!S_->HasRecord(velocity_key_)) {
    S_->Require<CV_t, CVS_t>(velocity_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 2);
    AddDefaultPrimaryEvaluator(S_, velocity_key_);
  }

  // -- discharge
  if (!S_->HasRecord(discharge_key_)) {
    S_->Require<CV_t, CVS_t>(discharge_key_, Tags::DEFAULT, discharge_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 2);

    Teuchos::ParameterList elist(discharge_key_);
    elist.set<std::string>("my key", discharge_key_).set<std::string>("tag", "");
    auto eval = Teuchos::rcp(new DischargeEvaluator(elist));
    S_->SetEvaluator(discharge_key_, Tags::DEFAULT, eval);
  }

  //-------------------------------
  // secondary fields
  //-------------------------------

  // -- pipe drain
  if (!source_key_.empty()) {
    if (!S_->HasRecord(source_key_)) {
      S_->Require<CV_t, CVS_t>(source_key_, Tags::DEFAULT, source_key_)
        .SetMesh(mesh_)
        ->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);

      S_->RequireEvaluator(source_key_, Tags::DEFAULT);
    }
  }

  // -- bathymetry
  if (!S_->HasRecord(bathymetry_key_)) {
    std::vector<std::string> names({ "cell", "node", "cmax" });
    std::vector<int> ndofs({ 1, 1, 1 });
    std::vector<AmanziMesh::Entity_kind> locations({ AmanziMesh::Entity_kind::CELL,
                                                     AmanziMesh::Entity_kind::NODE,
                                                     AmanziMesh::Entity_kind::CELL });

    S_->Require<CV_t, CVS_t>(bathymetry_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);
  }

  // -- hydrostatic pressure
  if (!S_->HasRecord(hydrostatic_pressure_key_)) {
    S_->Require<CV_t, CVS_t>(hydrostatic_pressure_key_, Tags::DEFAULT, hydrostatic_pressure_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

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
      ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  }

  // -- total depth
  if (!S_->HasRecord(total_depth_key_)) {
    S_->Require<CV_t, CVS_t>(total_depth_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // -- previous state of primary variable (for coupling)
  if (!S_->HasRecord(prev_primary_variable_key_)) {
    S_->Require<CV_t, CVS_t>(prev_primary_variable_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S_->GetRecordW(prev_primary_variable_key_, passwd_).set_io_vis(false);
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

  //primary variable BC is required on nodes
  // primary variable BC
  SetPrimaryVariableBC(bc_list);

  // velocity BC is required on faces
  // -- velocity BC
  if (bc_list->isSublist("velocity")) {
    PK_DomainFunctionFactory<ShallowWaterBoundaryFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("velocity");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);

        bc = bc_factory.Create(
          spec, "velocity", AmanziMesh::Entity_kind::FACE, Teuchos::null, Tags::DEFAULT, true);
        bc->set_bc_name("velocity");
        bc->set_type(WhetStone::DOF_Type::VECTOR);
        bcs_.push_back(bc);
      }
    }
  }

  // -- discharge BC
  if (bc_list->isSublist("discharge")) {
    PK_DomainFunctionFactory<ShallowWaterBoundaryFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("discharge");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);

        bc = bc_factory.Create(spec, "discharge", AmanziMesh::Entity_kind::FACE, Teuchos::null);
        bc->set_bc_name("discharge");
        bc->set_type(WhetStone::DOF_Type::VECTOR);
        bcs_.push_back(bc);
      }
    }
  }


  if (bc_list->isSublist("outward discharge")) {
    // in the xml, a negative sign on the BC will mean inflow of water
    PK_DomainFunctionFactory<ShallowWaterBoundaryFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("outward discharge");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);

        bc = bc_factory.Create(
          spec, "outward discharge", AmanziMesh::Entity_kind::FACE, Teuchos::null);
        bc->set_bc_name("outward discharge");
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

        srcs_.push_back(
          factory.Create(spec, "source", AmanziMesh::Entity_kind::CELL, Teuchos::null));
      }
    }
  }

  // gravity
  g_ = norm(S_->Get<AmanziGeometry::Point>("gravity"));

  // numerical flux
  Teuchos::ParameterList model_list;
  model_list
    .set<std::string>("numerical flux",
                      sw_list_->get<std::string>("numerical flux", "central upwind"))
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
  int ncells_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  if (!S_->GetRecord(bathymetry_key_).initialized()) {
    InitializeCVField(S_, *vo_, bathymetry_key_, Tags::DEFAULT, passwd_, 0.0);
  }

  const auto& B_n = *S_->Get<CV_t>(bathymetry_key_).ViewComponent("node");
  auto& B_c = *S_->GetW<CV_t>(bathymetry_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");
  auto& B_max = *S_->GetW<CV_t>(bathymetry_key_, Tags::DEFAULT, passwd_).ViewComponent("cmax");

  // compute B_c from B_n for well balanced scheme (Beljadid et. al. 2016)
  S_->Get<CV_t>(bathymetry_key_).ScatterMasterToGhosted("node");

  double cell_area_max = 0.0;
  for (int c = 0; c < ncells_owned; ++c) {
    const Amanzi::AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);

    cell_area_max = std::max(cell_area_max, mesh_->getCellVolume(c));
    auto cfaces = mesh_->getCellFaces(c);
    int nfaces_cell = cfaces.size();

    // compute cell averaged bathymery (Bc)
    Amanzi::AmanziGeometry::Point x0, x1;
    double tmp(0.0);
    for (int f = 0; f < nfaces_cell; ++f) {
      int edge = cfaces[f];

      auto face_nodes = mesh_->getFaceNodes(edge);
      int n0 = face_nodes[0];
      int n1 = face_nodes[1];

      x0 = mesh_->getNodeCoordinate(n0);
      x1 = mesh_->getNodeCoordinate(n1);

      AmanziGeometry::Point area = (xc - x0) ^ (xc - x1);
      tmp += norm(area) * (B_n[0][n0] + B_n[0][n1]) / 4;

      B_max[0][c] = std::max(B_max[0][c], B_n[0][n0]);
      B_max[0][c] = std::max(B_max[0][c], B_n[0][n1]);
    }
    B_c[0][c] = tmp / mesh_->getCellVolume(c);
  }
  // redistribute the result
  S_->Get<CV_t>(bathymetry_key_).ScatterMasterToGhosted("cell");
  S_->Get<CV_t>(bathymetry_key_).ScatterMasterToGhosted("cmax");

  // calculate cell area square (used as a tolerance)
  cell_area2_max_ = cell_area_max * cell_area_max;

  InitializeFields();

  InitializeCVField(S_, *vo_, velocity_key_, Tags::DEFAULT, passwd_, 0.0);
  InitializeCVField(S_, *vo_, discharge_key_, Tags::DEFAULT, discharge_key_, 0.0);

  // secondary fields
  S_->GetEvaluator(hydrostatic_pressure_key_).Update(*S_, passwd_);

  InitializeCVField(S_, *vo_, riemann_flux_key_, Tags::DEFAULT, passwd_, 0.0);
  InitializeCVFieldFromCVField(
    S_, *vo_, prev_primary_variable_key_, primary_variable_key_, passwd_);

  // soln_ is the TreeVector of conservative variables [h hu hv]
  Teuchos::RCP<TreeVector> tmp_h = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> tmp_q = Teuchos::rcp(new TreeVector());

  soln_->PushBack(tmp_h);
  soln_->PushBack(tmp_q);

  auto soln_h = S_->GetPtrW<CV_t>(primary_variable_key_, Tags::DEFAULT, passwd_);
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
// Advance conservative variables: (h, hu, hv)
//--------------------------------------------------------------
bool
ShallowWater_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  double dt = t_new - t_old;
  iters_++;

  int ncells_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  S_->GetEvaluator(discharge_key_).Update(*S_, passwd_);

  // distribute data to ghost cells
  S_->Get<CV_t>(total_depth_key_).ScatterMasterToGhosted("cell");
  S_->Get<CV_t>(primary_variable_key_).ScatterMasterToGhosted("cell");
  S_->Get<CV_t>(velocity_key_).ScatterMasterToGhosted("cell");
  S_->Get<CV_t>(discharge_key_).ScatterMasterToGhosted("cell");
  ScatterMasterToGhostedExtraEvaluators();

  // save a copy of primary and conservative fields
  auto& B_c = *S_->GetW<CV_t>(bathymetry_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
  auto& h_c =
    *S_->GetW<CV_t>(primary_variable_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
  auto& ht_c =
    *S_->GetW<CV_t>(total_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
  auto& vel_c = *S_->GetW<CV_t>(velocity_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);

  S_->GetEvaluator(discharge_key_).Update(*S_, passwd_);
  auto& q_c =
    *S_->GetW<CV_t>(discharge_key_, Tags::DEFAULT, discharge_key_).ViewComponent("cell", true);

  // create copies of primary fields
  std::vector<std::string> fields({ primary_variable_key_, velocity_key_, discharge_key_ });
  StateArchive archive(S_, vo_);
  archive.Add(fields, Tags::DEFAULT);

  Epetra_MultiVector& h_old = *soln_->SubVector(0)->Data()->ViewComponent("cell");
  Epetra_MultiVector& q_old = *soln_->SubVector(1)->Data()->ViewComponent("cell");

  for (int c = 0; c < ncells_owned; ++c) {
    double factor = inverse_with_tolerance(h_old[0][c], cell_area2_max_);
    h_c[0][c] = h_old[0][c];
    q_c[0][c] = q_old[0][c];
    q_c[1][c] = q_old[1][c];
    vel_c[0][c] = factor * q_old[0][c];
    vel_c[1][c] = factor * q_old[1][c];
  }
  UpdateSecondaryFields();

  // update source (external) terms
  for (int i = 0; i < srcs_.size(); ++i) { srcs_[i]->Compute(t_old, t_new); }

  // compute total source value for each time step
  total_source_ = 0.0;
  for (int i = 0; i < srcs_.size(); ++i) {
    for (auto it = srcs_[i]->begin(); it != srcs_[i]->end(); ++it) {
      int c = it->first;
      total_source_ += it->second[0] * mesh_->getCellVolume(c) * dt; // data unit is [m^3]
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
    cfl_ *= 0.75;
    return true;
  }

  // For consistency with other flow models, we need to track previous h
  // which was placed earlier in the archive.
  S_->GetW<CV_t>(prev_primary_variable_key_, Tags::DEFAULT, passwd_) =
    archive.get(primary_variable_key_);

  return false;
}


//--------------------------------------------------------------
// Advance conservative variables: (h, hu, hv)
//--------------------------------------------------------------
void
ShallowWater_PK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  int ncells_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  int nfaces_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);

  auto& B_c = *S_->GetW<CV_t>(bathymetry_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
  auto& h_c =
    *S_->GetW<CV_t>(primary_variable_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
  auto& ht_c =
    *S_->GetW<CV_t>(total_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
  auto& q_c =
    *S_->GetW<CV_t>(discharge_key_, Tags::DEFAULT, discharge_key_).ViewComponent("cell", true);
  auto& vel_c = *S_->GetW<CV_t>(velocity_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);

  Epetra_MultiVector& h_new = *soln_->SubVector(0)->Data()->ViewComponent("cell");
  Epetra_MultiVector& q_new = *soln_->SubVector(1)->Data()->ViewComponent("cell");

  // update solution
  for (int c = 0; c < ncells_owned; ++c) {
    double factor = inverse_with_tolerance(h_new[0][c], cell_area2_max_);
    h_c[0][c] = h_new[0][c];
    q_c[0][c] = q_new[0][c];
    q_c[1][c] = q_new[1][c];
    vel_c[0][c] = factor * q_new[0][c];
    vel_c[1][c] = factor * q_new[1][c];
  }
  UpdateSecondaryFields();

  if (!source_key_.empty()) {
    S_->GetEvaluator(source_key_).Update(*S_, passwd_); // in this evaluator when running pipe flow
  }

  Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CV_t, CVS_t>>(
    S_->GetEvaluatorPtr(primary_variable_key_, Tags::DEFAULT))
    ->SetChanged();

  Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CV_t, CVS_t>>(
    S_->GetEvaluatorPtr(velocity_key_, Tags::DEFAULT))
    ->SetChanged();

  S_->GetEvaluator(hydrostatic_pressure_key_).Update(*S_, passwd_);

  UpdateExtraEvaluators();

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
      double flux = riemann_f[0][f] / mesh_->getFaceArea(f);
      qmin = std::min(qmin, flux);
      qmax = std::max(qmax, flux);
    }

    double outmin[2], inmin[2] = { hmin, qmin };
    double outmax[2], inmax[2] = { hmax, qmax };
    mesh_->getComm()->MinAll(inmin, outmin, 2);
    mesh_->getComm()->MaxAll(inmax, outmax, 2);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "min/max(h): " << outmin[0] << "/" << outmax[0]
               << ", min/max(flux): " << outmin[1] << "/" << outmax[1] << std::endl;
  }
}


//--------------------------------------------------------------
// Total Depth ht = h + B (Evaluate value at edge midpoint for a polygonal cell)
// Reconstruct if necessary for positivity
//--------------------------------------------------------------
double
ShallowWater_PK::TotalDepthEdgeValue(int c,
                                     int e,
                                     double htc,
                                     double Bc,
                                     double Bmax,
                                     const Epetra_MultiVector& B_n)
{
  double ht_edge(0.0); // value to return

  auto& ht_grad = *total_depth_grad_->data()->ViewComponent("cell", true);

  const auto& xc = mesh_->getCellCentroid(c);
  const auto& xf = mesh_->getFaceCentroid(e);

  bool cell_is_dry, cell_is_fully_flooded, cell_is_partially_wet;
  cell_is_partially_wet = false;
  cell_is_dry = false;
  cell_is_fully_flooded = false;

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
  } else if (cell_is_dry) {
    ht_edge = BathymetryEdgeValue(e, B_n);
  } else if (cell_is_partially_wet) {
    auto cfaces = mesh_->getCellFaces(c);

    double mu_eps_sum = 0.0;

    for (int f = 0; f < cfaces.size(); ++f) {
      Amanzi::AmanziGeometry::Point x0, x1;
      int edge = cfaces[f];

      auto face_nodes = mesh_->getFaceNodes(edge);
      int n0 = face_nodes[0], n1 = face_nodes[1];

      x0 = mesh_->getNodeCoordinate(n0);
      x1 = mesh_->getNodeCoordinate(n1);

      double area = norm((xc - x0) ^ (xc - x1)) / 2.0;

      double epsilon;

      if ((htc < B_n[0][n0]) && (htc < B_n[0][n1])) {
        epsilon = 0.0;
      } else if ((htc >= B_n[0][n0]) && (htc >= B_n[0][n1])) {
        epsilon = 1.0;
      } else {
        epsilon = 0.5;
      }

      mu_eps_sum += (area / mesh_->getCellVolume(c)) * (epsilon);
    }

    auto face_nodes = mesh_->getFaceNodes(e);

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
// Discretization of the bed slope source term (well-balanced for lake at rest)
// To be used for the shallow water model only
//--------------------------------------------------------------------
std::vector<double>
ShallowWater_PK::NumericalSourceBedSlope(int c,
                                         double htc,
                                         double Bc,
                                         double Bmax,
                                         const Epetra_MultiVector& B_n)
{
  auto cfaces = mesh_->getCellFaces(c);
  auto cnodes = mesh_->getCellNodes(c);

  int orientation;
  double S1(0.0), S2(0.0);
  double vol = mesh_->getCellVolume(c);

  for (int n = 0; n < cfaces.size(); ++n) {
    int f = cfaces[n];
    const auto& normal = mesh_->getFaceNormal(f, c, &orientation);

    double ht_rec = TotalDepthEdgeValue(c, f, htc, Bc, Bmax, B_n);
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
  int dir;
  double d, d_min = 1.e10, vn, dt = 1.e10, dt_dry = 1.e-1;
  double h, vx, vy;
  double vnMax = 0.0;
  double hMax = 0.0;

  const auto& h_c = *S_->Get<CV_t>(primary_variable_key_).ViewComponent("cell", true);
  const auto& vel_c = *S_->Get<CV_t>(velocity_key_).ViewComponent("cell", true);

  int ncells_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    const Amanzi::AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);

    auto cfaces = mesh_->getCellFaces(c);

    for (int n = 0; n < cfaces.size(); ++n) {
      int f = cfaces[n];
      double farea = mesh_->getFaceArea(f);
      const auto& xf = mesh_->getFaceCentroid(f);
      const auto& normal = mesh_->getFaceNormal(f);

      h = h_c[0][c];
      vx = vel_c[0][c];
      vy = vel_c[1][c];

      // computing local (cell, face) time step using Kurganov's estimate d / (2a)
      vn = (vx * normal[0] + vy * normal[1]) / farea;
      d = norm(xc - xf);
      d_min = std::min(d_min, d);

      dt = std::min(d / std::max((2.0 * (std::abs(vn) + std::sqrt(g_ * h))), 1.e-12), dt);
      vnMax = std::max(std::abs(vn), vnMax);
      hMax = std::max(h, hMax);
    }
  }

  // reduce dt_min when dt is too large for completely dry conditions (h = 0, qx = 0, qy = 0)
  if (dt >= d_min * 1.e8) { dt = d_min * dt_dry; }

  double dt_min;
  mesh_->getComm()->MinAll(&dt, &dt_min, 1);

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "stable dt = " << dt_min << ", cfl = " << cfl_ << std::endl;
    *vo_->os() << "max abs vel normal to face = " << vnMax << ", max primary variable = " << hMax
               << std::endl;
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
void
ShallowWater_PK::VerifySolution_(TreeVector& u)
{
  int ncells_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  const auto& h_c = *u.SubVector(0)->Data()->ViewComponent("cell");

  int ierr(0);
  for (int c = 0; c < ncells_owned; ++c) {
    if (h_c[0][c] < 0.0) {
      ierr = -1;
      break;
    }
  }

  int ierr_tmp(ierr);
  mesh_->getComm()->MinAll(&ierr_tmp, &ierr, 1);
  if (ierr < 0) {
    Errors::Message msg;
    msg << "Negative primary variable.\n";
    Exceptions::amanzi_throw(msg);
  }
}


//--------------------------------------------------------------
// Bathymetry (evaluate value at edge midpoint for a polygonal cell)
//--------------------------------------------------------------
double
ShallowWater_PK::BathymetryEdgeValue(int e, const Epetra_MultiVector& B_n)
{
  auto nodes = mesh_->getFaceNodes(e);

  return (B_n[0][nodes[0]] + B_n[0][nodes[1]]) / 2.0;
}

//--------------------------------------------------------------
// Setup primary variable keys
//--------------------------------------------------------------
void
ShallowWater_PK::SetupPrimaryVariableKeys()
{
  primary_variable_key_ = Keys::getKey(domain_, "ponded_depth");
  prev_primary_variable_key_ = Keys::getKey(domain_, "prev_ponded_depth");
}

//--------------------------------------------------------------
// Set primary variable BC
//--------------------------------------------------------------
void
ShallowWater_PK::SetPrimaryVariableBC(Teuchos::RCP<Teuchos::ParameterList>& bc_list)
{
  // -- ponded depth BC
  if (bc_list->isSublist("ponded depth")) {
    PK_DomainFunctionFactory<ShallowWaterBoundaryFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("ponded depth");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);

        Teuchos::RCP<ShallowWaterBoundaryFunction> bc = bc_factory.Create(
          spec, "ponded depth", AmanziMesh::Entity_kind::NODE, Teuchos::null, Tags::DEFAULT, true);
        bc->set_bc_name("ponded depth");
        bc->set_type(WhetStone::DOF_Type::SCALAR);
        bcs_.push_back(bc);
      }
    }
  }
}

//--------------------------------------------------------------------
// Update total depth for shallow water model
//--------------------------------------------------------------------
void
ShallowWater_PK::UpdateSecondaryFields()
{
  auto& h_c =
    *S_->GetW<CV_t>(primary_variable_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
  auto& ht_c =
    *S_->GetW<CV_t>(total_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
  auto& B_c = *S_->GetW<CV_t>(bathymetry_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
  int ncells_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  for (int c = 0; c < ncells_owned; ++c) { ht_c[0][c] = ComputeTotalDepth(h_c[0][c], B_c[0][c]); }
}

//--------------------------------------------------------------
// Initialize fields
//--------------------------------------------------------------
void
ShallowWater_PK::InitializeFields()
{
  int ncells_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  auto& B_c = *S_->GetW<CV_t>(bathymetry_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");

  if (!S_->GetRecord(primary_variable_key_, Tags::DEFAULT).initialized()) {
    auto& h_c =
      *S_->GetW<CV_t>(primary_variable_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");
    auto& ht_c = *S_->GetW<CV_t>(total_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");

    for (int c = 0; c < ncells_owned; c++) { h_c[0][c] = ht_c[0][c] - B_c[0][c]; }

    S_->GetRecordW(primary_variable_key_, Tags::DEFAULT, passwd_).set_initialized();
  }

  if (!S_->GetRecord(total_depth_key_).initialized()) {
    const auto& h_c = *S_->Get<CV_t>(primary_variable_key_).ViewComponent("cell");
    auto& ht_c = *S_->GetW<CV_t>(total_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");

    for (int c = 0; c < ncells_owned; c++) { ht_c[0][c] = ComputeTotalDepth(h_c[0][c], B_c[0][c]); }

    S_->GetRecordW(total_depth_key_, Tags::DEFAULT, passwd_).set_initialized();
  }
}

//--------------------------------------------------------------
// Compute external forcing on cells
//--------------------------------------------------------------
void
ShallowWater_PK::ComputeExternalForcingOnCells(std::vector<double>& forcing)
{
  for (int i = 0; i < srcs_.size(); ++i) {
    for (auto it = srcs_[i]->begin(); it != srcs_[i]->end(); ++it) {
      int c = it->first;
      forcing[c] = it->second[0] / mesh_->getCellVolume(c); // [m/s] for SW
    }
  }
}

//--------------------------------------------------------------------
// Compute ponded depth at edge location
//--------------------------------------------------------------------
double
ShallowWater_PK::ComputeFieldOnEdge(int c,
                                    int e,
                                    double htc,
                                    double Bc,
                                    double Bmax,
                                    const Epetra_MultiVector& B_n)
{
  double ht_rec = TotalDepthEdgeValue(c, e, htc, Bc, Bmax, B_n);
  double B_rec = BathymetryEdgeValue(e, B_n);
  double V_rec = ht_rec - B_rec;

  return V_rec;
}

//--------------------------------------------------------------
// Inversion operation protected for small values
//--------------------------------------------------------------
double
ShallowWater_PK::inverse_with_tolerance(double h, double tol)
{
  double eps(velocity_desingularization_eps_);

  if (h > eps) return 1.0 / h;

  double h2 = h * h;
  return 2 * h / (h2 + std::max(h2, tol));
}

} // namespace ShallowWater
} // namespace Amanzi

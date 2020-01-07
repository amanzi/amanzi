/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Multiphase multi-component flow. We assume that one component (water)
  is always in liquid form, i,e does not evaporate. This allows us to
  compare results with other codes. Later, this PK will be split into
  a base and derived PKs depending on physical models.
*/


// TPLs
#include "Teuchos_RCP.hpp"

// Amanzi
#include "LinearOperatorFactory.hh"
#include "PDE_Accumulation.hh"
#include "PDE_AdvectionUpwind.hh"
#include "PDE_DiffusionFVwithGravity.hh"
#include "PK_DomainFunctionFactory.hh"
#include "primary_variable_field_evaluator.hh"
#include "RelPermEvaluator.hh"

// Multiphase
#include "ModelMeshPartition.hh"
#include "Multiphase_PK.hh"
#include "Multiphase_Utils.hh"
#include "MultiphaseTypeDefs.hh"
#include "PressureGasEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Standard constructor
****************************************************************** */
Multiphase_PK::Multiphase_PK(Teuchos::ParameterList& pk_tree,
                             const Teuchos::RCP<Teuchos::ParameterList>& glist,
                             const Teuchos::RCP<State>& S,
                             const Teuchos::RCP<TreeVector>& soln) :
  glist_(glist),
  S_(S),
  soln_(soln),
  passwd_("multiphase"),
  num_phases_(2),
  num_itrs_(0)
{
  std::string pk_name = pk_tree.name();
  auto found = pk_name.rfind("->");
  if (found != std::string::npos) pk_name.erase(0, found + 2);

  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  mp_list_ = Teuchos::sublist(pk_list, pk_name, true);

  // We also need miscaleneous sublists
  pc_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(mp_list_, "time integrator", true);

  // computational domain
  domain_ = mp_list_->template get<std::string>("domain name", "domain");
}

  
/* ******************************************************************
* Setup
****************************************************************** */
void Multiphase_PK::Setup(const Teuchos::Ptr<State>& S)
{
  double t_ini = S->time(); 
  mesh_ = S->GetMesh(domain_);
  dim_ = mesh_->space_dimension();

  // keys
  // -- primary unknowns
  pressure_liquid_key_ = Keys::getKey(domain_, "pressure_liquid"); 
  xl_key_ = Keys::getKey(domain_, "molar_fraction_liquid"); 
  saturation_liquid_key_ = Keys::getKey(domain_, "saturation_liquid"); 

  // -- other fields
  permeability_key_ = Keys::getKey(domain_, "permeability"); 
  relperm_liquid_key_ = Keys::getKey(domain_, "rel_permeability_liquid"); 
  relperm_gas_key_ = Keys::getKey(domain_, "rel_permeability_gas"); 
  porosity_key_ = Keys::getKey(domain_, "porosity"); 

  pressure_gas_key_ = Keys::getKey(domain_, "pressure_gas"); 
  temperature_key_ = Keys::getKey(domain_, "temperature"); 

  // register non-standard fields
  if (!S->HasField("gravity"))
      S->RequireConstantVector("gravity", passwd_, dim_);

  if (!S->HasField("fluid_density"))
      S->RequireScalar("fluid_density", passwd_);

  if (!S->HasField("fluid_viscosity"))
      S->RequireScalar("fluid_viscosity", passwd_);

  if (!S->HasField("gas_viscosity"))
      S->RequireScalar("gas_viscosity", passwd_);

  // pressure as the primary solution
  if (!S->HasField(pressure_liquid_key_)) {
    S->RequireField(pressure_liquid_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", pressure_liquid_key_);
    auto eval = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S->SetFieldEvaluator(pressure_liquid_key_, eval);
  }

  // liquid mole fraction is the primary solution
  component_names_ = mp_list_->template get<Teuchos::Array<std::string> >("primary component names").toVector();
  std::vector<std::vector<std::string> > subfield_names(1);
  subfield_names[0] = component_names_;

  S->RequireField(xl_key_, passwd_, subfield_names)
    ->SetMesh(mesh_)->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::CELL, component_names_.size());

  // water saturation is the primary solution
  if (!S->HasField(saturation_liquid_key_)) {
    S->RequireField(saturation_liquid_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", saturation_liquid_key_);
    saturation_liquid_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S->SetFieldEvaluator(saturation_liquid_key_, saturation_liquid_eval_);
  }

  // pressure gas
  auto wrm_list = Teuchos::sublist(mp_list_, "water retention models", true);
  wrm_ = CreateModelPartition<WRMmp>(mesh_, wrm_list, "water retention model");

  if (!S->HasField(pressure_gas_key_)) {
    S->RequireField(pressure_gas_key_, pressure_gas_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("pressure liquid key", pressure_liquid_key_)
         .set<std::string>("saturation liquid key", saturation_liquid_key_);

    auto eval = Teuchos::rcp(new PressureGasEvaluator(elist, wrm_));
    S->SetFieldEvaluator(pressure_gas_key_, eval);
  }

  // relative permeability of liquid phase
  S->RequireField(relperm_liquid_key_, relperm_liquid_key_)->SetMesh(mesh_)->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1)
    ->AddComponent("face", AmanziMesh::FACE, 1)
    ->AddComponent("dirichlet_faces", AmanziMesh::BOUNDARY_FACE, 1);

  Teuchos::ParameterList elist;
  elist.set<std::string>("my key", relperm_liquid_key_)
       .set<std::string>("saturation liquid key", saturation_liquid_key_)
       .set<std::string>("phase name", "liquid");

  auto eval = Teuchos::rcp(new RelPermEvaluator(elist, wrm_));
  S->SetFieldEvaluator(relperm_liquid_key_, eval);

  // relative permeability of gas phase
  S->RequireField(relperm_gas_key_, relperm_gas_key_)->SetMesh(mesh_)->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1)
    ->AddComponent("face", AmanziMesh::FACE, 1)
    ->AddComponent("dirichlet_faces", AmanziMesh::BOUNDARY_FACE, 1);

  elist.set<std::string>("my key", relperm_gas_key_)
       .set<std::string>("saturation liquid key", saturation_liquid_key_)
       .set<std::string>("phase name", "gas");

  eval = Teuchos::rcp(new RelPermEvaluator(elist, wrm_));
  S->SetFieldEvaluator(relperm_gas_key_, eval);

  // material properties
  if (!S->HasField(permeability_key_)) {
    S->RequireField(permeability_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim_);
  }

  if (!S->HasField(porosity_key_)) {
    S->RequireField(porosity_key_, porosity_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // other fields

  // fields from previous time step
  if (!S->HasField("prev_total_water_storage")) {
    S->RequireField("prev_total_water_storage", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->GetField("prev_total_water_storage", passwd_)->set_io_vis(false);
  }
  if (!S->HasField("prev_total_component_storage")) {
    S->RequireField("prev_total_component_storage", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->GetField("prev_total_component_storage", passwd_)->set_io_vis(false);
  }
}


/* ******************************************************************
* Standard constructor
****************************************************************** */
void Multiphase_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  ncells_owned_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncells_wghost_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  ncp_ = mp_list_->get<std::string>("NCP function", "min");
  smooth_mu_ = mp_list_->get<double>("smoothing parameter mu", 0.0);
  // H_ = mp_list_->get<double>("Henry law constants", 7.65e-6);
  num_aqueous_ = mp_list_->get<int>("number aqueous components");
  num_gaseous_ = mp_list_->get<int>("number gaseous components");
  num_primary_ = component_names_.size();

  // verbose object must go first to support initialization reports
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = mp_list_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("Multiphase", vlist));

  // fundamental physical quantities
  double* gravity_data;
  S->GetConstantVectorData("gravity")->ExtractView(&gravity_data);
  gravity_.set(dim_, &(gravity_data[0]));  // do it in complicated way because we
                                           // are not sure if gravity_data is an
                                           // array or vector
  g_ = fabs(gravity_[dim_ - 1]);

  rho_l_ = *S->GetScalarData("fluid_density");
  mu_l_ = *S->GetScalarData("fluid_viscosity");
  mu_g_ = *S->GetScalarData("gas_viscosity");

  // process CPR list
  cpr_enhanced_ = mp_list_->isSublist("CPR enhancement");
  if (cpr_enhanced_) {
    auto& cpr_list_ = mp_list_->sublist("CPR parameters");
    std::vector<int> correction_blocks = cpr_list_.get<Teuchos::Array<int> >("correction blocks").toVector();
    auto block_names_ = cpr_list_.get<Teuchos::Array<std::string> >("preconditioner").toVector();
    for (int i = 0; i < block_names_.size(); i++)
      AMANZI_ASSERT(pc_list_->isSublist(block_names_[i]));
  }

  // vector of primary variables
  InitSolutionVector();

  // boundary conditions
  Teuchos::RCP<MultiphaseBoundaryFunction> bc;
  auto& bc_list = mp_list_->sublist("boundary conditions");

  bcs_.clear();

  // -- pressure 
  if (bc_list.isSublist("pressure")) {
    PK_DomainFunctionFactory<MultiphaseBoundaryFunction> bc_factory(mesh_);

    Teuchos::ParameterList& tmp_list = bc_list.sublist("pressure");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);
        bc = bc_factory.Create(spec, "boundary pressure", AmanziMesh::FACE, Teuchos::null);
        bc->set_bc_name("pressure");
        bcs_.push_back(bc);
      }
    }
  }

  // -- saturation
  if (bc_list.isSublist("saturation")) {
    PK_DomainFunctionFactory<MultiphaseBoundaryFunction> bc_factory(mesh_);

    Teuchos::ParameterList& tmp_list = bc_list.sublist("saturation");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);
        bc = bc_factory.Create(spec, "boundary saturation", AmanziMesh::FACE, Teuchos::null);
        bc->set_bc_name("saturation");
        bcs_.push_back(bc);
      }
    }
  }

  // allocate memory and populate boundary conditions
  op_bcs_.clear();

  for (int i = 0; i < soln_names_.size(); ++i) {
    auto op_bc = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
    op_bcs_.push_back(op_bc);
  }

  double t_ini = S->time();
  for (int i = 0; i < bcs_.size(); i++) {
    bcs_[i]->Compute(t_ini, t_ini);
    bcs_[i]->ComputeSubmodel(mesh_);
  }

  PopulateBCs();

  // matrix is used to eleviate calculation of residual. It is the mathematical
  // representation of differential operators
  auto& opd_list = mp_list_->sublist("operators").sublist("diffusion operator").sublist("matrix");
  auto& opa_list = mp_list_->sublist("operators").sublist("advection operator").sublist("matrix");

  op_matrix_ = Teuchos::rcp(new Operators::TreeOperator(Teuchos::rcpFromRef(soln_->Map())));
  op_preconditioner_ = Teuchos::rcp(new Operators::TreeOperator(Teuchos::rcpFromRef(soln_->Map())));

  // -- liquid pressure is always go first
  auto opd = Teuchos::rcp(new Operators::PDE_DiffusionFVwithGravity(opd_list, mesh_, 1.0, gravity_));
  pde_matrix_diff_.push_back(opd);
  opd->SetBCs(op_bcs_[0], op_bcs_[0]);
  op_matrix_->SetOperatorBlock(0, 0, opd->global_operator());

  // -- chemical components go second
  for (int i = 0; i < num_primary_; ++i) {
    auto opa = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(opa_list, mesh_));
    opa->SetBCs(op_bcs_[i + 1], op_bcs_[i + 1]);

    opd = Teuchos::rcp(new Operators::PDE_DiffusionFVwithGravity(opd_list, opa->global_operator()));
    pde_matrix_diff_.push_back(opd);
    op_matrix_->SetOperatorBlock(i + 1, i + 1, opa->global_operator());
  }

  // preconditioner is reflects equations of states and is model-specific
  InitPreconditioner();

  for (int i = 0; i < num_primary_ + 1; ++i) {
    auto op = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_)); 
    if (eval_adv_[i] != "") {
      auto tmp = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(opa_list, op->global_operator()));
      tmp->SetBCs(op_bcs_[i], op_bcs_[i]);
    }
    if (eval_diff_[i] != "") {
      auto tmp = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(opd_list, op->global_operator()));
      tmp->SetBCs(op_bcs_[i], op_bcs_[i]);
    }
    op_preconditioner_->SetOperatorBlock(i, i, op->global_operator());
  }

  for (int i = 0; i < num_primary_ + 2; ++i) {
    auto op = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_)); 
    op->SetBCs(op_bcs_[2], op_bcs_[i]);
    op_preconditioner_->SetOperatorBlock(num_primary_, i, op->global_operator());
  }
  
  // initialize time integrator
  std::string ti_method_name = ti_list_->get<std::string>("time integration method", "none");
  AMANZI_ASSERT(ti_method_name == "BDF1");
  Teuchos::ParameterList& bdf1_list = ti_list_->sublist("BDF1");

  if (! bdf1_list.isSublist("verbose object"))
      bdf1_list.sublist("verbose object") = mp_list_->sublist("verbose object");

  bdf1_dae_ = Teuchos::rcp(new BDF1_TI<TreeVector, TreeVectorSpace>(*this, bdf1_list, soln_));

  // linear solver
  pc_name_ = mp_list_->get<std::string>("preconditioner");
  solver_name_ = mp_list_->get<std::string>("linear solver");
  AMANZI_ASSERT(pc_list_->isSublist(pc_name_));

  AmanziSolvers::LinearOperatorFactory<Operators::TreeOperator, TreeVector, TreeVectorSpace> factory;
  op_pc_solver_ = factory.Create(solver_name_, *linear_operator_list_, op_preconditioner_);

  // absolute permeability
  ncells_wghost_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  K_.resize(ncells_wghost_);
  ConvertFieldToTensor(S_, dim_, "permeability", K_);

  // upwind operator with a face model (FIXME)
  auto upw_list = mp_list_->sublist("operators")
                           .sublist("diffusion operator")
                           .sublist("upwind");

  upwind_ = Teuchos::rcp(new Operators::UpwindFlux<int>(mesh_, Teuchos::null));
  upwind_->Init(upw_list);

  // initialize other fields
  InitializeFields_();
}


/* ****************************************************************
* This completes initialization of fields left out by state.
**************************************************************** */
void Multiphase_PK::InitializeFields_()
{
  Teuchos::OSTab tab = vo_->getOSTab();

  InitializeFieldFromField_("prev_total_water_storage", "total_water_storage", true);
  InitializeFieldFromField_("prev_total_component_storage", "total_component_storage", true);
}


/* ****************************************************************
* Auxiliary initialization technique.
**************************************************************** */
void Multiphase_PK::InitializeFieldFromField_(
    const std::string& field0, const std::string& field1, bool call_evaluator)
{
  if (S_->HasField(field0)) {
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


/* ******************************************************************* 
* Performs one time step from time t_old to time t_new either for
* steady-state or transient simulation. If reinit=true, enforce 
* p-lambda constraints.
******************************************************************* */
bool Multiphase_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  dt_ = t_new - t_old;

  // make a copy of primary and conservative fields

  // initialization
  if (num_itrs_ == 0) {
    Teuchos::RCP<TreeVector> udot = Teuchos::rcp(new TreeVector(*soln_));
    udot->PutScalar(0.0);
    bdf1_dae_->SetInitialState(t_old, soln_, udot);

    UpdatePreconditioner(t_old, soln_, dt_);
    num_itrs_++;
  }

  // trying to make a step
  bool failed(false);
  failed = bdf1_dae_->TimeStep(dt_, dt_next_, soln_);
  if (failed) {
    dt_ = dt_next_;
    return failed;
  }

  dt_ = dt_next_;
  num_itrs_++;

  return failed;
}


/* ******************************************************************
* Push data to the state
****************************************************************** */
void Multiphase_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  *S_->GetFieldData(pressure_liquid_key_, passwd_) = *soln_->SubVector(0)->Data();
  *S_->GetFieldData(xl_key_, passwd_) = *soln_->SubVector(1)->Data();
  *S_->GetFieldData(saturation_liquid_key_, passwd_) = *soln_->SubVector(2)->Data();
}

}  // namespace Multiphase
}  // namespace Amanzi

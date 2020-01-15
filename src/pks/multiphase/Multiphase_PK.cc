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
#include "MolarMobilityEvaluator.hh"
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
  num_itrs_(0),
  op_pc_assembled_(false)
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
  x_liquid_key_ = Keys::getKey(domain_, "molar_fraction_liquid"); 
  saturation_liquid_key_ = Keys::getKey(domain_, "saturation_liquid"); 

  // -- other fields
  permeability_key_ = Keys::getKey(domain_, "permeability"); 
  porosity_key_ = Keys::getKey(domain_, "porosity"); 

  relperm_liquid_key_ = Keys::getKey(domain_, "rel_permeability_liquid"); 
  relperm_gas_key_ = Keys::getKey(domain_, "rel_permeability_gas"); 
  molar_mobility_liquid_key_ = Keys::getKey(domain_, "molar_mobility_liquid"); 
  molar_mobility_gas_key_ = Keys::getKey(domain_, "molar_mobility_gas"); 
  viscosity_liquid_key_ = Keys::getKey(domain_, "viscosity_liquid"); 
  viscosity_gas_key_ = Keys::getKey(domain_, "viscosity_gas"); 

  molar_density_liquid_key_ = Keys::getKey(domain_, "molar_density_liquid"); 
  molar_density_gas_key_ = Keys::getKey(domain_, "molar_density_gas"); 

  pressure_gas_key_ = Keys::getKey(domain_, "pressure_gas"); 
  temperature_key_ = Keys::getKey(domain_, "temperature"); 
  darcy_flux_liquid_key_ = Keys::getKey(domain_, "darcy_flux_liquid"); 
  darcy_flux_gas_key_ = Keys::getKey(domain_, "darcy_flux_gas"); 

  tws_key_ = Keys::getKey(domain_, "total_water_storage");
  tcs_key_ = Keys::getKey(domain_, "total_component_storage");

  // register non-standard fields
  if (!S->HasField("gravity"))
      S->RequireConstantVector("gravity", passwd_, dim_);

  if (!S->HasField("const_fluid_density"))
      S->RequireScalar("const_fluid_density", passwd_);

  if (!S->HasField("const_fluid_viscosity"))
      S->RequireScalar("const_fluid_viscosity", passwd_);

  if (!S->HasField("const_gas_viscosity"))
      S->RequireScalar("const_gas_viscosity", passwd_);

  // pressure is the primary solution
  if (!S->HasField(pressure_liquid_key_)) {
    S->RequireField(pressure_liquid_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", pressure_liquid_key_);
    pressure_liquid_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S->SetFieldEvaluator(pressure_liquid_key_, pressure_liquid_eval_);
  }

  // liquid mole fraction is the primary solution
  if (!S->HasField(x_liquid_key_)) {
    component_names_ = mp_list_->sublist("molecular diffusion")
       .get<Teuchos::Array<std::string> >("primary component names").toVector();
    std::vector<std::vector<std::string> > subfield_names(1);
    subfield_names[0] = component_names_;

    S->RequireField(x_liquid_key_, passwd_, subfield_names)
      ->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, component_names_.size());

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", x_liquid_key_);
    x_liquid_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S->SetFieldEvaluator(x_liquid_key_, x_liquid_eval_);
  }

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

  // Darcy velocities
  if (!S->HasField(darcy_flux_liquid_key_)) {
    S->RequireField(darcy_flux_liquid_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);
  }

  if (!S->HasField(darcy_flux_gas_key_)) {
    S->RequireField(darcy_flux_gas_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);
  }

  // relative permeability of liquid phase
  S->RequireField(relperm_liquid_key_, relperm_liquid_key_)->SetMesh(mesh_)->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  Teuchos::ParameterList elist;
  elist.set<std::string>("my key", relperm_liquid_key_)
       .set<std::string>("saturation liquid key", saturation_liquid_key_)
       .set<std::string>("phase name", "liquid");

  auto eval = Teuchos::rcp(new RelPermEvaluator(elist, wrm_));
  S->SetFieldEvaluator(relperm_liquid_key_, eval);

  // relative permeability of gas phase
  S->RequireField(relperm_gas_key_, relperm_gas_key_)->SetMesh(mesh_)->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  elist.set<std::string>("my key", relperm_gas_key_)
       .set<std::string>("saturation liquid key", saturation_liquid_key_)
       .set<std::string>("phase name", "gas");

  eval = Teuchos::rcp(new RelPermEvaluator(elist, wrm_));
  S->SetFieldEvaluator(relperm_gas_key_, eval);

  // molar mobility of liquid phase
  {
    S->RequireField(molar_mobility_liquid_key_, molar_mobility_liquid_key_)
      ->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

    elist.set<std::string>("my key", molar_mobility_liquid_key_)
         .set<std::string>("viscosity key", viscosity_liquid_key_)
         .set<std::string>("molar density key", molar_density_liquid_key_)
         .set<std::string>("phase name", "liquid");

    auto eval = Teuchos::rcp(new MolarMobilityEvaluator(elist, wrm_));
    S->SetFieldEvaluator(molar_mobility_liquid_key_, eval);
  }

  // molar mobility of gas phase
  {
    S->RequireField(molar_mobility_gas_key_, molar_mobility_gas_key_)
      ->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

    elist.set<std::string>("my key", molar_mobility_gas_key_)
         .set<std::string>("viscosity key", viscosity_gas_key_)
         .set<std::string>("molar density key", molar_density_gas_key_)
         .set<std::string>("phase name", "gas");

    auto eval = Teuchos::rcp(new MolarMobilityEvaluator(elist, wrm_));
    S->SetFieldEvaluator(molar_mobility_gas_key_, eval);
  }

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
* Initialize various PK objects
****************************************************************** */
void Multiphase_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  ncells_owned_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncells_wghost_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  nfaces_owned_ = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nfaces_wghost_ = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  ncp_ = mp_list_->get<std::string>("NCP function", "min");
  smooth_mu_ = mp_list_->get<double>("smoothing parameter mu", 0.0);
  num_aqueous_ = mp_list_->get<int>("number aqueous components");
  num_primary_ = component_names_.size();

  auto tmp_list = mp_list_->sublist("molecular diffusion");
  mol_diff_l_ = tmp_list.get<Teuchos::Array<double> >("aqueous values").toVector();
  mol_diff_g_ = tmp_list.get<Teuchos::Array<double> >("gaseous values").toVector();
  kH_ = tmp_list.get<Teuchos::Array<double> >("Henry dimensionless constants").toVector();

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

  rho_l_ = *S->GetScalarData("const_fluid_density");
  eta_l_ = rho_l_ / CommonDefs::MOLAR_MASS_H2O;
  mu_l_ = *S->GetScalarData("const_fluid_viscosity");
  mu_g_ = *S->GetScalarData("const_gas_viscosity");

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
  InitMPSolutionVector();

  // boundary conditions
  Teuchos::RCP<MultiphaseBoundaryFunction> bc;
  auto& bc_list = mp_list_->sublist("boundary conditions");

  bcs_.clear();

  // -- pressure 
  if (bc_list.isSublist("pressure liquid")) {
    PK_DomainFunctionFactory<MultiphaseBoundaryFunction> bc_factory(mesh_);

    Teuchos::ParameterList& tmp_list = bc_list.sublist("pressure liquid");
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

  // -- total injected mass flux
  if (bc_list.isSublist("mass flux total")) {
    PK_DomainFunctionFactory<MultiphaseBoundaryFunction> bc_factory(mesh_);

    Teuchos::ParameterList& tmp_list = bc_list.sublist("mass flux total");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      if (it->second.isList()) {
        Teuchos::ParameterList spec = Teuchos::getValue<Teuchos::ParameterList>(it->second);
        bc = bc_factory.Create(spec, "outward mass flux", AmanziMesh::FACE, Teuchos::null);
        bc->set_bc_name("flux");
        bc->SetComponentId(component_names_);
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

  op_bc_pg_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

  PopulateBCs(0);

  // matrix is used to eleviate calculation of residual. It is the mathematical
  // representation of differential operators
  auto& ddf_list = mp_list_->sublist("operators").sublist("diffusion operator").sublist("matrix");
  auto& mdf_list = mp_list_->sublist("operators").sublist("molecular diffusion operator").sublist("matrix");
  auto& adv_list = mp_list_->sublist("operators").sublist("advection operator").sublist("matrix");

  // -- liquid pressure is always go first
  pde_diff_K_ = Teuchos::rcp(new Operators::PDE_DiffusionFVwithGravity(ddf_list, mesh_, 1.0, gravity_));
  pde_diff_K_->SetBCs(op_bcs_[0], op_bcs_[0]);

  // -- chemical components go second (Darcy law + molecular diffusion)
  pde_diff_D_ = Teuchos::rcp(new Operators::PDE_DiffusionFV(mdf_list, mesh_));
  pde_diff_D_->SetBCs(op_bcs_[1], op_bcs_[1]);

  // preconditioner incorporates equations of states and is model-specific.
  // It is created in the scope of global assembly to reduce memory footprint.
  // Here we just initialize required basis supporting structures.
  op_preconditioner_ = Teuchos::rcp(new Operators::TreeOperator(Teuchos::rcpFromRef(soln_->Map())));
  InitMPPreconditioner();

  std::string pc_name = mp_list_->get<std::string>("preconditioner");
  AMANZI_ASSERT(pc_list_->isSublist(pc_name));
  auto tmp = pc_list_->sublist(pc_name);
  op_preconditioner_->InitializePreconditioner(tmp);
  
  // linear solver
  op_pc_solver_ = op_preconditioner_;
  std::string solver_name = mp_list_->get<std::string>("linear solver");
  AmanziSolvers::LinearOperatorFactory<Operators::TreeOperator, TreeVector, TreeVectorSpace> sfactory;
  op_pc_solver_ = sfactory.Create(solver_name, *linear_operator_list_, op_preconditioner_);

  // initialize time integrator
  std::string ti_method_name = ti_list_->get<std::string>("time integration method", "none");
  AMANZI_ASSERT(ti_method_name == "BDF1");
  Teuchos::ParameterList& bdf1_list = ti_list_->sublist("BDF1");

  if (! bdf1_list.isSublist("verbose object"))
      bdf1_list.sublist("verbose object") = mp_list_->sublist("verbose object");

  bdf1_dae_ = Teuchos::rcp(new BDF1_TI<TreeVector, TreeVectorSpace>(*this, bdf1_list, soln_));

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

  InitializeField(S_.ptr(), passwd_, darcy_flux_liquid_key_, 0.0);
  InitializeField(S_.ptr(), passwd_, darcy_flux_gas_key_, 0.0);
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
  // -- primary field
  CompositeVector pl_copy(*S_->GetFieldData(pressure_liquid_key_));
  CompositeVector sl_copy(*S_->GetFieldData(saturation_liquid_key_));
  CompositeVector xl_copy(*S_->GetFieldData(x_liquid_key_));
  CompositeVector ql_copy(*S_->GetFieldData(darcy_flux_liquid_key_));
  CompositeVector qg_copy(*S_->GetFieldData(darcy_flux_gas_key_));

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

    // recover the original fields
    *S_->GetFieldData(pressure_liquid_key_, passwd_) = pl_copy;
    pressure_liquid_eval_->SetFieldAsChanged(S_.ptr());

    *S_->GetFieldData(saturation_liquid_key_, passwd_) = sl_copy;
    saturation_liquid_eval_->SetFieldAsChanged(S_.ptr());

    *S_->GetFieldData(x_liquid_key_, passwd_) = xl_copy;
    x_liquid_eval_->SetFieldAsChanged(S_.ptr());

    *S_->GetFieldData(darcy_flux_liquid_key_, passwd_) = ql_copy;
    *S_->GetFieldData(darcy_flux_gas_key_, passwd_) = qg_copy;

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
  auto psol = PressureToSolution();
  auto ssol = SaturationToSolution();
  auto csol = ComponentToSolution(num_primary_);
  auto pl = soln_->SubVector(psol.first)->Data();

  *S_->GetFieldData(pressure_liquid_key_, passwd_) = *pl;
  *S_->GetFieldData(saturation_liquid_key_, passwd_) = *soln_->SubVector(ssol.first)->Data();
  *S_->GetFieldData(x_liquid_key_, passwd_) = *soln_->SubVector(csol.first)->Data();

  *S_->GetFieldData("prev_total_component_storage", passwd_) = *S_->GetFieldData(tcs_key_);
  *S_->GetFieldData("prev_total_water_storage", passwd_) = *S_->GetFieldData(tws_key_);

  // update time integration history
  bdf1_dae_->CommitSolution(t_new - t_old, soln_);
  ChangedSolution();

  // update Darcy fluxes
  // -- fluxes
  auto flux_liquid = S_->GetFieldData(darcy_flux_liquid_key_, passwd_);
  auto flux_gas = S_->GetFieldData(darcy_flux_gas_key_, passwd_);

  // -- other fields
  S_->GetFieldEvaluator(molar_mobility_gas_key_)->HasFieldChanged(S_.ptr(), passwd_);
  const auto& mobility_gc = *S_->GetFieldData(molar_mobility_gas_key_)->ViewComponent("cell");

  S_->GetFieldEvaluator(molar_density_gas_key_)->HasFieldChanged(S_.ptr(), passwd_);
  auto eta_g = S_->GetFieldData(molar_density_gas_key_);
  const auto& eta_gc = *eta_g->ViewComponent("cell");

  S_->GetFieldEvaluator(pressure_gas_key_)->HasFieldChanged(S_.ptr(), passwd_);
  auto pg = S_->GetFieldData(pressure_gas_key_);

  // -- work memory
  auto kr = CreateCVforUpwind(mesh_);
  auto& kr_c = *kr->ViewComponent("cell");

 Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);

  // -- liquid (replace rho_l_ with mass density vector FIXME)
  auto& pdeK = pde_diff_K_;
  pdeK->Setup(Kptr, kr, Teuchos::null, rho_l_, gravity_);
  pdeK->global_operator()->Init();
  pdeK->UpdateMatrices(Teuchos::null, Teuchos::null);
  pdeK->UpdateFlux(pl.ptr(), flux_liquid.ptr());
  flux_liquid->Scale(1.0 / eta_l_);

  // -- gas (replace eta_g with mass density vector FIXME)
  for (int c = 0; c < ncells_owned_; ++c) {
    kr_c[0][c] = mobility_gc[0][c] / eta_gc[0][c];
  }
  kr->ViewComponent("dirichlet_faces")->PutScalar(0.0);  // FIXME
  upwind_->Compute(*flux_gas, *kr, op_bcs_[0]->bc_model(), *kr);

  pdeK = pde_diff_K_;
  pdeK->Setup(Kptr, kr, Teuchos::null, eta_g, gravity_);
  pdeK->global_operator()->Init();
  pdeK->UpdateMatrices(Teuchos::null, Teuchos::null);
  pdeK->UpdateFlux(pg.ptr(), flux_gas.ptr());
}

}  // namespace Multiphase
}  // namespace Amanzi

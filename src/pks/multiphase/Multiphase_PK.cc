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

// Multiphase
#include "ModelMeshPartition.hh"
#include "Multiphase_PK.hh"
#include "Multiphase_Utils.hh"
#include "MultiphaseTypeDefs.hh"

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
  pressure_liquid_key_ = Keys::getKey(domain_, "pressure_liquid"); 
  xl_key_ = Keys::getKey(domain_, "mole_fraction_liquid"); 
  saturation_liquid_key_ = Keys::getKey(domain_, "saturation_liquid"); 

  permeability_key_ = Keys::getKey(domain_, "permeability"); 

  // register non-standard fields
  if (!S->HasField("gravity")) {
    S->RequireConstantVector("gravity", passwd_, dim_);
  } 

  // pressure as the primary solution
  S->RequireField(pressure_liquid_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::CELL, 1);

  // liquid mole fraction is the primary solution
  component_names_ = mp_list_->template get<Teuchos::Array<std::string> >("primary component names").toVector();
  std::vector<std::vector<std::string> > subfield_names(1);
  subfield_names[0] = component_names_;

  S->RequireField(xl_key_, passwd_, subfield_names)
    ->SetMesh(mesh_)->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::CELL, component_names_.size());

  // water saturation is the primary solution
  S->RequireField(saturation_liquid_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::CELL, 1);

  // material properties
  if (!S->HasField(permeability_key_)) {
    S->RequireField(permeability_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim_);
  }
}


/* ******************************************************************
* Standard constructor
****************************************************************** */
void Multiphase_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  ncp_ = mp_list_->get<std::string>("NCP function", "min");
  mu_ = mp_list_->get<double>("smoothing parameter mu", 0.0);
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

  // -- hydrogen density
  if (bc_list.isSublist("hydrogen density")) {
    PK_DomainFunctionFactory<MultiphaseBoundaryFunction> bc_factory(mesh_);

    Teuchos::ParameterList& tmp_list = bc_list.sublist("hydrogen density");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);
        bc = bc_factory.Create(spec, "boundary hydrogen density", AmanziMesh::FACE, Teuchos::null);
        bc->set_bc_name("hydrogen density");
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
  opd->SetBCs(op_bcs_[0], op_bcs_[0]);
  op_matrix_->SetOperatorBlock(0, 0, opd->global_operator());

  // -- chemical components go second
  for (int i = 0; i < num_primary_; ++i) {
    auto opa = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(opa_list, mesh_));
    opa->SetBCs(op_bcs_[i + 1], op_bcs_[i + 1]);

    opd = Teuchos::rcp(new Operators::PDE_DiffusionFVwithGravity(opd_list, opa->global_operator()));
    op_matrix_->SetOperatorBlock(i + 1, i + 1, opa->global_operator());
  }

  // preconditioner is reflects equations of states and is model-specific
  InitPreconditioner();

  for (int i = 0; i < num_primary_ + 1; ++i) {
    auto op = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_)); 
    if (eval_adv_[i] != Teuchos::null) {
      auto tmp = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(opa_list, op->global_operator()));
      tmp->SetBCs(op_bcs_[i], op_bcs_[i]);
    }
    if (eval_diff_[i] != Teuchos::null) {
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

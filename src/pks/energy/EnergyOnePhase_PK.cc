/*
  This is the Energy component of the Amanzi code.
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Process kernel for the single-phase energy equation.
*/

// Amanzi
#include "EOSEvaluator.hh"
#include "PDE_DiffusionFactory.hh"

// Amanzi::Energy
#include "EnergyOnePhase_PK.hh"
#include "EnthalpyEvaluator.hh"
#include "IEMEvaluator.hh"
#include "TCMEvaluator_OnePhase.hh"
#include "TotalEnergyEvaluator.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Default constructor.
****************************************************************** */
EnergyOnePhase_PK::EnergyOnePhase_PK(
                   Teuchos::ParameterList& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& glist,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& soln) :
    Energy_PK(glist, S),
    soln_(soln)
{
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  ep_list_ = Teuchos::sublist(pk_list, "energy", true);

  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  ti_list_ = Teuchos::sublist(ep_list_, "time integrator");
   
  // domain name
  domain_ = ep_list_->get<std::string>("domain name", "domain");
};


/* ******************************************************************
* Create the physical evaluators for energy, enthalpy, thermal
* conductivity, and any sources.
****************************************************************** */
void EnergyOnePhase_PK::Setup(const Teuchos::Ptr<State>& S)
{
  // basic class setup
  Energy_PK::Setup(S);

  // Get data and evaluators needed by the PK
  // -- energy, the conserved quantity
  S->RequireField(energy_key_)->SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  Teuchos::ParameterList ee_list = ep_list_->sublist("energy evaluator");
  ee_list.set<std::string>("energy key", energy_key_)
         .set<bool>("vapor diffusion", false);
  Teuchos::RCP<TotalEnergyEvaluator> ee = Teuchos::rcp(new TotalEnergyEvaluator(ee_list));
  S->SetFieldEvaluator(energy_key_, ee);

  // -- advection of enthalpy
  S->RequireField(enthalpy_key_)->SetMesh(mesh_)
    ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);

  if (!S->HasFieldEvaluator(enthalpy_key_)) {
    Teuchos::ParameterList enth_plist = ep_list_->sublist("enthalpy evaluator");
    enth_plist.set("enthalpy key", enthalpy_key_);
    Teuchos::RCP<EnthalpyEvaluator> enth = Teuchos::rcp(new EnthalpyEvaluator(enth_plist));
    S->SetFieldEvaluator(enthalpy_key_, enth);
  }

  // -- thermal conductivity
  S->RequireField(conductivity_key_)->SetMesh(mesh_)
    ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList tcm_plist = ep_list_->sublist("thermal conductivity evaluator");
  Teuchos::RCP<TCMEvaluator_OnePhase> tcm = Teuchos::rcp(new TCMEvaluator_OnePhase(tcm_plist));
  S->SetFieldEvaluator(conductivity_key_, tcm);
}


/* ******************************************************************
* Initialize the needed models to plug in enthalpy.
****************************************************************** */
void EnergyOnePhase_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  // create verbosity object
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = ep_list_->sublist("verbose object");
  vo_ =  Teuchos::rcp(new VerboseObject("EnergyPK::2Phase", vlist)); 

  // Call the base class initialize.
  Energy_PK::Initialize(S);

  // Create pointers to the primary flow field pressure.
  solution = S->GetFieldData(temperature_key_, passwd_);
  soln_->SetData(solution); 

  // Create local evaluators. Initialize local fields.
  InitializeFields_();

  // initialize independent operators: diffusion and advection 
  Teuchos::ParameterList tmp_list = ep_list_->sublist("operators").sublist("diffusion operator");
  Teuchos::ParameterList oplist_matrix = tmp_list.sublist("matrix");
  Teuchos::ParameterList oplist_pc = tmp_list.sublist("preconditioner");

  Operators::PDE_DiffusionFactory opfactory;
  op_matrix_diff_ = opfactory.Create(oplist_matrix, mesh_, op_bc_);
  op_matrix_diff_->SetBCs(op_bc_, op_bc_);
  op_matrix_ = op_matrix_diff_->global_operator();
  op_matrix_->Init();
  op_matrix_diff_->SetScalarCoefficient(S->GetFieldData(conductivity_key_), Teuchos::null);

  Teuchos::ParameterList oplist_adv = ep_list_->sublist("operators").sublist("advection operator");
  op_matrix_advection_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(oplist_adv, mesh_));

  const CompositeVector& flux = *S->GetFieldData("darcy_flux");
  op_matrix_advection_->Setup(flux);
  op_advection_ = op_matrix_advection_->global_operator();

  // initialize coupled operators: diffusion + advection + accumulation
  op_preconditioner_diff_ = opfactory.Create(oplist_pc, mesh_, op_bc_);
  op_preconditioner_diff_->SetBCs(op_bc_, op_bc_);
  op_preconditioner_ = op_preconditioner_diff_->global_operator();
  op_preconditioner_->Init();
  op_preconditioner_diff_->SetScalarCoefficient(S->GetFieldData(conductivity_key_), Teuchos::null);

  op_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op_preconditioner_));
  op_preconditioner_advection_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(oplist_adv, op_preconditioner_));
  op_preconditioner_->SymbolicAssembleMatrix();

  // initialize preconditioner
  AMANZI_ASSERT(ti_list_->isParameter("preconditioner"));
  std::string name = ti_list_->get<std::string>("preconditioner");
  Teuchos::ParameterList slist = preconditioner_list_->sublist(name);
  op_preconditioner_->InitializePreconditioner(slist);

  // initialize time integrator
  std::string ti_method_name = ti_list_->get<std::string>("time integration method", "none");
  if (ti_method_name == "BDF1") {
    Teuchos::ParameterList& bdf1_list = ti_list_->sublist("BDF1");

    if (! bdf1_list.isSublist("verbose object"))
        bdf1_list.sublist("verbose object") = ep_list_->sublist("verbose object");

    bdf1_dae_ = Teuchos::rcp(new BDF1_TI<TreeVector, TreeVectorSpace>(*this, bdf1_list, soln_));
  }

  // output of initialization header
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << std::endl << vo_->color("green")
               << "Initialization of TI period is complete." << vo_->reset() << std::endl;
  }
}


/* ****************************************************************
* This completes initialization of missed fields in the state.
**************************************************************** */
void EnergyOnePhase_PK::InitializeFields_()
{
  Teuchos::OSTab tab = vo_->getOSTab();

  if (S_->HasField(prev_energy_key_)) {
    if (!S_->GetField(prev_energy_key_, passwd_)->initialized()) {
      temperature_eval_->SetFieldAsChanged(S_.ptr());
      S_->GetFieldEvaluator(energy_key_)->HasFieldChanged(S_.ptr(), passwd_);

      const CompositeVector& e1 = *S_->GetFieldData(energy_key_);
      CompositeVector& e0 = *S_->GetFieldData(prev_energy_key_, passwd_);
      e0 = e1;

      S_->GetField(prev_energy_key_, passwd_)->set_initialized();

      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
          *vo_->os() << "initialized prev_energy to previous energy" << std::endl;  
    }
  }
}


/* ******************************************************************* 
* Performs one time step of size dt_ either for steady-state or 
* transient sumulation.
******************************************************************* */
bool EnergyOnePhase_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  dt_ = t_new - t_old;

  // save a copy of pressure
  CompositeVector temperature_copy(*S_->GetFieldData(temperature_key_, passwd_));

  // swap conserved field (i.e., energy) and save
  S_->GetFieldEvaluator(energy_key_)->HasFieldChanged(S_.ptr(), passwd_);
  const CompositeVector& e = *S_->GetFieldData(energy_key_);
  CompositeVector& e_prev = *S_->GetFieldData(prev_energy_key_, passwd_);

  CompositeVector e_prev_copy(e_prev);
  e_prev = e;

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

    // restore the original primary solution, temperature
    *S_->GetFieldData(temperature_key_, passwd_) = temperature_copy;
    temperature_eval_->SetFieldAsChanged(S_.ptr());

    // restore the original fields
    *S_->GetFieldData(prev_energy_key_, passwd_) = e_prev_copy;

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Step failed. Restored temperature, prev_energy." << std::endl;

    return failed;
  }

  // commit solution (should we do it here ?)
  bdf1_dae_->CommitSolution(dt_, soln_);
  temperature_eval_->SetFieldAsChanged(S_.ptr());

  num_itrs_++;
  dt_ = dt_next_;
  
  return failed;
}


/* ******************************************************************
* TBW 
****************************************************************** */
void EnergyOnePhase_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  dt_ = dt_next_;
}

}  // namespace Energy
}  // namespace Amanzi

/*
  This is the Energy component of the Amanzi code.
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Process kernel for thermal Richards' flow.
*/

#include "OperatorDiffusionFactory.hh"

#include "EnergyTwoPhase_PK.hh"
#include "eos_evaluator.hh"
#include "EnthalpyEvaluator.hh"
#include "iem_evaluator.hh"
#include "TwoPhaseEnergyEvaluator.hh"
#include "twophase_thermal_conductivity_evaluator.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Default constructor for Thermal Richrads PK.
****************************************************************** */
EnergyTwoPhase_PK::EnergyTwoPhase_PK(
                   Teuchos::ParameterList& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& glist,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& soln) :
    Energy_PK(glist, S),
    soln_(soln)
    
{
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  ep_list_ = Teuchos::sublist(pk_list, "Energy", true);

  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "Preconditioners", true);
  ti_list_ = Teuchos::sublist(ep_list_, "time integrator");
};


/* ******************************************************************
* Create the physical evaluators for energy, enthalpy, thermal
* conductivity, and any sources.
****************************************************************** */
void EnergyTwoPhase_PK::Setup()
{
  // basic class setup
  Energy_PK::Setup();

  // Get data and evaluators needed by the PK
  // -- energy, the conserved quantity
  S_->RequireField(energy_key_)->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  Teuchos::ParameterList ee_list = glist_->sublist("PKs").sublist("Energy").sublist("energy evaluator");
  ee_list.set("energy key", energy_key_);
  Teuchos::RCP<TwoPhaseEnergyEvaluator> ee = Teuchos::rcp(new TwoPhaseEnergyEvaluator(ee_list));
  S_->SetFieldEvaluator(energy_key_, ee);

  // advection of enthalpy
  S_->RequireField(enthalpy_key_)->SetMesh(mesh_)
    ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);

  Teuchos::ParameterList enth_plist = glist_->sublist("PKs").sublist("Energy").sublist("enthalpy evaluator");
  enth_plist.set("enthalpy key", enthalpy_key_);
  Teuchos::RCP<EnthalpyEvaluator> enth = Teuchos::rcp(new EnthalpyEvaluator(enth_plist));
  S_->SetFieldEvaluator(enthalpy_key_, enth);

  // thermal conductivity
  S_->RequireField(conductivity_key_)->SetMesh(mesh_)
    ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList tcm_plist = glist_->sublist("PKs")
                                            .sublist("Energy")
                                            .sublist("thermal conductivity evaluator");
  Teuchos::RCP<ThermalConductivityTwoPhaseEvaluator> tcm =
      Teuchos::rcp(new ThermalConductivityTwoPhaseEvaluator(tcm_plist));
  S_->SetFieldEvaluator(conductivity_key_, tcm);
}


/* ******************************************************************
* Initialize the needed models to plug in enthalpy.
****************************************************************** */
void EnergyTwoPhase_PK::Initialize()
{
  // create verbosity object
  Teuchos::ParameterList vlist;
  vlist.sublist("VerboseObject") = ep_list_->sublist("VerboseObject");
  vo_ = new VerboseObject("EnergyPK::2Phase", vlist); 

  // Create a scalar tensor so far
  K.resize(ncells_owned);
  for (int c = 0; c < ncells_owned; c++) {
    K[c].Init(dim, 1);
    K[c](0, 0) = 1.0;
  }

  // Call the base class initialize.
  Energy_PK::Initialize();

  // Create pointers to the primary flow field pressure.
  solution = S_->GetFieldData("temperature", passwd_);
  soln_->SetData(solution); 

  // create evaluators
  Teuchos::RCP<FieldEvaluator> eos_fe = S_->GetFieldEvaluator("molar_density_liquid");
  Teuchos::RCP<Relations::EOSEvaluator> eos_eval = Teuchos::rcp_dynamic_cast<Relations::EOSEvaluator>(eos_fe);
  ASSERT(eos_eval != Teuchos::null);
  eos_liquid_ = eos_eval->get_EOS();

  Teuchos::RCP<FieldEvaluator> iem_fe = S_->GetFieldEvaluator("internal_energy_liquid");
  Teuchos::RCP<IEMEvaluator> iem_eval = Teuchos::rcp_dynamic_cast<IEMEvaluator>(iem_fe);
  ASSERT(iem_eval != Teuchos::null);
  iem_liquid_ = iem_eval->get_IEM();

  // Create local evaluators. Initialize local fields.
  InitializeFields_();

  // initialize matrix and preconditioner operators
  Teuchos::ParameterList tmp_list = ep_list_->sublist("operators").sublist("diffusion operator");
  Teuchos::ParameterList oplist_matrix = tmp_list.sublist("matrix");
  Teuchos::ParameterList oplist_pc = tmp_list.sublist("preconditioner");

  AmanziGeometry::Point g(dim);

  Operators::OperatorDiffusionFactory opfactory;
  op_matrix_diff_ = opfactory.Create(mesh_, op_bc_, oplist_matrix, g, 0);
  op_matrix_diff_->SetBCs(op_bc_);
  op_matrix_ = op_matrix_diff_->global_operator();
  op_matrix_->Init();
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K);
  op_matrix_diff_->Setup(Kptr, Teuchos::null, Teuchos::null, 1.0, 1.0);

  Teuchos::ParameterList oplist_adv = ep_list_->sublist("operators").sublist("advection operator");
  op_advection_ = Teuchos::rcp(new Operators::OperatorAdvection(oplist_adv, op_matrix_));

  const CompositeVector& flux = *S_->GetFieldData("darcy_flux");
  op_advection_->Setup(flux);

  op_preconditioner_diff_ = opfactory.Create(mesh_, op_bc_, oplist_pc, g, 0);
  op_preconditioner_diff_->SetBCs(op_bc_);
  op_preconditioner_ = op_preconditioner_diff_->global_operator();
  op_preconditioner_->Init();
  op_preconditioner_diff_->Setup(Kptr, Teuchos::null, Teuchos::null, 1.0, 1.0);
  op_preconditioner_->SymbolicAssembleMatrix();

  op_acc_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, op_preconditioner_));

  // preconditioner and optional linear solver
  ASSERT(ti_list_->isParameter("preconditioner"));
  preconditioner_name_ = ti_list_->get<std::string>("preconditioner");

  // output of initializa header
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << std::endl << vo_->color("green")
               << "Initalization of TI period is complete." << vo_->reset() << std::endl;
  }
}


/* ******************************************************************
 * * Transfer part of the internal data needed by transport to the 
 * * flow state FS_MPC. MPC may request to populate the original FS.
 * * The consistency condition is improved by adjusting saturation while
 * * preserving its LED property.
 * ****************************************************************** */
void EnergyTwoPhase_PK::CommitStep(double t_old, double t_new)
{
}

}  // namespace Energy
}  // namespace Amanzi

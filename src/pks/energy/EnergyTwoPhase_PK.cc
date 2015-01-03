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

#include "EnergyTwoPhase_PK.hh"
#include "eos_evaluator.hh"
#include "enthalpy_evaluator.hh"
#include "iem_evaluator.hh"
#include "twophase_energy_evaluator.hh"
#include "twophase_thermal_conductivity_evaluator.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Default constructor for Thermal Richrads PK.
****************************************************************** */
EnergyTwoPhase_PK::EnergyTwoPhase_PK(
    Teuchos::RCP<const Teuchos::ParameterList>& glist, Teuchos::RCP<State>& S)
    : Energy_PK(glist, S) {};


/* ******************************************************************
* Create the physical evaluators for energy, enthalpy, thermal
* conductivity, and any sources.
****************************************************************** */
void EnergyTwoPhase_PK::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S)
{
  // Get data and evaluators needed by the PK
  // -- energy, the conserved quantity
  S->RequireField(energy_key_)->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  Teuchos::ParameterList ee_list = glist_->sublist("Energy").sublist("energy evaluator");
  ee_list.set("energy key", energy_key_);
  Teuchos::RCP<TwoPhaseEnergyEvaluator> ee = Teuchos::rcp(new TwoPhaseEnergyEvaluator(ee_list));
  S->SetFieldEvaluator(energy_key_, ee);

  // -- advection of enthalpy
  S->RequireField(enthalpy_key_)->SetMesh(mesh_)
      ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);

  Teuchos::ParameterList enth_plist = glist_->sublist("energy").sublist("enthalpy evaluator");
  enth_plist.set("enthalpy key", enthalpy_key_);
  Teuchos::RCP<EnthalpyEvaluator> enth = Teuchos::rcp(new EnthalpyEvaluator(enth_plist));
  S->SetFieldEvaluator(enthalpy_key_, enth);

  // -- thermal conductivity
  S->RequireField(conductivity_key_)->SetMesh(mesh_)
      ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList tcm_plist = glist_->sublist("Energy").sublist("thermal conductivity evaluator");
  Teuchos::RCP<ThermalConductivityTwoPhaseEvaluator> tcm =
      Teuchos::rcp(new ThermalConductivityTwoPhaseEvaluator(tcm_plist));
  S->SetFieldEvaluator(conductivity_key_, tcm);
}


/* ******************************************************************
* Initialize the needed models to plug in enthalpy.
****************************************************************** */
void EnergyTwoPhase_PK::Initialize()
{
  // Call the base class's initialize.
  Energy_PK::Initialize();

  // For the boundary conditions, we currently hack in the enthalpy to
  // the boundary faces to correctly advect in a Dirichlet temperature BC.
  // This requires density and internal energy, which in turn require 
  // a model based on p,T. This will be removed once boundary faces 
  // are implemented.
  Teuchos::RCP<FieldEvaluator> eos_fe = S_->GetFieldEvaluator("molar_density_liquid");
  Teuchos::RCP<Relations::EOSEvaluator> eos_eval = Teuchos::rcp_dynamic_cast<Relations::EOSEvaluator>(eos_fe);
  ASSERT(eos_eval != Teuchos::null);
  eos_liquid_ = eos_eval->get_EOS();

  Teuchos::RCP<FieldEvaluator> iem_fe = S_->GetFieldEvaluator("internal_energy_liquid");
  Teuchos::RCP<IEMEvaluator> iem_eval = Teuchos::rcp_dynamic_cast<IEMEvaluator>(iem_fe);
  ASSERT(iem_eval != Teuchos::null);
  iem_liquid_ = iem_eval->get_IEM();
}

}  // namespace Energy
}  // namespace Amanzi

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
    Teuchos::ParameterList& pk_tree,
    const Teuchos::RCP<Teuchos::ParameterList>& glist,
    const Teuchos::RCP<State>& S,
    const Teuchos::RCP<TreeVector>& soln) : Energy_PK(glist, S) 
{
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  ep_list_ = Teuchos::sublist(pk_list, "Energy", true);
};


/* ******************************************************************
* Create the physical evaluators for energy, enthalpy, thermal
* conductivity, and any sources.
****************************************************************** */
void EnergyTwoPhase_PK::Setup()
{
  Energy_PK::Setup();

  // Get data and evaluators needed by the PK
  // -- energy, the conserved quantity
  S_->RequireField(energy_key_)->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  Teuchos::ParameterList ee_list = glist_->sublist("PKs").sublist("Energy").sublist("energy evaluator");
  ee_list.set("energy key", energy_key_);
  Teuchos::RCP<TwoPhaseEnergyEvaluator> ee = Teuchos::rcp(new TwoPhaseEnergyEvaluator(ee_list));
  S_->SetFieldEvaluator(energy_key_, ee);

  // -- advection of enthalpy
  S_->RequireField(enthalpy_key_)->SetMesh(mesh_)
      ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);

  Teuchos::ParameterList enth_plist = glist_->sublist("PKs").sublist("Energy").sublist("enthalpy evaluator");
  enth_plist.set("enthalpy key", enthalpy_key_);
  Teuchos::RCP<EnthalpyEvaluator> enth = Teuchos::rcp(new EnthalpyEvaluator(enth_plist));
  S_->SetFieldEvaluator(enthalpy_key_, enth);

  // -- thermal conductivity
  S_->RequireField(conductivity_key_)->SetMesh(mesh_)
      ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList tcm_plist = glist_->sublist("PKs").sublist("Energy").sublist("thermal conductivity evaluator");
  Teuchos::RCP<ThermalConductivityTwoPhaseEvaluator> tcm =
      Teuchos::rcp(new ThermalConductivityTwoPhaseEvaluator(tcm_plist));
  S_->SetFieldEvaluator(conductivity_key_, tcm);
}


/* ******************************************************************
* Initialize the needed models to plug in enthalpy.
****************************************************************** */
void EnergyTwoPhase_PK::Initialize()
{
  // Call the base class's initialize.
  Energy_PK::Initialize();

  // create verbosity object
  Teuchos::ParameterList vlist;
  vlist.sublist("VerboseObject") = ep_list_->sublist("VerboseObject");
  vo_ = new VerboseObject("EnergyPK::2Phase", vlist); 

  // create evaluators
  Teuchos::RCP<FieldEvaluator> eos_fe = S_->GetFieldEvaluator("molar_density_liquid");
  Teuchos::RCP<Relations::EOSEvaluator> eos_eval = Teuchos::rcp_dynamic_cast<Relations::EOSEvaluator>(eos_fe);
  ASSERT(eos_eval != Teuchos::null);
  eos_liquid_ = eos_eval->get_EOS();

  Teuchos::RCP<FieldEvaluator> iem_fe = S_->GetFieldEvaluator("internal_energy_liquid");
  Teuchos::RCP<IEMEvaluator> iem_eval = Teuchos::rcp_dynamic_cast<IEMEvaluator>(iem_fe);
  ASSERT(iem_eval != Teuchos::null);
  iem_liquid_ = iem_eval->get_IEM();

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << std::endl 
        << vo_->color("green") << "Initalization of TI period is complete." << vo_->reset() << std::endl;
  }
}

}  // namespace Energy
}  // namespace Amanzi

/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Variables: pressure liquid, mole gas fraction, saturation liquid.
*/

// TPLs
#include "Teuchos_RCP.hpp"

// Amanzi::Energy
#include "TCMEvaluator_TwoPhase.hh"
#include "EnthalpyEvaluator.hh"

// Multiphase
#include "MoleFractionLiquid.hh"
#include "MultiphaseModel1_PK.hh"
#include "NCP_F.hh"
#include "NCP_MoleFractions.hh"
#include "ProductEvaluator.hh"
#include "SaturationGasEvaluator.hh"
#include "TotalComponentStorage.hh"
#include "TotalEnergyEvaluator.hh"
#include "TotalWaterStorage.hh"
#include "VaporPressureEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* Standard constructor
****************************************************************** */
MultiphaseModel1_PK::MultiphaseModel1_PK(Teuchos::ParameterList& pk_tree,
                                         const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                         const Teuchos::RCP<State>& S,
                                         const Teuchos::RCP<TreeVector>& soln)
  : Multiphase_PK(pk_tree, glist, S, soln) {};

  
/* ******************************************************************
* Setup
****************************************************************** */
void MultiphaseModel1_PK::Setup()
{
  Multiphase_PK::Setup();

  advection_water_key_ = Keys::getKey(domain_, "advection_water"); 
  pressure_vapor_key_ = Keys::getKey(domain_, "pressure_vapor"); 
  x_vapor_key_ = Keys::getKey(domain_, "mole_fraction_vapor"); 

  molecular_diff_liquid_key_ = Keys::getKey(domain_, "molecular_diff_liquid"); 
  molecular_diff_gas_key_ = Keys::getKey(domain_, "molecular_diff_gas"); 

  diffusion_liquid_key_ = Keys::getKey(domain_, "diffusion_liquid"); 
  diffusion_vapor_key_ = Keys::getKey(domain_, "diffusion_vapor"); 
  diffusion_gas_key_ = Keys::getKey(domain_, "diffusion_gas"); 

  ie_rock_key_ = Keys::getKey(domain_, "internal_energy_rock");
  ie_liquid_key_ = Keys::getKey(domain_, "internal_energy_liquid");
  ie_gas_key_ = Keys::getKey(domain_, "internal_energy_gas");
  conductivity_key_ = Keys::getKey(domain_, "thermal_conductivity");
  particle_density_key_ = Keys::getKey(domain_, "particle_density");
  enthalpy_liquid_key_ = Keys::getKey(domain_, "enthalpy_liquid");
  enthalpy_gas_key_ = Keys::getKey(domain_, "enthalpy_gas");

  advection_enthalpy_liquid_key_ = Keys::getKey(domain_, "advection_enthalpy_liquid");
  advection_enthalpy_gas_key_ = Keys::getKey(domain_, "advection_enthalpy_gas");

  // gas mole fraction is the primary solution
  if (!S_->HasRecord(x_gas_key_)) {
    component_names_ = mp_list_->sublist("molecular diffusion")
       .get<Teuchos::Array<std::string> >("aqueous names").toVector();

    S_->Require<CV_t, CVS_t>(x_gas_key_, Tags::DEFAULT, passwd_, component_names_)
      .SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, component_names_.size());

    Teuchos::ParameterList elist(x_gas_key_);
    elist.set<std::string>("tag", "");
    x_gas_eval_ = Teuchos::rcp(new EvaluatorPrimary<CV_t, CVS_t>(elist));
    S_->SetEvaluator(x_gas_key_, Tags::DEFAULT, x_gas_eval_);

    S_->RequireDerivative<CV_t, CVS_t>(x_gas_key_, Tags::DEFAULT,
                                       x_gas_key_, Tags::DEFAULT, x_gas_key_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // conserved quantities
  // -- total water storage
  if (!S_->HasRecord(tws_key_)) {
    auto elist = MyRequire_(tws_key_, tws_key_);
    elist.set<std::string>("molar density liquid key", mol_density_liquid_key_)
         .set<std::string>("molar density gas key", mol_density_gas_key_)
         .set<std::string>("porosity key", porosity_key_)
         .set<std::string>("saturation liquid key", saturation_liquid_key_)
         .set<std::string>("mole fraction vapor key", x_vapor_key_);

    eval_tws_ = Teuchos::rcp(new TotalWaterStorage(elist));
    S_->SetEvaluator(tws_key_, Tags::DEFAULT, eval_tws_);

    S_->RequireDerivative<CV_t, CVS_t>(tws_key_, Tags::DEFAULT,
                                       pressure_liquid_key_, Tags::DEFAULT, tws_key_);
    S_->RequireDerivative<CV_t, CVS_t>(tws_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, tws_key_);
  }

  // -- total component storage (for one component)
  if (!S_->HasRecord(tcs_key_)) {
    auto elist = MyRequire_(tcs_key_, tcs_key_);
    elist.set<std::string>("saturation liquid key", saturation_liquid_key_)
         .set<std::string>("porosity key", porosity_key_)
         .set<std::string>("molar density liquid key", mol_density_liquid_key_)
         .set<std::string>("molar density gas key", mol_density_gas_key_)
         .set<std::string>("mole fraction liquid key", x_liquid_key_)
         .set<std::string>("mole fraction gas key", x_gas_key_);

    eval_tcs_ = Teuchos::rcp(new TotalComponentStorage(elist));
    S_->SetEvaluator(tcs_key_, Tags::DEFAULT, eval_tcs_);

    S_->RequireDerivative<CV_t, CVS_t>(tcs_key_, Tags::DEFAULT,
                                       pressure_liquid_key_, Tags::DEFAULT, tcs_key_);
    S_->RequireDerivative<CV_t, CVS_t>(tcs_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, tcs_key_);
    S_->RequireDerivative<CV_t, CVS_t>(tcs_key_, Tags::DEFAULT,
                                       x_gas_key_, Tags::DEFAULT, tcs_key_);
  }

  // -- total energy
  non_isothermal_ = mp_list_->sublist("system").isParameter("energy eqn");

  if (non_isothermal_) {
    if (!S_->HasRecord(energy_key_)) {
      S_->Require<CV_t, CVS_t>(energy_key_, Tags::DEFAULT, energy_key_)
        .SetMesh(mesh_)->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);

      Teuchos::ParameterList elist = mp_list_->sublist("energy evaluator");
      elist.set<std::string>("energy key", energy_key_)
           .set<std::string>("tag", "")
           .set<bool>("vapor diffusion", true)
           .set<std::string>("particle density key", particle_density_key_)
           .set<std::string>("internal energy rock key", ie_rock_key_);
      elist.setName(energy_key_);

      auto ee = Teuchos::rcp(new Energy::TotalEnergyEvaluator(elist));
      S_->SetEvaluator(energy_key_, Tags::DEFAULT, ee);

      S_->RequireDerivative<CV_t, CVS_t>(energy_key_, Tags::DEFAULT,
                                        temperature_key_, Tags::DEFAULT, energy_key_);
    }

    // previous conserved quantities
    if (!S_->HasRecord(prev_energy_key_)) {
      S_->Require<CV_t, CVS_t>(prev_energy_key_, Tags::DEFAULT, passwd_)
        .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
      S_->GetRecordW(prev_energy_key_, passwd_).set_io_vis(false);
    }

    // internal energies 
    MyRequire_(ie_liquid_key_, ie_liquid_key_);
    S_->RequireEvaluator(ie_liquid_key_, Tags::DEFAULT);

    S_->RequireDerivative<CV_t, CVS_t>(ie_liquid_key_, Tags::DEFAULT,
                                       temperature_key_, Tags::DEFAULT, ie_liquid_key_);

    MyRequire_(ie_gas_key_, ie_gas_key_);
    S_->RequireEvaluator(ie_gas_key_, Tags::DEFAULT);

    S_->RequireDerivative<CV_t, CVS_t>(ie_gas_key_, Tags::DEFAULT,
                                       temperature_key_, Tags::DEFAULT, ie_gas_key_);

    MyRequire_(ie_rock_key_, ie_rock_key_);
    S_->RequireEvaluator(ie_rock_key_, Tags::DEFAULT);

    S_->RequireDerivative<CV_t, CVS_t>(ie_rock_key_, Tags::DEFAULT,
                                       temperature_key_, Tags::DEFAULT, ie_rock_key_);

    // -- thermal conductivity
    if (!S_->HasRecord(conductivity_key_)) {
      S_->Require<CV_t, CVS_t>(conductivity_key_, Tags::DEFAULT, conductivity_key_)
        .SetMesh(mesh_)->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);

      Teuchos::ParameterList elist = mp_list_->sublist("thermal conductivity evaluator");
      elist.set("thermal conductivity key", conductivity_key_)
           .set<std::string>("tag", "");
      elist.setName(conductivity_key_);

      auto tcm = Teuchos::rcp(new Energy::TCMEvaluator_TwoPhase(elist));
      S_->SetEvaluator(conductivity_key_, Tags::DEFAULT, tcm);
    }

    // -- advection of enthalpy
    if (!S_->HasRecord(enthalpy_liquid_key_)) {
      S_->Require<CV_t, CVS_t>(enthalpy_liquid_key_, Tags::DEFAULT, enthalpy_liquid_key_)
        .SetMesh(mesh_)->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);

      Teuchos::ParameterList elist = mp_list_->sublist("enthalpy evaluator");
      elist.set("enthalpy key", enthalpy_liquid_key_)
           .set<std::string>("tag", "")
           .set<std::string>("internal energy key", ie_liquid_key_)
           .set<bool>("include work term", false)
           .set<std::string>("pressure key", pressure_liquid_key_)
           .set<std::string>("molar density key", mol_density_liquid_key_);
      elist.setName(enthalpy_liquid_key_);

      auto enth = Teuchos::rcp(new Energy::EnthalpyEvaluator(elist));
      S_->SetEvaluator(enthalpy_liquid_key_, Tags::DEFAULT, enth);

      S_->RequireDerivative<CV_t, CVS_t>(enthalpy_liquid_key_, Tags::DEFAULT,
                                         temperature_key_, Tags::DEFAULT, enthalpy_liquid_key_);
    }

    // -- coefficient for enthalpy advection operator in liquid phase (4 fields)
    if (!S_->HasRecord(advection_enthalpy_liquid_key_)) {
      auto elist = MyRequire_(advection_enthalpy_liquid_key_, advection_enthalpy_liquid_key_);

      Teuchos::Array<int> dep_powers(4, 1);
      Teuchos::Array<std::string> dep_names;

      dep_names.push_back(mass_density_liquid_key_);
      dep_names.push_back(enthalpy_liquid_key_);
      dep_names.push_back(relperm_liquid_key_);
      dep_names.push_back(viscosity_liquid_key_);
      dep_powers[3] = -1;

      elist.set<Teuchos::Array<std::string> >("dependencies", dep_names) 
           .set<bool>("dependency tags are my tag", true)
           .set<Teuchos::Array<int> >("powers", dep_powers);

      auto eval = Teuchos::rcp(new ProductEvaluator(elist));
      S_->SetEvaluator(advection_enthalpy_liquid_key_, Tags::DEFAULT, eval);
    }
  }

  // liquid (water) densities
  if (!S_->HasRecord(mol_density_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(mol_density_liquid_key_, Tags::DEFAULT, mol_density_liquid_key_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireEvaluator(mol_density_liquid_key_, Tags::DEFAULT);
  }

  if (!S_->HasRecord(mass_density_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(mass_density_liquid_key_, Tags::DEFAULT, mass_density_liquid_key_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireEvaluator(mass_density_liquid_key_, Tags::DEFAULT);
  }

  /*
  if (!S_->HasRecord(mass_density_gas_key_)) {
    S_->Require<CV_t, CVS_t>(mass_density_gas_key_, Tags::DEFAULT, mass_density_gas_key_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  */

  // saturation
  if (!S_->HasRecord(saturation_gas_key_)) {
    auto elist = MyRequire_(saturation_gas_key_, saturation_gas_key_);
    elist.set<std::string>("saturation liquid key", saturation_liquid_key_);
    auto eval = Teuchos::rcp(new SaturationGasEvaluator(elist));
    S_->SetEvaluator(saturation_gas_key_, Tags::DEFAULT, eval);
  }

  // water vapor
  // -- vapor pressure
  if (!S_->HasRecord(pressure_vapor_key_)) {
    auto elist = MyRequire_(pressure_vapor_key_, pressure_vapor_key_);
    elist.set<std::string>("temperature key", temperature_key_)
         .set<std::string>("molar density liquid key", mol_density_liquid_key_)
         .set<std::string>("saturation liquid key", saturation_liquid_key_)
         .set<std::string>("eos type", "water vapor over water/ice");
    auto eval = Teuchos::rcp(new VaporPressureEvaluator(elist, wrm_));
    S_->SetEvaluator(pressure_vapor_key_, Tags::DEFAULT, eval);
  }

  // -- coefficient for water vapor diffusion operator in liquid phase
  if (!S_->HasRecord(diffusion_vapor_key_)) {
    auto elist = MyRequire_(diffusion_vapor_key_, diffusion_vapor_key_);

    Teuchos::Array<int> dep_powers(4, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(mol_density_gas_key_);
    dep_names.push_back(molecular_diff_gas_key_);
    dep_names.push_back(porosity_key_);
    dep_names.push_back(saturation_gas_key_);

    elist.set<Teuchos::Array<std::string> >("dependencies", dep_names) 
         .set<bool>("dependency tags are my tag", true)
         .set<Teuchos::Array<int> >("powers", dep_powers);

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S_->SetEvaluator(diffusion_vapor_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(diffusion_vapor_key_, Tags::DEFAULT,
                                       pressure_liquid_key_, Tags::DEFAULT, diffusion_vapor_key_);
    S_->RequireDerivative<CV_t, CVS_t>(diffusion_vapor_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, diffusion_vapor_key_);
  }

  // liquid mole fraction (for current component)
  if (!S_->HasRecord(x_liquid_key_)) {
    auto elist = MyRequire_(x_liquid_key_, x_liquid_key_);
    elist.set<std::string>("pressure gas key", pressure_gas_key_)
         .set<std::string>("mole fraction gas key", x_gas_key_);
    auto eval = Teuchos::rcp(new MoleFractionLiquid(elist));
    S_->SetEvaluator(x_liquid_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(x_liquid_key_, Tags::DEFAULT,
                                       x_gas_key_, Tags::DEFAULT, x_liquid_key_);
    S_->RequireDerivative<CV_t, CVS_t>(x_liquid_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, x_liquid_key_);
  }

  // product evaluators
  // -- coefficient for advection operator (div eta_l q_l) in liquid phase 
  if (!S_->HasRecord(advection_liquid_key_)) {
    auto elist = MyRequire_(advection_liquid_key_, advection_liquid_key_);

    Teuchos::Array<int> dep_powers(4, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(mol_density_liquid_key_);
    dep_names.push_back(x_liquid_key_);
    dep_names.push_back(relperm_liquid_key_);
    dep_names.push_back(viscosity_liquid_key_);
    dep_powers[3] = -1;

    elist.set<Teuchos::Array<std::string> >("dependencies", dep_names) 
         .set<bool>("dependency tags are my tag", true)
         .set<Teuchos::Array<int> >("powers", dep_powers);
    // elist.sublist("verbose object").set<std::string>("verbosity level", "extreme");

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S_->SetEvaluator(advection_liquid_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(advection_liquid_key_, Tags::DEFAULT,
                                       pressure_liquid_key_, Tags::DEFAULT, advection_liquid_key_);
    S_->RequireDerivative<CV_t, CVS_t>(advection_liquid_key_, Tags::DEFAULT,
                                       x_gas_key_, Tags::DEFAULT, advection_liquid_key_);
    S_->RequireDerivative<CV_t, CVS_t>(advection_liquid_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, advection_liquid_key_);
  }

  if (!S_->HasRecord(advection_water_key_)) {
    auto elist = MyRequire_(advection_water_key_, advection_water_key_);

    Teuchos::Array<int> dep_powers(3, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(mol_density_liquid_key_);
    dep_names.push_back(relperm_liquid_key_);
    dep_names.push_back(viscosity_liquid_key_);
    dep_powers[2] = -1;

    elist.set<Teuchos::Array<std::string> >("dependencies", dep_names) 
         .set<bool>("dependency tags are my tag", true)
         .set<Teuchos::Array<int> >("powers", dep_powers);

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S_->SetEvaluator(advection_water_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(advection_water_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, advection_water_key_);
  }

  // -- coefficient for advection operator in gas phase (4 fields)
  if (!S_->HasRecord(advection_gas_key_)) {
    auto elist = MyRequire_(advection_gas_key_, advection_gas_key_);

    Teuchos::Array<int> dep_powers(4, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(mol_density_gas_key_);
    dep_names.push_back(x_gas_key_);
    dep_names.push_back(relperm_gas_key_);
    dep_names.push_back(viscosity_gas_key_);
    dep_powers[3] = -1;

    elist.set<Teuchos::Array<std::string> >("dependencies", dep_names) 
         .set<bool>("dependency tags are my tag", true)
         .set<Teuchos::Array<int> >("powers", dep_powers);

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S_->SetEvaluator(advection_gas_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(advection_gas_key_, Tags::DEFAULT,
                                       pressure_liquid_key_, Tags::DEFAULT, advection_gas_key_);
    S_->RequireDerivative<CV_t, CVS_t>(advection_gas_key_, Tags::DEFAULT,
                                       x_gas_key_, Tags::DEFAULT, advection_gas_key_);
    S_->RequireDerivative<CV_t, CVS_t>(advection_gas_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, advection_gas_key_);
  }

  // -- coefficient for diffusion operator in liquid phase
  if (!S_->HasRecord(diffusion_liquid_key_)) {
    auto elist = MyRequire_(diffusion_liquid_key_, diffusion_liquid_key_);

    Teuchos::Array<int> dep_powers(3, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(molecular_diff_liquid_key_);
    dep_names.push_back(mol_density_liquid_key_);
    dep_names.push_back(saturation_liquid_key_);

    elist.set<Teuchos::Array<std::string> >("dependencies", dep_names) 
         .set<bool>("dependency tags are my tag", true)
         .set<Teuchos::Array<int> >("powers", dep_powers);

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S_->SetEvaluator(diffusion_liquid_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(diffusion_liquid_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, diffusion_liquid_key_);
  }

  // -- coefficient for diffusion operator in gas phase
  if (!S_->HasRecord(diffusion_gas_key_)) {
    auto elist = MyRequire_(diffusion_gas_key_, diffusion_gas_key_);

    Teuchos::Array<int> dep_powers(3, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(molecular_diff_gas_key_);
    dep_names.push_back(mol_density_gas_key_);
    dep_names.push_back(saturation_gas_key_);

    elist.set<Teuchos::Array<std::string> >("dependencies", dep_names) 
         .set<bool>("dependency tags are my tag", true)
         .set<Teuchos::Array<int> >("powers", dep_powers);

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S_->SetEvaluator(diffusion_gas_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(diffusion_gas_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, diffusion_gas_key_);
  }

  // mole fraction of water vapor
  if (!S_->HasRecord(x_vapor_key_)) {
    auto elist = MyRequire_(x_vapor_key_, x_vapor_key_);

    Teuchos::Array<int> dep_powers(2, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(pressure_gas_key_);
    dep_names.push_back(pressure_vapor_key_);
    dep_powers[0] = -1;

    elist.set<Teuchos::Array<std::string> >("dependencies", dep_names) 
         .set<bool>("dependency tags are my tag", true)
         .set<Teuchos::Array<int> >("powers", dep_powers);

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S_->SetEvaluator(x_vapor_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(x_vapor_key_, Tags::DEFAULT,
                                       pressure_liquid_key_, Tags::DEFAULT, x_vapor_key_);
    S_->RequireDerivative<CV_t, CVS_t>(x_vapor_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, x_vapor_key_);
  }

  // nonlinear complimentary problem
  if (!S_->HasRecord(ncp_f_key_)) {
    auto elist = MyRequire_(ncp_f_key_, ncp_f_key_);
    elist.set<std::string>("saturation liquid key", saturation_liquid_key_);
    auto eval = Teuchos::rcp(new NCP_F(elist));
    S_->SetEvaluator(ncp_f_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(ncp_f_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, ncp_f_key_);
  }

  if (!S_->HasRecord(ncp_g_key_)) {
    auto elist = MyRequire_(ncp_g_key_, ncp_g_key_);
    elist.set<std::string>("mole fraction vapor key", x_vapor_key_)
         .set<std::string>("mole fraction gas key", x_gas_key_);
    auto eval = Teuchos::rcp(new NCP_MoleFractions(elist));
    S_->SetEvaluator(ncp_g_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(ncp_g_key_, Tags::DEFAULT,
                                       pressure_liquid_key_, Tags::DEFAULT, ncp_g_key_);
    S_->RequireDerivative<CV_t, CVS_t>(ncp_g_key_, Tags::DEFAULT,
                                       x_gas_key_, Tags::DEFAULT, ncp_g_key_);
    S_->RequireDerivative<CV_t, CVS_t>(ncp_g_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, ncp_g_key_);
  }
}


/* ******************************************************************
* Initialize various PK objects
****************************************************************** */
void MultiphaseModel1_PK::Initialize()
{
  Multiphase_PK::Initialize();
  if (non_isothermal_) {
    InitializeFieldFromField_(prev_energy_key_, energy_key_, true);
  }
}


/* ******************************************************************
* Push data to the state
****************************************************************** */
void MultiphaseModel1_PK::CommitStep(
    double t_old, double t_new, const Tag& tag)
{
  Multiphase_PK::CommitStep(t_old, t_new, tag);

  // miscalleneous fields
  S_->GetEvaluator(ncp_fg_key_).Update(*S_, passwd_);
}


/********************************************************************
* Modifies nonlinear update du using .. TBW
****************************************************************** */
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
MultiphaseModel1_PK::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                                      Teuchos::RCP<const TreeVector> u,
                                      Teuchos::RCP<TreeVector> du)
{
  // clip liquid pressure to stays positive
  const auto& u0c = *u->SubVector(0)->Data()->ViewComponent("cell");
  auto& du0c = *du->SubVector(0)->Data()->ViewComponent("cell");

  for (int c = 0; c < ncells_owned_; ++c) {
    du0c[0][c] = std::min(du0c[0][c], u0c[0][c]);
  }

  // clip mole fraction to range [0; 1]
  /*
  const auto& u1c = *u->SubVector(1)->Data()->ViewComponent("cell");
  auto& du1c = *du->SubVector(1)->Data()->ViewComponent("cell");

  for (int i = 0; i < u1c.NumVectors(); ++i) {
    for (int c = 0; c < ncells_owned_; ++c) {
      du1c[i][c] = std::min(du1c[i][c], u1c[i][c]);
      du1c[i][c] = std::max(du1c[i][c], u1c[i][c] - 1.0);
    }    
  }
  */

  // clip saturation (residual saturation is missing, FIXME)
  /*
  const auto& u2c = *u->SubVector(2)->Data()->ViewComponent("cell");
  auto& du2c = *du->SubVector(2)->Data()->ViewComponent("cell");

  for (int c = 0; c < ncells_owned_; ++c) {
    du2c[0][c] = std::min(du2c[0][c], u2c[0][c]);
    du2c[0][c] = std::max(du2c[0][c], u2c[0][c] - 1.0);
  }
  */

  return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED;
}


/* ******************************************************************* 
* Tweak evaluators.
******************************************************************* */
void MultiphaseModel1_PK::ModifyEvaluators(int neqn)
{
  int n0 = (non_isothermal_) ? 1 : 0;
 
  if (neqn > n0) {
    int ifield(0), n(neqn - 1 - n0);
    // ifield is dummy here
    Teuchos::rcp_dynamic_cast<TotalComponentStorage>(eval_tcs_)->set_subvector(ifield, n, kH_[n]);
    Teuchos::rcp_dynamic_cast<TotalComponentStorage>(eval_tcs_)->Update(*S_, passwd_, true);

    auto eval = S_->GetEvaluatorPtr(x_liquid_key_, Tags::DEFAULT);
    Teuchos::rcp_dynamic_cast<MoleFractionLiquid>(eval)->set_subvector(ifield, n, kH_[n]);
    Teuchos::rcp_dynamic_cast<MoleFractionLiquid>(eval)->Update(*S_, passwd_, true);
  }
}

}  // namespace Multiphase
}  // namespace Amanzi

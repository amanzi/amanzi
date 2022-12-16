/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  MultiPhase PK

  Multiphase multi-component flow. We assume that one component (water)
  is always in liquid form, i,e does not evaporate. This allows us to
  compare results with other codes. Later, this PK will be split into
  a base and derived PKs depending on physical models.
*/


// TPLs
#include "Teuchos_RCP.hpp"

// Amanzi
#include "EnthalpyEvaluator.hh"
#include "EvaluatorPrimary.hh"
#include "InverseFactory.hh"
#include "IO.hh"
#include "PDE_Accumulation.hh"
#include "PDE_AdvectionUpwind.hh"
#include "PDE_DiffusionFactory.hh"
#include "PK_DomainFunctionFactory.hh"
#include "RelPermEvaluator.hh"
#include "TCMEvaluator_TwoPhase.hh"
#include "UpwindFactory.hh"

// Multiphase
#include "ModelMeshPartition.hh"
#include "MoleFractionLiquid.hh"
#include "Multiphase_PK.hh"
#include "Multiphase_Utils.hh"
#include "MultiphaseTypeDefs.hh"
#include "NCP_F.hh"
#include "NCP_HenryLaw.hh"
#include "PressureGasEvaluator.hh"
#include "ProductEvaluator.hh"
#include "TotalComponentStorage.hh"
#include "TotalEnergyEvaluator.hh"
#include "VaporPressureEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* Standard constructor
****************************************************************** */
Multiphase_PK::Multiphase_PK(Teuchos::ParameterList& pk_tree,
                             const Teuchos::RCP<Teuchos::ParameterList>& glist,
                             const Teuchos::RCP<State>& S,
                             const Teuchos::RCP<TreeVector>& soln)
  : passwd_(""), soln_(soln), num_phases_(2), op_pc_assembled_(false), glist_(glist), num_itrs_(0)
{
  S_ = S;

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

  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = mp_list_->sublist("verbose object");
  std::string ioname = "MultiphasePK" + ((domain_ != "domain") ? "-" + domain_ : "");
  vo_ = Teuchos::rcp(new VerboseObject(ioname, vlist));
}


/* ******************************************************************
* Setup
****************************************************************** */
void
Multiphase_PK::Setup()
{
  mesh_ = S_->GetMesh(domain_);
  dim_ = mesh_->space_dimension();

  // keys
  // -- fields
  pressure_liquid_key_ = Keys::getKey(domain_, "pressure_liquid");
  pressure_gas_key_ = Keys::getKey(domain_, "pressure_gas");
  pressure_vapor_key_ = Keys::getKey(domain_, "pressure_vapor");

  saturation_liquid_key_ = Keys::getKey(domain_, "saturation_liquid");
  saturation_gas_key_ = Keys::getKey(domain_, "saturation_gas");

  temperature_key_ = Keys::getKey(domain_, "temperature");
  energy_key_ = Keys::getKey(domain_, "energy");
  prev_energy_key_ = Keys::getKey(domain_, "prev_energy");

  x_liquid_key_ = Keys::getKey(domain_, "mole_fraction_liquid");
  x_gas_key_ = Keys::getKey(domain_, "mole_fraction_gas");

  tcc_liquid_key_ = Keys::getKey(domain_, "total_component_concentration_liquid");
  tcc_gas_key_ = Keys::getKey(domain_, "total_component_concentration_gas");

  permeability_key_ = Keys::getKey(domain_, "permeability");
  porosity_key_ = Keys::getKey(domain_, "porosity");

  relperm_liquid_key_ = Keys::getKey(domain_, "rel_permeability_liquid");
  relperm_gas_key_ = Keys::getKey(domain_, "rel_permeability_gas");
  viscosity_liquid_key_ = Keys::getKey(domain_, "viscosity_liquid");
  viscosity_gas_key_ = Keys::getKey(domain_, "viscosity_gas");

  advection_liquid_key_ = Keys::getKey(domain_, "advection_liquid");
  advection_gas_key_ = Keys::getKey(domain_, "advection_gas");
  mol_density_liquid_key_ = Keys::getKey(domain_, "molar_density_liquid");
  mol_density_gas_key_ = Keys::getKey(domain_, "molar_density_gas");
  mass_density_liquid_key_ = Keys::getKey(domain_, "mass_density_liquid");
  mass_density_gas_key_ = Keys::getKey(domain_, "mass_density_gas");

  vol_flowrate_liquid_key_ = Keys::getKey(domain_, "volumetric_flow_rate_liquid");
  vol_flowrate_gas_key_ = Keys::getKey(domain_, "volumetric_flow_rate_gas");

  tws_key_ = Keys::getKey(domain_, "total_water_storage");
  tcs_key_ = Keys::getKey(domain_, "total_component_storage");
  prev_tws_key_ = Keys::getKey(domain_, "prev_total_water_storage");
  prev_tcs_key_ = Keys::getKey(domain_, "prev_total_component_storage");

  ncp_f_key_ = Keys::getKey(domain_, "ncp_f");
  ncp_g_key_ = Keys::getKey(domain_, "ncp_g");
  ncp_fg_key_ = Keys::getKey(domain_, "ncp_fg");

  ie_rock_key_ = Keys::getKey(domain_, "internal_energy_rock");
  ie_liquid_key_ = Keys::getKey(domain_, "internal_energy_liquid");
  ie_gas_key_ = Keys::getKey(domain_, "internal_energy_gas");
  conductivity_key_ = Keys::getKey(domain_, "thermal_conductivity");
  particle_density_key_ = Keys::getKey(domain_, "particle_density");
  enthalpy_liquid_key_ = Keys::getKey(domain_, "enthalpy_liquid");
  enthalpy_gas_key_ = Keys::getKey(domain_, "enthalpy_gas");

  advection_enthalpy_liquid_key_ = Keys::getKey(domain_, "advection_enthalpy_liquid");
  advection_enthalpy_gas_key_ = Keys::getKey(domain_, "advection_enthalpy_gas");

  // model assumptions
  auto physical_models = Teuchos::sublist(mp_list_, "physical models and assumptions");
  flow_on_manifold_ = physical_models->get<bool>("flow and transport in fractures", false);

  // extract information about nonlinear system
  component_names_ = mp_list_->sublist("molecular diffusion")
                       .get<Teuchos::Array<std::string>>("aqueous names")
                       .toVector();
  num_primary_ = component_names_.size();

  int id0(0), id1, id2, id3, id4;
  id1 = InitMPSystem_("pressure eqn", id0, 1);
  id2 = InitMPSystem_("energy eqn", id1, 1);
  id3 = InitMPSystem_("solute eqn", id2, num_primary_);
  id4 = InitMPSystem_("constraint eqn", id3, 1);

  // register non-standard fields
  if (!S_->HasRecord("gravity"))
    S_->Require<AmanziGeometry::Point>("gravity", Tags::DEFAULT, "state");

  if (!S_->HasRecord("const_fluid_density"))
    S_->Require<double>("const_fluid_density", Tags::DEFAULT, "state");

  if (!S_->HasRecord("const_fluid_viscosity"))
    S_->Require<double>("const_fluid_viscosity", Tags::DEFAULT, "state");

  if (!S_->HasRecord("const_gas_viscosity"))
    S_->Require<double>("const_gas_viscosity", Tags::DEFAULT, "state");

  if (!S_->HasRecord("atmospheric_pressure"))
    S_->Require<double>("atmospheric_pressure", Tags::DEFAULT, "state");

  // primary unknowns
  // -- pressure liquid
  if (!S_->HasRecord(pressure_liquid_key_)) {
    auto elist = MyRequire_(pressure_liquid_key_, passwd_);
    auto eval = Teuchos::rcp(new EvaluatorPrimary<CV_t, CVS_t>(elist));
    S_->SetEvaluator(pressure_liquid_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(pressure_liquid_key_,
                                       Tags::DEFAULT,
                                       pressure_liquid_key_,
                                       Tags::DEFAULT,
                                       pressure_liquid_key_)
      .SetGhosted();
  }

  // -- temperature
  if (!S_->HasRecord(temperature_key_)) {
    if (mp_list_->sublist("system").isSublist("energy eqn")) {
      auto elist = MyRequire_(temperature_key_, passwd_);
      auto eval = Teuchos::rcp(new EvaluatorPrimary<CV_t, CVS_t>(elist));
      S_->SetEvaluator(temperature_key_, Tags::DEFAULT, eval);
    } else {
      MyRequire_(temperature_key_, temperature_key_);
      S_->RequireEvaluator(temperature_key_, Tags::DEFAULT);
    }
  }

  // -- other primary unknowns
  std::string primary_key =
    mp_list_->sublist("system").sublist("solute eqn").get<std::string>("primary unknown");

  if (!S_->HasRecord(primary_key)) {
    S_->Require<CV_t, CVS_t>(primary_key, Tags::DEFAULT, passwd_, component_names_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, component_names_.size());

    Teuchos::ParameterList elist(primary_key);
    elist.set<std::string>("name", primary_key).set<std::string>("tag", "");
    auto eval = Teuchos::rcp(new EvaluatorPrimary<CV_t, CVS_t>(elist));
    S_->SetEvaluator(primary_key, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(
        primary_key, Tags::DEFAULT, primary_key, Tags::DEFAULT, primary_key)
      .SetGhosted();
  }

  // conserved quantities
  // -- total water storage

  // water saturation is the primary solution
  if (!S_->HasRecord(saturation_liquid_key_)) {
    auto elist = MyRequire_(saturation_liquid_key_, passwd_);
    saturation_liquid_eval_ = Teuchos::rcp(new EvaluatorPrimary<CV_t, CVS_t>(elist));
    S_->SetEvaluator(saturation_liquid_key_, Tags::DEFAULT, saturation_liquid_eval_);

    S_->RequireDerivative<CV_t, CVS_t>(saturation_liquid_key_,
                                       Tags::DEFAULT,
                                       saturation_liquid_key_,
                                       Tags::DEFAULT,
                                       saturation_liquid_key_)
      .SetGhosted();
  }

  // total pressure gas
  auto wrm_list = Teuchos::sublist(mp_list_, "water retention models", true);
  wrm_ = CreateModelPartition<WRMmp>(mesh_, wrm_list, "water retention model");

  if (!S_->HasRecord(pressure_gas_key_)) {
    auto elist = MyRequire_(pressure_gas_key_, pressure_gas_key_);
    elist.set<std::string>("pressure liquid key", pressure_liquid_key_)
      .set<std::string>("saturation liquid key", saturation_liquid_key_);

    auto eval = Teuchos::rcp(new PressureGasEvaluator(elist, wrm_));
    S_->SetEvaluator(pressure_gas_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(
        pressure_gas_key_, Tags::DEFAULT, pressure_liquid_key_, Tags::DEFAULT, pressure_gas_key_)
      .SetGhosted();
    S_->RequireDerivative<CV_t, CVS_t>(
        pressure_gas_key_, Tags::DEFAULT, pressure_liquid_key_, Tags::DEFAULT, pressure_gas_key_)
      .SetGhosted();
    S_->RequireDerivative<CV_t, CVS_t>(
        pressure_gas_key_, Tags::DEFAULT, saturation_liquid_key_, Tags::DEFAULT, pressure_gas_key_)
      .SetGhosted();
  }

  // Darcy volume fluxes
  if (!S_->HasRecord(vol_flowrate_liquid_key_)) {
    if (!flow_on_manifold_) {
      S_->Require<CV_t, CVS_t>(vol_flowrate_liquid_key_, Tags::DEFAULT, passwd_)
        .SetMesh(mesh_)
        ->SetGhosted(true)
        ->SetComponent("face", AmanziMesh::FACE, 1);
    } else {
      auto cvs = Operators::CreateManifoldCVS(mesh_);
      *S_->Require<CV_t, CVS_t>(vol_flowrate_liquid_key_, Tags::DEFAULT, passwd_)
         .SetMesh(mesh_)
         ->SetGhosted(true) = *cvs;
    }
 }

  if (!S_->HasRecord(vol_flowrate_gas_key_)) {
    if (!flow_on_manifold_) {
      S_->Require<CV_t, CVS_t>(vol_flowrate_gas_key_, Tags::DEFAULT, passwd_)
        .SetMesh(mesh_)
        ->SetGhosted(true)
        ->SetComponent("face", AmanziMesh::FACE, 1);
    } else {
      auto cvs = Operators::CreateManifoldCVS(mesh_);
      *S_->Require<CV_t, CVS_t>(vol_flowrate_gas_key_, Tags::DEFAULT, passwd_)
         .SetMesh(mesh_)
         ->SetGhosted(true) = *cvs;
    }
  }

  // densities
  if (!S_->HasRecord(mass_density_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(mass_density_liquid_key_, Tags::DEFAULT, mass_density_liquid_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  }

  // viscosities
  if (!S_->HasRecord(viscosity_liquid_key_)) {
    MyRequire_(viscosity_liquid_key_, viscosity_liquid_key_);
    S_->RequireEvaluator(viscosity_liquid_key_, Tags::DEFAULT);
  }

  if (!S_->HasRecord(viscosity_gas_key_)) {
    MyRequire_(viscosity_gas_key_, viscosity_gas_key_);
    S_->RequireEvaluator(viscosity_gas_key_, Tags::DEFAULT);
  }

  // relative permeability of liquid phase
  if (!S_->HasRecord(relperm_liquid_key_)) {
    auto elist = MyRequire_(relperm_liquid_key_, relperm_liquid_key_);
    elist.set<std::string>("saturation liquid key", saturation_liquid_key_)
      .set<std::string>("phase name", "liquid");

    auto eval = Teuchos::rcp(new RelPermEvaluator(elist, wrm_));
    S_->SetEvaluator(relperm_liquid_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(relperm_liquid_key_,
                                       Tags::DEFAULT,
                                       saturation_liquid_key_,
                                       Tags::DEFAULT,
                                       relperm_liquid_key_)
      .SetGhosted();
  }

  // relative permeability of gas phase
  if (!S_->HasRecord(relperm_gas_key_)) {
    auto elist = MyRequire_(relperm_gas_key_, relperm_gas_key_);
    elist.set<std::string>("saturation liquid key", saturation_liquid_key_)
      .set<std::string>("phase name", "gas");

    auto eval = Teuchos::rcp(new RelPermEvaluator(elist, wrm_));
    S_->SetEvaluator(relperm_gas_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(
        relperm_gas_key_, Tags::DEFAULT, saturation_liquid_key_, Tags::DEFAULT, relperm_gas_key_)
      .SetGhosted();
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

  // material properties
  if (!S_->HasRecord(permeability_key_)) {
    S_->Require<CV_t, CVS_t>(permeability_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim_);
  }

  if (!S_->HasRecord(porosity_key_)) {
    MyRequire_(porosity_key_, porosity_key_);
    S_->RequireEvaluator(porosity_key_, Tags::DEFAULT);
  }

  // non-isothermal model
  if (system_["energy eqn"]) {
    if (!S_->HasRecord(energy_key_)) {
      S_->Require<CV_t, CVS_t>(energy_key_, Tags::DEFAULT, energy_key_)
        .SetMesh(mesh_)
        ->SetGhosted()
        ->AddComponent("cell", AmanziMesh::CELL, 1);

      Teuchos::ParameterList elist = mp_list_->sublist("energy evaluator");
      elist.set<std::string>("energy key", energy_key_)
        .set<std::string>("tag", "")
        .set<bool>("vapor diffusion", true)
        .set<std::string>("particle density key", particle_density_key_)
        .set<std::string>("internal energy rock key", ie_rock_key_);
      elist.setName(energy_key_);

      auto ee = Teuchos::rcp(new Energy::TotalEnergyEvaluator(elist));
      S_->SetEvaluator(energy_key_, Tags::DEFAULT, ee);

      S_->RequireDerivative<CV_t, CVS_t>(
          energy_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT, energy_key_)
        .SetGhosted();
    }

    // previous conserved quantities
    if (!S_->HasRecord(prev_energy_key_)) {
      S_->Require<CV_t, CVS_t>(prev_energy_key_, Tags::DEFAULT, passwd_)
        .SetMesh(mesh_)
        ->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
      S_->GetRecordW(prev_energy_key_, passwd_).set_io_vis(false);
    }

    // internal energies
    MyRequire_(ie_liquid_key_, ie_liquid_key_);
    S_->RequireEvaluator(ie_liquid_key_, Tags::DEFAULT);

    S_->RequireDerivative<CV_t, CVS_t>(
        ie_liquid_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT, ie_liquid_key_)
      .SetGhosted();

    MyRequire_(ie_gas_key_, ie_gas_key_);
    S_->RequireEvaluator(ie_gas_key_, Tags::DEFAULT);

    S_->RequireDerivative<CV_t, CVS_t>(
        ie_gas_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT, ie_gas_key_)
      .SetGhosted();

    MyRequire_(ie_rock_key_, ie_rock_key_);
    S_->RequireEvaluator(ie_rock_key_, Tags::DEFAULT);

    S_->RequireDerivative<CV_t, CVS_t>(
        ie_rock_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT, ie_rock_key_)
      .SetGhosted();

    // -- thermal conductivity
    if (!S_->HasRecord(conductivity_key_)) {
      S_->Require<CV_t, CVS_t>(conductivity_key_, Tags::DEFAULT, conductivity_key_)
        .SetMesh(mesh_)
        ->SetGhosted()
        ->AddComponent("cell", AmanziMesh::CELL, 1);

      Teuchos::ParameterList elist = mp_list_->sublist("thermal conductivity evaluator");
      elist.set("thermal conductivity key", conductivity_key_).set<std::string>("tag", "");
      elist.setName(conductivity_key_);

      auto tcm = Teuchos::rcp(new Energy::TCMEvaluator_TwoPhase(elist));
      S_->SetEvaluator(conductivity_key_, Tags::DEFAULT, tcm);
    }

    // -- advection of enthalpy
    if (!S_->HasRecord(enthalpy_liquid_key_)) {
      S_->Require<CV_t, CVS_t>(enthalpy_liquid_key_, Tags::DEFAULT, enthalpy_liquid_key_)
        .SetMesh(mesh_)
        ->SetGhosted()
        ->AddComponent("cell", AmanziMesh::CELL, 1);

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

      S_->RequireDerivative<CV_t, CVS_t>(enthalpy_liquid_key_,
                                         Tags::DEFAULT,
                                         temperature_key_,
                                         Tags::DEFAULT,
                                         enthalpy_liquid_key_)
        .SetGhosted();
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

      elist.set<Teuchos::Array<std::string>>("dependencies", dep_names)
        .set<bool>("dependency tags are my tag", true)
        .set<Teuchos::Array<int>>("powers", dep_powers);

      auto eval = Teuchos::rcp(new ProductEvaluator(elist));
      S_->SetEvaluator(advection_enthalpy_liquid_key_, Tags::DEFAULT, eval);
    }
  }

  // fields from previous time step
  if (!S_->HasRecord(prev_tws_key_)) {
    S_->Require<CV_t, CVS_t>(prev_tws_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->GetRecordW(prev_tws_key_, passwd_).set_io_vis(false);
  }
  if (!S_->HasRecord(prev_tcs_key_)) {
    S_->Require<CV_t, CVS_t>(prev_tcs_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->GetRecordW(prev_tcs_key_, passwd_).set_io_vis(false);
  }

  // nonlinear complimentary problem fields
  if (id4 > id3) {
    if (!S_->HasRecord(ncp_f_key_)) {
      S_->Require<CV_t, CVS_t>(ncp_f_key_, Tags::DEFAULT, ncp_f_key_)
        .SetMesh(mesh_)
        ->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);

      Teuchos::ParameterList elist(ncp_f_key_);
      elist.set<std::string>("my key", ncp_f_key_)
        .set<std::string>("saturation liquid key", saturation_liquid_key_)
        .set<std::string>("tag", "");
      auto eval = Teuchos::rcp(new NCP_F(elist));
      S_->SetEvaluator(ncp_f_key_, Tags::DEFAULT, eval);

      S_->RequireDerivative<CV_t, CVS_t>(
          ncp_f_key_, Tags::DEFAULT, saturation_liquid_key_, Tags::DEFAULT, ncp_f_key_)
        .SetGhosted();
    }

    if (!S_->HasRecord(ncp_fg_key_)) {
      auto elist = MyRequire_(ncp_fg_key_, ncp_fg_key_);

      Teuchos::Array<int> dep_powers(2, 1);
      Teuchos::Array<std::string> dep_names;

      dep_names.push_back(ncp_f_key_);
      dep_names.push_back(ncp_g_key_);

      elist.set<Teuchos::Array<std::string>>("dependencies", dep_names)
        .set<bool>("dependency tags are my tag", true)
        .set<Teuchos::Array<int>>("powers", dep_powers);

      auto eval = Teuchos::rcp(new ProductEvaluator(elist));
      S_->SetEvaluator(ncp_fg_key_, Tags::DEFAULT, eval);
    }
  }

  // additional evaluators
  if (mp_list_->isParameter("evaluators")) {
    auto evals = mp_list_->get<Teuchos::Array<std::string>>("evaluators").toVector();
    for (auto& e : evals) {
      S_->Require<CV_t, CVS_t>(e, Tags::DEFAULT, e)
        .SetMesh(mesh_)
        ->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
      S_->RequireEvaluator(e, Tags::DEFAULT);
    }
  }

  // derivatives with respect to primary
  for (const auto& key : eval_flattened_) {
    // require dependent evaluators to enable queering of dependency tree
    S_->GetEvaluator(key).EnsureEvaluators(*S_);

    for (const auto& name : soln_names_) {
      if (S_->GetEvaluator(key).IsDependency(*S_, name, Tags::DEFAULT)) {
        S_->RequireDerivative<CV_t, CVS_t>(key, Tags::DEFAULT, name, Tags::DEFAULT, key)
          .SetGhosted();
      }
    }
  }
}


/* ******************************************************************
* Initialize various PK objects
****************************************************************** */
void
Multiphase_PK::Initialize()
{
  ncells_owned_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncells_wghost_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  nfaces_owned_ = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nfaces_wghost_ = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  ncp_ = mp_list_->get<std::string>("NCP function", "min");
  smooth_mu_ = mp_list_->get<double>("smoothing parameter mu", 0.0);

  // some defaults
  flux_names_ = { vol_flowrate_liquid_key_, vol_flowrate_gas_key_ };

  auto aux_list = mp_list_->sublist("molecular diffusion");
  mol_diff_l_ = aux_list.get<Teuchos::Array<double>>("aqueous values").toVector();
  mol_diff_g_ = aux_list.get<Teuchos::Array<double>>("gaseous values").toVector();
  mol_mass_ = aux_list.get<Teuchos::Array<double>>("molar masses").toVector();
  kH_ = aux_list.get<Teuchos::Array<double>>("Henry dimensionless constants").toVector();

  mol_mass_H2O_ = mp_list_->get<double>("molar mass of water");

  // fundamental physical quantities
  gravity_ = S_->Get<AmanziGeometry::Point>("gravity");
  g_ = fabs(gravity_[dim_ - 1]);

  rho_l_ = S_->Get<double>("const_fluid_density");
  eta_l_ = rho_l_ / CommonDefs::MOLAR_MASS_H2O;
  mu_l_ = S_->Get<double>("const_fluid_viscosity");
  mu_g_ = S_->Get<double>("const_gas_viscosity");

  // process CPR list
  cpr_enhanced_ = mp_list_->isSublist("CPR enhancement");
  if (cpr_enhanced_) {
    auto& cpr_list_ = mp_list_->sublist("CPR parameters");
    std::vector<int> correction_blocks =
      cpr_list_.get<Teuchos::Array<int>>("correction blocks").toVector();
    auto block_names_ = cpr_list_.get<Teuchos::Array<std::string>>("preconditioner").toVector();
    for (int i = 0; i < block_names_.size(); i++)
      AMANZI_ASSERT(pc_list_->isSublist(block_names_[i]));
  }

  // create system structure using pairs of evaluators and scalar factors
  // create solution vector
  for (const auto& primary_name : soln_names_) {
    auto field = Teuchos::rcp(new TreeVector());
    soln_->PushBack(field);
    field->SetData(S_->GetPtrW<CompositeVector>(primary_name, Tags::DEFAULT, passwd_));
  }

  // boundary conditions
  Teuchos::RCP<MultiphaseBoundaryFunction> bc;
  auto& bc_list = mp_list_->sublist("boundary conditions");

  bcs_.clear();

  // -- pressure
  if (bc_list.isSublist("pressure liquid")) {
    PK_DomainFunctionFactory<MultiphaseBoundaryFunction> bc_factory(mesh_, S_);

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
    PK_DomainFunctionFactory<MultiphaseBoundaryFunction> bc_factory(mesh_, S_);

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
    PK_DomainFunctionFactory<MultiphaseBoundaryFunction> bc_factory(mesh_, S_);

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

  // -- concentration
  if (bc_list.isSublist("concentration")) {
    PK_DomainFunctionFactory<MultiphaseBoundaryFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list.sublist("concentration");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);
        bc = bc_factory.Create(spec, "boundary concentration", AmanziMesh::FACE, Teuchos::null);
        bc->set_bc_name("concentration");
        bc->SetComponentId(component_names_);
        bcs_.push_back(bc);

        AMANZI_ASSERT(bc->component_id() >= 0);
      }
    }
  }

  // allocate metadata and populate boundary conditions
  op_bcs_.clear();

  for (const auto& name : soln_names_) {
    auto op_bc =
      Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
    op_bcs_[name] = op_bc;
  }
  for (const auto& name : secondary_names_) {
    auto op_bc =
      Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
    op_bcs_[name] = op_bc;
  }

  double t_ini = S_->get_time();
  for (int i = 0; i < bcs_.size(); i++) {
    bcs_[i]->Compute(t_ini, t_ini);
    bcs_[i]->ComputeSubmodel(mesh_);
  }

  PopulateBCs(0, true);

  // matrix is used to simplify calculation of residual
  // -- absolute permeability
  ncells_wghost_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  K_.resize(ncells_wghost_);
  ConvertFieldToTensor(S_, dim_, permeability_key_, K_);

  // -- pre-process full tensor once to re-use it later for flux calculation
  auto& ddf_list = mp_list_->sublist("operators").sublist("diffusion operator").sublist("matrix");
  Teuchos::RCP<std::vector<WhetStone::Tensor>> Kptr = Teuchos::rcpFromRef(K_);

  Operators::PDE_DiffusionFactory opfactory(ddf_list, mesh_);
  opfactory.SetVariableTensorCoefficient(Kptr);
  opfactory.SetConstantGravitationalTerm(gravity_, rho_l_);

  pde_diff_K_.resize(2);
  pde_diff_K_[0] = opfactory.Create();
  pde_diff_K_[0]->SetBCs(op_bcs_[pressure_liquid_key_], op_bcs_[pressure_liquid_key_]);
  pde_diff_K_[0]->UpdateMatrices(Teuchos::null, Teuchos::null);

  auto eta_g = S_->GetPtr<CV_t>(mol_density_gas_key_, Tags::DEFAULT);
  opfactory.SetVariableGravitationalTerm(gravity_, eta_g);

  pde_diff_K_[1] = opfactory.Create();
  pde_diff_K_[1]->SetBCs(op_bcs_[pressure_gas_key_], op_bcs_[pressure_gas_key_]);
  pde_diff_K_[1]->UpdateMatrices(Teuchos::null, Teuchos::null);

  auto& mdf_list =
    mp_list_->sublist("operators").sublist("molecular diffusion operator").sublist("matrix");
  pde_diff_D_ = Teuchos::rcp(new Operators::PDE_DiffusionFV(mdf_list, mesh_));
  pde_diff_D_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // preconditioner is model-specific. It is created in the scope of global assembly to
  // reduce memory footprint for large number of components.
  std::string pc_name = mp_list_->get<std::string>("preconditioner");
  std::string ls_name = mp_list_->get<std::string>("linear solver");
  auto inv_list = AmanziSolvers::mergePreconditionerSolverLists(
    pc_name, *pc_list_, ls_name, *linear_operator_list_, true);
  inv_list.setName(pc_name);

  op_preconditioner_ =
    Teuchos::rcp(new Operators::FlattenedTreeOperator(Teuchos::rcpFromRef(soln_->Map())));
  op_preconditioner_->AddColoring(inv_list);
  op_pc_solver_ = AmanziSolvers::createInverse(inv_list, op_preconditioner_);

  // initialize time integrator
  std::string ti_method_name = ti_list_->get<std::string>("time integration method", "none");
  AMANZI_ASSERT(ti_method_name == "BDF1");
  Teuchos::ParameterList& bdf1_list = ti_list_->sublist("BDF1");

  if (!bdf1_list.isSublist("verbose object"))
    bdf1_list.sublist("verbose object") = mp_list_->sublist("verbose object");

  bdf1_dae_ = Teuchos::rcp(new BDF1_TI<TreeVector, TreeVectorSpace>(*this, bdf1_list, soln_));

  // upwind operator with a face model (FIXME)
  Operators::UpwindFactory upwfact;
  auto upw_list = mp_list_->sublist("operators").sublist("diffusion operator").sublist("upwind");
  upwind_ = upwfact.Create(mesh_, upw_list);

  // initialize other fields and evaluators
  S_->GetEvaluator(relperm_liquid_key_).Update(*S_, passwd_);
  S_->GetEvaluator(relperm_gas_key_).Update(*S_, passwd_);

  InitializeFields_();

  // io
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    for (const auto& name : soln_names_) { *vo_->os() << "unknown: \"" << name << "\"\n\n"; }
    for (int i = 0; i < bcs_.size(); i++) {
      *vo_->os() << "bc \"" << bcs_[i]->keyword() << "\" has " << bcs_[i]->size() << " entities"
                 << std::endl;
    }
    UpdatePreconditioner(0.0, soln_, 1.0);
    *vo_->os() << "preconditioner:" << std::endl
               << op_preconditioner_->PrintDiagnostics() << std::endl;
    *vo_->os() << vo_->color("green")
               << "Initialization of PK is complete, T=" << units_.OutputTime(S_->get_time())
               << vo_->reset() << std::endl
               << std::endl;
  }
}


/* ****************************************************************
* This completes initialization of fields left out by state.
**************************************************************** */
void
Multiphase_PK::InitializeFields_()
{
  Teuchos::OSTab tab = vo_->getOSTab();

  InitializeFieldFromField_(prev_tws_key_, tws_key_, true);
  InitializeFieldFromField_(prev_tcs_key_, tcs_key_, true);
  if (system_["energy eqn"]) InitializeFieldFromField_(prev_energy_key_, energy_key_, true);

  InitializeCVField(S_, *vo_, vol_flowrate_liquid_key_, Tags::DEFAULT, passwd_, 0.0);
  InitializeCVField(S_, *vo_, vol_flowrate_gas_key_, Tags::DEFAULT, passwd_, 0.0);
}


/* ****************************************************************
* Auxiliary initialization technique.
**************************************************************** */
void
Multiphase_PK::InitializeFieldFromField_(const std::string& field0,
                                         const std::string& field1,
                                         bool call_evaluator)
{
  if (S_->HasRecord(field0)) {
    if (!S_->GetRecord(field0).initialized()) {
      if (call_evaluator) S_->GetEvaluator(field1).Update(*S_, passwd_);

      const auto& f1 = S_->Get<CV_t>(field1);
      auto& f0 = S_->GetW<CV_t>(field0, Tags::DEFAULT, passwd_);
      f0 = f1;

      S_->GetRecordW(field0, Tags::DEFAULT, passwd_).set_initialized();

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
bool
Multiphase_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  dt_ = t_new - t_old;

  // make a copy of primary and conservative fields
  auto copy_names = soln_names_;
  copy_names.push_back(prev_tws_key_);
  copy_names.push_back(prev_tcs_key_);
  copy_names.push_back(vol_flowrate_liquid_key_);
  copy_names.push_back(vol_flowrate_gas_key_);
  if (system_["energy eqn"]) copy_names.push_back(prev_energy_key_);

  int ncopy = copy_names.size();
  std::vector<CompositeVector> copies;
  for (int i = 0; i < ncopy; ++i) { copies.push_back(S_->Get<CV_t>(copy_names[i])); }

  // initialization
  if (num_itrs_ == 0) {
    auto udot = Teuchos::rcp(new TreeVector(*soln_));
    udot->PutScalar(0.0);
    bdf1_dae_->SetInitialState(t_old, soln_, udot);
    num_itrs_++;
  }

  // update fields from previous time step
  S_->GetEvaluator(tws_key_).Update(*S_, passwd_);
  S_->GetW<CV_t>(prev_tws_key_, Tags::DEFAULT, passwd_) = S_->Get<CV_t>(tws_key_);

  S_->GetEvaluator(tcs_key_).Update(*S_, passwd_);
  S_->GetW<CV_t>(prev_tcs_key_, Tags::DEFAULT, passwd_) = S_->Get<CV_t>(tcs_key_);

  if (system_["energy eqn"]) {
    S_->GetEvaluator(energy_key_).Update(*S_, passwd_);
    S_->GetW<CV_t>(prev_energy_key_, Tags::DEFAULT, passwd_) = S_->Get<CV_t>(energy_key_);
  }

  // trying to make a step
  bool failed(false);
  failed = bdf1_dae_->TimeStep(dt_, dt_next_, soln_);
  if (failed) {
    dt_ = dt_next_;

    // recover the original fields
    for (int i = 0; i < ncopy; ++i) {
      S_->GetW<CV_t>(copy_names[i], Tags::DEFAULT, passwd_) = copies[i];
    }

    ChangedSolution();
    return failed;
  }

  bdf1_dae_->CommitSolution(t_new - t_old, soln_);
  ChangedSolution();

  dt_ = dt_next_;
  num_itrs_++;

  return failed;
}


/* ******************************************************************
* Push data to the state
****************************************************************** */
void
Multiphase_PK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  // no BC for upwind algorithms
  std::vector<int> bcnone(nfaces_wghost_, Operators::OPERATOR_BC_NONE);

  S_->GetEvaluator(advection_gas_key_).Update(*S_, passwd_);

  // work memory
  auto kr = CreateCVforUpwind_();
  auto& kr_c = *kr->ViewComponent("cell");

  Teuchos::RCP<std::vector<WhetStone::Tensor>> Kptr = Teuchos::rcpFromRef(K_);
  PopulateBCs(0, true);

  std::vector<std::string> relperm_name{ advection_liquid_key_, advection_gas_key_ };
  std::vector<std::string> viscosity_name{ viscosity_liquid_key_, viscosity_gas_key_ };
  std::vector<std::string> varp_name{ pressure_liquid_key_, pressure_gas_key_ };
  std::vector<std::string> flux_name{ vol_flowrate_liquid_key_, vol_flowrate_gas_key_ };

  for (int phase = 0; phase < 2; ++phase) {
    S_->GetEvaluator(relperm_name[phase]).Update(*S_, passwd_);

    const auto& relperm_c = *S_->Get<CV_t>(relperm_name[phase]).ViewComponent("cell");
    const auto& viscosity_c = *S_->Get<CV_t>(viscosity_name[phase]).ViewComponent("cell");
    auto flux = S_->GetPtrW<CV_t>(flux_name[phase], Tags::DEFAULT, passwd_);

    for (int c = 0; c < ncells_owned_; ++c) { kr_c[0][c] = relperm_c[0][c] / viscosity_c[0][c]; }
    upwind_->Compute(*flux, bcnone, *kr);

    auto pdeK = pde_diff_K_[phase];
    pdeK->Setup(Kptr, kr, Teuchos::null);
    pdeK->SetBCs(op_bcs_[varp_name[phase]], op_bcs_[varp_name[phase]]);
    pdeK->global_operator()->Init();
    pdeK->UpdateMatrices(Teuchos::null, Teuchos::null);

    auto var = S_->GetPtr<CV_t>(varp_name[phase], Tags::DEFAULT);
    pdeK->UpdateFlux(var.ptr(), flux.ptr());
  }

  // other fields
  if (S_->HasEvaluator(ncp_fg_key_, Tags::DEFAULT)) {
    S_->GetEvaluator(ncp_fg_key_).Update(*S_, passwd_);
  }
}


/********************************************************************
* Modifies nonlinear update du using .. TBW
****************************************************************** */
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
Multiphase_PK::ModifyCorrection(double h,
                                Teuchos::RCP<const TreeVector> res,
                                Teuchos::RCP<const TreeVector> u,
                                Teuchos::RCP<TreeVector> du)
{
  // clip molar density to range [0; +\infty]
  const auto& u1c = *u->SubVector(1)->Data()->ViewComponent("cell");
  auto& du1c = *du->SubVector(1)->Data()->ViewComponent("cell");

  for (int i = 0; i < u1c.NumVectors(); ++i) {
    for (int c = 0; c < ncells_owned_; ++c) {
      du1c[i][c] = std::min(du1c[i][c], u1c[i][c]);
      // du1c[i][c] = std::max(du1c[i][c], u1c[i][c] - 1.0);
    }
  }

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
* Parse system structure
******************************************************************* */
int
Multiphase_PK::InitMPSystem_(const std::string& eqn_name, int eqn_id, int eqn_num)
{
  system_[eqn_name] = false;
  const auto& plist = mp_list_->sublist("system");
  if (!plist.isSublist(eqn_name)) return eqn_id;

  auto slist = plist.sublist(eqn_name);
  system_[eqn_name] = true;

  std::string name, primary_name;
  std::vector<std::string> names;

  primary_name = slist.get<std::string>("primary unknown");
  soln_names_.push_back(primary_name);

  for (int i = 0; i < eqn_num; ++i) {
    int n = eqn_id + i;
    eqns_.resize(n + 1);
    eqns_flattened_.resize(n + 1);

    eqns_flattened_[n].resize(3);
    eqns_flattened_[n][0] = soln_names_.size() - 1;
    eqns_flattened_[n][1] = i;
    eqns_flattened_[n][2] = (eqn_num > 1) ? n : -1; // dEval_dSoln=0 for all eqns except this

    // advection evaluators
    if (slist.isParameter("advection factors"))
      eqns_[n].adv_factors = slist.get<Teuchos::Array<double>>("advection factors").toVector();

    if (slist.isParameter("advection liquid")) {
      names = slist.get<Teuchos::Array<std::string>>("advection liquid").toVector();
      eqns_[n].advection.push_back(std::make_pair(names[0], names[1]));
      eval_flattened_.push_back(names[0]);
      secondary_names_.insert(names[1]);
    } else {
      eqns_[n].advection.push_back(std::make_pair("", "")); // no liquid phase
    }

    if (slist.isParameter("advection gas")) {
      names = slist.get<Teuchos::Array<std::string>>("advection gas").toVector();
      eqns_[n].advection.push_back(std::make_pair(names[0], names[1]));
      eval_flattened_.push_back(names[0]);
      secondary_names_.insert(names[1]);
    } else {
      eqns_[n].advection.push_back(std::make_pair("", "")); // no gas phase
    }

    // diffusion evaluators
    if (slist.isParameter("diffusion factors"))
      eqns_[n].diff_factors = slist.get<Teuchos::Array<double>>("diffusion factors").toVector();

    if (slist.isParameter("diffusion liquid")) {
      names = slist.get<Teuchos::Array<std::string>>("diffusion liquid").toVector();
      eqns_[n].diffusion.push_back(std::make_pair(names[0], names[1]));
      eval_flattened_.push_back(names[0]);
      secondary_names_.insert(names[1]);
    } else {
      eqns_[n].diffusion.push_back(std::make_pair("", ""));
    }

    if (slist.isParameter("diffusion gas")) {
      names = slist.get<Teuchos::Array<std::string>>("diffusion gas").toVector();
      eqns_[n].diffusion.push_back(std::make_pair(names[0], names[1]));
      eval_flattened_.push_back(names[0]);
      secondary_names_.insert(names[1]);
    } else {
      eqns_[n].diffusion.push_back(std::make_pair("", ""));
    }

    // storage
    if (slist.isParameter("accumulation")) {
      name = slist.get<std::string>("accumulation");
      eqns_[n].storage = name;
      eval_flattened_.push_back(name);
    }

    // constraint
    if (slist.isParameter("ncp evaluators")) {
      names = slist.get<Teuchos::Array<std::string>>("ncp evaluators").toVector();
      eqns_[n].constraint = std::make_pair(names[0], names[1]);
      eval_flattened_.push_back(names[0]);
      eval_flattened_.push_back(names[1]);
    }
  }

  // remove solution vector from the list of secondary names
  for (const auto& soln : soln_names_) { secondary_names_.erase(soln); }

  return eqn_id + eqn_num;
}


/* *******************************************************************
* Populate boundary conditions for various bc types
******************************************************************* */
Teuchos::ParameterList
Multiphase_PK::MyRequire_(const Key& key, const std::string& owner)
{
  S_->Require<CV_t, CVS_t>(key, Tags::DEFAULT, owner)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::CELL, 1);

  Teuchos::ParameterList elist(key);
  elist.set<std::string>("my key", key).set<std::string>("tag", Tags::DEFAULT.get());

  return elist;
}


/* *******************************************************************
* Tweak evaluators.
******************************************************************* */
void
Multiphase_PK::ModifyEvaluators(int neqn)
{
  int n0 = (system_["energy eqn"]) ? 1 : 0;
  if (neqn > n0) {
    int ifield(0), n(neqn - 1 - n0);

    // ifield is dummy here
    if (S_->HasEvaluator(tcs_key_, Tags::DEFAULT)) {
      auto eval = S_->GetEvaluatorPtr(tcs_key_, Tags::DEFAULT);
      Teuchos::rcp_dynamic_cast<MultiphaseBaseEvaluator>(eval)->set_subvector(ifield, n, kH_[n]);
      Teuchos::rcp_dynamic_cast<MultiphaseBaseEvaluator>(eval)->Update(*S_, passwd_, true);
    }

    if (S_->HasEvaluator(ncp_g_key_, Tags::DEFAULT)) {
      auto eval = S_->GetEvaluatorPtr(ncp_g_key_, Tags::DEFAULT);
      Teuchos::rcp_dynamic_cast<MultiphaseBaseEvaluator>(eval)->set_subvector(ifield, n, kH_[n]);
      Teuchos::rcp_dynamic_cast<MultiphaseBaseEvaluator>(eval)->Update(*S_, passwd_, true);
    }

    if (S_->HasEvaluator(x_liquid_key_, Tags::DEFAULT)) {
      auto eval = S_->GetEvaluatorPtr(x_liquid_key_, Tags::DEFAULT);
      Teuchos::rcp_dynamic_cast<MoleFractionLiquid>(eval)->set_subvector(ifield, n, kH_[n]);
      Teuchos::rcp_dynamic_cast<MoleFractionLiquid>(eval)->Update(*S_, passwd_, true);
    }

    if (S_->HasEvaluator(tcc_liquid_key_, Tags::DEFAULT)) {
      auto eval = S_->GetEvaluatorPtr(tcc_liquid_key_, Tags::DEFAULT);
      if (eval->get_type() != EvaluatorType::PRIMARY) {
        Teuchos::rcp_dynamic_cast<MultiphaseBaseEvaluator>(eval)->set_subvector(ifield, n, kH_[n]);
        Teuchos::rcp_dynamic_cast<MultiphaseBaseEvaluator>(eval)->Update(*S_, passwd_, true);
      }
    }

    if (S_->HasEvaluator(tcc_gas_key_, Tags::DEFAULT)) {
      auto eval = S_->GetEvaluatorPtr(tcc_gas_key_, Tags::DEFAULT);
      if (eval->get_type() != EvaluatorType::PRIMARY) {
        Teuchos::rcp_dynamic_cast<MultiphaseBaseEvaluator>(eval)->set_subvector(ifield, n, kH_[n]);
        Teuchos::rcp_dynamic_cast<MultiphaseBaseEvaluator>(eval)->Update(*S_, passwd_, true);
      }
    }
  }
}


/* ******************************************************************* 
* Special cell-face structure for upwind field.
******************************************************************* */
Teuchos::RCP<CompositeVector> Multiphase_PK::CreateCVforUpwind_()
{
  auto cvs = S_->Get<CV_t>(vol_flowrate_liquid_key_, Tags::DEFAULT).Map();
  cvs.SetOwned(false);
  cvs.AddComponent("cell", AmanziMesh::CELL, 1);
  return Teuchos::rcp(new CompositeVector(cvs));
}

} // namespace Multiphase
} // namespace Amanzi

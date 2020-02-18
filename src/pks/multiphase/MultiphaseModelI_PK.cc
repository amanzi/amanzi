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

// Multiphase
#include "MoleFractionLiquid.hh"
#include "MultiphaseModelI_PK.hh"
#include "NCP_F.hh"
#include "NCP_MoleFractions.hh"
#include "ProductEvaluator.hh"
#include "SaturationGasEvaluator.hh"
#include "TotalComponentStorage.hh"
#include "TotalWaterStorage.hh"
#include "VaporPressureEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Standard constructor
****************************************************************** */
MultiphaseModelI_PK::MultiphaseModelI_PK(Teuchos::ParameterList& pk_tree,
                                         const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                         const Teuchos::RCP<State>& S,
                                         const Teuchos::RCP<TreeVector>& soln)
  : Multiphase_PK(pk_tree, glist, S, soln) {};

  
/* ******************************************************************
* Setup
****************************************************************** */
void MultiphaseModelI_PK::Setup(const Teuchos::Ptr<State>& S)
{
  Multiphase_PK::Setup(S);

  advection_water_key_ = Keys::getKey(domain_, "advection_water"); 
  pressure_vapor_key_ = Keys::getKey(domain_, "pressure_vapor"); 
  x_vapor_key_ = Keys::getKey(domain_, "mole_fraction_vapor"); 

  molecular_diff_liquid_key_ = Keys::getKey(domain_, "molecular_diff_liquid"); 
  molecular_diff_gas_key_ = Keys::getKey(domain_, "molecular_diff_gas"); 

  diffusion_liquid_key_ = Keys::getKey(domain_, "diffusion_liquid"); 
  diffusion_vapor_key_ = Keys::getKey(domain_, "diffusion_vapor"); 

  // gas mole fraction is the primary solution
  if (!S->HasField(x_gas_key_)) {
    component_names_ = mp_list_->sublist("molecular diffusion")
       .get<Teuchos::Array<std::string> >("aqueous names").toVector();
    std::vector<std::vector<std::string> > subfield_names(1);
    subfield_names[0] = component_names_;

    S->RequireField(x_gas_key_, passwd_, subfield_names)
      ->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, component_names_.size());

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", x_gas_key_);
    x_gas_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S->SetFieldEvaluator(x_gas_key_, x_gas_eval_);
  }

  // conserved quantities
  // -- total water storage
  if (!S->HasField(tws_key_)) {
    S->RequireField(tws_key_, tws_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", tws_key_)
         .set<std::string>("molar density liquid key", molar_density_liquid_key_)
         .set<std::string>("molar density gas key", molar_density_gas_key_)
         .set<std::string>("porosity key", porosity_key_)
         .set<std::string>("saturation liquid key", saturation_liquid_key_)
         .set<std::string>("mole fraction vapor key", x_vapor_key_);

    eval_tws_ = Teuchos::rcp(new TotalWaterStorage(elist));
    S->SetFieldEvaluator(tws_key_, eval_tws_);
  }

  // -- total component storage (for one component)
  if (!S->HasField(tcs_key_)) {
    S->RequireField(tcs_key_, tcs_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", tcs_key_)
         .set<std::string>("saturation liquid key", saturation_liquid_key_)
         .set<std::string>("porosity key", porosity_key_);

    eval_tcs_ = Teuchos::rcp(new TotalComponentStorage(elist));
    S->SetFieldEvaluator(tcs_key_, eval_tcs_);
  }

  // saturation
  if (!S->HasField(saturation_gas_key_)) {
    S->RequireField(saturation_gas_key_, saturation_gas_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", saturation_gas_key_)
         .set<std::string>("saturation liquid key", saturation_liquid_key_);
    auto eval = Teuchos::rcp(new SaturationGasEvaluator(elist));
    S->SetFieldEvaluator(saturation_gas_key_, eval);
  }

  // water vapor
  // -- vapor pressure
  if (!S->HasField(pressure_vapor_key_)) {
    S->RequireField(pressure_vapor_key_, pressure_vapor_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", pressure_vapor_key_)
         .set<std::string>("temperature key", temperature_key_)
         .set<std::string>("molar density liquid key", molar_density_liquid_key_)
         .set<std::string>("saturation liquid key", saturation_liquid_key_)
         .set<std::string>("vapor pressure model type", "water vapor over water/ice");
    auto eval = Teuchos::rcp(new VaporPressureEvaluator(elist, wrm_));
    S->SetFieldEvaluator(pressure_vapor_key_, eval);
  }

  // -- coefficient for diffusion operator in liquid phase
  if (!S->HasField(diffusion_vapor_key_)) {
    S->RequireField(diffusion_vapor_key_, diffusion_vapor_key_)
      ->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::Array<int> dep_powers(4, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(molar_density_gas_key_);
    dep_names.push_back(molecular_diff_gas_key_);
    dep_names.push_back(porosity_key_);
    dep_names.push_back(saturation_gas_key_);

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", diffusion_vapor_key_)
         .set<Teuchos::Array<std::string> >("evaluator dependencies", dep_names) 
         .set<Teuchos::Array<int> >("powers", dep_powers);

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S->SetFieldEvaluator(diffusion_vapor_key_, eval);
  }

  // liquid mole fraction (for current component)
  if (!S->HasField(x_liquid_key_)) {
    S->RequireField(x_liquid_key_, x_liquid_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", x_liquid_key_)
         .set<std::string>("pressure gas key", pressure_gas_key_)
         .set<std::string>("mole fraction gas key", x_gas_key_);
    auto eval = Teuchos::rcp(new MoleFractionLiquid(elist));
    S->SetFieldEvaluator(x_liquid_key_, eval);
  }

  // product evaluators
  // -- coefficient for advection operator (div eta_l q_l) in liquid phase 
  if (!S->HasField(advection_liquid_key_)) {
    S->RequireField(advection_liquid_key_, advection_liquid_key_)
      ->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::Array<int> dep_powers(4, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(molar_density_liquid_key_);
    dep_names.push_back(x_liquid_key_);
    dep_names.push_back(relperm_liquid_key_);
    dep_names.push_back(viscosity_liquid_key_);
    dep_powers[3] = -1;

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", advection_liquid_key_)
         .set<Teuchos::Array<std::string> >("evaluator dependencies", dep_names) 
         .set<Teuchos::Array<int> >("powers", dep_powers);
    // elist.sublist("verbose object").set<std::string>("verbosity level", "extreme");

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S->SetFieldEvaluator(advection_liquid_key_, eval);
  }

  if (!S->HasField(advection_water_key_)) {
    S->RequireField(advection_water_key_, advection_water_key_)
      ->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::Array<int> dep_powers(3, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(molar_density_liquid_key_);
    dep_names.push_back(relperm_liquid_key_);
    dep_names.push_back(viscosity_liquid_key_);
    dep_powers[2] = -1;

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", advection_water_key_)
         .set<Teuchos::Array<std::string> >("evaluator dependencies", dep_names) 
         .set<Teuchos::Array<int> >("powers", dep_powers);

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S->SetFieldEvaluator(advection_water_key_, eval);
  }

  // -- coefficient for advection operator in gas phase (4 fields)
  if (!S->HasField(advection_gas_key_)) {
    S->RequireField(advection_gas_key_, advection_gas_key_)
      ->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::Array<int> dep_powers(4, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(molar_density_gas_key_);
    dep_names.push_back(x_gas_key_);
    dep_names.push_back(relperm_gas_key_);
    dep_names.push_back(viscosity_gas_key_);
    dep_powers[3] = -1;

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", advection_gas_key_)
         .set<Teuchos::Array<std::string> >("evaluator dependencies", dep_names) 
         .set<Teuchos::Array<int> >("powers", dep_powers);

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S->SetFieldEvaluator(advection_gas_key_, eval);
  }

  // -- coefficient for diffusion operator in liquid phase
  if (!S->HasField(diffusion_liquid_key_)) {
    S->RequireField(diffusion_liquid_key_, diffusion_liquid_key_)
      ->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::Array<int> dep_powers(2, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(molecular_diff_liquid_key_);
    dep_names.push_back(molar_density_liquid_key_);

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", diffusion_liquid_key_)
         .set<Teuchos::Array<std::string> >("evaluator dependencies", dep_names) 
         .set<Teuchos::Array<int> >("powers", dep_powers);

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S->SetFieldEvaluator(diffusion_liquid_key_, eval);
  }

  // mole fraction of water vapor
  if (!S->HasField(x_vapor_key_)) {
    S->RequireField(x_vapor_key_, x_vapor_key_)
      ->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::Array<int> dep_powers(2, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(pressure_gas_key_);
    dep_names.push_back(pressure_vapor_key_);
    dep_powers[0] = -1;

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", x_vapor_key_)
         .set<Teuchos::Array<std::string> >("evaluator dependencies", dep_names) 
         .set<Teuchos::Array<int> >("powers", dep_powers);

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S->SetFieldEvaluator(x_vapor_key_, eval);
  }

  // nonlinear complimentary problem
  if (!S->HasField(ncp_f_key_)) {
    S->RequireField(ncp_f_key_, ncp_f_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", ncp_f_key_)
         .set<std::string>("saturation liquid key", saturation_liquid_key_);
    auto eval = Teuchos::rcp(new NCP_F(elist));
    S->SetFieldEvaluator(ncp_f_key_, eval);
  }

  if (!S->HasField(ncp_g_key_)) {
    S->RequireField(ncp_g_key_, ncp_g_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", ncp_g_key_)
         .set<std::string>("mole fraction vapor key", x_vapor_key_)
         .set<std::string>("mole fraction gas key", x_gas_key_);
    auto eval = Teuchos::rcp(new NCP_MoleFractions(elist));
    S->SetFieldEvaluator(ncp_g_key_, eval);
  }
}


/* ******************************************************************
* Push data to the state
****************************************************************** */
void MultiphaseModelI_PK::CommitStep(
    double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  Multiphase_PK::CommitStep(t_old, t_new, S);

  // miscalleneous fields
  S_->GetFieldEvaluator(ncp_fg_key_)->HasFieldChanged(S_.ptr(), passwd_);
}


/********************************************************************
* Modifies nonlinear update du using .. TBW
****************************************************************** */
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
MultiphaseModelI_PK::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                                      Teuchos::RCP<const TreeVector> u,
                                      Teuchos::RCP<TreeVector> du)
{
  // clip mole fraction to range [0; 1]
  const auto& u1c = *u->SubVector(1)->Data()->ViewComponent("cell");
  auto& du1c = *du->SubVector(1)->Data()->ViewComponent("cell");

  for (int i = 0; i < u1c.NumVectors(); ++i) {
    for (int c = 0; c < ncells_owned_; ++c) {
      // du1c[i][c] = std::min(du1c[i][c], u1c[i][c]);
      // du1c[i][c] = std::max(du1c[i][c], u1c[i][c] - 1.0);
    }    
  }

  // clip saturation (residual saturation is missing, FIXME)
  const auto& u2c = *u->SubVector(2)->Data()->ViewComponent("cell");
  auto& du2c = *du->SubVector(2)->Data()->ViewComponent("cell");

  for (int c = 0; c < ncells_owned_; ++c) {
    // du2c[0][c] = std::min(du2c[0][c], u2c[0][c]);
    // du2c[0][c] = std::max(du2c[0][c], u2c[0][c] - 1.0);
  }

  return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED;
}


/* ******************************************************************* 
* Create vector of solutions
******************************************************************* */
void MultiphaseModelI_PK::InitMPSolutionVector()
{
  soln_names_.push_back(pressure_liquid_key_);
  soln_names_.push_back(x_gas_key_);
  soln_names_.push_back(saturation_liquid_key_);

  for (int i = 0; i < soln_names_.size(); ++i) {
    auto field = Teuchos::rcp(new TreeVector());
    soln_->PushBack(field);
    field->SetData(S_->GetFieldData(soln_names_[i], passwd_));
  }
}


/* ******************************************************************* 
* Populate structure of equations with names of evqluators.
******************************************************************* */
void MultiphaseModelI_PK::InitMPPreconditioner()
{
  eqns_.resize(num_primary_ + 2);

  eqns_[0].adv_factors.resize(2, 1.0);
  eqns_[0].advection.push_back(std::make_pair(advection_water_key_, pressure_liquid_key_));
  eqns_[0].advection.push_back(std::make_pair("", ""));  // no gas phase

  eqns_[0].diff_factors.resize(2, 1.0);
  eqns_[0].diffusion.push_back(std::make_pair(diffusion_vapor_key_, x_vapor_key_));
  eqns_[0].diffusion.push_back(std::make_pair("", ""));  // no gas phase flux

  eqns_[0].storage = tws_key_;

  for (int i = 0; i < num_primary_; ++i) {
    int n = i + 1;
    eqns_[n].adv_factors.resize(2, 1.0);
    eqns_[n].advection.push_back(std::make_pair(advection_liquid_key_, pressure_liquid_key_));
    eqns_[n].advection.push_back(std::make_pair(advection_gas_key_, pressure_gas_key_));

    eqns_[n].diff_factors.resize(2, 1.0);
    eqns_[n].diffusion.push_back(std::make_pair(diffusion_liquid_key_, x_liquid_key_));
    eqns_[n].diffusion.push_back(std::make_pair(diffusion_gas_key_, x_gas_key_));

    eqns_[n].storage = tcs_key_;
  }

  // the last equaiton is special
  int n = num_primary_ + 1;
  eqns_[n].constraint = std::make_pair(ncp_f_key_, ncp_g_key_);
}


/* ******************************************************************* 
* Populate boundary conditions for various bc types
******************************************************************* */
void MultiphaseModelI_PK::PopulateBCs(int icomp, bool flag)
{
  // calculus of variations produces zero for fixed BCs
  double factor = (flag) ? 1.0 : 0.0;

  auto& bc_model_p = op_bcs_[0]->bc_model();
  auto& bc_value_p = op_bcs_[0]->bc_value();

  auto& bc_model_x = op_bcs_[1]->bc_model();
  auto& bc_value_x = op_bcs_[1]->bc_value();

  auto& bc_model_s = op_bcs_[2]->bc_model();
  auto& bc_value_s = op_bcs_[2]->bc_value();

  // initialize to zero
  int nfaces = bc_model_p.size();
  for (int n = 0; n < nfaces; ++n) {
    bc_model_p[n] = Operators::OPERATOR_BC_NONE;
    bc_value_p[n] = 0.0;

    bc_model_x[n] = Operators::OPERATOR_BC_NONE;
    bc_value_x[n] = 0.0;

    bc_model_s[n] = Operators::OPERATOR_BC_NONE;
    bc_value_s[n] = 0.0;
  }

  // populate boundary conditions
  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->bc_name() == "pressure") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        bc_model_p[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value_p[f] = it->second[0];
      }
    }

    if (bcs_[i]->bc_name() == "flux") {
      if (bcs_[i]->component_name() == "water") { 
        for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
          int f = it->first;
          bc_model_p[f] = Operators::OPERATOR_BC_NEUMANN;
          bc_value_p[f] = it->second[0] * factor;
        }
      } else if (bcs_[i]->component_id() == icomp) {
        for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
          int f = it->first;
          bc_model_x[f] = Operators::OPERATOR_BC_NEUMANN;
          bc_value_x[f] = it->second[0] * factor;
        }
      }
    }

    if (bcs_[i]->bc_name() == "saturation") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        bc_model_s[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value_s[f] = it->second[0];

        bc_model_x[f] = Operators::OPERATOR_BC_DIRICHLET;  // a huck
        bc_value_x[f] = 0.0;
      }
    }
  }

  // mark missing boundary conditions as zero flux conditions 
  AmanziMesh::Entity_ID_List cells;
  missed_bc_faces_ = 0;

  for (int f = 0; f < nfaces_owned_; f++) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

    if (cells.size() == 1) {
      if (bc_model_p[f] == Operators::OPERATOR_BC_NONE) {
        bc_model_p[f] = Operators::OPERATOR_BC_NEUMANN;
        missed_bc_faces_++;
      }

      if (bc_model_x[f] == Operators::OPERATOR_BC_NONE)
        bc_model_x[f] = Operators::OPERATOR_BC_NEUMANN;

      if (bc_model_s[f] == Operators::OPERATOR_BC_NONE)
        bc_model_s[f] = Operators::OPERATOR_BC_NEUMANN;
    }
  }

  // boundary conditions for derived fields 
  auto& bc_model_pg = op_bc_pg_->bc_model();
  auto& bc_value_pg = op_bc_pg_->bc_value();

  bc_model_pg = bc_model_p;

  for (int f = 0; f != nfaces_owned_; ++f) {
    if (bc_model_pg[f] == Operators::OPERATOR_BC_DIRICHLET) {
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int c = cells[0];

      bc_value_pg[f] = bc_value_p[f] + wrm_->second[(*wrm_->first)[c]]->capillaryPressure(bc_value_s[f]);
    }
    else if (bc_model_pg[f] == Operators::OPERATOR_BC_NEUMANN) {
      bc_value_pg[f] = 0.0;
    }
  }
}


/* ******************************************************************* 
* Map for indices
******************************************************************* */
SolutionStructure MultiphaseModelI_PK::EquationToSolution(int neqn)
{
  SolutionStructure soln(neqn, 0, -1);
  if (neqn == 0) return soln;

  if (neqn == num_primary_ + 1) {
    soln.var = 2;
    return soln;
  }
 
  soln.var = 1;
  soln.comp = neqn - 1;
  soln.matching_eqn = neqn;  // dEval_dSoln=0 for all eqns except this
  return soln;
}


/* ******************************************************************* 
* Tweak evaluators.
******************************************************************* */
void MultiphaseModelI_PK::ModifyEvaluators(int neqn)
{
  if (neqn > 0) {
    int ifield(0), n(neqn - 1);
    // ifield is dummy here
    Teuchos::rcp_dynamic_cast<TotalComponentStorage>(eval_tcs_)->set_subvector(ifield, n, kH_[n]);
    Teuchos::rcp_dynamic_cast<TotalComponentStorage>(eval_tcs_)->HasFieldChanged(S_.ptr(), passwd_, true);

    auto eval = S_->GetFieldEvaluator(x_liquid_key_);
    Teuchos::rcp_dynamic_cast<MoleFractionLiquid>(eval)->set_subvector(ifield, n, kH_[n]);
    Teuchos::rcp_dynamic_cast<MoleFractionLiquid>(eval)->HasFieldChanged(S_.ptr(), passwd_, true);
  }
}

}  // namespace Multiphase
}  // namespace Amanzi

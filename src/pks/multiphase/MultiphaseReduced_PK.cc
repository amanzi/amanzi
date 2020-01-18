/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Reduced multiphase model: H2O is only in liquid phase.
*/


// TPLs
#include "Teuchos_RCP.hpp"

// Multiphase
#include "MultiphaseReduced_PK.hh"
#include "ProductEvaluator.hh"
#include "TotalComponentStorage.hh"
#include "TotalWaterStorage.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Standard constructor
****************************************************************** */
MultiphaseReduced_PK::MultiphaseReduced_PK(Teuchos::ParameterList& pk_tree,
                                           const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                           const Teuchos::RCP<State>& S,
                                           const Teuchos::RCP<TreeVector>& soln) :
    Multiphase_PK(pk_tree, glist, S, soln) {};

  
/* ******************************************************************
* Setup
****************************************************************** */
void MultiphaseReduced_PK::Setup(const Teuchos::Ptr<State>& S)
{
  Multiphase_PK::Setup(S);

  // conserved quantities
  // -- total water storage
  S->RequireField("total_water_storage", "total_water_storage")->SetMesh(mesh_)->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::CELL, 1);

  Teuchos::ParameterList plist;
  plist.set<std::string>("my key", tws_key_)
       .set<std::string>("saturation liquid key", saturation_liquid_key_)
       .set<std::string>("porosity key", porosity_key_);

  eval_tws_ = Teuchos::rcp(new TotalWaterStorage(plist));
  S->SetFieldEvaluator(tws_key_, eval_tws_);

  // -- total component storage is allocated for component set with a modifier function
  S->RequireField(tcs_key_, tcs_key_)->SetMesh(mesh_)->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::CELL, 1);

  plist.set<std::string>("my key", tcs_key_);

  eval_tcs_ = Teuchos::rcp(new TotalComponentStorage(plist));
  S->SetFieldEvaluator(tcs_key_, eval_tcs_);

  // product evaluators
  {
    S->RequireField(advection_liquid_key_, advection_liquid_key_)
      ->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::Array<int> dep_factors(4, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(molar_density_liquid_key_);
    dep_names.push_back(relperm_liquid_key_);
    dep_names.push_back(viscosity_liquid_key_);
    dep_names.push_back(molar_fraction_liquid_key_);
    dep_factors[2] = -1;

    elist.set<std::string>("my key", advection_liquid_key_)
         .set<Teuchos::Array<std::string> >("evaluator dependencies", dep_names) 
         .set<Teuchos::Array<int> >("factors", dep_factors);

    auto eval = Teuchos::rcp(new MolarMobilityEvaluator(elist));
    S->SetFieldEvaluator(advection_liquid_key_, eval);
  }

  // molar mobility of gas phase
  {
    S->RequireField(advection_gas_key_, advection_gas_key_)
      ->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::Array<int> dep_factors(3, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(molar_density_gas_key_);
    dep_names.push_back(relperm_gas_key_);
    dep_names.push_back(viscosity_gas_key_);
    dep_names.push_back(molar_fraction_gas_key_);
    dep_factors[2] = -1;

    elist.set<std::string>("my key", advection_gas_key_)
         .set<Teuchos::Array<std::string> >("evaluator dependencies", dep_names) 
         .set<Teuchos::Array<int> >("factors", dep_factors);

    auto eval = Teuchos::rcp(new MolarMobilityEvaluator(elist));
    S->SetFieldEvaluator(advection_gas_key_, eval);
  }
}


/* ******************************************************************* 
* Create vector of solutions
******************************************************************* */
void MultiphaseReduced_PK::InitMPSolutionVector()
{
  soln_names_.push_back(pressure_liquid_key_);
  soln_names_.push_back(saturation_liquid_key_);
  soln_names_.push_back(x_liquid_key_);

  for (int i = 0; i < soln_names_.size(); ++i) {
    auto field = Teuchos::rcp(new TreeVector());
    soln_->PushBack(field);
    field->SetData(S_->GetFieldData(soln_names_[i], passwd_));
  }
}


/* ******************************************************************* 
* Create matrix structure of operators
******************************************************************* */
void MultiphaseReduced_PK::InitMPPreconditioner()
{
  eval_eqns_.push_back({advection_liquid_key_, "", "", " ", tws_key_});
  for (int i = 0; i < num_primary_; ++i) {
    eval_eqns_.push_back({advection_liquid_key_, advection_gas_key_, "any", "any", twc_key_});
  }
}


/* ******************************************************************* 
* Populate boundary conditions for various bc types
******************************************************************* */
void MultiphaseReduced_PK::PopulateBCs(int icomp)
{
  auto& bc_model_p = op_bcs_[0]->bc_model();
  auto& bc_value_p = op_bcs_[0]->bc_value();

  auto& bc_model_s = op_bcs_[1]->bc_model();
  auto& bc_value_s = op_bcs_[1]->bc_value();

  auto& bc_model_x = op_bcs_[2]->bc_model();
  auto& bc_value_x = op_bcs_[2]->bc_value();

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
          bc_value_p[f] = it->second[0] * (eta_l_ / rho_l_);
        }
      } else if (bcs_[i]->component_id() == icomp) {
        if (icomp > 0) {
          for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
            int f = it->first;
            bc_model_x[f] = Operators::OPERATOR_BC_NEUMANN;
            bc_value_x[f] = it->second[0] * (eta_l_ / rho_l_);
          }
        } else {
          for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
            int f = it->first;
            bc_model_s[f] = Operators::OPERATOR_BC_NEUMANN;
            bc_value_s[f] = it->second[0] * (eta_l_ / rho_l_);
          }
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

      if (bc_model_s[f] == Operators::OPERATOR_BC_NONE)
        bc_model_s[f] = Operators::OPERATOR_BC_NEUMANN;

      if (bc_model_x[f] == Operators::OPERATOR_BC_NONE)
        bc_model_x[f] = Operators::OPERATOR_BC_NEUMANN;
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
std::pair<Key, int> MultiphaseReduced_PK::EquationToSolution(int neqn) {
  if (neqn == 0) return std::make_pair(pressure_liquid_key_, 0);
  if (neqn == 1) return std::make_pair(saturation_liquid_key_, 0);
  return std::make_pair(x_liquid_key, neqn - 2);
}


/* ******************************************************************* 
* Remove water molar fraction (= 1) from this evalautor
******************************************************************* */
void MultiphaseReduced_PK::ModifyEvaluator(
  int neqn, int pos, const Teuchos::RCP<MultiphaseBaseEvaluator>& eval)
{
  if (neqn == 0) {
    eval->set_subvector(2, -1);
    if (pos == 4)
      Teuchos::rcp_dynamic_cast<TotalComponentStorage>(eval)->set_kH(kH_[neqn - 1]);
  }
}

}  // namespace Multiphase
}  // namespace Amanzi

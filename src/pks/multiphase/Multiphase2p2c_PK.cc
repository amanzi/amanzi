/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Reduced multiphase model for water and hydrogen.
*/

// TPLs
#include "Teuchos_RCP.hpp"

// Multiphase
#include "MoleFractionLiquid.hh"
#include "Multiphase2p2c_PK.hh"
#include "NCP_F.hh"
#include "NCP_HenryLaw.hh"
#include "ProductEvaluator.hh"
#include "TotalComponentStorageTest.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Standard constructor
****************************************************************** */
Multiphase2p2c_PK::Multiphase2p2c_PK(Teuchos::ParameterList& pk_tree,
                                     const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                     const Teuchos::RCP<State>& S,
                                     const Teuchos::RCP<TreeVector>& soln) :
    Multiphase_PK(pk_tree, glist, S, soln) {};

  
/* ******************************************************************
* Setup
****************************************************************** */
void Multiphase2p2c_PK::Setup(const Teuchos::Ptr<State>& S)
{
  Multiphase_PK::Setup(S);

  molar_density_water_key_ = Keys::getKey(domain_, "molar_density_water"); 
  advection_liquid_reduced_key_ = Keys::getKey(domain_, "advection_liquid_reduced"); 

  ncp_f_key_ = Keys::getKey(domain_, "ncp_f"); 
  ncp_g_key_ = Keys::getKey(domain_, "ncp_g"); 
  ncp_fg_key_ = Keys::getKey(domain_, "ncp_fg"); 

  // liquid molar fraction is the primary solution
  if (!S->HasField(molar_density_liquid_key_)) {
    component_names_ = mp_list_->sublist("molecular diffusion")
       .get<Teuchos::Array<std::string> >("primary component names").toVector();
    std::vector<std::vector<std::string> > subfield_names(1);
    subfield_names[0] = component_names_;

    S->RequireField(molar_density_liquid_key_, passwd_, subfield_names)
      ->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, component_names_.size());

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", molar_density_liquid_key_);
    auto eval = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S->SetFieldEvaluator(molar_density_liquid_key_, eval);
  }

  // conserved quantities
  // -- total water storage
  {
    S->RequireField(tws_key_, tws_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::Array<int> dep_powers(3, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(molar_density_water_key_);
    dep_names.push_back(porosity_key_);
    dep_names.push_back(saturation_liquid_key_);

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", tws_key_)
         .set<Teuchos::Array<std::string> >("evaluator dependencies", dep_names) 
         .set<Teuchos::Array<int> >("powers", dep_powers);

    eval_tws_ = Teuchos::rcp(new ProductEvaluator(elist));
    S->SetFieldEvaluator(tws_key_, eval_tws_);
  }

  // -- total component storage
  {
    S->RequireField(tcs_key_, tcs_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", tcs_key_)
         .set<std::string>("saturation liquid key", saturation_liquid_key_)
         .set<std::string>("porosity key", porosity_key_);

    eval_tcs_ = Teuchos::rcp(new TotalComponentStorageTest(elist));
    S->SetFieldEvaluator(tcs_key_, eval_tcs_);
  }

  // product evaluators
  // -- coefficient for advection operator (div eta_l q_l) in liquid phase
  {
    S->RequireField(advection_liquid_key_, advection_liquid_key_)
      ->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::Array<int> dep_powers(3, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(molar_density_liquid_key_);
    dep_names.push_back(relperm_liquid_key_);
    dep_names.push_back(viscosity_liquid_key_);
    dep_powers[2] = -1;

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", advection_liquid_key_)
         .set<Teuchos::Array<std::string> >("evaluator dependencies", dep_names) 
         .set<Teuchos::Array<int> >("powers", dep_powers);

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S->SetFieldEvaluator(advection_liquid_key_, eval);
  }

  {
    S->RequireField(advection_liquid_reduced_key_, advection_liquid_reduced_key_)
      ->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::Array<int> dep_powers(3, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(molar_density_water_key_);
    dep_names.push_back(relperm_liquid_key_);
    dep_names.push_back(viscosity_liquid_key_);
    dep_powers[2] = -1;

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", advection_liquid_reduced_key_)
         .set<Teuchos::Array<std::string> >("evaluator dependencies", dep_names) 
         .set<Teuchos::Array<int> >("powers", dep_powers);

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S->SetFieldEvaluator(advection_liquid_reduced_key_, eval);
  }

  // -- coefficient for advection operator in gas phase
  {
    S->RequireField(advection_gas_key_, advection_gas_key_)
      ->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::Array<int> dep_powers(3, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(molar_density_gas_key_);
    dep_names.push_back(relperm_gas_key_);
    dep_names.push_back(viscosity_gas_key_);
    dep_powers[2] = -1;

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", advection_gas_key_)
         .set<Teuchos::Array<std::string> >("evaluator dependencies", dep_names) 
         .set<Teuchos::Array<int> >("powers", dep_powers);

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S->SetFieldEvaluator(advection_gas_key_, eval);
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
         .set<std::string>("pressure gas key", pressure_gas_key_)
         .set<std::string>("molar density liquid key", molar_density_liquid_key_);
    auto eval = Teuchos::rcp(new NCP_HenryLaw(elist));
    S->SetFieldEvaluator(ncp_g_key_, eval);
  }

  if (!S->HasField(ncp_fg_key_)) {
    S->RequireField(ncp_fg_key_, ncp_fg_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::Array<int> dep_powers(2, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(ncp_f_key_);
    dep_names.push_back(ncp_g_key_);

    Teuchos::ParameterList elist;
    elist.set<std::string>("my key", ncp_fg_key_)
         .set<Teuchos::Array<std::string> >("evaluator dependencies", dep_names) 
         .set<Teuchos::Array<int> >("powers", dep_powers);

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S->SetFieldEvaluator(ncp_fg_key_, eval);
  }
}


/* ******************************************************************
* Push data to the state
****************************************************************** */
void Multiphase2p2c_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  Multiphase_PK::CommitStep(t_old, t_new, S);

  // miscalleneous fields
  S_->GetFieldEvaluator(ncp_fg_key_)->HasFieldChanged(S_.ptr(), passwd_);
}


/* ******************************************************************* 
* Create vector of solutions
******************************************************************* */
void Multiphase2p2c_PK::InitMPSolutionVector()
{
  soln_names_.push_back(pressure_liquid_key_);
  soln_names_.push_back(saturation_liquid_key_);
  soln_names_.push_back(molar_density_liquid_key_);

  for (int i = 0; i < soln_names_.size(); ++i) {
    auto field = Teuchos::rcp(new TreeVector());
    soln_->PushBack(field);
    field->SetData(S_->GetFieldData(soln_names_[i], passwd_));
  }

  // create names used by molecular diffusion
  varx_name_ = {molar_density_liquid_key_, molar_density_gas_key_};
}


/* ******************************************************************* 
* Create matrix structure of pairs of evaluators and scalar factors
******************************************************************* */
void Multiphase2p2c_PK::InitMPPreconditioner()
{
  eval_eqns_.push_back({{advection_liquid_reduced_key_, 1.0}, {"", 0.0},
                        {"any", -1.0}, {"any", -1.0}, {tws_key_, 1.0}});

  for (int i = 0; i < num_primary_; ++i) {
    eval_eqns_.push_back({{advection_liquid_key_, 1.0}, {advection_gas_key_, 1.0},
                          {"any", 1.0}, {"any", 1.0}, {tcs_key_, 1.0}});
  }

  eval_eqns_.push_back({{ncp_f_key_, 1.0}, {ncp_g_key_, 1.0}});
}


/* ******************************************************************* 
* Populate boundary conditions for various bc types
******************************************************************* */
void Multiphase2p2c_PK::PopulateBCs(int icomp, bool flag)
{
  // calculus of variations produces zero for fixed BCs
  double factor = (flag) ? 1.0 : 0.0;

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
          bc_value_p[f] = it->second[0] * factor;
        }
      } else if (bcs_[i]->component_id() == icomp) {
        if (icomp > 0) {
          for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
            int f = it->first;
            bc_model_x[f] = Operators::OPERATOR_BC_NEUMANN;
            bc_value_x[f] = it->second[0] * factor;
          }
        } else {
          for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
            int f = it->first;
            bc_model_s[f] = Operators::OPERATOR_BC_NEUMANN;
            bc_value_s[f] = it->second[0] * factor;
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
std::pair<int, int> Multiphase2p2c_PK::EquationToSolution(int neqn) {
  if (neqn == 0) return std::make_pair(0, 0);
  if (neqn == 1) return std::make_pair(1, 0);
  return std::make_pair(2, neqn - 2);
}


/* ******************************************************************* 
* Tweak evaluators.
******************************************************************* */
void Multiphase2p2c_PK::ModifyEvaluators(int neqn)
{
  if (neqn > 0) {
    int ifield(0), n(neqn - 1);
    // ifield is dummy here
    Teuchos::rcp_dynamic_cast<TotalComponentStorageTest>(eval_tcs_)->set_subvector(ifield, n, kH_[n]);

    // mole fraction is second in the dependencies set
    auto eval = S_->GetFieldEvaluator(advection_liquid_key_);
    Teuchos::rcp_dynamic_cast<ProductEvaluator>(eval)->set_subvector(1, n, kH_[n]);

    eval = S_->GetFieldEvaluator(ncp_g_key_);
    Teuchos::rcp_dynamic_cast<NCP_HenryLaw>(eval)->set_subvector(ifield, n, kH_[n]);
  }
}

}  // namespace Multiphase
}  // namespace Amanzi

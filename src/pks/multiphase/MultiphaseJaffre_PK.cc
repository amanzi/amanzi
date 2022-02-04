/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Two-phase two-component system with water and hydrogen. 
  It is used currently for verification purpose.
*/

// TPLs
#include "Teuchos_RCP.hpp"

// Multiphase
#include "MultiphaseJaffre_PK.hh"
#include "NCP_F.hh"
#include "NCP_HenryLaw.hh"
#include "ProductEvaluator.hh"
#include "TotalComponentStorage_MolarDensity.hh"

namespace Amanzi {
namespace Multiphase {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* Standard constructor
****************************************************************** */
MultiphaseJaffre_PK::MultiphaseJaffre_PK(
    Teuchos::ParameterList& pk_tree,
    const Teuchos::RCP<Teuchos::ParameterList>& glist,
    const Teuchos::RCP<State>& S,
    const Teuchos::RCP<TreeVector>& soln)
  : Multiphase_PK(pk_tree, glist, S, soln) {};

  
/* ******************************************************************
* Setup
****************************************************************** */
void MultiphaseJaffre_PK::Setup()
{
  Multiphase_PK::Setup();

  advection_water_key_ = Keys::getKey(domain_, "advection_water"); 
  molar_density_water_key_ = Keys::getKey(domain_, "molar_density_water"); 

  diffusion_liquid_key_ = Keys::getKey(domain_, "diffusion_liquid"); 

  molecular_diff_liquid_key_ = Keys::getKey(domain_, "molecular_diff_liquid"); 
  molecular_diff_gas_key_ = Keys::getKey(domain_, "molecular_diff_gas"); 

  // liquid molar fraction is the primary solution
  if (!S_->HasRecord(molar_density_liquid_key_)) {
    component_names_ = mp_list_->sublist("molecular diffusion")
       .get<Teuchos::Array<std::string> >("aqueous names").toVector();

    S_->Require<CV_t, CVS_t>(molar_density_liquid_key_, Tags::DEFAULT, passwd_, component_names_)
      .SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, component_names_.size());

    Teuchos::ParameterList elist(molar_density_liquid_key_);
    elist.set<std::string>("name", molar_density_liquid_key_)
         .set<std::string>("tag", "");
    auto eval = Teuchos::rcp(new EvaluatorPrimary<CV_t, CVS_t>(elist));
    S_->SetEvaluator(molar_density_liquid_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(molar_density_liquid_key_, Tags::DEFAULT,
                                       molar_density_liquid_key_, Tags::DEFAULT, molar_density_liquid_key_);
  }

  // conserved quantities
  // -- total water storage
  if (!S_->HasRecord(tws_key_)) {
    S_->Require<CV_t, CVS_t>(tws_key_, Tags::DEFAULT, tws_key_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::Array<int> dep_powers(3, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(molar_density_water_key_);
    dep_names.push_back(porosity_key_);
    dep_names.push_back(saturation_liquid_key_);

    Teuchos::ParameterList elist(tws_key_);
    elist.set<std::string>("my key", tws_key_)
         .set<std::string>("tag", "")
         .set<Teuchos::Array<std::string> >("dependencies", dep_names) 
         .set<bool>("dependency tags are my tag", true)
         .set<Teuchos::Array<int> >("powers", dep_powers);

    eval_tws_ = Teuchos::rcp(new ProductEvaluator(elist));
    S_->SetEvaluator(tws_key_, Tags::DEFAULT, eval_tws_);

    S_->RequireDerivative<CV_t, CVS_t>(tws_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, tws_key_);
  }

  // -- total component storage
  if (!S_->HasRecord(tcs_key_)) {
    S_->Require<CV_t, CVS_t>(tcs_key_, Tags::DEFAULT, tcs_key_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist(tcs_key_);
    elist.set<std::string>("my key", tcs_key_)
         .set<std::string>("tag", "")
         .set<std::string>("saturation liquid key", saturation_liquid_key_)
         .set<std::string>("porosity key", porosity_key_)
         .set<std::string>("molar density liquid key", molar_density_liquid_key_)
         .set<std::string>("molar density gas key", molar_density_gas_key_);

    eval_tcs_ = Teuchos::rcp(new TotalComponentStorage_MolarDensity(elist));
    S_->SetEvaluator(tcs_key_, Tags::DEFAULT, eval_tcs_);

    S_->RequireDerivative<CV_t, CVS_t>(tcs_key_, Tags::DEFAULT,
                                       pressure_liquid_key_, Tags::DEFAULT, tcs_key_);
    S_->RequireDerivative<CV_t, CVS_t>(tcs_key_, Tags::DEFAULT,
                                       molar_density_liquid_key_, Tags::DEFAULT, tcs_key_);
    S_->RequireDerivative<CV_t, CVS_t>(tcs_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, tcs_key_);
  }

  // product evaluators
  // -- coefficient for advection operator (div eta_l q_l) in liquid phase
  if (!S_->HasRecord(advection_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(advection_liquid_key_, Tags::DEFAULT, advection_liquid_key_)
      .SetMesh(mesh_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::Array<int> dep_powers(3, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(molar_density_liquid_key_);
    dep_names.push_back(relperm_liquid_key_);
    dep_names.push_back(viscosity_liquid_key_);
    dep_powers[2] = -1;

    Teuchos::ParameterList elist(advection_liquid_key_);
    elist.set<std::string>("my key", advection_liquid_key_)
         .set<std::string>("tag", "")
         .set<Teuchos::Array<std::string> >("dependencies", dep_names) 
         .set<bool>("dependency tags are my tag", true)
         .set<Teuchos::Array<int> >("powers", dep_powers);

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S_->SetEvaluator(advection_liquid_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(advection_liquid_key_, Tags::DEFAULT,
                                       molar_density_liquid_key_, Tags::DEFAULT, advection_liquid_key_);
    S_->RequireDerivative<CV_t, CVS_t>(advection_liquid_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, advection_liquid_key_);
  }

  if (!S_->HasRecord(advection_water_key_)) {
    S_->Require<CV_t, CVS_t>(advection_water_key_, Tags::DEFAULT, advection_water_key_)
      .SetMesh(mesh_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::Array<int> dep_powers(3, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(molar_density_water_key_);
    dep_names.push_back(relperm_liquid_key_);
    dep_names.push_back(viscosity_liquid_key_);
    dep_powers[2] = -1;

    Teuchos::ParameterList elist(advection_water_key_);
    elist.set<std::string>("my key", advection_water_key_)
         .set<std::string>("tag", "")
         .set<Teuchos::Array<std::string> >("dependencies", dep_names) 
         .set<bool>("dependency tags are my tag", true)
         .set<Teuchos::Array<int> >("powers", dep_powers);

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S_->SetEvaluator(advection_water_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(advection_water_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, advection_water_key_);
  }

  // -- coefficient for advection operator in gas phase
  if (!S_->HasRecord(advection_gas_key_)) {
    S_->Require<CV_t, CVS_t>(advection_gas_key_, Tags::DEFAULT, advection_gas_key_)
      .SetMesh(mesh_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::Array<int> dep_powers(3, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(molar_density_gas_key_);
    dep_names.push_back(relperm_gas_key_);
    dep_names.push_back(viscosity_gas_key_);
    dep_powers[2] = -1;

    Teuchos::ParameterList elist(advection_gas_key_);
    elist.set<std::string>("my key", advection_gas_key_)
         .set<std::string>("tag", "")
         .set<Teuchos::Array<std::string> >("dependencies", dep_names) 
         .set<bool>("dependency tags are my tag", true)
         .set<Teuchos::Array<int> >("powers", dep_powers);

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S_->SetEvaluator(advection_gas_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(advection_gas_key_, Tags::DEFAULT,
                                       pressure_liquid_key_, Tags::DEFAULT, advection_gas_key_);
    S_->RequireDerivative<CV_t, CVS_t>(advection_gas_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, advection_gas_key_);
  }

  // -- coefficient for diffusion operator in liquid phase
  if (!S_->HasRecord(diffusion_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(diffusion_liquid_key_, Tags::DEFAULT, diffusion_liquid_key_)
      .SetMesh(mesh_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::Array<int> dep_powers(3, 1);
    Teuchos::Array<std::string> dep_names;

    dep_names.push_back(molecular_diff_liquid_key_);
    dep_names.push_back(porosity_key_);
    dep_names.push_back(saturation_liquid_key_);

    Teuchos::ParameterList elist(diffusion_liquid_key_);
    elist.set<std::string>("my key", diffusion_liquid_key_)
         .set<std::string>("tag", "")
         .set<Teuchos::Array<std::string> >("dependencies", dep_names) 
         .set<bool>("dependency tags are my tag", true)
         .set<Teuchos::Array<int> >("powers", dep_powers);

    auto eval = Teuchos::rcp(new ProductEvaluator(elist));
    S_->SetEvaluator(diffusion_liquid_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(diffusion_liquid_key_, Tags::DEFAULT,
                                       molar_density_liquid_key_, Tags::DEFAULT, diffusion_liquid_key_);
    S_->RequireDerivative<CV_t, CVS_t>(diffusion_liquid_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, diffusion_liquid_key_);
  }

  // -- coefficient for diffusion operator in gas phase

  // nonlinear complimentary problem
  if (!S_->HasRecord(ncp_f_key_)) {
    S_->Require<CV_t, CVS_t>(ncp_f_key_, Tags::DEFAULT, ncp_f_key_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist(ncp_f_key_);
    elist.set<std::string>("my key", ncp_f_key_)
         .set<std::string>("saturation liquid key", saturation_liquid_key_)
         .set<std::string>("tag", "");
    auto eval = Teuchos::rcp(new NCP_F(elist));
    S_->SetEvaluator(ncp_f_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(ncp_f_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, ncp_f_key_);
  }

  if (!S_->HasRecord(ncp_g_key_)) {
    S_->Require<CV_t, CVS_t>(ncp_g_key_, Tags::DEFAULT, ncp_g_key_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist(ncp_g_key_);
    elist.set<std::string>("my key", ncp_g_key_)
         .set<std::string>("tag", "")
         .set<std::string>("pressure gas key", pressure_gas_key_)
         .set<std::string>("molar density liquid key", molar_density_liquid_key_);
    // elist.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    auto eval = Teuchos::rcp(new NCP_HenryLaw(elist));
    S_->SetEvaluator(ncp_g_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(ncp_g_key_, Tags::DEFAULT,
                                       pressure_liquid_key_, Tags::DEFAULT, ncp_g_key_);
    S_->RequireDerivative<CV_t, CVS_t>(ncp_g_key_, Tags::DEFAULT,
                                       molar_density_liquid_key_, Tags::DEFAULT, ncp_g_key_);
    S_->RequireDerivative<CV_t, CVS_t>(ncp_g_key_, Tags::DEFAULT,
                                       saturation_liquid_key_, Tags::DEFAULT, ncp_g_key_);
  }
}


/* ******************************************************************
* Push data to the state
****************************************************************** */
void MultiphaseJaffre_PK::CommitStep(
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
MultiphaseJaffre_PK::ModifyCorrection(
    double h, Teuchos::RCP<const TreeVector> res,
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
void MultiphaseJaffre_PK::InitMPSolutionVector()
{
  soln_names_.push_back(pressure_liquid_key_);
  soln_names_.push_back(molar_density_liquid_key_);
  soln_names_.push_back(saturation_liquid_key_);

  for (int i = 0; i < soln_names_.size(); ++i) {
    auto field = Teuchos::rcp(new TreeVector());
    soln_->PushBack(field);
    field->SetData(S_->GetPtrW<CompositeVector>(soln_names_[i], Tags::DEFAULT, passwd_));
  }
}


/* ******************************************************************* 
* Create matrix structure of pairs of evaluators and scalar factors
******************************************************************* */
void MultiphaseJaffre_PK::InitMPPreconditioner()
{
  eqns_.resize(num_primary_ + 2);

  eqns_[0].adv_factors.resize(2, 1.0);
  eqns_[0].advection.push_back(std::make_pair(advection_water_key_, pressure_liquid_key_));
  eqns_[0].advection.push_back(std::make_pair("", ""));  // no gas phase

  eqns_[0].diff_factors.resize(2, -mol_mass_[0] / mol_mass_H2O_);
  eqns_[0].diffusion.push_back(std::make_pair(diffusion_liquid_key_, molar_density_liquid_key_));
  eqns_[0].diffusion.push_back(std::make_pair("", ""));  // no gas phase flux

  eqns_[0].storage = tws_key_;

  for (int i = 0; i < num_primary_; ++i) {
    int n = i + 1;
    eqns_[n].adv_factors.resize(2, 1.0);
    eqns_[n].advection.push_back(std::make_pair(advection_liquid_key_, pressure_liquid_key_));
    eqns_[n].advection.push_back(std::make_pair(advection_gas_key_, pressure_gas_key_));

    eqns_[n].diff_factors.resize(2, 1.0);
    eqns_[n].diffusion.push_back(std::make_pair(diffusion_liquid_key_, molar_density_liquid_key_));
    eqns_[n].diffusion.push_back(std::make_pair("", ""));

    eqns_[n].storage = tcs_key_;
  }

  // the last equaiton is special
  int n = num_primary_ + 1;
  eqns_[n].constraint = std::make_pair(ncp_f_key_, ncp_g_key_);
}


/* ******************************************************************* 
* Populate boundary conditions for various bc types
******************************************************************* */
void MultiphaseJaffre_PK::PopulateBCs(int icomp, bool flag)
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

    bc_model_s[n] = Operators::OPERATOR_BC_NONE;
    bc_value_s[n] = 0.0;

    bc_model_x[n] = Operators::OPERATOR_BC_NONE;
    bc_value_x[n] = 0.0;
  }

  // populate boundary conditions
  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->get_bc_name() == "pressure") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        bc_model_p[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value_p[f] = it->second[0];
      }
    }

    if (bcs_[i]->get_bc_name() == "flux") {
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

    if (bcs_[i]->get_bc_name() == "saturation") {
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
SolutionStructure MultiphaseJaffre_PK::EquationToSolution(int neqn)
{
  SolutionStructure soln(neqn, 0, -1);
  return soln;
}


/* ******************************************************************* 
* Tweak evaluators.
******************************************************************* */
void MultiphaseJaffre_PK::ModifyEvaluators(int neqn)
{
  if (neqn > 0) {
    int ifield(0), n(neqn - 1);
    // ifield is dummy here
    Teuchos::rcp_dynamic_cast<TotalComponentStorage_MolarDensity>(eval_tcs_)->set_subvector(ifield, n, kH_[n]);
    Teuchos::rcp_dynamic_cast<TotalComponentStorage_MolarDensity>(eval_tcs_)->Update(*S_, passwd_, true);

    // mole fraction is second in the dependencies set
    auto eval = S_->GetEvaluatorPtr(advection_liquid_key_, Tags::DEFAULT);
    Teuchos::rcp_dynamic_cast<ProductEvaluator>(eval)->set_subvector(1, n, kH_[n]);
    Teuchos::rcp_dynamic_cast<ProductEvaluator>(eval)->Update(*S_, passwd_, true);

    eval = S_->GetEvaluatorPtr(ncp_g_key_, Tags::DEFAULT);
    Teuchos::rcp_dynamic_cast<NCP_HenryLaw>(eval)->set_subvector(ifield, n, kH_[n]);
    Teuchos::rcp_dynamic_cast<NCP_HenryLaw>(eval)->Update(*S_, passwd_, true);
  }
}

}  // namespace Multiphase
}  // namespace Amanzi

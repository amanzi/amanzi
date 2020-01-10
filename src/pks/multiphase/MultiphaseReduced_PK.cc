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
}


/* ******************************************************************* 
* Create vector of solutions
******************************************************************* */
void MultiphaseReduced_PK::InitMPSolutionVector()
{
  soln_names_.push_back(pressure_liquid_key_);
  soln_names_.push_back(x_liquid_key_);
  soln_names_.push_back(saturation_liquid_key_);

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
  eval_acc_.push_back("total_water_storage");
  eval_adv_.push_back("");
  eval_diff_.push_back("");

  for (int i = 0; i < num_primary_; ++i) {
    eval_acc_.push_back("total_component_storage");
  }
}


/* ******************************************************************* 
* Populate boundary conditions for various bc types
******************************************************************* */
void MultiphaseReduced_PK::PopulateBCs(int icomp)
{
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
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        bc_model_p[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value_p[f] = it->second[0];
      }
    }

    if (bcs_[i]->bc_name() == "saturation") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        bc_value_s[f] = it->second[0];
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
        bc_value_p[f] = 0.0;
        bc_value_s[f] = 0.0;
        bc_value_x[f] = 0.0;
      }

      missed_bc_faces_++;
    }
  }

  bc_model_s = bc_model_p;
  bc_model_x = bc_model_p;

  // boundary conditions for derived fields 
  auto& bc_model_pg = op_bc_pg_->bc_model();
  auto& bc_value_pg = op_bc_pg_->bc_value();

  for (int f = 0; f != nfaces_owned_; ++f) {
    if (bc_model_p[f] == Operators::OPERATOR_BC_DIRICHLET) {
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int c = cells[0];

      bc_value_pg[f] = bc_value_p[f] + wrm_->second[(*wrm_->first)[c]]->capillaryPressure(bc_value_s[f]);
    }
  }

  bc_model_pg = bc_model_p;
}

}  // namespace Multiphase
}  // namespace Amanzi

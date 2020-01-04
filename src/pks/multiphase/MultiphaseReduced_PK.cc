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

  
/* ******************************************************************* 
* Create vector of solutions
******************************************************************* */
void MultiphaseReduced_PK::InitializeSolution()
{
  soln_names_.push_back(pressure_liquid_key_);
  soln_names_.push_back(tcc_key_);
  soln_names_.push_back(saturation_liquid_key_);
  for (int i = 0; i < soln_names_.size(); ++i) {
    auto field = Teuchos::rcp(new TreeVector());
    soln_->PushBack(field);
    field->SetData(S_->GetFieldData(soln_names_[i], passwd_));
  }
}


/* ******************************************************************* 
* Populate boundary conditions for various bc types
******************************************************************* */
void MultiphaseReduced_PK::PopulateBCs()
{
  auto& bc_model_p = op_bcs_[0]->bc_model();
  auto& bc_value_p = op_bcs_[0]->bc_value();

  auto& bc_model_tcc = op_bcs_[1]->bc_model();
  auto& bc_value_tcc = op_bcs_[1]->bc_value_vector(num_primary_);

  auto& bc_model_s = op_bcs_[2]->bc_model();
  auto& bc_value_s = op_bcs_[2]->bc_value();

  // initialize to zero
  int nfaces = bc_model_p.size();
  for (int n = 0; n < nfaces; ++n) {
    bc_model_p[n] = Operators::OPERATOR_BC_NONE;
    bc_value_p[n] = 0.0;

    bc_model_tcc[n] = Operators::OPERATOR_BC_NONE;
    for (int i = 0; i < num_primary_; ++i) bc_value_tcc[n][i] = 0.0;

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

    if (bcs_[i]->bc_name() == "saturation") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        bc_value_s[f] = it->second[0];
      }
    }

    if (bcs_[i]->bc_name() == "hydrogen density") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        bc_value_tcc[f][0] = it->second[0];
      }
    }
  }

  // mark missing boundary conditions as zero flux conditions 
  AmanziMesh::Entity_ID_List cells;
  missed_bc_faces_ = 0;
  nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  for (int f = 0; f < nfaces; f++) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

    if (cells.size() == 1) {
      if (bc_model_p[f] == Operators::OPERATOR_BC_NONE) {
        bc_model_p[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value_p[f] = 0.0;
        bc_value_s[f] = 0.0;
        bc_value_tcc[f][0] = 0.0;
      }

      missed_bc_faces_++;
    }
  }

  bc_model_s = bc_model_p;
  bc_model_tcc = bc_model_p;
}

}  // namespace Multiphase
}  // namespace Amanzi

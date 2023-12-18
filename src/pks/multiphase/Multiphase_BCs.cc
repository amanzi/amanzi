/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  MultiPhase PK

  Support of boundary conditions. We may have multiple cases and the list
  below will grow in time. Let EQN(u, v) = 0 be an equation for u involving
  another primary unknown v and secondary fields w:

        \Sum_i div[ fi(u,v) \grad gi(u,v) ] = 0

  1. Total flux. Only one operator in this equation may impose this condition.
  The other operators must use zero flux boundary conditions.

  2. Dirichlet for u and v. We have to compute essential BCs for gi(u, v) using
  the dependency tree and copy data to corresponding boundary operators.
  Currently this generality is not supported, since it is not clear how to handle
  provided BCs for secondary fields. Hence, we assume that either gi(u,v) = u or
  gi(u,v) = v.

  3. Dirichlet for u and flux for v. We assume that fi \grad gi = flux u.
*/

// TPLs
#include "Teuchos_RCP.hpp"

// Amanzi
#include "Mesh_Algorithms.hh"

// Multiphase
#include "Multiphase_PK.hh"
#include "Multiphase_Utils.hh"

namespace Amanzi {
namespace Multiphase {

/* *******************************************************************
* Populate boundary conditions for various bc types
******************************************************************* */
void
Multiphase_PK::PopulateBCs(int icomp, bool flag)
{
  int n0 = (system_["energy eqn"]) ? 2 : 1;

  // initialize primary and secondary fields to no-BCs
  auto all_names = soln_names_;
  all_names.insert(all_names.end(), secondary_names_.begin(), secondary_names_.end());

  for (int i = 0; i < all_names.size(); ++i) {
    Key name = all_names[i];
    auto& bc_model = op_bcs_[name]->bc_model();
    auto& bc_value = op_bcs_[name]->bc_value();

    int nfaces = bc_model.size();
    for (int n = 0; n < nfaces; ++n) {
      bc_model[n] = Operators::OPERATOR_BC_NONE;
      bc_value[n] = 0.0;
    }
  }

  // calculus of variations produces zero for fixed BCs
  double factor = (flag) ? 1.0 : 0.0;

  // populate boundary conditions
  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->get_bc_name() == "pressure") {
      auto& bc_model = op_bcs_[pressure_liquid_key_]->bc_model();
      auto& bc_value = op_bcs_[pressure_liquid_key_]->bc_value();

      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value[f] = it->second[0];
      }
    }

    if (bcs_[i]->get_bc_name() == "flux") {
      if (bcs_[i]->component_name() == "water") {
        auto& bc_model = op_bcs_[pressure_liquid_key_]->bc_model();
        auto& bc_value = op_bcs_[pressure_liquid_key_]->bc_value();

        for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
          int f = it->first;
          bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
          bc_value[f] = it->second[0] * factor;
        }
      } else if (bcs_[i]->component_id() == icomp) {
        Key x_key_base = splitPhase(soln_names_[n0]).first;
        Key x_key = mergePhase(x_key_base, bcs_[i]->component_phase());

        auto& bc_model = op_bcs_[x_key]->bc_model();
        auto& bc_value = op_bcs_[x_key]->bc_value();

        for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
          int f = it->first;
          bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
          bc_value[f] = it->second[0] * factor;
        }
      }
    }

    if (bcs_[i]->get_bc_name() == "concentration" && bcs_[i]->component_id() == icomp) {
      Key x_key_base = splitPhase(soln_names_[n0]).first;
      Key x_key = mergePhase(x_key_base, bcs_[i]->component_phase());

      auto& bc_model = op_bcs_[x_key]->bc_model();
      auto& bc_value = op_bcs_[x_key]->bc_value();

      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value[f] = it->second[0];
      }
    }

    if (bcs_[i]->get_bc_name() == "saturation") {
      auto& bc_model_s = op_bcs_[saturation_liquid_key_]->bc_model();
      auto& bc_value_s = op_bcs_[saturation_liquid_key_]->bc_value();

      Key x_key_base = splitPhase(soln_names_[n0]).first;
      Key x_key = mergePhase(x_key_base, MULTIPHASE_PHASE_LIQUID);
      auto& bc_model_x = op_bcs_[x_key]->bc_model();
      auto& bc_value_x = op_bcs_[x_key]->bc_value();

      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        bc_model_s[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value_s[f] = it->second[0];

        bc_model_x[f] = Operators::OPERATOR_BC_DIRICHLET; // a huck
        bc_value_x[f] = 0.0;
      }
    }
  }

  // mark missing boundary conditions as zero flux conditions
  missed_bc_faces_ = 0;

  for (int f = 0; f < nfaces_owned_; f++) {
    auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);

    if (cells.size() == 1) {
      for (auto it = op_bcs_.begin(); it != op_bcs_.end(); ++it) {
        auto& bc_model = it->second->bc_model();

        if (bc_model[f] == Operators::OPERATOR_BC_NONE) {
          bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
          // if (i == 0) missed_bc_faces_++;
        }
      }
    }
  }

  PopulateSecondaryBCs_();
}


/* *******************************************************************
* Populate boundary conditions for derived fields
******************************************************************* */
void
Multiphase_PK::CheckCompatibilityBCs(const Key& keyr, const Key& gname)
{
  auto& bc_model_r = op_bcs_[keyr]->bc_model();
  auto& bc_model_g = op_bcs_[gname]->bc_model();

  for (int f = 0; f != nfaces_owned_; ++f) {
    if (bc_model_r[f] == Operators::OPERATOR_BC_DIRICHLET &&
        bc_model_g[f] == Operators::OPERATOR_BC_NEUMANN)
      AMANZI_ASSERT(false);
  }
}


/* *******************************************************************
* Populate boundary conditions for derived fields
******************************************************************* */
void
Multiphase_PK::PopulateSecondaryBCs_()
{
  auto& bc_model_pg = op_bcs_[pressure_gas_key_]->bc_model();
  auto& bc_value_pg = op_bcs_[pressure_gas_key_]->bc_value();

  auto& bc_model_pl = op_bcs_[pressure_liquid_key_]->bc_model();
  auto& bc_value_pl = op_bcs_[pressure_liquid_key_]->bc_value();

  auto& bc_model_sl = op_bcs_[saturation_liquid_key_]->bc_model();
  auto& bc_value_sl = op_bcs_[saturation_liquid_key_]->bc_value();

  auto& bc_model_advl = op_bcs_[advection_liquid_key_]->bc_model();
  auto& bc_value_advl = op_bcs_[advection_liquid_key_]->bc_value();

  auto& bc_model_advg = op_bcs_[advection_gas_key_]->bc_model();
  auto& bc_value_advg = op_bcs_[advection_gas_key_]->bc_value();

  auto& bc_model_etal = bc_model_pl;
  auto& bc_value_etal = bc_value_pl;

  for (int i = 0; i < bcs_.size(); ++i) {
     if (bcs_[i]->get_bc_name() == "concentration" && bcs_[i]->component_id() == 0) { // FIX ME; icomp == 0 assumption; i.e., one additional component only apart from water
       Key x_key_base = splitPhase(soln_names_[1]).first;
       Key x_key = mergePhase(x_key_base, bcs_[i]->component_phase());

       bc_model_etal = op_bcs_[x_key]->bc_model();
       bc_value_etal = op_bcs_[x_key]->bc_value();
    }
  }

  bc_model_pg = bc_model_pl;
  bc_model_advl = bc_model_pl;
  bc_model_advg = bc_model_pg;

  const auto& sl_c =
    *S_->Get<CompositeVector>(saturation_liquid_key_, Tags::DEFAULT).ViewComponent("cell");

  for (int f = 0; f != nfaces_wghost_; ++f) {
    if (bc_model_pl[f] == Operators::OPERATOR_BC_DIRICHLET) {
      int c = getFaceOnBoundaryInternalCell(*mesh_, f);

      //bc_value_pg[f] =
      //  bc_value_pl[f] + wrm_->second[(*wrm_->first)[c]]->capillaryPressure(sl_c[0][c]);
      
      // debugging
      /*
      std::cout<<"wrm_->first = "<<c<<"; "<<(*wrm_->first)[c]<<std::endl;
      std::cout<<"wrm_->second = "<<wrm_->second[(*wrm_->first)[c]]<<std::endl;
      */

      bc_value_pg[f] = bc_value_pl[f] + wrm_->second[(*wrm_->first)[c]]->capillaryPressure(bc_value_sl[f]);

      // FIX ME: NEEDS TRANSMISSIBILITY (and below as well)
      bc_value_advl[f] =  wrm_->second[(*wrm_->first)[c]]->k_relative(bc_value_sl[f], MULTIPHASE_PHASE_LIQUID);
      bc_value_advl[f] *= bc_value_etal[f] * (1.0 / mu_l_) * 5e-20;

      //if (bc_value_sl[f] < 1.0 - 1.e-12) { // i.e., if gas phase is present
        double bc_sg = 1.0 - bc_value_sl[f];

        //const auto& T = *S_->Get<CompositeVector)(temperature_key_, Tags::DEFAULT).ViewComponent("cell");
        //double R = eos_plist_.get<double>("ideal gas constant", 8.31446261815324);
        double bc_etag = bc_value_pg[f] / (303.0 * 8.31446261815324); // FIX ME temperature T, gas constant R needs to be extracted from xml file

        bc_value_advg[f] =  wrm_->second[(*wrm_->first)[c]]->k_relative(bc_sg, MULTIPHASE_PHASE_GAS);
        bc_value_advg[f] *= bc_etag * (1.0 / mu_g_) * 5e-20;
    //}

    } else if (bc_model_pl[f] == Operators::OPERATOR_BC_NEUMANN) {
      bc_value_pg[f] = 0.0;
      bc_value_advl[f] = 0.0;
      bc_value_advg[f] = 0.0;
    }
  }
}

} // namespace Multiphase
} // namespace Amanzi

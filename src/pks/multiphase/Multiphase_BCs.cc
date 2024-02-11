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
#include "MeshAlgorithms.hh"

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
      auto phase = bcs_[i]->component_phase();
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
    auto cells = mesh_->getFaceCells(f);

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

  // current advection BC implementation for one component only
  auto& bc_model_advl = op_bcs_[advection_liquid_key_]->bc_model();
  auto& bc_value_advl = op_bcs_[advection_liquid_key_]->bc_value();

  auto& bc_model_advg = op_bcs_[advection_gas_key_]->bc_model();
  auto& bc_value_advg = op_bcs_[advection_gas_key_]->bc_value();

  auto& bc_model_etal = op_bcs_[mol_density_liquid_key_]->bc_model();
  auto& bc_value_etal = op_bcs_[mol_density_liquid_key_]->bc_value();

  bc_model_pg = bc_model_pl;
  bc_model_advl = bc_model_pl;
  bc_model_advg = bc_model_pg;

  double bc_value_sg, bc_value_etag;

  // ideal gas law parameters and variables
  auto T = *S_->GetPtr<CompositeVector>(temperature_key_, Tags::DEFAULT);
  const Epetra_MultiVector& T_c = *T.ViewComponent("cell", true);

  auto eos_parameter_list = glist_->sublist("state").sublist("evaluators").sublist("molar_density_gas").sublist("EOS parameters");
  double R;
  if(eos_parameter_list.isParameter("ideal gas constant")) {
    R = eos_parameter_list.get<double>("ideal gas constant");
  } else {
    R = 8.31446261815324;
  }

  for (int f = 0; f != nfaces_wghost_; ++f) {
    if (bc_model_pl[f] == Operators::OPERATOR_BC_DIRICHLET) {
      int c = getFaceOnBoundaryInternalCell(*mesh_, f);

      // calculate gas phase BCs from liquid saturation/pressure
      bc_value_pg[f] = bc_value_pl[f] + wrm_->second[(*wrm_->first)[c]]->capillaryPressure(bc_value_sl[f]);
      bc_value_sg = 1.0 - bc_value_sl[f];

      bc_value_advl[f] =  wrm_->second[(*wrm_->first)[c]]->k_relative(bc_value_sl[f], MULTIPHASE_PHASE_LIQUID);
      bc_value_advl[f] *= bc_value_etal[f] * (1.0 / mu_l_);

      bc_value_advg[f] = 0.0; // dummy value for now

      // need more info on the phase: for example, need to specify gas BC concentration, but currently bc->component_phase is either 1 or 2, for liquid or gas.  
      /*
      bc_value_sg = 1.0 - bc_value_sl[f];
      bc_value_etag = bc_value_pg[f] / (R * T_c[0][c]);
      
      bc_value_advg[f] =  wrm_->second[(*wrm_->first)[c]]->k_relative(bc_value_sl[f], MULTIPHASE_PHASE_GAS);
      bc_value_advg[f] *= bc_value_etag * (1.0 / mu_g_);
      */

    } else if (bc_model_pl[f] == Operators::OPERATOR_BC_NEUMANN) {
      bc_value_pg[f] = 0.0;
      bc_value_advl[f] = 0.0;
      bc_value_advg[f] = 0.0;
    }
  }
}

} // namespace Multiphase
} // namespace Amanzi

/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>

#include "ReconstructionCell.hh"
#include "OperatorDefs.hh"
#include "TransportExplicit_PK.hh"
#include "UniqueLocalIndex.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************* 
* Routine takes a parallel vector C and returns parallel vector F(C). 
****************************************************************** */
void TransportExplicit_PK::FunctionalTimeDerivative(
    double t, const Epetra_Vector& component,
    Epetra_Vector& f_component)
{
  auto darcy_flux = S_->GetFieldData(darcy_flux_key_)->ViewComponent("face", true);

  // distribute vector
  Epetra_Vector component_tmp(component);
  component_tmp.Import(component, *tcc->importer("cell"), Insert);

  // transport routines need an RCP pointer
  Teuchos::RCP<const Epetra_Vector> component_rcp(&component_tmp, false);

  Teuchos::ParameterList plist = tp_list_->sublist("reconstruction");
  lifting_->Init(plist);
  lifting_->ComputeGradient(component_rcp);

  // extract boundary conditions for the current component
  std::vector<int> bc_model(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost);

  for (int m = 0; m < bcs_.size(); m++) {
    std::vector<int>& tcc_index = bcs_[m]->tcc_index();
    int ncomp = tcc_index.size();

    for (int i = 0; i < ncomp; i++) {
      if (current_component_ == tcc_index[i]) {
        for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
          int f = it->first;

          std::vector<double>& values = it->second;

          bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_value[f] = values[i];
        }
      }
    }
  }

  limiter_->Init(plist, darcy_flux);
  limiter_->ApplyLimiter(component_rcp, 0, lifting_->gradient(), bc_model, bc_value);
  limiter_->gradient()->ScatterMasterToGhosted("cell");

  // ADVECTIVE FLUXES
  // We assume that limiters made their job up to round-off errors.
  // Min-max condition will enforce robustness w.r.t. these errors.
  int c1, c2;
  double u, u1, u2, umin, umax, upwind_tcc, tcc_flux;

  f_component.PutScalar(0.0);
  for (int f = 0; f < nfaces_wghost; f++) {
    c1 = (upwind_cells_[f].size() > 0) ? upwind_cells_[f][0] : -1;
    c2 = (downwind_cells_[f].size() > 0) ? downwind_cells_[f][0] : -1;

    if (c1 >= 0 && c2 >= 0) {
      u1 = component[c1];
      u2 = component[c2];
      umin = std::min(u1, u2);
      umax = std::max(u1, u2);
    } else if (c1 >= 0) {
      u1 = u2 = umin = umax = component[c1];
    } else if (c2 >= 0) {
      u1 = u2 = umin = umax = component[c2];
    }

    u = fabs((*darcy_flux)[0][f]);
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

    if (c1 >= 0 && c1 < ncells_owned && c2 >= 0 && c2 < ncells_owned) {
      upwind_tcc = limiter_->getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c1] -= tcc_flux;
      f_component[c2] += tcc_flux;
    } else if (c1 >= 0 && c1 < ncells_owned && (c2 >= ncells_owned || c2 < 0)) {
      upwind_tcc = limiter_->getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c1] -= tcc_flux;
    } else if (c1 >= ncells_owned && c2 >= 0 && c2 < ncells_owned) {
      upwind_tcc = limiter_->getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c2] += tcc_flux;
    }
  }

  // process external sources
  if (srcs_.size() != 0) {
    ComputeSources_(t, 1.0, f_component, 
                    *component_rcp, current_component_, current_component_);
  }

  for (int c = 0; c < ncells_owned; c++) {  // calculate conservative quantatity
    double a = t / dt_;
    double ws = a * (*ws_end)[0][c] + (1.0 - a) * (*ws_start)[0][c]; 
    double vol_phi_ws = mesh_->cell_volume(c) * (*phi)[0][c] * ws;
    f_component[c] /= vol_phi_ws;
  }

  // BOUNDARY CONDITIONS for ADVECTION
  const auto& flux_map = S_->GetFieldData(darcy_flux_key_)->Map().Map("face", true);

  for (int m = 0; m < bcs_.size(); m++) {
    std::vector<int>& tcc_index = bcs_[m]->tcc_index();
    int ncomp = tcc_index.size();

    for (int i = 0; i < ncomp; i++) {
      if (current_component_ == tcc_index[i]) {
        for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
          int f = it->first;

          std::vector<double>& values = it->second;

          if (downwind_cells_[f].size() > 0 && f < nfaces_owned) {
            for (int k = 0; k < downwind_cells_[f].size(); k++) {
              c2 = downwind_cells_[f][k];
              int g = flux_map->FirstPointInElement(f);
              g += Operators::UniqueIndexFaceToCells(*mesh_, f, c2);

              if (c2 >=0) {
                u = fabs((*darcy_flux)[0][g]);
                double vol_phi_ws = mesh_->cell_volume(c2) * (*phi)[0][c2] * (*ws_start)[0][c2];
                tcc_flux = u * values[i];
                f_component[c2] += tcc_flux / vol_phi_ws;
              }
            }
          }
        }
      }
    }
  }
}  


/* ******************************************************************* 
* Routine takes a parallel overlapping vector C and returns parallel
* overlapping vector F(C). Old version.
****************************************************************** */
void TransportExplicit_PK::DudtOld(double t,
                           const Epetra_Vector& component,
                           Epetra_Vector& f_component)
{
  auto darcy_flux = S_->GetFieldData(darcy_flux_key_)->ViewComponent("face", true);

  // transport routines need an RCP pointer
  Teuchos::RCP<const Epetra_Vector> component_rcp(&component, false);

  Teuchos::ParameterList plist = tp_list_->sublist("reconstruction");
  lifting_->Init(plist);
  lifting_->ComputeGradient(component_rcp);

  // extract boundary conditions for the current component
  std::vector<int> bc_model(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost);

  for (int m = 0; m < bcs_.size(); m++) {
    std::vector<int>& tcc_index = bcs_[m]->tcc_index();
    int ncomp = tcc_index.size();

    for (int i = 0; i < ncomp; i++) {
      if (current_component_ == tcc_index[i]) {
        for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
          int f = it->first;

          std::vector<double>& values = it->second;

          bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_value[f] = values[i];
        }
      }
    }
  }

  limiter_->Init(plist, darcy_flux);
  limiter_->ApplyLimiter(component_rcp, 0, lifting_->gradient(), bc_model, bc_value);
  limiter_->gradient()->ScatterMasterToGhosted("cell");

  // ADVECTIVE FLUXES
  // We assume that limiters made their job up to round-off errors.
  // Min-max condition will enforce robustness w.r.t. these errors.
  int c1, c2;
  double u, u1, u2, umin, umax, upwind_tcc, tcc_flux;

  f_component.PutScalar(0.0);
  for (int f = 0; f < nfaces_wghost; f++) {  // loop over master and slave faces
    c1 = (upwind_cells_[f].size() > 0) ? upwind_cells_[f][0] : -1;
    c2 = (downwind_cells_[f].size() > 0) ? downwind_cells_[f][0] : -1;

    if (c1 >= 0 && c2 >= 0) {
      u1 = component[c1];
      u2 = component[c2];
      umin = std::min(u1, u2);
      umax = std::max(u1, u2);
    } else if (c1 >= 0) {
      u1 = u2 = umin = umax = component[c1];
    } else if (c2 >= 0) {
      u1 = u2 = umin = umax = component[c2];
    }

    u = fabs((*darcy_flux)[0][f]);
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

    if (c1 >= 0 && c1 < ncells_owned && c2 >= 0 && c2 < ncells_owned) {
      upwind_tcc = limiter_->getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c1] -= tcc_flux;
      f_component[c2] += tcc_flux;
    } else if (c1 >= 0 && c1 < ncells_owned && (c2 >= ncells_owned || c2 < 0)) {
      upwind_tcc = limiter_->getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c1] -= tcc_flux;
    } else if (c1 >= ncells_owned && c2 >= 0 && c2 < ncells_owned) {
      upwind_tcc = limiter_->getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c2] += tcc_flux;
    }
  }

  // process external sources
  if (srcs_.size() != 0) {
    ComputeSources_(t, 1.0, f_component, 
                    *component_rcp, current_component_, current_component_);
  }

  for (int c = 0; c < ncells_owned; c++) {  // calculate conservative quantatity
    double vol_phi_ws = mesh_->cell_volume(c) * (*phi)[0][c] * (*ws_start)[0][c];
    f_component[c] /= vol_phi_ws;
  }

  // BOUNDARY CONDITIONS for ADVECTION
  for (int m = 0; m < bcs_.size(); m++) {
    std::vector<int>& tcc_index = bcs_[m]->tcc_index();
    int ncomp = tcc_index.size();

    for (int i = 0; i < ncomp; i++) {
      if (current_component_ == tcc_index[i]) {
        for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
          int f = it->first;

          std::vector<double>& values = it->second;

          if (downwind_cells_[f].size() > 0 && f < nfaces_owned) {
            for (int j = 0; j < downwind_cells_[f].size(); j++) {
              c2 = downwind_cells_[f][0];
              if (c2 < 0) continue;
              u = fabs((*darcy_flux)[0][f]);
              double vol_phi_ws = mesh_->cell_volume(c2) * (*phi)[0][c2] * (*ws_start)[0][c2];
              tcc_flux = u * values[i];
              f_component[c2] += tcc_flux / vol_phi_ws;
            }
          }
        }
      }
    }
  }
}




}  // namespace Transport
}  // namespace Amanzi



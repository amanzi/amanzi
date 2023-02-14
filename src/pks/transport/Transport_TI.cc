/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

*/

#include <algorithm>

#include "Mesh_Algorithms.hh"
#include "OperatorDefs.hh"
#include "TransportExplicit_PK.hh"
#include "UniqueLocalIndex.hh"

namespace Amanzi {
namespace Transport {

/* *******************************************************************
* Monotonic Upstream-centered Scheme for Conservation Laws
****************************************************************** */
void
Transport_PK::FunctionalTimeDerivative_MUSCL_(double t,
                                              const CompositeVector& component,
                                              CompositeVector& func,
                                              bool scale)
{
  auto flowrate = S_->Get<CompositeVector>(vol_flowrate_key_).ViewComponent("face", true);

  // distribute vector
  component.ScatterMasterToGhosted("cell");
  const auto& component_c = *component.ViewComponent("cell", true);
  auto& f_c = *func.ViewComponent("cell", true);

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
          bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_value[f] = it->second[i];
        }
      }
    }
  }

  // transport routines need an RCP pointer
  Teuchos::RCP<const Epetra_MultiVector> component_rcp(&component_c, false);

  Teuchos::ParameterList plist = tp_list_->sublist("reconstruction");
  lifting_->Init(plist);
  lifting_->Compute(component_rcp);

  limiter_->Init(plist, flowrate);
  limiter_->ApplyLimiter(component_rcp, 0, lifting_, bc_model, bc_value);
  lifting_->data()->ScatterMasterToGhosted("cell");

  // ADVECTIVE FLUXES
  // We assume that limiters made their job up to round-off errors.
  // Min-max condition will enforce robustness w.r.t. these errors.
  int c1, c2;
  double u, u1, u2, umin, umax, upwind_tcc, tcc_flux;

  const auto& flux_map = S_->Get<CompositeVector>(vol_flowrate_key_).Map().Map("face", true);

  func.PutScalar(0.0);
  if (mesh_->getSpaceDimension() == mesh_->getManifoldDimension()) {
    for (int f = 0; f < nfaces_wghost; f++) {
      c1 = (upwind_cells_[f].size() > 0) ? upwind_cells_[f][0] : -1;
      c2 = (downwind_cells_[f].size() > 0) ? downwind_cells_[f][0] : -1;

      if (c1 >= 0 && c2 >= 0) {
        u1 = component_c[0][c1];
        u2 = component_c[0][c2];
        umin = std::min(u1, u2);
        umax = std::max(u1, u2);
      } else if (c1 >= 0) {
        u1 = u2 = umin = umax = component_c[0][c1];
      } else if (c2 >= 0) {
        u1 = u2 = umin = umax = component_c[0][c2];
      }

      int g = flux_map->FirstPointInElement(f);
      u = fabs((*flowrate)[0][g]);
      const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);

      if (c1 >= 0 && c1 < ncells_owned && c2 >= 0 && c2 < ncells_owned) {
        upwind_tcc = limiter_->getValue(c1, xf);
        upwind_tcc = std::max(upwind_tcc, umin);
        upwind_tcc = std::min(upwind_tcc, umax);

        tcc_flux = u * upwind_tcc;
        f_c[0][c1] -= tcc_flux;
        f_c[0][c2] += tcc_flux;
      } else if (c1 >= 0 && c1 < ncells_owned && (c2 >= ncells_owned || c2 < 0)) {
        upwind_tcc = limiter_->getValue(c1, xf);
        upwind_tcc = std::max(upwind_tcc, umin);
        upwind_tcc = std::min(upwind_tcc, umax);

        tcc_flux = u * upwind_tcc;
        f_c[0][c1] -= tcc_flux;
      } else if (c1 >= ncells_owned && c2 >= 0 && c2 < ncells_owned) {
        upwind_tcc = limiter_->getValue(c1, xf);
        upwind_tcc = std::max(upwind_tcc, umin);
        upwind_tcc = std::min(upwind_tcc, umax);

        tcc_flux = u * upwind_tcc;
        f_c[0][c2] += tcc_flux;
      }
    }
  } else {
    for (int f = 0; f < nfaces_wghost; f++) {
      double flux_in(0.0), tcc_out(0.0);

      for (int n = 0; n < upwind_cells_[f].size(); ++n) {
        int c = upwind_cells_[f][n];
        u = upwind_flux_[f][n];
        tcc_out += u * component_c[0][c];
      }

      for (int n = 0; n < downwind_cells_[f].size(); ++n) { flux_in -= downwind_flux_[f][n]; }
      if (flux_in == 0.0) flux_in = 1e-12;

      // update solutes
      for (int n = 0; n < upwind_cells_[f].size(); ++n) {
        int c = upwind_cells_[f][n];
        u = upwind_flux_[f][n];

        if (c < ncells_owned) { f_c[0][c] -= u * component_c[0][c]; }
      }

      for (int n = 0; n < downwind_cells_[f].size(); ++n) {
        int c = downwind_cells_[f][n];
        u = downwind_flux_[f][n];

        if (c < ncells_owned) {
          double tmp = u / flux_in;
          f_c[0][c] -= tmp * tcc_out;
        }
      }
    }
  }

  // BOUNDARY CONDITIONS for ADVECTION
  func.PutScalarGhosted(0.0);
  int flag(0);

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
              if (c2 >= ncells_owned) flag = 1;

              int g = flux_map->FirstPointInElement(f);
              g += Operators::UniqueIndexFaceToCells(*mesh_, f, c2);

              if (c2 >= 0) {
                u = fabs((*flowrate)[0][g]);
                tcc_flux = u * values[i];
                f_c[0][c2] += tcc_flux;
              }
            }
          }
        }
      }
    }
  }

  int flag_tmp(flag);
  mesh_->getComm()->MaxAll(&flag_tmp, &flag, 1);
  if (flag == 1) func.GatherGhostedToMaster("cell");

  // process external sources
  if (srcs_.size() != 0) {
    ComputeSources_(t, 1.0, f_c, *component_rcp, current_component_, current_component_);
  }

  // optional convertion to the primary variable: dC/dt = F / wc
  if (scale) {
    for (int c = 0; c < ncells_owned; c++) {
      double a = t / dt_;
      double wc = a * (*wc_end)[0][c] + (1.0 - a) * (*wc_start)[0][c];
      double vol_wc = mesh_->getCellVolume(c) * wc;
      f_c[0][c] /= vol_wc;
    }
  }

  if (vo_->getVerbLevel() > Teuchos::VERB_MEDIUM) {
    limiter_->limiter()->MeanValue(&limiter_mean_);
  }
}


/* *******************************************************************
* Flux corrected transport
****************************************************************** */
void
Transport_PK::FunctionalTimeDerivative_FCT_(double t,
                                            const CompositeVector& component,
                                            CompositeVector& func)
{
  auto flowrate = S_->Get<CompositeVector>(vol_flowrate_key_).ViewComponent("face", true);

  S_->GetEvaluator(wc_key_).Update(*S_, "transport");
  auto weight1 = S_->Get<CompositeVector>(wc_key_).ViewComponent("cell", true);
  auto weight0 = S_->Get<CompositeVector>(prev_wc_key_).ViewComponent("cell", true);

  // distribute vector
  // distribute vector
  component.ScatterMasterToGhosted("cell");
  const auto& component_c = *component.ViewComponent("cell", true);
  auto& f_c = *func.ViewComponent("cell", true);

  // extract boundary conditions for the current component
  auto bcs = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bcs->bc_model();
  std::vector<double>& bc_value = bcs->bc_value();

  for (int m = 0; m < bcs_.size(); m++) {
    std::vector<int>& tcc_index = bcs_[m]->tcc_index();
    int ncomp = tcc_index.size();

    for (int i = 0; i < ncomp; i++) {
      if (current_component_ == tcc_index[i]) {
        for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
          int f = it->first;
          bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_value[f] = it->second[i];
        }
      }
    }
  }

  // transport routines need an RCP pointer
  int c1, c2;
  double u, tcc_flux;

  Teuchos::RCP<const Epetra_MultiVector> component_rcp(&component_c, false);

  Teuchos::ParameterList plist = tp_list_->sublist("reconstruction");
  lifting_->Init(plist);
  lifting_->Compute(component_rcp);
  lifting_->data()->ScatterMasterToGhosted("cell");

  // low-order and high-order fluxes
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_)->SetGhosted(false)->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  CompositeVector flux_lo(cvs), flux_ho(cvs), flux_numer(cvs);

  auto& flux_lo_f = *flux_lo.ViewComponent("face");
  auto& flux_ho_f = *flux_ho.ViewComponent("face");
  auto& flux_numer_f = *flux_numer.ViewComponent("face");

  for (int f = 0; f < nfaces_owned; f++) {
    auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    int ncells = cells.size();

    c1 = (upwind_cells_[f].size() > 0) ? upwind_cells_[f][0] : -1;
    c2 = (downwind_cells_[f].size() > 0) ? downwind_cells_[f][0] : -1;

    u = dt_ * fabs((*flowrate)[0][f]);
    const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);

    if (c1 == cells[0]) {
      flux_lo_f[0][f] = u * component_c[0][c1];
      flux_ho_f[0][f] = u * lifting_->getValue(c1, xf);
    } else if (c1 >= 0 && ncells == 2) {
      flux_lo_f[0][f] = -u * component_c[0][c1];
      flux_ho_f[0][f] = -u * lifting_->getValue(c1, xf);
    } else if (c1 < 0) {
      flux_lo_f[0][f] = -u * component_c[0][c2];
      flux_ho_f[0][f] = -u * lifting_->getValue(c2, xf);
    } else if (c2 < 0) {
      flux_lo_f[0][f] = u * component_c[0][c1];
      flux_ho_f[0][f] = u * lifting_->getValue(c1, xf);
    }
  }

  // boundary fluxes
  const auto& flux_map = S_->Get<CompositeVector>(vol_flowrate_key_).Map().Map("face", true);

  for (int m = 0; m < bcs_.size(); m++) {
    std::vector<int>& tcc_index = bcs_[m]->tcc_index();
    int ncomp = tcc_index.size();

    for (int i = 0; i < ncomp; i++) {
      if (current_component_ == tcc_index[i]) {
        for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
          int f = it->first;

          if (downwind_cells_[f].size() > 0 && f < nfaces_owned) {
            for (int k = 0; k < downwind_cells_[f].size(); k++) {
              c2 = downwind_cells_[f][k];
              int g = flux_map->FirstPointInElement(f);
              g += Operators::UniqueIndexFaceToCells(*mesh_, f, c2);

              if (c2 >= 0) {
                u = dt_ * fabs((*flowrate)[0][g]);
                tcc_flux = u * it->second[i];
                flux_lo_f[0][f] = -tcc_flux;
                flux_ho_f[0][f] = -tcc_flux;
              }
            }
          }
        }
      }
    }
  }

  // compute bounds-preserving flux
  fct_->Init(component_rcp, 0, weight0, weight1);
  fct_->Compute(flux_lo, flux_ho, *bcs, flux_numer);

  // update functional for time integrator
  CompositeVectorSpace cvs2;
  cvs2.SetMesh(mesh_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  CompositeVector residual(cvs2);
  auto& residual_c = *residual.ViewComponent("cell", true);

  for (int f = 0; f < nfaces_owned; ++f) {
    tcc_flux = flux_numer_f[0][f] / dt_;

    auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    int ncells = cells.size();

    if (ncells == 2) {
      residual_c[0][cells[0]] -= tcc_flux;
      residual_c[0][cells[1]] += tcc_flux;
    } else {
      residual_c[0][cells[0]] -= tcc_flux;
    }
  }
  residual.GatherGhostedToMaster();

  for (int c = 0; c < ncells_owned; ++c) { f_c[0][c] = residual_c[0][c]; }

  // process external sources
  if (srcs_.size() != 0) {
    ComputeSources_(t, 1.0, f_c, *component_rcp, current_component_, current_component_);
  }

  for (int c = 0; c < ncells_owned; c++) { // calculate conservative quantity
    double a = t / dt_;
    double wc = a * (*wc_end)[0][c] + (1.0 - a) * (*wc_start)[0][c];
    double vol_wc = mesh_->getCellVolume(c) * wc;
    f_c[0][c] /= vol_wc;
  }
}

} // namespace Transport
} // namespace Amanzi

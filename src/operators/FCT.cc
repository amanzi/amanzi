/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Kontantin Lipnikov (lipnikov@lanl.gov)

  Flux Corrected Transport.
*/

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"

#include "FCT.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Flux is from 1st to 2nd cell in the list face->cells
****************************************************************** */
void FCT::Compute(const CompositeVector& flux_lo,
                  const CompositeVector& flux_ho,
                  const BCs& bc,
                  CompositeVector& flux)
{
  int nfaces_owned = mesh0_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int ncells_owned = mesh0_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh0_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  AmanziMesh::Entity_ID_List cells, faces;

  const auto& flux_lo_f = *flux_lo.ViewComponent("face"); 
  const auto& flux_ho_f = *flux_ho.ViewComponent("face"); 
  auto& flux_f = *flux.ViewComponent("face"); 

  // allocate memory
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh0_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1);
  CompositeVector dlo(cvs), pos(cvs), neg(cvs);

  auto& dlo_c = *dlo.ViewComponent("cell", true);
  auto& pos_c = *pos.ViewComponent("cell", true);
  auto& neg_c = *neg.ViewComponent("cell", true);

  // collect positive and negative fluxes in each mesh cell
  int dir;

  for (int f = 0; f < nfaces_owned; ++f) {
    mesh0_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    mesh0_->face_normal(f, false, cells[0], &dir);
    double tmp0 = -flux_lo_f[0][f];
    double dtmp = -flux_ho_f[0][f] + flux_lo_f[0][f];

    for (int i = 0; i < ncells; ++i) {
      int c = cells[i];
      pos_c[0][c] += std::max(dtmp, 0.0);
      neg_c[0][c] += std::min(dtmp, 0.0);
      dlo_c[0][c] += tmp0;
      dtmp = -dtmp;
      tmp0 = -tmp0;
    }
  }

  dlo.GatherGhostedToMaster();
  pos.GatherGhostedToMaster();
  neg.GatherGhostedToMaster();

  // collect cell-limiters for positive and negative fluxes
  const auto& bc_model = bc.bc_model();
  const auto& bc_value = bc.bc_value();

  Teuchos::RCP<const CompositeVector> bounds;
  if (limiter_->get_external_bounds())
    bounds = limiter_->get_bounds();
  else
    bounds = limiter_->BoundsForCells(*field_, bc_model, bc_value, OPERATOR_LIMITER_STENCIL_C2C_ALL);
  auto& bounds_c = *bounds->ViewComponent("cell");

  double Qmin, Qmax;
  std::vector<double> alpha(nfaces_owned, 1.0);

  for (int c = 0; c < ncells_owned; ++c) {
    double vol0 = mesh0_->cell_volume(c);
    double vol1 = mesh1_->cell_volume(c);

    if (weight0_ != Teuchos::null) vol0 *= (*weight0_)[0][c];
    if (weight0_ != Teuchos::null) vol1 *= (*weight1_)[0][c];

    Qmin = Qmax = 1.0;
    if (neg_c[0][c] != 0.0) 
      Qmin = std::min(0.0, (vol1 * bounds_c[0][c] - vol0 * (*field_)[0][c] - dlo_c[0][c])) / neg_c[0][c];

    if (pos_c[0][c] != 0.0) 
      Qmax = std::max(0.0, (vol1 * bounds_c[1][c] - vol0 * (*field_)[0][c] - dlo_c[0][c])) / pos_c[0][c];

    neg_c[0][c] = std::fabs(Qmin);  // re-using allocated memory
    pos_c[0][c] = std::fabs(Qmax); 
  }

  pos.ScatterMasterToGhosted();
  neg.ScatterMasterToGhosted();

  // move cell-limiters to face limiters
  for (int f = 0; f < nfaces_owned; ++f) {
    mesh0_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

    if (cells.size() == 2) {
      double tmp = -flux_ho_f[0][f] + flux_lo_f[0][f];
    
      if (tmp > 0.0)
        alpha[f] = std::min({ alpha[f], pos_c[0][cells[0]], neg_c[0][cells[1]] });
      else
        alpha[f] = std::min({ alpha[f], neg_c[0][cells[0]], pos_c[0][cells[1]] });
    }
  }

  alpha_mean_ = 0.0;
  for (int f = 0; f < nfaces_owned; ++f) {
    alpha_mean_ += alpha[f];
    flux_f[0][f] = flux_lo_f[0][f] + alpha[f] * (flux_ho_f[0][f] - flux_lo_f[0][f]);
  }

  double tmp(alpha_mean_);
  mesh0_->get_comm()->SumAll(&tmp, &alpha_mean_, 1);
  alpha_mean_ /= flux_f.GlobalLength();
}

}  // namespace Operators
}  // namespace Amanzi


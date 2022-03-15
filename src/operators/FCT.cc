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
  int ncells_owned = mesh0_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_owned = mesh0_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  AmanziMesh::Entity_ID_List cells, faces;

  const auto& flux_lo_f = *flux_lo.ViewComponent("face"); 
  const auto& flux_ho_f = *flux_ho.ViewComponent("face"); 
  auto& flux_f = *flux.ViewComponent("face"); 

  int dir;
  std::vector<double> dlo(ncells_owned, 0.0);
  std::vector<double> pos(ncells_owned, 0.0), neg(ncells_owned, 0.0);

  for (int f = 0; f < nfaces_owned; ++f) {
    mesh0_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    mesh0_->face_normal(f, false, cells[0], &dir);
    double tmp0 = -flux_lo_f[0][f];
    double dtmp = -flux_ho_f[0][f] + flux_lo_f[0][f];

    for (int i = 0; i < ncells; ++i) {
      int c = cells[i];
      pos[c] += std::max(dtmp, 0.0);
      neg[c] += std::min(dtmp, 0.0);
      dlo[c] += tmp0;
      dtmp = -dtmp;
      tmp0 = -tmp0;
    }
  }

  const auto& bc_model = bc.bc_model();
  const auto& bc_value = bc.bc_value();
  auto bounds = lifting_->BoundsForCells(*field_, bc_model, bc_value, OPERATOR_LIMITER_STENCIL_C2C_ALL);
  auto& bounds_c = *bounds->ViewComponent("cell");

  double Qmin, Qmax;
  std::vector<double> alpha(nfaces_owned, 1.0);

  for (int c = 0; c < ncells_owned; ++c) {
    double vol0 = mesh0_->cell_volume(c);
    double vol1 = mesh1_->cell_volume(c);

    Qmin = Qmax = 1.0;
    if (neg[c] != 0.0) 
      Qmin = (vol1 * bounds_c[0][c] - vol0 * (*field_)[0][c] - dlo[c]) / neg[c];

    if (pos[c] != 0.0) 
      Qmax = (vol1 * bounds_c[1][c] - vol0 * (*field_)[0][c] - dlo[c]) / pos[c];

    neg[c] = std::fabs(Qmin);  // re-using allocated memory
    pos[c] = std::fabs(Qmax); 
  }

  for (int f = 0; f < nfaces_owned; ++f) {
    mesh0_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

    if (cells.size() == 2) {
      double tmp = -flux_ho_f[0][f] + flux_lo_f[0][f];
    
      if (tmp > 0.0)
        alpha[f] = std::min({ alpha[f], pos[cells[0]], neg[cells[1]] });
      else
        alpha[f] = std::min({ alpha[f], neg[cells[0]], pos[cells[1]] });
    }
  }

  alpha_mean_ = 0.0;
  for (int f = 0; f < nfaces_owned; ++f) {
    alpha_mean_ += alpha[f];
    flux_f[0][f] = flux_lo_f[0][f] + alpha[f] * (flux_ho_f[0][f] - flux_lo_f[0][f]);
  }

  alpha_mean_ /= nfaces_owned;
}

}  // namespace Operators
}  // namespace Amanzi


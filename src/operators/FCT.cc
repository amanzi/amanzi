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
*
****************************************************************** */
void FCT::Compute(const CompositeVector& flux_lo,
                  const CompositeVector& flux_ho, CompositeVector& flux)
{
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  AmanziMesh::Entity_ID_List cells, faces;

  const auto& flux_lo_f = *flux_lo.ViewComponent("face"); 
  const auto& flux_ho_f = *flux_ho.ViewComponent("face"); 
  auto& flux_f = *flux.ViewComponent("face"); 

  int dir;
  std::vector<double> dlo(ncells_owned, 0.0);
  std::vector<double> pos(ncells_owned, 0.0), neg(ncells_owned, 0.0);

  for (int f = 0; f < nfaces_owned; ++f) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    mesh_->face_normal(f, false, cells[0], &dir);
    double tmp = (flux_ho_f[0][f] - flux_lo_f[0][f]) * dir;

    for (int i = 0; i < ncells; ++i) {
      int c = cells[i];
      pos[c] += std::max(tmp, 0.0);
      neg[c] += std::min(tmp, 0.0);
      dlo[c] += tmp;
      tmp = -tmp;
    }
  }

  double Qmin, Qmax;
  const auto& bounds = *lifting_->bounds()->ViewComponent("cell");
  std::vector<double> alpha(nfaces_owned, 1.0);

  for (int c = 0; c < ncells_owned; ++c) {
    double vol = mesh_->cell_volume(c);

    Qmin = Qmax = 1.0;
    if (neg[c] != 0.0) 
      Qmin = vol * (bounds[0][c] - (*field_)[0][c] - dlo[c]) / neg[c];

    if (pos[c] != 0.0) 
      Qmax = vol * (bounds[1][c] - (*field_)[0][c] - dlo[c]) / pos[c];

    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int i = 0; i < nfaces; ++i) {
      int f = faces[i];
      if (flux_lo_f[0][f] > 0.0)
        alpha[f] = std::min(alpha[f], Qmax);
      else
        alpha[f] = std::min(alpha[f], Qmin);
    }
  }

  for (int f = 0; f < nfaces_owned; ++f) {
    flux_f[0][f] = flux_lo_f[0][f] + alpha[f] * (flux_ho_f[0][f] - flux_lo_f[0][f]);
  }
}

}  // namespace Operators
}  // namespace Amanzi


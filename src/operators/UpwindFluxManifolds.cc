/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Upwind a cell-centered field defined on a network of manifolds
  using a given face-based flux.
*/

// TPLs
#include "Epetra_IntVector.h"
#include "Teuchos_RCP.hpp"

// Operators
#include "UpwindFluxManifolds.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Init method is not used.
****************************************************************** */
void
UpwindFluxManifolds::Init(Teuchos::ParameterList& plist)
{
  method_ = Operators::OPERATOR_UPWIND_FLUX_MANIFOLDS;
  tolerance_ = plist.get<double>("tolerance", OPERATOR_UPWIND_RELATIVE_TOLERANCE);
}


/* ******************************************************************
* Upwind cells -> faces.
****************************************************************** */
void
UpwindFluxManifolds::Compute(const CompositeVector& flux,
                             const std::vector<int>& bc_model,
                             CompositeVector& field)
{
  AMANZI_ASSERT(field.hasComponent("cell"));
  AMANZI_ASSERT(field.hasComponent("face"));

  flux.scatterMasterToGhosted("face");
  field.scatterMasterToGhosted("cell");

  const auto& flux_f = *flux.viewComponent("face", true);
  const auto& field_c = *field.viewComponent("cell", true);
  const auto& field_bf = *field.viewComponent("boundary_face", true);
  auto& field_f = *field.viewComponent("face", true);

  double flxmin, flxmax, tol;
  flux_f.MinValue(&flxmin);
  flux_f.MaxValue(&flxmax);
  tol = tolerance_ * std::max(fabs(flxmin), fabs(flxmax));

  int nfaces_wghost = mesh_->getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::ALL);

  // multiple DOFs on faces require usage of block map
  const auto& fmap = *flux.ComponentMap("face", true);

  for (int f = 0; f < nfaces_wghost; ++f) {
    auto cells = mesh_->getFaceCells(f);
    int ncells = cells.size();

    int g = fmap.FirstPointInElement(f);

    // volume-weigthed average over downwind cells
    if (ncells > 1) {
      int dir;
      double uvol(0.0), umean(0.0), tvol(0.0), tmean(0.0), vol;

      for (int i = 0; i < ncells; ++i) {
        int c = cells[i];
        mesh_->getFaceNormal(f, c, &dir);
        vol = mesh_->getCellVolume(c);

        tvol += vol;
        tmean += vol * field_c[0][c];

        if (flux_f[0][g + i] * dir >= tol) {
          uvol += vol;
          umean += vol * field_c[0][c];
        }
      }

      // average upwind cell values and use them for all downwind faces
      if (uvol > 0.0) {
        umean /= uvol;
        for (int i = 0; i < ncells; ++i) {
          int c = cells[i];
          mesh_->getFaceNormal(f, c, &dir);
          if (flux_f[0][g + i] * dir >= tol)
            field_f[0][g + i] = field_c[0][c];
          else
            field_f[0][g + i] = umean;
        }
        // flow is negligent, use average all cell values
      } else {
        tmean /= tvol;
        for (int i = 0; i < ncells; ++i) field_f[0][g + i] = tmean;
      }

      // upwind only on inflow Dirichlet faces
    } else {
      field_f[0][g] = field_c[0][cells[0]];
      if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
        int bf = getFaceOnBoundaryBoundaryFace(*mesh_, f);
        field_f[0][g] = field_bf[0][bf];
      }
    }
  }
}

} // namespace Operators
} // namespace Amanzi

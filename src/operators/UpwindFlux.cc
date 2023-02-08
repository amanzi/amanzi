/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Upwind a cell-centered field (e.g. rel perm) using a given
  face-based flux (e.g. Darcy flux).
*/

// TPLs
#include "Epetra_IntVector.h"
#include "Teuchos_RCP.hpp"

// Operators
#include "UniqueLocalIndex.hh"
#include "UpwindFlux.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Public init method. It is not yet used.
****************************************************************** */
void
UpwindFlux::Init(Teuchos::ParameterList& plist)
{
  method_ = Operators::OPERATOR_UPWIND_FLUX;
  tolerance_ = plist.get<double>("tolerance", OPERATOR_UPWIND_RELATIVE_TOLERANCE);
  order_ = plist.get<int>("polynomial order", 1);
}


/* ******************************************************************
* Upwind field uses flux. The result is placed in field.
* Upwinded field must be calculated on all faces of the owned cells.
****************************************************************** */
void
UpwindFlux::Compute(const CompositeVector& flux,
                    const std::vector<int>& bc_model,
                    CompositeVector& field)
{
  AMANZI_ASSERT(field.HasComponent("cell"));
  AMANZI_ASSERT(field.HasComponent(face_comp_));

  flux.ScatterMasterToGhosted("face");
  field.ScatterMasterToGhosted("cell");

  const Epetra_MultiVector& flux_f = *flux.ViewComponent("face", true);

  const Epetra_MultiVector& field_c = *field.ViewComponent("cell", true);
  const Epetra_MultiVector& field_bf = *field.ViewComponent("boundary_face", true);
  Epetra_MultiVector& field_f = *field.ViewComponent(face_comp_, true);

  double flxmin, flxmax, tol;
  flux_f.MinValue(&flxmin);
  flux_f.MaxValue(&flxmax);
  tol = tolerance_ * std::max(fabs(flxmin), fabs(flxmax));

  int nfaces_wghost =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  // multiple DOFs on faces require usage of block map
  const auto& fmap = *flux.ComponentMap("face", true);

  int c1, c2, dir;
  double kc1, kc2;
  for (int f = 0; f < nfaces_wghost; ++f) {
    auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    int ncells = cells.size();

    int g = fmap.FirstPointInElement(f);
    int ndofs = fmap.ElementSize(f);

    c1 = cells[0];
    kc1 = field_c[0][c1];

    mesh_->getFaceNormal(f, c1, &dir);
    bool flag = (flux_f[0][g] * dir <= -tol); // upwind flag

    // average field on almost vertical faces
    if (ncells == 2 && ndofs == 1) {
      c2 = cells[1];
      kc2 = field_c[0][c2];

      if (fabs(flux_f[0][g]) <= tol) {
        double v1 = mesh_->getCellVolume(c1);
        double v2 = mesh_->getCellVolume(c2);

        double tmp = v2 / (v1 + v2);
        field_f[0][g] = kc1 * tmp + kc2 * (1.0 - tmp);
      } else {
        field_f[0][g] = (flag) ? kc2 : kc1;
      }

      // copy cell value on fractures (matrix only)
    } else if (ncells == 2 && ndofs == 2) {
      c2 = cells[1];
      kc2 = field_c[0][c2];

      int k = UniqueIndexFaceToCells(*mesh_, f, c1);
      field_f[0][g + k] = kc1;
      field_f[0][g + 1 - k] = kc2;

      // upwind only on inflow dirichlet faces
    } else {
      field_f[0][g] = kc1;
      if (bc_model[f] == OPERATOR_BC_DIRICHLET && flag) {
        int bf = getFaceOnBoundaryBoundaryFace(*mesh_, f);
        field_f[0][g] = field_bf[0][bf];
      }
    }
  }
}

} // namespace Operators
} // namespace Amanzi

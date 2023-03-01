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
  constant velocity (e.g. gravity).
*/

#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "Mesh.hh"
#include "Mesh_Algorithms.hh"
#include "Point.hh"

// Operators
#include "UpwindGravity.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Public init method. It is not yet used.
****************************************************************** */
void
UpwindGravity::Init(Teuchos::ParameterList& plist)
{
  method_ = Operators::OPERATOR_UPWIND_GRAVITY;
  tolerance_ = plist.get<double>("tolerance", OPERATOR_UPWIND_RELATIVE_TOLERANCE);
  order_ = plist.get<int>("polynomial", 1);

  int dim = mesh_->getSpaceDimension();
  g_[dim - 1] = -1.0;
}


/* ******************************************************************
* Upwind field uses gravity and places the result in field.
* Upwinded field must be calculated on all faces of the owned cells.
****************************************************************** */
void
UpwindGravity::Compute(const CompositeVector& flux,
                       const CompositeVector& solution,
                       const std::vector<int>& bc_model,
                       CompositeVector& field)
{
  AMANZI_ASSERT(field.HasComponent("cell"));
  AMANZI_ASSERT(field.HasComponent(face_comp_));

  field.ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& field_c = *field.ViewComponent("cell", true);
  const Epetra_MultiVector& field_bf = *field.ViewComponent("boundary_face", true);
  Epetra_MultiVector& field_f = *field.ViewComponent(face_comp_, true);

  int nfaces_wghost = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  int c1, c2, dir;
  double kc1, kc2;
  for (int f = 0; f < nfaces_wghost; ++f) {
    auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    int ncells = cells.size();

    c1 = cells[0];
    kc1 = field_c[0][c1];

    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f, c1, &dir);
    double flx_face = g_ * normal;
    bool flag = (flx_face <= -tolerance_); // upwind flag

    if (ncells == 2) {
      c2 = cells[1];
      kc2 = field_c[0][c2];

      // We average field on almost vertical faces.
      if (fabs(flx_face) <= tolerance_) {
        double v1 = mesh_->getCellVolume(c1);
        double v2 = mesh_->getCellVolume(c2);

        double tmp = v2 / (v1 + v2);
        field_f[0][f] = kc1 * tmp + kc2 * (1.0 - tmp);
      } else {
        field_f[0][f] = (flag) ? kc2 : kc1;
      }

      // We upwind only on inflow dirichlet faces.
    } else {
      field_f[0][f] = kc1;
      if (bc_model[f] == OPERATOR_BC_DIRICHLET && flag) {
        int bf = getFaceOnBoundaryBoundaryFace(*mesh_, f);
        field_f[0][f] = field_bf[0][bf];
      }
    }
  }
}

} // namespace Operators
} // namespace Amanzi

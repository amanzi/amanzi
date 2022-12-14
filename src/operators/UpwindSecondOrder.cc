/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

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

// Operators
#include "UpwindSecondOrder.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Public init method. It is not yet used.
****************************************************************** */
void
UpwindSecondOrder::Init(Teuchos::ParameterList& plist)
{
  method_ = Operators::OPERATOR_UPWIND_FLUX_SECOND_ORDER;
  tolerance_ = plist.get<double>("tolerance", OPERATOR_UPWIND_RELATIVE_TOLERANCE);
  order_ = plist.get<int>("polynomial order", 2);
}


/* ******************************************************************
* Flux-based upwind consistent with mimetic discretization.
****************************************************************** */
void
UpwindSecondOrder::Compute(const CompositeVector& flux,
                           const CompositeVector& solution,
                           const std::vector<int>& bc_model,
                           CompositeVector& field)
{
  AMANZI_ASSERT(field.HasComponent("cell"));
  AMANZI_ASSERT(field.HasComponent("grad"));
  AMANZI_ASSERT(field.HasComponent(face_comp_));

  field.ScatterMasterToGhosted("cell");
  flux.ScatterMasterToGhosted("face");

  const Epetra_MultiVector& flx_face = *flux.ViewComponent("face", true);
  // const Epetra_MultiVector& sol_face = *solution.ViewComponent("face", true);

  const Epetra_MultiVector& fld_cell = *field.ViewComponent("cell", true);
  const Epetra_MultiVector& fld_grad = *field.ViewComponent("grad", true);
  const Epetra_MultiVector& fld_boundary = *field.ViewComponent("boundary_face", true);
  const Epetra_Map& ext_face_map = mesh_->exterior_face_map(true);
  const Epetra_Map& face_map = mesh_->face_map(true);
  Epetra_MultiVector& upw_face = *field.ViewComponent(face_comp_, true);
  upw_face.PutScalar(0.0);

  double flxmin, flxmax;
  flx_face.MinValue(&flxmin);
  flx_face.MaxValue(&flxmax);
  double tol = tolerance_ * std::max(fabs(flxmin), fabs(flxmax));

  int dim = mesh_->space_dimension();
  AmanziGeometry::Point grad(dim);

  int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  for (int c = 0; c < ncells_wghost; c++) {
    const auto& faces = mesh_->cell_get_faces(c);
    const auto& dirs = mesh_->cell_get_face_dirs(c);
    int nfaces = faces.size();

    double kc(fld_cell[0][c]);
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    for (int i = 0; i < dim; i++) grad[i] = fld_grad[i][c];

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      bool flag = (flx_face[0][f] * dirs[n] <= -tol); // upwind flag

      // Internal faces. We average field on almost vertical faces.
      if (bc_model[f] == OPERATOR_BC_NONE && fabs(flx_face[0][f]) <= tol) {
        double tmp(0.5);
        int c2 = cell_get_face_adj_cell(*mesh_, c, f);
        if (c2 >= 0) {
          double v1 = mesh_->cell_volume(c);
          double v2 = mesh_->cell_volume(c2);
          tmp = v2 / (v1 + v2);
        }
        const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
        upw_face[0][f] += (kc + grad * (xf - xc)) * tmp;
        // Boundary faces. We upwind only on inflow dirichlet faces.
      } else if (bc_model[f] == OPERATOR_BC_DIRICHLET && flag) {
        upw_face[0][f] = fld_boundary[0][ext_face_map.LID(face_map.GID(f))];
      } else if (bc_model[f] == OPERATOR_BC_NEUMANN && flag) {
        upw_face[0][f] = kc;
      } else if (bc_model[f] == OPERATOR_BC_MIXED && flag) {
        upw_face[0][f] = kc;
        // Internal and boundary faces.
      } else if (!flag) {
        const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
        upw_face[0][f] = kc + grad * (xf - xc);
      }
    }
  }
}

} // namespace Operators
} // namespace Amanzi

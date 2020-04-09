/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_GRAVITY_HH_
#define AMANZI_GRAVITY_HH_

#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "Mesh.hh"
#include "Point.hh"

// Operators
#include "Upwind.hh"

namespace Amanzi {
namespace Operators {

template <class Model>
class UpwindGravity : public Upwind<Model> {
 public:
  UpwindGravity(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                Teuchos::RCP<const Model> model)
    : Upwind<Model>(mesh, model), g_(mesh->space_dimension()){};
  ~UpwindGravity(){};

  // main methods
  void Init(Teuchos::ParameterList& plist);

  void Compute(const CompositeVector& flux, const CompositeVector& solution,
               const std::vector<int>& bc_model, CompositeVector& field);

 private:
  using Upwind<Model>::mesh_;
  using Upwind<Model>::model_;
  using Upwind<Model>::face_comp_;

 private:
  int method_, order_;
  double tolerance_;
  AmanziGeometry::Point g_;
};


/* ******************************************************************
 * Public init method. It is not yet used.
 ****************************************************************** */
template <class Model>
void
UpwindGravity<Model>::Init(Teuchos::ParameterList& plist)
{
  method_ = Operators::OPERATOR_UPWIND_GRAVITY;
  tolerance_ =
    plist.get<double>("tolerance", OPERATOR_UPWIND_RELATIVE_TOLERANCE);
  order_ = plist.get<int>("polynomial", 1);

  int dim = mesh_->space_dimension();
  g_[dim - 1] = -1.0;
}


/* ******************************************************************
 * Upwind field uses gravity and places the result in field.
 * Upwinded field must be calculated on all faces of the owned cells.
 ****************************************************************** */
template <class Model>
void
UpwindGravity<Model>::Compute(const CompositeVector& flux,
                              const CompositeVector& solution,
                              const std::vector<int>& bc_model,
                              CompositeVector& field)
{
  AMANZI_ASSERT(field.HasComponent("cell"));
  AMANZI_ASSERT(field.HasComponent(face_comp_));

  field.ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& fld_cell = *field.ViewComponent("cell", true);
  const Epetra_MultiVector& fld_boundary =
    *field.ViewComponent("dirichlet_faces", true);
  const Epetra_Map& ext_face_map = mesh_->exterior_face_map(true);
  const Epetra_Map& face_map = mesh_->face_map(true);
  Epetra_MultiVector& upw_face = *field.ViewComponent(face_comp_, true);
  // const Epetra_MultiVector& sol_face = *solution.ViewComponent("face", true);

  int nfaces_wghost =
    mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  AmanziMesh::Entity_ID_List cells;

  int c1, c2, dir;
  double kc1, kc2;
  for (int f = 0; f < nfaces_wghost; ++f) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    c1 = cells[0];
    kc1 = fld_cell[0][c1];

    const AmanziGeometry::Point& normal =
      mesh_->face_normal(f, false, c1, &dir);
    double flx_face = g_ * normal;
    bool flag = (flx_face <= -tolerance_); // upwind flag

    if (ncells == 2) {
      c2 = cells[1];
      kc2 = fld_cell[0][c2];

      // We average field on almost vertical faces.
      if (fabs(flx_face) <= tolerance_) {
        double v1 = mesh_->cell_volume(c1, false);
        double v2 = mesh_->cell_volume(c2, false);

        double tmp = v2 / (v1 + v2);
        upw_face[0][f] = kc1 * tmp + kc2 * (1.0 - tmp);
      } else {
        upw_face[0][f] = (flag) ? kc2 : kc1;
      }

      // We upwind only on inflow dirichlet faces.
    } else {
      upw_face[0][f] = kc1;
      if (bc_model[f] == OPERATOR_BC_DIRICHLET && flag) {
        upw_face[0][f] = fld_boundary[0][ext_face_map.LID(face_map.GID(f))];
      }
      // if (bc_model[f] == OPERATOR_BC_NEUMANN) {
      //   upw_face[0][f] = ((*model_).*Value)(c, sol_face[0][f]);
      // }
    }
  }
}

} // namespace Operators
} // namespace Amanzi

#endif

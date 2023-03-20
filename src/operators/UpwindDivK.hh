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

#ifndef AMANZI_UPWIND_DIVK_HH_
#define AMANZI_UPWIND_DIVK_HH_

#include <string>
#include <vector>

#include "Epetra_IntVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "CompositeVector.hh"
#include "Mesh.hh"
#include "Mesh_Algorithms.hh"
#include "VerboseObject.hh"

#include "Upwind.hh"

namespace Amanzi {
namespace Operators {

class UpwindDivK : public Upwind {
 public:
  UpwindDivK(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : Upwind(mesh){};
  ~UpwindDivK(){};

  // main methods
  void Init(Teuchos::ParameterList& plist);

  void
  Compute(const CompositeVector& flux, const std::vector<int>& bc_model, CompositeVector& field);

 private:
  int method_, order_;
  double tolerance_;
};


/* ******************************************************************
* Public init method. It is not yet used.
****************************************************************** */
inline void
UpwindDivK::Init(Teuchos::ParameterList& plist)
{
  method_ = Operators::OPERATOR_UPWIND_DIVK;
  tolerance_ = plist.get<double>("tolerance", OPERATOR_UPWIND_RELATIVE_TOLERANCE);
  order_ = plist.get<int>("polynomial order", 1);
}


/* ******************************************************************
* Flux-based upwind consistent with mimetic discretization.
****************************************************************** */
inline void
UpwindDivK::Compute(const CompositeVector& flux,
                    const std::vector<int>& bc_model,
                    CompositeVector& field)
{
  AMANZI_ASSERT(field.HasComponent("cell"));
  AMANZI_ASSERT(field.HasComponent(face_comp_));

  field.ScatterMasterToGhosted("cell");
  flux.ScatterMasterToGhosted("face");

  const Epetra_MultiVector& flx_face = *flux.ViewComponent("face", true);

  const Epetra_MultiVector& fld_cell = *field.ViewComponent("cell", true);
  const Epetra_MultiVector& fld_boundary = *field.ViewComponent("boundary_face", true);
  const Epetra_Map& ext_face_map = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,true);
  const Epetra_Map& face_map = mesh_->getMap(AmanziMesh::Entity_kind::FACE,true);
  Epetra_MultiVector& upw_face = *field.ViewComponent(face_comp_, true);
  upw_face.PutScalar(0.0);

  double flxmin, flxmax;
  flx_face.MinValue(&flxmin);
  flx_face.MaxValue(&flxmax);
  double tol = tolerance_ * std::max(fabs(flxmin), fabs(flxmax));

  int ncells_wghost = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  for (int c = 0; c < ncells_wghost; c++) {
    auto [faces, dirs] = mesh_->getCellFacesAndDirections(c);
    int nfaces = faces.size();
    double kc(fld_cell[0][c]);

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      bool flag = (flx_face[0][f] * dirs[n] <= -tol); // upwind flag

      // Internal faces. We average field on almost vertical faces.
      if (bc_model[f] == OPERATOR_BC_NONE && fabs(flx_face[0][f]) <= tol) {
        double tmp(0.5);
        int c2 = AmanziMesh::MeshAlgorithms::getFaceAdjacentCell(*mesh_, c, f);
        if (c2 >= 0) {
          double v1 = mesh_->getCellVolume(c);
          double v2 = mesh_->getCellVolume(c2);
          tmp = v2 / (v1 + v2);
        }
        upw_face[0][f] += kc * tmp;
        // Boundary faces. We upwind only on inflow dirichlet faces.
      } else if (bc_model[f] == OPERATOR_BC_DIRICHLET && flag) {
        upw_face[0][f] = fld_boundary[0][ext_face_map.LID(face_map.GID(f))];
      } else if (bc_model[f] == OPERATOR_BC_NEUMANN && flag) {
        upw_face[0][f] = kc;
      } else if (bc_model[f] == OPERATOR_BC_MIXED && flag) {
        upw_face[0][f] = kc;
        // Internal and boundary faces.
      } else if (!flag) {
        int c2 = AmanziMesh::MeshAlgorithms::getFaceAdjacentCell(*mesh_, c, f);
        if (c2 >= 0) {
          double kc2(fld_cell[0][c2]);
          upw_face[0][f] = std::pow(kc * (kc + kc2) / 2, 0.5);
        } else {
          upw_face[0][f] = kc;
        }
      }
    }
  }
}

} // namespace Operators
} // namespace Amanzi

#endif

/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

// -----------------------------------------------------------------------------
// ATS
//
// License: see $ATS_DIR/COPYRIGHT
// Author: Ethan Coon (ecoon@lanl.gov)
//
// Scheme for taking coefficients for div-grad operators from cells to
// faces.
// -----------------------------------------------------------------------------

#include "tensor.hpp"
#include "upwind_gravity_flux.hh"

namespace Amanzi {
namespace Operators {

UpwindGravityFlux::UpwindGravityFlux(std::string pkname,
        std::string cell_coef,
        std::string face_coef,
        const Teuchos::RCP<std::vector<WhetStone::Tensor> > K) :
    pkname_(pkname),
    cell_coef_(cell_coef),
    face_coef_(face_coef),
    K_(K) {};


void UpwindGravityFlux::Update(const Teuchos::Ptr<State>& S) {
  Teuchos::RCP<const CompositeVector> cell = S->GetFieldData(cell_coef_);
  Teuchos::RCP<const Epetra_Vector> g_vec = S->GetConstantVectorData("gravity");
  Teuchos::RCP<CompositeVector> face = S->GetFieldData(face_coef_, pkname_);
  CalculateCoefficientsOnFaces(*cell, *g_vec, face.ptr());
};


void UpwindGravityFlux::CalculateCoefficientsOnFaces(
        const CompositeVector& cell_coef,
        const Epetra_Vector& g_vec,
        const Teuchos::Ptr<CompositeVector>& face_coef) {

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  double flow_eps = 1.e-10;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = face_coef->mesh();

  // set up gravity
  AmanziGeometry::Point gravity(g_vec.MyLength());
  for (int i=0; i!=g_vec.MyLength(); ++i) gravity[i] = g_vec[i];

  face_coef->ViewComponent("cell")->PutScalar(1.0);

  // communicate cell values to ghosted
  cell_coef.ScatterMasterToGhosted("cell");

  int c_used = cell_coef.size("cell", true);
  for (int c=0; c!=c_used; ++c) {
    mesh->cell_get_faces_and_dirs(c, &faces, &dirs);
    AmanziGeometry::Point Kgravity = (*K_)[c] * gravity;

    for (int n=0; n!=faces.size(); ++n) {
      int f = faces[n];

      const AmanziGeometry::Point& normal = mesh->face_normal(f);

      if ((normal * Kgravity) * dirs[n] >= flow_eps) {
        (*face_coef)("face",f) = cell_coef("cell",c);
      } else if (std::abs((normal * Kgravity) * dirs[n]) < flow_eps) {
        (*face_coef)("face",f) += cell_coef("cell",c) / 2.0;
      }
    }
  }
};


} //namespace
} //namespace

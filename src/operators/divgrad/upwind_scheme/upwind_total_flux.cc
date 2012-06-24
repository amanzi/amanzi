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

#include "upwind_total_flux.hh"

namespace Amanzi {
namespace Operators {

UpwindTotalFlux::UpwindTotalFlux(std::string pkname,
        std::string cell_coef,
        std::string face_coef,
        std::string flux) :
    pkname_(pkname),
    cell_coef_(cell_coef),
    face_coef_(face_coef),
    flux_(flux) {};


void UpwindTotalFlux::Update(const Teuchos::Ptr<State>& S) {
  Teuchos::RCP<const CompositeVector> cell = S->GetFieldData(cell_coef_);
  Teuchos::RCP<const CompositeVector> flux = S->GetFieldData(flux_);
  Teuchos::RCP<CompositeVector> face = S->GetFieldData(face_coef_, pkname_);
  CalculateCoefficientsOnFaces(*cell, *flux, face.ptr());
};


void UpwindTotalFlux::CalculateCoefficientsOnFaces(
        const CompositeVector& cell_coef,
        const CompositeVector& flux,
        const Teuchos::Ptr<CompositeVector>& face_coef) {

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = face_coef->mesh();
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  double flow_eps = 1.e-10;

  face_coef->ViewComponent("cell")->PutScalar(1.0);

  cell_coef.ScatterMasterToGhosted("cell");

  int c_used = cell_coef.size("cell", true);
  for (int c=0; c!=c_used; ++c) {
    mesh->cell_get_faces_and_dirs(c, &faces, &dirs);

    for (int n=0; n!=faces.size(); ++n) {
      int f = faces[n];

      if ((flux("face",f) * dirs[n] >= flow_eps)) {
        (*face_coef)("face",f) = cell_coef("cell",c);
      } else if (std::abs(flux("face",f)) < flow_eps) {
        (*face_coef)("face",f) += cell_coef("cell",c) / 2.0;
      }
    }
  }
};


} //namespace
} //namespace

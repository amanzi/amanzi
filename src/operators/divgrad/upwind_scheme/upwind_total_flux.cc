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

  // initialize the face coefficients
  face_coef->ViewComponent("face",true)->PutScalar(0.0);
  if (face_coef->has_component("cell")) {
    face_coef->ViewComponent("cell",true)->PutScalar(1.0);
  }

  // Note that by scattering, and then looping over all USED cells, we
  // end up getting the correct upwind values in all faces (owned or
  // not) bordering an owned cell.  This is the necessary data for
  // making the local matrices in MFD, so there is no need to
  // communicate the resulting face coeficients.

  // communicate ghosted cells
  cell_coef.ScatterMasterToGhosted("cell");
  flux.ScatterMasterToGhosted("face");

  for (int c=0; c!=cell_coef.size("cell", true); ++c) {
    mesh->cell_get_faces_and_dirs(c, &faces, &dirs);
    for (int n=0; n!=faces.size(); ++n) {
      int f = faces[n];

      if ((flux("face",f) * dirs[n] >= flow_eps)) {
        (*face_coef)("face",f) = cell_coef("cell",c);
      } else if (std::abs(flux("face",f) * dirs[n]) < flow_eps) {
        (*face_coef)("face",f) += cell_coef("cell",c) / 2.0;
      }

      //      std::cout << "UPWIND (" << c << "," << f << "): flux = " << flux("face",f) << " dir = " << dirs[n] << " cell coef = " << cell_coef("cell",c) << " face coef = " << (*face_coef)("face",f) << std::endl;

    }
  }
};


} //namespace
} //namespace

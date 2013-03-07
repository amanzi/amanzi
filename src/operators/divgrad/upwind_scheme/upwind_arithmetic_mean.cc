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

#include "composite_vector.hh"
#include "state.hh"
#include "upwind_arithmetic_mean.hh"

namespace Amanzi {
namespace Operators {

UpwindArithmeticMean::UpwindArithmeticMean(std::string pkname,
        std::string cell_coef,
        std::string face_coef) :
    pkname_(pkname),
    cell_coef_(cell_coef),
    face_coef_(face_coef) {};


void UpwindArithmeticMean::Update(const Teuchos::Ptr<State>& S) {
  Teuchos::RCP<const CompositeVector> cell = S->GetFieldData(cell_coef_);
  Teuchos::RCP<CompositeVector> face = S->GetFieldData(face_coef_, pkname_);
  CalculateCoefficientsOnFaces(*cell, face.ptr());
};


void UpwindArithmeticMean::CalculateCoefficientsOnFaces(
        const CompositeVector& cell_coef,
        const Teuchos::Ptr<CompositeVector>& face_coef) {

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = face_coef->mesh();
  AmanziMesh::Entity_ID_List faces;

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

  Epetra_MultiVector& face_coef_f = *face_coef->ViewComponent("face",true);
  const Epetra_MultiVector& cell_coef_c = *cell_coef.ViewComponent("cell",true);

  int c_used = cell_coef.size("cell", true);
  for (int c=0; c!=c_used; ++c) {
    std::vector<int> fdirs;
    mesh->cell_get_faces_and_dirs(c, &faces, &fdirs);

    for (int n=0; n!=faces.size(); ++n) {
      int f = faces[n];
      face_coef_f[0][f] += cell_coef_c[0][c] / 2.0;
    }
  }
};


} //namespace
} //namespace

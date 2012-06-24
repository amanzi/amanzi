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

#include "Mesh.hh"

#include "upwind_cell_centered.hh"

namespace Amanzi {
namespace Operators {

UpwindCellCentered::UpwindCellCentered(std::string pkname,
        std::string cell_coef,
        std::string face_coef) :
    pkname_(pkname),
    cell_coef_(cell_coef),
    face_coef_(face_coef) {};


void UpwindCellCentered::Update(const Teuchos::Ptr<State>& S) {
  Teuchos::RCP<const CompositeVector> cell = S->GetFieldData(cell_coef_);
  Teuchos::RCP<CompositeVector> face = S->GetFieldData(face_coef_, pkname_);
  CalculateCoefficientsOnFaces(*cell, face.ptr());
};


void UpwindCellCentered::CalculateCoefficientsOnFaces(
        const CompositeVector& cell_coef,
        const Teuchos::Ptr<CompositeVector>& face_coef) {

  *face_coef->ViewComponent("cell") = *cell_coef.ViewComponent("cell");
  if (face_coef->has_component("face")) {
    face_coef->ViewComponent("face")->PutScalar(1.0);
  }
};


} //namespace
} //namespace

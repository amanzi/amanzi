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
#include "composite_vector.hh"
#include "state.hh"
#include "upwind_total_flux.hh"
#include "Epetra_IntVector.h"

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

  // initialize the face coefficients
  face_coef->ViewComponent("face",true)->PutScalar(0.0);
  if (face_coef->has_component("cell")) {
    face_coef->ViewComponent("cell",true)->PutScalar(1.0);
  }

  // communicate ghosted cells
  cell_coef.ScatterMasterToGhosted("cell");
  flux.ScatterMasterToGhosted("face");

  // identify upwind/downwind cells for each face
  Epetra_IntVector upwind_cell(*face_coef->map("face",true));
  upwind_cell.PutValue(-1);
  Epetra_IntVector downwind_cell(*face_coef->map("face",true));
  downwind_cell.PutValue(-1);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;
  for (int c=0; c!=cell_coef.size("cell",false); ++c) {
    mesh->cell_get_faces_and_dirs(c, &faces, &fdirs);

    for (int n=0; n!=faces.size(); ++n) {
      int f = faces[n];
      if (flux("face",f) * fdirs[n] >= 0) {
        upwind_cell[f] = c;
      } else {
        downwind_cell[f] = c;
      }
    }
  }

  // Determine the face coefficient
  AmanziMesh::Entity_ID_List cells;
  double eps = 1.e-15;
  double flow_eps_factor = 1.e-8;

  for (int f=0; f!=face_coef->size("face",false); ++f) {
    mesh->face_get_cells(f, AmanziMesh::USED, &cells);

    if (cells.size() == 1) {
      (*face_coef)("face",f) = cell_coef("cell",cells[0]);
    } else {
      // Determine the size of the overlap region, a smooth transition region
      // near zero flux
      double flow_eps = 0.0;
      if ((cell_coef("cell",cells[0]) > 0) || (cell_coef("cell",cells[1]) > 0)) {
        flow_eps = (cell_coef("cell",cells[0]) * cell_coef("cell",cells[1]))
            / (cell_coef("cell",cells[0]) + cell_coef("cell",cells[1]));
      }

      // these parameters may be key!
      flow_eps = std::max(flow_eps*flow_eps_factor, eps);

      // Determine the coefficient
      if (abs(flux("face",f)) >= flow_eps) {
        (*face_coef)("face",f) = cell_coef("cell", upwind_cell[f]);
      } else {
        // Parameterization of a linear scaling between upwind and downwind.
        double param;
        if (flow_eps < 2*eps) {
          param = 0.5;
        } else {
          param = abs(flux("face",f)) / (2*flow_eps) + 0.5;
        }
        ASSERT(param >= 0.0);
        ASSERT(param <= 1.0);
        (*face_coef)("face",f) = cell_coef("cell", upwind_cell[f]) * param
            + cell_coef("cell", downwind_cell[f]) * (1. - param);
      }
    }
  }
};

} //namespace
} //namespace

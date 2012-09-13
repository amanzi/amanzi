/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

// -----------------------------------------------------------------------------
// ATS
//
// License: see $ATS_DIR/COPYRIGHT
// Author: Ethan Coon (ecoon@lanl.gov)
//
// Scheme for taking coefficients for div-grad operators from cells to faces.
// Upwinds based upon a potential vector, with an overlap region size
// determined by an (optional) other field.
// -----------------------------------------------------------------------------

#include "composite_vector.hh"
#include "state.hh"
#include "upwind_potential_difference.hh"

namespace Amanzi {
namespace Operators {

UpwindPotentialDifference::UpwindPotentialDifference(std::string pkname,
        std::string cell_coef,
        std::string face_coef,
        std::string potential,
        std::string overlap) :
    pkname_(pkname),
    cell_coef_(cell_coef),
    face_coef_(face_coef),
    potential_(potential),
    overlap_(overlap) {
  if (overlap_ == std::string("")) {
    overlap_ = potential_;
  }
};


void UpwindPotentialDifference::Update(const Teuchos::Ptr<State>& S) {
  Teuchos::RCP<const CompositeVector> cell = S->GetFieldData(cell_coef_);
  Teuchos::RCP<const CompositeVector> potential = S->GetFieldData(potential_);
  Teuchos::RCP<const CompositeVector> overlap = S->GetFieldData(overlap_);
  Teuchos::RCP<CompositeVector> face = S->GetFieldData(face_coef_, pkname_);
  CalculateCoefficientsOnFaces(*cell, *potential, *overlap, face.ptr());
};


void UpwindPotentialDifference::CalculateCoefficientsOnFaces(
        const CompositeVector& cell_coef,
        const CompositeVector& potential,
        const CompositeVector& overlap,
        const Teuchos::Ptr<CompositeVector>& face_coef) {

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = face_coef->mesh();
  AmanziMesh::Entity_ID_List cells;
  std::vector<int> dirs;
  double eps = 1.e-16;

  // communicate ghosted cells
  cell_coef.ScatterMasterToGhosted("cell");
  potential.ScatterMasterToGhosted("cell");
  overlap.ScatterMasterToGhosted("cell");

  for (int f=0; f!=face_coef->size("face",false); ++f) {
    mesh->face_get_cells(f, AmanziMesh::USED, &cells);

    if (cells.size() == 1) {
      if (potential("cell", cells[0]) >= potential("face", f)) {
        (*face_coef)("face",f) = cell_coef("cell", cells[0]);
      }
    } else {
      // Determine the size of the overlap region, a smooth transition region
      // near zero potential difference.
      double flow_eps = 0.0;
      if ((overlap("cell",cells[0]) > 0) || (overlap("cell",cells[1]) > 0)) {
        flow_eps = (overlap("cell",cells[0]) * overlap("cell",cells[1]))
            / (overlap("cell",cells[0]) + overlap("cell",cells[1]));
      }
      flow_eps = std::max(flow_eps, eps);

      // Determine the coefficient.
      if (potential("cell",cells[0]) - potential("cell", cells[1]) > flow_eps) {
        (*face_coef)("face",f) = cell_coef("cell", cells[0]);
      } else if (potential("cell",cells[1]) - potential("cell", cells[0]) > flow_eps) {
        (*face_coef)("face",f) = cell_coef("cell", cells[1]);
      } else {
        // Parameterization of a linear scaling between upwind and downwind.
        double param;
        if (flow_eps < 2*eps) {
          param = 0.5;
        } else {
          param = (potential("cell",cells[1]) - potential("cell", cells[0]))
              / (2*flow_eps) + 0.5;
        }
        ASSERT(param >= 0.0);
        ASSERT(param <= 1.0);
        (*face_coef)("face",f) = cell_coef("cell", cells[1]) * param
            + cell_coef("cell", cells[0]) * (1. - param);
      }

      // if ((0.0 < (*face_coef)("face",f)) && ((*face_coef)("face",f) < 1.0e-40)) {
      //   std::cout << "Near null space rel perm calculated: face = " << f << ", rel perm = " <<  (*face_coef)("face",f) << std::endl;
      //   std::cout << "  pres+elev: " << potential("cell", cells[0]) << "  "
      //             << potential("cell", cells[0]) << std::endl;
      //   std::cout << "  cell perms: " << cell_coef("cell", cells[0]) << "  "
      //             << cell_coef("cell", cells[0]) << std::endl;
      //   std::cout << "  flow_eps: " << flow_eps << std::endl;
      //   std::cout << "  if branches: " <<  (potential("cell",cells[0]) - potential("cell", cells[1]) > flow_eps)
      //             << (potential("cell",cells[1]) - potential("cell", cells[0]) > flow_eps)
      //             << (flow_eps == 0.0)
      //             << (potential("cell",cells[1]) - potential("cell", cells[0]))
      //       / (2*flow_eps) + 0.5 << std::endl;

      // }
    }
  }
};


} //namespace
} //namespace

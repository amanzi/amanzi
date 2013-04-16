/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  An evaluator for pulling the darcy flux, at the surface, from the
  subsurface field and putting it into a surface field.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Mesh_MSTK.hh"
#include "matrix_mfd.hh"
#include "nonlinear_source_from_subsurface_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

NonlinearSourceFromSubsurfaceEvaluator::NonlinearSourceFromSubsurfaceEvaluator(
        Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  my_key_ = plist_.get<std::string>("source key",
          "overland_source_from_subsurface");
  setLinePrefix(my_key_+std::string(" evaluator"));

  height_key_ = plist_.get<std::string>("subsurface height key", "ponded_depth");
  dependencies_.insert(height_key_);

  density_key_ = plist_.get<std::string>("molar density key", "molar_density_liquid");
  dependencies_.insert(density_key_);

  pressure_key_ = plist_.get<std::string>("subsurface pressure key", "pressure");
  dependencies_.insert(pressure_key_);

  surface_mesh_key_ = plist_.get<std::string>("surface mesh key", "surface");
  subsurface_mesh_key_ = plist_.get<std::string>("subsurface mesh key", "domain");
}

NonlinearSourceFromSubsurfaceEvaluator::NonlinearSourceFromSubsurfaceEvaluator(
        const NonlinearSourceFromSubsurfaceEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    height_key_(other.height_key_),
    density_key_(other.density_key_),
    pressure_key_(other.pressure_key_),
    surface_mesh_key_(other.surface_mesh_key_),
    subsurface_mesh_key_(other.subsurface_mesh_key_),
    face_and_dirs_(other.face_and_dirs_),
    op_(other.op_) {}

Teuchos::RCP<FieldEvaluator> NonlinearSourceFromSubsurfaceEvaluator::Clone() const {
  return Teuchos::rcp(new NonlinearSourceFromSubsurfaceEvaluator(*this));
}


void NonlinearSourceFromSubsurfaceEvaluator::IdentifyFaceAndDirection_(
        const Teuchos::Ptr<State>& S) {
  // grab the meshes
  Teuchos::RCP<const AmanziMesh::Mesh> subsurface = S->GetMesh(subsurface_mesh_key_);
  Teuchos::RCP<const AmanziMesh::Mesh_MSTK> surface =
      Teuchos::rcp_static_cast<const AmanziMesh::Mesh_MSTK>(S->GetMesh("surface"));

  // allocate space for face IDs and directions
  int ncells = surface->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  face_and_dirs_ = Teuchos::rcp(new std::vector<FaceDir>(ncells));

  for (int c=0; c!=ncells; ++c) {
    // Get the face on the subsurface mesh corresponding to the cell
    // of the surface mesh.
    AmanziMesh::Entity_ID domain_face =
      surface->entity_get_parent(AmanziMesh::CELL, c);

    // Get the direction corresponding to that face wrt its only cell.
    // -- get the cell
    AmanziMesh::Entity_ID_List cells;
    subsurface->face_get_cells(domain_face, AmanziMesh::OWNED, &cells);
    ASSERT(cells.size() == 1);

    // -- Get directions
    AmanziMesh::Entity_ID_List faces;
    std::vector<int> fdirs;
    subsurface->cell_get_faces_and_dirs(cells[0], &faces, &fdirs);
    int index = std::find(faces.begin(), faces.end(), domain_face) - faces.begin();
    int direction = fdirs[index];

    // Put (face,dir) into cached data.
    (*face_and_dirs_)[c] = std::make_pair(domain_face, direction);
  }
}

// Required methods from SecondaryVariableFieldEvaluator
void
NonlinearSourceFromSubsurfaceEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  if (face_and_dirs_ == Teuchos::null) {
    IdentifyFaceAndDirection_(S);
  }

  double eps = 1.e-5;
  Teuchos::RCP<const AmanziMesh::Mesh> subsurface = S->GetMesh(subsurface_mesh_key_);

  Teuchos::RCP<const CompositeVector> height = S->GetFieldData(height_key_);
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pressure_key_);
  Teuchos::RCP<const CompositeVector> n_liq = S->GetFieldData(density_key_);
  Teuchos::RCP<const CompositeVector> rho_liq = S->GetFieldData("mass_density_liquid");
  Teuchos::RCP<const double> p_atm = S->GetScalarData("atmospheric_pressure");
  Teuchos::RCP<const Epetra_Vector> gvec = S->GetConstantVectorData("gravity");

  std::vector<Teuchos::SerialDenseMatrix<int, double> > Aff_cells = op_->Aff_cells();

  for (int c=0; c!=result->size("cell",false); ++c) {
    // subsurface face
    int ss_f = (*face_and_dirs_)[c].first;

    // subsurface cell
    AmanziMesh::Entity_ID_List cells;
    subsurface->face_get_cells(ss_f, AmanziMesh::OWNED, &cells);
    ASSERT(cells.size() == 1);
    int ss_c = cells[0];

    // calculate the surface pressure given height
    double surface_p = (*p_atm) - (*rho_liq)("cell",ss_c) * (*gvec)[2] * (*height)("cell",c);

    // calculate the flux out the surface
    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;
    subsurface->cell_get_faces_and_dirs(ss_c, &faces, &dirs);

    int ss_nfaces = faces.size();
    std::vector<double> dp(ss_nfaces);
    for (int n=0; n!=ss_nfaces; ++n) {
      int f = faces[n];
      if (f == (*face_and_dirs_)[c].first) {
        dp[n] = surface_p - (*pres)("cell", ss_c);
      } else {
        dp[n] = (*pres)("face",f) - (*pres)("cell",ss_c);
      }
    }

    double s = 0.0;
    for (int m=0; m!=ss_nfaces; ++m) s += Aff_cells[ss_c](ss_f,m)*dp[m];

    // note sign: positive direction = outward height, so positive
    // source to surface.  Height is in mol/(m^2-s), and our source is
    // in m/s.
    (*result)("cell",c) = (*face_and_dirs_)[c].second * s / (*n_liq)("cell", ss_c);
  }
}

void
NonlinearSourceFromSubsurfaceEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  ASSERT(0);
  // this would require differentiating height wrt pressure, which we
  // don't do for now.
}

} // namespace
} // namespace
} // namespace

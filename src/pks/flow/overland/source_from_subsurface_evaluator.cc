/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  An evaluator for pulling the darcy flux, at the surface, from the
  subsurface field and putting it into a surface field.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Mesh_MSTK.hh"
#include "source_from_subsurface_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

SourceFromSubsurfaceEvaluator::SourceFromSubsurfaceEvaluator(
        Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  my_key_ = plist_.get<std::string>("source key",
          "overland_source_from_subsurface");
  setLinePrefix(my_key_+std::string(" evaluator"));

  flux_key_ = plist_.get<std::string>("subsurface flux key", "darcy_flux");
  dependencies_.insert(flux_key_);

  density_key_ = plist_.get<std::string>("molar density key", "molar_density_liquid");
  dependencies_.insert(density_key_);

  surface_mesh_key_ = plist_.get<std::string>("surface mesh key", "surface");
  subsurface_mesh_key_ = plist_.get<std::string>("subsurface mesh key", "domain");
}

SourceFromSubsurfaceEvaluator::SourceFromSubsurfaceEvaluator(
        const SourceFromSubsurfaceEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    flux_key_(other.flux_key_),
    density_key_(other.density_key_),
    surface_mesh_key_(other.surface_mesh_key_),
    subsurface_mesh_key_(other.subsurface_mesh_key_),
    face_and_dirs_(other.face_and_dirs_) {}

Teuchos::RCP<FieldEvaluator> SourceFromSubsurfaceEvaluator::Clone() const {
  return Teuchos::rcp(new SourceFromSubsurfaceEvaluator(*this));
}


void SourceFromSubsurfaceEvaluator::IdentifyFaceAndDirection_(
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
void SourceFromSubsurfaceEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  if (face_and_dirs_ == Teuchos::null) {
    IdentifyFaceAndDirection_(S);
  }

  Teuchos::RCP<const AmanziMesh::Mesh> subsurface = S->GetMesh(subsurface_mesh_key_);
  Teuchos::RCP<const CompositeVector> flux = S->GetFieldData(flux_key_);
  Teuchos::RCP<const CompositeVector> n_liq = S->GetFieldData(density_key_);


  for (int c=0; c!=result->size("cell",false); ++c) {
    // this will hopefully change to be the density of the surface water?
    AmanziMesh::Entity_ID_List cells;
    subsurface->face_get_cells((*face_and_dirs_)[c].first, AmanziMesh::OWNED, &cells);
    ASSERT(cells.size() == 1);

    // note sign: positive direction = outward flux, so positive
    // source to surface.  Flux is in mol/(m^2-s), and our source is
    // in m/s.
    (*result)("cell",c) = (*face_and_dirs_)[c].second *
      (*flux)("face", (*face_and_dirs_)[c].first) / (*n_liq)("cell", cells[0]);
  }
}

void SourceFromSubsurfaceEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  ASSERT(0);
  // this would require differentiating flux wrt pressure, which we
  // don't do for now.
}

} // namespace
} // namespace
} // namespace

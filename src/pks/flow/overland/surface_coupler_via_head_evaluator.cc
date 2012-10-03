/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  An evaluator for pulling the darcy flux, at the surface, from the
  subsurface field and putting it into a surface field.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Mesh_MSTK.hh"
#include "surface_coupler_via_head_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

Utils::RegisteredFactory<FieldEvaluator,SurfaceCouplerViaHeadEvaluator>
SurfaceCouplerViaHeadEvaluator::fac_("surface_coupler_via_head");

SurfaceCouplerViaHeadEvaluator::SurfaceCouplerViaHeadEvaluator(
        Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  my_key_ = plist_.get<std::string>("source key",
          "overland_source_from_subsurface");
  setLinePrefix(my_key_+std::string(" evaluator"));

  flux_key_ = plist_.get<std::string>("flux key", "darcy_flux");

  // since we cannot have flux as a dependency (it has no model), we have to
  // use pressure as a proxy.
  Key pres_key = plist_.get<std::string>("pressure key", "pressure");
  dependencies_.insert(pres_key);

  dens_key_ = plist_.get<std::string>("molar density key", "molar_density_liquid");
  dependencies_.insert(dens_key_);

  surface_mesh_key_ = plist_.get<std::string>("surface mesh key", "surface");
  subsurface_mesh_key_ = plist_.get<std::string>("subsurface mesh key", "domain");
}

SurfaceCouplerViaHeadEvaluator::SurfaceCouplerViaHeadEvaluator(
        const SurfaceCouplerViaHeadEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    flux_key_(other.flux_key_),
    dens_key_(other.dens_key_),
    surface_mesh_key_(other.surface_mesh_key_),
    subsurface_mesh_key_(other.subsurface_mesh_key_),
    face_and_dirs_(other.face_and_dirs_) {}

Teuchos::RCP<FieldEvaluator> SurfaceCouplerViaHeadEvaluator::Clone() const {
  return Teuchos::rcp(new SurfaceCouplerViaHeadEvaluator(*this));
}


void SurfaceCouplerViaHeadEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // for now just passing... might do something later here?
  S->RequireField(my_key_, my_key_);
  S->RequireField(flux_key_);
}


void SurfaceCouplerViaHeadEvaluator::IdentifyFaceAndDirection_(
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
    // note this direction does not include face surface area!
    int direction = fdirs[index] > 0 ? 1 : -1;

    // Put (face,dir) into cached data.
    (*face_and_dirs_)[c] = std::make_pair(domain_face, direction);
  }
}

// Required methods from SecondaryVariableFieldEvaluator
void SurfaceCouplerViaHeadEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  if (face_and_dirs_ == Teuchos::null) {
    IdentifyFaceAndDirection_(S);
  }

  Teuchos::RCP<const AmanziMesh::Mesh> subsurface = S->GetMesh(subsurface_mesh_key_);
  const Epetra_MultiVector& flux = *S->GetFieldData(flux_key_)->ViewComponent("face",false);
  const Epetra_MultiVector& dens = *S->GetFieldData(dens_key_)->ViewComponent("cell",false);

  int ncells = result->size("cell",false);
  for (int c=0; c!=ncells; ++c) {
    AmanziMesh::Entity_ID_List cells;
    subsurface->face_get_cells((*face_and_dirs_)[c].first, AmanziMesh::OWNED, &cells);
    ASSERT(cells.size() == 1);

    // upwind head in the "rel perm" term, k = head*coef
    (*result)("cell",c) = flux[0][(*face_and_dirs_)[c].first]*(*face_and_dirs_)[c].second
        / dens[0][cells[0]];

    Teuchos::OSTab tab = getOSTab();
    if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
      *out_ << "Source, sub -> surface = " << (*result)("cell",c) << std::endl;
    }

  }
}

void SurfaceCouplerViaHeadEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  ASSERT(0);
  // this would require differentiating flux wrt pressure, which we
  // don't do for now.
}

} // namespace
} // namespace
} // namespace

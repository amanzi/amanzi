/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  An evaluator for pulling the darcy flux, at the surface, from the
  subsurface field and putting it into a surface field.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Mesh_MSTK.hh"
#include "surface_coupler_via_source_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

Utils::RegisteredFactory<FieldEvaluator,SurfaceCouplerViaSourceEvaluator>
SurfaceCouplerViaSourceEvaluator::fac_("surface_coupler_via_source");

SurfaceCouplerViaSourceEvaluator::SurfaceCouplerViaSourceEvaluator(
        Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  my_key_ = plist_.get<std::string>("source key",
          "overland_source_from_subsurface");
  setLinePrefix(my_key_+std::string(" evaluator"));

  pres_key_ = plist_.get<std::string>("pressure key", "pressure");
  dependencies_.insert(pres_key_);

  density_key_ = plist_.get<std::string>("mass density key", "mass_density_liquid");
  dependencies_.insert(density_key_);

  ponded_depth_key_ = plist_.get<std::string>("ponded depth key", "ponded_depth");
  dependencies_.insert(ponded_depth_key_);

  //  surface_density_key_ = plist_.get<std::string>("surface mass density key", "surface_mass_density_liquid");
  //  dependencies_.insert(surface_density_key_);

  surface_mesh_key_ = plist_.get<std::string>("surface mesh key", "surface");
  subsurface_mesh_key_ = plist_.get<std::string>("subsurface mesh key", "domain");
}

SurfaceCouplerViaSourceEvaluator::SurfaceCouplerViaSourceEvaluator(
        const SurfaceCouplerViaSourceEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    pres_key_(other.pres_key_),
    density_key_(other.density_key_),
    surface_density_key_(other.surface_density_key_),
    ponded_depth_key_(other.ponded_depth_key_),
    surface_mesh_key_(other.surface_mesh_key_),
    subsurface_mesh_key_(other.subsurface_mesh_key_),
    face_and_dirs_(other.face_and_dirs_) {}

Teuchos::RCP<FieldEvaluator> SurfaceCouplerViaSourceEvaluator::Clone() const {
  return Teuchos::rcp(new SurfaceCouplerViaSourceEvaluator(*this));
}


void SurfaceCouplerViaSourceEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // for now just passing... might do something later here?
  S->RequireField(my_key_, my_key_);
}


void SurfaceCouplerViaSourceEvaluator::IdentifyFaceAndDirection_(
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
    int direction = fdirs[index] > 0 ? 1 : -1;

    // Put (face,dir) into cached data.
    (*face_and_dirs_)[c] = std::make_pair(domain_face, direction);
  }
}

// Required methods from SecondaryVariableFieldEvaluator
void SurfaceCouplerViaSourceEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  if (face_and_dirs_ == Teuchos::null) {
    IdentifyFaceAndDirection_(S);
  }

  Teuchos::RCP<const AmanziMesh::Mesh> subsurface = S->GetMesh(subsurface_mesh_key_);

  //  Teuchos::RCP<const CompositeVector> surface_rho = S->GetFieldData(surface_density_key_);
  Teuchos::RCP<const CompositeVector> surface_depth = S->GetFieldData(ponded_depth_key_);

  Teuchos::RCP<const CompositeVector> rho = S->GetFieldData(density_key_);
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);
  const double& p_atm = *(S->GetScalarData("atmospheric_pressure"));
  const Epetra_Vector& gvec = *(S->GetConstantVectorData("gravity"));

  int ncells = result->size("cell",false);
  for (int c=0; c!=ncells; ++c) {
    AmanziMesh::Entity_ID_List cells;
    subsurface->face_get_cells((*face_and_dirs_)[c].first, AmanziMesh::OWNED, &cells);
    ASSERT(cells.size() == 1);

    // surface head
    // -- this will change to be the surface water density when that exists...
    double surface_head = (*surface_depth)("cell",c);

    // subsurface head
    // -- this will change to be the boundary water density when that exists...
    double subsurface_head = ((*pres)("face", (*face_and_dirs_)[c].first) - p_atm) /
        ((*rho)("cell",cells[0])*(-gvec[2]));

    // determine flux
    double eps = 1.e-6;
    double L = 0.01;
    double coef = 1.e-5;

    // upwind head in the "rel perm" term, k = head*coef
    if (surface_head > subsurface_head) {
      (*result)("cell",c) = (subsurface_head - surface_head) / L * std::pow<double>(surface_head, 5./3.) * coef;
    } else {
      (*result)("cell",c) = (subsurface_head - surface_head) / L * std::pow<double>(subsurface_head, 5./3.) * coef;
    }

    Teuchos::OSTab tab = getOSTab();
    if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
      *out_ << "Source, sub -> surface = " << (*result)("cell",c) << std::endl;
    }

  }
}

void SurfaceCouplerViaSourceEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  ASSERT(0);
  // this would require differentiating flux wrt pressure, which we
  // don't do for now.
}

} // namespace
} // namespace
} // namespace

/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $AMANZI_DIR/COPYRIGHT
Author Ethan Coon

Function applied to a mesh component, along with meta-data to store the values
of this function in a ComposteVector.

------------------------------------------------------------------------- */

#include "errors.hh"
#include "composite_vector_function.hh"

namespace Amanzi {
namespace Functions {

CompositeVectorFunction::CompositeVectorFunction(
        const Teuchos::RCP<const MeshFunction>& func,
        const std::vector<std::string>& names) :
    func_(func) {

  ASSERT(names.size() == func->size());

  // TODO: This is horrid.  Google won't tell me if STL has python's zip()
  // functionality, so this is probably bad C++.

  // Zip the two iterators, adding the resulting "tuple" to the container of pairs.
  std::vector<std::string>::const_iterator name = names.begin();
  for (MeshFunction::spec_iterator spec=func->begin();
       spec!=func->end(); ++spec) {
    cv_spec_list_.push_back(Teuchos::rcp(new CompositeVectorSpec(*name,*spec)));
    ++name;
  }
}

void CompositeVectorFunction::Compute(double time,
        const Teuchos::Ptr<CompositeVector>& cv) {
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = func_->mesh();
  cv->PutScalar(0.0);

  // create the input tuple
  int dim = mesh->space_dimension();
  double *args = new double[1+dim];
  double *xargs = args+1;
  args[0] = time;

  // loop over the name/spec pair
  for (CompositeVectorSpecList::const_iterator cv_spec=cv_spec_list_.begin();
       cv_spec!=cv_spec_list_.end(); ++cv_spec) {
    std::string compname = (*cv_spec)->first;
    Teuchos::RCP<MeshFunction::Spec> spec = (*cv_spec)->second;

    AmanziMesh::Entity_kind kind = spec->first->second;

    // loop over all regions in the spec
    for (MeshFunction::RegionList::const_iterator region=spec->first->first.begin();
         region!=spec->first->first.end(); ++region) {

      // special case for BOUNDARY_FACE
      if (kind == AmanziMesh::BOUNDARY_FACE) {
        if (mesh->valid_set_name(*region, AmanziMesh::FACE)) {
          // get the indices of the domain.
          AmanziMesh::Entity_ID_List id_list;
          mesh->get_set_entities(*region, AmanziMesh::FACE, AmanziMesh::USED, &id_list);

          const Epetra_Map& vandelay_map = mesh->exterior_face_epetra_map();

          // loop over indices
          AmanziMesh::Entity_ID_List cells;
          for (AmanziMesh::Entity_ID_List::const_iterator id=id_list.begin();
               id!=id_list.end(); ++id) {
            mesh->face_get_cells(*id, AmanziMesh::USED, &cells);
            if (cells.size() == 1) {
              AmanziMesh::Entity_ID bf = vandelay_map.LID(*id);

              // get the coordinate
              AmanziGeometry::Point xf = mesh->face_centroid(*id);
              for (int i=0; i!=dim; ++i) xargs[i] = xf[i];

              // evaluate the functions and stuff the result into the CV
              double *value = (*spec->second)(args);
              for (int i=0; i!=(*spec->second).size(); ++i) {
                (*cv)(compname, i, bf) = value[i];
              }
            }
          }
        } else {
          std::stringstream m;
          m << "unknown region: \"" << *region << "\"";
          Errors::Message message(m.str());
          Exceptions::amanzi_throw(message);
        }

      } else {
        if (mesh->valid_set_name(*region, kind)) {
          // get the indices of the domain.
          AmanziMesh::Entity_ID_List id_list;
          mesh->get_set_entities(*region, kind, AmanziMesh::USED, &id_list);

          // loop over indices
          for (AmanziMesh::Entity_ID_List::const_iterator id=id_list.begin();
               id!=id_list.end(); ++id) {
            // get the coordinate
            AmanziGeometry::Point xc;
            if (kind == AmanziMesh::CELL) {
              xc = mesh->cell_centroid(*id);
            } else if (kind == AmanziMesh::FACE) {
              xc = mesh->face_centroid(*id);
            } else if (kind == AmanziMesh::NODE) {
              mesh->node_get_coordinates(*id, &xc);
            } else {
              ASSERT(0);
            }
            for (int i=0; i!=dim; ++i) xargs[i] = xc[i];

            // evaluate the functions and stuff the result into the CV
            double *value = (*spec->second)(args);
            for (int i=0; i!=(*spec->second).size(); ++i) {
              (*cv)(compname, i, *id) = value[i];
            }
          }
        } else {
          std::stringstream m;
          m << "unknown region: \"" << *region << "\"";
          Errors::Message message(m.str());
          Exceptions::amanzi_throw(message);
        }
      }
    }
  }
  delete [] args;
}

} // namespace
} // namespace


/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!
  Function applied to a mesh component, along with meta-data to store the values
  of this function in a ComposteVector.
*/

#include "errors.hh"
#include "VerboseObject.hh"
#include "MeshDefs.hh"
#include "CompositeVectorFunction.hh"

namespace Amanzi {
namespace Functions {

CompositeVectorFunction::CompositeVectorFunction(const Teuchos::RCP<const MeshFunction>& func,
                                                 const std::vector<std::string>& names)
  : func_(func)
{
  AMANZI_ASSERT(names.size() == func->size());

  // Zip the two iterators, adding the resulting "tuple" to the container of pairs.
  std::vector<std::string>::const_iterator name = names.begin();
  for (MeshFunction::spec_iterator spec = func->begin(); spec != func->end(); ++spec) {
    cv_spec_list_.push_back(Teuchos::rcp(new CompositeVectorSpec(*name, *spec)));
    ++name;
  }
}

void
CompositeVectorFunction::Compute(double time,
                                 const Teuchos::Ptr<CompositeVector>& cv,
                                 const VerboseObject* vo)
{
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = func_->mesh();
  cv->PutScalar(0.);

#ifdef ENSURE_INITIALIZED_CVFUNCS
  // ensure all components are touched
  std::map<std::string, bool> done;
  for (auto compname : *cv) { done[compname] = false; }
#endif

  // create the input tuple
  int dim = mesh->getSpaceDimension();
  std::vector<double> args(1 + dim, 0.);
  args[0] = time;

  // loop over the name/spec pair
  for (CompositeVectorSpecList::const_iterator cv_spec = cv_spec_list_.begin();
       cv_spec != cv_spec_list_.end();
       ++cv_spec) {
    std::string compname = (*cv_spec)->first;
#ifdef ENSURE_INITIALIZED_CVFUNCS
    done[compname] = true;
#endif

    Epetra_MultiVector& compvec = *cv->ViewComponent(compname, false);
    Teuchos::RCP<MeshFunction::Spec> spec = (*cv_spec)->second;

    AmanziMesh::Entity_kind kind = spec->first->second;
    if (vo && vo->os_OK(Teuchos::VERB_EXTREME)) {
      *vo->os() << "Writing function on component: " << compname << " and entity "
                << AmanziMesh::to_string(kind) << std::endl;
    }

    // loop over all regions in the spec
    for (MeshFunction::RegionList::const_iterator region = spec->first->first.begin();
         region != spec->first->first.end();
         ++region) {
      // special case for BOUNDARY_FACE because BOUNDARY_FACE cannot currently
      // be used in sets, see #558.  Hopefully this will be fixed soon.
      if (kind == AmanziMesh::Entity_kind::BOUNDARY_FACE) {
        if (mesh->isValidSetName(*region, AmanziMesh::Entity_kind::FACE)) {
          // get the indices of the domain.
          auto id_list = mesh->getSetEntities(
            *region, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);

          const Epetra_Map& face_map = mesh->getMap(AmanziMesh::Entity_kind::FACE, false);
          const Epetra_Map& vandelay_map =
            mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false);

          // loop over indices
          AmanziMesh::Entity_ID_List cells;
          for (AmanziMesh::Entity_ID_List::const_iterator id = id_list.begin(); id != id_list.end();
               ++id) {
            cells = mesh->getFaceCells(*id, AmanziMesh::Parallel_type::ALL);
            if (cells.size() == 1) {
              AmanziMesh::Entity_ID bf = vandelay_map.LID(face_map.GID(*id));
              AMANZI_ASSERT(bf >= 0);

              // get the coordinate
              AmanziGeometry::Point xf = mesh->getFaceCentroid(*id);
              for (int i = 0; i != dim; ++i) args[i + 1] = xf[i];

              // evaluate the functions and stuff the result into the CV
              double* value = (*spec->second)(args);
              for (int i = 0; i != (*spec->second).size(); ++i) { compvec[i][bf] = value[i]; }
            }
          }
        } else {
          Errors::Message message;
          message << "CV: unknown region \"" << *region << "\"";
          Exceptions::amanzi_throw(message);
        }
      } else {
        bool valid = mesh->isValidSetName(*region, kind);
        if (vo && vo->os_OK(Teuchos::VERB_EXTREME)) {
          if (!valid) *vo->os() << "  region: " << *region << " not valid!" << std::endl;
        }
        if (valid) {
          // get the indices of the domain.
          auto id_list = mesh->getSetEntities(*region, kind, AmanziMesh::Parallel_type::OWNED);
          auto map = cv->Map().Map(compname, false);

          if (vo && vo->os_OK(Teuchos::VERB_EXTREME)) {
            *vo->os() << "  region: " << *region << " contains " << id_list.size()
                      << " local entities" << std::endl;
          }

          // loop over indices
          for (AmanziMesh::Entity_ID_List::const_iterator id = id_list.begin(); id != id_list.end();
               ++id) {
            // get the coordinate
            AmanziGeometry::Point xc;
            if (kind == AmanziMesh::Entity_kind::CELL) {
              xc = mesh->getCellCentroid(*id);
            } else if (kind == AmanziMesh::Entity_kind::FACE) {
              xc = mesh->getFaceCentroid(*id);
            } else if (kind == AmanziMesh::Entity_kind::NODE) {
              xc = mesh->getNodeCoordinate(*id);
            } else {
              AMANZI_ASSERT(0);
            }
            for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

            // evaluate the functions and stuff the result into the CV
            double* value = (*spec->second)(args);

            int g = map->FirstPointInElement(*id);
            int ndofs = map->ElementSize(*id);
            for (int n = 0; n < ndofs; ++n) {
              for (int i = 0; i != (*spec->second).size(); ++i) { compvec[i][g + n] = value[i]; }
            }
          }
        } else {
          Errors::Message message;
          message << "CV: unknown region \"" << *region << "\"";
          Exceptions::amanzi_throw(message);
        }
      }
    }
  }


#ifdef ENSURE_INITIALIZED_CVFUNCS
  for (auto compname : *cv) {
    if (!done[compname]) {
      Errors::Message message;
      message << "CV: component \"" << compname << "\" was not set";
      Exceptions::amanzi_throw(message);
    }
  }
#endif
}

} // namespace Functions
} // namespace Amanzi

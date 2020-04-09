/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//! <MISSING_ONELINE_DOCSTRING>

#include "errors.hh"
#include "MeshDefs.hh"
#include "CompositeVectorFunction.hh"

namespace Amanzi {
namespace Functions {

CompositeVectorFunction::CompositeVectorFunction(
  const Teuchos::RCP<const MeshFunction>& func,
  const std::vector<std::string>& names)
  : func_(func)
{
  AMANZI_ASSERT(names.size() == func->size());

  // Zip the two iterators, adding the resulting "tuple" to the container of
  // pairs.
  std::vector<std::string>::const_iterator name = names.begin();
  for (MeshFunction::spec_iterator spec = func->begin(); spec != func->end();
       ++spec) {
    cv_spec_list_.push_back(
      Teuchos::rcp(new CompositeVectorSpec(*name, *spec)));
    ++name;
  }
}

void
CompositeVectorFunction::Compute(double time, CompositeVector& cv)
{
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = func_->mesh();

  cv.putScalar(0.);

#ifdef ENSURE_INITIALIZED_CVFUNCS
  // ensure all components are touched
  std::map<std::string, bool> done;
  for (auto compname : cv) { done[compname] = false; }
#endif

  // create the input tuple
  int dim = mesh->space_dimension();
  Kokkos::View<double*> args("args", 1 + dim);
  for (int i = 0; i < args.extent(0); ++i) args(i) = 0;
  args[0] = time;

  // loop over the name/spec pair
  for (CompositeVectorSpecList::const_iterator cv_spec = cv_spec_list_.begin();
       cv_spec != cv_spec_list_.end();
       ++cv_spec) {
    std::string compname = (*cv_spec)->first;
#ifdef ENSURE_INITIALIZED_CVFUNCS
    done[compname] = true;
#endif

    auto compvec = cv.ViewComponent<AmanziDefaultHost>(compname, false);
    Teuchos::RCP<MeshFunction::Spec> spec = (*cv_spec)->second;

    AmanziMesh::Entity_kind kind = spec->first->second;

    // loop over all regions in the spec
    for (MeshFunction::RegionList::const_iterator region =
           spec->first->first.begin();
         region != spec->first->first.end();
         ++region) {
      // special case for BOUNDARY_FACE
      if (kind == AmanziMesh::BOUNDARY_FACE) {
        if (mesh->valid_set_name(*region, AmanziMesh::FACE)) {
          // get the indices of the domain.
          Kokkos::View<AmanziMesh::Entity_ID*> id_list;
          mesh->get_set_entities(*region,
                                 AmanziMesh::FACE,
                                 AmanziMesh::Parallel_type::OWNED,
                                 id_list);

          const auto& face_map = *mesh->face_map(false);
          const auto& vandelay_map = *mesh->exterior_face_map(false);

          // loop over indices
          Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cells;
          for (int id = 0; id < id_list.extent(0); ++id) {
            mesh->face_get_cells(
              id_list(id), AmanziMesh::Parallel_type::ALL, cells);
            if (cells.extent(0) == 1) {
              AmanziMesh::Entity_ID bf = vandelay_map.getLocalElement(
                face_map.getGlobalElement(id_list(id)));
              AMANZI_ASSERT(bf >= 0);

              // get the coordinate
              AmanziGeometry::Point xf = mesh->face_centroid(id_list(id));
              for (int i = 0; i != dim; ++i) args[i + 1] = xf[i];

              // evaluate the functions and stuff the result into the CV
              Kokkos::View<double*> value = (*spec->second)(args);
              for (int i = 0; i != (*spec->second).size(); ++i) {
                compvec(i, bf) = value[i];
              }
            }
          }
        } else {
          std::stringstream m;
          m << "CV: unknown boundary region: \"" << *region << "\"";
          Errors::Message message(m.str());
          Exceptions::amanzi_throw(message);
        }

      } else {
        if (mesh->valid_set_name(*region, kind)) {
          // get the indices of the domain.
          Kokkos::View<AmanziMesh::Entity_ID*> id_list;
          mesh->get_set_entities(
            *region, kind, AmanziMesh::Parallel_type::OWNED, id_list);

          // convert
          Kokkos::View<double**> txyz("txyz", dim + 1, id_list.extent(0));
          // Feed the array with data
          for (int id = 0; id < id_list.extent(0); ++id) {
            // get the coordinate
            AmanziGeometry::Point xc;
            if (kind == AmanziMesh::CELL) {
              xc = mesh->cell_centroid(id_list(id));
            } else if (kind == AmanziMesh::FACE) {
              xc = mesh->face_centroid(id_list(id));
            } else if (kind == AmanziMesh::NODE) {
              mesh->node_get_coordinates(id_list(id), &xc);
            } else {
              AMANZI_ASSERT(0);
            }
            txyz(0, id) = time;
            for (int i = 0; i < dim; ++i) txyz(i + 1, id) = xc[i];
          } // for

          Kokkos::View<double**, Kokkos::LayoutLeft> result(
            "result", id_list.extent(0), (*spec->second).size());
          spec->second->apply(txyz, result);

          assert(id_list.extent(0) == result.extent(0));
          assert((*spec->second).size() == result.extent(1));

          for (int id = 0; id < id_list.extent(0); ++id) {
            for (int i = 0; i < (*spec->second).size(); ++i) {
              compvec(id_list(id), i) = result(id, i);
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

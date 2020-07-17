/*
  Mesh Functions

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author Ethan Coon

  Function applied to a mesh component, along with meta-data to store the values
  of this function in a ComposteVector.
*/

#include "errors.hh"
#include "MeshDefs.hh"
#include "CompositeVectorFunction.hh"

namespace Amanzi {
namespace Functions {

CompositeVectorFunction::CompositeVectorFunction(
    const Teuchos::RCP<const MeshFunction>& func,
    const std::vector<std::string>& names) :
    func_(func) {

  AMANZI_ASSERT(names.size() == func->size());

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

  cv->PutScalar(0.);

#ifdef ENSURE_INITIALIZED_CVFUNCS  
  // ensure all components are touched
  std::map<std::string,bool> done;
  for (auto compname : *cv) {
    done[compname] = false;
  }
#endif

  // create the input tuple
  int dim = mesh->space_dimension();
  std::vector<double> args(1+dim, 0.);
  args[0] = time;

  // loop over the name/spec pair
  for (CompositeVectorSpecList::const_iterator cv_spec=cv_spec_list_.begin();
       cv_spec!=cv_spec_list_.end(); ++cv_spec) {
    std::string compname = (*cv_spec)->first;
#ifdef ENSURE_INITIALIZED_CVFUNCS  
    done[compname] = true;
#endif
    
    Epetra_MultiVector& compvec = *cv->ViewComponent(compname,false);
    Teuchos::RCP<MeshFunction::Spec> spec = (*cv_spec)->second;

    AmanziMesh::Entity_kind kind = spec->first->second;

    // loop over all regions in the spec
    for (MeshFunction::RegionList::const_iterator region=spec->first->first.begin();
         region!=spec->first->first.end(); ++region) {

      // region ENTIRE_MESH_REGION
      if (*region == std::string("ENTIRE_MESH_REGION")) {
        if (kind == AmanziMesh::BOUNDARY_FACE) {
          unsigned int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
          const Epetra_Map& vandelay_map = mesh->exterior_face_map(false);
          const Epetra_Map& face_map = mesh->face_map(false);

          // loop over indices
          AmanziMesh::Entity_ID_List cells;
          for (AmanziMesh::Entity_ID id=0; id!=nfaces; ++id) {
            mesh->face_get_cells(id, AmanziMesh::Parallel_type::ALL, &cells);
            if (cells.size() == 1) {
              AmanziMesh::Entity_ID bf = vandelay_map.LID(face_map.GID(id));
              AMANZI_ASSERT(bf >= 0);

              // get the coordinate
              AmanziGeometry::Point xf = mesh->face_centroid(id);
              for (int i=0; i!=dim; ++i) args[i+1] = xf[i];

              // evaluate the functions and stuff the result into the CV
              double *value = (*spec->second)(args);
              for (int i=0; i!=(*spec->second).size(); ++i) {
                compvec[i][bf] = value[i];
              }
            }
          }
        } else {
          unsigned int nentities = mesh->num_entities(kind, AmanziMesh::Parallel_type::OWNED);
          for (AmanziMesh::Entity_ID id=0; id!=nentities; ++id) {
            AmanziGeometry::Point xc;
            if (kind == AmanziMesh::CELL) {
              xc = mesh->cell_centroid(id);
            } else if (kind == AmanziMesh::FACE) {
              xc = mesh->face_centroid(id);
            } else if (kind == AmanziMesh::NODE) {
              mesh->node_get_coordinates(id, &xc);
            } else {
              Errors::Message msg;
              msg << "In CompositeVectorFunction: unknown mesh kind: \"" << kind << "\"";
              Exceptions::amanzi_throw(msg);
            }
            for (int i=0; i!=dim; ++i) args[i+1] = xc[i];

            // evaluate the functions and stuff the result into the CV
            double *value = (*spec->second)(args);
            for (int i=0; i!=(*spec->second).size(); ++i) {
              compvec[i][id] = value[i];
            }
          }
        }

      } else {
        // special case for BOUNDARY_FACE
        if (kind == AmanziMesh::BOUNDARY_FACE) {
          if (mesh->valid_set_name(*region, AmanziMesh::FACE)) {
            // get the indices of the domain.
            AmanziMesh::Entity_ID_List id_list;
            mesh->get_set_entities(*region, AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED, &id_list);

            const Epetra_Map& face_map = mesh->face_map(false);
            const Epetra_Map& vandelay_map = mesh->exterior_face_map(false);

            // loop over indices
            AmanziMesh::Entity_ID_List cells;
            for (AmanziMesh::Entity_ID_List::const_iterator id=id_list.begin();
                 id!=id_list.end(); ++id) {
              mesh->face_get_cells(*id, AmanziMesh::Parallel_type::ALL, &cells);
              if (cells.size() == 1) {
                AmanziMesh::Entity_ID bf = vandelay_map.LID(face_map.GID(*id));
                AMANZI_ASSERT(bf >= 0);


                // get the coordinate
                AmanziGeometry::Point xf = mesh->face_centroid(*id);
                for (int i=0; i!=dim; ++i) args[i+1] = xf[i];

                // evaluate the functions and stuff the result into the CV
                double *value = (*spec->second)(args);
                for (int i=0; i!=(*spec->second).size(); ++i) {
                  compvec[i][bf] = value[i];
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
            AmanziMesh::Entity_ID_List id_list;
            mesh->get_set_entities(*region, kind, AmanziMesh::Parallel_type::OWNED, &id_list);

	    auto map = cv->Map().Map(compname, false);

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
                AMANZI_ASSERT(0);
              }
              for (int i=0; i!=dim; ++i) args[i+1] = xc[i];

              // evaluate the functions and stuff the result into the CV
              double *value = (*spec->second)(args);

              int g = map->FirstPointInElement(*id);
              int ndofs = map->ElementSize(*id);
              for (int n = 0; n < ndofs; ++n) {
                for (int i=0; i!=(*spec->second).size(); ++i) {
                  compvec[i][g + n] = value[i];
                }
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

} // namespace
} // namespace

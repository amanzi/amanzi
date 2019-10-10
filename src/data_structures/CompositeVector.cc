/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon  
*/


//! <MISSING_ONELINE_DOCSTRING>

#include "Mesh.hh"
#include "CompositeVector.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Interpolate pressure ICs on cells to ICs for lambda (faces).
// -----------------------------------------------------------------------------
void DeriveFaceValuesFromCellValues(CompositeVector& cv) {
  if (cv.HasComponent("face")) {
    cv.ScatterMasterToGhosted("cell");

    const auto cv_c = cv.ViewComponent("cell",true);
    auto cv_f = cv.ViewComponent("face",false);
    const AmanziMesh::Mesh* mesh = &*cv.Mesh();

    Kokkos::parallel_for(cv_f.extent(0),
                         KOKKOS_LAMBDA(decltype(cv_f)::size_type f) {
      AmanziMesh::Entity_ID_View cells;
      mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
      int ncells = cells.size();

      double face_value = 0.0;
      for (int n=0; n!=ncells; ++n) {
        face_value += cv_c(cells[n], 0);
      }
      cv_f(f,0) = face_value / ncells;
    });
  }
  else if (cv.HasComponent("boundary_face")) {
    const auto cv_c = cv.ViewComponent("cell",true);
    auto cv_f = cv.ViewComponent("boundary_face",false);
    const AmanziMesh::Mesh* mesh = &*cv.Mesh();

    Map_ptr_type fb_map = cv.Mesh()->exterior_face_map(false);
    Map_ptr_type f_map = cv.Mesh()->face_map(false);

    Kokkos::parallel_for(cv_f.extent(0),
                         KOKKOS_LAMBDA(decltype(cv_f)::size_type fb) {
      AmanziMesh::Entity_ID_View cells;
      int f_gid = fb_map->getGlobalElement(fb);
      int f_lid = f_map->getLocalElement(f_gid);
      
      mesh->face_get_cells(f_lid, AmanziMesh::Parallel_type::ALL, cells);
      int ncells = cells.size();
      assert((ncells==1));
      cv_f(fb,0) = cv_c(cells[0],0);
    });
  }
}

} // namespace


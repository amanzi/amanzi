/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/
// Helpers for copying data between various data structures.

#include "MeshHelpers.hh"
#include "DataStructuresHelpers.hh"

namespace Amanzi {


// -----------------------------------------------------------------------------
// Interpolate pressure ICs on cells to ICs for lambda (faces).
// -----------------------------------------------------------------------------
void
DeriveFaceValuesFromCellValues(CompositeVector& cv)
{
  if (cv.hasComponent("face")) {
    cv.scatterMasterToGhosted("cell");

    const auto cv_c = cv.viewComponent("cell", true);
    auto cv_f = cv.viewComponent("face", false);
    const AmanziMesh::Mesh* mesh = &*cv.getMesh();

    Kokkos::parallel_for(
      "CompositeVector::DeriveFaceValuesFromCellValues loop 1",
      cv_f.extent(0),
      KOKKOS_LAMBDA(decltype(cv_f)::size_type f) {
        int ncells = mesh->getFaceNumCells(f, AmanziMesh::Parallel_kind::ALL);
        double face_value = 0.0;
        for (int n = 0; n != ncells; ++n) { face_value += cv_c(mesh->getFaceCell(f,n), 0); }
        cv_f(f, 0) = face_value / ncells;
      });
  } else if (cv.hasComponent("boundary_face")) {
    AmanziMesh::copyCellsToBoundaryFaces(*cv.getMesh(), *cv.getComponent("cell", false), *cv.getComponent("boundary_face", false));
  }
}


void
copyMeshCoordinatesToVector(const AmanziMesh::Mesh& mesh,
                            AmanziMesh::Entity_kind kind,
                            CompositeVector& vec)
{
  auto view = vec.viewComponent(to_string(kind), false);
  int ndim = mesh.getSpaceDimension();
  Kokkos::parallel_for("copyMeshCoordinatesToVector", view.extent(0),
                       KOKKOS_LAMBDA(const int& i) {
                         auto nc = mesh.getCentroid(kind, i);
                         for (int j = 0; j != ndim; ++j) view(i,j) = nc[j];
                       });
}


void
copyVectorToMeshCoordinates(const CompositeVector& vec, AmanziMesh::Mesh& mesh)
{
  const auto nodes = vec.viewComponent<DefaultHostMemorySpace>("node", true);
  int ndim = mesh.getSpaceDimension();

  AmanziMesh::Mesh::Entity_ID_View node_ids("node ids", nodes.extent(0));
  AmanziMesh::Mesh::Point_View new_positions("new pos", nodes.extent(0));
  Kokkos::parallel_for("copyVectorToMeshCoordinates", nodes.extent(0),
                       KOKKOS_LAMBDA(const int& n) {
                         node_ids[n] = n;
                         if (mesh.getSpaceDimension() == 2) {
                           new_positions[n] = Amanzi::AmanziGeometry::Point{ nodes(n,0), nodes(n,1) };
                         } else {
                           new_positions[n] = Amanzi::AmanziGeometry::Point{ nodes(n,0), nodes(n,1), nodes(n,2) };
                         }
                       });
  AmanziMesh::deform(mesh, node_ids, new_positions);
}


} // namespace Amanzi

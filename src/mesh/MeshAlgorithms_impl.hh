/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
#pragma once

#include "MeshAlgorithms_decl.hh"
#include "MeshCache_impl.hh"

namespace Amanzi {
namespace AmanziMesh {


// -----------------------------------------------------------------------------
// Given a boundary face ID, get the corresponding face ID
// -----------------------------------------------------------------------------
template <class Mesh_type>
KOKKOS_INLINE_FUNCTION Entity_ID
getBoundaryFaceFace(const Mesh_type& mesh, Entity_ID bf)
{
  return mesh.getBoundaryFaces()(bf);
}


// -----------------------------------------------------------------------------
// Given a boundary face ID, get the cell internal to that face.
// -----------------------------------------------------------------------------
template <class Mesh_type>
KOKKOS_INLINE_FUNCTION Entity_ID
getBoundaryFaceInternalCell(const Mesh_type& mesh, Entity_ID bf)
{
  return getFaceOnBoundaryInternalCell(mesh, getBoundaryFaceFace(mesh, bf));
}

// -----------------------------------------------------------------------------
// Given a face ID, and assuming it is a boundary face, get the cell internal.
// -----------------------------------------------------------------------------
template <class Mesh_type>
KOKKOS_INLINE_FUNCTION Entity_ID
getFaceOnBoundaryInternalCell(const Mesh_type& mesh, Entity_ID f)
{
  if constexpr (std::is_same<Mesh_type, Mesh>::value) {
    assert(mesh.getFaceNumCells(f, Parallel_kind::ALL) == 1 &&
           "getFaceOnBoundaryInternalCell() requires a face on the boundary.");
  } else {
    if (mesh.getFaceNumCells(f, Parallel_kind::ALL) != 1) {
      AmanziGeometry::Point fc = mesh.getFaceCentroid(f);
      std::stringstream msgs;
      msgs << "getFaceOnBoundaryInternalCell called with non-internal face GID "
           << mesh.getMap(AmanziMesh::Entity_kind::FACE, true)->getGlobalElement(f) << " at " << fc;
      Errors::Message msg(msgs.str());
      Exceptions::amanzi_throw(msg);
    }
  }
  return mesh.getFaceCell(f, 0);
}


// -----------------------------------------------------------------------------
// Exterior boundary normal: dir = 0 for internal face
// -----------------------------------------------------------------------------
template <class Mesh_type>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
getFaceNormalExterior(const Mesh_type& mesh, int f, int* dir)
{
  auto cells = mesh.getFaceCells(f);
  auto normal = mesh.getFaceNormal(f, cells[0], dir);
  if (cells.size() > 1) *dir = 0;
  return normal;
}


// -----------------------------------------------------------------------------
// Given a face ID, get the corresponding boundary face ID (assuming it is a bf)
// -----------------------------------------------------------------------------
template <MemSpace_kind MEM>
Entity_ID
getFaceOnBoundaryBoundaryFace(const MeshCache<MEM>& mesh, Entity_ID f)
{
  // should this be deprecated?  It seems likely to not perform well
  const auto& fmap = mesh.getMap(AmanziMesh::Entity_kind::FACE, true);
  const auto& bfmap = mesh.getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, true);
  return bfmap.LID(fmap.GID(f));
}

// -----------------------------------------------------------------------------
// Given a boundary face ID, get the cell internal to that face.
// -----------------------------------------------------------------------------
template <MemSpace_kind MEM>
Entity_ID
getBoundaryFaceInternalCell(const MeshCache<MEM>& mesh, Entity_ID bf)
{
  return getFaceOnBoundaryInternalCell(mesh, getBoundaryFaceFace(mesh, bf));
}


// -----------------------------------------------------------------------------
// Given a face ID, and assuming it is a boundary face, get the cell internal.
// -----------------------------------------------------------------------------
template <MemSpace_kind MEM>
Entity_ID
getFaceOnBoundaryInternalCell(const MeshCache<MEM>& mesh, Entity_ID f)
{
  auto cells = mesh.getFaceCells(f);
  assert(cells.size() == 1 && "getFaceOnBoundaryInternalCell() called with non-internal face GID");
  return cells[0];
}

// -----------------------------------------------------------------------------
// Given a vector on faces, import to vector on boundary faces
// -----------------------------------------------------------------------------
template <MemSpace_kind MEM>
void
copyFacesToBoundaryFaces(const MeshCache<MEM>& mesh,
                         const Epetra_MultiVector& faces,
                         Epetra_MultiVector& boundary_faces)
{
  int ierr = boundary_faces.Import(faces, mesh.getBoundaryFaceImporter(), Insert);
  AMANZI_ASSERT(!ierr);
}

// -----------------------------------------------------------------------------
// Given a vector on faces, import to vector on boundary faces
// -----------------------------------------------------------------------------
template <MemSpace_kind MEM>
void
copyBoundaryFacesToFaces(const MeshCache<MEM>& mesh,
                         const Epetra_MultiVector& boundary_faces,
                         Epetra_MultiVector& faces)
{
  int ierr = faces.Export(boundary_faces, mesh.getBoundaryFaceImporter(), Insert);
  AMANZI_ASSERT(!ierr);
}

// -----------------------------------------------------------------------------
// Given a vector on cells, set the boundary_face entries by their internal cell
// -----------------------------------------------------------------------------
template <MemSpace_kind MEM>
void
copyCellsToBoundaryFaces(const MeshCache<MEM>& mesh,
                         const Epetra_MultiVector& cells,
                         Epetra_MultiVector& boundary_faces)
{
  AMANZI_ASSERT(cells.NumVectors() == boundary_faces.NumVectors());
  for (Entity_ID bf = 0; bf != boundary_faces.MyLength(); ++bf) {
    Entity_ID c = getBoundaryFaceInternalCell(mesh, bf);
    for (int i = 0; i != boundary_faces.NumVectors(); ++i) { boundary_faces[i][bf] = cells[i][c]; }
  }
}


//
// This class provides default, virtual algorithms for computing geometric
// quantities given nodal coordinates and topological information.
//
// Split into two classes to aid in deletion of the MeshFramework class, while
// keeping the MeshAlgorithms class around for use by the Cache.
//
// lumped things for more efficient calculation
// KOKKOS_INLINE_FUNCTION
// std::pair<double, AmanziGeometry::Point>
// MeshAlgorithms::computeCellGeometry(const Mesh& mesh, const Entity_ID c) const {
//   return Impl::computeCellGeometry(mesh, c);
// }

std::pair<double, AmanziGeometry::Point>
MeshAlgorithms::computeCellGeometry(const Mesh& mesh, const Entity_ID c) const
{
  return Impl::computeCellGeometry(mesh, c);
}

// KOKKOS_INLINE_FUNCTION
// std::tuple<double, AmanziGeometry::Point, typename Mesh::cPoint_View>
// MeshAlgorithms::computeFaceGeometry(const Mesh& mesh, const Entity_ID f) const {
//   return Impl::computeFaceGeometry(mesh, f);
// }

std::tuple<double, AmanziGeometry::Point, typename Mesh::cPoint_View>
MeshAlgorithms::computeFaceGeometry(const Mesh& mesh, const Entity_ID f) const
{
  return Impl::computeFaceGeometry(mesh, f);
}

// KOKKOS_INLINE_FUNCTION
// std::pair<AmanziGeometry::Point, AmanziGeometry::Point>
// MeshAlgorithms::computeEdgeGeometry(const Mesh& mesh, const Entity_ID e) const {
//   return Impl::computeEdgeGeometry(mesh, e);
// }

std::pair<AmanziGeometry::Point, AmanziGeometry::Point>
MeshAlgorithms::computeEdgeGeometry(const Mesh& mesh, const Entity_ID e) const
{
  return Impl::computeEdgeGeometry(mesh, e);
}

// KOKKOS_INLINE_FUNCTION
// double
// MeshAlgorithms::computeCellVolume(const Mesh& mesh, const Entity_ID c) const {
//   return Impl::computeCellGeometry(mesh, c).first;
// }

double
MeshAlgorithms::computeCellVolume(const Mesh& mesh, const Entity_ID c) const
{
  return Impl::computeCellGeometry(mesh, c).first;
}

// KOKKOS_INLINE_FUNCTION
// AmanziGeometry::Point
// MeshAlgorithms::computeCellCentroid(const Mesh& mesh, const Entity_ID c) const {
//   return Impl::computeCellGeometry(mesh, c).second;
// }

AmanziGeometry::Point
MeshAlgorithms::computeCellCentroid(const Mesh& mesh, const Entity_ID c) const
{
  return Impl::computeCellGeometry(mesh, c).second;
}

// KOKKOS_INLINE_FUNCTION
// double
// MeshAlgorithms::computeFaceArea(const Mesh& mesh, const Entity_ID f) const {
//   return std::get<0>(Impl::computeFaceGeometry(mesh, f));
// }

double
MeshAlgorithms::computeFaceArea(const Mesh& mesh, const Entity_ID f) const
{
  return std::get<0>(Impl::computeFaceGeometry(mesh, f));
}

// KOKKOS_INLINE_FUNCTION
// AmanziGeometry::Point
// MeshAlgorithms::computeFaceCentroid(const Mesh& mesh, const Entity_ID f) const {
//   return std::get<1>(Impl::computeFaceGeometry(mesh, f));
// }

AmanziGeometry::Point
MeshAlgorithms::computeFaceCentroid(const Mesh& mesh, const Entity_ID f) const
{
  return std::get<1>(Impl::computeFaceGeometry(mesh, f));
}

// KOKKOS_INLINE_FUNCTION
// AmanziGeometry::Point
// MeshAlgorithms::computeFaceNormal(
//   const Mesh& mesh, const Entity_ID f, const Entity_ID c, int * const orientation) const {
//   assert(false && "computeFaceNormal() not implemented on device");
// }

AmanziGeometry::Point
MeshAlgorithms::computeFaceNormal(const Mesh& mesh,
                                  const Entity_ID f,
                                  const Entity_ID c,
                                  int* const orientation) const
{
  auto geom = Impl::computeFaceGeometry(mesh, f);

  Mesh::cEntity_ID_View fcells;
  mesh.getFaceCells(f, fcells);
  if (orientation) *orientation = 0;

  Entity_ID cc;
  std::size_t i;
  if (c < 0) {
    cc = fcells[0];
    i = 0;
  } else {
    cc = c;
    auto ncells = fcells.size();
    for (i = 0; i != ncells; ++i)
      if (fcells[i] == cc) break;
  }
  AmanziGeometry::Point normal = std::get<2>(geom)[i];

  if (mesh.getSpaceDimension() == mesh.getManifoldDimension()) {
    if (c < 0) {
      assert(orientation == nullptr);
      normal *= Impl::getFaceDirectionInCell(mesh, f, cc);
    } else if (orientation) {
      *orientation = Impl::getFaceDirectionInCell(mesh, f, cc);
    }
  } else {
    // manifold case
    if (c < 0) {
      assert(orientation == nullptr);

      if (fcells.size() != 2) {
        normal *= Impl::getFaceDirectionInCell(mesh, f, cc);
      } else {
        // average normals oriented from lower to higher GIDs
        int pos_i = mesh.getEntityGID(Entity_kind::CELL, fcells[0]) >
                        mesh.getEntityGID(Entity_kind::CELL, fcells[1]) ?
                      0 :
                      1;
        normal = (std::get<2>(geom)[1 - pos_i] - std::get<2>(geom)[pos_i]) / 2;
      }
    } else if (orientation) {
      *orientation = Impl::getFaceDirectionInCell(mesh, f, cc);
    }
  }

  if (orientation) assert(*orientation != 0);
  return normal;
}

// KOKKOS_INLINE_FUNCTION
// double
// MeshAlgorithms::computeEdgeLength(const Mesh& mesh, const Entity_ID e) const {
//   return AmanziGeometry::norm(Impl::computeEdgeGeometry(mesh, e).first);
// }

double
MeshAlgorithms::computeEdgeLength(const Mesh& mesh, const Entity_ID e) const
{
  return AmanziGeometry::norm(Impl::computeEdgeGeometry(mesh, e).first);
}

// KOKKOS_INLINE_FUNCTION
// AmanziGeometry::Point
// MeshAlgorithms::computeEdgeVector(const Mesh& mesh, const Entity_ID e, const Entity_ID n, int * const orientation) const {
//   auto geom = Impl::computeEdgeGeometry(mesh, e);
//   if (n >= 0) {
//     auto nodes = mesh.getEdgeNodes(e);
//     if (n == nodes[0]) {
//       if (orientation) *orientation = 1;
//       return geom.first;
//     } else if (n == nodes[1]) {
//       if (orientation) *orientation = -1;
//       return -geom.first;
//     } else {
//       AMANZI_ASSERT(0);
//     }
//   }
//   return geom.first;
// }

AmanziGeometry::Point
MeshAlgorithms::computeEdgeVector(const Mesh& mesh,
                                  const Entity_ID e,
                                  const Entity_ID n,
                                  int* const orientation) const
{
  auto geom = Impl::computeEdgeGeometry(mesh, e);
  if (n >= 0) {
    auto nodes = mesh.getEdgeNodes(e);
    if (n == nodes[0]) {
      if (orientation) *orientation = 1;
      return geom.first;
    } else if (n == nodes[1]) {
      if (orientation) *orientation = -1;
      return -geom.first;
    } else {
      AMANZI_ASSERT(0);
    }
  }
  return geom.first;
}

// KOKKOS_INLINE_FUNCTION
// AmanziGeometry::Point
// MeshAlgorithms::computeEdgeCentroid(const Mesh& mesh, const Entity_ID e) const {
//   return Impl::computeEdgeGeometry(mesh, e).second;
// }

AmanziGeometry::Point
MeshAlgorithms::computeEdgeCentroid(const Mesh& mesh, const Entity_ID e) const
{
  return Impl::computeEdgeGeometry(mesh, e).second;
}

// KOKKOS_INLINE_FUNCTION
// void
// MeshAlgorithms::computeCellFacesAndBisectors(const Mesh& mesh, const Entity_ID cellid,
//         typename Mesh::cEntity_ID_View& faceids, typename Mesh::cPoint_View * const bisectors) const {
//   mesh.getCellFaces(cellid, faceids);
//   if (bisectors)
//     *bisectors = Impl::computeBisectors(mesh, cellid, faceids);
// }

void
MeshAlgorithms::computeCellFacesAndBisectors(const Mesh& mesh,
                                             const Entity_ID cellid,
                                             typename Mesh::cEntity_ID_View& faceids,
                                             typename Mesh::cPoint_View* const bisectors) const
{
  mesh.getCellFaces(cellid, faceids);
  if (bisectors) *bisectors = Impl::computeBisectors(mesh, cellid, faceids);
}


// //
// // Creates a pair of maps, <ALL, OWNED>, for a given entity_kind.
// // Uses GIDs provided by the Mesh object.
// //
// template <class Mesh_type>
// std::pair<Map_ptr_type, Map_ptr_type>
// createMapsFromMeshGIDs(const Mesh_type& mesh, const Entity_kind kind)
// {
//   Entity_ID num_owned = mesh.getNumEntities(kind, Parallel_kind::OWNED);
//   auto gids = mesh.getEntityGIDs(kind, true);
//   for (auto gid : gids) AMANZI_ASSERT(gid >= 0);
//   auto gids_owned =
//     Kokkos::subview(gids, Kokkos::make_pair((std::size_t)0, (std::size_t)num_owned));
//   return std::make_pair(Teuchos::rcp(new Map_type(-1, gids, 0, mesh.getComm())),
//                         Teuchos::rcp(new Map_type(-1, gids_owned, 0, mesh.getComm())));
// }

// //
// // Given discontinuous maps, make continuous maps
// //
// inline std::pair<Map_ptr_type, Map_ptr_type>
// createContiguousMaps(const Map_ptr_type& ghosted, const Map_ptr_type& owned)
// {
//   // create the owned contiguous map.
//   auto owned_contiguous = Teuchos::rcp(
//     new Map_type(owned->getGlobalNumElements(), owned->getLocalNumElements(), 0, owned->getComm()));

//   // communicated owned to ghosted using the mesh's maps
//   Vector_type_<GO> owned_contiguous_vec(owned);
//   {
//     auto data = owned_contiguous_vec.getLocalViewHost(Tpetra::Access::ReadWrite);
//     for (int i = 0; i != data.size(); ++i) { data(i, 0) = owned_contiguous->getGlobalElement(i); }
//   }
//   Import_type importer(owned, ghosted);
//   Vector_type_<GO> all_contiguous_vec(ghosted);
//   all_contiguous_vec.doImport(owned_contiguous_vec, importer, Tpetra::CombineMode::INSERT);

//   auto ghosted_contiguous =
//     Teuchos::rcp(new Map_type(-1, all_contiguous_vec.getData()(), 0, ghosted->getComm()));
//   return std::make_pair(ghosted_contiguous, owned_contiguous);
// }


// //
// // Creates a pair of maps, <ALL, OWNED>, for a given entity_kind.
// // Uses a natural ordering of GIDs, proc 0 == 0...n, proc 1 = n..., etc.
// //
// template <class Mesh_type>
// std::pair<Map_ptr_type, Map_ptr_type>
// createMapsFromContiguousGIDs(const Mesh_type& mesh, const Entity_kind kind)
// {
//   auto [ghosted_mesh_map, owned_mesh_map] = createMapsFromMeshGIDs(mesh, kind);
//   return createContiguousMaps(ghosted_mesh_map, owned_mesh_map);
// }


template <MemSpace_kind MEM>
void
cacheAll(MeshCache<MEM>& mesh)
{
  std::cout << "############# CacheAll" << std::endl;
  mesh.cacheCellNodes();
  mesh.cacheCellCoordinates();
  mesh.cacheFaceCoordinates();
  mesh.cacheNodeCells();
  mesh.cacheNodeFaces();
  if (mesh.hasEdges()) {
    mesh.cacheCellEdges();
    mesh.cacheEdgeCells();
    mesh.cacheNodeEdges();
    mesh.cacheEdgeNodes();
    mesh.cacheEdgeCoordinates();
  }
}

template <MemSpace_kind MEM>
void
cacheDefault(MeshCache<MEM>& mesh)
{
  std::cout << "############# CacheDefault ";
  std::cout << " hasNodes: " << mesh.hasNodes() << " hasEdges: " << mesh.hasEdges() << std::endl;
  // caches what the developers currently think is best
  if (mesh.hasNodes()) { mesh.cacheNodeCoordinates(); }
  mesh.cacheCellFaces();
  mesh.cacheFaceCells();
  if (mesh.hasNodes()) { mesh.cacheFaceNodes(); }
  mesh.cacheCellGeometry();
  mesh.cacheFaceGeometry();
  if (mesh.hasEdges()) {
    mesh.cacheFaceEdges();
    mesh.cacheEdgeFaces();
    mesh.cacheEdgeGeometry();
  }
}


template <class Mesh_type>
void
deform(Mesh_type& mesh,
       const typename Mesh_type::cEntity_ID_View& nodeids,
       const typename Mesh_type::cPoint_View& newpos)
{
  mesh.setNodeCoordinates(nodeids, newpos);
}


} // namespace AmanziMesh
} // namespace Amanzi

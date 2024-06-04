/*
   2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
*/

#pragma once

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>

#include "AmanziMap.hh"
#include "AmanziVector.hh"
#include "ViewUtils.hh"
#include "Geometry.hh"
#include "CellTopology.hh"

#include "MeshAudit_decl.hh"

namespace Amanzi {
namespace AmanziMesh {
namespace Impl {

#define MAX_OUT 5

template <class Mesh_type>
MeshAudit_Base<Mesh_type>::MeshAudit_Base(const Teuchos::RCP<const Mesh_type>& mesh,
                                          std::ostream& os)
  : mesh_(mesh),
    comm_(mesh_->getComm()),
    os_(os),
    nnodes_all_(mesh_->getNumEntities(Entity_kind::NODE, Parallel_kind::ALL)),
    nfaces_all_(mesh_->getNumEntities(Entity_kind::FACE, Parallel_kind::ALL)),
    ncells_all_(mesh_->getNumEntities(Entity_kind::CELL, Parallel_kind::ALL)),
    nnodes_owned_(mesh_->getNumEntities(Entity_kind::NODE, Parallel_kind::OWNED)),
    nfaces_owned_(mesh_->getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED)),
    ncells_owned_(mesh_->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED))
{}


template <class MeshAudit_type>
void
createTestDependencies_Base(MeshAudit_type& audit, AuditDirectedGraph<MeshAudit_type>& graph)
{}


// This creates the dependency graph for the individual tests.  Adding a new
// test to the graph is a simple matter of adding a new stanza of the form
//
//   Graph::vertex_descriptor my_test_handle = AddVertex_(g);
//   g_[my_test_handle].name = "description of my test";
//   g_[my_test_handle].test = &MeshAudit::my_test;
//   addEdge_(other_test_handle, my_test_handle);
//
// The last line specifies that other_test_handle is a pre-requisite for
// my_test_handle.  There may be multiple pre-requisites or none.
template <class MeshAudit_type>
void
createTestDependencies_Geometry(MeshAudit_type& audit, AuditDirectedGraph<MeshAudit_type>& graph)
{
  createTestDependencies_Base(audit, graph);

  // Cell_to_nodes tests
  auto test02 = graph.AddVertex("cell_to_nodes", &MeshAudit_type::checkCellToNodes);

  auto test03 = graph.AddVertex("node references by cells", &MeshAudit_type::checkNodeRefsByCells);
  graph.AddEdge(test02, test03);

  // Cell_to_faces tests
  auto test04 = graph.AddVertex("cell_to_faces", &MeshAudit_type::checkCellToFaces);

  auto test05 = graph.AddVertex("face references by cells", &MeshAudit_type::checkFaceRefsByCells);
  graph.AddEdge(test04, test05);

  // face_to_nodes tests
  auto test06 = graph.AddVertex("face_to_nodes", &MeshAudit_type::checkFaceToNodes);

  auto test08 = graph.AddVertex("node references by faces", &MeshAudit_type::checkNodeRefsByFaces);
  graph.AddEdge(test06, test08);

  // cell_to_face_dirs tests
  auto test09 = graph.AddVertex("cell_to_face_dirs", &MeshAudit_type::checkCellToFaceDirs);

  // cell degeneracy test
  auto test10 =
    graph.AddVertex("topological non-degeneracy of cells", &MeshAudit_type::checkCellDegeneracy);
  graph.AddEdge(test02, test10);

  // Consistency between the various mesh connectivity data.
  auto test11 = graph.AddVertex("consistency of mesh connectivity data",
                                &MeshAudit_type::checkCellToFacesToNodes);
  graph.AddEdge(test02, test11);
  graph.AddEdge(test04, test11);
  graph.AddEdge(test06, test11);
  graph.AddEdge(test09, test11);
  graph.AddEdge(test10, test11);

  // node_to_coordinates tests
  auto test12 = graph.AddVertex("node_to_coordinates", &MeshAudit_type::checkNodeToCoordinates);

  // cell_to_coordinates tests
  auto test13 = graph.AddVertex("cell_to_coordinates", &MeshAudit_type::checkCellToCoordinates);
  graph.AddEdge(test02, test13);
  graph.AddEdge(test12, test13);

  // face_to_coordinates tests
  auto test14 = graph.AddVertex("face_to_coordinates", &MeshAudit_type::checkFaceToCoordinates);
  graph.AddEdge(test06, test14);
  graph.AddEdge(test12, test14);

  // cell topology/geometry test
  auto test15 = graph.AddVertex("cell geometry", &MeshAudit_type::checkCellGeometry);
  graph.AddEdge(test10, test15);
  graph.AddEdge(test13, test15);

  // cell-to-face adjacencies and orientations
  auto test16 = graph.AddVertex("consistency of face-cell adjacencies",
                                &MeshAudit_type::checkFaceCellAdjacencyConsistency);
  graph.AddEdge(test04, test16);
  graph.AddEdge(test09, test16);

  // face normals relto cell
  auto test17 =
    graph.AddVertex("face normals are outward", &MeshAudit_type::checkFaceNormalReltoCell);
  graph.AddEdge(test16, test17);

  // face normals
  auto test18 =
    graph.AddVertex("face normals are oriented", &MeshAudit_type::checkFaceNormalOrientation);
  graph.AddEdge(test16, test18);
}


template <class MeshAudit_type>
void
createTestDependencies_Maps(MeshAudit_type& audit, AuditDirectedGraph<MeshAudit_type>& graph)
{
  createTestDependencies_Geometry(audit, graph);

  // Entity_counts tests
  auto test01 = graph.AddVertex("entity_counts", &MeshAudit_type::checkEntityCounts);

  // map tests
  auto test16 = graph.AddVertex("owned and overlap node maps", &MeshAudit_type::checkNodeMaps);

  auto test17 = graph.AddVertex("owned and overlap face maps", &MeshAudit_type::checkFaceMaps);

  auto test18 = graph.AddVertex("owned and overlap cell maps", &MeshAudit_type::checkCellMaps);

  // ghost data tests
  auto test19 = graph.AddVertex("node_to_coordinates ghost data",
                                &MeshAudit_type::checkNodeToCoordinatesGhostData);
  graph.AddEdge("node_to_coordinates", test19);
  graph.AddEdge("owned and overlap node maps", test19);

  auto test20 =
    graph.AddVertex("face_to_nodes ghost data", &MeshAudit_type::checkFaceToNodesGhostData);
  graph.AddEdge("face_to_nodes", test20);
  graph.AddEdge("owned and overlap node maps", test20);
  graph.AddEdge("owned and overlap face maps", test20);

  auto test21 =
    graph.AddVertex("cell_to_nodes ghost data", &MeshAudit_type::checkCellToNodesGhostData);
  graph.AddEdge("cell_to_nodes", test21);
  graph.AddEdge("owned and overlap cell maps", test21);
  graph.AddEdge("owned and overlap node maps", test21);

  auto test22 =
    graph.AddVertex("cell_to_faces ghost data", &MeshAudit_type::checkCellToFacesGhostData);
  graph.AddEdge("cell_to_faces", test22);
  graph.AddEdge("owned and overlap cell maps", test22);
  graph.AddEdge("owned and overlap face maps", test22);

  // partition tests
  auto test32 = graph.AddVertex("face partition", &MeshAudit_type::checkFacePartition);
  graph.AddEdge("owned and overlap face maps", test32);
  graph.AddEdge("owned and overlap cell maps", test32);
  graph.AddEdge(test18, test32);

  auto test33 = graph.AddVertex("node partition", &MeshAudit_type::checkNodePartition);
  graph.AddEdge("owned and overlap node maps", test33);
  graph.AddEdge("owned and overlap cell maps", test33);
};

template <class MeshAudit_type>
void
createTestDependencies_Sets(MeshAudit_type& audit, AuditDirectedGraph<MeshAudit_type>& graph)
{
  createTestDependencies_Maps(audit, graph);

  // node set data tests
  auto test23 = graph.AddVertex("node set IDs", &MeshAudit_type::checkNodeSetIDs);

  auto test24 = graph.AddVertex("node sets", &MeshAudit_type::checkNodeSets);
  graph.AddEdge("owned and overlap node maps", test24);

  auto test25 = graph.AddVertex("valid node set IDs", &MeshAudit_type::checkValidNodeSetID);
  graph.AddEdge("node set IDs", test25);

  // face set data tests
  auto test26 = graph.AddVertex("face set IDs", &MeshAudit_type::checkFaceSetIDs);

  auto test27 = graph.AddVertex("face sets", &MeshAudit_type::checkFaceSets);
  graph.AddEdge("owned and overlap face maps", test27);
  graph.AddEdge("face set IDs", test27);

  auto test28 = graph.AddVertex("valid face set IDs", &MeshAudit_type::checkValidFaceSetID);
  graph.AddEdge("face set IDs", test28);

  // cell set data tests
  auto test29 = graph.AddVertex("cell set IDs", &MeshAudit_type::checkCellSetIDs);

  auto test30 = graph.AddVertex("cell sets", &MeshAudit_type::checkCellSets);
  graph.AddEdge("owned and overlap cell maps", test30);
  graph.AddEdge("cell set IDs", test30);

  auto test31 = graph.AddVertex("valid cell set IDs", &MeshAudit_type::checkValidCellSetID);
  graph.AddEdge("cell set IDs", test31);
}

////////////////////////////////////////////////////////////////////////////////
//
// Tests (and auxillary functions) follow.  Tests must be const functions
// that take no arguments and that return a bool result: true if an error
// was found, otherwise false.  Tests must be considered collective procedures
// when in a parallel context, and so return a comm_on collective result.
// This is easily done using the function globalAny_(bool) which returns true
// on all processes if the argument is true on any of the processes.
//
////////////////////////////////////////////////////////////////////////////////

// The count_entities method should return values that match the number of
// elements in the corresponding Map_types.  This applies to nodes, faces
// and cells, including ghosts and not.  Here the maps are considered to be
// the authoritative source of information.  A positive value is returned
// if any discrepancy is found, but it is safe to perform other tests as
// they do not use these methods.
template <class Mesh_type>
bool
MeshAudit_Geometry<Mesh_type>::checkEntityCounts() const
{
  int n, nref;
  bool error = false;
  // Check the number of owned nodes.
  n = mesh_->getNumEntities(Entity_kind::NODE, Parallel_kind::OWNED);
  nref = mesh_->getMap(AmanziMesh::Entity_kind::NODE, false)->getLocalNumElements();
  if (n != nref) {
    os_ << ": ERROR: getNumEntities(NODE,OWNED)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  // Check the number of used nodes.
  n = mesh_->getNumEntities(Entity_kind::NODE, Parallel_kind::ALL);
  nref = mesh_->getMap(AmanziMesh::Entity_kind::NODE, true)->getLocalNumElements();
  if (n != nref) {
    os_ << "ERROR: getNumEntities(NODE,ALL)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  // Check the number of owned faces.
  n = mesh_->getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
  nref = mesh_->getMap(AmanziMesh::Entity_kind::FACE, false)->getLocalNumElements();
  if (n != nref) {
    os_ << "ERROR: getNumEntities(FACE,OWNED)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  // Check the number of used faces.
  n = mesh_->getNumEntities(Entity_kind::FACE, Parallel_kind::ALL);
  nref = mesh_->getMap(AmanziMesh::Entity_kind::FACE, true)->getLocalNumElements();
  if (n != nref) {
    os_ << "ERROR: getNumEntities(FACE,ALL)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  // Check the number of owned cells.
  n = mesh_->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  nref = mesh_->getMap(AmanziMesh::Entity_kind::CELL, false)->getLocalNumElements();
  if (n != nref) {
    os_ << "ERROR: getNumEntities(CELL,OWNED)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  // Check the number of used cells.
  n = mesh_->getNumEntities(Entity_kind::CELL, Parallel_kind::ALL);
  nref = mesh_->getMap(AmanziMesh::Entity_kind::CELL, true)->getLocalNumElements();
  if (n != nref) {
    os_ << "ERROR: getNumEntities(CELL,ALL)=" << n << "; should be " << nref << std::endl;
    error = true;
  }
  return globalAny_(error);
}


// Check that cell_to_nodes successfully returns valid references to local
// nodes.  A nonzero return value signals an error, and further tests using its
// data should be avoided.
DISABLE_CUDA_WARNING
template <class Mesh_type>
bool
MeshAudit_Geometry<Mesh_type>::checkCellToNodes() const
{
  Entity_ID_View bad_cells("bad cells", ncells_all_);

  const Mesh_type& m = *mesh_;
  int nnodes_all = nnodes_all_;
  using Range_type = Kokkos::RangePolicy<LO, typename cEntity_ID_View::execution_space>;
  Kokkos::parallel_for(
    "MeshAudit::checkCellToNodes", Range_type(0, ncells_all_), KOKKOS_LAMBDA(const LO c) {
      auto cnode = m.getCellNodes(c);
      for (auto& n : cnode) {
        if (n >= nnodes_all) bad_cells(c) = 1;
      }
    });

  bool error = checkErrorList(bad_cells, "invalid nodes referenced by cells", this->os_);
  return globalAny_(error);
}


// Check that every node is referenced by at least one cell.  This assumes
// that cell_to_nodes have been verified to return valid data.  A nonzero
// return value indicates that one or more nodes are not attached to any cell.
DISABLE_CUDA_WARNING
template <class Mesh_type>
bool
MeshAudit_Geometry<Mesh_type>::checkNodeRefsByCells() const
{
  Entity_ID_View free_nodes("MeshAudit::checkNodeRefsByCells", nnodes_all_);
  Kokkos::deep_copy(free_nodes, 1);

  const Mesh_type& m = *mesh_;
  using Range_type = Kokkos::RangePolicy<LO, typename cEntity_ID_View::execution_space>;
  Kokkos::parallel_for(
    "MeshAudit::checkNodeRefsByCells",
    Range_type(0, ncells_all_),
    KOKKOS_LAMBDA(const Entity_ID c) {
      cEntity_ID_View cnode;
      m.getCellNodes(c, cnode);
      for (auto& n : cnode) Kokkos::atomic_store(&free_nodes(n), 0);
    });

  bool error = checkErrorList(free_nodes, "found unreferenced nodes", this->os_);
  return globalAny_(error);
}


// Check that cell_to_faces successfully returns valid references to local
// faces.  A nonzero return value
// signals an error, and further tests using its data should be avoided.
template <class Mesh_type>
bool
MeshAudit_Geometry<Mesh_type>::checkCellToFaces() const
{
  Entity_ID_View bad_cells("bad cells", ncells_all_);

  using Range_type = Kokkos::RangePolicy<LO, typename cEntity_ID_View::execution_space>;
  const Mesh_type& m = *mesh_;
  int nfaces_all = nfaces_all_;
  Kokkos::parallel_for(
    "MeshAudit::checkCellToFaces", Range_type(0, ncells_all_), KOKKOS_LAMBDA(const Entity_ID c) {
      auto cface = m.getCellFaces(c);
      for (auto& f : cface) {
        if (f >= nfaces_all) bad_cells(c) = 1;
      }
    });

  bool error = checkErrorList(bad_cells, "invalid faces referenced by cells", this->os_);
  return globalAny_(error);
}


// Check that every face is referenced by exactly one or two cells.  This
// assumes that cell_to_faces have been verified to return valid data.  A
// nonzero return value indicates that faces were found that do not belong
// to any cell or belong to more than two cells (bad topology).
template <class Mesh_type>
bool
MeshAudit_Geometry<Mesh_type>::checkFaceRefsByCells() const
{
  Entity_ID_View face_counts("MeshAudit::checkFaceRefsByCells", nfaces_all_);
  Entity_ID_View free_faces("MeshAudit::checkFaceRefsByCells", nfaces_all_);

  const Mesh_type& m = *mesh_;
  using Range_type = Kokkos::RangePolicy<LO, typename cEntity_ID_View::execution_space>;
  Kokkos::parallel_for(
    "MeshAudit::checkNodeRefsByCells",
    Range_type(0, ncells_all_),
    KOKKOS_LAMBDA(const Entity_ID c) {
      cEntity_ID_View cface;
      m.getCellFaces(c, cface);
      for (auto& n : cface) Kokkos::atomic_add(&face_counts(n), 1);
    });


  Kokkos::parallel_for(
    "MeshAudit::checkNodeRefsByCells2",
    Range_type(0, nfaces_all_),
    KOKKOS_LAMBDA(const Entity_ID f) {
      if (face_counts(f) == 0) free_faces(f) = 1;
      if (face_counts(f) <= 2) face_counts(f) = 0;
    });

  bool error = checkErrorList(free_faces, "found unreferenced faces", this->os_);

  // note that 3D space meshes of 2D manifolds may allow for more than 2 cells on a face
  if (mesh_->getSpaceDimension() == mesh_->getManifoldDimension()) {
    error |= checkErrorList(face_counts, "found faces shared by more than two cells", this->os_);
  }
  return globalAny_(error);
}


// Check that face_to_nodes successfully returns valid references to local
// nodes.  A nonzero return value signals an error, and further tests using its
// data should be avoided.
template <class Mesh_type>
bool
MeshAudit_Geometry<Mesh_type>::checkFaceToNodes() const
{
  Entity_ID_View bad_faces("bad faces", nfaces_all_);

  using Range_type = Kokkos::RangePolicy<LO, typename cEntity_ID_View::execution_space>;
  const Mesh_type& m = *mesh_;
  int nnodes_all = nnodes_all_;
  Kokkos::parallel_for(
    "MeshAudit::checkFaceToNodes", Range_type(0, nfaces_all_), KOKKOS_LAMBDA(const LO c) {
      auto cnode = m.getFaceNodes(c);
      for (auto& n : cnode) {
        if (n >= nnodes_all) bad_faces(c) = 1;
      }
    });

  bool error = checkErrorList(bad_faces, "invalid nodes referenced by faces", this->os_);
  return globalAny_(error);
}


// Check that every node is referenced by at least one face.  This assumes
// that face_to_nodes has been verified to return valid data.  A nonzero
// return value indicates that one or more nodes are not attached to any face.
template <class Mesh_type>
bool
MeshAudit_Geometry<Mesh_type>::checkNodeRefsByFaces() const
{
  Entity_ID_View free_nodes("MeshAudit::checkNodeRefsByFaces", nnodes_all_);
  Kokkos::deep_copy(free_nodes, 1);

  const Mesh_type& m = *mesh_;
  using Range_type = Kokkos::RangePolicy<LO, typename cEntity_ID_View::execution_space>;
  Kokkos::parallel_for(
    "MeshAudit::checkNodeRefsByFaces",
    Range_type(0, nfaces_all_),
    KOKKOS_LAMBDA(const Entity_ID c) {
      cEntity_ID_View cnode;
      m.getFaceNodes(c, cnode);
      for (auto& n : cnode) Kokkos::atomic_store(&free_nodes(n), 0);
    });

  bool error = checkErrorList(free_nodes, "found unreferenced nodes", this->os_);
  return globalAny_(error);
}

// Check that cell_to_face_dirs successfully returns data for all cells and
// that the values are either +1 or -1. If this test fails, further tests using
// this data should be avoided.
template <class Mesh_type>
bool
MeshAudit_Geometry<Mesh_type>::checkCellToFaceDirs() const
{
  Entity_ID_View bad_cells("bad cells", ncells_all_);

  using Range_type = Kokkos::RangePolicy<LO, typename cEntity_ID_View::execution_space>;
  const Mesh_type& m = *mesh_;
  Kokkos::parallel_for(
    "MeshAudit::checkCellToFaces", Range_type(0, ncells_all_), KOKKOS_LAMBDA(const Entity_ID c) {
      auto [cface, fdirs] = m.getCellFacesAndDirections(c);
      for (auto& fd : fdirs) {
        if (fd != 1 && fd != -1) bad_cells(c) = 1;
      }
    });

  bool error = checkErrorList(bad_cells, "inadmissable or no data for cells", this->os_);
  return globalAny_(error);
}


// Check that cells are not topologically degenerate (repeated node index).
// If this test fails, many further tests involving the cell_to_nodes data
// should be avoided.
//
template <class Mesh_type>
bool
MeshAudit_Geometry<Mesh_type>::checkCellDegeneracy() const
{
  os_ << "Checking cells for topological degeneracy ..." << std::endl;

  Entity_ID_View bad_cells("bad cells", ncells_all_);

  const Mesh_type& m = *mesh_;
  using Range_type = Kokkos::RangePolicy<LO, typename cEntity_ID_View::execution_space>;
  Kokkos::parallel_for(
    "MeshAudit::checkCellDegeneracy",
    Range_type(0, ncells_all_),
    KOKKOS_CLASS_LAMBDA(const Entity_ID c) {
      auto cnodes = m.getCellNodes(c); // should not fail
      if (!this->areDistinctValues_(cnodes)) bad_cells(c) = 1;
    });

  bool error = checkErrorList(bad_cells, "found topologically degenerate cells", this->os_);
  return globalAny_(error);
}

// Check that cell_to_faces and face_to_nodes are returning the correct
// values by composing those maps and comparing the result against the
// result returned by cell_to_nodes and the local face numbering convention
// described by cell_topology::HexFaceVert.  Also check that the relative
// orientation value returned by cell_to_face_dirs is correct.  If this test
// fails, one or more of the methods cell_to_faces, face_to_nodes, and
// cell_to_face_dirs are returning incorrect results and that further tests
// using their values should be avoided.
template <class Mesh_type>
bool
MeshAudit_Geometry<Mesh_type>::checkCellToFacesToNodes() const
{
  // This can't be done on device unless we move Topology onto device...
  if constexpr (std::is_same_v<Mesh_type, AmanziMesh::Mesh>) {
    return false;
  } else {
    // unordered meshes cannot pass this test.
    if (!mesh_->isOrdered()) {
      os_ << "  Skipped (unordered mesh)" << std::endl;
      return false;
    }

    Entity_ID_View bad_cells("bad cells", ncells_all_);
    Entity_ID_View bad_cells1("bad cells1", ncells_all_);

    for (LO c = 0; c != ncells_all_; ++c) {
      Cell_kind ctype_enum = mesh_->getCellType(c);

      // If this is a general, non-standard element there is nothing to
      // to check against
      if (ctype_enum == Cell_kind::UNKNOWN || ctype_enum == Cell_kind::POLYGON ||
          ctype_enum == Cell_kind::POLYHED)
        continue;

      int ctype = (int)ctype_enum;
      auto cnode = mesh_->getCellNodes(c);
      auto [cface, fdirs] = mesh_->getCellFacesAndDirections(c);

      bool bad_face = false;
      bool bad_dir = false;

      if (cface.size() != Topology::nface_std[ctype]) {
        bad_face = true;
      } else {
        for (int k = 0; k < cface.size(); ++k) {
          auto fnode = mesh_->getFaceNodes(cface(k)); // this should not fail
          int nfn = Topology::nfnodes_std[ctype][k];

          if (fnode.size() != nfn) {
            bad_face = true;
            break;
          }

          typename Mesh_type::Entity_ID_View fnode_ref("fnode_red", nfn);
          for (int i = 0; i < nfn; ++i) {
            int nodenum = Topology::fnodes_std[ctype][k][i];
            fnode_ref[i] = cnode(nodenum);
          }

          int dir = this->isSameFace_(fnode, fnode_ref); // should be the same face
          if (dir == 0) {                                // wrong face
            bad_face = true;
            break;
          } else if (dir != fdirs(k)) { // right face but wrong dir value
            bad_dir = true;
            break;
          }
        }
      }
      if (bad_face) bad_cells(c) = 1;
      if (bad_dir) bad_cells1(c) = 1;
    }

    bool error = checkErrorList(bad_cells, "bad getCellFaces values for cells", this->os_);

    // algorithms on a non-manifold use multiple normals and special continuity equations
    // for fluxes, so that orientation does not play role. This may change in the future.
    if (mesh_->getSpaceDimension() == mesh_->getManifoldDimension()) {
      error |= checkErrorList(bad_cells1, "bad cell_get_face_dirs values for cells", this->os_);
    }
    return globalAny_(error);
  }
}


// Check that getNodeCoordinate successfully returns data for all nodes
// A negative return value signals a terminal error
//
template <class Mesh_type>
bool
MeshAudit_Geometry<Mesh_type>::checkNodeToCoordinates() const
{
  int spdim = mesh_->getSpaceDimension();
  Entity_ID_View bad_nodes("bad nodes", nnodes_all_);

  using Range_type = Kokkos::RangePolicy<LO, typename cEntity_ID_View::execution_space>;
  const Mesh_type& m = *mesh_;
  Kokkos::parallel_for(
    "MeshAudit::checkNodeToCoordinates", Range_type(0, nnodes_all_), KOKKOS_LAMBDA(const LO n) {
      auto x0 = m.getNodeCoordinate(n); // this may fail
      if (x0.dim() != spdim) bad_nodes(n) = 1;
    });

  bool error = checkErrorList(bad_nodes, "missing coordinate data for nodes", this->os_);
  return globalAny_(error);
}


// Check that cell_get_coordinates successfully returns data for all cells, and
// that the data is identical to that returned by composing with cell_to_nodes.
// If this test fails, further tests using cell_to_coordinates data should be
// avoided.
//
template <class Mesh_type>
bool
MeshAudit_Geometry<Mesh_type>::checkCellToCoordinates() const
{
  int spdim = mesh_->getSpaceDimension();
  Entity_ID_View bad_cells("bad cells", ncells_all_);

  using Range_type = Kokkos::RangePolicy<LO, typename cEntity_ID_View::execution_space>;
  const Mesh_type& m = *mesh_;
  Kokkos::parallel_for(
    "MeshAudit::checkNodeToCoordinates", Range_type(0, ncells_all_), KOKKOS_LAMBDA(const LO c) {
      auto x0 = m.getCellCoordinates(c); // this may fail
      cEntity_ID_View cnode;
      m.getCellNodes(c, cnode); // this should not fail
      for (int k = 0; k < cnode.size(); ++k) {
        bool bad_data = false;
        auto xref = m.getNodeCoordinate(cnode(k)); // this should not fail
        for (int i = 0; i < spdim; i++)
          if (x0(k)[i] != xref[i]) {
            bad_data = true;
            bad_cells(c) = 1;
            break;
          }
        if (bad_data) break;
      }
    });

  bool error = checkErrorList(bad_cells, "bad cell_to_coordinates data for cells", this->os_);
  return globalAny_(error);
}


// Check that face_get_coordinates successfully returns data for all faces, and
// that the data is identical to that returned by composing with face_to_nodes.
// If this test fails, further tests using face_to_coordinates data should be
// avoided.
//
template <class Mesh_type>
bool
MeshAudit_Geometry<Mesh_type>::checkFaceToCoordinates() const
{
  int spdim = mesh_->getSpaceDimension();
  Entity_ID_View bad_faces("bad faces", nfaces_all_);

  using Range_type = Kokkos::RangePolicy<LO, typename cEntity_ID_View::execution_space>;
  const Mesh_type& m = *mesh_;
  Kokkos::parallel_for(
    "MeshAudit::checkNodeToCoordinates", Range_type(0, nfaces_all_), KOKKOS_LAMBDA(const LO c) {
      auto x0 = m.getFaceCoordinates(c); // this may fail
      cEntity_ID_View cnode;
      m.getFaceNodes(c, cnode); // this should not fail
      for (int k = 0; k < cnode.size(); ++k) {
        bool bad_data = false;
        auto xref = m.getNodeCoordinate(cnode(k)); // this should not fail
        for (int i = 0; i < spdim; i++)
          if (x0(k)[i] != xref[i]) {
            bad_data = true;
            bad_faces(c) = 1;
            break;
          }
        if (bad_data) break;
      }
    });

  bool error = checkErrorList(bad_faces, "bad face_to_coordinates data for faces", this->os_);
  return globalAny_(error);
}


//
// Check that upward and downward adjacencies and orientations are consistent
// between cells and faces.
//
template <class Mesh_type>
bool
MeshAudit_Geometry<Mesh_type>::checkFaceCellAdjacencyConsistency() const
{
  Entity_ID_View bad_faces("bad faces", nfaces_all_);
  Entity_ID_View bad_faces1("bad faces1", nfaces_all_);

  const Mesh_type& m = *mesh_;
  using Range_type = Kokkos::RangePolicy<LO, typename cEntity_ID_View::execution_space>;
  Kokkos::parallel_for(
    "MeshAudit::checkFaceCellAdjacencyConsistency",
    Range_type(0, nfaces_all_),
    KOKKOS_LAMBDA(const LO f) {
      cEntity_ID_View fcells;
      m.getFaceCells(f, fcells);
      for (const auto& c : fcells) {
        cEntity_ID_View cfaces;
        typename Mesh_type::cDirection_View cfdirs;
        m.getCellFacesAndDirs(c, cfaces, &cfdirs);

        // find the match
        int i = 0;
        for (; i != cfaces.size(); ++i) {
          if (f == cfaces[i]) { break; }
        }
        if (i == cfaces.size()) {
          bad_faces[f] = 1;
          break;
        }

        // check directions
        int orientation = 0;
        auto fnormal = m.getFaceNormal(f, c, &orientation);
        if (orientation == 0) {
          if (m.getSpaceDimension() == m.getManifoldDimension()) {
            // should be an orientation but there isn't
            bad_faces1[f] = 1;
            break;
          }
        } else {
          if (orientation != cfdirs[i]) {
            // orientation is there, not equivalent to the dir from the cell
            bad_faces1[f] = 1;
            break;
          }
        }
      }
    });

  bool error = checkErrorList(bad_faces, "inconsistent adjacencies", this->os_);
  error |= checkErrorList(bad_faces1, "inconsistent orientations", this->os_);
  return globalAny_(error);
}


//
// Check that a face normal, relative to a cell, is outward.  Note this
// requires star-convex meshes.
//
template <class Mesh_type>
bool
MeshAudit_Geometry<Mesh_type>::checkFaceNormalReltoCell() const
{
  Entity_ID_View bad_faces("bad faces", nfaces_all_);

  const Mesh_type& m = *mesh_;
  using Range_type = Kokkos::RangePolicy<LO, typename cEntity_ID_View::execution_space>;
  Kokkos::parallel_for(
    "MeshAudit::checkFaceNormalReltoCell",
    Range_type(0, nfaces_all_),
    KOKKOS_CLASS_LAMBDA(const LO f) {
      auto fc = m.getFaceCentroid(f);

      cEntity_ID_View fcells;
      m.getFaceCells(f, fcells);
      for (int k = 0; k < fcells.size(); ++k) {
        AmanziGeometry::Point cc = m.getCellCentroid(fcells[k]);
        AmanziGeometry::Point fnormal = m.getFaceNormal(f, fcells[k]);
        if ((fc - cc) * fnormal < 0) {
          bad_faces[f] = 1;
          break;
        }
      }
    });

  bool error = checkErrorList(bad_faces, "bad face_normal, not outward", this->os_);
  return globalAny_(error);
}


//
// Check that a face normal has the correct orientation.
//
template <class Mesh_type>
bool
MeshAudit_Geometry<Mesh_type>::checkFaceNormalOrientation() const
{
  // this test is only for meshes which have a natural normal
  if (mesh_->getSpaceDimension() != mesh_->getManifoldDimension()) return false;

  Entity_ID_View bad_faces("bad faces", nfaces_all_);

  const Mesh_type& m = *mesh_;
  using Range_type = Kokkos::RangePolicy<LO, typename cEntity_ID_View::execution_space>;
  Kokkos::parallel_for(
    "MeshAudit::checkFaceNormalOrientation", Range_type(0, nfaces_all_), KOKKOS_LAMBDA(const LO f) {
      AmanziGeometry::Point fnormal = m.getFaceNormal(f);

      cEntity_ID_View fcells;
      m.getFaceCells(f, fcells);
      for (int k = 0; k < fcells.size(); ++k) {
        int orientation = 0;
        AmanziGeometry::Point fnormalc = m.getFaceNormal(f, fcells[k], &orientation);
        if (orientation != 0) {
          if (AmanziGeometry::norm((orientation * fnormal) - fnormalc) > 1.e-16) {
            bad_faces[f] = 1;
            break;
          }
        }
      }
    });

  bool error = checkErrorList(bad_faces, "bad face_normal, inconsistent orientation", this->os_);
  return globalAny_(error);
}


// The cells must not be degenerate, either topologically (repeated
// node index) or geometrically (with coincident nodes).
//
template <class Mesh_type>
bool
MeshAudit_Geometry<Mesh_type>::checkCellGeometry() const
{
  os_ << "Checking cell geometry ..." << std::endl;
  Entity_ID_View bad_cells("bad cells", ncells_all_);

  const Mesh_type& m = *mesh_;
  using Range_type = Kokkos::RangePolicy<LO, typename cEntity_ID_View::execution_space>;
  Kokkos::parallel_for(
    "MeshAudit::checkCellGeometry", Range_type(0, ncells_all_), KOKKOS_LAMBDA(const LO c) {
      double hvol = m.getCellVolume(c);
      if (hvol <= 0.0) bad_cells[c] = 1;
    });
  bool error = checkErrorList(bad_cells, "found cells with non-positive volumes", this->os_);
  return globalAny_(error);
}

// We assume the maps returned by the map accessors are well-formed as
// Epetra_maps.  Here we are checking that they have the required
// characteristics, especially for the relationship between the owned
// and overlapped version of the maps.  Specifically, we require that
// the index base for all maps is 0, and that the owned maps (ones without
// ghost) are 1-1 maps.  In serial the owned and overlap (one with ghosts)
// maps should be the same.  In parallel, the overlap version should extend
// the owned version (owned.GID() = overlapped.GID() for all owned LIDs),
// and that the remaining overlapped LIDs (if any) refer to GIDs owned
// by other processes.  In addition there should be no duplicate GIDs on
// any processes map.
//
template <class Mesh_type>
bool
MeshAudit_Maps<Mesh_type>::checkNodeMaps() const
{
  return checkMaps_(*mesh_->getMap(AmanziMesh::Entity_kind::NODE, false),
                    *mesh_->getMap(AmanziMesh::Entity_kind::NODE, true));
}

template <class Mesh_type>
bool
MeshAudit_Maps<Mesh_type>::checkFaceMaps() const
{
  return checkMaps_(*mesh_->getMap(AmanziMesh::Entity_kind::FACE, false),
                    *mesh_->getMap(AmanziMesh::Entity_kind::FACE, true));
}

template <class Mesh_type>
bool
MeshAudit_Maps<Mesh_type>::checkCellMaps() const
{
  return checkMaps_(*mesh_->getMap(AmanziMesh::Entity_kind::CELL, false),
                    *mesh_->getMap(AmanziMesh::Entity_kind::CELL, true));
}

template <class Mesh_type>
bool
MeshAudit_Maps<Mesh_type>::checkMaps_(const Map_type& map_own, const Map_type& map_use) const
{
  bool error = false;

  // Local index should start at 0.
  if (map_own.getIndexBase() != 0) {
    os_ << "ERROR: the owned map's index base is not 0." << std::endl;
    error = true;
  }
  if (map_use.getIndexBase() != 0) {
    os_ << "ERROR: the overlap map's index base is not 0." << std::endl;
    error = true;
  }

  // // Check that the owned map is 1-1.
  // if (!map_own.UniqueGIDs()) {
  //   os_ << "ERROR: owned map is not 1-to-1" << std::endl;
  //   error = true;
  // }

  // Check that the global ID space is contiguously divided (but not
  // necessarily uniformly) across all processers
  // This is not actually required, particularly on e.g. extraced meshes.
  // std::vector<int> owned_GIDs(map_own.getLocalNumElements());
  // for (int i = 0; i < map_own.getLocalNumElements(); i++)
  //   owned_GIDs[i] = map_own.getGlobalElement(i);
  // std::sort(owned_GIDs.begin(), owned_GIDs.end());

  // for (int i = 0; i < ((int)map_own.getLocalNumElements() - 1); i++) {
  //   int diff = owned_GIDs[i + 1] - owned_GIDs[i];
  //   if (diff > 1) {
  //     os_ << "ERROR: owned map is not contiguous" << std::endl;
  //     os_ << "Global IDs jump from " << owned_GIDs[i] << " to " << owned_GIDs[i + 1] << std::endl;
  //     error = true;
  //   }
  // }

  error = globalAny_(error);
  if (error) return error;

  if (this->comm_->getSize() == 1) {
    // Serial or 1-process MPI
    if (!map_use.isSameAs(map_own)) {
      os_ << "ERROR: the overlap map differs from the owned map (single process)." << std::endl;
      error = true;
    }
    return error;

  } else {
    // Multi-process MPI
    int num_own = map_own.getLocalNumElements();
    int num_use = map_use.getLocalNumElements();
    int num_ovl = num_use - num_own;

    // Verify that the used map extends the owned map.
    bool bad_map = false;
    if (num_ovl < 0)
      bad_map = true;
    else {
      for (int j = 0; j < num_own; ++j)
        if (map_use.getGlobalElement(j) != map_own.getGlobalElement(j)) bad_map = true;
    }
    if (bad_map) {
      os_ << "ERROR: overlap map does not extend the owned map." << std::endl;
      error = true;
    }

    error = globalAny_(error);
    if (error) return error;

    // Verify that the overlap indices are owned by other processes.
    Teuchos::Array<Entity_GID> gids(num_ovl);
    Teuchos::Array<Entity_ID> pids(num_ovl);
    Teuchos::Array<Entity_ID> lids(num_ovl);
    for (int j = 0; j < num_ovl; ++j) gids[j] = map_use.getGlobalElement(j + num_own);
    map_own.getRemoteIndexList(gids, pids, lids);
    bad_map = false;
    for (int j = 0; j < num_ovl; ++j)
      if (pids[j] < 0 || pids[j] == this->comm_->getRank()) bad_map = true;
    if (bad_map) {
      os_ << "ERROR: invalid ghosts in overlap map." << std::endl;
      error = true;
    }

    // Look for duplicates among the overlap indices.
    std::sort(gids.begin(), gids.end());
    if (adjacent_find(gids.begin(), gids.end()) != gids.end()) {
      os_ << "ERROR: duplicate ghosts in overlap map." << std::endl;
      error = true;
    }
    return globalAny_(error);
  }
}

// Check that ghost nodes are exact copies of their master.
// This simply means that they have the same coordinates.
//
template <class Mesh_type>
bool
MeshAudit_Maps<Mesh_type>::checkNodeToCoordinatesGhostData() const
{
  int spdim = mesh_->getSpaceDimension();

  const auto& node_map_own = mesh_->getMap(AmanziMesh::Entity_kind::NODE, false);
  const auto& node_map_use = mesh_->getMap(AmanziMesh::Entity_kind::NODE, true);

  int nnode_own = node_map_own->getLocalNumElements();
  int nnode_use = node_map_use->getLocalNumElements();

  MultiVector_type coord_use(node_map_use, spdim);
  auto coord_own = coord_use.offsetViewNonConst(node_map_own, 0);
  if constexpr (std::is_same_v<Mesh_type, Mesh>) {
    auto coord_own_view = coord_own->getLocalViewDevice(Tpetra::Access::ReadWrite);
    const Mesh_type& m = *mesh_;
    Kokkos::parallel_for(
      "MeshAudit::check_node_to_coordinates_ghost_data", nnode_own, KOKKOS_LAMBDA(const int n) {
        auto coord = m.getNodeCoordinate(n);
        for (int k = 0; k < spdim; ++k) coord_own_view(n, k) = coord[k];
      });
  } else {
    auto coord_own_view = coord_own->getLocalViewHost(Tpetra::Access::ReadWrite);
    for (size_t n = 0; n != nnode_own; ++n) {
      auto coord = mesh_->getNodeCoordinate(n);
      for (int k = 0; k < spdim; ++k) coord_own_view(n, k) = coord[k];
    }
  }

  Import_type importer(node_map_own, node_map_use);
  coord_use.doImport(*coord_own, importer, Tpetra::INSERT);

  Entity_ID_View bad_nodes("bad nodes", nnode_use - nnode_own);
  if constexpr (std::is_same_v<Mesh_type, Mesh>) {
    auto coord_use_view = coord_use.getLocalViewDevice(Tpetra::Access::ReadOnly);
    const Mesh_type& m = *mesh_;
    Kokkos::parallel_for(
      "MeshAudit::check_node_to_coordinates_ghost_data",
      nnode_use - nnode_own,
      KOKKOS_LAMBDA(const int n) {
        auto coord = m.getNodeCoordinate(n + nnode_own);
        for (int k = 0; k < spdim; ++k) {
          if (coord[k] != coord_use_view(n + nnode_own, k)) { bad_nodes[n] = 1; }
        }
      });
  } else {
    auto coord_use_view = coord_use.getLocalViewHost(Tpetra::Access::ReadOnly);
    for (size_t n = 0; n != (nnode_use - nnode_own); ++n) {
      auto coord = mesh_->getNodeCoordinate(n + nnode_own);
      for (int k = 0; k < spdim; ++k) {
        if (coord[k] != coord_use_view(n + nnode_own, k)) { bad_nodes[n] = 1; }
      }
    }
  }

  bool error = checkErrorList(bad_nodes, "found ghost nodes with incorrect coordinates", this->os_);
  return globalAny_(error);
}

// Check that ghost faces are exact copies of their master.  This means
// that the the GIDs of the nodes defining the face are the same, including
// their order (face orientation).
//
template <class Mesh_type>
bool
MeshAudit_Maps<Mesh_type>::checkFaceToNodesGhostData() const
{
  os_ << "WARNING: check_face_to_nodes_ghost_data() test disabled in tpetra" << std::endl;
  return false;
  // const Map_type &node_map = mesh_->getMap(AmanziMesh::Entity_kind::NODE, true);
  // const Map_type &face_map_own = mesh_->getMap(AmanziMesh::Entity_kind::FACE, false);
  // const Map_type &face_map_use = mesh_->getMap(AmanziMesh::Entity_kind::FACE, true);

  // int nface_own = face_map_own->getLocalNumElements();
  // int nface_use = face_map_use->getLocalNumElements();

  // Entity_ID_List fnode;
  // Entity_ID_List bad_faces, bad_faces1, bad_faces2;

  // // Create a matrix of the GIDs for all owned faces.

  // int maxnodes = 0;
  // for (Entity_ID j = 0; j < nface_own; ++j) {
  //   mesh_->getFaceNodes(j, fnode);
  //   maxnodes = fnode.size() > maxnodes ? fnode.size() : maxnodes;
  // }

  // // parallel calls require uniformity
  // int maxnodes_tmp(maxnodes);
  // this->comm_->MaxAll(&maxnodes_tmp, &maxnodes, 1);

  // Epetra_IntSerialDenseMatrix gids(nface_use,maxnodes); // no Epetra_IntMultiVector :(
  // for (Entity_ID j = 0; j < nface_own; ++j) {
  //   mesh_->getFaceNodes(j, fnode);
  //   for (int k = 0; k < fnode.size(); ++k)
  //     gids(j,k) = node_map.getGlobalElement(fnode[k]);
  //   for (int k = fnode.size(); k < maxnodes; ++k)
  //     gids(j, k) = 0;
  // }

  // // Import these GIDs to all used faces; sets values on ghost faces.
  // Import_type importer(face_map_use, face_map_own);
  // for (int k = 0; k < maxnodes; ++k) {
  //   Epetra_IntVector kgids_own(View, face_map_own, gids[k]);
  //   Epetra_IntVector kgids_use(View, face_map_use, gids[k]);
  //   kgids_use.Import(kgids_own, importer, Insert);
  // }

  // // Compare the ghost face GIDs against the reference values just computed.
  // for (Entity_ID j = nface_own; j < nface_use; ++j) {
  //   mesh_->getFaceNodes(j, fnode);
  //   bool bad_data = false;
  //   for (int k = 0; k < fnode.size(); ++k)
  //     if (node_map.getGlobalElement(fnode[k]) != gids(j,k)) {
  //       bad_data = true;
  //     }
  //   if (bad_data) {
  //     // Determine just how bad the data is.
  //     Entity_ID_List fnode_ref(maxnodes);
  //     for (int k = 0; k < maxnodes; ++k) {
  //       fnode[k] = node_map.getGlobalElement(fnode[k]);
  //       fnode_ref[k] = gids(j,k);
  //     }
  //     int n = this->isSameFace_(fnode, fnode_ref);
  //     switch (n) {
  //       case 0: // completely bad -- different face
  //         bad_faces.push_back(j);

  //         os_ << "P " << this->comm_->getRank() << ": ghost face " << j << " (GID "
  //                   << face_map_use.getGlobalElement(j) << ")," <<
  //           " has different nodes than its master " << std::endl;
  //         os_ << "ghost face nodes (GIDs): ";
  //         for (int k = 0; k < fnode.size(); ++k)
  //           os_ << fnode[k] << " ";
  //         os_ << std::endl;
  //         os_ << "master face nodes (GIDs): ";
  //         for (int k = 0; k < fnode_ref.size(); ++k)
  //           os_ << fnode_ref[k] << " ";
  //         os_ << std::endl;
  //         break;
  //       case -1: // very bad -- same face but wrong orientation
  //         bad_faces1.push_back(j);

  //         os_ << "P " << this->comm_->getRank() << ": ghost face " << j << ", has different orientation than its master " << std::endl;
  //         os_ << "ghost face nodes (GIDs): ";
  //         for (int k = 0; k < fnode.size(); ++k)
  //           os_ << fnode[k];
  //         os_ << std::endl;
  //         os_ << "master face nodes (GIDs): ";
  //         for (int k = 0; k < fnode_ref.size(); ++k)
  //           os_ << fnode_ref[k];
  //         os_ << std::endl;
  //         break;
  //       case 1:
  //         // This is fine because there is no way to ensure this for
  //         // general meshes (think of building a tet mesh and defining
  //         // the faces such that they return the same vertices in the
  //         // same order (unpermuted) no matter which way the elements
  //         // match up. You can only ensure this for hexahedral meshes
  //         // As long as the faces have the same vertices and have the
  //         // same direction, but the vertices are permuted, its fine
  //         break;
  //     }
  //   }
  // }

  // bool error = false;

  // if (!bad_faces.empty()) {
  //   os_ << "ERROR: found bad data for ghost faces:";
  //   checkErrorList(bad_faces, this->os_);
  //   error = true;
  // }

  // if (!bad_faces1.empty()) {
  //   os_ << "ERROR: found mis-oriented ghost faces:";
  //   checkErrorList(bad_faces1, this->os_);
  //   error = true;
  // }

  // return globalAny_(error);
}

// Check that ghost cells are exact copies of their master.  This means
// that the the GIDs of the nodes defining the cell are the same, including
// their order (face orientation).
//
template <class Mesh_type>
bool
MeshAudit_Maps<Mesh_type>::checkCellToNodesGhostData() const
{
  os_ << "WARNING: check_cell_to_nodes_ghost_data() test disabled in tpetra" << std::endl;
  return false;
  // const Map_type &node_map = mesh_->getMap(AmanziMesh::Entity_kind::NODE, true);
  // const Map_type &cell_map_own = mesh_->getMap(AmanziMesh::Entity_kind::CELL, false);
  // const Map_type &cell_map_use = mesh_->getMap(AmanziMesh::Entity_kind::CELL, true);

  // int ncell_own = cell_map_own->getLocalNumElements();
  // int ncell_use = cell_map_use->getLocalNumElements();

  // Entity_ID_List cnode;
  // Entity_ID_List bad_cells;

  // int maxnodes = 0;
  // for (Entity_ID j = 0; j < ncell_own; ++j) {
  //   mesh_->getCellNodes(j, cnode);
  //   maxnodes = (cnode.size() > maxnodes) ? cnode.size() : maxnodes;
  // }

  // // parallel calls require uniformity
  // int maxnodes_tmp(maxnodes);
  // this->comm_->MaxAll(&maxnodes_tmp, &maxnodes, 1);

  // Epetra_IntSerialDenseMatrix gids(ncell_use,maxnodes); // no Epetra_IntMultiVector :(

  // // Create a matrix of the GIDs for all owned cells.
  // for (Entity_ID j = 0; j < ncell_own; ++j) {
  //   mesh_->getCellNodes(j, cnode);
  //   for (int k = 0; k < cnode.size(); ++k)
  //     gids(j,k) = node_map.getGlobalElement(cnode[k]);
  //   for (int k = cnode.size(); k < maxnodes; ++k)
  //     gids(j,k) = 0;
  // }

  // // Import these GIDs to all used cells; sets values on ghost cells.
  // Import_type importer(cell_map_use, cell_map_own);
  // for (int k = 0; k < maxnodes; ++k) {
  //   Epetra_IntVector kgids_own(View, cell_map_own, gids[k]);
  //   Epetra_IntVector kgids_use(View, cell_map_use, gids[k]);
  //   kgids_use.Import(kgids_own, importer, Insert);
  // }

  // // Compare the ghost cell GIDs against the reference values just computed.
  // for (Entity_ID j = ncell_own; j < ncell_use; ++j) {
  //   mesh_->getCellNodes(j, cnode);
  //   bool bad_data = false;
  //   for (int k = 0; k < cnode.size(); ++k)
  //     if (node_map.getGlobalElement(cnode[k]) != gids(j,k)) bad_data = true;
  //   if (bad_data) {
  //     for (int k = 0; k < cnode.size(); ++k) {
  //       bool found = false;
  //       for (int l = 0; l < cnode.size(); ++l)
  //         if (node_map.getGlobalElement(cnode[k]) == gids(j,l)) {
  //           found = true;
  //           break;
  //         }
  //       if (!found) {
  //         bad_cells(j) = 1;
  //         break;
  //       }
  //     }
  //   }
  // }

  // bool error = false;

  // if (!bad_cells.empty()) {
  //   os_ << "ERROR: found bad data for ghost cells:";
  //   checkErrorList(bad_cells, this->os_);
  //   os_ << "       The ghost cells don't have the same nodes as their respective masters." << std::endl;
  //   error = true;
  // }

  // return globalAny_(error);
}


// Check that ghost cells reference the same faces, and in the same
// order, as their master.  This is part of a ghost being an exact
// copies the master.  Even if the preceding ghost checks pass it is
// still possible for this one to fail.  In this case the ghost would
// reference a different face (by GID) but that face would reference
// the same nodes as the correct face.  So the two faces would be
// geometrically identical, including orientation, but be distinct.
//
template <class Mesh_type>
bool
MeshAudit_Maps<Mesh_type>::checkCellToFacesGhostData() const
{
  os_ << "WARNING: check_face_to_faces_ghost_data() test disabled in tpetra" << std::endl;
  return false;

  // const Map_type &face_map = mesh_->getMap(AmanziMesh::Entity_kind::FACE, true);
  // const Map_type &cell_map_own = mesh_->getMap(AmanziMesh::Entity_kind::CELL, false);
  // const Map_type &cell_map_use = mesh_->getMap(AmanziMesh::Entity_kind::CELL, true);

  // int ncell_own = cell_map_own->getLocalNumElements();
  // int ncell_use = cell_map_use->getLocalNumElements();

  // Entity_ID_List cface;
  // Entity_ID_List bad_cells;

  // // Create a matrix of the GIDs for all owned cells.
  // int maxfaces = 0;
  // for (Entity_ID j = 0; j < ncell_own; ++j) {
  //   mesh_->getCellFaces(j, cface);
  //   maxfaces = (cface.size() > maxfaces) ? cface.size() : maxfaces;
  // }

  // // parallel calls require uniformity
  // int maxfaces_tmp(maxfaces);
  // this->comm_->MaxAll(&maxfaces_tmp, &maxfaces, 1);

  // Epetra_IntSerialDenseMatrix gids(ncell_use,maxfaces); // no Epetra_IntMultiVector :(

  // for (Entity_ID j = 0; j < ncell_own; ++j) {
  //   mesh_->getCellFaces(j, cface);
  //   for (int k = 0; k < cface.size(); ++k)
  //     gids(j,k) = face_map.getGlobalElement(cface[k]);
  //   for (int k = cface.size(); k < maxfaces; ++k)
  //     gids(j,k) = 0;
  // }

  // // Import these GIDs to all used cells; sets values on ghost cells.
  // Import_type importer(cell_map_use, cell_map_own);
  // for (int k = 0; k < maxfaces; ++k) {
  //   Epetra_IntVector kgids_own(View, cell_map_own, gids[k]);
  //   Epetra_IntVector kgids_use(View, cell_map_use, gids[k]);
  //   kgids_use.Import(kgids_own, importer, Insert);
  // }

  // // Compare the ghost cell GIDs against the reference values just computed.
  // for (Entity_ID j = ncell_own; j < ncell_use; ++j) {
  //   mesh_->getCellFaces(j, cface);
  //   bool bad_data = false;
  //   for (int k = 0; k < cface.size(); ++k)
  //     if (face_map.getGlobalElement(cface[k]) != gids(j,k)) bad_data = true;
  //   if (bad_data) bad_cells(j) = 1;
  // }

  // bool error = false;

  // if (!bad_cells.empty()) {
  //   os_ << "ERROR: found bad data for ghost cells:";
  //   checkErrorList(bad_cells, this->os_);
  //   os_ << "       The ghost cells are not exact copies of their master." << std::endl;
  //   error = true;
  // }

  // return globalAny_(error);
}

////////////////////////////////////////////////////////////////////////////////
//
// TESTS OF SET DATA
//
////////////////////////////////////////////////////////////////////////////////


// Check that get_set_ids successfully returns the vector of set IDs, without
// duplicates, and that each process gets the exact same vector of set IDs.
// This is a collective test, returning a collective pass/fail result.
//
template <class Mesh_type>
bool
MeshAudit_Sets<Mesh_type>::checkNodeSetIDs() const
{
  return checkGetSetIDs_(Entity_kind::NODE);
}

template <class Mesh_type>
bool
MeshAudit_Sets<Mesh_type>::checkFaceSetIDs() const
{
  return checkGetSetIDs_(Entity_kind::FACE);
}

template <class Mesh_type>
bool
MeshAudit_Sets<Mesh_type>::checkCellSetIDs() const
{
  return checkGetSetIDs_(Entity_kind::CELL);
}

template <class Mesh_type>
bool
MeshAudit_Sets<Mesh_type>::checkGetSetIDs_(Entity_kind kind) const
{
  // Previously we checked if set entities were valid, with no duplicate entries.
  //
  // Since we lazily resolve sets on demand, there is no real way to know what
  // combination of region and mesh the user actually intends to use as a set.
  return false;
}

// Check that valid_set_id() returns the correct results.
// This is a collective test, returning a collective pass/fail result.
//
template <class Mesh_type>
bool
MeshAudit_Sets<Mesh_type>::checkValidNodeSetID() const
{
  return checkValidSetID_(Entity_kind::NODE);
}

template <class Mesh_type>
bool
MeshAudit_Sets<Mesh_type>::checkValidFaceSetID() const
{
  return checkValidSetID_(Entity_kind::FACE);
}

template <class Mesh_type>
bool
MeshAudit_Sets<Mesh_type>::checkValidCellSetID() const
{
  return checkValidSetID_(Entity_kind::CELL);
}

template <class Mesh_type>
bool
MeshAudit_Sets<Mesh_type>::checkValidSetID_(Entity_kind kind) const
{
  // This needs some work.  Since we lazily resolve sets on demand, there is no
  // real way to know what sets the user actually intends to use.
  return false;

  // // Get the list of set IDs.
  // int nset = mesh_->num_sets(kind); // this should not fail
  // std::vector<Set_ID> sids(nset);
  // mesh_->get_set_ids(kind, &sids); // this should not fail

  // std::vector<Set_ID> bad_sids;
  // int max_id = 0;
  // for (int j = 0; j < nset; ++j)
  //   if (max_id < sids[j]) max_id = sids[j];
  // std::vector<bool> valid(max_id+2, false);
  // for (int j = 0; j < nset; ++j)
  //   valid[sids[j]] = true;

  // std::vector<Set_ID> bad_sids1, bad_sids2;
  // for (int n = 0; n < valid.size(); ++n) {
  //   if (valid[n] && !mesh_->valid_set_id(n, kind)) bad_sids1.push_back(n);
  //   if (!valid[n] && mesh_->valid_set_id(n, kind)) bad_sids2.push_back(n);
  // }

  // bool error = false;
  // if (!bad_sids1.empty()) {
  //   os_ << "ERROR: valid_set_id() returned false for valid set IDs:";
  //   checkErrorList(bad_sids1, this->os_);
  //   error = true;
  // }
  // if (!bad_sids2.empty()) {
  //   os_ << "ERROR: valid_set_id() returned true for invalid set IDs:";
  //   checkErrorList(bad_sids2, this->os_);
  //   error = true;
  // }
  // return globalAny_(error);
}

// For each set, check that get_set successfully returns valid references to
// local entities, without duplicates, and that the used set is consistent
// with the owned set.
//
template <class Mesh_type>
bool
MeshAudit_Sets<Mesh_type>::checkNodeSets() const
{
  return checkSets_(Entity_kind::NODE,
                    *mesh_->getMap(AmanziMesh::Entity_kind::NODE, false),
                    *mesh_->getMap(AmanziMesh::Entity_kind::NODE, true));
}

template <class Mesh_type>
bool
MeshAudit_Sets<Mesh_type>::checkFaceSets() const
{
  return checkSets_(Entity_kind::FACE,
                    *mesh_->getMap(AmanziMesh::Entity_kind::FACE, false),
                    *mesh_->getMap(AmanziMesh::Entity_kind::FACE, true));
}

template <class Mesh_type>
bool
MeshAudit_Sets<Mesh_type>::checkCellSets() const
{
  return checkSets_(Entity_kind::CELL,
                    *mesh_->getMap(AmanziMesh::Entity_kind::CELL, false),
                    *mesh_->getMap(AmanziMesh::Entity_kind::CELL, true));
}

template <class Mesh_type>
bool
MeshAudit_Sets<Mesh_type>::checkSets_(Entity_kind kind,
                                      const Map_type& map_own,
                                      const Map_type& map_use) const
{
  // This needs some work.  Since we lazily resolve sets on demand, there is no
  // real way to know what sets the user actually intends to use.
  return false;

  // bool error = false;
  // // Get the list of set IDs.
  // int nset = mesh_->num_sets(kind);
  // std::vector<Set_ID> sids(nset);
  // mesh_->get_set_ids(kind, &sids);

  // for (int n = 0; n < sids.size(); ++n) {
  //   os_ << "  Checking set ID=" << sids[n] << " ..." << std::endl;

  //   // Basic sanity checks of the owned and used sets.
  //   bool bad_set = check_get_set(sids[n], kind, Parallel_kind::OWNED,
  //       			 map_own) ||
  //                  check_get_set(sids[n], kind, Parallel_kind::ALL,
  //       			 map_use);
  //   bad_set = globalAny_(bad_set);

  //   // Verify the used set relates correctly to the owned set.
  //   if (!bad_set) bad_set = check_used_set(sids[n], kind, map_own, map_use);

  //   // OUGHT TO DO TESTING OF THE GHOST SETS
  //   if (bad_set) error = true;
  // }
  // return error;
}

// Basic sanity check on set values: no duplicates, and all LID values belong
// to the map.  This test runs independently on each process and returns a
// per-process pass/fail result.
//
template <class Mesh_type>
bool
MeshAudit_Sets<Mesh_type>::checkGetSet_(Set_ID sid,
                                        Entity_kind kind,
                                        Parallel_kind ptype,
                                        const Map_type& map) const
{
  os_ << "WARNING: Checks on sets disabled until MeshAudit handles new set specification methods "
         "(Tkt #686)"
      << std::endl;
  return false;

  // // Get the size of the set.
  // try {
  //   std::string set_name = mesh_->geometric_model()->FindRegion(sid)->name();
  //   mesh_->getSetSize(set_name, kind, ptype); // this may fail
  // } catch (...) {
  //   os_ << "  ERROR: caught exception from getSetSize()" << std::endl;
  //   return true;
  // }

  // // Get the set.
  // cEntity_ID_View set;
  // try {
  //   std::string set_name = mesh_->geometric_model()->FindRegion(sid)->name();
  //   mesh_->getSetEntities(set_name, kind, ptype, &set); // this may fail
  // } catch (...) {
  //   os_ << "  ERROR: caught exception from get_set()" << std::endl;
  //   return true;
  // }

  // // Check that all values were assigned.
  // bool bad_data = false;
  // for (int j = 0; j < set.size(); ++j)
  //   if (set[j] == std::numeric_limits<unsigned int>::max()) bad_data = true;
  // if (bad_data) {
  //   os_ << "  ERROR: not all values assigned by get_set()" << std::endl;
  //   return true;
  // }

  // Check that the LIDs in the set belong to the map.
  // os_ << "  WARNING: LIDs in the set check is not implemented yet for Tpetra" << std::endl;
  // Entity_ID_List bad_LIDs;
  // for (int j = 0; j < set.size(); ++j)
  //   if (!map.MyLID(set[j])) bad_LIDs.push_back(set[j]);
  // if (!bad_LIDs.empty()) {
  //   os_ << "  ERROR: set contains invalid LIDs:";
  //   checkErrorList(bad_LIDs, this->os_);
  //   return true;
  // }

  // // Check that there are no duplicates in the set.
  // if (!this->areDistinctValues_(set)) {
  //   os_ << "  ERROR: set contains duplicate LIDs." << std::endl;
  //   // it would be nice to output the duplicates
  //   return true;
  // }

  // return false;
}

// The correct used set is completely determined by the owned set.  This test
// verifies that the used set is what it should be, considering the owned set
// as definitive.  Note that we do not require the vector of used set LIDs to
// extend the vector of owned set LIDs; the values in each list can be
// presented in any order.  This is a collective test, returning a collective
// pass/fail result.
//
template <class Mesh_type>
bool
MeshAudit_Sets<Mesh_type>::checkUsedSet_(Set_ID sid,
                                         Entity_kind kind,
                                         const Map_type& map_own,
                                         const Map_type& map_use) const
{
  os_ << "WARNING: Checks on sets disabled until MeshAudit handles new set specification methods "
         "(Tkt #686)"
      << std::endl;
  return false;

  // std::string set_name = mesh_->geometric_model()->FindRegion(sid)->name();

  // if (this->comm_->getSize() == 1) {
  //   // In serial, the owned and used sets should be identical.

  // int n = mesh_->getSetSize(set_name, kind, Parallel_kind::OWNED);
  // Entity_ID_List set_own;
  // mesh_->getSetEntities(set_name, kind, Parallel_kind::OWNED, &set_own);

  // // Set sizes had better be the same.
  // if (mesh_->getSetSize(set_name, kind, Parallel_kind::ALL) !=
  //     set_own.size()) {
  //   os_ << "  ERROR: owned and used set sizes differ" << std::endl;
  //   return true;
  // }

  // // Verify that the two sets are identical.
  // auto set_use = mesh_->getSetEntities(set_name, kind, Parallel_kind::ALL);
  // bool bad_data = false;
  // for (int j = 0; j < n; ++j)
  //   if (set_use[j] != set_own[j]) bad_data = true;
  // if (bad_data) {
  //   os_ << "  ERROR: owned and used sets differ" << std::endl;
  //   return true;
  // }


  // } else {

  // auto set_own = mesh_->getSetEntities(set_name, kind, Parallel_kind::OWNED);
  // auto set_use = mesh_->getSetEntities(set_name, kind, Parallel_kind::ALL);

  //   // Tag all LIDs in the used map that should belong to the used set;
  //   // the owned set LIDs are taken as definitive.
  //   Epetra_IntVector tag_use(map_use); // fills with zero values
  //   int *tag_data;
  //   tag_use.ExtractView(&tag_data);
  //   Epetra_IntVector tag_own(View, map_own, tag_data);
  //   for (int j = 0; j < set_own.size(); ++j) tag_own[set_own[j]] = 1;
  //   Import_type importer(map_use, map_own);
  //   tag_use.Import(tag_own, importer, Insert);

  //   // Now untag all the LIDs that belong to the used set.  If things
  //   // are correct, the tag vector will be filled with zeros afterwards.
  //   for (int j = 0; j < set_use.size(); ++j) --tag_use[set_use[j]];

  //   bool error = false;

  //   // Check for negative tag values;
  //   // these mark used LIDs that shouldn't be in the set but are.
  //   Entity_ID_List bad_LIDs;
  //   for (int j = 0; j < set_use.size(); ++j)
  //     if (tag_use[j] < 0) bad_LIDs.push_back(j);
  //   if (!bad_LIDs.empty()) {
  //     os_ << "  ERROR: found used LIDs that belong to the set but shouldn't:";
  //     checkErrorList(bad_LIDs, this->os_);
  //     error = true;
  //   }

  //   // Check for positive tag values;
  //   // these mark used LIDs that should be in the set but aren't.
  //   bad_LIDs.resize(0);
  //   for (int j = 0; j < set_own.size(); ++j)
  //     if (tag_use[j] > 0) bad_LIDs.push_back(j);
  //   if (!bad_LIDs.empty()) {
  //     os_ << "  ERROR: found used LIDs that should belong to set but don't:";
  //     checkErrorList(bad_LIDs, this->os_);
  //     error = true;
  //   }

  //   return globalAny_(error);
  // }
}


// SANE PARTITIONING CHECKS.  In parallel the cells, faces and nodes of the
// mesh are each partitioned across the processes, and while in principle
// these partitionings may be completely independent of each other, practical
// considerations lead to certain conditions that reasonable partitions should
// satisfy.  Taking the partitioning of the cells as given, one property that
// should be satisfied by the face and node partitionings is that the process
// that owns a particular face (node) must also own one of the cells containing
// the face (node).
//
template <class Mesh_type>
bool
MeshAudit_Maps<Mesh_type>::checkFacePartition() const
{
  // Mark all the faces contained by owned cells.
  Entity_ID_View all_faces("owned faces", nfaces_all_);
  Kokkos::deep_copy(all_faces, 1);

  const Mesh_type& m = *mesh_;
  using Range_type = Kokkos::RangePolicy<LO, typename cEntity_ID_View::execution_space>;
  Kokkos::parallel_for(
    "MeshAudit::checkFacePartition", Range_type(0, ncells_owned_), KOKKOS_LAMBDA(const LO c) {
      cEntity_ID_View cface;
      m.getCellFaces(c, cface);
      for (auto& f : cface) Kokkos::atomic_store(&all_faces(f), 0);
    });

  auto owned_faces = Kokkos::subview(all_faces, Kokkos::pair{ 0, nfaces_owned_ });
  bool error = checkErrorList(owned_faces, "found orphaned owned faces", this->os_);
  return globalAny_(error);
}


template <class Mesh_type>
bool
MeshAudit_Maps<Mesh_type>::checkNodePartition() const
{
  // Mark all the faces contained by owned cells.
  Entity_ID_View all_nodes("owned nodes", nnodes_all_);
  Kokkos::deep_copy(all_nodes, 1);

  const Mesh_type& m = *mesh_;
  using Range_type = Kokkos::RangePolicy<LO, typename cEntity_ID_View::execution_space>;
  Kokkos::parallel_for(
    "MeshAudit::checkNodePartition", Range_type(0, ncells_owned_), KOKKOS_LAMBDA(const LO c) {
      cEntity_ID_View cnode;
      m.getCellNodes(c, cnode);
      for (auto& f : cnode) Kokkos::atomic_store(&all_nodes(f), 0);
    });

  auto owned_nodes = Kokkos::subview(all_nodes, Kokkos::pair{ 0, nnodes_owned_ });
  bool error = checkErrorList(owned_nodes, "found orphaned owned nodes", this->os_);
  return globalAny_(error);
}


// Returns true if the values in the list are distinct -- no repeats.
template <class Mesh_type>
bool
MeshAudit_Base<Mesh_type>::areDistinctValues_(const cEntity_ID_View& list) const
{
  for (size_t i = 0; i != (list.extent(0) - 1); ++i) {
    for (size_t j = (i + 1); j != list.extent(0); ++j) {
      if (list(i) == list(j)) return false;
    }
  }
  return true;
}


// Returns 1 if the face node lists fnode1 and fnode2 describe the same face
// with the same orientation.  Returns -1 if the lists describe the same face
// but with opposite orientations.  Returns 0 if the lists describe different
// faces.  Implicitly assumes non-degenerate faces; the results are not
// reliable otherwise.
template <class Mesh_type>
template <class View1_type, class View2_type>
int
MeshAudit_Base<Mesh_type>::isSameFace_(const View1_type& fnode1, const View2_type& fnode2) const
{
  int nn = fnode1.size();
  if (nn != fnode2.size()) return 0;

  // Locate position in fnode1 of fnode2(0).
  int i, n;
  for (i = 0, n = -1; i < nn; ++i)
    if (fnode1(i) == fnode2(0)) {
      n = i;
      break;
    }
  if (n == -1) return 0; // did not find it -- different faces

  if (nn == 2) {
    // These are edges in a 2D mesh

    if (n == 0 && fnode1(1) == fnode2(1)) return 1;
    if (n == 1 && fnode1(0) == fnode2(1)) return -1;

    if (n == 0 && fnode1(1) == fnode2(1)) return 1;
    if (n == 1 && fnode1(0) == fnode2(1)) return -1;

  } else {
    for (i = 1; i < nn; ++i)
      if (fnode1[(n + i) % nn] != fnode2(i)) break;
    if (i == nn) return 1; // they match

    // Modify the permutation to reverse the orientation of fnode1.
    for (i = 1; i < nn; ++i)
      if (fnode1[(n - i + nn) % nn] != fnode2(i)) break;
    if (i == nn) return -1; // matched nodes but orientation is reversed
  }
  return 0; // different faces
}


template <class View_type>
bool
checkErrorList(const View_type& list, const std::string& msg, std::ostream& os)
{
  using Range_type = Kokkos::RangePolicy<LO, typename View_type::execution_space>;
  int count = 0;
  Kokkos::parallel_reduce(
    "MeshAudit::checkErrorList()",
    Range_type(0, list.extent(0)),
    KOKKOS_LAMBDA(const Entity_ID i, int& lcount) { lcount += list(i) > 0 ? 1 : 0; },
    count);

  if (count > 0) {
    if constexpr (Kokkos::SpaceAccessibility<typename View_type::execution_space,
                                             typename Kokkos::HostSpace>::accessible) {
      printErrorList(list, msg, count, os);
    } else {
      Kokkos::View<typename View_type::data_type, Amanzi::DefaultHost> list_host(
        list.label() + "_h", list.extent(0));
      Kokkos::deep_copy(list_host, list);
      printErrorList(list_host, msg, count, os);
    }
  }
  return count > 0;
}


template <class View_type>
void
printErrorList(const View_type& list, const std::string& msg, int count, std::ostream& os)
{
  static_assert(Kokkos::SpaceAccessibility<typename View_type::execution_space,
                                           typename Kokkos::HostSpace>::accessible);
  os << "ERROR: " << msg << ":";
  size_t lcount = 0;
  for (size_t i = 0; i != list.size(); ++i) {
    if (list[i] > 0) {
      os << " " << i;
      ++lcount;
      if (lcount >= MAX_OUT) break;
    }
  }
  if (lcount >= MAX_OUT) os << " [" << count - lcount << " items omitted]";
  os << std::endl;
}


template <class Mesh_type>
bool
MeshAudit_Base<Mesh_type>::globalAny_(bool value) const
{
  return Impl::globalAny(*this->comm_, value);
}


} // namespace Impl
} // namespace AmanziMesh
} // namespace Amanzi

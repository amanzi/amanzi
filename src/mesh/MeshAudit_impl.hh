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

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"

#include "Geometry.hh"
#include "CellTopology.hh"

#include "MeshAudit_decl.hh"

namespace Amanzi {
namespace AmanziMesh {
namespace Impl {

template<class Mesh_type>
MeshAudit_Base<Mesh_type>::MeshAudit_Base(
  const Teuchos::RCP<const Mesh_type>& mesh,
  std::ostream& os)
  : mesh_(mesh),
    comm_(mesh_->getComm()),
    os_(os),
    nnodes_all_(mesh_->getNumEntities(Entity_kind::NODE, Parallel_type::ALL)),
    nfaces_all_(mesh_->getNumEntities(Entity_kind::FACE, Parallel_type::ALL)),
    ncells_all_(mesh_->getNumEntities(Entity_kind::CELL, Parallel_type::ALL)),
    MAX_OUT(5)
{}


template<class MeshAudit_type>
void createTestDependencies_Base(MeshAudit_type& audit, AuditGraph<MeshAudit_type>& graph) {}


// This creates the dependency graph for the individual tests.  Adding a new
// test to the graph is a simple matter of adding a new stanza of the form
//
//   Graph::vertex_descriptor my_test_handle = addVertex_(g);
//   g_[my_test_handle].name = "description of my test";
//   g_[my_test_handle].test = &MeshAudit::my_test;
//   addEdge_(other_test_handle, my_test_handle);
//
// The last line specifies that other_test_handle is a pre-requisite for
// my_test_handle.  There may be multiple pre-requisites or none.
template<class MeshAudit_type>
void createTestDependencies_Geometry(MeshAudit_type& audit, AuditGraph<MeshAudit_type>& graph)
{
  createTestDependencies_Base(audit, graph);

  // Cell_to_nodes tests
  auto test02 = graph.addVertex();
  graph.getVertex(test02).name = "cell_to_nodes";
  graph.getVertex(test02).test = &MeshAudit_type::check_cell_to_nodes;

  auto test03 = graph.addVertex();
  graph.getVertex(test03).name = "node references by cells";
  graph.getVertex(test03).test = &MeshAudit_type::check_node_refs_by_cells;
  graph.addEdge(test02, test03);

  // Cell_to_faces tests
  auto test04 = graph.addVertex();
  graph.getVertex(test04).name = "cell_to_faces";
  graph.getVertex(test04).test = &MeshAudit_type::check_cell_to_faces;

  auto test05 = graph.addVertex();
  graph.getVertex(test05).name = "face references by cells";
  graph.getVertex(test05).test = &MeshAudit_type::check_face_refs_by_cells;
  graph.addEdge(test04, test05);

  // face_to_nodes tests
  auto test06 = graph.addVertex();
  graph.getVertex(test06).name = "face_to_nodes";
  graph.getVertex(test06).test = &MeshAudit_type::check_face_to_nodes;

  auto test08 = graph.addVertex();
  graph.getVertex(test08).name = "node references by faces";
  graph.getVertex(test08).test = &MeshAudit_type::check_node_refs_by_faces;
  graph.addEdge(test06, test08);

  // cell_to_face_dirs tests
  auto test09 = graph.addVertex();
  graph.getVertex(test09).name = "cell_to_face_dirs";
  graph.getVertex(test09).test = &MeshAudit_type::check_cell_to_face_dirs;

  // cell degeneracy test
  auto test10 = graph.addVertex();
  graph.getVertex(test10).name = "topological non-degeneracy of cells";
  graph.getVertex(test10).test = &MeshAudit_type::check_cell_degeneracy;
  graph.addEdge(test02, test10);

  // Consistency between the various mesh connectivity data.
  auto test11 = graph.addVertex();
  graph.getVertex(test11).name = "consistency of mesh connectivity data";
  graph.getVertex(test11).test = &MeshAudit_type::check_cell_to_faces_to_nodes;
  graph.addEdge(test02, test11);
  graph.addEdge(test04, test11);
  graph.addEdge(test06, test11);
  graph.addEdge(test09, test11);
  graph.addEdge(test10, test11);

  // node_to_coordinates tests
  auto test12 = graph.addVertex();
  graph.getVertex(test12).name = "node_to_coordinates";
  graph.getVertex(test12).test = &MeshAudit_type::check_node_to_coordinates;

  // cell_to_coordinates tests
  auto test13 = graph.addVertex();
  graph.getVertex(test13).name = "cell_to_coordinates";
  graph.getVertex(test13).test = &MeshAudit_type::check_cell_to_coordinates;
  graph.addEdge(test02, test13);
  graph.addEdge(test12, test13);

  // face_to_coordinates tests
  auto test14 = graph.addVertex();
  graph.getVertex(test14).name = "face_to_coordinates";
  graph.getVertex(test14).test = &MeshAudit_type::check_face_to_coordinates;
  graph.addEdge(test06, test14);
  graph.addEdge(test12, test14);

  // cell topology/geometry test
  auto test15 = graph.addVertex();
  graph.getVertex(test15).name = "cell geometry";
  graph.getVertex(test15).test = &MeshAudit_type::check_cell_geometry;
  graph.addEdge(test10, test15);
  graph.addEdge(test13, test15);

  // cell-to-face adjacencies and orientations
  auto test16 = graph.addVertex();
  graph.getVertex(test16).name = "consistency of face-cell adjacencies";
  graph.getVertex(test16).test = &MeshAudit_type::check_face_cell_adjacency_consistency;
  graph.addEdge(test04, test16);
  graph.addEdge(test09, test16);

  // face normals relto cell
  auto test17 = graph.addVertex();
  graph.getVertex(test17).name = "face normals are outward";
  graph.getVertex(test17).test = &MeshAudit_type::check_face_normal_relto_cell;
  graph.addEdge(test16, test17);

  // face normals
  auto test18 = graph.addVertex();
  graph.getVertex(test18).name = "face normals are oriented";
  graph.getVertex(test18).test = &MeshAudit_type::check_face_normal_orientation;
  graph.addEdge(test16, test18);

}


template<class MeshAudit_type>
void createTestDependencies_Maps(MeshAudit_type& audit, AuditGraph<MeshAudit_type>& graph)
{
  createTestDependencies_Geometry(audit, graph);

  // Entity_counts tests
  auto test01 = graph.addVertex();
  graph.getVertex(test01).name = "entity_counts";
  graph.getVertex(test01).test = &MeshAudit_type::check_entity_counts;

  // map tests
  auto test16 = graph.addVertex();
  graph.getVertex(test16).name = "owned and overlap node maps";
  graph.getVertex(test16).test = &MeshAudit_type::check_node_maps;

  auto test17 = graph.addVertex();
  graph.getVertex(test17).name = "owned and overlap face maps";
  graph.getVertex(test17).test = &MeshAudit_type::check_face_maps;

  auto test18 = graph.addVertex();
  graph.getVertex(test18).name = "owned and overlap cell maps";
  graph.getVertex(test18).test = &MeshAudit_type::check_cell_maps;

  // ghost data tests
  auto test19 = graph.addVertex();
  graph.getVertex(test19).name = "node_to_coordinates ghost data";
  graph.getVertex(test19).test = &MeshAudit_type::check_node_to_coordinates_ghost_data;
  graph.addEdge("node_to_coordinates", test19);
  graph.addEdge("owned and overlap node maps", test19);

  auto test20 = graph.addVertex();
  graph.getVertex(test20).name = "face_to_nodes ghost data";
  graph.getVertex(test20).test = &MeshAudit_type::check_face_to_nodes_ghost_data;
  graph.addEdge("face_to_nodes", test20);
  graph.addEdge("owned and overlap node maps", test20);
  graph.addEdge("owned and overlap face maps", test20);

  auto test21 = graph.addVertex();
  graph.getVertex(test21).name = "cell_to_nodes ghost data";
  graph.getVertex(test21).test = &MeshAudit_type::check_cell_to_nodes_ghost_data;
  graph.addEdge("cell_to_nodes", test21);
  graph.addEdge("owned and overlap cell maps", test21);
  graph.addEdge("owned and overlap node maps", test21);

  auto test22 = graph.addVertex();
  graph.getVertex(test22).name = "cell_to_faces ghost data";
  graph.getVertex(test22).test = &MeshAudit_type::check_cell_to_faces_ghost_data;
  graph.addEdge("cell_to_faces", test22);
  graph.addEdge("owned and overlap cell maps", test22);
  graph.addEdge("owned and overlap face maps", test22);

  // partition tests
  auto test32 = graph.addVertex();
  graph.getVertex(test32).name = "face partition";
  graph.getVertex(test32).test = &MeshAudit_type::check_face_partition;
  graph.addEdge("owned and overlap face maps", test32);
  graph.addEdge("owned and overlap cell maps", test32);
  graph.addEdge(test18, test32);

  auto test33 = graph.addVertex();
  graph.getVertex(test33).name = "node partition";
  graph.getVertex(test33).test = &MeshAudit_type::check_node_partition;
  graph.addEdge("owned and overlap node maps", test33);
  graph.addEdge("owned and overlap cell maps", test33);
};

template<class MeshAudit_type>
void createTestDependencies_Sets(MeshAudit_type& audit, AuditGraph<MeshAudit_type>& graph)
{
  createTestDependencies_Maps(audit, graph);

  // node set data tests
  auto test23 = graph.addVertex();
  graph.getVertex(test23).name = "node set IDs";
  graph.getVertex(test23).test = &MeshAudit_type::check_node_set_ids;

  auto test24 = graph.addVertex();
  graph.getVertex(test24).name = "node sets";
  graph.getVertex(test24).test = &MeshAudit_type::check_node_sets;
  graph.addEdge("owned and overlap node maps", test24);

  auto test25 = graph.addVertex();
  graph.getVertex(test25).name = "valid node set IDs";
  graph.getVertex(test25).test = &MeshAudit_type::check_valid_node_set_id;
  graph.addEdge("node set IDs", test25);

  // face set data tests
  auto test26 = graph.addVertex();
  graph.getVertex(test26).name = "face set IDs";
  graph.getVertex(test26).test = &MeshAudit_type::check_face_set_ids;

  auto test27 = graph.addVertex();
  graph.getVertex(test27).name = "face sets";
  graph.getVertex(test27).test = &MeshAudit_type::check_face_sets;
  graph.addEdge("owned and overlap face maps", test27);
  graph.addEdge("face set IDs", test27);

  auto test28 = graph.addVertex();
  graph.getVertex(test28).name = "valid face set IDs";
  graph.getVertex(test28).test = &MeshAudit_type::check_valid_face_set_id;
  graph.addEdge("face set IDs", test28);

  // cell set data tests
  auto test29 = graph.addVertex();
  graph.getVertex(test29).name = "cell set IDs";
  graph.getVertex(test29).test = &MeshAudit_type::check_cell_set_ids;

  auto test30 = graph.addVertex();
  graph.getVertex(test30).name = "cell sets";
  graph.getVertex(test30).test = &MeshAudit_type::check_cell_sets;
  graph.addEdge("owned and overlap cell maps", test30);
  graph.addEdge("cell set IDs", test30);

  auto test31 = graph.addVertex();
  graph.getVertex(test31).name = "valid cell set IDs";
  graph.getVertex(test31).test = &MeshAudit_type::check_valid_cell_set_id;
  graph.addEdge("cell set IDs", test31);
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
// elements in the corresponding Epetra_Maps.  This applies to nodes, faces
// and cells, including ghosts and not.  Here the maps are considered to be
// the authoritative source of information.  A positive value is returned
// if any discrepancy is found, but it is safe to perform other tests as
// they do not use these methods.
template<class Mesh_type>
bool MeshAudit_Geometry<Mesh_type>::check_entity_counts() const
{
  int n, nref;
  bool error = false;
  // Check the number of owned nodes.
  n = mesh_->getNumEntities(Entity_kind::NODE,Parallel_type::OWNED);
  nref = mesh_->getMap(AmanziMesh::Entity_kind::NODE, false).NumMyElements();
  if (n != nref) {
    os_ << ": ERROR: getNumEntities(NODE,OWNED)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  // Check the number of used nodes.
  n = mesh_->getNumEntities(Entity_kind::NODE,Parallel_type::ALL);
  nref = mesh_->getMap(AmanziMesh::Entity_kind::NODE, true).NumMyElements();
  if (n != nref) {
    os_ << "ERROR: getNumEntities(NODE,ALL)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  // Check the number of owned faces.
  n = mesh_->getNumEntities(Entity_kind::FACE,Parallel_type::OWNED);
  nref = mesh_->getMap(AmanziMesh::Entity_kind::FACE, false).NumMyElements();
  if (n != nref) {
    os_ << "ERROR: getNumEntities(FACE,OWNED)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  // Check the number of used faces.
  n = mesh_->getNumEntities(Entity_kind::FACE,Parallel_type::ALL);
  nref = mesh_->getMap(AmanziMesh::Entity_kind::FACE, true).NumMyElements();
  if (n != nref) {
    os_ << "ERROR: getNumEntities(FACE,ALL)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  // Check the number of owned cells.
  n = mesh_->getNumEntities(Entity_kind::CELL,Parallel_type::OWNED);
  nref = mesh_->getMap(AmanziMesh::Entity_kind::CELL, false).NumMyElements();
  if (n != nref) {
    os_ << "ERROR: getNumEntities(CELL,OWNED)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  // Check the number of used cells.
  n = mesh_->getNumEntities(Entity_kind::CELL,Parallel_type::ALL);
  nref = mesh_->getMap(AmanziMesh::Entity_kind::CELL, true).NumMyElements();
  if (n != nref) {
    os_ << "ERROR: getNumEntities(CELL,ALL)=" << n << "; should be " << nref << std::endl;
    error = true;
  }
  return globalAny_(error);
}


// Check that cell_to_nodes successfully returns valid references to local
// nodes.  Here the std::vector accessor is used.  The consistency of
// the alternative accessors is checked elsewhere.  A nonzero return value
// signals an error, and further tests using its data should be avoided.
template<class Mesh_type>
bool MeshAudit_Geometry<Mesh_type>::check_cell_to_nodes() const
{
  std::vector<Entity_ID> bad_cells, bad_cells1;
  Entity_ID_List free_nodes;
  Entity_ID_List cnode;

  for (Entity_ID j = 0; j < ncells_all_; ++j) {
    try {
      mesh_->getCellNodes(j, cnode); // this may fail
      bool invalid_refs = false;
      for (int k = 0; k < cnode.size(); ++k) {
        if (cnode[k] >= nnodes_all_) invalid_refs = true;
      }
      if (invalid_refs) bad_cells.push_back(j);
    } catch (...) {
      bad_cells1.push_back(j);
    }
  }

  bool error = false;

  if (!bad_cells.empty()) {
    os_ << "ERROR: invalid nodes referenced by cells:";
    writeList_(bad_cells);
    error = true;
  }

  if (!bad_cells1.empty()) {
    os_ << "ERROR: caught exception for cells:";
    writeList_(bad_cells1);
    error = true;
  }
  return globalAny_(error);
}

// Check that every node is referenced by at least one cell.  This assumes
// that cell_to_nodes have been verified to return valid data.  A nonzero
// return value indicates that one or more nodes are not attached to any cell.
template<class Mesh_type>
bool MeshAudit_Geometry<Mesh_type>::check_node_refs_by_cells() const
{
  Entity_ID_List cnode;
  std::vector<bool> ref(nnodes_all_, false);

  for (Entity_ID j = 0; j < ncells_all_; ++j) {
    mesh_->getCellNodes(j, cnode); // this should not fail
    for (int k = 0; k < cnode.size(); ++k) ref[cnode[k]] = true;
  }

  std::vector<Entity_ID> free_nodes;
  for (int j = 0; j < nnodes_all_; ++j) {
    if (!ref[j]) free_nodes.push_back(j);
  }

  bool error = false;

  if (!free_nodes.empty()) {
    os_ << "ERROR: found unreferenced nodes:";
    writeList_(free_nodes);
    error = true;
  }

  return globalAny_(error);
}

// Check that cell_to_faces successfully returns valid references to local
// faces.  Here the std::vector accessor is used.  The consistency of
// the alternative accessors is checked elsewhere.  A nonzero return value
// signals an error, and further tests using its data should be avoided.
template<class Mesh_type>
bool MeshAudit_Geometry<Mesh_type>::check_cell_to_faces() const
{
  std::vector<Entity_ID> bad_cells, bad_cells1;
  Entity_ID_List cface;

  for (Entity_ID j = 0; j < ncells_all_; ++j) {
    try {
      mesh_->getCellFaces(j, cface); // this may fail
      bool invalid_refs = false;
      for (int k = 0; k < cface.size(); ++k) {
        if (cface[k] >= nfaces_all_) invalid_refs = true;
      }
      if (invalid_refs) bad_cells.push_back(j);
    } catch (...) {
      bad_cells1.push_back(j);
    }
  }

  bool error = false;

  if (!bad_cells.empty()) {
    os_ << "ERROR: invalid faces referenced by cells:";
    writeList_(bad_cells);
    error = true;
  }

  if (!bad_cells1.empty()) {
    os_ << "ERROR: caught exception for cells:";
    writeList_(bad_cells1);
    error = true;
  }

  return globalAny_(error);
}

// Check that every face is referenced by exactly one or two cells.  This
// assumes that cell_to_faces have been verified to return valid data.  A
// nonzero return value indicates that faces were found that do not belong
// to any cell or belong to more than two cells (bad topology).
template<class Mesh_type>
bool MeshAudit_Geometry<Mesh_type>::check_face_refs_by_cells() const
{
  Entity_ID_List cface;
  std::vector<unsigned int> refs(nfaces_all_, 0);

  for (Entity_ID j = 0; j < ncells_all_; ++j) {
    mesh_->getCellFaces(j, cface);
    for (int k = 0; k < cface.size(); ++k) (refs[cface[k]])++;
  }

  std::vector<Entity_ID> free_faces;
  std::vector<Entity_ID> bad_faces;

  for (int j = 0; j < nfaces_all_; ++j) {
    if (refs[j] == 0)
      free_faces.push_back(j);
    else if (refs[j] > 2)
      bad_faces.push_back(j);
  }

  bool error = false;

  if (!free_faces.empty()) {
    os_ << "ERROR: found unreferenced faces:";
    writeList_(free_faces);
    error = true;
  }

  // note that 3D space meshes of 2D manifolds may allow for more than 2 cells on a face
  if (!bad_faces.empty() && mesh_->getSpaceDimension() == mesh_->getManifoldDimension()) {
    os_ << "ERROR: found faces shared by more than two cells:";
    writeList_(bad_faces);
    error = true;
  }

  return globalAny_(error);
}

// Check that face_to_nodes successfully returns valid references to local
// nodes.  Here the std::vector accessor is used.  The consistency of
// the alternative accessors is checked elsewhere.  A nonzero return value
// signals an error, and further tests using its data should be avoided.
template<class Mesh_type>
bool MeshAudit_Geometry<Mesh_type>::check_face_to_nodes() const
{
  std::vector<Entity_ID> bad_faces, bad_faces1;
  Entity_ID_List fnode;

  for (Entity_ID j = 0; j < nfaces_all_; ++j) {
    try {
      mesh_->getFaceNodes(j, fnode); // this may fail
      bool invalid_refs = false;
      for (int k = 0; k < fnode.size(); ++k) {
        if (fnode[k] >= nnodes_all_) invalid_refs = true;
      }
      if (invalid_refs) bad_faces.push_back(j);
    } catch (...) {
      bad_faces1.push_back(j);
    }
  }

  bool error = false;

  if (!bad_faces.empty()) {
    os_ << "ERROR: invalid nodes referenced by faces:";
    writeList_(bad_faces);
    error = true;
  }

  if (!bad_faces1.empty()) {
    os_ << "ERROR: caught exception for faces:";
    writeList_(bad_faces1);
    error = true;
  }
  return globalAny_(error);
}

// Check that every node is referenced by at least one face.  This assumes
// that face_to_nodes has been verified to return valid data.  A nonzero
// return value indicates that one or more nodes are not attached to any face.
template<class Mesh_type>
bool MeshAudit_Geometry<Mesh_type>::check_node_refs_by_faces() const
{
  Entity_ID_List fnode;
  std::vector<bool> ref(nnodes_all_, false);

  for (Entity_ID j = 0; j < nfaces_all_; ++j) {
    mesh_->getFaceNodes(j, fnode);
    for (int k = 0; k < fnode.size(); ++k) ref[fnode[k]] = true;
  }

  std::vector<Entity_ID> free_nodes;
  for (int j = 0; j < nnodes_all_; ++j) {
    if (!ref[j]) free_nodes.push_back(j);
  }

  bool error = false;

  if (!free_nodes.empty()) {
    os_ << "ERROR: found unreferenced nodes:";
    writeList_(free_nodes);
    error = true;
  }

  return globalAny_(error);
}

// Check that cell_to_face_dirs successfully returns data for all cells and
// that the values are either +1 or -1. The std::vector-based method is used;
// the consistency of the alternative methods is checked elsewhere.  If this
// test fails, further tests using this data should be avoided.
template<class Mesh_type>
bool MeshAudit_Geometry<Mesh_type>::check_cell_to_face_dirs() const
{
  Entity_ID_List faces;
  Entity_Direction_List fdirs;
  std::vector<Entity_ID> bad_cells, bad_cells_exc;

  for (Entity_ID j = 0; j < ncells_all_; ++j) {
    //Kokkos::resize(fdirs, 6); 
    //initView(fdirs, std::numeric_limits<int>::max());
    try {
      mesh_->getCellFacesAndDirs(j, faces, &fdirs);  // this may fail
      bool bad_data = false;
      for (int k = 0; k < fdirs.size(); ++k)
        if (fdirs[k] != -1 && fdirs[k] != 1) bad_data = true;
      if (bad_data) bad_cells.push_back(j);
    } catch (...) {
      bad_cells_exc.push_back(j);
    }
  }

  bool error = false;

  if (!bad_cells.empty()) {
    os_ << "ERROR: inadmissable or no data for cells:";
    writeList_(bad_cells);
    error = true;
  }

  if (!bad_cells_exc.empty()) {
    os_ << "ERROR: caught exception for cells:";
    writeList_(bad_cells_exc);
    error = true;
  }

  return globalAny_(error);
}


// Check that cells are not topologically degenerate (repeated node index).
// If this test fails, many further tests involving the cell_to_nodes data
// should be avoided.
//
template<class Mesh_type>
bool MeshAudit_Geometry<Mesh_type>::check_cell_degeneracy() const
{
  os_ << "Checking cells for topological degeneracy ..." << std::endl;

  Entity_ID_List cnode;
  std::vector<Entity_ID> bad_cells;

  for (Entity_ID j = 0; j < ncells_all_; ++j) {
    mesh_->getCellNodes(j, cnode); // should not fail
    if (!this->areDistinctValues_(cnode)) bad_cells.push_back(j);
  }

  bool error = false;

  if (!bad_cells.empty()) {
    os_ << "ERROR: found topologically degenerate cells:";
    writeList_(bad_cells);
    error = true;
  }

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
template<class Mesh_type>
bool MeshAudit_Geometry<Mesh_type>::check_cell_to_faces_to_nodes() const
{
  Entity_ID_List cnode;
  Entity_ID_List cface;
  Entity_ID_List fnode_ref;
  Entity_ID_List fnode;
  Entity_Direction_List fdirs;
  std::vector<Entity_ID> bad_cells0;
  std::vector<Entity_ID> bad_cells1;

  // unordered meshes cannot pass this test.
  if (!mesh_->isOrdered()) {
    os_ << "  Skipped (unordered mesh)" << std::endl;
    return false;
  }

  for (Entity_ID j = 0; j < ncells_all_; ++j) {
    Cell_type ctype_enum = mesh_->getCellType(j);

    // If this is a general, non-standard element there is nothing to
    // to check against
    if (ctype_enum == Cell_type::UNKNOWN ||
        ctype_enum == Cell_type::POLYGON ||
        ctype_enum == Cell_type::POLYHED)
      continue;

    int ctype = (int) ctype_enum;
    mesh_->getCellNodes(j, cnode); // this should not fail
    mesh_->getCellFacesAndDirs(j, cface, &fdirs); // this should not fail

    bool bad_face = false;
    bool bad_dir  = false;

    if (cface.size() != Topology::nface_std[ctype])
      bad_face = true;
    else {

      for (int k = 0; k < cface.size(); ++k) {

        mesh_->getFaceNodes(cface[k], fnode); // this should not fail

        int nfn = Topology::nfnodes_std[ctype][k];

        if (fnode.size() != nfn) {
          bad_face = true;
          break;
        }

        Kokkos::resize(fnode_ref, nfn); 
        for (int i = 0; i < nfn; ++i) {
          int nodenum = Topology::fnodes_std[ctype][k][i];
          fnode_ref[i] = cnode[nodenum];
        }

        int dir = this->isSameFace_(fnode, fnode_ref); // should be the same face
        if (dir == 0) { // wrong face
          bad_face = true;
          break;
        }
        else if (dir != fdirs[k]) { // right face but wrong dir value
          bad_dir = true;
          break;
        }
      }

    }
    if (bad_face) bad_cells0.push_back(j);
    if (bad_dir)  bad_cells1.push_back(j);
  }

  bool error = false;

  if (!bad_cells0.empty()) {
    os_ << "ERROR: bad getCellFaces values for cells:";
    writeList_(bad_cells0);
    error = true;
  }

  // algorithms on a non-manifold use multiple normals and special continuity equations
  // for fluxes, so that orientation does not play role. This may change in the future.
  if (!bad_cells1.empty() && mesh_->getSpaceDimension() == mesh_->getManifoldDimension()) {
    os_ << "ERROR: bad cell_get_face_dirs values for cells:";
    writeList_(bad_cells1);
    error = true;
  }

  return globalAny_(error);
}

// Check that getNodeCoordinate successfully returns data for all nodes
// A negative return value signals a terminal error
//
template<class Mesh_type>
bool MeshAudit_Geometry<Mesh_type>::check_node_to_coordinates() const
{
  std::vector<double> x(3);
  std::vector<Entity_ID> bad_nodes, bad_nodes_exc;
  int spdim = mesh_->getSpaceDimension();
  std::vector<Entity_ID> bad_nodes0;

  for (Entity_ID j = 0; j < nnodes_all_; ++j) {
    AmanziGeometry::Point x0(spdim);
    x0.set(std::numeric_limits<double>::max());
    try {
      x0 = mesh_->getNodeCoordinate(j); // this may fail
      bool bad_data = false;
      for (int k = 0; k < spdim; ++k)
        if (x0[k] == std::numeric_limits<double>::max()) bad_data = true;
      if (bad_data) bad_nodes.push_back(j);
    } catch (...) {
      bad_nodes_exc.push_back(j);
    }
  }

  bool error = false;

  if (!bad_nodes.empty()) {
    os_ << "ERROR: missing coordinate data for nodes:";
    writeList_(bad_nodes);
    error = true;
  }

  if (!bad_nodes_exc.empty()) {
    os_ << "ERROR: caught exception for nodes:";
    writeList_(bad_nodes_exc);
    error = true;
  }

  return globalAny_(error);
}

// Check that cell_get_coordinates successfully returns data for all
// cells, and that the data is identical to that returned by composing
// with cell_to_nodes.  The consistency of the alternative accessors is checked
// elsewhere.  If this test fails, further tests using cell_to_coordinates data
// should be avoided.
//
template<class Mesh_type>
bool MeshAudit_Geometry<Mesh_type>::check_cell_to_coordinates() const
{
  int spdim = mesh_->getSpaceDimension();
  Entity_ID_List cnode;
  std::vector<Entity_ID> bad_cells, bad_cells_exc;

  for (Entity_ID j = 0; j < ncells_all_; ++j) {
    try {
      auto x0 = mesh_->getCellCoordinates(j); // this may fail
      mesh_->getCellNodes(j, cnode); // this should not fail
      for (int k = 0; k < cnode.size(); ++k){
        bool bad_data = false;
        auto xref = mesh_->getNodeCoordinate(cnode[k]); // this should not fail
        for (int i = 0; i < spdim; i++)
          if (x0[k][i] != xref[i]) {
            bad_data = true;
            bad_cells.push_back(j);
            break;
          }
        if (bad_data)
           break;
      }
    } catch (...) {
      bad_cells_exc.push_back(j);
    }
  }

  bool error = false;

  if (!bad_cells.empty()) {
    os_ << "ERROR: bad cell_to_coordinates data for cells:";
    writeList_(bad_cells);
    error = true;
  }

  if (!bad_cells_exc.empty()) {
    os_ << "ERROR: caught exception for cells:";
    writeList_(bad_cells_exc);
    error = true;
  }

  return globalAny_(error);
}


// Check that face_get_coordinates successfully returns data for all
// faces, and that this data is identical to that returned by
// composing getNodeCoordinate with getFaceNodes.  This assumes
// that the std::vector form of getFaceNodes and pointer-based form
// of getNodeCoordinate have been verified to return valid data. A
// negative return value signals a terminal error.
//
template<class Mesh_type>
bool MeshAudit_Geometry<Mesh_type>::check_face_to_coordinates() const
{
  int spdim = mesh_->getSpaceDimension();
  Entity_ID_List fnode;
  std::vector<Entity_ID> bad_faces, bad_faces_exc;

  for (Entity_ID j = 0; j < nfaces_all_; ++j) {
    try {
      bool bad_data = false;
      auto x0 = mesh_->getFaceCoordinates(j); // this may fail
      mesh_->getFaceNodes(j, fnode); // this should not fail
      for (int k = 0; k < fnode.size(); ++k) {
        AmanziGeometry::Point xref(spdim);
        xref = mesh_->getNodeCoordinate(fnode[k]); // this should not fail
        for (int i = 0; i < spdim; i++)
          if (x0[k][i] != xref[i]) {
            bad_data = true;
            bad_faces.push_back(j);
            break;
          }
        if (bad_data)
          break;
      }
    } catch (...) {
      bad_faces_exc.push_back(j);
    }
  }

  bool error = false;

  if (!bad_faces.empty()) {
    os_ << "ERROR: bad face_to_coordinates data for faces:";
    writeList_(bad_faces);
    error = true;
  }

  if (!bad_faces_exc.empty()) {
    os_ << "ERROR: caught exception for faces:";
    writeList_(bad_faces_exc);
    error = true;
  }

  return globalAny_(error);
}

//
// Check that upward and downward adjacencies and orientations are consistent
// between cells and faces.
//
template<class Mesh_type>
bool MeshAudit_Geometry<Mesh_type>::check_face_cell_adjacency_consistency() const
{
  std::vector<Entity_ID> bad_faces1, bad_faces2, bad_faces_exc;

  for (Entity_ID f=0; f!=nfaces_all_; ++f) {
    try {
      bool bad_data = false;

      Entity_ID_List fcells;
      mesh_->getFaceCells(f, Parallel_type::ALL, fcells);
      for (const auto& c : fcells) {
        Entity_ID_List cfaces;
        Entity_Direction_List cfdirs;
        mesh_->getCellFacesAndDirs(c, cfaces, &cfdirs);

        // find the match
        int i = 0;
        for (; i != cfaces.size(); ++i)
          if (f == cfaces[i]) break;
        if (i == cfaces.size()) {
          bad_faces1.push_back(f);
          bad_data = true;
          break;
        }

        // check directions
        int orientation = 0;
        auto fnormal = mesh_->getFaceNormal(f, c, &orientation);
        if (orientation == 0) {
          if (mesh_->getSpaceDimension() == mesh_->getManifoldDimension()) {
            // should be an orientation but there isn't
            bad_faces2.push_back(f);
            bad_data = true;
            break;
          }
        } else {
          if (orientation != cfdirs[i]) {
            // orientation is there, not equivalent to the dir from the cell
            bad_faces2.push_back(f);
            bad_data = true;
            break;
          }
        }
      }
      if (bad_data) break;
    } catch (...) {
      bad_faces_exc.push_back(f);
    }
  }

  bool error = false;
  if (!bad_faces1.empty()) {
    os_ << "ERROR: , inconsistent adjacencies:";
    writeList_(bad_faces1);
    error = true;
  }
  if (!bad_faces2.empty()) {
    os_ << "ERROR: , inconsistent orientations:";
    writeList_(bad_faces2);
    error = true;
  }
  if (!bad_faces_exc.empty()) {
    os_ << "ERROR: caught exception:";
    writeList_(bad_faces_exc);
    error = true;
  }
  return globalAny_(error);
}


//
// Check that a face normal, relative to a cell, is outward.  Note this
// requires star-convex meshes.
//
template<class Mesh_type>
bool MeshAudit_Geometry<Mesh_type>::check_face_normal_relto_cell() const
{
  std::vector<Entity_ID> bad_faces, bad_faces_exc;

  for (Entity_ID j = 0; j < nfaces_all_; ++j) {
    try {
      auto fc = mesh_->getFaceCentroid(j);

      Entity_ID_List fcells;
      mesh_->getFaceCells(j, Parallel_type::ALL, fcells);
      for (int k = 0; k < fcells.size(); ++k) {
        AmanziGeometry::Point cc = mesh_->getCellCentroid(fcells[k]);
        AmanziGeometry::Point fnormal = mesh_->getFaceNormal(j, fcells[k]);
        if ((fc - cc) * fnormal < 0) {
          bad_faces.push_back(j);
          break;
        }
      }
    } catch (...) {
      bad_faces_exc.push_back(j);
    }
  }

  bool error = false;
  if (!bad_faces.empty()) {
    os_ << "ERROR: bad face_normal, not outward:";
    writeList_(bad_faces);
    error = true;
  }
  if (!bad_faces_exc.empty()) {
    os_ << "ERROR: caught exception for getFaceNormal:";
    writeList_(bad_faces_exc);
    error = true;
  }
  return globalAny_(error);
}


//
// Check that a face normal has the correct orientation.
//
template<class Mesh_type>
bool MeshAudit_Geometry<Mesh_type>::check_face_normal_orientation() const
{
  // this test is only for meshes which have a natural normal
  if (mesh_->getSpaceDimension() != mesh_->getManifoldDimension()) return false;
  std::vector<Entity_ID> bad_faces, bad_faces_exc;

  for (Entity_ID j = 0; j < nfaces_all_; ++j) {
    try {
      AmanziGeometry::Point fnormal = mesh_->getFaceNormal(j);

      Entity_ID_List fcells;
      mesh_->getFaceCells(j, Parallel_type::ALL, fcells);
      for (int k = 0; k < fcells.size(); ++k) {
        int orientation = 0;
        AmanziGeometry::Point fnormalc = mesh_->getFaceNormal(j, fcells[k], &orientation);
        if (orientation != 0) {
          if (AmanziGeometry::norm((orientation * fnormal) - fnormalc) > 1.e-16) {
            bad_faces.push_back(j);
            break;
          }
        }
      }
    } catch (...) {
      bad_faces_exc.push_back(j);
    }
  }

  bool error = false;
  if (!bad_faces.empty()) {
    os_ << "ERROR: bad face_normal, inconsistent orientation:";
    writeList_(bad_faces);
    error = true;
  }
  if (!bad_faces_exc.empty()) {
    os_ << "ERROR: caught exception for getFaceNormal with orientation:";
    writeList_(bad_faces_exc);
    error = true;
  }
  return globalAny_(error);
}



// The cells must not be degenerate, either topologically (repeated
// node index) or geometrically (with coincident nodes).
//
template<class Mesh_type>
bool MeshAudit_Geometry<Mesh_type>::check_cell_geometry() const
{
  os_ << "Checking cell geometry ..." << std::endl;
  std::vector<Entity_ID> bad_cells;

  for (Entity_ID c=0; c!=ncells_all_; ++c) {
    double hvol = mesh_->getCellVolume(c);
    if (hvol <= 0.0) bad_cells.push_back(c);
  }

  bool error = false;
  if (!bad_cells.empty()) {
    os_ << "ERROR: found cells with non-positive volumes:";
    writeList_(bad_cells);
    os_ << "       Cells either geometrically degenerate or have bad topology.";
    error = true;
  }
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
template<class Mesh_type>
bool MeshAudit_Maps<Mesh_type>::check_node_maps() const
{
  return check_maps(mesh_->getMap(AmanziMesh::Entity_kind::NODE, false),
                    mesh_->getMap(AmanziMesh::Entity_kind::NODE, true));
}

template<class Mesh_type>
bool MeshAudit_Maps<Mesh_type>::check_face_maps() const
{
  return check_maps(mesh_->getMap(AmanziMesh::Entity_kind::FACE, false),
                    mesh_->getMap(AmanziMesh::Entity_kind::FACE, true));
}

template<class Mesh_type>
bool MeshAudit_Maps<Mesh_type>::check_cell_maps() const
{
  return check_maps(mesh_->getMap(AmanziMesh::Entity_kind::CELL, false),
                    mesh_->getMap(AmanziMesh::Entity_kind::CELL, true));
}

template<class Mesh_type>
bool MeshAudit_Maps<Mesh_type>::check_maps(const Epetra_Map &map_own, const Epetra_Map &map_use) const
{
  bool error = false;

  // Local index should start at 0.
  if (map_own.IndexBase() != 0) {
    os_ << "ERROR: the owned map's index base is not 0." << std::endl;
    error = true;
  }
  if (map_use.IndexBase() != 0) {
    os_ << "ERROR: the overlap map's index base is not 0." << std::endl;
    error = true;
  }

  // Check that the owned map is 1-1.
  if (!map_own.UniqueGIDs()) {
    os_ << "ERROR: owned map is not 1-to-1" << std::endl;
    error = true;
  }

  // Check that the global ID space is contiguously divided (but not
  // necessarily uniformly) across all processers

  std::vector<int> owned_GIDs(map_own.NumMyElements());
  for (int i = 0; i < map_own.NumMyElements(); i++)
    owned_GIDs[i] = map_own.GID(i);
  std::sort(owned_GIDs.begin(), owned_GIDs.end());

  for (int i = 0; i < map_own.NumMyElements()-1; i++) {
    int diff = owned_GIDs[i+1]-owned_GIDs[i];
    if (diff > 1) {
      os_ << "ERROR: owned map is not contiguous" << std::endl;
      os_ << "Global IDs jump from " << owned_GIDs[i] << " to " <<
        owned_GIDs[i+1] << std::endl;
      error = true;
    }
  }

  error = globalAny_(error);
  if (error) return error;

  if (this->comm_->NumProc() == 1) {
    // Serial or 1-process MPI
    if (!map_use.SameAs(map_own)) {
      os_ << "ERROR: the overlap map differs from the owned map (single process)." << std::endl;
      error = true;
    }
    return globalAny_(error);

  } else {
    // Multi-process MPI
    int num_own = map_own.NumMyElements();
    int num_use = map_use.NumMyElements();
    int num_ovl = num_use - num_own;

    // Verify that the used map extends the owned map.
    bool bad_map = false;
    if (num_ovl < 0)
      bad_map = true;
    else {
      for (int j = 0; j < num_own; ++j)
        if (map_use.GID(j) != map_own.GID(j)) bad_map = true;
    }
    if (bad_map) {
      os_ << "ERROR: overlap map does not extend the owned map." << std::endl;
      error = true;
    }

    error = globalAny_(error);
    if (error) return error;

    // Verify that the overlap indices are owned by other processes.
    int *gids = new int[num_ovl];
    int *pids = new int[num_ovl];
    int *lids = new int[num_ovl];
    for (int j = 0; j < num_ovl; ++j) gids[j] = map_use.GID(j+num_own);
    map_own.RemoteIDList(num_ovl, gids, pids, lids);
    bad_map = false;
    for (int j = 0; j < num_ovl; ++j)
      if (pids[j] < 0 || pids[j] == this->comm_->MyPID()) bad_map = true;
    if (bad_map) {
      os_ << "ERROR: invalid ghosts in overlap map." << std::endl;
      error = true;
    }

    // Look for duplicates among the overlap indices.
    std::vector<int> ovl_gids(gids, gids+num_ovl);
    sort(ovl_gids.begin(), ovl_gids.end());
    if (adjacent_find(ovl_gids.begin(),ovl_gids.end()) != ovl_gids.end()) {
      os_ << "ERROR: duplicate ghosts in overlap map." << std::endl;
      error = true;
    }

    delete [] lids;
    delete [] pids;
    delete [] gids;
    return globalAny_(error);
  }
}

// Check that ghost nodes are exact copies of their master.
// This simply means that they have the same coordinates.
//
template<class Mesh_type>
bool MeshAudit_Maps<Mesh_type>::check_node_to_coordinates_ghost_data() const
{
  int spdim = mesh_->getSpaceDimension();

  const Epetra_Map &node_map_own = mesh_->getMap(AmanziMesh::Entity_kind::NODE, false);
  const Epetra_Map &node_map_use = mesh_->getMap(AmanziMesh::Entity_kind::NODE, true);

  int nnode_own = node_map_own.NumMyElements();
  int nnode_use = node_map_use.NumMyElements();

  AmanziGeometry::Point coord(spdim);
  std::vector<Entity_ID> bad_nodes;

  Epetra_MultiVector coord_use(node_map_use,spdim);
  double **data;
  coord_use.ExtractView(&data);
  Epetra_MultiVector coord_own(View, node_map_own, data, spdim);

  for (Entity_ID j = 0; j < nnode_own; ++j) {
    coord = mesh_->getNodeCoordinate(j);
    for (int k = 0; k < spdim; ++k) coord_own[k][j] = coord[k];
  }

  Epetra_Import importer(node_map_use, node_map_own);
  coord_use.Import(coord_own, importer, Insert);

  for (Entity_ID j = nnode_own; j < nnode_use; ++j) {
    coord = mesh_->getNodeCoordinate(j);
    bool bad_data = false;
    for (int k = 0; k < spdim; ++k)
      if (coord[k] != coord_use[k][j]) bad_data = true;
    if (bad_data) bad_nodes.push_back(j);
  }

  bool error = false;

  if (!bad_nodes.empty()) {
    os_ << "ERROR: found ghost nodes with incorrect coordinates:";
    writeList_(bad_nodes);
    error = true;
  }

  return globalAny_(error);
}

// Check that ghost faces are exact copies of their master.  This means
// that the the GIDs of the nodes defining the face are the same, including
// their order (face orientation).
//
template<class Mesh_type>
bool MeshAudit_Maps<Mesh_type>::check_face_to_nodes_ghost_data() const
{
  const Epetra_Map &node_map = mesh_->getMap(AmanziMesh::Entity_kind::NODE, true);
  const Epetra_Map &face_map_own = mesh_->getMap(AmanziMesh::Entity_kind::FACE, false);
  const Epetra_Map &face_map_use = mesh_->getMap(AmanziMesh::Entity_kind::FACE, true);

  int nface_own = face_map_own.NumMyElements();
  int nface_use = face_map_use.NumMyElements();

  Entity_ID_List fnode;
  std::vector<Entity_ID> bad_faces, bad_faces1, bad_faces2;

  // Create a matrix of the GIDs for all owned faces.

  int maxnodes = 0;
  for (Entity_ID j = 0; j < nface_own; ++j) {
    mesh_->getFaceNodes(j, fnode);
    maxnodes = fnode.size() > maxnodes ? fnode.size() : maxnodes;
  }

  // parallel calls require uniformity
  int maxnodes_tmp(maxnodes);
  this->comm_->MaxAll(&maxnodes_tmp, &maxnodes, 1);

  Epetra_IntSerialDenseMatrix gids(nface_use,maxnodes); // no Epetra_IntMultiVector :(
  for (Entity_ID j = 0; j < nface_own; ++j) {
    mesh_->getFaceNodes(j, fnode);
    for (int k = 0; k < fnode.size(); ++k)
      gids(j,k) = node_map.GID(fnode[k]);
    for (int k = fnode.size(); k < maxnodes; ++k)
      gids(j, k) = 0;
  }

  // Import these GIDs to all used faces; sets values on ghost faces.
  Epetra_Import importer(face_map_use, face_map_own);
  for (int k = 0; k < maxnodes; ++k) {
    Epetra_IntVector kgids_own(View, face_map_own, gids[k]);
    Epetra_IntVector kgids_use(View, face_map_use, gids[k]);
    kgids_use.Import(kgids_own, importer, Insert);
  }

  // Compare the ghost face GIDs against the reference values just computed.
  for (Entity_ID j = nface_own; j < nface_use; ++j) {
    mesh_->getFaceNodes(j, fnode);
    bool bad_data = false;
    for (int k = 0; k < fnode.size(); ++k)
      if (node_map.GID(fnode[k]) != gids(j,k)) {
        bad_data = true;
      }
    if (bad_data) {
      // Determine just how bad the data is.
      Entity_ID_List fnode_ref("fnode_ref", maxnodes);
      for (int k = 0; k < maxnodes; ++k) {
        fnode[k] = node_map.GID(fnode[k]);
        fnode_ref[k] = gids(j,k);
      }
      int n = this->isSameFace_(fnode, fnode_ref);
      switch (n) {
        case 0: // completely bad -- different face
          bad_faces.push_back(j);

          os_ << "P " << this->comm_->MyPID() << ": ghost face " << j << " (GID "
                    << face_map_use.GID(j) << ")," <<
            " has different nodes than its master " << std::endl;
          os_ << "ghost face nodes (GIDs): ";
          for (int k = 0; k < fnode.size(); ++k)
            os_ << fnode[k] << " ";
          os_ << std::endl;
          os_ << "master face nodes (GIDs): ";
          for (int k = 0; k < fnode_ref.size(); ++k)
            os_ << fnode_ref[k] << " ";
          os_ << std::endl;
          break;
        case -1: // very bad -- same face but wrong orientation
          bad_faces1.push_back(j);

          os_ << "P " << this->comm_->MyPID() << ": ghost face " << j << ", has different orientation than its master " << std::endl;
          os_ << "ghost face nodes (GIDs): ";
          for (int k = 0; k < fnode.size(); ++k)
            os_ << fnode[k];
          os_ << std::endl;
          os_ << "master face nodes (GIDs): ";
          for (int k = 0; k < fnode_ref.size(); ++k)
            os_ << fnode_ref[k];
          os_ << std::endl;
          break;
        case 1:
          // This is fine because there is no way to ensure this for
          // general meshes (think of building a tet mesh and defining
          // the faces such that they return the same vertices in the
          // same order (unpermuted) no matter which way the elements
          // match up. You can only ensure this for hexahedral meshes
          // As long as the faces have the same vertices and have the
          // same direction, but the vertices are permuted, its fine
          break;
      }
    }
  }

  bool error = false;

  if (!bad_faces.empty()) {
    os_ << "ERROR: found bad data for ghost faces:";
    writeList_(bad_faces);
    error = true;
  }

  if (!bad_faces1.empty()) {
    os_ << "ERROR: found mis-oriented ghost faces:";
    writeList_(bad_faces1);
    error = true;
  }

  return globalAny_(error);
}

// Check that ghost cells are exact copies of their master.  This means
// that the the GIDs of the nodes defining the cell are the same, including
// their order (face orientation).
//
template<class Mesh_type>
bool MeshAudit_Maps<Mesh_type>::check_cell_to_nodes_ghost_data() const
{
  const Epetra_Map &node_map = mesh_->getMap(AmanziMesh::Entity_kind::NODE, true);
  const Epetra_Map &cell_map_own = mesh_->getMap(AmanziMesh::Entity_kind::CELL, false);
  const Epetra_Map &cell_map_use = mesh_->getMap(AmanziMesh::Entity_kind::CELL, true);

  int ncell_own = cell_map_own.NumMyElements();
  int ncell_use = cell_map_use.NumMyElements();

  Entity_ID_List cnode;
  std::vector<Entity_ID> bad_cells;

  int maxnodes = 0;
  for (Entity_ID j = 0; j < ncell_own; ++j) {
    mesh_->getCellNodes(j, cnode);
    maxnodes = (cnode.size() > maxnodes) ? cnode.size() : maxnodes;
  }

  // parallel calls require uniformity
  int maxnodes_tmp(maxnodes);
  this->comm_->MaxAll(&maxnodes_tmp, &maxnodes, 1);

  Epetra_IntSerialDenseMatrix gids(ncell_use,maxnodes); // no Epetra_IntMultiVector :(

  // Create a matrix of the GIDs for all owned cells.
  for (Entity_ID j = 0; j < ncell_own; ++j) {
    mesh_->getCellNodes(j, cnode);
    for (int k = 0; k < cnode.size(); ++k)
      gids(j,k) = node_map.GID(cnode[k]);
    for (int k = cnode.size(); k < maxnodes; ++k)
      gids(j,k) = 0;
  }

  // Import these GIDs to all used cells; sets values on ghost cells.
  Epetra_Import importer(cell_map_use, cell_map_own);
  for (int k = 0; k < maxnodes; ++k) {
    Epetra_IntVector kgids_own(View, cell_map_own, gids[k]);
    Epetra_IntVector kgids_use(View, cell_map_use, gids[k]);
    kgids_use.Import(kgids_own, importer, Insert);
  }

  // Compare the ghost cell GIDs against the reference values just computed.
  for (Entity_ID j = ncell_own; j < ncell_use; ++j) {
    mesh_->getCellNodes(j, cnode);
    bool bad_data = false;
    for (int k = 0; k < cnode.size(); ++k)
      if (node_map.GID(cnode[k]) != gids(j,k)) bad_data = true;
    if (bad_data) {
      for (int k = 0; k < cnode.size(); ++k) {
        bool found = false;
        for (int l = 0; l < cnode.size(); ++l)
          if (node_map.GID(cnode[k]) == gids(j,l)) {
            found = true;
            break;
          }
        if (!found) {
          bad_cells.push_back(j);
          break;
        }
      }
    }
  }

  bool error = false;

  if (!bad_cells.empty()) {
    os_ << "ERROR: found bad data for ghost cells:";
    writeList_(bad_cells);
    os_ << "       The ghost cells don't have the same nodes as their respective masters." << std::endl;
    error = true;
  }

  return globalAny_(error);
}


// Check that ghost cells reference the same faces, and in the same
// order, as their master.  This is part of a ghost being an exact
// copies the master.  Even if the preceding ghost checks pass it is
// still possible for this one to fail.  In this case the ghost would
// reference a different face (by GID) but that face would reference
// the same nodes as the correct face.  So the two faces would be
// geometrically identical, including orientation, but be distinct.
//
template<class Mesh_type>
bool MeshAudit_Maps<Mesh_type>::check_cell_to_faces_ghost_data() const
{
  const Epetra_Map &face_map = mesh_->getMap(AmanziMesh::Entity_kind::FACE, true);
  const Epetra_Map &cell_map_own = mesh_->getMap(AmanziMesh::Entity_kind::CELL, false);
  const Epetra_Map &cell_map_use = mesh_->getMap(AmanziMesh::Entity_kind::CELL, true);

  int ncell_own = cell_map_own.NumMyElements();
  int ncell_use = cell_map_use.NumMyElements();

  Entity_ID_List cface;
  std::vector<Entity_ID> bad_cells;

  // Create a matrix of the GIDs for all owned cells.
  int maxfaces = 0;
  for (Entity_ID j = 0; j < ncell_own; ++j) {
    mesh_->getCellFaces(j, cface);
    maxfaces = (cface.size() > maxfaces) ? cface.size() : maxfaces;
  }

  // parallel calls require uniformity
  int maxfaces_tmp(maxfaces);
  this->comm_->MaxAll(&maxfaces_tmp, &maxfaces, 1);

  Epetra_IntSerialDenseMatrix gids(ncell_use,maxfaces); // no Epetra_IntMultiVector :(

  for (Entity_ID j = 0; j < ncell_own; ++j) {
    mesh_->getCellFaces(j, cface);
    for (int k = 0; k < cface.size(); ++k)
      gids(j,k) = face_map.GID(cface[k]);
    for (int k = cface.size(); k < maxfaces; ++k)
      gids(j,k) = 0;
  }

  // Import these GIDs to all used cells; sets values on ghost cells.
  Epetra_Import importer(cell_map_use, cell_map_own);
  for (int k = 0; k < maxfaces; ++k) {
    Epetra_IntVector kgids_own(View, cell_map_own, gids[k]);
    Epetra_IntVector kgids_use(View, cell_map_use, gids[k]);
    kgids_use.Import(kgids_own, importer, Insert);
  }

  // Compare the ghost cell GIDs against the reference values just computed.
  for (Entity_ID j = ncell_own; j < ncell_use; ++j) {
    mesh_->getCellFaces(j, cface);
    bool bad_data = false;
    for (int k = 0; k < cface.size(); ++k)
      if (face_map.GID(cface[k]) != gids(j,k)) bad_data = true;
    if (bad_data) bad_cells.push_back(j);
  }

  bool error = false;

  if (!bad_cells.empty()) {
    os_ << "ERROR: found bad data for ghost cells:";
    writeList_(bad_cells);
    os_ << "       The ghost cells are not exact copies of their master." << std::endl;
    error = true;
  }

  return globalAny_(error);
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
template<class Mesh_type>
bool MeshAudit_Sets<Mesh_type>::check_node_set_ids() const
{
  return check_get_set_ids(Entity_kind::NODE);
}

template<class Mesh_type>
bool MeshAudit_Sets<Mesh_type>::check_face_set_ids() const
{
  return check_get_set_ids(Entity_kind::FACE);
}

template<class Mesh_type>
bool MeshAudit_Sets<Mesh_type>::check_cell_set_ids() const
{
  return check_get_set_ids(Entity_kind::CELL);
}

template<class Mesh_type>
bool MeshAudit_Sets<Mesh_type>::check_get_set_ids(Entity_kind kind) const
{
  // This needs some work.  Since we lazily resolve sets on demand, there is no
  // real way to know what sets the user actually intends to use.
  return false;

  // bool error = false;

  // // Get the number of sets.
  // int nset;
  // try {
  //   nset = mesh_->num_sets(kind); // this may fail
  // } catch (...) {
  //   os_ << "ERROR: caught exception from num_sets()" << std::endl;
  //   error = true;
  // }
  // error = globalAny_(error);
  // if (error) return error;

  // // Get the vector of set IDs.
  // std::vector<Set_ID> sids(nset, std::numeric_limits<unsigned int>::max());
  // try {
  //   mesh_->get_set_ids(kind, &sids); // this may fail
  // } catch (...) {
  //   os_ << "ERROR: caught exception from get_set_ids()" << std::endl;
  //   error = true;
  // }
  // error = globalAny_(error);
  // if (error) return error;

  // // Check to see that set ID values were actually assigned.  This assumes
  // // UINT_MAX is not a valid set ID.  This is a little iffy; perhaps 0 should
  // // be declared as an invalid set ID instead (the case for ExodusII), or
  // // perhaps we should just skip this check.
  // bool bad_data = false;
  // for (int j = 0; j < nset; ++j)
  //   if (sids[j] == std::numeric_limits<unsigned int>::max()) bad_data = true;
  // if (bad_data) {
  //   os_ << "ERROR: get_set_ids() failed to set all values" << std::endl;
  //   error = true;
  // }
  // error = globalAny_(error);
  // if (error) return error;

  // // Verify that the vector of set IDs contains no duplicates.
  // if (!this->areDistinctValues_(sids)) {
  //   os_ << "ERROR: get_set_ids() returned duplicate IDs" << std::endl;
  //   // it would be nice to output the duplicates
  //   error = true;
  // }
  // error = globalAny_(error);
  // if (error) return error;

  // // In parallel, verify that each process returns the exact same result.
  // if (this->comm_->NumProc() > 1) {
  //   // Check the number of sets are the same.
  //   this->comm_->Broadcast(&nset, 1, 0);
  //   if (nset != mesh_->num_sets(kind)) {
  //     os_ << "ERROR: inconsistent num_sets() value" << std::endl;
  //     error = true;
  //   }
  //   error = globalAny_(error);

  //   if (!error) {
  //     // Broadcast the set IDs on processor 0.
  //   std::vector<Set_ID> sids(nset, std::numeric_limits<unsigned int>::max());
  //   mesh_->get_set_ids(kind, &sids);
  //     int *sids0 = new int[nset];
  //     for (int j = 0; j < nset; ++j) sids0[j] = sids[j];
  //     this->comm_->Broadcast(sids0, nset, 0);

  //     // Check the set IDs, using the vector on process 0 as the reference.
  //     bool bad_data = false;
  //     for (int j = 0; j < nset; ++j)
  //       if (sids[j] != sids0[j]) bad_data = true;
  //     if (bad_data) {
  //       os_ << "ERROR: get_set_ids() returned inconsistent values" << std::endl;
  //       error = true;
  //     }
  //     delete [] sids0;
  //   }
  // }

  // return globalAny_(error);
}

// Check that valid_set_id() returns the correct results.
// This is a collective test, returning a collective pass/fail result.
//
template<class Mesh_type>
bool MeshAudit_Sets<Mesh_type>::check_valid_node_set_id() const
{
  return check_valid_set_id(Entity_kind::NODE);
}

template<class Mesh_type>
bool MeshAudit_Sets<Mesh_type>::check_valid_face_set_id() const
{
  return check_valid_set_id(Entity_kind::FACE);
}

template<class Mesh_type>
bool MeshAudit_Sets<Mesh_type>::check_valid_cell_set_id() const
{
  return check_valid_set_id(Entity_kind::CELL);
}

template<class Mesh_type>
bool MeshAudit_Sets<Mesh_type>::check_valid_set_id(Entity_kind kind) const
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
  //   writeList_(bad_sids1);
  //   error = true;
  // }
  // if (!bad_sids2.empty()) {
  //   os_ << "ERROR: valid_set_id() returned true for invalid set IDs:";
  //   writeList_(bad_sids2);
  //   error = true;
  // }
  // return globalAny_(error);
}

// For each set, check that get_set successfully returns valid references to
// local entities, without duplicates, and that the used set is consistent
// with the owned set.
//
template<class Mesh_type>
bool MeshAudit_Sets<Mesh_type>::check_node_sets() const
{
  return check_sets(Entity_kind::NODE,
                    mesh_->getMap(AmanziMesh::Entity_kind::NODE, false),
                    mesh_->getMap(AmanziMesh::Entity_kind::NODE, true));
}

template<class Mesh_type>
bool MeshAudit_Sets<Mesh_type>::check_face_sets() const
{
  return check_sets(Entity_kind::FACE,
                    mesh_->getMap(AmanziMesh::Entity_kind::FACE, false),
                    mesh_->getMap(AmanziMesh::Entity_kind::FACE, true));
}

template<class Mesh_type>
bool MeshAudit_Sets<Mesh_type>::check_cell_sets() const
{
  return check_sets(Entity_kind::CELL,
                    mesh_->getMap(AmanziMesh::Entity_kind::CELL, false),
                    mesh_->getMap(AmanziMesh::Entity_kind::CELL, true));
}

template<class Mesh_type>
bool MeshAudit_Sets<Mesh_type>::check_sets(Entity_kind kind,
                           const Epetra_Map &map_own, const Epetra_Map &map_use) const
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
  //   bool bad_set = check_get_set(sids[n], kind, Parallel_type::OWNED,
  //       			 map_own) ||
  //                  check_get_set(sids[n], kind, Parallel_type::ALL,
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
template<class Mesh_type>
bool MeshAudit_Sets<Mesh_type>::check_get_set(Set_ID sid,
        Entity_kind kind,
        Parallel_type ptype,
        const Epetra_Map &map) const
{
  os_ << "WARNING: Checks on sets disabled until MeshAudit handles new set specification methods (Tkt #686)" << std::endl;
  return false;

  // Get the size of the set.
  try {
    std::string set_name = mesh_->geometric_model()->FindRegion(sid)->name();
    mesh_->getSetSize(set_name, kind, ptype); // this may fail
  } catch (...) {
    os_ << "  ERROR: caught exception from getSetSize()" << std::endl;
    return true;
  }

  // Get the set.
  Entity_ID_List set;
  try {
    std::string set_name = mesh_->geometric_model()->FindRegion(sid)->name();
    mesh_->getSetEntities(set_name, kind, ptype, &set);  // this may fail
  } catch (...) {
    os_ << "  ERROR: caught exception from get_set()" << std::endl;
    return true;
  }

  // Check that all values were assigned.
  bool bad_data = false;
  for (int j = 0; j < set.size(); ++j)
    if (set[j] == std::numeric_limits<unsigned int>::max()) bad_data = true;
  if (bad_data) {
    os_ << "  ERROR: not all values assigned by get_set()" << std::endl;
    return true;
  }

  // Check that the LIDs in the set belong to the map.
  std::vector<Entity_ID> bad_LIDs;
  for (int j = 0; j < set.size(); ++j)
    if (!map.MyLID(set[j])) bad_LIDs.push_back(set[j]);
  if (!bad_LIDs.empty()) {
    os_ << "  ERROR: set contains invalid LIDs:";
    writeList_(bad_LIDs);
    return true;
  }

  // Check that there are no duplicates in the set.
  if (!this->areDistinctValues_(set)) {
    os_ << "  ERROR: set contains duplicate LIDs." << std::endl;
    // it would be nice to output the duplicates
    return true;
  }

  return false;
}

// The correct used set is completely determined by the owned set.  This test
// verifies that the used set is what it should be, considering the owned set
// as definitive.  Note that we do not require the vector of used set LIDs to
// extend the vector of owned set LIDs; the values in each list can be
// presented in any order.  This is a collective test, returning a collective
// pass/fail result.
//
template<class Mesh_type>
bool MeshAudit_Sets<Mesh_type>::check_used_set(Set_ID sid,
        Entity_kind kind,
        const Epetra_Map &map_own,
        const Epetra_Map &map_use) const
{
  os_ << "WARNING: Checks on sets disabled until MeshAudit handles new set specification methods (Tkt #686)" << std::endl;
  return false;

  std::string set_name = mesh_->geometric_model()->FindRegion(sid)->name();

  if (this->comm_->NumProc() == 1) {
    // In serial, the owned and used sets should be identical.

    int n = mesh_->getSetSize(set_name, kind, Parallel_type::OWNED);
    Entity_ID_List set_own;
    mesh_->getSetEntities(set_name, kind, Parallel_type::OWNED, &set_own);

    // Set sizes had better be the same.
    if (mesh_->getSetSize(set_name, kind, Parallel_type::ALL) !=
        set_own.size()) {
      os_ << "  ERROR: owned and used set sizes differ" << std::endl;
      return true;
    }

    // Verify that the two sets are identical.
    auto set_use = mesh_->getSetEntities(set_name, kind, Parallel_type::ALL);
    bool bad_data = false;
    for (int j = 0; j < n; ++j)
      if (set_use[j] != set_own[j]) bad_data = true;
    if (bad_data) {
      os_ << "  ERROR: owned and used sets differ" << std::endl;
      return true;
    }

    return false;

  } else {

    auto set_own = mesh_->getSetEntities(set_name, kind, Parallel_type::OWNED);
    auto set_use = mesh_->getSetEntities(set_name, kind, Parallel_type::ALL);

    // Tag all LIDs in the used map that should belong to the used set;
    // the owned set LIDs are taken as definitive.
    Epetra_IntVector tag_use(map_use); // fills with zero values
    int *tag_data;
    tag_use.ExtractView(&tag_data);
    Epetra_IntVector tag_own(View, map_own, tag_data);
    for (int j = 0; j < set_own.size(); ++j) tag_own[set_own[j]] = 1;
    Epetra_Import importer(map_use, map_own);
    tag_use.Import(tag_own, importer, Insert);

    // Now untag all the LIDs that belong to the used set.  If things
    // are correct, the tag vector will be filled with zeros afterwards.
    for (int j = 0; j < set_use.size(); ++j) --tag_use[set_use[j]];

    bool error = false;

    // Check for negative tag values;
    // these mark used LIDs that shouldn't be in the set but are.
    std::vector<Entity_ID> bad_LIDs;
    for (int j = 0; j < set_use.size(); ++j)
      if (tag_use[j] < 0) bad_LIDs.push_back(j);
    if (!bad_LIDs.empty()) {
      os_ << "  ERROR: found used LIDs that belong to the set but shouldn't:";
      writeList_(bad_LIDs);
      error = true;
    }

    // Check for positive tag values;
    // these mark used LIDs that should be in the set but aren't.
    bad_LIDs.resize(0);
    for (int j = 0; j < set_own.size(); ++j)
      if (tag_use[j] > 0) bad_LIDs.push_back(j);
    if (!bad_LIDs.empty()) {
      os_ << "  ERROR: found used LIDs that should belong to set but don't:";
      writeList_(bad_LIDs);
      error = true;
    }

    return globalAny_(error);
  }
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
template<class Mesh_type>
bool MeshAudit_Maps<Mesh_type>::check_face_partition() const
{
  // Mark all the faces contained by owned cells.
  bool owned[nfaces_all_];
  for (int j = 0; j < nfaces_all_; ++j) owned[j] = false;
  Entity_ID_List cface;
  for (Entity_ID j = 0; j < mesh_->getMap(AmanziMesh::Entity_kind::CELL, false).NumMyElements(); ++j) {
    mesh_->getCellFaces(j, cface);
    for (int k = 0; k < cface.size(); ++k) owned[cface[k]] = true;
  }

  // Verify that every owned face has been marked as belonging to an owned cell.
  std::vector<Entity_ID> bad_faces;
  for (Entity_ID j = 0; j < mesh_->getMap(AmanziMesh::Entity_kind::FACE, false).NumMyElements(); ++j)
    if (!owned[j]) bad_faces.push_back(j);

  if (!bad_faces.empty()) {
    os_ << "ERROR: found orphaned owned faces:";
    writeList_(bad_faces);
    os_ << "       Process doesn't own either of the cells sharing the face." << std::endl;
    return true;
  }

  return false;
}


template<class Mesh_type>
bool MeshAudit_Maps<Mesh_type>::check_node_partition() const
{
  // Mark all the nodes contained by owned cells.
  bool owned[nnodes_all_];
  for (int j = 0; j < nnodes_all_; ++j) owned[j] = false;
  Entity_ID_List cnode;
  for (Entity_ID j = 0; j < mesh_->getMap(AmanziMesh::Entity_kind::CELL, false).NumMyElements(); ++j) {
    mesh_->getCellNodes(j, cnode);
    for (int k = 0; k < cnode.size(); ++k) owned[cnode[k]] = true;
  }

  // Verify that every owned node has been marked as belonging to an owned cell.
  std::vector<Entity_ID> bad_nodes;
  for (Entity_ID j = 0; j < mesh_->getMap(AmanziMesh::Entity_kind::NODE, false).NumMyElements(); ++j)
    if (!owned[j]) bad_nodes.push_back(j);

  if (!bad_nodes.empty()) {
    os_ << "ERROR: found orphaned owned nodes:";
    writeList_(bad_nodes);
    os_ << "       Process doesn't own any of the cells containing the node." << std::endl;
    return true;
  }

  return false;
}


// Returns true if the values in the list are distinct -- no repeats.
template<class Mesh_type>
bool MeshAudit_Base<Mesh_type>::areDistinctValues_(const Entity_ID_List &list) const
{
  std::vector<Entity_ID> copy(list.data(), list.data()+list.size());
  sort(copy.begin(), copy.end());
  return (adjacent_find(copy.begin(),copy.end()) == copy.end());
}


// Returns 1 if the face node lists fnode1 and fnode2 describe the same face
// with the same orientation.  Returns -1 if the lists describe the same face
// but with opposite orientations.  Returns 0 if the lists describe different
// faces.  Implicitly assumes non-degenerate faces; the results are not
// reliable otherwise.
template<class Mesh_type>
int MeshAudit_Base<Mesh_type>::isSameFace_(const Entity_ID_List fnode1, const Entity_ID_List fnode2) const
{
  int nn = fnode1.size();

  if (nn != fnode2.size()) return 0;

  // Locate position in fnode1 of fnode2[0].
  int i, n;
  for (i = 0, n = -1; i < nn; ++i)
    if (fnode1[i] == fnode2[0]) {
      n = i;
      break;
    }
  if (n == -1) return 0; // did not find it -- different faces

  if (nn == 2) {
    // These are edges in a 2D mesh

    if (n == 0 && fnode1[1] == fnode2[1]) return 1;
    if (n == 1 && fnode1[0] == fnode2[1]) return -1;

    if (n == 0 && fnode1[1] == fnode2[1]) return 1;
    if (n == 1 && fnode1[0] == fnode2[1]) return -1;

  } else {
    for (i = 1; i < nn; ++i)
      if (fnode1[(n+i)%nn] != fnode2[i]) break;
    if (i == nn) return 1;  // they match

    // Modify the permutation to reverse the orientation of fnode1.
    for (i = 1; i < nn; ++i)
      if (fnode1[(n-i+nn)%nn] != fnode2[i]) break;
    if (i == nn) return -1;   // matched nodes but orientation is reversed
  }
  return 0; // different faces
}

template<class Mesh_type>
void MeshAudit_Base<Mesh_type>::writeList_(const std::vector<Entity_ID>& list) const
{
  int num_out = std::min((unsigned int) list.size(), this->MAX_OUT);
  for (int i = 0; i < num_out; ++i) os_ << " " << list[i];
  if (num_out < list.size()) os_ << " [" << list.size()-num_out << " items omitted]";
  os_ << std::endl;
}

template<class Mesh_type>
bool MeshAudit_Base<Mesh_type>::globalAny_(bool value) const
{
  int lval=value, gval;
  this->comm_->MaxAll(&lval, &gval, 1);
  return gval;
}


} // namespace Impl
} // namespace AmanziMesh
} // namespace Amanzi

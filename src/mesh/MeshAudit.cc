#include "MeshAudit.hh"

#include <algorithm>
#include <cfloat>

#include <boost/graph/topological_sort.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"

#include "Geometry.hh"

#include <iostream>
#include <iomanip>

using namespace std;
using namespace boost;

namespace Amanzi
{

  MeshAudit:: MeshAudit(Teuchos::RCP<AmanziMesh::MeshFramework> &mesh_, std::ostream& os_) :
      mesh(mesh_),
      comm_(mesh_->get_comm()),
      MyPID(mesh_->get_comm()->MyPID()),
      os(os_),
      nnode(mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_type::ALL)),
      nface(mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL)),
      ncell(mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL)),
      MAX_OUT(5)
    { create_test_dependencies(); }

// Verify runs all the tests in an order that respects the dependencies
// between tests.  If a particular test fails, all other tests that have it
// as a pre-requisite are skipped.  It is important that each test return
// a collective fail/pass result in parallel, so that all processes proceed
// through the tests in lockstep.

int MeshAudit::Verify() const
{
  int status = 0;

  typedef Graph::vertex_descriptor GraphVertex;
  std::list<GraphVertex> run_order;
  topological_sort(g, std::front_inserter(run_order));

  mark_do_not_run vis;

  for (auto itr = run_order.begin(); itr != run_order.end(); ++itr) {
    if (g[*itr].run) {
      os << "Checking " << g[*itr].name << " ..." << std::endl;
      if (((*this).*(g[*itr].test))()) {
        status = 1;
        os << "  Test failed!" << std::endl;
        breadth_first_search(g, *itr, visitor(vis));
      }
    } else {
      os << "Skipping " << g[*itr].name << " check because of previous failures." << std::endl;
    }
  }

  return status;
}

// This creates the dependency graph for the individual tests.  Adding a new
// test to the graph is a simple matter of adding a new stanza of the form
//
//   Graph::vertex_descriptor my_test_handle = add_vertex(g);
//   g[my_test_handle].name = "description of my test";
//   g[my_test_handle].test = &MeshAudit::my_test;
//   add_edge(other_test_handle, my_test_handle, g);
//
// The last line specifies that other_test_handle is a pre-requisite for
// my_test_handle.  There may be multiple pre-requisites or none.

void MeshAudit::create_test_dependencies()
{
  // Entity_counts tests
  Graph::vertex_descriptor test01 = add_vertex(g);
  g[test01].name = "entity_counts";
  g[test01].test = &MeshAudit::check_entity_counts;

  // Cell_to_nodes tests
  Graph::vertex_descriptor test02 = add_vertex(g);
  g[test02].name = "cell_to_nodes";
  g[test02].test = &MeshAudit::check_cell_to_nodes;

  Graph::vertex_descriptor test03 = add_vertex(g);
  g[test03].name = "node references by cells";
  g[test03].test = &MeshAudit::check_node_refs_by_cells;
  add_edge(test02, test03, g);

  // Cell_to_faces tests
  Graph::vertex_descriptor test04 = add_vertex(g);
  g[test04].name = "cell_to_faces";
  g[test04].test = &MeshAudit::check_cell_to_faces;

  Graph::vertex_descriptor test05 = add_vertex(g);
  g[test05].name = "face references by cells";
  g[test05].test = &MeshAudit::check_face_refs_by_cells;
  add_edge(test04, test05, g);

  // face_to_nodes tests
  Graph::vertex_descriptor test06 = add_vertex(g);
  g[test06].name = "face_to_nodes";
  g[test06].test = &MeshAudit::check_face_to_nodes;

  Graph::vertex_descriptor test08 = add_vertex(g);
  g[test08].name = "node references by faces";
  g[test08].test = &MeshAudit::check_node_refs_by_faces;
  add_edge(test06, test08, g);

  // cell_to_face_dirs tests
  Graph::vertex_descriptor test09 = add_vertex(g);
  g[test09].name = "cell_to_face_dirs";
  g[test09].test = &MeshAudit::check_cell_to_face_dirs;

  // cell degeneracy test
  Graph::vertex_descriptor test10 = add_vertex(g);
  g[test10].name = "topological non-degeneracy of cells";
  g[test10].test = &MeshAudit::check_cell_degeneracy;
  add_edge(test02, test10, g);

  // Consistency between the various mesh connectivity data.
  Graph::vertex_descriptor test11 = add_vertex(g);
  g[test11].name = "consistency of mesh connectivity data";
  g[test11].test = &MeshAudit::check_cell_to_faces_to_nodes;
  add_edge(test02, test11, g);
  add_edge(test04, test11, g);
  add_edge(test06, test11, g);
  add_edge(test09, test11, g);
  add_edge(test10, test11, g);

  // node_to_coordinates tests
  Graph::vertex_descriptor test12 = add_vertex(g);
  g[test12].name = "node_to_coordinates";
  g[test12].test = &MeshAudit::check_node_to_coordinates;

  // cell_to_coordinates tests
  Graph::vertex_descriptor test13 = add_vertex(g);
  g[test13].name = "cell_to_coordinates";
  g[test13].test = &MeshAudit::check_cell_to_coordinates;
  add_edge(test02, test13, g);
  add_edge(test12, test13, g);


  // face_to_coordinates tests
  Graph::vertex_descriptor test14 = add_vertex(g);
  g[test14].name = "face_to_coordinates";
  g[test14].test = &MeshAudit::check_face_to_coordinates;
  add_edge(test06, test14, g);
  add_edge(test12, test14, g);

  // cell topology/geometry test
  Graph::vertex_descriptor test15 = add_vertex(g);
  g[test15].name = "cell geometry";
  g[test15].test = &MeshAudit::check_cell_geometry;
  add_edge(test10, test15, g);
  add_edge(test13, test15, g);

#if 0
  // map tests
  Graph::vertex_descriptor test16 = add_vertex(g);
  g[test16].name = "owned and overlap node maps";
  g[test16].test = &MeshAudit::check_node_maps;

  Graph::vertex_descriptor test17 = add_vertex(g);
  g[test17].name = "owned and overlap face maps";
  g[test17].test = &MeshAudit::check_face_maps;

  Graph::vertex_descriptor test18 = add_vertex(g);
  g[test18].name = "owned and overlap cell maps";
  g[test18].test = &MeshAudit::check_cell_maps;

  // ghost data tests
  Graph::vertex_descriptor test19 = add_vertex(g);
  g[test19].name = "node_to_coordinates ghost data";
  g[test19].test = &MeshAudit::check_node_to_coordinates_ghost_data;
  add_edge(test12, test19, g);
  add_edge(test16, test19, g);

  Graph::vertex_descriptor test20 = add_vertex(g);
  g[test20].name = "face_to_nodes ghost data";
  g[test20].test = &MeshAudit::check_face_to_nodes_ghost_data;
  add_edge(test06, test20, g);
  add_edge(test16, test20, g);
  add_edge(test17, test20, g);

  Graph::vertex_descriptor test21 = add_vertex(g);
  g[test21].name = "cell_to_nodes ghost data";
  g[test21].test = &MeshAudit::check_cell_to_nodes_ghost_data;
  add_edge(test02, test21, g);
  add_edge(test16, test21, g);
  add_edge(test18, test21, g);

  Graph::vertex_descriptor test22 = add_vertex(g);
  g[test22].name = "cell_to_faces ghost data";
  g[test22].test = &MeshAudit::check_cell_to_faces_ghost_data;
  add_edge(test04, test22, g);
  add_edge(test17, test22, g);
  add_edge(test18, test22, g);

  // node set data tests
  Graph::vertex_descriptor test23 = add_vertex(g);
  g[test23].name = "node set IDs";
  g[test23].test = &MeshAudit::check_node_set_ids;

  Graph::vertex_descriptor test24 = add_vertex(g);
  g[test24].name = "node sets";
  g[test24].test = &MeshAudit::check_node_sets;
  add_edge(test16, test24, g);
  add_edge(test23, test24, g);

  Graph::vertex_descriptor test25 = add_vertex(g);
  g[test25].name = "valid node set IDs";
  g[test25].test = &MeshAudit::check_valid_node_set_id;
  add_edge(test23, test25, g);

  // face set data tests
  Graph::vertex_descriptor test26 = add_vertex(g);
  g[test26].name = "face set IDs";
  g[test26].test = &MeshAudit::check_face_set_ids;

  Graph::vertex_descriptor test27 = add_vertex(g);
  g[test27].name = "face sets";
  g[test27].test = &MeshAudit::check_face_sets;
  add_edge(test17, test27, g);
  add_edge(test26, test27, g);

  Graph::vertex_descriptor test28 = add_vertex(g);
  g[test28].name = "valid face set IDs";
  g[test28].test = &MeshAudit::check_valid_face_set_id;
  add_edge(test26, test28, g);


  // cell set data tests
  Graph::vertex_descriptor test29 = add_vertex(g);
  g[test29].name = "cell set IDs";
  g[test29].test = &MeshAudit::check_cell_set_ids;

  Graph::vertex_descriptor test30 = add_vertex(g);
  g[test30].name = "cell sets";
  g[test30].test = &MeshAudit::check_cell_sets;
  add_edge(test18, test30, g);
  add_edge(test29, test30, g);

  Graph::vertex_descriptor test31 = add_vertex(g);
  g[test31].name = "valid cell set IDs";
  g[test31].test = &MeshAudit::check_valid_cell_set_id;
  add_edge(test29, test31, g);

  // partition tests
  Graph::vertex_descriptor test32 = add_vertex(g);
  g[test32].name = "face partition";
  g[test32].test = &MeshAudit::check_face_partition;
  add_edge(test17, test32, g);
  add_edge(test18, test32, g);

  Graph::vertex_descriptor test33 = add_vertex(g);
  g[test33].name = "node partition";
  g[test33].test = &MeshAudit::check_node_partition;
  add_edge(test16, test33, g);
  add_edge(test18, test33, g);
#endif
}

////////////////////////////////////////////////////////////////////////////////
//
// Tests (and auxillary functions) follow.  Tests must be const functions
// that take no arguments and that return a bool result: true if an error
// was found, otherwise false.  Tests must be considered collective procedures
// when in a parallel context, and so return a comm_on collective result.
// This is easily done using the function global_any(bool) which returns true
// on all processes if the argument is true on any of the processes.
//
////////////////////////////////////////////////////////////////////////////////

// The count_entities method should return values that match the number of
// elements in the corresponding Epetra_Maps.  This applies to nodes, faces
// and cells, including ghosts and not.  Here the maps are considered to be
// the authoritative source of information.  A positive value is returned
// if any discrepancy is found, but it is safe to perform other tests as
// they do not use these methods.

/*
bool MeshAudit::check_entity_counts() const
{
  int n, nref;
  bool error = false;

  // Check the number of owned nodes.
  n = mesh->getNumEntities(AmanziMesh::Entity_kind::NODE,AmanziMesh::Parallel_type::OWNED);
  nref = mesh->node_map(false).NumMyElements();
  if (n != nref) {
    os << ": ERROR: getNumEntities(NODE,OWNED)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  // Check the number of used nodes.
  n = mesh->getNumEntities(AmanziMesh::Entity_kind::NODE,AmanziMesh::Parallel_type::ALL);
  nref = mesh->node_map(true).NumMyElements();
  if (n != nref) {
    os << "ERROR: getNumEntities(NODE,ALL)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  // Check the number of owned faces.
  n = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE,AmanziMesh::Parallel_type::OWNED);
  nref = mesh->face_map(false).NumMyElements();
  if (n != nref) {
    os << "ERROR: getNumEntities(FACE,OWNED)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  // Check the number of used faces.
  n = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE,AmanziMesh::Parallel_type::ALL);
  nref = mesh->face_map(true).NumMyElements();
  if (n != nref) {
    os << "ERROR: getNumEntities(FACE,ALL)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  // Check the number of owned cells.
  n = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL,AmanziMesh::Parallel_type::OWNED);
  nref = mesh->cell_map(false).NumMyElements();
  if (n != nref) {
    os << "ERROR: getNumEntities(CELL,OWNED)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  // Check the number of used cells.
  n = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL,AmanziMesh::Parallel_type::ALL);
  nref = mesh->cell_map(true).NumMyElements();
  if (n != nref) {
    os << "ERROR: getNumEntities(CELL,ALL)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  return global_any(error);
}
*/

// Check that cell_to_nodes successfully returns valid references to local
// nodes.  Here the std::vector accessor is used.  The consistency of
// the alternative accessors is checked elsewhere.  A nonzero return value
// signals an error, and further tests using its data should be avoided.

bool MeshAudit::check_cell_to_nodes() const
{
  AmanziMesh::Entity_ID_List bad_cells, bad_cells1;
  AmanziMesh::Entity_ID_List free_nodes;
  AmanziMesh::Entity_ID_List cnode;

  for (AmanziMesh::Entity_ID j = 0; j < ncell; ++j) {
    try {
      mesh->getCellNodes(j, &cnode); // this may fail
      bool invalid_refs = false;
      for (int k = 0; k < cnode.size(); ++k) {
        if (cnode[k] >= nnode) invalid_refs = true;
      }
      if (invalid_refs) bad_cells.push_back(j);
    } catch (...) {
      bad_cells1.push_back(j);
    }
  }

  bool error = false;

  if (!bad_cells.empty()) {
    os << "ERROR: invalid nodes referenced by cells:";
    write_list(bad_cells, MAX_OUT);
    error = true;
  }

  if (!bad_cells1.empty()) {
    os << "ERROR: caught exception for cells:";
    write_list(bad_cells1, MAX_OUT);
    error = true;
  }

  return global_any(error);
}

// Check that every node is referenced by at least one cell.  This assumes
// that cell_to_nodes have been verified to return valid data.  A nonzero
// return value indicates that one or more nodes are not attached to any cell.

bool MeshAudit::check_node_refs_by_cells() const
{
  AmanziMesh::Entity_ID_List cnode;
  vector<bool> ref(nnode, false);

  for (AmanziMesh::Entity_ID j = 0; j < ncell; ++j) {
    mesh->getCellNodes(j, &cnode); // this should not fail
    for (int k = 0; k < cnode.size(); ++k) ref[cnode[k]] = true;
  }

  AmanziMesh::Entity_ID_List free_nodes;
  for (int j = 0; j < nnode; ++j) {
    if (!ref[j]) free_nodes.push_back(j);
  }

  bool error = false;

  if (!free_nodes.empty()) {
    os << "ERROR: found unreferenced nodes:";
    write_list(free_nodes, MAX_OUT);
    error = true;
  }

  return global_any(error);
}

// Check that cell_to_faces successfully returns valid references to local
// faces.  Here the std::vector accessor is used.  The consistency of
// the alternative accessors is checked elsewhere.  A nonzero return value
// signals an error, and further tests using its data should be avoided.

bool MeshAudit::check_cell_to_faces() const
{
  AmanziMesh::Entity_ID_List bad_cells, bad_cells1;
  AmanziMesh::Entity_ID_List bad_faces;
  AmanziMesh::Entity_ID_List free_faces;
  AmanziMesh::Entity_ID_List cface;

  for (AmanziMesh::Entity_ID j = 0; j < ncell; ++j) {
    try {
      mesh->getCellFaces(j, &cface); // this may fail
      bool invalid_refs = false;
      for (int k = 0; k < cface.size(); ++k) {
        if (cface[k] >= nface) invalid_refs = true;
      }
      if (invalid_refs) bad_cells.push_back(j);
    } catch (...) {
      bad_cells1.push_back(j);
    }
  }

  bool error = false;

  if (!bad_cells.empty()) {
    os << "ERROR: invalid faces referenced by cells:";
    write_list(bad_cells, MAX_OUT);
    error = true;
  }

  if (!bad_cells1.empty()) {
    os << "ERROR: caught exception for cells:";
    write_list(bad_cells1, MAX_OUT);
    error = true;
  }

  return global_any(error);
}

// Check that every face is referenced by exactly one or two cells.  This
// assumes that cell_to_faces have been verified to return valid data.  A
// nonzero return value indicates that faces were found that do not belong
// to any cell or belong to more than two cells (bad topology).

bool MeshAudit::check_face_refs_by_cells() const
{
  AmanziMesh::Entity_ID_List cface;
  vector<unsigned int> refs(nface, 0);

  for (AmanziMesh::Entity_ID j = 0; j < ncell; ++j) {
    mesh->getCellFaces(j, &cface);
    for (int k = 0; k < cface.size(); ++k) (refs[cface[k]])++;
  }

  AmanziMesh::Entity_ID_List free_faces;
  AmanziMesh::Entity_ID_List bad_faces;

  for (int j = 0; j < nface; ++j) {
    if (refs[j] == 0)
      free_faces.push_back(j);
    else if (refs[j] > 2)
      bad_faces.push_back(j);
  }

  bool error = false;

  if (!free_faces.empty()) {
    os << "ERROR: found unreferenced faces:";
    write_list(free_faces, MAX_OUT);
    error = true;
  }

  if (!bad_faces.empty() && mesh->get_space_dimension() == mesh->get_manifold_dimension()) {
    os << "ERROR: found faces shared by more than two cells:";
    write_list(bad_faces, MAX_OUT);
    error = true;
  }

  return global_any(error);
}

// Check that face_to_nodes successfully returns valid references to local
// nodes.  Here the std::vector accessor is used.  The consistency of
// the alternative accessors is checked elsewhere.  A nonzero return value
// signals an error, and further tests using its data should be avoided.

bool MeshAudit::check_face_to_nodes() const
{
  AmanziMesh::Entity_ID_List bad_faces, bad_faces1;
  AmanziMesh::Entity_ID_List fnode;

  for (AmanziMesh::Entity_ID j = 0; j < nface; ++j) {
    try {
      mesh->getFaceNodes(j, &fnode); // this may fail
      bool invalid_refs = false;
      for (int k = 0; k < fnode.size(); ++k) {
        if (fnode[k] >= nnode) invalid_refs = true;
      }
      if (invalid_refs) bad_faces.push_back(j);
    } catch (...) {
      bad_faces1.push_back(j);
    }
  }

  bool error = false;

  if (!bad_faces.empty()) {
    os << "ERROR: invalid nodes referenced by faces:";
    write_list(bad_faces, MAX_OUT);
    error = true;
  }

  if (!bad_faces1.empty()) {
    os << "ERROR: caught exception for faces:";
    write_list(bad_faces1, MAX_OUT);
    error = true;
  }

  return global_any(error);
}

// Check that every node is referenced by at least one face.  This assumes
// that face_to_nodes has been verified to return valid data.  A nonzero
// return value indicates that one or more nodes are not attached to any face.

bool MeshAudit::check_node_refs_by_faces() const
{
  AmanziMesh::Entity_ID_List fnode;
  vector<bool> ref(nnode, false);

  for (AmanziMesh::Entity_ID j = 0; j < nface; ++j) {
    mesh->getFaceNodes(j, &fnode);
    for (int k = 0; k < fnode.size(); ++k) ref[fnode[k]] = true;
  }

  AmanziMesh::Entity_ID_List free_nodes;
  for (int j = 0; j < nnode; ++j) {
    if (!ref[j]) free_nodes.push_back(j);
  }

  bool error = false;

  if (!free_nodes.empty()) {
    os << "ERROR: found unreferenced nodes:";
    write_list(free_nodes, MAX_OUT);
    error = true;
  }

  return global_any(error);
}

// Check that cell_to_face_dirs successfully returns data for all cells and
// that the values are either +1 or -1. The std::vector-based method is used;
// the consistency of the alternative methods is checked elsewhere.  If this
// test fails, further tests using this data should be avoided.

bool MeshAudit::check_cell_to_face_dirs() const
{
  AmanziMesh::Entity_ID_List faces;
  vector<int> fdirs;
  AmanziMesh::Entity_ID_List bad_cells, bad_cells_exc;

  for (AmanziMesh::Entity_ID j = 0; j < ncell; ++j) {
    fdirs.assign(6, INT_MAX);
    try {
      mesh->getCellFacesAndDirs(j, &faces, &fdirs);  // this may fail
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
    os << "ERROR: inadmissable or no data for cells:";
    write_list(bad_cells, MAX_OUT);
    error = true;
  }

  if (!bad_cells_exc.empty()) {
    os << "ERROR: caught exception for cells:";
    write_list(bad_cells_exc, MAX_OUT);
    error = true;
  }

  return global_any(error);
}


// Check that cells are not topologically degenerate (repeated node index).
// If this test fails, many further tests involving the cell_to_nodes data
// should be avoided.

bool MeshAudit::check_cell_degeneracy() const
{
  os << "Checking cells for topological degeneracy ..." << std::endl;

  AmanziMesh::Entity_ID_List cnode;
  AmanziMesh::Entity_ID_List bad_cells;

  for (AmanziMesh::Entity_ID j = 0; j < ncell; ++j) {
    mesh->getCellNodes(j, &cnode); // should not fail
    if (!distinct_values(cnode)) bad_cells.push_back(j);
  }

  bool error = false;

  if (!bad_cells.empty()) {
    os << "ERROR: found topologically degenerate cells:";
    write_list(bad_cells, MAX_OUT);
    error = true;
  }

    return global_any(error);
}

// Check that cell_to_faces and face_to_nodes are returning the correct
// values by composing those maps and comparing the result against the
// result returned by cell_to_nodes and the local face numbering convention
// described by cell_topology::HexFaceVert.  Also check that the relative
// orientation value returned by cell_to_face_dirs is correct.  If this test
// fails, one or more of the methods cell_to_faces, face_to_nodes, and
// cell_to_face_dirs are returning incorrect results and that further tests
// using their values should be avoided.

bool MeshAudit::check_cell_to_faces_to_nodes() const
{
  AmanziMesh::Entity_ID_List cnode;
  AmanziMesh::Entity_ID_List cface;
  AmanziMesh::Entity_ID_List fnode_ref;
  AmanziMesh::Entity_ID_List fnode;
  vector<int> fdirs;
  AmanziMesh::Entity_ID_List bad_cells0;
  AmanziMesh::Entity_ID_List bad_cells1;

  for (AmanziMesh::Entity_ID j = 0; j < ncell; ++j) {
    AmanziMesh::Cell_type ctype = mesh->cell_get_type(j);

    // If this is a general, non-standard element there is nothing to
    // to check against

    if (ctype == AmanziMesh::Entity_kind::CELLTYPE_UNKNOWN || 
	ctype == AmanziMesh::POLYGON || 
	ctype == AmanziMesh::POLYHED)
      continue;

    mesh->getCellNodes(j, &cnode); // this should not fail

    mesh->getCellFacesAndDirs(j, &cface, &fdirs, true); // this should not fail

    bool bad_face = false;
    bool bad_dir  = false;

    if (cface.size() != AmanziMesh::nface_std[ctype])
      bad_face = true;
    else {

      for (int k = 0; k < cface.size(); ++k) {

	mesh->getFaceNodes(cface[k], &fnode); // this should not fail

	int nfn = AmanziMesh::nfnodes_std[ctype][k];

	if (fnode.size() != nfn) {
	  bad_face = true;
	  break;
	}

	fnode_ref.clear();
	for (int i = 0; i < nfn; ++i) {
	  int nodenum = AmanziMesh::fnodes_std[ctype][k][i];
	  fnode_ref.push_back(cnode[nodenum]);
	}

	int dir = same_face(fnode, fnode_ref); // should be the same face
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
    os << "ERROR: bad getCellFaces values for cells:";
    write_list(bad_cells0, MAX_OUT);
    error = true;
  }

  // algorithms on a non-manifold use multiple normals and special continuity equations
  // for fluxes, so that orientation does not play role. This may change in the future.
  if (!bad_cells1.empty() && mesh->get_space_dimension() == mesh->get_manifold_dimension()) {
    os << "ERROR: bad cell_get_face_dirs values for cells:";
    write_list(bad_cells1, MAX_OUT);
    error = true;
  }

  return global_any(error);
}

// Check that getNodeCoordinate successfully returns data for all nodes
// A negative return value signals a terminal error

bool MeshAudit::check_node_to_coordinates() const
{
  vector<double> x(3);
  AmanziMesh::Entity_ID_List bad_nodes, bad_nodes_exc;
  int spdim = mesh->get_space_dimension();
  AmanziMesh::Entity_ID_List bad_nodes0;

  for (AmanziMesh::Entity_ID j = 0; j < nnode; ++j) {
    AmanziGeometry::Point x0(spdim);
    x0.set(DBL_MAX);
    try {
      mesh->getNodeCoordinate(j, &x0); // this may fail
      bool bad_data = false;
      for (int k = 0; k < spdim; ++k)
        if (x0[k] == DBL_MAX) bad_data = true;
      if (bad_data) bad_nodes.push_back(j);
    } catch (...) {
      bad_nodes_exc.push_back(j);
    }
  }

  bool error = false;

  if (!bad_nodes.empty()) {
    os << "ERROR: missing coordinate data for nodes:";
    write_list(bad_nodes, MAX_OUT);
    error = true;
  }

  if (!bad_nodes_exc.empty()) {
    os << "ERROR: caught exception for nodes:";
    write_list(bad_nodes_exc, MAX_OUT);
    error = true;
  }

  return global_any(error);
}

// Check that cell_get_coordinates successfully returns data for all
// cells, and that the data is identical to that returned by composing
// with cell_to_nodes.  The consistency of the alternative accessors is checked
// elsewhere.  If this test fails, further tests using cell_to_coordinates data
// should be avoided.

bool MeshAudit::check_cell_to_coordinates() const
{
  int spdim = mesh->get_space_dimension();
  AmanziMesh::Entity_ID_List cnode;
  AmanziMesh::Entity_ID_List bad_cells, bad_cells_exc;

  for (AmanziMesh::Entity_ID j = 0; j < ncell; ++j) {
    vector<AmanziGeometry::Point> x0;
    try {
      mesh->cell_get_coordinates(j, &x0); // this may fail
      mesh->getCellNodes(j, &cnode); // this should not fail
      for (int k = 0; k < cnode.size(); ++k) {
	AmanziGeometry::Point xref(spdim);
	bool bad_data = false;
        mesh->getNodeCoordinate(cnode[k], &xref); // this should not fail
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
    os << "ERROR: bad cell_to_coordinates data for cells:";
    write_list(bad_cells, MAX_OUT);
    error = true;
  }

  if (!bad_cells_exc.empty()) {
    os << "ERROR: caught exception for cells:";
    write_list(bad_cells_exc, MAX_OUT);
    error = true;
  }

  return global_any(error);
}


// Check that face_get_coordinates successfully returns data for all
// faces, and that this data is identical to that returned by
// composing getNodeCoordinate with getFaceNodes.  This assumes
// that the std::vector form of getFaceNodes and pointer-based form
// of getNodeCoordinate have been verified to return valid data. A
// negative return value signals a terminal error.


bool MeshAudit::check_face_to_coordinates() const
{
  int spdim = mesh->get_space_dimension();
  AmanziMesh::Entity_ID_List fnode;
  AmanziMesh::Entity_ID_List bad_faces, bad_faces_exc;

  for (AmanziMesh::Entity_ID j = 0; j < nface; ++j) {
    try {
      vector<AmanziGeometry::Point> x0;
      bool bad_data = false;
      mesh->face_get_coordinates(j, &x0); // this may fail
      mesh->getFaceNodes(j, &fnode); // this should not fail
      for (int k = 0; k < fnode.size(); ++k) {
	AmanziGeometry::Point xref(spdim);
        mesh->getNodeCoordinate(fnode[k], &xref); // this should not fail
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
    os << "ERROR: bad face_to_coordinates data for faces:";
    write_list(bad_faces, MAX_OUT);
    error = true;
  }

  if (!bad_faces_exc.empty()) {
    os << "ERROR: caught exception for faces:";
    write_list(bad_faces_exc, MAX_OUT);
    error = true;
  }

  return global_any(error);
}


// The cells must not be degenerate, either topologically (repeated
// node index) or geometrically (with coincident nodes).


bool MeshAudit::check_cell_geometry() const
{
  os << "Checking cell geometry ..." << std::endl;
  AmanziGeometry::Point centroid;
  double hvol;
  AmanziMesh::Entity_ID_List bad_cells;

  for (AmanziMesh::Entity_ID j = 0; j < ncell; ++j) {
    hvol = mesh->getCellVolume(j);
      
    if (hvol <= 0.0) bad_cells.push_back(j);
  }

  bool error = false;

  if (!bad_cells.empty()) {
    os << "ERROR: found cells with non-positive volumes:";
    write_list(bad_cells, MAX_OUT);
    os << "       Cells either geometrically degenerate or have bad topology.";
    error = true;
  }

  return global_any(error);
}

#if 0

// We assume the maps returned by the map accessors are well-formed as
// Epetra_maps.  Here we are checking that they have the required
// characteristics, especially for the relationship between the owned
// and overlapped version of the maps.  Specifically, we require that
// the index base for all maps is 0, and that the owned maps (ones without
// ghost) are 1-1 maps.  In serial the owned and overlap (one with ghosts)
// maps should be the same.  In parallel, the overlap version should extend
// the owned version (owned.getEntityGID() = overlapped.getEntityGID() for all owned LIDs),
// and that the remaining overlapped LIDs (if any) refer to getEntityGIDs owned
// by other processes.  In addition there should be no duplicate getEntityGIDs on
// any processes map.

bool MeshAudit::check_node_maps() const
{
  return check_maps(mesh->node_map(false), mesh->node_map(true));
}

bool MeshAudit::check_face_maps() const
{
  return check_maps(mesh->face_map(false), mesh->face_map(true));
}

bool MeshAudit::check_cell_maps() const
{
  return check_maps(mesh->cell_map(false), mesh->cell_map(true));
}

bool MeshAudit::check_maps(const Epetra_Map &map_own, const Epetra_Map &map_use) const
{
  bool error = false;

  // Local index should start at 0.
  if (map_own.IndexBase() != 0) {
    os << "ERROR: the owned map's index base is not 0." << std::endl;
    error = true;
  }
  if (map_use.IndexBase() != 0) {
    os << "ERROR: the overlap map's index base is not 0." << std::endl;
    error = true;
  }

  // Check that the owned map is 1-1.
  if (!map_own.UniquegetEntityGIDs()) {
    os << "ERROR: owned map is not 1-to-1" << std::endl;
    error = true;
  }

  // Check that the global ID space is contiguously divided (but not
  // necessarily uniformly) across all processers

  std::vector<int> owned_getEntityGIDs(map_own.NumMyElements());
  for (int i = 0; i < map_own.NumMyElements(); i++)
    owned_getEntityGIDs[i] = map_own.getEntityGID(i);
  std::sort(owned_getEntityGIDs.begin(), owned_getEntityGIDs.end());

  for (int i = 0; i < map_own.NumMyElements()-1; i++) {
    int diff = owned_getEntityGIDs[i+1]-owned_getEntityGIDs[i];
    if (diff > 1) {
      os << "ERROR: owned map is not contiguous" << std::endl;
      os << "Global IDs jump from " << owned_getEntityGIDs[i] << " to " <<
          owned_getEntityGIDs[i+1] << std::endl;
      error = true;
    }
  }

  error = global_any(error);
  if (error) return error;

  if (comm_->NumProc() == 1)
  {

    // Serial or 1-process MPI

    if (!map_use.SameAs(map_own)) {
      os << "ERROR: the overlap map differs from the owned map (single process)." << std::endl;
      error = true;
    }

    return global_any(error);

  }
  else
  {

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
        if (map_use.getEntityGID(j) != map_own.getEntityGID(j)) bad_map = true;
    }
    if (bad_map) {
      os << "ERROR: overlap map does not extend the owned map." << std::endl;
      error = true;
    }

    error = global_any(error);
    if (error) return error;

    // Verify that the overlap indices are owned by other processes.
    int *gids = new int[num_ovl];
    int *pids = new int[num_ovl];
    int *lids = new int[num_ovl];
    for (int j = 0; j < num_ovl; ++j) gids[j] = map_use.getEntityGID(j+num_own);
    map_own.RemoteIDList(num_ovl, gids, pids, lids);
    bad_map = false;
    for (int j = 0; j < num_ovl; ++j)
      if (pids[j] < 0 || pids[j] == comm_->MyPID()) bad_map = true;
    if (bad_map) {
      os << "ERROR: invalid ghosts in overlap map." << std::endl;
      error = true;
    }

    // Look for duplicates among the overlap indices.
    vector<int> ovl_gids(gids, gids+num_ovl);
    sort(ovl_gids.begin(), ovl_gids.end());
    if (adjacent_find(ovl_gids.begin(),ovl_gids.end()) != ovl_gids.end()) {
      os << "ERROR: duplicate ghosts in overlap map." << std::endl;
      error = true;
    }

    delete [] lids;
    delete [] pids;
    delete [] gids;

    return global_any(error);
  }
}

// Check that ghost nodes are exact copies of their master.
// This simply means that they have the same coordinates.

bool MeshAudit::check_node_to_coordinates_ghost_data() const
{
  int spdim = mesh->get_space_dimension();

  const Epetra_Map &node_map_own = mesh->node_map(false);
  const Epetra_Map &node_map_use = mesh->node_map(true);

  int nnode_own = node_map_own.NumMyElements();
  int nnode_use = node_map_use.NumMyElements();

  AmanziGeometry::Point coord(spdim);
  AmanziMesh::Entity_ID_List bad_nodes;

  Epetra_MultiVector coord_use(node_map_use,spdim);
  double **data;
  coord_use.ExtractView(&data);
  Epetra_MultiVector coord_own(View, node_map_own, data, spdim);

  for (AmanziMesh::Entity_ID j = 0; j < nnode_own; ++j) {
    mesh->getNodeCoordinate(j, &coord);
    for (int k = 0; k < spdim; ++k) coord_own[k][j] = coord[k];
  }

  Epetra_Import importer(node_map_use, node_map_own);
  coord_use.Import(coord_own, importer, Insert);

  for (AmanziMesh::Entity_ID j = nnode_own; j < nnode_use; ++j) {
    mesh->getNodeCoordinate(j, &coord);
    bool bad_data = false;
    for (int k = 0; k < spdim; ++k)
      if (coord[k] != coord_use[k][j]) bad_data = true;
    if (bad_data) bad_nodes.push_back(j);
  }

  bool error = false;

  if (!bad_nodes.empty()) {
    os << "ERROR: found ghost nodes with incorrect coordinates:";
    write_list(bad_nodes, MAX_OUT);
    error = true;
  }

  return global_any(error);
}

// Check that ghost faces are exact copies of their master.  This means
// that the the getEntityGIDs of the nodes defining the face are the same, including
// their order (face orientation).

bool MeshAudit::check_face_to_nodes_ghost_data() const
{
  const Epetra_Map &node_map = mesh->node_map(true);
  const Epetra_Map &face_map_own = mesh->face_map(false);
  const Epetra_Map &face_map_use = mesh->face_map(true);

  int nface_own = face_map_own.NumMyElements();
  int nface_use = face_map_use.NumMyElements();

  AmanziMesh::Entity_ID_List fnode;
  AmanziMesh::Entity_ID_List bad_faces, bad_faces1, bad_faces2;

  // Create a matrix of the getEntityGIDs for all owned faces.

  int maxnodes = 0;
  for (AmanziMesh::Entity_ID j = 0; j < nface_own; ++j) {
    mesh->getFaceNodes(j, &fnode);
    maxnodes = fnode.size() > maxnodes ? fnode.size() : maxnodes;
  }

  // parallel calls require uniformity
  int maxnodes_tmp(maxnodes);
  comm_->MaxAll(&maxnodes_tmp, &maxnodes, 1);

  Epetra_IntSerialDenseMatrix gids(nface_use,maxnodes); // no Epetra_IntMultiVector :(
  for (AmanziMesh::Entity_ID j = 0; j < nface_own; ++j) {
    mesh->getFaceNodes(j, &fnode);
    for (int k = 0; k < fnode.size(); ++k)
      gids(j,k) = node_map.getEntityGID(fnode[k]);
    for (int k = fnode.size(); k < maxnodes; ++k)
      gids(j, k) = 0;
  }

  // Import these getEntityGIDs to all used faces; sets values on ghost faces.
  Epetra_Import importer(face_map_use, face_map_own);
  for (int k = 0; k < maxnodes; ++k) {
    Epetra_IntVector kgids_own(View, face_map_own, gids[k]);
    Epetra_IntVector kgids_use(View, face_map_use, gids[k]);
    kgids_use.Import(kgids_own, importer, Insert);
  }

  // Compare the ghost face getEntityGIDs against the reference values just computed.
  for (AmanziMesh::Entity_ID j = nface_own; j < nface_use; ++j) {
    mesh->getFaceNodes(j, &fnode);
    bool bad_data = false;
    for (int k = 0; k < fnode.size(); ++k)
      if (node_map.getEntityGID(fnode[k]) != gids(j,k)) { 
        bad_data = true;
      }
    if (bad_data) {
      // Determine just how bad the data is.
      AmanziMesh::Entity_ID_List fnode_ref(maxnodes);
      for (int k = 0; k < maxnodes; ++k) {
        fnode[k] = node_map.getEntityGID(fnode[k]);
        fnode_ref[k] = gids(j,k);
      }
      int n = same_face(fnode, fnode_ref);
      switch (n) {
      case 0: // completely bad -- different face
        bad_faces.push_back(j);

        std::cerr << "P " << comm_->MyPID() << ": ghost face " << j << " (getEntityGID "
                  << face_map_use.getEntityGID(j) << ")," <<
            " has different nodes than its master " << std::endl;
        std::cerr << "ghost face nodes (getEntityGIDs): ";
        for (int k = 0; k < fnode.size(); ++k)
          std::cerr << fnode[k] << " ";
        std::cerr << std::endl;
        std::cerr << "master face nodes (getEntityGIDs): ";
        for (int k = 0; k < fnode_ref.size(); ++k)
          std::cerr << fnode_ref[k] << " ";
        std::cerr << std::endl;
        break;
      case -1: // very bad -- same face but wrong orientation
        bad_faces1.push_back(j);

        std::cerr << "P " << comm_->MyPID() << ": ghost face " << j << ", has different orientation than its master " << std::endl;
        std::cerr << "ghost face nodes (getEntityGIDs): ";
        for (int k = 0; k < fnode.size(); ++k)
          std::cerr << fnode[k];
        std::cerr << std::endl;
        std::cerr << "master face nodes (getEntityGIDs): ";
        for (int k = 0; k < fnode_ref.size(); ++k)
          std::cerr << fnode_ref[k];
        std::cerr << std::endl;
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
    os << "ERROR: found bad data for ghost faces:";
    write_list(bad_faces, MAX_OUT);
    error = true;
  }

  if (!bad_faces1.empty()) {
    os << "ERROR: found mis-oriented ghost faces:";
    write_list(bad_faces1, MAX_OUT);
    error = true;
  }

  return global_any(error);
}

// Check that ghost cells are exact copies of their master.  This means
// that the the getEntityGIDs of the nodes defining the cell are the same, including
// their order (face orientation).

bool MeshAudit::check_cell_to_nodes_ghost_data() const
{
  const Epetra_Map &node_map = mesh->node_map(true);
  const Epetra_Map &cell_map_own = mesh->cell_map(false);
  const Epetra_Map &cell_map_use = mesh->cell_map(true);

  int ncell_own = cell_map_own.NumMyElements();
  int ncell_use = cell_map_use.NumMyElements();

  AmanziMesh::Entity_ID_List cnode;
  AmanziMesh::Entity_ID_List bad_cells;

  int maxnodes = 0;
  for (AmanziMesh::Entity_ID j = 0; j < ncell_own; ++j) {
    mesh->getCellNodes(j, &cnode);
    maxnodes = (cnode.size() > maxnodes) ? cnode.size() : maxnodes;
  }

  // parallel calls require uniformity
  int maxnodes_tmp(maxnodes);
  comm_->MaxAll(&maxnodes_tmp, &maxnodes, 1);

  Epetra_IntSerialDenseMatrix gids(ncell_use,maxnodes); // no Epetra_IntMultiVector :(

  // Create a matrix of the getEntityGIDs for all owned cells.
  for (AmanziMesh::Entity_ID j = 0; j < ncell_own; ++j) {
    mesh->getCellNodes(j, &cnode);
    for (int k = 0; k < cnode.size(); ++k)
      gids(j,k) = node_map.getEntityGID(cnode[k]);
    for (int k = cnode.size(); k < maxnodes; ++k)
      gids(j,k) = 0;
  }

  // Import these getEntityGIDs to all used cells; sets values on ghost cells.
  Epetra_Import importer(cell_map_use, cell_map_own);
  for (int k = 0; k < maxnodes; ++k) {
    Epetra_IntVector kgids_own(View, cell_map_own, gids[k]);
    Epetra_IntVector kgids_use(View, cell_map_use, gids[k]);
    kgids_use.Import(kgids_own, importer, Insert);
  }

  // Compare the ghost cell getEntityGIDs against the reference values just computed.
  for (AmanziMesh::Entity_ID j = ncell_own; j < ncell_use; ++j) {
    mesh->getCellNodes(j, &cnode);
    bool bad_data = false;
    for (int k = 0; k < cnode.size(); ++k)
      if (node_map.getEntityGID(cnode[k]) != gids(j,k)) bad_data = true;
    if (bad_data) {
      for (int k = 0; k < cnode.size(); ++k) {
        bool found = false;
        for (int l = 0; l < cnode.size(); ++l)
          if (node_map.getEntityGID(cnode[k]) == gids(j,l)) {
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
    os << "ERROR: found bad data for ghost cells:";
    write_list(bad_cells, MAX_OUT);
    os << "       The ghost cells don't have the same nodes as their respective masters." << std::endl;
    error = true;
  }

  return global_any(error);
}


// Check that ghost cells reference the same faces, and in the same
// order, as their master.  This is part of a ghost being an exact
// copies the master.  Even if the preceding ghost checks pass it is
// still possible for this one to fail.  In this case the ghost would
// reference a different face (by getEntityGID) but that face would reference
// the same nodes as the correct face.  So the two faces would be
// geometrically identical, including orientation, but be distinct.

bool MeshAudit::check_cell_to_faces_ghost_data() const
{
  const Epetra_Map &face_map = mesh->face_map(true);
  const Epetra_Map &cell_map_own = mesh->cell_map(false);
  const Epetra_Map &cell_map_use = mesh->cell_map(true);

  int ncell_own = cell_map_own.NumMyElements();
  int ncell_use = cell_map_use.NumMyElements();

  AmanziMesh::Entity_ID_List cface;
  AmanziMesh::Entity_ID_List bad_cells;

  // Create a matrix of the getEntityGIDs for all owned cells.
  int maxfaces = 0;
  for (AmanziMesh::Entity_ID j = 0; j < ncell_own; ++j) {
    mesh->getCellFaces(j, &cface);
    maxfaces = (cface.size() > maxfaces) ? cface.size() : maxfaces;
  }

  // parallel calls require uniformity
  int maxfaces_tmp(maxfaces);
  comm_->MaxAll(&maxfaces_tmp, &maxfaces, 1);

  Epetra_IntSerialDenseMatrix gids(ncell_use,maxfaces); // no Epetra_IntMultiVector :(

  for (AmanziMesh::Entity_ID j = 0; j < ncell_own; ++j) {
    mesh->getCellFaces(j, &cface);
    for (int k = 0; k < cface.size(); ++k)
      gids(j,k) = face_map.getEntityGID(cface[k]);
    for (int k = cface.size(); k < maxfaces; ++k)
      gids(j,k) = 0;
  }

  // Import these getEntityGIDs to all used cells; sets values on ghost cells.
  Epetra_Import importer(cell_map_use, cell_map_own);
  for (int k = 0; k < maxfaces; ++k) {
    Epetra_IntVector kgids_own(View, cell_map_own, gids[k]);
    Epetra_IntVector kgids_use(View, cell_map_use, gids[k]);
    kgids_use.Import(kgids_own, importer, Insert);
  }

  // Compare the ghost cell getEntityGIDs against the reference values just computed.
  for (AmanziMesh::Entity_ID j = ncell_own; j < ncell_use; ++j) {
    mesh->getCellFaces(j, &cface);
    bool bad_data = false;
    for (int k = 0; k < cface.size(); ++k)
      if (face_map.getEntityGID(cface[k]) != gids(j,k)) bad_data = true;
    if (bad_data) bad_cells.push_back(j);
  }

  bool error = false;

  if (!bad_cells.empty()) {
    os << "ERROR: found bad data for ghost cells:";
    write_list(bad_cells, MAX_OUT);
    os << "       The ghost cells are not exact copies of their master." << std::endl;
    error = true;
  }

  return global_any(error);
}

////////////////////////////////////////////////////////////////////////////////
//
// TESTS OF SET DATA
//
////////////////////////////////////////////////////////////////////////////////


// Check that get_set_ids successfully returns the vector of set IDs, without
// duplicates, and that each process gets the exact same vector of set IDs.
// This is a collective test, returning a collective pass/fail result.

bool MeshAudit::check_node_set_ids() const
{
  return check_get_set_ids(AmanziMesh::Entity_kind::NODE);
}

bool MeshAudit::check_face_set_ids() const
{
  return check_get_set_ids(AmanziMesh::Entity_kind::FACE);
}

bool MeshAudit::check_cell_set_ids() const
{
  return check_get_set_ids(AmanziMesh::Entity_kind::CELL);
}

bool MeshAudit::check_get_set_ids(AmanziMesh::Entity_kind kind) const
{
  os << "WARNING: Checks on sets disabled until MeshAudit handles new set specification methods (Tkt #686)" << std::endl;
  return false;

  bool error = false;

  // // Get the number of sets.
  // int nset;
  // try {
  //   nset = mesh->num_sets(kind); // this may fail
  // } catch (...) {
  //   os << "ERROR: caught exception from num_sets()" << std::endl;
  //   error = true;
  // }
  // error = global_any(error);
  // if (error) return error;

  // // Get the vector of set IDs.
  // AmanziMesh::Set_ID_List sids(nset, UINT_MAX);
  // try {
  //   mesh->get_set_ids(kind, &sids); // this may fail
  // } catch (...) {
  //   os << "ERROR: caught exception from get_set_ids()" << std::endl;
  //   error = true;
  // }
  // error = global_any(error);
  // if (error) return error;

  // // Check to see that set ID values were actually assigned.  This assumes
  // // UINT_MAX is not a valid set ID.  This is a little iffy; perhaps 0 should
  // // be declared as an invalid set ID instead (the case for ExodusII), or
  // // perhaps we should just skip this check.
  // bool bad_data = false;
  // for (int j = 0; j < nset; ++j)
  //   if (sids[j] == UINT_MAX) bad_data = true;
  // if (bad_data) {
  //   os << "ERROR: get_set_ids() failed to set all values" << std::endl;
  //   error = true;
  // }
  // error = global_any(error);
  // if (error) return error;

  // // Verify that the vector of set IDs contains no duplicates.
  // if (!distinct_values(sids)) {
  //   os << "ERROR: get_set_ids() returned duplicate IDs" << std::endl;
  //   // it would be nice to output the duplicates
  //   error = true;
  // }
  // error = global_any(error);
  // if (error) return error;

  // // In parallel, verify that each process returns the exact same result.
  // if (comm_->NumProc() > 1) {
  //   // Check the number of sets are the same.
  //   comm_->Broadcast(&nset, 1, 0);
  //   if (nset != mesh->num_sets(kind)) {
  //     os << "ERROR: inconsistent num_sets() value" << std::endl;
  //     error = true;
  //   }
  //   error = global_any(error);

  //   if (!error) {
  //     // Broadcast the set IDs on processor 0.
  //   AmanziMesh::Set_ID_List sids(nset, UINT_MAX);
  //   mesh->get_set_ids(kind, &sids);
  //     int *sids0 = new int[nset];
  //     for (int j = 0; j < nset; ++j) sids0[j] = sids[j];
  //     comm_->Broadcast(sids0, nset, 0);

  //     // Check the set IDs, using the vector on process 0 as the reference.
  //     bool bad_data = false;
  //     for (int j = 0; j < nset; ++j)
  //       if (sids[j] != sids0[j]) bad_data = true;
  //     if (bad_data) {
  //       os << "ERROR: get_set_ids() returned inconsistent values" << std::endl;
  //       error = true;
  //     }
  //     delete [] sids0;
  //   }
  // }

  return global_any(error);
}

// Check that valid_set_id() returns the correct results.
// This is a collective test, returning a collective pass/fail result.

bool MeshAudit::check_valid_node_set_id() const
{
  return check_valid_set_id(AmanziMesh::Entity_kind::NODE);
}

bool MeshAudit::check_valid_face_set_id() const
{
  return check_valid_set_id(AmanziMesh::Entity_kind::FACE);
}

bool MeshAudit::check_valid_cell_set_id() const
{
  return check_valid_set_id(AmanziMesh::Entity_kind::CELL);
}

bool MeshAudit::check_valid_set_id(AmanziMesh::Entity_kind kind) const
{
  os << "WARNING: Checks on sets disabled until MeshAudit handles new set specification methods (Tkt #686)" << std::endl;
  return false;

  // // Get the list of set IDs.
  // int nset = mesh->num_sets(kind); // this should not fail
  // AmanziMesh::Set_ID_List sids(nset);
  // mesh->get_set_ids(kind, &sids); // this should not fail

  // AmanziMesh::Set_ID_List bad_sids;
  // int max_id = 0;
  // for (int j = 0; j < nset; ++j)
  //   if (max_id < sids[j]) max_id = sids[j];
  // vector<bool> valid(max_id+2, false);
  // for (int j = 0; j < nset; ++j)
  //   valid[sids[j]] = true;

  // AmanziMesh::Set_ID_List bad_sids1, bad_sids2;
  // for (int n = 0; n < valid.size(); ++n) {
  //   if (valid[n] && !mesh->valid_set_id(n, kind)) bad_sids1.push_back(n);
  //   if (!valid[n] && mesh->valid_set_id(n, kind)) bad_sids2.push_back(n);
  // }

  bool error = false;

  // if (!bad_sids1.empty()) {
  //   os << "ERROR: valid_set_id() returned false for valid set IDs:";
  //   write_list(bad_sids1, MAX_OUT);
  //   error = true;
  // }

  // if (!bad_sids2.empty()) {
  //   os << "ERROR: valid_set_id() returned true for invalid set IDs:";
  //   write_list(bad_sids2, MAX_OUT);
  //   error = true;
  // }

  return global_any(error);
}

// For each set, check that get_set successfully returns valid references to
// local entities, without duplicates, and that the used set is consistent
// with the owned set.

bool MeshAudit::check_node_sets() const
{
  return check_sets(AmanziMesh::Entity_kind::NODE, mesh->node_map(false), mesh->node_map(true));
}

bool MeshAudit::check_face_sets() const
{
  return check_sets(AmanziMesh::Entity_kind::FACE, mesh->face_map(false), mesh->face_map(true));
}

bool MeshAudit::check_cell_sets() const
{
  return check_sets(AmanziMesh::Entity_kind::CELL, mesh->cell_map(false), mesh->cell_map(true));
}

bool MeshAudit::check_sets(AmanziMesh::Entity_kind kind,
                          const Epetra_Map &map_own, const Epetra_Map &map_use) const
{
  os << "WARNING: Checks on sets disabled until MeshAudit handles new set specification methods (Tkt #686)" << std::endl;
  return false;

  bool error = false;

  // // Get the list of set IDs.
  // int nset = mesh->num_sets(kind);
  // AmanziMesh::Set_ID_List sids(nset);
  // mesh->get_set_ids(kind, &sids);

  // for (int n = 0; n < sids.size(); ++n) {
  //   os << "  Checking set ID=" << sids[n] << " ..." << std::endl;

  //   // Basic sanity checks of the owned and used sets.
  //   bool bad_set = check_get_set(sids[n], kind, AmanziMesh::Parallel_type::OWNED, 
  //       			 map_own) ||
  //                  check_get_set(sids[n], kind, AmanziMesh::Parallel_type::ALL,  
  //       			 map_use);
  //   bad_set = global_any(bad_set);

  //   // Verify the used set relates correctly to the owned set.
  //   if (!bad_set) bad_set = check_used_set(sids[n], kind, map_own, map_use);

  //   // OUGHT TO DO TESTING OF THE GHOST SETS

  //   if (bad_set) error = true;
  // }

  return error;
}

// Basic sanity check on set values: no duplicates, and all LID values belong
// to the map.  This test runs independently on each process and returns a
// per-process pass/fail result.

bool MeshAudit::check_get_set(AmanziMesh::Set_ID sid, 
			      AmanziMesh::Entity_kind kind,
			      AmanziMesh::Parallel_type ptype, 
			      const Epetra_Map &map) const
{
  os << "WARNING: Checks on sets disabled until MeshAudit handles new set specification methods (Tkt #686)" << std::endl;
  return false;

  // Get the size of the set.
  try {
    std::string set_name = mesh->geometric_model()->FindRegion(sid)->name();
    mesh->get_set_size(set_name, kind, ptype); // this may fail
  } catch (...) {
    os << "  ERROR: caught exception from get_set_size()" << std::endl;
    return true;
  }

  // Get the set.
  AmanziMesh::Entity_ID_List set;
  try {
    std::string set_name = mesh->geometric_model()->FindRegion(sid)->name();
    mesh->getSetEntities(set_name, kind, ptype, &set);  // this may fail
  } catch (...) {
    os << "  ERROR: caught exception from get_set()" << std::endl;
    return true;
  }

  // Check that all values were assigned.
  bool bad_data = false;
  for (int j = 0; j < set.size(); ++j)
    if (set[j] == UINT_MAX) bad_data = true;
  if (bad_data) {
    os << "  ERROR: not all values assigned by get_set()" << std::endl;
    return true;
  }

  // Check that the LIDs in the set belong to the map.
  AmanziMesh::Entity_ID_List bad_LIDs;
  for (int j = 0; j < set.size(); ++j)
    if (!map.MyLID(set[j])) bad_LIDs.push_back(set[j]);
  if (!bad_LIDs.empty()) {
    os << "  ERROR: set contains invalid LIDs:";
    write_list(bad_LIDs, MAX_OUT);
    return true;
  }

  // Check that there are no duplicates in the set.
  if (!distinct_values(set)) {
    os << "  ERROR: set contains duplicate LIDs." << std::endl;
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

bool MeshAudit::check_used_set(AmanziMesh::Set_ID sid, 
			       AmanziMesh::Entity_kind kind,
                               const Epetra_Map &map_own, 
			       const Epetra_Map &map_use) const
{
  os << "WARNING: Checks on sets disabled until MeshAudit handles new set specification methods (Tkt #686)" << std::endl;
  return false;

  std::string set_name = mesh->geometric_model()->FindRegion(sid)->name();

  if (comm_->NumProc() == 1) {
    // In serial, the owned and used sets should be identical.

    int n = mesh->get_set_size(set_name, kind, AmanziMesh::Parallel_type::OWNED);
    AmanziMesh::Entity_ID_List set_own;
    mesh->getSetEntities(set_name, kind, AmanziMesh::Parallel_type::OWNED, &set_own);

    // Set sizes had better be the same.
    if (mesh->get_set_size(set_name, kind, AmanziMesh::Parallel_type::ALL) != 
	set_own.size()) {
      os << "  ERROR: owned and used set sizes differ" << std::endl;
      return true;
    }

    // Verify that the two sets are identical.
    AmanziMesh::Entity_ID_List set_use;
    mesh->getSetEntities(set_name, kind, AmanziMesh::Parallel_type::ALL, &set_use);
    bool bad_data = false;
    for (int j = 0; j < n; ++j)
      if (set_use[j] != set_own[j]) bad_data = true;
    if (bad_data) {
      os << "  ERROR: owned and used sets differ" << std::endl;
      return true;
    }

    return false;

  } else {

    int n = mesh->get_set_size(set_name, kind, AmanziMesh::Parallel_type::OWNED);
    AmanziMesh::Entity_ID_List set_own;
    mesh->getSetEntities(set_name, kind, AmanziMesh::Parallel_type::OWNED, &set_own);

    n = mesh->get_set_size(set_name, kind, AmanziMesh::Parallel_type::ALL);
    AmanziMesh::Entity_ID_List set_use(n);
    mesh->getSetEntities(set_name, kind, AmanziMesh::Parallel_type::ALL, &set_use);

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
    AmanziMesh::Entity_ID_List bad_LIDs;
    for (int j = 0; j < set_use.size(); ++j)
      if (tag_use[j] < 0) bad_LIDs.push_back(j);
    if (!bad_LIDs.empty()) {
      os << "  ERROR: found used LIDs that belong to the set but shouldn't:";
      write_list(bad_LIDs, MAX_OUT);
      error = true;
    }

    // Check for positive tag values;
    // these mark used LIDs that should be in the set but aren't.
    bad_LIDs.resize(0);
    for (int j = 0; j < set_own.size(); ++j)
      if (tag_use[j] > 0) bad_LIDs.push_back(j);
    if (!bad_LIDs.empty()) {
      os << "  ERROR: found used LIDs that should belong to set but don't:";
      write_list(bad_LIDs, MAX_OUT);
      error = true;
    }

    return global_any(error);
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

bool MeshAudit::check_face_partition() const
{
  // Mark all the faces contained by owned cells.
  bool owned[nface];
  for (int j = 0; j < nface; ++j) owned[j] = false;
  AmanziMesh::Entity_ID_List cface;
  for (AmanziMesh::Entity_ID j = 0; j < mesh->cell_map(false).NumMyElements(); ++j) {
    mesh->getCellFaces(j, &cface);
    for (int k = 0; k < cface.size(); ++k) owned[cface[k]] = true;
  }

  // Verify that every owned face has been marked as belonging to an owned cell.
  AmanziMesh::Entity_ID_List bad_faces;
  for (AmanziMesh::Entity_ID j = 0; j < mesh->face_map(false).NumMyElements(); ++j)
    if (!owned[j]) bad_faces.push_back(j);

  if (!bad_faces.empty()) {
    os << "ERROR: found orphaned owned faces:";
    write_list(bad_faces, MAX_OUT);
    os << "       Process doesn't own either of the cells sharing the face." << std::endl;
    return true;
  }

  return false;
}


bool MeshAudit::check_node_partition() const
{
  // Mark all the nodes contained by owned cells.
  bool owned[nnode];
  for (int j = 0; j < nnode; ++j) owned[j] = false;
  AmanziMesh::Entity_ID_List cnode;
  for (AmanziMesh::Entity_ID j = 0; j < mesh->cell_map(false).NumMyElements(); ++j) {
    mesh->getCellNodes(j, &cnode);
    for (int k = 0; k < cnode.size(); ++k) owned[cnode[k]] = true;
  }

  // Verify that every owned node has been marked as belonging to an owned cell.
  AmanziMesh::Entity_ID_List bad_nodes;
  for (AmanziMesh::Entity_ID j = 0; j < mesh->node_map(false).NumMyElements(); ++j)
    if (!owned[j]) bad_nodes.push_back(j);

  if (!bad_nodes.empty()) {
    os << "ERROR: found orphaned owned nodes:";
    write_list(bad_nodes, MAX_OUT);
    os << "       Process doesn't own any of the cells containing the node." << std::endl;
    return true;
  }

  return false;
}

#endif

// Returns true if the values in the list are distinct -- no repeats.

bool MeshAudit::distinct_values(const AmanziMesh::Entity_ID_List &list) const
{
  AmanziMesh::Entity_ID_List copy(list);
  sort(copy.begin(), copy.end());
  return (adjacent_find(copy.begin(),copy.end()) == copy.end());
}


// Returns 1 if the face node lists fnode1 and fnode2 describe the same face
// with the same orientation.  Returns -1 if the lists describe the same face
// but with opposite orientations.  Returns 0 if the lists describe different
// faces.  Implicitly assumes non-degenerate faces; the results are not
// reliable otherwise.

int MeshAudit::same_face(const AmanziMesh::Entity_ID_List fnode1, const AmanziMesh::Entity_ID_List fnode2) const
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

  }
  else {
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


void MeshAudit::write_list(const AmanziMesh::Entity_ID_List &list, unsigned int max_out) const
{
  int num_out = min((unsigned int) list.size(), max_out);
  for (int i = 0; i < num_out; ++i) os << " " << list[i];
  if (num_out < list.size()) os << " [" << list.size()-num_out << " items omitted]";
  os << std::endl;
}

bool MeshAudit::global_any(bool value) const
{
  int lval=value, gval;
  comm_->MaxAll(&lval, &gval, 1);
  return gval;
}

} // close namespace Amanzi

#include "MeshAuditOld.hh"

#include <algorithm>
#include <cfloat>

#include <boost/graph/topological_sort.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"

#include "Cell_topology.hh"
#include "cell_geometry.hh"

#include <iostream>
#include <iomanip>

using namespace std;
using namespace boost;
using namespace Amanzi;
using namespace AmanziMesh;

MeshAuditOld:: MeshAuditOld(Teuchos::RCP<Mesh> &mesh_, ostream& os_) :
      mesh(mesh_), comm(*(mesh_->get_comm())), MyPID(mesh_->get_comm()->MyPID()),
      os(os_),
      nnode(mesh_->node_map(true).NumMyElements()),
      nface(mesh_->face_map(true).NumMyElements()),
      ncell(mesh_->cell_map(true).NumMyElements()),
      MAX_OUT(5)
    { create_test_dependencies(); }

// Verify runs all the tests in an order that respects the dependencies
// between tests.  If a particular test fails, all other tests that have it
// as a pre-requisite are skipped.  It is important that each test return
// a collective fail/pass result in parallel, so that all processes proceed
// through the tests in lockstep.

int MeshAuditOld::Verify() const
{
  int status = 0;

  typedef Graph::vertex_descriptor Vertex;
  std::list<Vertex> run_order;
  topological_sort(g, std::front_inserter(run_order));

  mark_do_not_run vis;

  for (std::list<Vertex>::iterator itr = run_order.begin(); itr != run_order.end(); ++itr) {
    if (g[*itr].run){
      os << "Checking " << g[*itr].name << " ..." << endl;
      if (((*this).*(g[*itr].test))()) {
        status = 1;
        os << "  Test failed!" << endl;
        breadth_first_search(g, *itr, visitor(vis));
      }
    } else {
      os << "Skipping " << g[*itr].name << " check because of previous failures." << endl;
    }
  }

  return status;
}

// This creates the dependency graph for the individual tests.  Adding a new
// test to the graph is a simple matter of adding a new stanza of the form
//
//   Graph::vertex_descriptor my_test_handle = add_vertex(g);
//   g[my_test_handle].name = "description of my test";
//   g[my_test_handle].test = &MeshAuditOld::my_test;
//   add_edge(other_test_handle, my_test_handle, g);
//
// The last line specifies that other_test_handle is a pre-requisite for
// my_test_handle.  There may be multiple pre-requisites or none.

void MeshAuditOld::create_test_dependencies()
{
  // Entity_counts tests
  Graph::vertex_descriptor test01 = add_vertex(g);
  g[test01].name = "entity_counts";
  g[test01].test = &MeshAuditOld::check_entity_counts;

  // Cell_to_nodes tests
  Graph::vertex_descriptor test02 = add_vertex(g);
  g[test02].name = "cell_to_nodes";
  g[test02].test = &MeshAuditOld::check_cell_to_nodes;

  Graph::vertex_descriptor test03 = add_vertex(g);
  g[test03].name = "consistency of cell_to_nodes methods";
  g[test03].test = &MeshAuditOld::check_cell_to_nodes_consistency;
  add_edge(test02, test03, g);

  Graph::vertex_descriptor test04 = add_vertex(g);
  g[test04].name = "node references by cells";
  g[test04].test = &MeshAuditOld::check_node_refs_by_cells;
  add_edge(test02, test04, g);

  // Cell_to_faces tests
  Graph::vertex_descriptor test05 = add_vertex(g);
  g[test05].name = "cell_to_faces";
  g[test05].test = &MeshAuditOld::check_cell_to_faces;

  Graph::vertex_descriptor test06 = add_vertex(g);
  g[test06].name = "consistency of cell_to_faces methods";
  g[test06].test = &MeshAuditOld::check_cell_to_faces_consistency;
  add_edge(test05, test06, g);

  Graph::vertex_descriptor test07 = add_vertex(g);
  g[test07].name = "face references by cells";
  g[test07].test = &MeshAuditOld::check_face_refs_by_cells;
  add_edge(test05, test07, g);

  // face_to_nodes tests
  Graph::vertex_descriptor test08 = add_vertex(g);
  g[test08].name = "face_to_nodes";
  g[test08].test = &MeshAuditOld::check_face_to_nodes;

  Graph::vertex_descriptor test09 = add_vertex(g);
  g[test09].name = "consistency of face_to_nodes methods";
  g[test09].test = &MeshAuditOld::check_face_to_nodes_consistency;
  add_edge(test08, test09, g);

  Graph::vertex_descriptor test10 = add_vertex(g);
  g[test10].name = "node references by faces";
  g[test10].test = &MeshAuditOld::check_node_refs_by_faces;
  add_edge(test08, test10, g);

  // cell_to_face_dirs tests
  Graph::vertex_descriptor test11 = add_vertex(g);
  g[test11].name = "cell_to_face_dirs";
  g[test11].test = &MeshAuditOld::check_cell_to_face_dirs;

  Graph::vertex_descriptor test12 = add_vertex(g);
  g[test12].name = "consistency of cell_to_face_dirs methods";
  g[test12].test = &MeshAuditOld::check_cell_to_face_dirs_consistency;
  add_edge(test11, test12, g);

  // cell degeneracy test
  Graph::vertex_descriptor test13 = add_vertex(g);
  g[test13].name = "topological non-degeneracy of cells";
  g[test13].test = &MeshAuditOld::check_cell_degeneracy;
  add_edge(test02, test13, g);

  // Consistency between the various mesh connectivity data.
  Graph::vertex_descriptor test14 = add_vertex(g);
  g[test14].name = "consistency of mesh connectivity data";
  g[test14].test = &MeshAuditOld::check_cell_to_faces_to_nodes;
  add_edge(test02, test14, g);
  add_edge(test05, test14, g);
  add_edge(test08, test14, g);
  add_edge(test11, test14, g);
  add_edge(test13, test14, g);

  // node_to_coordinates tests
  Graph::vertex_descriptor test15 = add_vertex(g);
  g[test15].name = "node_to_coordinates";
  g[test15].test = &MeshAuditOld::check_node_to_coordinates;

  Graph::vertex_descriptor test16 = add_vertex(g);
  g[test16].name = "consistency of node_to_coordinates methods";
  g[test16].test = &MeshAuditOld::check_node_to_coordinates_alt;
  add_edge(test15, test16, g);

  // cell_to_coordinates tests
  Graph::vertex_descriptor test17 = add_vertex(g);
  g[test17].name = "cell_to_coordinates";
  g[test17].test = &MeshAuditOld::check_cell_to_coordinates;
  add_edge(test02, test17, g);
  add_edge(test15, test17, g);

  Graph::vertex_descriptor test18 = add_vertex(g);
  g[test18].name = "consistency of cell_to_coordinates methods";
  g[test18].test = &MeshAuditOld::check_cell_to_coordinates_alt;
  add_edge(test17, test18, g);

  // face_to_coordinates tests
  Graph::vertex_descriptor test19 = add_vertex(g);
  g[test19].name = "face_to_coordinates";
  g[test19].test = &MeshAuditOld::check_face_to_coordinates;
  add_edge(test08, test19, g);
  add_edge(test15, test19, g);

  Graph::vertex_descriptor test20 = add_vertex(g);
  g[test20].name = "consistency of face_to_coordinates methods";
  g[test20].test = &MeshAuditOld::check_face_to_coordinates_alt;
  add_edge(test19, test20, g);

  // cell topology/geometry test
  Graph::vertex_descriptor test21 = add_vertex(g);
  g[test21].name = "cell topology/geometry";
  g[test21].test = &MeshAuditOld::check_cell_topology;
  add_edge(test13, test21, g);
  add_edge(test18, test21, g);

  // map tests
  Graph::vertex_descriptor test22 = add_vertex(g);
  g[test22].name = "owned and overlap node maps";
  g[test22].test = &MeshAuditOld::check_node_maps;

  Graph::vertex_descriptor test23 = add_vertex(g);
  g[test23].name = "owned and overlap face maps";
  g[test23].test = &MeshAuditOld::check_face_maps;

  Graph::vertex_descriptor test24 = add_vertex(g);
  g[test24].name = "owned and overlap cell maps";
  g[test24].test = &MeshAuditOld::check_cell_maps;

  // ghost data tests
  Graph::vertex_descriptor test25 = add_vertex(g);
  g[test25].name = "node_to_coordinates ghost data";
  g[test25].test = &MeshAuditOld::check_node_to_coordinates_ghost_data;
  add_edge(test15, test25, g);
  add_edge(test22, test25, g);

  Graph::vertex_descriptor test26 = add_vertex(g);
  g[test26].name = "face_to_nodes ghost data";
  g[test26].test = &MeshAuditOld::check_face_to_nodes_ghost_data;
  add_edge(test08, test26, g);
  add_edge(test22, test26, g);
  add_edge(test23, test26, g);

  Graph::vertex_descriptor test27 = add_vertex(g);
  g[test27].name = "cell_to_nodes ghost data";
  g[test27].test = &MeshAuditOld::check_cell_to_nodes_ghost_data;
  add_edge(test02, test27, g);
  add_edge(test22, test27, g);
  add_edge(test24, test27, g);

  Graph::vertex_descriptor test28 = add_vertex(g);
  g[test28].name = "cell_to_faces ghost data";
  g[test28].test = &MeshAuditOld::check_cell_to_faces_ghost_data;
  add_edge(test05, test28, g);
  add_edge(test23, test28, g);
  add_edge(test24, test28, g);

  // node set data tests
  Graph::vertex_descriptor test29 = add_vertex(g);
  g[test29].name = "node set IDs";
  g[test29].test = &MeshAuditOld::check_node_set_ids;

  Graph::vertex_descriptor test30 = add_vertex(g);
  g[test30].name = "node sets";
  g[test30].test = &MeshAuditOld::check_node_sets;
  add_edge(test22, test30, g);
  add_edge(test29, test30, g);

  Graph::vertex_descriptor test31 = add_vertex(g);
  g[test31].name = "valid node set IDs";
  g[test31].test = &MeshAuditOld::check_valid_node_set_id;
  add_edge(test29, test31, g);

  Graph::vertex_descriptor test32 = add_vertex(g);
  g[test32].name = "alternative node set methods";
  g[test32].test = &MeshAuditOld::check_node_sets_alt;
  add_edge(test29, test32, g);
  add_edge(test30, test32, g);

  // face set data tests
  Graph::vertex_descriptor test33 = add_vertex(g);
  g[test33].name = "face set IDs";
  g[test33].test = &MeshAuditOld::check_face_set_ids;

  Graph::vertex_descriptor test34 = add_vertex(g);
  g[test34].name = "face sets";
  g[test34].test = &MeshAuditOld::check_face_sets;
  add_edge(test23, test34, g);
  add_edge(test33, test34, g);

  Graph::vertex_descriptor test35 = add_vertex(g);
  g[test35].name = "valid face set IDs";
  g[test35].test = &MeshAuditOld::check_valid_face_set_id;
  add_edge(test33, test35, g);

  Graph::vertex_descriptor test36 = add_vertex(g);
  g[test36].name = "alternative face set methods";
  g[test36].test = &MeshAuditOld::check_face_sets_alt;
  add_edge(test33, test36, g);
  add_edge(test34, test36, g);

  // cell set data tests
  Graph::vertex_descriptor test37 = add_vertex(g);
  g[test37].name = "cell set IDs";
  g[test37].test = &MeshAuditOld::check_cell_set_ids;

  Graph::vertex_descriptor test38 = add_vertex(g);
  g[test38].name = "cell sets";
  g[test38].test = &MeshAuditOld::check_cell_sets;
  add_edge(test24, test38, g);
  add_edge(test37, test38, g);

  Graph::vertex_descriptor test39 = add_vertex(g);
  g[test39].name = "valid cell set IDs";
  g[test39].test = &MeshAuditOld::check_valid_cell_set_id;
  add_edge(test37, test39, g);

  Graph::vertex_descriptor test40 = add_vertex(g);
  g[test40].name = "alternative cell set methods";
  g[test40].test = &MeshAuditOld::check_cell_sets_alt;
  add_edge(test37, test40, g);
  add_edge(test38, test40, g);

  // partition tests
  Graph::vertex_descriptor test41 = add_vertex(g);
  g[test41].name = "face partition";
  g[test41].test = &MeshAuditOld::check_face_partition;
  add_edge(test23, test41, g);
  add_edge(test24, test41, g);

  Graph::vertex_descriptor test42 = add_vertex(g);
  g[test42].name = "node partition";
  g[test42].test = &MeshAuditOld::check_node_partition;
  add_edge(test22, test42, g);
  add_edge(test24, test42, g);
}

////////////////////////////////////////////////////////////////////////////////
//
// Tests (and auxillary functions) follow.  Tests must be const functions
// that take no arguments and that return a bool result: true if an error
// was found, otherwise false.  Tests must be considered collective procedures
// when in a parallel context, and so return a common collective result.
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

bool MeshAuditOld::check_entity_counts() const
{
  int n, nref;
  bool error = false;

  // Check the number of owned nodes.
  n = mesh->count_entities(NODE,OWNED);
  nref = mesh->node_map(false).NumMyElements();
  if (n != nref) {
    os << ": ERROR: count_entities(NODE,OWNED)=" << n << "; should be " << nref << endl;
    error = true;
  }

  // Check the number of used nodes.
  n = mesh->count_entities(NODE,USED);
  nref = mesh->node_map(true).NumMyElements();
  if (n != nref) {
    os << "ERROR: count_entities(NODE,USED)=" << n << "; should be " << nref << endl;
    error = true;
  }

  // Check the number of owned faces.
  n = mesh->count_entities(FACE,OWNED);
  nref = mesh->face_map(false).NumMyElements();
  if (n != nref) {
    os << "ERROR: count_entities(FACE,OWNED)=" << n << "; should be " << nref << endl;
    error = true;
  }

  // Check the number of used faces.
  n = mesh->count_entities(FACE,USED);
  nref = mesh->face_map(true).NumMyElements();
  if (n != nref) {
    os << "ERROR: count_entities(FACE,USED)=" << n << "; should be " << nref << endl;
    error = true;
  }

  // Check the number of owned cells.
  n = mesh->count_entities(CELL,OWNED);
  nref = mesh->cell_map(false).NumMyElements();
  if (n != nref) {
    os << "ERROR: count_entities(CELL,OWNED)=" << n << "; should be " << nref << endl;
    error = true;
  }

  // Check the number of used cells.
  n = mesh->count_entities(CELL,USED);
  nref = mesh->cell_map(true).NumMyElements();
  if (n != nref) {
    os << "ERROR: count_entities(CELL,USED)=" << n << "; should be " << nref << endl;
    error = true;
  }

  return global_any(error);
}

// Check that cell_to_nodes successfully returns valid references to local
// nodes.  Here the std::vector accessor is used.  The consistency of
// the alternative accessors is checked elsewhere.  A nonzero return value
// signals an error, and further tests using its data should be avoided.

bool MeshAuditOld::check_cell_to_nodes() const
{
  vector<unsigned int> bad_cells, bad_cells1;
  vector<unsigned int> cnode(8);

  for (unsigned int j = 0; j < ncell; ++j) {
    cnode.assign(8, UINT_MAX);
    try {
      mesh->cell_to_nodes(j, cnode.begin(), cnode.end()); // this may fail
      bool invalid_refs = false;
      for (int k = 0; k < cnode.size(); ++k) {
        if (cnode[k] < 0 || cnode[k] >= nnode) invalid_refs = true;
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

bool MeshAuditOld::check_node_refs_by_cells() const
{
  vector<unsigned int> cnode(8);
  vector<bool> ref(nnode, false);

  for (unsigned int j = 0; j < ncell; ++j) {
    mesh->cell_to_nodes(j, cnode.begin(), cnode.end()); // this should not fail
    for (int k = 0; k < cnode.size(); ++k) ref[cnode[k]] = true;
  }

  vector<unsigned int> free_nodes;
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

bool MeshAuditOld::check_cell_to_faces() const
{
  vector<unsigned int> bad_cells, bad_cells1;
  vector<unsigned int> cface(6);

  for (unsigned int j = 0; j < ncell; ++j) {
    cface.assign(6, UINT_MAX);
    try {
      mesh->cell_to_faces(j, cface.begin(), cface.end()); // this may fail
      bool invalid_refs = false;
      for (int k = 0; k < cface.size(); ++k) {
        if (cface[k] < 0 || cface[k] >= nface) invalid_refs = true;
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

bool MeshAuditOld::check_face_refs_by_cells() const
{
  vector<unsigned int> cface(6);
  vector<unsigned int> refs(nface, 0);

  for (unsigned int j = 0; j < ncell; ++j) {
    mesh->cell_to_faces(j, cface.begin(), cface.end());
    for (int k = 0; k < cface.size(); ++k) (refs[cface[k]])++;
  }

  vector<unsigned int> free_faces;
  vector<unsigned int> bad_faces;

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

  if (!bad_faces.empty()) {
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

bool MeshAuditOld::check_face_to_nodes() const
{
  vector<unsigned int> bad_faces, bad_faces1;
  vector<unsigned int> fnode(4);

  for (unsigned int j = 0; j < nface; ++j) {
    fnode.assign(4, UINT_MAX);
    try {
      mesh->face_to_nodes(j, fnode.begin(), fnode.end()); // this may fail
      bool invalid_refs = false;
      for (int k = 0; k < fnode.size(); ++k) {
        if (fnode[k] < 0 || fnode[k] >= nnode) invalid_refs = true;
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

bool MeshAuditOld::check_node_refs_by_faces() const
{
  vector<unsigned int> fnode(4);
  vector<bool> ref(nnode, false);

  for (unsigned int j = 0; j < nface; ++j) {
    mesh->face_to_nodes(j, fnode.begin(), fnode.end());
    for (int k = 0; k < fnode.size(); ++k) ref[fnode[k]] = true;
  }

  vector<unsigned int> free_nodes;
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

bool MeshAuditOld::check_cell_to_face_dirs() const
{
  vector<int> fdirs(6);
  vector<unsigned int> bad_cells, bad_cells_exc;

  for (unsigned int j = 0; j < ncell; ++j) {
    fdirs.assign(6, INT_MAX);
    try {
      mesh->cell_to_face_dirs(j, fdirs.begin(), fdirs.end());  // this may fail
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

// The mesh interface provides alternative methods for accessing
// connectivity data.  These alternates should return the same data.
// Here the accessors using std::vector are considered normative,
// and it is assumed that they have been verified to return valid
// (though perhaps not correct) results.  A positive value is
// returned if any inconsistency is found.  Further tests, which
// use only the std::vector-based accessors, are safe.

bool MeshAuditOld::check_cell_to_nodes_consistency() const
{
  vector<unsigned int> cnode(8);
  vector<unsigned int> bad_cells1, bad_cells1_exc;

  for (unsigned int j = 0; j < ncell; ++j) {
    mesh->cell_to_nodes(j, cnode.begin(), cnode.end()); // this should not fail
    try {
      unsigned int cnode1[8] = { UINT_MAX };
      mesh->cell_to_nodes(j, cnode1, cnode1+8); // this may fail
      bool bad_data = false;
      for (int i = 0; i < 8; ++i)
        if (cnode1[i] != cnode[i]) bad_data = true;
      if (bad_data) bad_cells1.push_back(j);
    } catch (...) {
      bad_cells1_exc.push_back(j);
    }
  }

  bool error = false;

  if (!bad_cells1.empty()) {
    os << "ERROR: bad values from pointer-based accessor for cells:";
    write_list(bad_cells1, MAX_OUT);
    error = true;
  }

  if (!bad_cells1_exc.empty()) {
    os << "ERROR: caught exception from pointer-based accessor for cells:";
    write_list(bad_cells1_exc, MAX_OUT);
    error = true;
  }

  return global_any(error);
}


bool MeshAuditOld::check_face_to_nodes_consistency() const
{
  vector<unsigned int> fnode(4);
  vector<unsigned int> bad_faces1, bad_faces1_exc;

  for (unsigned int j = 0; j < nface; ++j) {
    mesh->face_to_nodes(j, fnode.begin(), fnode.end()); // this should not fail
    try {
      unsigned int fnode1[4] = { UINT_MAX };
      mesh->face_to_nodes(j, fnode1, fnode1+4); // this may fail
      bool bad_data = false;
      for (int i = 0; i < 4; ++i)
        if (fnode1[i] != fnode[i]) bad_data = true;
      if (bad_data) bad_faces1.push_back(j);
    } catch (...) {
      bad_faces1_exc.push_back(j);
    }
  }

  bool error = false;

  if (!bad_faces1.empty()) {
    os << "ERROR: bad values from pointer-based accessor for faces:";
    write_list(bad_faces1, MAX_OUT);
    error = true;
  }

  if (!bad_faces1_exc.empty()) {
    os << "ERROR: caught exception from pointer-based accessor for faces:";
    write_list(bad_faces1_exc, MAX_OUT);
    error = true;
  }

  return global_any(error);
}


bool MeshAuditOld::check_cell_to_faces_consistency() const
{
  vector<unsigned int> cface(6);
  vector<unsigned int> bad_cells1, bad_cells1_exc;

  for (unsigned int j = 0; j < ncell; ++j) {
    mesh->cell_to_faces(j, cface.begin(), cface.end()); // this should not fail
    try {
      unsigned int cface1[6] = { UINT_MAX };
      mesh->cell_to_faces(j, cface1, cface1+6); // this may fail
      bool bad_data = false;
      for (int i = 0; i < 6; ++i)
        if (cface1[i] != cface[i]) bad_data = true;
      if (bad_data) bad_cells1.push_back(j);
    } catch (...) {
      bad_cells1_exc.push_back(j);
    }
  }

  bool error = false;

  if (!bad_cells1.empty()) {
    os << "ERROR: bad values from pointer-based accessor for cells:";
    write_list(bad_cells1, MAX_OUT);
    error = true;
  }

  if (!bad_cells1_exc.empty()) {
    os << "ERROR: caught exception from pointer-based accessor for cells:";
    write_list(bad_cells1_exc, MAX_OUT);
    error = true;
  }

  return global_any(error);
}


bool MeshAuditOld::check_cell_to_face_dirs_consistency() const
{
  vector<int> fdirs(6);
  vector<unsigned int> bad_cells1, bad_cells1_exc;

  for (unsigned int j = 0; j < ncell; ++j) {
    mesh->cell_to_face_dirs(j, fdirs.begin(), fdirs.end());  // this should not fail
    try {
      int fdirs1[6] = { INT_MAX };
      mesh->cell_to_face_dirs(j, fdirs1, fdirs1+6); // this may fail
      bool bad_data = false;
      for (int k = 0; k < 6; ++k)
        if (fdirs1[k] != fdirs[k]) bad_data = true;
      if (bad_data) bad_cells1.push_back(j);
    } catch (...) {
      bad_cells1_exc.push_back(j);
    }
  }

  bool error = false;

  if (!bad_cells1.empty()) {
    os << "ERROR: pointer-based accessor returned inconsistent values for cells:";
    write_list(bad_cells1, MAX_OUT);
    error = true;
  }

  if (!bad_cells1_exc.empty()) {
    os << "ERROR: caught exception from pointer-based accessor for cells:";
    write_list(bad_cells1_exc, MAX_OUT);
    error = true;
  }

  return global_any(error);
}

// Check that cells are not topologically degenerate (repeated node index).
// If this test fails, many further tests involving the cell_to_nodes data
// should be avoided.

bool MeshAuditOld::check_cell_degeneracy() const
{
  os << "Checking cells for topological degeneracy ..." << endl;

  vector<unsigned int> cnode(8);
  vector<unsigned int> bad_cells;

  for (unsigned int j = 0; j < ncell; ++j) {
    mesh->cell_to_nodes(j, cnode.begin(), cnode.end()); // should not fail
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
// using there values should be avoided.

bool MeshAuditOld::check_cell_to_faces_to_nodes() const
{
  vector<unsigned int> cnode(8);
  vector<unsigned int> cface(6);
  vector<unsigned int> fnode_ref(4);
  vector<unsigned int> fnode(4);
  vector<int> fdirs(6);
  vector<unsigned int> bad_cells0;
  vector<unsigned int> bad_cells1;

  for (unsigned int j = 0; j < ncell; ++j) {
    mesh->cell_to_nodes(j, cnode.begin(), cnode.end()); // this should not fail
    mesh->cell_to_faces(j, cface.begin(), cface.end()); // this should not fail
    mesh->cell_to_face_dirs(j, fdirs.begin(), fdirs.end()); // this should not fail
    bool bad_face = false;
    bool bad_dir  = false;
    for (int k = 0; k < cface.size(); ++k) {
      for (int i = 0; i < fnode_ref.size(); ++i)
        fnode_ref[i] = cnode[Amanzi::AmanziMesh::HexFaceVert[k][i]];
      mesh->face_to_nodes(cface[k], fnode.begin(), fnode.end()); // this should not fail
      int dir = same_face(fnode, fnode_ref); // should be the same face
      if (dir == 0) // wrong face
        bad_face = true;
      else if (dir != fdirs[k]) // right face but wrong dir value
        bad_dir = true;
    }
    if (bad_face) bad_cells0.push_back(j);
    if (bad_dir)  bad_cells1.push_back(j);
  }

  bool error = false;

  if (!bad_cells0.empty()) {
    os << "ERROR: bad cell_to_faces values for cells:";
    write_list(bad_cells0, MAX_OUT);
    error = true;
  }

  if (!bad_cells1.empty()) {
    os << "ERROR: bad cell_to_face_dirs values for cells:";
    write_list(bad_cells1, MAX_OUT);
    error = true;
  }

  return global_any(error);
}

// Check that node_to_coordinates successfully returns data for all nodes.
// The std::vector-based accessor is used here.  The consistency of the
// alternative accessors is checked elsewhere.  If this test fails, further
// tests using its data should be avoided.

bool MeshAuditOld::check_node_to_coordinates() const
{
  vector<double> x(3);
  vector<unsigned int> bad_nodes, bad_nodes_exc;

  for (unsigned int j = 0; j < nnode; ++j) {
    try {
      x.assign(3,DBL_MAX);
      mesh->node_to_coordinates(j, x.begin(), x.end()); // this may fail
      bool bad_data = false;
      for (int k = 0; k < x.size(); ++k)
        if (x[k] == DBL_MAX) bad_data = true;
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

// Check that the alternative node_to_coordinates accessors successfully
// return the same data as returned by the std::vector accessor, which is
// considered normative.  If this test fails, further tests, which use the
// std::vector-based accessor, are safe.

bool MeshAuditOld::check_node_to_coordinates_alt() const
{
  vector<double> x(3);
  vector<unsigned int> bad_nodes1, bad_nodes1_exc;

  for (unsigned int j = 0; j < nnode; ++j) {
    mesh->node_to_coordinates(j, x.begin(), x.end()); // this should not fail
    try {
      double x1[3] = { DBL_MAX };
      mesh->node_to_coordinates(j, x1, x1+3); // this may fail
      bool bad_data = false;
      for (int k = 0; k < 3; ++k)
        if (x1[k] != x[k]) bad_data = true;
      if (bad_data) bad_nodes1.push_back(j);
    } catch (...) {
      bad_nodes1_exc.push_back(j);
    }
  }

  bool error = false;

  if (!bad_nodes1.empty()) {
    os << "ERROR: pointer-based accessor returned inconsistent values for nodes:";
    write_list(bad_nodes1, MAX_OUT);
    error = true;
  }

  if (!bad_nodes1_exc.empty()) {
    os << "ERROR: caught exception from pointer-based accessor for nodes:";
    write_list(bad_nodes1_exc, MAX_OUT);
    error = true;
  }

  return global_any(error);
}

// Check that cell_to_coordinates successfully returns data for all cells, and
// that this data is identical to that returned by composing node_to_coordinates
// with cell_to_nodes.  The consistency of the alternative accessors is checked
// elsewhere.  If this test fails, further tests using cell_to_coordinates data
// should be avoided.

bool MeshAuditOld::check_cell_to_coordinates() const
{
  vector<double> xref(24), x(24); // 3x8
  vector<unsigned int> cnode(8);
  vector<unsigned int> bad_cells, bad_cells_exc;

  for (unsigned int j = 0; j < ncell; ++j) {
    mesh->cell_to_nodes(j, cnode.begin(), cnode.end()); // this should not fail
    vector<double>::iterator xbeg = xref.begin();
    for (int k = 0; k < 8; ++k) {
      mesh->node_to_coordinates(cnode[k], xbeg, xbeg+3); // this should not fail
      xbeg += 3;
    }
    try {
      x.assign(24, DBL_MAX);
      mesh->cell_to_coordinates(j, x.begin(), x.end()); // this may fail
      bool bad_data = false;
      for (int k = 0; k < 24; ++k)
        if (x[k] != xref[k]) bad_data = true;
      if (bad_data) bad_cells.push_back(j);
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

// Check that the alternative cell_to_coordinates accessors successfully
// return the same data as returned by the std::vector accessor, which is
// considered normative.  If this test fails, further tests, which use the
// std::vector-based accessor, are safe.

bool MeshAuditOld::check_cell_to_coordinates_alt() const
{
  vector<double> x(24); // 3x8
  vector<unsigned int> bad_cells1, bad_cells1_exc;

  for (unsigned int j = 0; j < ncell; ++j) {
    mesh->cell_to_coordinates(j, x.begin(), x.end()); // this should not fail
    try {
      double x1[24] = { DBL_MAX };
      mesh->cell_to_coordinates(j, x1, x1+24); // this may fail
      bool bad_data = false;
      for (int k = 0; k < 24; ++k)
        if (x1[k] != x[k]) bad_data = true;
      if (bad_data) bad_cells1.push_back(j);
    } catch (...) {
      bad_cells1_exc.push_back(j);
    }
  }

  bool error = false;

  if (!bad_cells1.empty()) {
    os << "ERROR: bad values from pointer-based accessor for cells:";
    write_list(bad_cells1, MAX_OUT);
    error = true;
  }

  if (!bad_cells1_exc.empty()) {
    os << "ERROR: caught exception from pointer-based accessor for cells:";
    write_list(bad_cells1_exc, MAX_OUT);
    error = true;
  }

  return global_any(error);
}

// Check that face_to_coordinates successfully returns data for all faces, and
// that this data is identical to that returned by composing node_to_coordinates
// with face_to_nodes.  The consistency of the alternative accessors is checked
// elsewhere.  If this test fails, further tests using face_to_coordinates data
// should be avoided.

bool MeshAuditOld::check_face_to_coordinates() const
{
  vector<double> xref(12), x(12); // 3x4
  vector<unsigned int> fnode(4);
  vector<unsigned int> bad_faces, bad_faces_exc;

  for (unsigned int j = 0; j < nface; ++j) {
    mesh->face_to_nodes(j, fnode.begin(), fnode.end()); // this should not fail
    vector<double>::iterator xbeg = xref.begin();
    for (int k = 0; k < 4; ++k) {
      mesh->node_to_coordinates(fnode[k], xbeg, xbeg+3); // this should not fail
      xbeg += 3;
    }
    try {
      x.assign(12, DBL_MAX);
      mesh->face_to_coordinates(j, x.begin(), x.end()); // this may fail
      bool bad_data = false;
      for (int k = 0; k < x.size(); ++k)
        if (x[k] != xref[k]) bad_data = true;
      if (bad_data) bad_faces.push_back(j);
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

// Check that the alternative face_to_coordinates accessors successfully
// return the same data as returned by the std::vector accessor, which is
// considered normative.  If this test fails, further tests, which use the
// std::vector-based accessor, are safe.

bool MeshAuditOld::check_face_to_coordinates_alt() const
{
  vector<double> x(12); // 3x4
  vector<unsigned int> bad_faces1, bad_faces1_exc;

  for (unsigned int j = 0; j < nface; ++j) {
    mesh->face_to_coordinates(j, x.begin(), x.end()); // this should not fail
    try {
      double x1[12] = { DBL_MAX };
      mesh->face_to_coordinates(j, x1, x1+12); // this may fail
      bool bad_data = false;
      for (int k = 0; k < 12; ++k)
        if (x1[k] != x[k]) bad_data = true;
      if (bad_data) bad_faces1.push_back(j);
    } catch (...) {
      bad_faces1_exc.push_back(j);
    }
  }

  bool error = false;

  if (!bad_faces1.empty()) {
    os << "ERROR: bad values from pointer-based accessor for faces:";
    write_list(bad_faces1, MAX_OUT);
    error = true;
  }

  if (!bad_faces1_exc.empty()) {
    os << "ERROR: caught exception from pointer-based accessor for faces:";
    write_list(bad_faces1_exc, MAX_OUT);
    error = true;
  }

  return global_any(error);
}

// The hexahedral cells must not be degenerate, either topologically (repeated
// node index) or geometrically (with coincident nodes).  In addition, the
// cell vertices must be ordered in a pre-defined manner (see cell_topology).
// To detect geometric degeneracy or bad topology (vertices listed in the
// wrong order), the corner tet volumes of the hexahedron are evaluated and
// checked for positivity.  If any are negative this indicates bad topology
// (is this sufficient?  I think so.)  Geometric degeneracy is indicated
// by one or more zero corner volumes.

bool MeshAuditOld::check_cell_topology() const
{
  double x[24], hvol, cvol[8];
  Epetra_SerialDenseMatrix xmat(View, x, 3, 3, 8);
  vector<unsigned int> bad_cells;

  for (unsigned int j = 0; j < ncell; ++j) {
    mesh->cell_to_coordinates(j, x, x+24); // should not fail
    cell_geometry::compute_hex_volumes(xmat, hvol, cvol);
    bool bad_vol = (hvol <= 0.0);
    for (int k = 0; k < 8; ++k)
      if (cvol[k] <= 0.0) bad_vol = true;
    if (bad_vol) bad_cells.push_back(j);
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

bool MeshAuditOld::check_node_maps() const
{
  return check_maps(mesh->node_map(false), mesh->node_map(true));
}

bool MeshAuditOld::check_face_maps() const
{
  return check_maps(mesh->face_map(false), mesh->face_map(true));
}

bool MeshAuditOld::check_cell_maps() const
{
  return check_maps(mesh->cell_map(false), mesh->cell_map(true));
}

bool MeshAuditOld::check_maps(const Epetra_Map &map_own, const Epetra_Map &map_use) const
{
  bool error = false;

  // Local index should start at 0.
  if (map_own.IndexBase() != 0) {
    os << "ERROR: the owned map's index base is not 0." << endl;
    error = true;
  }
  if (map_use.IndexBase() != 0) {
    os << "ERROR: the overlap map's index base is not 0." << endl;
    error = true;
  }

  // Check that the owned map is 1-1.
  if (!map_own.UniqueGIDs()) {
    os << "ERROR: owned map is not 1-to-1" << endl;
    error = true;
  }

  error = global_any(error);
  if (error) return error;

  if (comm.NumProc() == 1)
  {

    // Serial or 1-process MPI

    if (!map_use.SameAs(map_own)) {
      os << "ERROR: the overlap map differs from the owned map (single process)." << endl;
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
        if (map_use.GID(j) != map_own.GID(j)) bad_map = true;
    }
    if (bad_map) {
      os << "ERROR: overlap map does not extend the owned map." << endl;
      error = true;
    }

    error = global_any(error);
    if (error) return error;

    // Verify that the overlap indices are owned by other processes.
    int *gids = new int[num_ovl];
    int *pids = new int[num_ovl];
    int *lids = new int[num_ovl];
    for (int j = 0; j < num_ovl; ++j) gids[j] = map_use.GID(j+num_own);
    map_own.RemoteIDList(num_ovl, gids, pids, lids);
    bad_map = false;
    for (int j = 0; j < num_ovl; ++j)
      if (pids[j] < 0 || pids[j] == comm.MyPID()) bad_map = true;
    if (bad_map) {
      os << "ERROR: invalid ghosts in overlap map." << endl;
      error = true;
    }

    // Look for duplicates among the overlap indices.
    vector<int> ovl_gids(gids, gids+num_ovl);
    sort(ovl_gids.begin(), ovl_gids.end());
    if (adjacent_find(ovl_gids.begin(),ovl_gids.end()) != ovl_gids.end()) {
      os << "ERROR: duplicate ghosts in overlap map." << endl;
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

bool MeshAuditOld::check_node_to_coordinates_ghost_data() const
{
  const Epetra_Map &node_map_own = mesh->node_map(false);
  const Epetra_Map &node_map_use = mesh->node_map(true);

  int nnode_own = node_map_own.NumMyElements();
  int nnode_use = node_map_use.NumMyElements();

  vector<double> coord(3);
  vector<unsigned int> bad_nodes;

  Epetra_MultiVector coord_use(node_map_use,3);
  double **data;
  coord_use.ExtractView(&data);
  Epetra_MultiVector coord_own(View, node_map_own, data, 3);

  for (unsigned int j = 0; j < nnode_own; ++j) {
    mesh->node_to_coordinates(j, coord.begin(), coord.end());
    for (int k = 0; k < 3; ++k) coord_own[k][j] = coord[k];
  }

  Epetra_Import importer(node_map_use, node_map_own);
  coord_use.Import(coord_own, importer, Insert);

  for (unsigned int j = nnode_own; j < nnode_use; ++j) {
    mesh->node_to_coordinates(j, coord.begin(), coord.end());
    bool bad_data = false;
    for (int k = 0; k < 3; ++k)
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
// that the the GIDs of the nodes defining the face are the same, including
// their order (face orientation).

bool MeshAuditOld::check_face_to_nodes_ghost_data() const
{
  const Epetra_Map &node_map = mesh->node_map(true);
  const Epetra_Map &face_map_own = mesh->face_map(false);
  const Epetra_Map &face_map_use = mesh->face_map(true);

  int nface_own = face_map_own.NumMyElements();
  int nface_use = face_map_use.NumMyElements();

  vector<unsigned int> fnode(4);
  vector<unsigned int> bad_faces, bad_faces1, bad_faces2;

  // Create a matrix of the GIDs for all owned faces.
  Epetra_IntSerialDenseMatrix gids(nface_use,4); // no Epetra_IntMultiVector :(
  for (unsigned int j = 0; j < nface_own; ++j) {
    mesh->face_to_nodes(j, fnode.begin(), fnode.end());
    for (int k = 0; k < 4; ++k)
      gids(j,k) = node_map.GID(fnode[k]);
  }

  // Import these GIDs to all used faces; sets values on ghost faces.
  Epetra_Import importer(face_map_use, face_map_own);
  for (int k = 0; k < 4; ++k) {
    Epetra_IntVector kgids_own(View, face_map_own, gids[k]);
    Epetra_IntVector kgids_use(View, face_map_use, gids[k]);
    kgids_use.Import(kgids_own, importer, Insert);
  }

  // Compare the ghost face GIDs against the reference values just computed.
  for (unsigned int j = nface_own; j < nface_use; ++j) {
    mesh->face_to_nodes(j, fnode.begin(), fnode.end());
    bool bad_data = false;
    for (int k = 0; k < 4; ++k)
      if (node_map.GID(fnode[k]) != gids(j,k)) bad_data = true;
    if (bad_data) {
      // Determine just how bad the data is.
      vector<unsigned int> fnode_ref(4);
      for (int k = 0; k < 4; ++k) {
        fnode[k] = node_map.GID(fnode[k]);
        fnode_ref[k] = gids(j,k);
      }
      int n = same_face(fnode, fnode_ref);
      switch (n) {
      case 0: // completely bad -- different face
        bad_faces.push_back(j);
        break;
      case -1: // very bad -- same face but wrong orientation
        bad_faces1.push_back(j);
        break;
      case 1:  // not good -- same face and orientation, but not an exact copy
        bad_faces2.push_back(j);
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

  if (!bad_faces2.empty()) {
    os << "ERROR: found ghost faces that are not exact copies of their master:";
    write_list(bad_faces2, MAX_OUT);
    error = true; // some controversy whether this should be considered an error
  }

  return global_any(error);
}

// Check that ghost cells are exact copies of their master.  This means
// that the the GIDs of the nodes defining the cell are the same, including
// their order (face orientation).

bool MeshAuditOld::check_cell_to_nodes_ghost_data() const
{
  const Epetra_Map &node_map = mesh->node_map(true);
  const Epetra_Map &cell_map_own = mesh->cell_map(false);
  const Epetra_Map &cell_map_use = mesh->cell_map(true);

  int ncell_own = cell_map_own.NumMyElements();
  int ncell_use = cell_map_use.NumMyElements();

  vector<unsigned int> cnode(8);
  vector<unsigned int> bad_cells;

  // Create a matrix of the GIDs for all owned cells.
  Epetra_IntSerialDenseMatrix gids(ncell_use,8); // no Epetra_IntMultiVector :(
  for (unsigned int j = 0; j < ncell_own; ++j) {
    mesh->cell_to_nodes(j, cnode.begin(), cnode.end());
    for (int k = 0; k < 8; ++k)
      gids(j,k) = node_map.GID(cnode[k]);
  }

  // Import these GIDs to all used cells; sets values on ghost cells.
  Epetra_Import importer(cell_map_use, cell_map_own);
  for (int k = 0; k < 8; ++k) {
    Epetra_IntVector kgids_own(View, cell_map_own, gids[k]);
    Epetra_IntVector kgids_use(View, cell_map_use, gids[k]);
    kgids_use.Import(kgids_own, importer, Insert);
  }

  // Compare the ghost cell GIDs against the reference values just computed.
  for (unsigned int j = ncell_own; j < ncell_use; ++j) {
    mesh->cell_to_nodes(j, cnode.begin(), cnode.end());
    bool bad_data = false;
    for (int k = 0; k < 8; ++k)
      if (node_map.GID(cnode[k]) != gids(j,k)) bad_data = true;
    if (bad_data) bad_cells.push_back(j);
    // for bad cells we could do a further check to see how bad they are.
    // For example, maybe the GIDs are in a different order but the cell
    // orientation is the same, etc.
  }

  bool error = false;

  if (!bad_cells.empty()) {
    os << "ERROR: found bad data for ghost cells:";
    write_list(bad_cells, MAX_OUT);
    os << "       The ghost cells are not exact copies of their master." << endl;
    error = true;
  }

  return global_any(error);
}

// Check that ghost cells reference the same faces, and in the same order, as their
// master.  This is part of a ghost being an exact copies the master.  Even if the
// preceding ghost checks pass it is still possible for this one to fail.  In this
// case the ghost would reference a different face (by GID) but that face would
// reference the same nodes as the correct face.  So the two faces would be
// geometrically identical, including orientation, but be distinct.

bool MeshAuditOld::check_cell_to_faces_ghost_data() const
{
  const Epetra_Map &face_map = mesh->face_map(true);
  const Epetra_Map &cell_map_own = mesh->cell_map(false);
  const Epetra_Map &cell_map_use = mesh->cell_map(true);

  int ncell_own = cell_map_own.NumMyElements();
  int ncell_use = cell_map_use.NumMyElements();

  vector<unsigned int> cface(6);
  vector<unsigned int> bad_cells;

  // Create a matrix of the GIDs for all owned cells.
  Epetra_IntSerialDenseMatrix gids(ncell_use,6); // no Epetra_IntMultiVector :(
  for (unsigned int j = 0; j < ncell_own; ++j) {
    mesh->cell_to_faces(j, cface.begin(), cface.end());
    for (int k = 0; k < 6; ++k)
      gids(j,k) = face_map.GID(cface[k]);
  }

  // Import these GIDs to all used cells; sets values on ghost cells.
  Epetra_Import importer(cell_map_use, cell_map_own);
  for (int k = 0; k < 6; ++k) {
    Epetra_IntVector kgids_own(View, cell_map_own, gids[k]);
    Epetra_IntVector kgids_use(View, cell_map_use, gids[k]);
    kgids_use.Import(kgids_own, importer, Insert);
  }

  // Compare the ghost cell GIDs against the reference values just computed.
  for (unsigned int j = ncell_own; j < ncell_use; ++j) {
    mesh->cell_to_faces(j, cface.begin(), cface.end());
    bool bad_data = false;
    for (int k = 0; k < 6; ++k)
      if (face_map.GID(cface[k]) != gids(j,k)) bad_data = true;
    if (bad_data) bad_cells.push_back(j);
  }

  bool error = false;

  if (!bad_cells.empty()) {
    os << "ERROR: found bad data for ghost cells:";
    write_list(bad_cells, MAX_OUT);
    os << "       The ghost cells are not exact copies of their master." << endl;
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

bool MeshAuditOld::check_node_set_ids() const
{
  return check_get_set_ids(NODE);
}

bool MeshAuditOld::check_face_set_ids() const
{
  return check_get_set_ids(FACE);
}

bool MeshAuditOld::check_cell_set_ids() const
{
  return check_get_set_ids(CELL);
}

bool MeshAuditOld::check_get_set_ids(Entity_kind kind) const
{
  bool error = false;

  // Get the number of sets.
  int nset;
  try {
    nset = mesh->num_sets(kind); // this may fail
  } catch (...) {
    os << "ERROR: caught exception from num_sets()" << endl;
    error = true;
  }
  error = global_any(error);
  if (error) return error;

  // Get the vector of set IDs.
  vector<unsigned int> sids(nset, UINT_MAX);
  try {
    mesh->get_set_ids(kind, sids.begin(), sids.end()); // this may fail
  } catch (...) {
    os << "ERROR: caught exception from get_set_ids()" << endl;
    error = true;
  }
  error = global_any(error);
  if (error) return error;

  // Check to see that set ID values were actually assigned.  This assumes
  // UINT_MAX is not a valid set ID.  This is a little iffy; perhaps 0 should
  // be declared as an invalid set ID instead (the case for ExodusII), or
  // perhaps we should just skip this check.
  bool bad_data = false;
  for (int j = 0; j < nset; ++j)
    if (sids[j] == UINT_MAX) bad_data = true;
  if (bad_data) {
    os << "ERROR: get_set_ids() failed to set all values" << endl;
    error = true;
  }
  error = global_any(error);
  if (error) return error;

  // Verify that the vector of set IDs contains no duplicates.
  if (!distinct_values(sids)) {
    os << "ERROR: get_set_ids() returned duplicate IDs" << endl;
    // it would be nice to output the duplicates
    error = true;
  }
  error = global_any(error);
  if (error) return error;

  // In parallel, verify that each process returns the exact same result.
  if (comm.NumProc() > 1) {
    // Check the number of sets are the same.
    comm.Broadcast(&nset, 1, 0);
    if (nset != mesh->num_sets(kind)) {
      os << "ERROR: inconsistent num_sets() value" << endl;
      error = true;
    }
    error = global_any(error);

    if (!error) {
      // Broadcast the set IDs on processor 0.
      vector<unsigned int> sids(nset);
      mesh->get_set_ids(kind, sids.begin(), sids.end());
      int *sids0 = new int[nset];
      for (int j = 0; j < nset; ++j) sids0[j] = sids[j];
      comm.Broadcast(sids0, nset, 0);

      // Check the set IDs, using the vector on process 0 as the reference.
      bool bad_data = false;
      for (int j = 0; j < nset; ++j)
        if (sids[j] != sids0[j]) bad_data = true;
      if (bad_data) {
        os << "ERROR: get_set_ids() returned inconsistent values" << endl;
        error = true;
      }
      delete [] sids0;
    }
  }

  return global_any(error);
}

// Check that valid_set_id() returns the correct results.
// This is a collective test, returning a collective pass/fail result.

bool MeshAuditOld::check_valid_node_set_id() const
{
  return check_valid_set_id(NODE);
}

bool MeshAuditOld::check_valid_face_set_id() const
{
  return check_valid_set_id(FACE);
}

bool MeshAuditOld::check_valid_cell_set_id() const
{
  return check_valid_set_id(CELL);
}

bool MeshAuditOld::check_valid_set_id(Entity_kind kind) const
{
  // Get the list of set IDs.
  int nset = mesh->num_sets(kind); // this should not fail
  vector<unsigned int> sids(nset);
  mesh->get_set_ids(kind, sids.begin(), sids.end()); // this should not fail

  // Generate a valid/invalid set ID reference array that encompasses the known ones.
  int max_id = 0;
  for (int j = 0; j < nset; ++j)
    if (max_id < sids[j]) max_id = sids[j];
  vector<bool> valid(max_id+2, false);
  for (int j = 0; j < nset; ++j)
    valid[sids[j]] = true;

  vector<unsigned int> bad_sids1, bad_sids2;
  for (int n = 0; n < valid.size(); ++n) {
    if (valid[n] && !mesh->valid_set_id(n, kind)) bad_sids1.push_back(n);
    if (!valid[n] && mesh->valid_set_id(n, kind)) bad_sids2.push_back(n);
  }

  bool error = false;

  if (!bad_sids1.empty()) {
    os << "ERROR: valid_set_id() returned false for valid set IDs:";
    write_list(bad_sids1, MAX_OUT);
    error = true;
  }

  if (!bad_sids2.empty()) {
    os << "ERROR: valid_set_id() returned true for invalid set IDs:";
    write_list(bad_sids2, MAX_OUT);
    error = true;
  }

  return global_any(error);
}

// For each set, check that get_set successfully returns valid references to
// local entities, without duplicates, and that the used set is consistent
// with the owned set.

bool MeshAuditOld::check_node_sets() const
{
  return check_sets(NODE, mesh->node_map(false), mesh->node_map(true));
}

bool MeshAuditOld::check_face_sets() const
{
  return check_sets(FACE, mesh->face_map(false), mesh->face_map(true));
}

bool MeshAuditOld::check_cell_sets() const
{
  return check_sets(CELL, mesh->cell_map(false), mesh->cell_map(true));
}

bool MeshAuditOld::check_sets(Entity_kind kind,
                          const Epetra_Map &map_own, const Epetra_Map &map_use) const
{
  bool error = false;

  // Get the list of set IDs.
  int nset = mesh->num_sets(kind);
  vector<unsigned int> sids(nset);
  mesh->get_set_ids(kind, sids.begin(), sids.end());

  for (int n = 0; n < sids.size(); ++n) {
    os << "  Checking set ID=" << sids[n] << " ..." << endl;

    // Basic sanity checks of the owned and used sets.
    bool bad_set = check_get_set(sids[n], kind, OWNED, map_own) ||
                   check_get_set(sids[n], kind, USED,  map_use);
    bad_set = global_any(bad_set);

    // Verify the used set relates correctly to the owned set.
    if (!bad_set) bad_set = check_used_set(sids[n], kind, map_own, map_use);

    // OUGHT TO DO TESTING OF THE GHOST SETS

    if (bad_set) error = true;
  }

  return error;
}

// Basic sanity check on set values: no duplicates, and all LID values belong
// to the map.  This test runs independently on each process and returns a
// per-process pass/fail result.

bool MeshAuditOld::check_get_set(unsigned int sid, Entity_kind kind,
                             Parallel_type ptype, const Epetra_Map &map) const
{
  // Get the size of the set.
  int n;
  try {
    n = mesh->get_set_size(sid, kind, ptype); // this may fail
  } catch (...) {
    os << "  ERROR: caught exception from get_set_size()" << endl;
    return true;
  }

  // Get the set.
  vector<unsigned int> set(n, UINT_MAX);
  try {
    mesh->get_set(sid, kind, ptype, set.begin(), set.end());  // this may fail
  } catch (...) {
    os << "  ERROR: caught exception from get_set()" << endl;
    return true;
  }

  // Check that all values were assigned.
  bool bad_data = false;
  for (int j = 0; j < set.size(); ++j)
    if (set[j] == UINT_MAX) bad_data = true;
  if (bad_data) {
    os << "  ERROR: not all values assigned by get_set()" << endl;
    return true;
  }

  // Check that the LIDs in the set belong to the map.
  vector<unsigned int> bad_LIDs;
  for (int j = 0; j < set.size(); ++j)
    if (!map.MyLID(set[j])) bad_LIDs.push_back(set[j]);
  if (!bad_LIDs.empty()) {
    os << "  ERROR: set contains invalid LIDs:";
    write_list(bad_LIDs, MAX_OUT);
    return true;
  }

  // Check that there are no duplicates in the set.
  if (!distinct_values(set)) {
    os << "  ERROR: set contains duplicate LIDs." << endl;
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

bool MeshAuditOld::check_used_set(unsigned int sid, Entity_kind kind,
                               const Epetra_Map &map_own, const Epetra_Map &map_use) const
{
  if (comm.NumProc() == 1) {

    // In serial, the owned and used sets should be identical.

    int n = mesh->get_set_size(sid, kind, OWNED);
    vector<unsigned int> set_own(n);
    mesh->get_set(sid, kind, OWNED, set_own.begin(), set_own.end());

    // Set sizes had better be the same.
    if (mesh->get_set_size(sid, kind, USED) != set_own.size()) {
      os << "  ERROR: owned and used set sizes differ" << endl;
      return true;
    }

    // Verify that the two sets are identical.
    vector<unsigned int> set_use(n);
    mesh->get_set(sid, kind, USED, set_use.begin(), set_use.end());
    bool bad_data = false;
    for (int j = 0; j < n; ++j)
      if (set_use[j] != set_own[j]) bad_data = true;
    if (bad_data) {
      os << "  ERROR: owned and used sets differ" << endl;
      return true;
    }

    return false;

  } else {

    int n = mesh->get_set_size(sid, kind, OWNED);
    vector<unsigned int> set_own(n);
    mesh->get_set(sid, kind, OWNED, set_own.begin(), set_own.end());

    n = mesh->get_set_size(sid, kind, USED);
    vector<unsigned int> set_use(n);
    mesh->get_set(sid, kind, USED,  set_use.begin(), set_use.end());

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
    vector<unsigned int> bad_LIDs;
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

// Check that the alternate pointer-based get_set_ids and get_set accessors
// give the same results as the normative std::vector-based accessors.  This
// is a collective test, returning a collective pass/fail result.

bool MeshAuditOld::check_node_sets_alt() const
{
  return check_sets_alt(NODE);
}

bool MeshAuditOld::check_face_sets_alt() const
{
  return check_sets_alt(FACE);
}

bool MeshAuditOld::check_cell_sets_alt() const
{
  return check_sets_alt(CELL);
}

bool MeshAuditOld::check_sets_alt(Entity_kind kind) const
{
  bool error = false;

  int nset = mesh->num_sets(kind); // this should not fail
  vector<unsigned int> sids(nset);
  mesh->get_set_ids(kind, sids.begin(), sids.end()); // this should not fail
  try {
    unsigned int *sids1 = new unsigned int[nset];
    for (int j = 0; j < nset; ++j) sids1[j] = UINT_MAX;
    mesh->get_set_ids(kind, sids1, sids1+nset); // this may fail
    bool bad_data = false;
    for (int j = 0; j < nset; ++j)
      if (sids1[j] != sids[j]) bad_data = true;
    if (bad_data) {
      os << "ERROR: pointer-based get_set_ids() returned inconsistent values." << endl;
      error = true;
    }
  } catch (...) {
    os << "ERROR: caught exception from pointer-based get_set_ids()." << endl;
    error = true;
  }
  error = global_any(error);
  if (error) return error;

  for (int n = 0; n < sids.size(); ++n) {
    os << "  Checking set ID=" << sids[n] << " ..." << endl;

    bool bad_set = check_get_set_alt(sids[n], kind, OWNED) ||
                   check_get_set_alt(sids[n], kind, USED);

    // OUGHT TO DO TESTING OF THE GHOST SETS

    bad_set = global_any(bad_set);
    if (bad_set) error = true;
  }

  return error;
}

// Check that the alternate pointer-based get_set() returns the same result as
// the normative std::vector-based method.  This test runs independently on each
// process and returns a per-process pass/fail result.

bool MeshAuditOld::check_get_set_alt(unsigned int sid, Entity_kind kind,
                                  Parallel_type ptype) const
{
  int n = mesh->get_set_size(sid, kind, ptype);
  vector<unsigned int> set0(n);
  mesh->get_set(sid, kind, ptype, set0.begin(), set0.end()); // this should not fail
  unsigned int *set1 = new unsigned int[n];
  try {
    mesh->get_set(sid, kind, ptype, set1, set1+n);
    bool bad_data = false;
    for (int j = 0; j < n; ++j)
      if (set1[j] != set0[j]) bad_data = true;
    if (bad_data) {
      os << "  ERROR: pointer-based get_set() returned inconsistent values." << endl;
      return true;
    }
  } catch (...) {
    os << "  ERROR: caught exception from pointer-based get_set()" << endl;
    return true;
  }
  return false;
}

// SANE PARTITIONING CHECKS.  In parallel the cells, faces and nodes of the
// mesh are each partitioned across the processes, and while in principle
// these partitionings may be completely independent of each other, practical
// considerations lead to certain conditions that reasonable partitions should
// satisfy.  Taking the partitioning of the cells as given, one property that
// should be satisfied by the face and node partitionings is that the process
// that owns a particular face (node) must also own one of the cells containing
// the face (node).

bool MeshAuditOld::check_face_partition() const
{
  // Mark all the faces contained by owned cells.
  bool owned[nface];
  for (int j = 0; j < nface; ++j) owned[j] = false;
  unsigned int cface[6];
  for (unsigned int j = 0; j < mesh->cell_map(false).NumMyElements(); ++j) {
    mesh->cell_to_faces(j, cface, cface+6);
    for (int k = 0; k < 6; ++k) owned[cface[k]] = true;
  }

  // Verify that every owned face has been marked as belonging to an owned cell.
  vector<unsigned int> bad_faces;
  for (unsigned int j = 0; j < mesh->face_map(false).NumMyElements(); ++j)
    if (!owned[j]) bad_faces.push_back(j);

  if (!bad_faces.empty()) {
    os << "ERROR: found orphaned owned faces:";
    write_list(bad_faces, MAX_OUT);
    os << "       Process doesn't own either of the cells sharing the face." << endl;
    return true;
  }

  return false;
}


bool MeshAuditOld::check_node_partition() const
{
  // Mark all the nodes contained by owned cells.
  bool owned[nnode];
  for (int j = 0; j < nnode; ++j) owned[j] = false;
  unsigned int cnode[8];
  for (unsigned int j = 0; j < mesh->cell_map(false).NumMyElements(); ++j) {
    mesh->cell_to_nodes(j, cnode, cnode+8);
    for (int k = 0; k < 8; ++k) owned[cnode[k]] = true;
  }

  // Verify that every owned node has been marked as belonging to an owned cell.
  vector<unsigned int> bad_nodes;
  for (unsigned int j = 0; j < mesh->node_map(false).NumMyElements(); ++j)
    if (!owned[j]) bad_nodes.push_back(j);

  if (!bad_nodes.empty()) {
    os << "ERROR: found orphaned owned nodes:";
    write_list(bad_nodes, MAX_OUT);
    os << "       Process doesn't own any of the cells containing the node." << endl;
    return true;
  }

  return false;
}

// Returns true if the values in the list are distinct -- no repeats.

bool MeshAuditOld::distinct_values(const vector<unsigned int> &list) const
{
  vector<unsigned int> copy(list);
  sort(copy.begin(), copy.end());
  return (adjacent_find(copy.begin(),copy.end()) == copy.end());
}

//bool MeshAuditOld::distinct_values (const std::vector<unsigned int> &list) const
//{
//  for (int i = 0; i < list.size(); ++i)
//    for (int j = i+1; j < list.size(); ++j)
//      if (list[i] == list[j]) return false;
//  return true;
//}

// Returns 1 if the face node lists fnode1 and fnode2 describe the same face
// with the same orientation.  Returns -1 if the lists describe the same face
// but with opposite orientations.  Returns 0 if the lists describe different
// faces.  Implicitly assumes non-degenerate faces; the results are not
// reliable otherwise.

int MeshAuditOld::same_face(const vector<unsigned int> fnode1, const vector<unsigned int> fnode2) const
{
  // Locate position in fnode1 of fnode2[0].
  int n;
  for (n = 0; n < 4; ++n)
    if (fnode1[n] == fnode2[0]) break;
  if (n == 4) return 0; // did not find it -- different faces

  // Cyclic permutation p such that p[0] = n.
  int p[4];
  for (int i = 0; i < 4; ++i) p[i] = (n + i) % 4;

  // Check if fnode1[p] equals fnode2.
  for (n = 1; n < 4; ++n)
    if (fnode1[p[n]] != fnode2[n]) break;
  if (n == 4) return 1;  // they match

  // Modify the permutation to reverse the orientation of fnode1.
  int s = p[1];
  p[1] = p[3];
  p[3] = s;

  // Check if fnode1[p] equals fnode2.
  for (n = 1; n < 4; ++n)
    if (fnode1[p[n]] != fnode2[n]) break;
  if (n == 4) return -1;

  return 0; // different faces
}


void MeshAuditOld::write_list(const vector<unsigned int> &list, unsigned int max_out) const
{
  int num_out = min((unsigned int) list.size(), max_out);
  for (int i = 0; i < num_out; ++i) os << " " << list[i];
  if (num_out < list.size()) os << " [" << list.size()-num_out << " items omitted]";
  os << endl;
}

bool MeshAuditOld::global_any(bool value) const
{
  int lval=value, gval;
  comm.MaxAll(&lval, &gval, 1);
  return gval;
}

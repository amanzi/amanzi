#include "MeshAudit.hh"

#include <algorithm>
#include <cfloat>

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"

#include "cell_geometry.hh"
#include "cell_topology.hh"

#include <iostream>
#include <iomanip>

using namespace std;

// Need checks whether illegal usage of the interface causes exceptions

MeshAudit:: MeshAudit(Teuchos::RCP<Mesh_maps_base> &mesh_, ostream& os_) :
      mesh(mesh_), comm(mesh_->get_comm()), MyPID(mesh_->get_comm()->MyPID()),
      os(os_),
      nnode(mesh_->node_map(true).NumMyElements()),
      nface(mesh_->face_map(true).NumMyElements()),
      ncell(mesh_->cell_map(true).NumMyElements()),
      MAX_OUT(5)
    {}


int MeshAudit::Verify() const
{
  int ierr;
  bool fail = false;
  bool abort = false;
  
  // What I'd love to do here is to describe the dependencies between tests
  // (i.e., what tests must pass before a given test can be run) and generate
  // the dependency graph (tree) ala make, and then walk the tree running all
  // the tests that can be run.  Instead I've got the horrible crufty hack
  // you find below.

  ierr = check_entity_counts();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  ierr = check_cell_to_nodes_refs();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  ierr = check_cell_to_nodes_consistency();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  ierr = check_cell_to_faces_refs();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  ierr = check_cell_to_faces_consistency();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  ierr = check_face_to_nodes_refs();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  ierr = check_face_to_nodes_consistency();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  ierr = check_cell_to_face_dirs_basic();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  ierr = check_node_to_coordinates();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  if (abort) return 1;

  ierr = check_cell_degeneracy();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  if (abort) return 1;

  ierr = check_cell_to_faces();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  if (abort) return 1;

  ierr = check_cell_to_coordinates();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  ierr = check_face_to_coordinates();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  if (abort) return 1;

  ierr = check_cell_topology();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  ierr = check_node_maps();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  ierr = check_face_maps();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  ierr = check_cell_maps();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;
  
  if (abort) return 1;

  // Need to check all the "set" accessors

  ierr = check_node_to_coordinates_ghost_data();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;
  
  ierr = check_face_to_nodes_ghost_data();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  ierr = check_cell_to_nodes_ghost_data();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  ierr = check_cell_to_faces_ghost_data();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  if (fail)
    return 0;
  else
    return 1;
}

// The count_entities method should return values that match the number of
// elements in the corresponding Epetra_Maps.  This applies to nodes, faces
// and cells, including ghosts and not.  Here the maps are considered to be
// the authoritative source of information.  A positive value is returned
// if any discrepancy is found, but it is safe to perform other tests as
// they do not use these methods.

int MeshAudit::check_entity_counts() const
{
  int n, nref;
  int status=0;

  os << "Checking entity_counts ..." << endl;

  // Check the number of owned nodes.
  n = mesh->count_entities(Mesh_data::NODE,OWNED);
  nref = mesh->node_map(false).NumMyElements();
  if (n != nref) {
    os << ": ERROR: count_entities(NODE,OWNED)=" << n << "; should be " << nref << endl;
    status = 1;
  }

  // Check the number of used nodes.
  n = mesh->count_entities(Mesh_data::NODE,USED);
  nref = mesh->node_map(true).NumMyElements();
  if (n != nref) {
    os << "ERROR: count_entities(NODE,USED)=" << n << "; should be " << nref << endl;
    status = 1;
  }

  // Check the number of owned faces.
  n = mesh->count_entities(Mesh_data::FACE,OWNED);
  nref = mesh->face_map(false).NumMyElements();
  if (n != nref) {
    os << "ERROR: count_entities(FACE,OWNED)=" << n << "; should be " << nref << endl;
    status = 1;
  }

  // Check the number of used faces.
  n = mesh->count_entities(Mesh_data::FACE,USED);
  nref = mesh->face_map(true).NumMyElements();
  if (n != nref) {
    os << "ERROR: count_entities(FACE,USED)=" << n << "; should be " << nref << endl;
    status = 1;
  }

  // Check the number of owned cells.
  n = mesh->count_entities(Mesh_data::CELL,OWNED);
  nref = mesh->cell_map(false).NumMyElements();
  if (n != nref) {
    os << "ERROR: count_entities(CELL,OWNED)=" << n << "; should be " << nref << endl;
    status = 1;
  }

  // Check the number of used cells.
  n = mesh->count_entities(Mesh_data::CELL,USED);
  nref = mesh->cell_map(true).NumMyElements();
  if (n != nref) {
    os << "ERROR: count_entities(CELL,USED)=" << n << "; should be " << nref << endl;
    status = 1;
  }

  return status;
}

// The data returned by cell_to_nodes should reference valid local nodes.
// In addition, every node should be referenced at least one cell.
// Here we use the std::vector interface to the accessor.  Its consistency
// with alternative accessors is checked elsewhere.  A negative return
// value signals a terminal error: cell_to_nodes references invalid nodes
// and further checks using its data should be avoided.  A positive return
// value also signals an error, but further checks using cell_to_nodes data
// are safe.

int MeshAudit::check_cell_to_nodes_refs() const
{
  vector<unsigned int> bad_cells, bad_cells1;
  vector<unsigned int> free_nodes;
  vector<unsigned int> cnode(8);
  vector<unsigned int> refs(nnode, 0);

  os << "Checking cell_to_nodes references ..." << endl;

  for (unsigned int j = 0; j < ncell; ++j) {
    cnode.assign(8, UINT_MAX);
    try {
      mesh->cell_to_nodes(j, cnode.begin(), cnode.end()); // this may fail
      bool invalid_refs = false;
      for (int k = 0; k < cnode.size(); ++k) {
        if (cnode[k] >= 0 && cnode[k] < nnode)
          (refs[cnode[k]])++;
        else
          invalid_refs = true;
      }
      if (invalid_refs) bad_cells.push_back(j);
    } catch (...) {
      bad_cells1.push_back(j);
    }
  }

  unsigned int degree_min = UINT_MAX;
  unsigned int degree_max = 0;
  for (int j = 0; j < nnode; ++j) {
    if (refs[j] > 0) {
      degree_min = min(refs[j], degree_min);
      degree_max = max(refs[j], degree_max);
    } else {
      free_nodes.push_back(j);
    }
  }

  int status = 0;

  if (!bad_cells.empty()) {
    os << "ERROR: invalid nodes referenced by cells:";
    write_list(bad_cells, MAX_OUT);
    status = -1;
  }

  if (!bad_cells1.empty()) {
    os << "ERROR: caught exception for cells:";
    write_list(bad_cells1, MAX_OUT);
    status = -1;
  }

  if (!free_nodes.empty()) {
    os << "WARNING: found unreferenced nodes:";
    write_list(free_nodes, MAX_OUT);
    if (status == 0) status = 1;
  }

  return status;
}

// The data returned by cell_to_faces should reference valid local faces.
// In addition, every face should be referenced by exactly one or two cells.
// Here we use the std::vector interface to the accessor.  Its consistency
// with alternative accessors is checked elsewhere.  A negative return
// value signals a terminal error: cell_to_faces references invalid faces
// and further checks using its data should be avoided.  A positive return
// value also signals an error, but further checks using cell_to_faces data
// are safe.

int MeshAudit::check_cell_to_faces_refs() const
{
  vector<unsigned int> bad_cells, bad_cells1;
  vector<unsigned int> bad_faces;
  vector<unsigned int> free_faces;
  vector<unsigned int> cface(6);
  vector<unsigned int> refs(nface, 0);

  os << "Checking cell_to_faces references ..." << endl;

  for (unsigned int j = 0; j < ncell; ++j) {
    cface.assign(6, UINT_MAX);
    try {
      mesh->cell_to_faces(j, cface.begin(), cface.end()); // this may fail
      bool invalid_refs = false;
      for (int k = 0; k < cface.size(); ++k) {
        if (cface[k] >= 0 && cface[k] < nface)
          (refs[cface[k]])++;
        else
          invalid_refs = true;
      }
      if (invalid_refs) bad_cells.push_back(j);
    } catch (...) {
      bad_cells1.push_back(j);
    }
  }

  for (int j = 0; j < nface; ++j) {
    if (refs[j] == 0)
      free_faces.push_back(j);
    else if (refs[j] > 2)
      bad_faces.push_back(j);
  }

  int status = 0;

  if (!bad_cells.empty()) {
    os << "ERROR: invalid faces referenced by cells:";
    write_list(bad_cells, MAX_OUT);
    status = -1;
  }

  if (!bad_cells1.empty()) {
    os << "ERROR: caught exception for cells:";
    write_list(bad_cells1, MAX_OUT);
    status = -1;
  }

  if (!free_faces.empty()) {
    os << "WARNING: found unreferenced faces:";
    write_list(free_faces, MAX_OUT);
    if (status == 0) status = 1;
  }

  if (!bad_faces.empty()) {
    os << "WARNING: found faces shared by more than two cells:";
    write_list(bad_faces, MAX_OUT);
    if (status == 0) status = 1;
  }

  return status;
}

// The data returned by face_to_nodes should reference valid local nodes.
// In addition, every node should be referenced at least one face.
// Here we use the std::vector interface to the accessor.  Its consistency
// with alternative accessors is checked elsewhere.  A negative return
// value signals a terminal error: face_to_nodes references invalid nodes
// and further checks using its data should be avoided.  A positive return
// value also signals an error, but further checks using face_to_nodes data
// are safe.

int MeshAudit::check_face_to_nodes_refs() const
{
  vector<unsigned int> bad_faces, bad_faces1;
  vector<unsigned int> free_nodes;
  vector<unsigned int> fnode(4);
  vector<unsigned int> refs(nnode, 0);

  os << "Checking face_to_nodes references ..." << endl;

  for (unsigned int j = 0; j < nface; ++j) {
    fnode.assign(4, UINT_MAX);
    try {
      mesh->face_to_nodes(j, fnode.begin(), fnode.end()); // this may fail
      bool invalid_refs = false;
      for (int k = 0; k < fnode.size(); ++k) {
        if (fnode[k] >= 0 && fnode[k] < nnode)
          (refs[fnode[k]])++;
        else
          invalid_refs = true;
      }
      if (invalid_refs) bad_faces.push_back(j);
    } catch (...) {
      bad_faces1.push_back(j);
    }
  }

  unsigned int degree_min = UINT_MAX;
  unsigned int degree_max = 0;
  for (int j = 0; j < nnode; ++j) {
    if (refs[j] > 0) {
      degree_min = min(refs[j], degree_min);
      degree_max = max(refs[j], degree_max);
    } else {
      free_nodes.push_back(j);
    }
  }

  int status = 0;

  if (!bad_faces.empty()) {
    os << "ERROR: invalid nodes referenced by faces:";
    write_list(bad_faces, MAX_OUT);
    status = -1;
  }

  if (!bad_faces1.empty()) {
    os << "ERROR: caught exception for faces:";
    write_list(bad_faces1, MAX_OUT);
    status = -1;
  }

  if (!free_nodes.empty()) {
    os << "WARNING: found unreferenced nodes:";
    write_list(free_nodes, MAX_OUT);
    if (status == 0) status = 1;
  }

  return status;
}

// The mesh interface provides alternative methods for accessing
// connectivity data.  These alternates should return the same data.
// Here the accessors using std::vector are considered normative,
// and it is assumed that they have been verified to return valid
// (though perhaps not correct) results.  A positive value is
// returned if any inconsistency is found.  Further tests, which
// use only the std::vector-based accessors, are safe.

int MeshAudit::check_cell_to_nodes_consistency() const
{
  os << "Checking consistency of cell_to_nodes methods ..." << endl;

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

  int status = 0;

  if (!bad_cells1.empty()) {
    os << "ERROR: bad values from pointer-based accessor for cells:";
    write_list(bad_cells1, MAX_OUT);
    status = 1;
  }

  if (!bad_cells1_exc.empty()) {
    os << "ERROR: caught exception from pointer-based accessor for cells:";
    write_list(bad_cells1_exc, MAX_OUT);
    status = 1;
  }

  return status;
}


int MeshAudit::check_face_to_nodes_consistency() const
{
  os << "Checking consistency of face_to_nodes methods ..." << endl;

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

  int status = 0;

  if (!bad_faces1.empty()) {
    os << "ERROR: bad values from pointer-based accessor for faces:";
    write_list(bad_faces1, MAX_OUT);
    status = 1;
  }

  if (!bad_faces1_exc.empty()) {
    os << "ERROR: caught exception from pointer-based accessor for faces:";
    write_list(bad_faces1_exc, MAX_OUT);
    status = 1;
  }

  return status;
}


int MeshAudit::check_cell_to_faces_consistency() const
{
  os << "Checking consistency of cell_to_faces methods ..." << endl;

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

  int status = 0;

  if (!bad_cells1.empty()) {
    os << "ERROR: bad values from pointer-based accessor for cells:";
    write_list(bad_cells1, MAX_OUT);
    status = 1;
  }

  if (!bad_cells1_exc.empty()) {
    os << "ERROR: caught exception from pointer-based accessor for cells:";
    write_list(bad_cells1_exc, MAX_OUT);
    status = 1;
  }

  return status;
}

// Check that cell_to_face_dirs successfully returns data for all cells
// and that the values are either +1 or -1.  Also check that the alternate
// accessor methods return the same data.  The std::vector-based method is
// considered normative here.  This basic check does not check the correctness
// of the values, only that the values are acceptable.

int MeshAudit::check_cell_to_face_dirs_basic() const
{
  os << "Checking cell_to_face_dirs (basic) ..." << endl;

  vector<int> fdirs0(6);
  vector<unsigned int> bad_cells0;
  vector<unsigned int> bad_cells1;

  for (unsigned int j = 0; j < ncell; ++j) {
    fdirs0.assign(6, INT_MAX);
    try {
      mesh->cell_to_face_dirs(j, fdirs0.begin(), fdirs0.end());  // this may fail
      bool bad_data = false;
      for (int k = 0; k < fdirs0.size(); ++k)
        if (fdirs0[k] != -1 && fdirs0[k] != 1) bad_data = true;
      if (bad_data)
        bad_cells0.push_back(j);
      else {
        int fdirs1[6] = { INT_MAX };
        try {
          mesh->cell_to_face_dirs(j, fdirs1, fdirs1+6); // this may fail
          bad_data = false;
          for (int k = 0; k < 6; ++k)
            if (fdirs1[k] != fdirs0[k]) bad_data = true;
          if (bad_data) bad_cells1.push_back(j);
        } catch (...) {
          bad_cells1.push_back(j);
        }
      }
    } catch (...) {
      bad_cells0.push_back(j);
    }
  }

  int status = 0;

  if (!bad_cells0.empty()) {
    os << "ERROR: inadmissable or no data for cells:";
    write_list(bad_cells0, MAX_OUT);
    status = -1;
  }

  if (!bad_cells1.empty()) {
    os << "ERROR: pointer-based accessor returned inconsistent values for cells:";
    write_list(bad_cells1, MAX_OUT);
    if (status == 0) status = -1;
  }

  return status;
}

// Checks that cells are not topologically degenerate (repeated node index).
// A negative value is returned if any such cells are found, and in this case
// many further tests involving the cell_to_nodes data should be avoided.

int MeshAudit::check_cell_degeneracy() const
{
  os << "Checking cells for topological degeneracy ..." << endl;

  vector<unsigned int> cnode(8);
  vector<unsigned int> bad_cells;

  for (unsigned int j = 0; j < ncell; ++j) {
    mesh->cell_to_nodes(j, cnode.begin(), cnode.end()); // should not fail
    if (!distinct_values(cnode)) bad_cells.push_back(j);
  }

  if (!bad_cells.empty()) {
    os << "ERROR: found topologically degenerate cells:";
    write_list(bad_cells, MAX_OUT);
    return -1;
  } else {
    return 0;
  }
}

// Check that cell_to_faces is returning the correct values by composing
// that map with face_to_nodes and comparing against the results returned
// by cell_to_nodes and the local face numbering convention described by
// cell_topology::HexFaceVert.  Also check that the relative orientation
// value returned by cell_to_face_dirs is correct.  A negative return
// value indicates cell_to_faces and/or cell_to_face_dirs return incorrect
// results and that further tests using their values should be avoided.

int MeshAudit::check_cell_to_faces() const
{
  os << "Checking correctness of cell_to_faces and cell_to_face_dirs ..." << endl;

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
        fnode_ref[i] = cnode[cell_topology::HexFaceVert[k][i]];
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

  int status = 0;

  if (!bad_cells0.empty()) {
    os << "ERROR: bad cell_to_faces values for cells:";
    write_list(bad_cells0, MAX_OUT);
    status = -1;
  }

  if (!bad_cells1.empty()) {
    os << "ERROR: bad cell_to_face_dirs values for cells:";
    write_list(bad_cells1, MAX_OUT);
    status = -1;
  }

  return status;
}

// Check that node_to_coordinates successfully returns data for all nodes
// and that the alternative accessor methods return the same data.
// The pointer-based method is considered normative here.  A negative return
// value signals a terminal error: node_to_coordinates did not return data
// for some nodes and further checks involving its data should be avoided.
// A postive return signals an error with the alternative accessor methods,
// but other checks, which use the pointer-based accessor, are safe.

int MeshAudit::check_node_to_coordinates() const
{
  os << "Checking node_to_coordinates ..." << endl;

  vector<double> x1(3);
  vector<unsigned int> bad_nodes0;
  vector<unsigned int> bad_nodes1;

  for (unsigned int j = 0; j < nnode; ++j) {
    double x0[3] = { DBL_MAX };
    try {
      mesh->node_to_coordinates(j, x0, x0+3); // this may fail
      bool bad_data = false;
      for (int k = 0; k < 3; ++k)
        if (x0[k] == DBL_MAX) bad_data = true;
      if (bad_data)
        bad_nodes0.push_back(j);
      else {
        // Check consistency of the other accessor.
        x1.assign(3,DBL_MAX);
        try {
          mesh->node_to_coordinates(j, x1.begin(), x1.end()); // this may fail
          bad_data = false;
          for (int k = 0; k < x1.size(); ++k)
            if (x1[k] != x0[k]) bad_data = true;
          if (bad_data) bad_nodes1.push_back(j);
        } catch (...) {
          bad_nodes1.push_back(j);
        }
      }
    } catch (...) {
      bad_nodes0.push_back(j);
    }
  }

  int status = 0;

  if (!bad_nodes0.empty()) {
    os << "ERROR: no coordinate data for nodes:";
    write_list(bad_nodes0, MAX_OUT);
    status = -1;
  }

  if (!bad_nodes1.empty()) {
    os << "ERROR: std::vector-based accessor returned inconsistent values for nodes:";
    write_list(bad_nodes1, MAX_OUT);
    if (status == 0) status = 1;
  }

  return status;
}

// Check that cell_to_coordinates successfully returns data for all cells, and
// that this data is identical to that returned by composing node_to_coordinates
// with cell_to_nodes.  This assumes that the std::vector form of cell_to_nodes
// and pointer-based form of node_to_coordinates have been verified to return
// valid data.  Also check that the alternative cell_to_coordinates accessors
// return the same data as the normative pointer-based method. A negative return
// value signals a terminal error: cell_to_coordinates returned incorrect data
// for some cells and further checks using its data should be avoided.  A
// positive return value signals an error with the alternative accessor methods,
// but other checks, which use the pointer-based method, are safe.

int MeshAudit::check_cell_to_coordinates() const
{
  os << "Checking cell_to_coordinates ..." << endl;

  double xref[24]; // 3x8
  vector<double> x1(24);
  vector<unsigned int> cnode(8);
  vector<unsigned int> bad_cells0;
  vector<unsigned int> bad_cells1;

  for (unsigned int j = 0; j < ncell; ++j) {
    double x0[24] = { DBL_MAX };
    try {
      mesh->cell_to_coordinates(j, x0, x0+24); // this may fail
      mesh->cell_to_nodes(j, cnode.begin(), cnode.end()); // this should not fail
      double *xbeg = xref;
      for (int k = 0; k < 8; ++k) {
        mesh->node_to_coordinates(cnode[k], xbeg, xbeg+3); // this should not fail
        xbeg += 3;
      }
      bool bad_data = false;
      for (int k = 0; k < 24; ++k)
        if (x0[k] != xref[k]) bad_data = true;
      if (bad_data)
        bad_cells0.push_back(j);
      else {
        // Check consistency of other accessors.
        x1.assign(24, DBL_MAX);
        try {
          mesh->cell_to_coordinates(j, x1.begin(), x1.end()); // this may fail
          bad_data = false;
          for (int k = 0; k < x1.size(); ++k)
            if (x1[k] != x0[k]) bad_data = true;
          if (bad_data) bad_cells1.push_back(j);
        } catch (...) {
          bad_cells1.push_back(j);
        }
      }
    } catch (...) {
      bad_cells0.push_back(j);
    }
  }

  int status = 0;

  if (!bad_cells0.empty()) {
    os << "ERROR: bad cell_to_coordinates data for cells:";
    write_list(bad_cells0, MAX_OUT);
    status = -1;
  }

  if (!bad_cells1.empty()) {
    os << "ERROR: std::vector-based accessor returned inconsistent values for cells:";
    write_list(bad_cells1, MAX_OUT);
    if (status == 0) status = 1;
  }

  return status;
}

// Check that face_to_coordinates successfully returns data for all faces, and
// that this data is identical to that returned by composing node_to_coordinates
// with face_to_nodes.  This assumes that the std::vector form of face_to_nodes
// and pointer-based form of node_to_coordinates have been verified to return
// valid data.  Also check that the alternative face_to_coordinates accessors
// return the same data as the normative pointer-based method. A negative return
// value signals a terminal error: face_to_coordinates returned incorrect data
// for some faces and further checks using its data should be avoided.  A
// positive return value signals an error with the alternative accessor methods,
// but other checks, which use the pointer-based accessor, are safe.

int MeshAudit::check_face_to_coordinates() const
{
  os << "Checking face_to_coordinates ..." << endl;

  double xref[12]; // 3x4
  vector<double> x1(12);
  vector<unsigned int> fnode(4);
  vector<unsigned int> bad_faces0;
  vector<unsigned int> bad_faces1;

  for (unsigned int j = 0; j < nface; ++j) {
    double x0[12] = { DBL_MAX };
    try {
      mesh->face_to_coordinates(j, x0, x0+12); // this may fail
      mesh->face_to_nodes(j, fnode.begin(), fnode.end()); // this should not fail
      double *xbeg = xref;
      for (int k = 0; k < 4; ++k) {
        mesh->node_to_coordinates(fnode[k], xbeg, xbeg+3); // this should not fail
        xbeg += 3;
      }
      bool bad_data = false;
      for (int k = 0; k < 12; ++k)
        if (x0[k] != xref[k]) bad_data = true;
      if (bad_data)
        bad_faces0.push_back(j);
      else {
        // Check consistency of other accessors.
        x1.assign(12, DBL_MAX);
        try {
          mesh->face_to_coordinates(j, x1.begin(), x1.end()); // this may fail
          bad_data = false;
          for (int k = 0; k < x1.size(); ++k)
            if (x1[k] != x0[k]) bad_data = true;
          if (bad_data) bad_faces1.push_back(j);
        } catch (...) {
          bad_faces1.push_back(j);
        }
      }
    } catch (...) {
      bad_faces0.push_back(j);
    }
  }

  int status = 0;

  if (!bad_faces0.empty()) {
    os << "ERROR: bad face_to_coordinates data for faces:";
    write_list(bad_faces0, MAX_OUT);
    status = -1;
  }

  if (!bad_faces1.empty()) {
    os << "ERROR: std::vector-based accessor returned inconsistent values for faces:";
    write_list(bad_faces1, MAX_OUT);
    if (status == 0) status = 1;
  }

  return status;
}

// The hexahedral cells must not be degenerate, either topologically (repeated
// node index) or geometrically (with coincident nodes).  In addition, the
// cell vertices must be ordered in a pre-defined manner (see cell_topology).
// To detect geometric degeneracy or bad topology (vertices listed in the
// wrong order), the corner tet volumes of the hexahedron are evaluated and
// checked for positivity.  If any are negative this indicates bad topology
// (is this sufficient?  I think so.)  Geometric degeneracy is indicated
// by one or more zero corner volumes.

int MeshAudit::check_cell_topology() const
{
  os << "Checking cell topology/geometry ..." << endl;

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

  if (!bad_cells.empty()) {
    os << "ERROR: found cells with non-positive volumes:";
    write_list(bad_cells, MAX_OUT);
    os << "       Cells either geometrically degenerate or have bad topology.";
    return 1;
  } else {
    return 0;
  }
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

int MeshAudit::check_node_maps() const
{
  os << "Checking owned and overlap node maps ..." << endl;
  return check_maps(mesh->node_map(false), mesh->node_map(true));
}

int MeshAudit::check_face_maps() const
{
  os << "Checking owned and overlap face maps ..." << endl;
  return check_maps(mesh->face_map(false), mesh->face_map(true));
}

int MeshAudit::check_cell_maps() const
{
  os << "Checking owned and overlap cell maps ..." << endl;
  return check_maps(mesh->cell_map(false), mesh->cell_map(true));
}

int MeshAudit::check_maps(const Epetra_Map &map_own, const Epetra_Map &map_use) const
{
  int status = 0;

  // Local index should start at 0.
  if (map_own.IndexBase() != 0) {
    os << "ERROR: the owned map's index base is not 0." << endl;
    status = -1;
  }
  if (map_use.IndexBase() != 0) {
    os << "ERROR: the overlap map's index base is not 0." << endl;
    status = -1;
  }

  // Check that the owned map is 1-1.
  if (!map_own.UniqueGIDs()) {
    os << "ERROR: owned map is not 1-to-1" << endl;
    status = -1;
  }

  if (status != 0) return status;

  if (comm->NumProc() == 1)
  {

    // Serial or 1-process MPI

    if (!map_use.SameAs(map_own)) {
      os << "ERROR: the overlap map differs from the owned map (single process)." << endl;
      return -1;
    } else {
      return 0;
    }

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
      status = -1;
    }

    // Verify that the overlap indices are owned by other processes.
    int *gids = new int[num_ovl];
    int *pids = new int[num_ovl];
    int *lids = new int[num_ovl];
    for (int j = 0; j < num_ovl; ++j) gids[j] = map_use.GID(j+num_own);
    map_own.RemoteIDList(num_ovl, gids, pids, lids);
    bad_map = false;
    for (int j = 0; j < num_ovl; ++j)
      if (pids[j] < 0 || pids[j] == comm->MyPID()) bad_map = true;
    if (bad_map) {
      os << "ERROR: invalid ghosts in overlap map." << endl;
      status = -1;
    }

    // Look for duplicates among the overlap indices.
    vector<int> ovl_gids(gids, gids+num_ovl);
    sort(ovl_gids.begin(), ovl_gids.end());
    if (adjacent_find(ovl_gids.begin(),ovl_gids.end()) != ovl_gids.end()) {
      os << "ERROR: duplicate ghosts in overlap map." << endl;
      status = -1;
    }

    delete [] lids;
    delete [] pids;
    delete [] gids;

    return status;
  }
}

// Check that ghost nodes are exact copies of their master.
// This simply means that they have the same coordinates.

int MeshAudit::check_node_to_coordinates_ghost_data() const
{
  os << "Checking node_to_coordinates ghost data ..." << endl;
  
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
  
  int status = 0;
  
  if (!bad_nodes.empty()) {
    os << "ERROR: found ghost nodes with incorrect coordinates:";
    write_list(bad_nodes, MAX_OUT);
    status = 1;
  }
  
  return status;
}

// Check that ghost faces are exact copies of their master.  This means
// that the the GIDs of the nodes defining the face are the same, including
// their order (face orientation).

int MeshAudit::check_face_to_nodes_ghost_data() const
{
  os << "Checking face_to_nodes ghost data ..." << endl;
  
  const Epetra_Map &node_map = mesh->node_map(true);
  const Epetra_Map &face_map_own = mesh->face_map(false);
  const Epetra_Map &face_map_use = mesh->face_map(true);
  
  int nface_own = face_map_own.NumMyElements();
  int nface_use = face_map_use.NumMyElements();
  
  vector<unsigned int> fnode(4);
  vector<unsigned int> bad_faces;
  
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
    if (bad_data) bad_faces.push_back(j);
    // for bad faces we could do a further check to see how bad they are.
    // For example, maybe the GIDs are in a different order but the face
    // orientation is the same, etc.
  }
  
  int status = 0;
  
  if (!bad_faces.empty()) {
    os << "ERROR: found bad data for ghost faces:";
    write_list(bad_faces, MAX_OUT);
    os << "       The ghost faces are not exact copies of their master." << endl;
    status = 1;
  }
  
  return status;
}

// Check that ghost cells are exact copies of their master.  This means
// that the the GIDs of the nodes defining the cell are the same, including
// their order (face orientation).

int MeshAudit::check_cell_to_nodes_ghost_data() const
{
  os << "Checking cell_to_nodes ghost data ..." << endl;
  
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
  
  int status = 0;
  
  if (!bad_cells.empty()) {
    os << "ERROR: found bad data for ghost cells:";
    write_list(bad_cells, MAX_OUT);
    os << "       The ghost cells are not exact copies of their master." << endl;
    status = 1;
  }
  
  return status;
}

// Check that ghost cells reference the same faces, and in the same order, as their
// master.  This is part of a ghost being an exact copies the master.  Even if the
// preceding ghost checks pass it is still possible for this one to fail.  In this
// case the ghost would reference a different face (by GID) but that face would
// reference the same nodes as the correct face.  So the two faces would be
// geometrically identical, including orientation, but be distinct.

int MeshAudit::check_cell_to_faces_ghost_data() const
{
  os << "Checking cell_to_face ghost data ..." << endl;
  
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
  
  int status = 0;
  
  if (!bad_cells.empty()) {
    os << "ERROR: found bad data for ghost cells:";
    write_list(bad_cells, MAX_OUT);
    os << "       The ghost cells are not exact copies of their master." << endl;
    status = 1;
  }
  
  return status;
}

// Check that num_sets successfully completes and that the returned value
// is the same across all processors.  Independent checks for NODE, FACE
// and CELL entities.

// Check that get_set_ids successfully completes and that the returned data is
// the same (including order) across all processors.  Check that the alternate
// methods return identical data.  Independent checks for NODE, FACE and CELL
// entities.

// Check that valid_set_id returns true for all set IDs, and false for some
// selection of invalid IDs.  Independent checks for NODE, FACE and CELL entities.

// Check that get_set_size successfully completes and that the value for
// for USED equals the sum of the values for OWNED and GHOST.  Check that
// get_set successfully completes for OWNED, GHOST and USED.
// Check that the lists (OWNED, GHOST, USED) have distinct values.
// Check that the OWNED list LIDs are valid LIDs from the owned map.
// Check that the USED list LIDs are valid LIDs from the used map.
// Check that if the OWNED list were exported f
// Check that the GHOST list LIDs are in fact ghost LIDs
// Check that the USED list equals the concatenation of the OWNED and GHOST
// lists.  (This is 

//int MeshAudit::check_set(const Epetra_Map &map_own; const vector<unsigned int> &set_own,
//                         const Epetra_Map &map_use; const vector<unsigned int> &set_use) const
//{
//  // Check that 
//  if (!distinct_values(set_use)) {
//  }
//  
//  bool bad_set = false;
//  for (int j = 0; j < set_own.size(); ++j)
//    if (!map_own.MyLID(set_own[j])) bad_set = true;
//  if (bad_set){
//  }
//  
//  // Tag all LIDs in the used map that should belong to the used set;
//  // the owned set LIDs are taken as definitive.
//  Epetra_IntVector tag_use(map_use); // fills with zero values
//  int *tag_data;
//  tag_use.ExtractView(&tag_data);
//  Epetra_IntVector tag_own(View, map_own, tag_data);
//  for (int j = 0; j < set_own.size(); ++j) tag_own[set_own[j]] = 1;
//  Epetra_Import importer(map_use, map_own);
//  tag_use.Import(tag_own, importer, Insert);
//  for (int j = 0; j < set_use.size(); ++j) --tag_use[set_use[j]];
//  
//  // Check for negative tag values;
//  // these indicate a ghost LID that shouldn't be in the set but is.
//  bad_set = false;
//  vector<unsigned int> bad_LIDs;
//  for (int j = 0; j < set_own.size(); ++j)
//    if (tag_use[j] < 0) bad_set.push_back(j);
//  if (!bad_LIDs.empty()) {
//    os << "ERROR: found LIDs that 
//  }
//  
//  // Positive values indicate a ghost LID that should be in the set but isn't.
//  
//  
//
//}

// Returns true if the values in the list are distinct -- no repeats.

bool MeshAudit::distinct_values(const vector<unsigned int> &list) const
{
  vector<unsigned int> copy(list);
  sort(copy.begin(), copy.end());
  return (adjacent_find(copy.begin(),copy.end()) == copy.end());
}

//bool MeshAudit::distinct_values (const std::vector<unsigned int> &list) const
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

int MeshAudit::same_face(const vector<unsigned int> fnode1, const vector<unsigned int> fnode2) const
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


void MeshAudit::write_list(const vector<unsigned int> &list, unsigned int max_out) const
{
  int num_out = min((unsigned int) list.size(), max_out);
  for (int i = 0; i < num_out; ++i) os << " " << list[i];
  if (num_out < list.size()) os << " [" << list.size()-num_out << " items omitted]";
  os << endl;
}


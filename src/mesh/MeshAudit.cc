#include "MeshAudit.hh"

#include <algorithm>
#include <cfloat>

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"

#include "Geometry.hh"

#include <iostream>
#include <iomanip>

using namespace std;
using namespace Amanzi;
using namespace AmanziMesh;
using namespace AmanziGeometry;

// Need checks whether illegal usage of the interface causes exceptions

MeshAudit:: MeshAudit(Teuchos::RCP<Mesh> &mesh_, ostream& os_) :
      mesh(mesh_), comm(*(mesh_->get_comm())), MyPID(mesh_->get_comm()->MyPID()),
      os(os_),
      nnode(mesh_->node_epetra_map(true).NumMyElements()),
      nface(mesh_->face_epetra_map(true).NumMyElements()),
      ncell(mesh_->cell_epetra_map(true).NumMyElements()),
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

  ierr = check_cell_to_faces_refs();
  if (ierr < 0) abort = true;
  if (ierr != 0) fail = true;

  ierr = check_face_to_nodes_refs();
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

  ierr = check_cell_geometry();
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

  if (abort) return 1;

  ierr = check_node_partition();
  if (ierr != 0) fail = true;
  
  ierr = check_face_partition();
  if (ierr != 0) fail = true;

  ierr = check_node_sets();
  if (ierr) fail = true;

  ierr = check_face_sets();
  if (ierr) fail = true;

  ierr = check_cell_sets();
  if (ierr) fail = true;

  if (fail)
    return 1;
  else
    return 0;
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
  n = mesh->num_entities(NODE,OWNED);
  nref = mesh->node_epetra_map(false).NumMyElements();
  if (n != nref) {
    os << ": ERROR: num_entities(NODE,OWNED)=" << n << "; should be " << nref << endl;
    status = 1;
  }

  // Check the number of used nodes.
  n = mesh->num_entities(NODE,USED);
  nref = mesh->node_epetra_map(true).NumMyElements();
  if (n != nref) {
    os << "ERROR: num_entities(NODE,USED)=" << n << "; should be " << nref << endl;
    status = 1;
  }

  // Check the number of owned faces.
  n = mesh->num_entities(FACE,OWNED);
  nref = mesh->face_epetra_map(false).NumMyElements();
  if (n != nref) {
    os << "ERROR: num_entities(FACE,OWNED)=" << n << "; should be " << nref << endl;
    status = 1;
  }

  // Check the number of used faces.
  n = mesh->num_entities(FACE,USED);
  nref = mesh->face_epetra_map(true).NumMyElements();
  if (n != nref) {
    os << "ERROR: num_entities(FACE,USED)=" << n << "; should be " << nref << endl;
    status = 1;
  }

  // Check the number of owned cells.
  n = mesh->num_entities(CELL,OWNED);
  nref = mesh->cell_epetra_map(false).NumMyElements();
  if (n != nref) {
    os << "ERROR: num_entities(CELL,OWNED)=" << n << "; should be " << nref << endl;
    status = 1;
  }

  // Check the number of used cells.
  n = mesh->num_entities(CELL,USED);
  nref = mesh->cell_epetra_map(true).NumMyElements();
  if (n != nref) {
    os << "ERROR: num_entities(CELL,USED)=" << n << "; should be " << nref << endl;
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
  Entity_ID_List bad_cells, bad_cells1;
  Entity_ID_List free_nodes;
  Entity_ID_List cnode;
  vector<unsigned int> refs(nnode, 0);

  os << "Checking cell_to_nodes references ..." << endl;

  for (unsigned int j = 0; j < ncell; ++j) {
    try {
      mesh->cell_get_nodes(j, &cnode); // this may fail
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

// The data returned by cell_to_faces should reference valid local
// faces.  In addition, every face should be referenced by exactly one
// or two cells.  Here we use the std::vector interface to the
// accessor.  A negative return value signals a terminal error:
// cell_to_faces references invalid faces and further checks using its
// data should be avoided.


int MeshAudit::check_cell_to_faces_refs() const
{
  Entity_ID_List bad_cells, bad_cells1;
  Entity_ID_List bad_faces;
  Entity_ID_List free_faces;
  Entity_ID_List cface;
  vector<unsigned int> refs(nface, 0);

  os << "Checking cell_to_faces references ..." << endl;

  for (unsigned int j = 0; j < ncell; ++j) {
    try {
      mesh->cell_get_faces(j, &cface); // this may fail
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
  Entity_ID_List bad_faces, bad_faces1;
  Entity_ID_List free_nodes;
  Entity_ID_List fnode;
  vector<unsigned int> refs(nnode, 0);

  os << "Checking face_to_nodes references ..." << endl;

  for (unsigned int j = 0; j < nface; ++j) {
    try {
      mesh->face_get_nodes(j, &fnode); // this may fail
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


// Check that cell_to_face_dirs successfully returns data for all
// cells and that the values are either +1 or -1. The
// std::vector-based method is considered normative here.  This basic
// check does not check the correctness of the values, only that the
// values are acceptable.

int MeshAudit::check_cell_to_face_dirs_basic() const
{
  os << "Checking cell_to_face_dirs (basic) ..." << endl;

  vector<int> fdirs0;
  Entity_ID_List bad_cells0;
  Entity_ID_List bad_cells1;

  for (unsigned int j = 0; j < ncell; ++j) {
    fdirs0.assign(6, INT_MAX);
    try {
      mesh->cell_get_face_dirs(j, &fdirs0);  // this may fail
      bool bad_data = false;
      for (int k = 0; k < fdirs0.size(); ++k)
        if (fdirs0[k] != -1 && fdirs0[k] != 1) bad_data = true;
      if (bad_data)
        bad_cells0.push_back(j);
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

  return status;
}

// Checks that cells are not topologically degenerate (repeated node index).
// A negative value is returned if any such cells are found, and in this case
// many further tests involving the cell_to_nodes data should be avoided.

int MeshAudit::check_cell_degeneracy() const
{
  os << "Checking cells for topological degeneracy ..." << endl;

  Entity_ID_List cnode;
  Entity_ID_List bad_cells;

  for (unsigned int j = 0; j < ncell; ++j) {
    mesh->cell_get_nodes(j, &cnode); // should not fail
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



// NEEDS TO BE REWORKED

// Check that cell_get_faces is returning the correct values by composing
// that map with face_get_nodes and comparing against the results returned
// by cell_get_nodes and the local face numbering convention described by
// cell_topology::HexFaceVert.  Also check that the relative orientation
// value returned by cell_get_face_dirs is correct.  A negative return
// value indicates cell_get_faces and/or cell_get_face_dirs return incorrect
// results and that further tests using their values should be avoided.

int MeshAudit::check_cell_to_faces() const
{
  os << "Checking correctness of cell_get_faces and cell_get_face_dirs ..." << endl;

  Entity_ID_List cnode;
  Entity_ID_List cface;
  Entity_ID_List fnode_ref;
  Entity_ID_List fnode;
  vector<int> fdirs;
  Entity_ID_List bad_cells0;
  Entity_ID_List bad_cells1;

  for (unsigned int j = 0; j < ncell; ++j) {
    Cell_type ctype = mesh->cell_get_type(j);

    // If this is a general, non-standard element there is nothing to
    // to check against

    if (ctype == UNKNOWN || ctype == POLYGON || ctype == POLYHED)
      continue;

    mesh->cell_get_nodes(j, &cnode); // this should not fail

    mesh->cell_get_faces(j, &cface); // this should not fail
    mesh->cell_get_face_dirs(j, &fdirs); // this should not fail

    bool bad_face = false;
    bool bad_dir  = false;

    if (cface.size() != nface_std[ctype])
      bad_face = true;
    else {

      for (int k = 0; k < cface.size(); ++k) {

	mesh->face_get_nodes(cface[k], &fnode); // this should not fail

	int nfn = nfnodes_std[ctype][k];

	if (fnode.size() != nfn) {
	  bad_face = true;
	  break;
	}

	fnode_ref.clear();
	for (int i = 0; i < nfn; ++i) {
	  int nodenum = fnodes_std[ctype][k][i];
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

  int status = 0;

  if (!bad_cells0.empty()) {
    os << "ERROR: bad cell_get_faces values for cells:";
    write_list(bad_cells0, MAX_OUT);
    status = -1;
  }

  if (!bad_cells1.empty()) {
    os << "ERROR: bad cell_get_face_dirs values for cells:";
    write_list(bad_cells1, MAX_OUT);
    status = -1;
  }

  return status;
}

// Check that node_get_coordinates successfully returns data for all nodes
// A negative return value signals a terminal error

int MeshAudit::check_node_to_coordinates() const
{
  os << "Checking node_get_coordinates ..." << endl;

  int spdim = mesh->space_dimension();
  Entity_ID_List bad_nodes0;

  for (unsigned int j = 0; j < nnode; ++j) {
    Point x0(spdim);
    x0.set(DBL_MAX);
    try {
      mesh->node_get_coordinates(j, &x0); // this may fail
      bool bad_data = false;
      for (int k = 0; k < spdim; ++k)
        if (x0[k] == DBL_MAX) bad_data = true;
      if (bad_data)
        bad_nodes0.push_back(j);
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

  return status;
}

// Check that cell_get_coordinates successfully returns data for all
// cells, and that the data is identical to that returned by composing
// node_get_coordinates with cell_get_nodes.  This assumes that the
// vector form of cell_get_nodes have been verified to return valid
// data. A negative return value signals a terminal error:
// cell_get_coordinates returned incorrect data for some cells and
// further checks using its data should be avoided.

int MeshAudit::check_cell_to_coordinates() const
{
  os << "Checking cell_get_coordinates ..." << endl;

  Entity_ID_List cnode;
  Entity_ID_List bad_cells0;
  int spdim = mesh->space_dimension();

  for (unsigned int j = 0; j < ncell; ++j) {
    vector<Point> x0;
    try {
      mesh->cell_get_coordinates(j, &x0); // this may fail
      mesh->cell_get_nodes(j, &cnode); // this should not fail
      for (int k = 0; k < cnode.size(); ++k) {
	Point xref(spdim);
	bool bad_data = false;
        mesh->node_get_coordinates(cnode[k], &xref); // this should not fail
	for (int i = 0; i < spdim; i++)
	  if (x0[k][i] != xref[i]) {
	    bad_data = true;
	    bad_cells0.push_back(j);
	    break;
	  }
	if (bad_data)
	  break;
      }
    } catch (...) {
      bad_cells0.push_back(j);
    }
  }

  int status = 0;

  if (!bad_cells0.empty()) {
    os << "ERROR: bad cell_get_coordinates data for cells:";
    write_list(bad_cells0, MAX_OUT);
    status = -1;
  }

  return status;
}

// Check that face_get_coordinates successfully returns data for all
// faces, and that this data is identical to that returned by
// composing node_get_coordinates with face_get_nodes.  This assumes
// that the std::vector form of face_get_nodes and pointer-based form
// of node_get_coordinates have been verified to return valid data. A
// negative return value signals a terminal error.


int MeshAudit::check_face_to_coordinates() const
{
  os << "Checking face_get_coordinates ..." << endl;

  Entity_ID_List fnode;
  Entity_ID_List bad_faces0;
  Entity_ID_List bad_faces1;
  int spdim = mesh->space_dimension();

  for (unsigned int j = 0; j < nface; ++j) {
    try {
      vector<Point> x0;
      bool bad_data = false;
      mesh->face_get_coordinates(j, &x0); // this may fail
      mesh->face_get_nodes(j, &fnode); // this should not fail
      for (int k = 0; k < fnode.size(); ++k) {
	Point xref(spdim);
        mesh->node_get_coordinates(fnode[k], &xref); // this should not fail
	for (int i = 0; i < spdim; i++)
	  if (x0[k][i] != xref[i]) {
	    bad_data = true;
	    bad_faces0.push_back(j);
	    break;
	  }
	if (bad_data)
	  break;
      }
    } catch (...) {
      bad_faces0.push_back(j);
    }
  }

  int status = 0;

  if (!bad_faces0.empty()) {
    os << "ERROR: bad face_get_coordinates data for faces:";
    write_list(bad_faces0, MAX_OUT);
    status = -1;
  }

  return status;
}

// The cells must not be degenerate, either topologically (repeated
// node index) or geometrically (with coincident nodes).


int MeshAudit::check_cell_geometry() const
{
  os << "Checking cell geometry ..." << endl;

  Point centroid;
  double hvol;
  Entity_ID_List bad_cells;

  for (Entity_ID j = 0; j < ncell; ++j) {
    hvol = mesh->cell_volume(j);
      
    if (hvol <= 0.0) bad_cells.push_back(j);
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
  return check_maps(mesh->node_epetra_map(false), mesh->node_epetra_map(true));
}

int MeshAudit::check_face_maps() const
{
  os << "Checking owned and overlap face maps ..." << endl;
  return check_maps(mesh->face_epetra_map(false), mesh->face_epetra_map(true));
}

int MeshAudit::check_cell_maps() const
{
  os << "Checking owned and overlap cell maps ..." << endl;
  return check_maps(mesh->cell_epetra_map(false), mesh->cell_epetra_map(true));
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

  if (comm.NumProc() == 1)
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
      if (pids[j] < 0 || pids[j] == comm.MyPID()) bad_map = true;
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
  os << "Checking node_get_coordinates ghost data ..." << endl;

  int spdim = mesh->space_dimension();

  const Epetra_Map &node_map_own = mesh->node_epetra_map(false);
  const Epetra_Map &node_map_use = mesh->node_epetra_map(true);

  int nnode_own = node_map_own.NumMyElements();
  int nnode_use = node_map_use.NumMyElements();

  Point coord(spdim);
  Entity_ID_List bad_nodes;

  Epetra_MultiVector coord_use(node_map_use,spdim);
  double **data;
  coord_use.ExtractView(&data);
  Epetra_MultiVector coord_own(View, node_map_own, data, spdim);

  for (unsigned int j = 0; j < nnode_own; ++j) {
    mesh->node_get_coordinates(j, &coord);
    for (int k = 0; k < spdim; ++k) coord_own[k][j] = coord[k];
  }

  Epetra_Import importer(node_map_use, node_map_own);
  coord_use.Import(coord_own, importer, Insert);

  for (unsigned int j = nnode_own; j < nnode_use; ++j) {
    mesh->node_get_coordinates(j, &coord);
    bool bad_data = false;
    for (int k = 0; k < spdim; ++k)
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
  os << "Checking face_get_nodes ghost data ..." << endl;

  const Epetra_Map &node_map = mesh->node_epetra_map(true);
  const Epetra_Map &face_map_own = mesh->face_epetra_map(false);
  const Epetra_Map &face_map_use = mesh->face_epetra_map(true);

  int nface_own = face_map_own.NumMyElements();
  int nface_use = face_map_use.NumMyElements();

  Entity_ID_List fnode;
  Entity_ID_List bad_faces, bad_faces1, bad_faces2;

  // Create a matrix of the GIDs for all owned faces.

  int maxnodes = 0;
  for (unsigned int j = 0; j < nface_own; ++j) {
    mesh->face_get_nodes(j, &fnode);
    maxnodes = fnode.size() > maxnodes ? fnode.size() : maxnodes;
  }

  Epetra_IntSerialDenseMatrix gids(nface_use,maxnodes); // no Epetra_IntMultiVector :(
  for (unsigned int j = 0; j < nface_own; ++j) {
    mesh->face_get_nodes(j, &fnode);
    for (int k = 0; k < fnode.size(); ++k)
      gids(j,k) = node_map.GID(fnode[k]);
    for (int k = fnode.size(); k < maxnodes; ++k)
      gids(j, k) = 0;
  }


  // ????????????????????????????
  // replaced 4 with maxnodes - Seems like it should work

  // Import these GIDs to all used faces; sets values on ghost faces.
  Epetra_Import importer(face_map_use, face_map_own);
  for (int k = 0; k < maxnodes; ++k) {
    Epetra_IntVector kgids_own(View, face_map_own, gids[k]);
    Epetra_IntVector kgids_use(View, face_map_use, gids[k]);
    kgids_use.Import(kgids_own, importer, Insert);
  }

  // Compare the ghost face GIDs against the reference values just computed.
  for (unsigned int j = nface_own; j < nface_use; ++j) {
    mesh->face_get_nodes(j, &fnode);
    bool bad_data = false;
    for (int k = 0; k < maxnodes; ++k)
      if (node_map.GID(fnode[k]) != gids(j,k)) bad_data = true;
    if (bad_data) {
      // Determine just how bad the data is.
      Entity_ID_List fnode_ref(maxnodes);
      for (int k = 0; k < maxnodes; ++k) {
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

  int status = 0;

  if (!bad_faces.empty()) {
    os << "ERROR: found bad data for ghost faces:";
    write_list(bad_faces, MAX_OUT);
    status = 1;
  }

  if (!bad_faces2.empty()) {
    os << "WARNING: found ghost faces that are not exact copies of their master:";
    write_list(bad_faces2, MAX_OUT);
    //status = 1; // don't flag this as an error for now
  }

  return status;
}

// Check that ghost cells are exact copies of their master.  This means
// that the the GIDs of the nodes defining the cell are the same, including
// their order (face orientation).

// CHECK CAREFULLY

int MeshAudit::check_cell_to_nodes_ghost_data() const
{
  os << "Checking cell_get_nodes ghost data ..." << endl;

  const Epetra_Map &node_map = mesh->node_epetra_map(true);
  const Epetra_Map &cell_map_own = mesh->cell_epetra_map(false);
  const Epetra_Map &cell_map_use = mesh->cell_epetra_map(true);

  int ncell_own = cell_map_own.NumMyElements();
  int ncell_use = cell_map_use.NumMyElements();

  Entity_ID_List cnode;
  Entity_ID_List bad_cells;

  int maxnodes = 0;
  for (unsigned int j = 0; j < ncell_own; ++j) {
    mesh->cell_get_nodes(j, &cnode);
    maxnodes = (cnode.size() > maxnodes) ? cnode.size() : maxnodes;
  }
  Epetra_IntSerialDenseMatrix gids(ncell_use,maxnodes); // no Epetra_IntMultiVector :(

  // Create a matrix of the GIDs for all owned cells.
  for (unsigned int j = 0; j < ncell_own; ++j) {
    mesh->cell_get_nodes(j, &cnode);
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
  for (unsigned int j = ncell_own; j < ncell_use; ++j) {
    mesh->cell_get_nodes(j, &cnode);
    bool bad_data = false;
    for (int k = 0; k < cnode.size(); ++k)
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


// Check that ghost cells reference the same faces, and in the same
// order, as their master.  This is part of a ghost being an exact
// copies the master.  Even if the preceding ghost checks pass it is
// still possible for this one to fail.  In this case the ghost would
// reference a different face (by GID) but that face would reference
// the same nodes as the correct face.  So the two faces would be
// geometrically identical, including orientation, but be distinct.

int MeshAudit::check_cell_to_faces_ghost_data() const
{
  os << "Checking cell_get_face ghost data ..." << endl;

  const Epetra_Map &face_map = mesh->face_epetra_map(true);
  const Epetra_Map &cell_map_own = mesh->cell_epetra_map(false);
  const Epetra_Map &cell_map_use = mesh->cell_epetra_map(true);

  int ncell_own = cell_map_own.NumMyElements();
  int ncell_use = cell_map_use.NumMyElements();

  Entity_ID_List cface;
  Entity_ID_List bad_cells;

  // Create a matrix of the GIDs for all owned cells.
  int maxfaces = 0;
  for (unsigned int j = 0; j < ncell_own; ++j) {
    mesh->cell_get_faces(j, &cface);
    maxfaces = (cface.size() > maxfaces) ? cface.size() : maxfaces;
  }

  Epetra_IntSerialDenseMatrix gids(ncell_use,maxfaces); // no Epetra_IntMultiVector :(

  for (unsigned int j = 0; j < ncell_own; ++j) {
    mesh->cell_get_faces(j, &cface);
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
  for (unsigned int j = ncell_own; j < ncell_use; ++j) {
    mesh->cell_get_faces(j, &cface);
    bool bad_data = false;
    for (int k = 0; k < cface.size(); ++k)
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

////////////////////////////////////////////////////////////////////////////////
//
// TESTS OF SET DATA
//

int MeshAudit::check_node_sets() const
{
  os << "Checking node sets ..." << endl;
  return check_sets(NODE, mesh->node_epetra_map(false), mesh->node_epetra_map(true));
}

int MeshAudit::check_face_sets() const
{
  os << "Checking face sets ..." << endl;
  return check_sets(FACE, mesh->face_epetra_map(false), mesh->face_epetra_map(true));
}

int MeshAudit::check_cell_sets() const
{
  os << "Checking cell sets ..." << endl;
  return check_sets(CELL, mesh->cell_epetra_map(false), mesh->cell_epetra_map(true));
}

int MeshAudit::check_sets(Entity_kind kind,
                          const Epetra_Map &map_own, const Epetra_Map &map_use) const
{
  int ierr(0);

  // Basic sanity check on set IDs.
  int ierr_loc = check_set_ids(kind);
  comm.MaxAll(&ierr_loc, &ierr, 1);
  if (ierr) return 1;

  // Check set IDs are same across all processes.
  ierr = check_set_ids_same(kind);
  if (ierr) return 1;

  int status_loc = 0, status = 0; // overall status of the test

  // Additional tests; if these fail we can still continue.
  ierr = check_valid_set_id(kind);
  if (ierr) status_loc = 1;

  // Now get the verified list of set IDs.
  int nset = mesh->num_sets(kind);
  Set_ID_List sids(nset);
  mesh->get_set_ids(kind, &sids);

  for (int n = 0; n < nset; ++n) {
    os << "  Checking set ID=" << sids[n] << " ..." << endl;

    // Basic sanity checks of the owned and used sets.
    int ierr1 = check_get_set(sids[n], kind, OWNED, map_own);
    int ierr2 = check_get_set(sids[n], kind, USED,  map_use);

    // If anyone failed, everyone bails on further tests of this set.
    if (ierr1 != 0 || ierr2 != 0) status_loc = ierr1 = 1;
    comm.MaxAll(&ierr1, &ierr, 1);
    if (ierr) continue;

    // Verify the used set relates correctly to the owned set.
    ierr = check_used_set(sids[n], kind, map_own, map_use);
    if (ierr) status_loc = 1;

    // OUGHT TO DO TESTING OF THE GHOST SETS
  }

  comm.MaxAll(&status_loc, &status, 1);
  return status;
}

// Basic sanity check for set IDs: check that each process can successfully
// get the vector of set IDs, and that the vector contains no duplicates.
// This test runs independently on each process, and returns a per-process
// pass/fail result.

int MeshAudit::check_set_ids(Entity_kind kind) const
{
  // Get the number of sets.
  int nset;
  try {
    nset = mesh->num_sets(kind); // this may fail
  } catch (...) {
    os << "ERROR: caught exception from num_sets()" << endl;
    return 1;
  }

  // Get the vector of set IDs.
  Set_ID_List sids(nset, UINT_MAX);
  try {
    mesh->get_set_ids(kind, &sids); // this may fail
  } catch (...) {
    os << "ERROR: caught exception from get_set_ids()" << endl;
    return 1;
  }

  // Check to see that set ID values were actually assigned.  This assumes
  // UINT_MAX is not a valid set ID.  This is a little iffy; perhaps 0 should
  // be declared as an invalid set ID instead (the case for ExodusII), or
  // perhaps we should just skip this check.
  bool bad_data = false;
  for (int j = 0; j < nset; ++j)
    if (sids[j] == UINT_MAX) bad_data = true;
  if (bad_data) {
    os << "ERROR: get_set_ids() failed to set all values" << endl;
    return 1;
  }

  // Verify that the vector of set IDs contains no duplicates.
  if (!distinct_values(sids)) {
    os << "ERROR: get_set_ids() returned duplicate IDs" << endl;
    // it would be nice to output the duplicates
    return 1;
  }

  return 0;
}

// Check that each process gets the exact same vector of set IDs.
// This is a collective test, returning a collective pass/fail result.

int MeshAudit::check_set_ids_same(Entity_kind kind) const
{
  if (comm.NumProc() > 1) {

    int status_loc = 0, status = 0;

    // Check the number of sets are the same,
    // using the number on process 0 as the reference.
    int nset = mesh->num_sets(kind);
    comm.Broadcast(&nset, 1, 0);
    if (nset != mesh->num_sets(kind)) {
      os << "ERROR: inconsistent num_sets() value" << endl;
      status_loc = 1;
    }
    comm.MaxAll(&status_loc, &status, 1);
    if (status != 0) return 1;

    // Broadcast the set IDs on processor 0.
    Set_ID_List sids(nset, UINT_MAX);
    mesh->get_set_ids(kind, &sids);
    int *sids0 = new int[nset];
    for (int j = 0; j < nset; ++j) sids0[j] = sids[j];
    comm.Broadcast(sids0, nset, 0);

    // Check the set IDs, using the vector on process 0 as the reference.
    bool bad_data = false;
    for (int j = 0; j < nset; ++j)
      if (sids[j] != sids0[j]) bad_data = true;
    if (bad_data) {
      os << "ERROR: get_set_ids() returned inconsistent values" << endl;
      status_loc = 1;
    }
    delete [] sids0;
    comm.MaxAll(&status_loc, &status, 1);
    if (status != 0) return 1;
  }

  return 0;
}

// Check that valid_set_id() returns correct results.  This test runs
// independently on each process and returns a per-process pass/fail result.

int MeshAudit::check_valid_set_id(Entity_kind kind) const
{
  int status = 0;
  int nset = mesh->num_sets(kind); // this should not fail
  Set_ID_List sids(nset);
  mesh->get_set_ids(kind, &sids); // this should not fail
  bool bad_data = false;
  Set_ID_List bad_sids;
  for (int j = 0; j < nset; ++j)
    if (!mesh->valid_set_id(sids[j], kind)) bad_sids.push_back(sids[j]);
  if (!bad_sids.empty()) {
    os << "ERROR: valid_set_id() returned false for valid set IDs:";
    write_list(bad_sids, MAX_OUT);
    status = 1;
  }
  // WE REALLY SHOULD ALSO CHECK THE RESULT FOR A FEW INVALID SET IDS.
  return status;
}


// Basic sanity check on set values: no duplicates, and all LID values belong
// to the map.  This test runs independently on each process and returns a
// per-process pass/fail result.

int MeshAudit::check_get_set(unsigned int sid, Entity_kind kind,
                             Parallel_type category, const Epetra_Map &map) const
{
  // Get the size of the set.
  int n;
  try {
    n = mesh->get_set_size(sid, kind, category); // this may fail
  } catch (...) {
    os << "  ERROR: caught exception from get_set_size()" << endl;
    return 1;
  }

  // Get the set.
  vector<unsigned int> set;
  try {
    mesh->get_set_entities(sid, kind, category, &set);  // this may fail
  } catch (...) {
    os << "  ERROR: caught exception from get_set()" << endl;
    return 1;
  }

  // Check that all values were assigned.
  bool bad_data = false;
  for (int j = 0; j < set.size(); ++j)
    if (set[j] == UINT_MAX) bad_data = true;
  if (bad_data) {
    os << "  ERROR: not all values assigned by get_set()" << endl;
    return 1;
  }

  // Check that the LIDs in the set belong to the map.
  vector<unsigned int> bad_LIDs;
  for (int j = 0; j < set.size(); ++j)
    if (!map.MyLID(set[j])) bad_LIDs.push_back(set[j]);
  if (!bad_LIDs.empty()) {
    os << "  ERROR: set contains invalid LIDs:";
    write_list(bad_LIDs, MAX_OUT);
    return 1;
  }

  // Check that there are no duplicates in the set.
  if (!distinct_values(set)) {
    os << "  ERROR: set contains duplicate LIDs." << endl;
    // it would be nice to output the duplicates
    return 1;
  }

  return 0;
}


// The correct used set is completely determined by the owned set.  This test
// verifies that the used set is what it should be, considering the owned set
// as definitive.  Note that we do not require the vector of used set LIDs to
// extend the vector of owned set LIDs.  The values in each list can be
// presented in any order.

int MeshAudit::check_used_set(unsigned int sid, Entity_kind kind,
                              const Epetra_Map &map_own, const Epetra_Map &map_use) const
{
  if (comm.NumProc() == 1) {

    // In serial, the owned and used sets should be identical.

    int n = mesh->get_set_size(sid, kind, OWNED);
    vector<unsigned int> set_own;
    mesh->get_set_entities(sid, kind, OWNED, &set_own);

    // Set sizes had better be the same.
    if (mesh->get_set_size(sid, kind, USED) != set_own.size()) {
      os << "  ERROR: owned and used set sizes differ" << endl;
      return 1;
    }

    // Verify that the two sets are identical.
    vector<unsigned int> set_use;
    mesh->get_set_entities(sid, kind, USED, &set_use);
    bool bad_data = false;
    for (int j = 0; j < n; ++j)
      if (set_use[j] != set_own[j]) bad_data = true;
    if (bad_data) {
      os << "  ERROR: owned and used sets differ" << endl;
      return 1;
    }

    return 0;

  } else {

    int n = mesh->get_set_size(sid, kind, OWNED);
    vector<unsigned int> set_own;
    mesh->get_set_entities(sid, kind, OWNED, &set_own);

    n = mesh->get_set_size(sid, kind, USED);
    vector<unsigned int> set_use(n);
    mesh->get_set_entities(sid, kind, USED,  &set_use);

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

    int status_loc = 0, status = 0;

    // Check for negative tag values;
    // these mark used LIDs that shouldn't be in the set but are.
    vector<unsigned int> bad_LIDs;
    for (int j = 0; j < set_use.size(); ++j)
      if (tag_use[j] < 0) bad_LIDs.push_back(j);
    if (!bad_LIDs.empty()) {
      os << "  ERROR: found used LIDs that belong to the set but shouldn't:";
      write_list(bad_LIDs, MAX_OUT);
      status_loc = 1;
    }

    // Check for positive tag values;
    // these mark used LIDs that should be in the set but aren't.
    bad_LIDs.resize(0);
    for (int j = 0; j < set_own.size(); ++j)
      if (tag_use[j] > 0) bad_LIDs.push_back(j);
    if (!bad_LIDs.empty()) {
      os << "  ERROR: found used LIDs that should belong to set but don't:";
      write_list(bad_LIDs, MAX_OUT);
      status_loc = 1;
    }

    comm.MaxAll(&status_loc, &status, 1);
    return status;
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

int MeshAudit::check_face_partition() const
{
  os << "Checking face partition ..." << endl;
  
  // Mark all the faces contained by owned cells.
  bool owned[nface];
  for (int j = 0; j < nface; ++j) owned[j] = false;
  Entity_ID_List cface;
  for (unsigned int j = 0; j < mesh->cell_epetra_map(false).NumMyElements(); ++j) {
    mesh->cell_get_faces(j, &cface);
    for (int k = 0; k < cface.size(); ++k) owned[cface[k]] = true;
  }
  
  // Verify that every owned face has been marked as belonging to an owned cell.
  Entity_ID_List bad_faces;
  for (unsigned int j = 0; j < mesh->face_epetra_map(false).NumMyElements(); ++j)
    if (!owned[j]) bad_faces.push_back(j);
  
  if (!bad_faces.empty()) {
    os << "ERROR: found orphaned owned faces:";
    write_list(bad_faces, MAX_OUT);
    os << "       Process doesn't own either of the cells sharing the face." << endl;
    return 1;
  }
  
  return 0;
}


int MeshAudit::check_node_partition() const
{
  os << "Checking node partition ..." << endl;
  
  // Mark all the nodes contained by owned cells.
  bool owned[nnode];
  for (int j = 0; j < nnode; ++j) owned[j] = false;
  Entity_ID_List cnode;
  for (unsigned int j = 0; j < mesh->cell_epetra_map(false).NumMyElements(); ++j) {
    mesh->cell_get_nodes(j, &cnode);
    for (int k = 0; k < cnode.size(); ++k) owned[cnode[k]] = true;
  }
  
  // Verify that every owned node has been marked as belonging to an owned cell.
  Entity_ID_List bad_nodes;
  for (unsigned int j = 0; j < mesh->node_epetra_map(false).NumMyElements(); ++j)
    if (!owned[j]) bad_nodes.push_back(j);
  
  if (!bad_nodes.empty()) {
    os << "ERROR: found orphaned owned nodes:";
    write_list(bad_nodes, MAX_OUT);
    os << "       Process doesn't own any of the cells containing the node." << endl;
    return 1;
  }
  
  return 0;
}

// Returns true if the values in the list are distinct -- no repeats.

bool MeshAudit::distinct_values(const Entity_ID_List &list) const
{
  Entity_ID_List copy(list);
  sort(copy.begin(), copy.end());
  return (adjacent_find(copy.begin(),copy.end()) == copy.end());
}


// Returns 1 if the face node lists fnode1 and fnode2 describe the same face
// with the same orientation.  Returns -1 if the lists describe the same face
// but with opposite orientations.  Returns 0 if the lists describe different
// faces.  Implicitly assumes non-degenerate faces; the results are not
// reliable otherwise.

int MeshAudit::same_face(const Entity_ID_List fnode1, const Entity_ID_List fnode2) const
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


  for (i = 1; i < nn; ++i)
    if (fnode1[(n+i)%nn] != fnode2[i]) break;
  if (i == nn) return 1;  // they match

  // Modify the permutation to reverse the orientation of fnode1.

  for (i = 1; i < nn; ++i)
    if (fnode1[(n-i+nn)%nn] != fnode2[i]) break;
  if (i == nn) return -1;   // matched nodes but orientation is reversed

  return 0; // different faces
}


void MeshAudit::write_list(const Entity_ID_List &list, unsigned int max_out) const
{
  int num_out = min((unsigned int) list.size(), max_out);
  for (int i = 0; i < num_out; ++i) os << " " << list[i];
  if (num_out < list.size()) os << " [" << list.size()-num_out << " items omitted]";
  os << endl;
}


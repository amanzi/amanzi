/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "MeshLogicalAudit.hh"

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

namespace Amanzi {
namespace AmanziMesh {


MeshLogicalAudit::MeshLogicalAudit(const Teuchos::RCP<const MeshHost>& mesh_, std::ostream& os_)
  : mesh(mesh_),
    comm_(mesh_->getComm()),
    MyPID(mesh_->getComm()->MyPID()),
    nface(mesh_->getMap(Entity_kind::FACE, true).NumMyElements()),
    ncell(mesh_->getMap(Entity_kind::CELL, true).NumMyElements()),
    os(os_),
    MAX_OUT(5)
{
  create_test_dependencies();
}

// Verify runs all the tests in an order that respects the dependencies
// between tests.  If a particular test fails, all other tests that have it
// as a pre-requisite are skipped.  It is important that each test return
// a collective fail/pass result in parallel, so that all processes proceed
// through the tests in lockstep.

int
MeshLogicalAudit::Verify() const
{
  int ok(0), v;

  while ((v = FindAnyRoot()) != -1) {
    os << "Checking " << vertices_[v].name << " ..." << std::endl;
    bool flag = (this->*(vertices_[v].test))();
    if (!flag) {
      vertices_[v].status = Status_kind::GOOD;
    } else {
      os << "  Test failed!" << std::endl;
      vertices_[v].status = Status_kind::FAIL;
      SkipBranches(v, os);
      ok = 1;
    }
  }
  return ok;
}

// This creates the dependency graph for the individual tests.  Adding a new
// test to the graph is a simple matter of adding a new vertex and edge.
void
MeshLogicalAudit::create_test_dependencies()
{
  AddVertex("entity_counts", &MeshLogicalAudit::check_entity_counts); // test01
  auto test04 = AddVertex("cell_to_faces", &MeshLogicalAudit::check_cell_to_faces);

  auto test05 = AddVertex("face references by cells", &MeshLogicalAudit::check_face_refs_by_cells);
  AddEdge(test04, test05);

  auto test09 = AddVertex("cell_to_face_dirs", &MeshLogicalAudit::check_cell_to_face_dirs);

  auto test08 =
    AddVertex("cell_to_face_consistency", &MeshLogicalAudit::check_faces_cell_consistency);
  AddEdge(test05, test08);

  // cell degeneracy test
  auto test10 =
    AddVertex("topological non-degeneracy of cells", &MeshLogicalAudit::check_cell_degeneracy);
  AddEdge(test04, test10);

  // cell topology/geometry test
  auto test15 = AddVertex("cell geometry", &MeshLogicalAudit::check_cell_geometry);
  AddEdge(test10, test15);

  auto test16 = AddVertex("face geometry", &MeshLogicalAudit::check_face_geometry);
  AddEdge(test09, test16);

  auto test14 =
    AddVertex("face-cell-bisector geometry", &MeshLogicalAudit::check_cell_face_bisector_geometry);
  AddEdge(test08, test14);

  auto test17 = AddVertex("owned and overlap face maps", &MeshLogicalAudit::check_face_maps);
  auto test18 = AddVertex("owned and overlap cell maps", &MeshLogicalAudit::check_cell_maps);
  auto test22 =
    AddVertex("cell_to_faces ghost data", &MeshLogicalAudit::check_cell_to_faces_ghost_data);
  AddEdge(test04, test22);
  AddEdge(test17, test22);
  AddEdge(test18, test22);

  // partition tests
  auto test32 = AddVertex("face partition", &MeshLogicalAudit::check_face_partition);
  AddEdge(test17, test32);
  AddEdge(test18, test32);
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

bool
MeshLogicalAudit::check_entity_counts() const
{
  int n, nref;
  bool error = false;

  // Check the number of owned faces.
  n = mesh->getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
  nref = mesh->getMap(Entity_kind::FACE, false).NumMyElements();
  if (n != nref) {
    os << "ERROR: num_entities(FACE,OWNED)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  // Check the number of used faces.
  n = mesh->getNumEntities(Entity_kind::FACE, Parallel_kind::ALL);
  nref = mesh->getMap(Entity_kind::FACE, true).NumMyElements();
  if (n != nref) {
    os << "ERROR: num_entities(FACE,ALL)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  // Check the number of owned cells.
  n = mesh->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  nref = mesh->getMap(Entity_kind::CELL, false).NumMyElements();
  if (n != nref) {
    os << "ERROR: num_entities(CELL,OWNED)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  // Check the number of used cells.
  n = mesh->getNumEntities(Entity_kind::CELL, Parallel_kind::ALL);
  nref = mesh->getMap(Entity_kind::CELL, true).NumMyElements();
  if (n != nref) {
    os << "ERROR: num_entities(CELL,ALL)=" << n << "; should be " << nref << std::endl;
    error = true;
  }

  return global_any(error);
}
// Returns true if the values in the list are distinct -- no repeats.

bool
MeshLogicalAudit::distinct_values(const MeshHost::cEntity_ID_View& list) const
{
  MeshHost::Entity_ID_View copy;
  copy.fromConst(list);
  sort(copy.begin(), copy.end());
  return (adjacent_find(copy.begin(), copy.end()) == copy.end());
}


// Check that cell_to_faces successfully returns valid references to local
// faces.  Here the std::vector accessor is used.  The consistency of
// the alternative accessors is checked elsewhere.  A nonzero return value
// signals an error, and further tests using its data should be avoided.

bool
MeshLogicalAudit::check_cell_to_faces() const
{
  Entity_ID_List bad_cells, bad_cells1;
  MeshHost::cEntity_ID_View cface;

  for (Entity_ID j = 0; j < ncell; ++j) {
    try {
      mesh->getCellFaces(j, cface); // this may fail
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

bool
MeshLogicalAudit::check_face_refs_by_cells() const
{
  MeshHost::cEntity_ID_View cface;
  vector<unsigned int> refs(nface, 0);

  for (Entity_ID j = 0; j < ncell; ++j) {
    mesh->getCellFaces(j, cface);
    for (int k = 0; k < cface.size(); ++k) (refs[cface[k]])++;
  }

  Entity_ID_List free_faces;
  Entity_ID_List bad_faces;

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


// Check the consistency of face_get_cells and cell_get_faces.
bool
MeshLogicalAudit::check_faces_cell_consistency() const
{
  MeshHost::cEntity_ID_View cface;
  MeshHost::cEntity_ID_View fcell;

  Entity_ID_List bad_cells;
  for (Entity_ID j = 0; j < ncell; ++j) {
    mesh->getCellFaces(j, cface);
    for (int k = 0; k < cface.size(); ++k) {
      mesh->getFaceCells(cface[k], fcell);
      if (std::find(fcell.begin(), fcell.end(), j) == fcell.end()) { bad_cells.push_back(j); }
    }
  }

  Entity_ID_List bad_faces;
  for (Entity_ID j = 0; j < nface; ++j) {
    mesh->getFaceCells(j, fcell);
    for (int k = 0; k < fcell.size(); ++k) {
      mesh->getCellFaces(fcell[k], cface);
      if (std::find(cface.begin(), cface.end(), j) == cface.end()) { bad_faces.push_back(j); }
    }
  }

  bool error = false;
  if (!bad_faces.empty()) {
    os << "ERROR: found inconsistent faces:";
    write_list(bad_faces, MAX_OUT);
    error = true;
  }

  if (!bad_cells.empty()) {
    os << "ERROR: found inconsistent cells:";
    write_list(bad_cells, MAX_OUT);
    error = true;
  }

  return global_any(error);
}

// Check that cell_to_face_dirs successfully returns data for all cells and
// that the values are either +1 or -1. The std::vector-based method is used;
// the consistency of the alternative methods is checked elsewhere.  If this
// test fails, further tests using this data should be avoided.

bool
MeshLogicalAudit::check_cell_to_face_dirs() const
{
  MeshHost::cEntity_ID_View faces;
  MeshHost::cDirection_View fdirs;
  Entity_ID_List bad_cells, bad_cells_exc;

  for (Entity_ID j = 0; j < ncell; ++j) {
    try {
      mesh->getCellFacesAndDirs(j, faces, &fdirs); // this may fail
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

bool
MeshLogicalAudit::check_cell_degeneracy() const
{
  os << "Checking cells for topological degeneracy ..." << std::endl;

  MeshHost::cEntity_ID_View cface;
  Entity_ID_List bad_cells;

  for (Entity_ID j = 0; j < ncell; ++j) {
    mesh->getCellFaces(j, cface); // should not fail
    if (!distinct_values(cface)) bad_cells.push_back(j);
  }

  bool error = false;

  if (!bad_cells.empty()) {
    os << "ERROR: found topologically degenerate cells:";
    write_list(bad_cells, MAX_OUT);
    error = true;
  }

  return global_any(error);
}

// The cells must not be degenerate, with zero volume
bool
MeshLogicalAudit::check_cell_geometry() const
{
  os << "Checking cell geometry ..." << std::endl;
  AmanziGeometry::Point centroid;
  double hvol;
  Entity_ID_List bad_cells;

  for (Entity_ID j = 0; j < ncell; ++j) {
    hvol = mesh->getCellVolume(j);
    if (hvol <= 1.e-10) bad_cells.push_back(j);
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

// The faces must not be degenerate
bool
MeshLogicalAudit::check_face_geometry() const
{
  os << "Checking face geometry ..." << std::endl;
  double hvol;
  Entity_ID_List bad_faces;

  for (Entity_ID j = 0; j < nface; ++j) {
    hvol = mesh->getFaceArea(j);
    if (hvol <= 0.0) bad_faces.push_back(j);
  }

  bool error = false;

  if (!bad_faces.empty()) {
    os << "ERROR: found faces with non-positive area:";
    write_list(bad_faces, MAX_OUT);
    error = true;
  }

  return global_any(error);
}


// The faces must not be degenerate
bool
MeshLogicalAudit::check_cell_face_bisector_geometry() const
{
  os << "Checking cell-to-face bisector geometry ..." << std::endl;
  Entity_ID_List bad_cells;

  MeshHost::cEntity_ID_View cface;
  MeshHost::cPoint_View bisectors;
  for (Entity_ID j = 0; j < ncell; ++j) {
    mesh->getCellFacesAndBisectors(j, cface, &bisectors);
    if (cface.size() != bisectors.size()) {
      bad_cells.push_back(j);
      continue;
    }
    for (int i = 0; i != cface.size(); ++i) {
      if (AmanziGeometry::norm(bisectors[i]) < 1.e-10) { bad_cells.push_back(j); }
    }
  }

  bool error = false;

  if (!bad_cells.empty()) {
    os << "ERROR: found cells with bad bisector data:";
    write_list(bad_cells, MAX_OUT);
    error = true;
  }

  return global_any(error);
}


bool
MeshLogicalAudit::check_face_maps() const
{
  return check_maps(mesh->getMap(Entity_kind::FACE, false), mesh->getMap(Entity_kind::FACE, true));
}

bool
MeshLogicalAudit::check_cell_maps() const
{
  return check_maps(mesh->getMap(Entity_kind::CELL, false), mesh->getMap(Entity_kind::CELL, true));
}

bool
MeshLogicalAudit::check_maps(const Epetra_Map& map_own, const Epetra_Map& map_use) const
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
  if (!map_own.UniqueGIDs()) {
    os << "ERROR: owned map is not 1-to-1" << std::endl;
    error = true;
  }

  // Check that the global ID space is contiguously divided (but not
  // necessarily uniformly) across all processers

  std::vector<int> owned_GIDs(map_own.NumMyElements());
  for (int i = 0; i < map_own.NumMyElements(); i++) owned_GIDs[i] = map_own.GID(i);
  std::sort(owned_GIDs.begin(), owned_GIDs.end());

  for (int i = 0; i < map_own.NumMyElements() - 1; i++) {
    int diff = owned_GIDs[i + 1] - owned_GIDs[i];
    if (diff > 1) {
      os << "ERROR: owned map is not contiguous" << std::endl;
      os << "Global IDs jump from " << owned_GIDs[i] << " to " << owned_GIDs[i + 1] << std::endl;
      error = true;
    }
  }

  error = global_any(error);
  if (error) return error;

  if (comm_->NumProc() == 1) {
    // Serial or 1-process MPI

    if (!map_use.SameAs(map_own)) {
      os << "ERROR: the overlap map differs from the owned map (single process)." << std::endl;
      error = true;
    }

    return global_any(error);

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
      os << "ERROR: overlap map does not extend the owned map." << std::endl;
      error = true;
    }

    error = global_any(error);
    if (error) return error;

    // Verify that the overlap indices are owned by other processes.
    int* gids = new int[num_ovl];
    int* pids = new int[num_ovl];
    int* lids = new int[num_ovl];
    for (int j = 0; j < num_ovl; ++j) gids[j] = map_use.GID(j + num_own);
    map_own.RemoteIDList(num_ovl, gids, pids, lids);
    bad_map = false;
    for (int j = 0; j < num_ovl; ++j)
      if (pids[j] < 0 || pids[j] == comm_->MyPID()) bad_map = true;
    if (bad_map) {
      os << "ERROR: invalid ghosts in overlap map." << std::endl;
      error = true;
    }

    // Look for duplicates among the overlap indices.
    vector<int> ovl_gids(gids, gids + num_ovl);
    sort(ovl_gids.begin(), ovl_gids.end());
    if (adjacent_find(ovl_gids.begin(), ovl_gids.end()) != ovl_gids.end()) {
      os << "ERROR: duplicate ghosts in overlap map." << std::endl;
      error = true;
    }

    delete[] lids;
    delete[] pids;
    delete[] gids;

    return global_any(error);
  }
}


// Check that ghost cells reference the same faces, and in the same
// order, as their master.  This is part of a ghost being an exact
// copies the master.  Even if the preceding ghost checks pass it is
// still possible for this one to fail.  In this case the ghost would
// reference a different face (by GID) but that face would reference
// the same nodes as the correct face.  So the two faces would be
// geometrically identical, including orientation, but be distinct.

bool
MeshLogicalAudit::check_cell_to_faces_ghost_data() const
{
  const Epetra_Map& face_map = mesh->getMap(Entity_kind::FACE, true);
  const Epetra_Map& cell_map_own = mesh->getMap(Entity_kind::CELL, false);
  const Epetra_Map& cell_map_use = mesh->getMap(Entity_kind::CELL, true);

  int ncell_own = cell_map_own.NumMyElements();
  int ncell_use = cell_map_use.NumMyElements();

  MeshHost::cEntity_ID_View cface;
  Entity_ID_List bad_cells;

  // Create a matrix of the GIDs for all owned cells.
  int maxfaces = 0;
  for (Entity_ID j = 0; j < ncell_own; ++j) {
    mesh->getCellFaces(j, cface);
    maxfaces = (cface.size() > maxfaces) ? cface.size() : maxfaces;
  }

  Epetra_IntSerialDenseMatrix gids(ncell_use, maxfaces); // no Epetra_IntMultiVector :(

  for (Entity_ID j = 0; j < ncell_own; ++j) {
    mesh->getCellFaces(j, cface);
    for (int k = 0; k < cface.size(); ++k) gids(j, k) = face_map.GID(cface[k]);
    for (int k = cface.size(); k < maxfaces; ++k) gids(j, k) = 0;
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
    mesh->getCellFaces(j, cface);
    bool bad_data = false;
    for (int k = 0; k < cface.size(); ++k)
      if (face_map.GID(cface[k]) != gids(j, k)) bad_data = true;
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


// SANE PARTITIONING CHECKS.  In parallel the cells, faces and nodes of the
// mesh are each partitioned across the processes, and while in principle
// these partitionings may be completely independent of each other, practical
// considerations lead to certain conditions that reasonable partitions should
// satisfy.  Taking the partitioning of the cells as given, one property that
// should be satisfied by the face and node partitionings is that the process
// that owns a particular face (node) must also own one of the cells containing
// the face (node).

bool
MeshLogicalAudit::check_face_partition() const
{
  // Mark all the faces contained by owned cells.
  std::vector<bool> owned(nface, false);
  MeshHost::cEntity_ID_View cface;
  for (Entity_ID j = 0; j < mesh->getMap(Entity_kind::CELL, false).NumMyElements(); ++j) {
    mesh->getCellFaces(j, cface);
    for (int k = 0; k < cface.size(); ++k) owned[cface[k]] = true;
  }

  // Verify that every owned face has been marked as belonging to an owned cell.
  Entity_ID_List bad_faces;
  for (Entity_ID j = 0; j < mesh->getMap(Entity_kind::FACE, false).NumMyElements(); ++j)
    if (!owned[j]) bad_faces.push_back(j);

  if (!bad_faces.empty()) {
    os << "ERROR: found orphaned owned faces:";
    write_list(bad_faces, MAX_OUT);
    os << "       Process doesn't own either of the cells sharing the face." << std::endl;
    return true;
  }

  return false;
}

// Returns 1 if the face node lists fnode1 and fnode2 describe the same face
// with the same orientation.  Returns -1 if the lists describe the same face
// but with opposite orientations.  Returns 0 if the lists describe different
// faces.  Implicitly assumes non-degenerate faces; the results are not
// reliable otherwise.

int
MeshLogicalAudit::same_face(const MeshHost::Entity_ID_View fnode1,
                            const MeshHost::Entity_ID_View fnode2) const
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
      if (fnode1[(n + i) % nn] != fnode2[i]) break;
    if (i == nn) return 1; // they match

    // Modify the permutation to reverse the orientation of fnode1.

    for (i = 1; i < nn; ++i)
      if (fnode1[(n - i + nn) % nn] != fnode2[i]) break;
    if (i == nn) return -1; // matched nodes but orientation is reversed
  }

  return 0; // different faces
}


void
MeshLogicalAudit::write_list(const Entity_ID_List& list, unsigned int max_out) const
{
  int num_out = min((unsigned int)list.size(), max_out);
  for (int i = 0; i < num_out; ++i) os << " " << list[i];
  if (num_out < list.size()) os << " [" << list.size() - num_out << " items omitted]";
  os << std::endl;
}

bool
MeshLogicalAudit::global_any(bool value) const
{
  int lval = value, gval;
  comm_->MaxAll(&lval, &gval, 1);
  return gval;
}

} // namespace AmanziMesh
} // namespace Amanzi

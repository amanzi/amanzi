/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
// Helper struct for resolving columnar structure of meshes.


#include "MeshColumns.hh"
#include "MeshCache.hh"

namespace Amanzi {
namespace AmanziMesh {


//
// Constructor that guesses how to infer columnar structure.
//
void
MeshColumns::initialize(const MeshCache& mesh)
{
  AMANZI_ASSERT(mesh.getSpaceDimension() == 3);
  AMANZI_ASSERT(mesh.getManifoldDimension() == 3);

  // loop over all boundary faces and look for those whose normal z component
  // is not zero and whose opposing face is below them.
  int col_id = 0;
  Entity_ID_List surface_faces;
  for (const Entity_ID f : mesh.getBoundaryFaces()) {
    auto f_normal = mesh.getFaceNormal(f);
    Entity_ID c = mesh.getFaceCells(f, Parallel_type::ALL)[0];
    if (Impl::orientFace(mesh, f, c) == 1) surface_faces.emplace_back(f);
  }

  initialize(mesh, surface_faces);
}

//
// Constructor that infers columnar structure from a mesh set.
//
void
MeshColumns::initialize(const MeshCache& mesh, const std::vector<std::string>& regions)
{
  AMANZI_ASSERT(mesh.getSpaceDimension() == 3);
  AMANZI_ASSERT(mesh.getManifoldDimension() == 3);

  if (regions.size() == 0) {
    Errors::Message msg("MeshColumns::initialize() called with empty list of regions.");
    Exceptions::amanzi_throw(msg);
  }

  // collect faces in all regions, keeping the collection sorted
  Entity_ID_List surface_faces;
  for (const auto& r : regions) {
    auto r_faces = mesh.getSetEntities(r, Entity_kind::FACE, Parallel_type::ALL);
    for (Entity_ID f : r_faces)
      surface_faces.insert(std::upper_bound(surface_faces.begin(),
              surface_faces.end(),f),f);
  }

  // build columns for each face
  initialize(mesh, surface_faces);
}

void
MeshColumns::initialize(const MeshCache& mesh, const Entity_ID_List& surface_faces)
{
  // figure out the correct size
  // Note, this is done to make life easier for Kokkos
  cells.resize(surface_faces.size());
  faces.resize(surface_faces.size());
  int i = 0;
  for (Entity_ID f : surface_faces) {
    auto size = Impl::countCellsInColumn(mesh, f);
    cells[i].resize(size);
    faces[i++].resize(size+1);
  }

  // fill
  i = 0;
  for (Entity_ID f : surface_faces) {
    buildColumn_(mesh, f, i++);
    if (f < mesh.nfaces_owned) num_columns_owned = i;
  }
  num_columns_all = i;
}


//
// Actually does the work of building the column, starting at f.
//
void
MeshColumns::buildColumn_(const MeshCache& mesh, Entity_ID f, int col)
{
  Parallel_type ptype = mesh.getParallelType(Entity_kind::FACE, f);

  Entity_ID c = mesh.getFaceCells(f, Parallel_type::ALL)[0]; // guaranteed size 1
  if (mesh.getParallelType(Entity_kind::CELL, c) != ptype) {
    Errors::Message msg("MeshColumns::buildColumn() mesh was not partitioned vertically using zoltan_rcb");
    Exceptions::amanzi_throw(msg);
  }
  bool done = false;
  int i = 0;
  while (!done) {
    cells[col][i] = c;
    faces[col][i] = f;

    Entity_ID f_opp = Impl::findDownFace(mesh, c);
    if (mesh.getParallelType(Entity_kind::FACE, f_opp) != ptype) {
      Errors::Message msg("MeshColumns::buildColumn() mesh was not partitioned vertically using zoltan_rcb");
      Exceptions::amanzi_throw(msg);
    }
    AMANZI_ASSERT(f_opp != f);

    Entity_ID c_opp = Impl::findOpposingCell(mesh, c, f_opp);
    if (c_opp < 0) {
      faces[col][i+1] = f_opp;
      done = true;
    } else {
      if (mesh.getParallelType(Entity_kind::CELL, c_opp) != ptype) {
        Errors::Message msg("MeshColumns::buildColumn() mesh was not partitioned vertically using zoltan_rcb");
        Exceptions::amanzi_throw(msg);
      }
      i++;
      f = f_opp;
      c = c_opp;
    }
  }
}

namespace Impl {

int
orientFace(const MeshCache& mesh, const Entity_ID f, const Entity_ID c)
{
  auto normal = mesh.getFaceNormal(f, c);
  normal /= AmanziGeometry::norm(normal);
  if (normal[normal.dim()-1] > 1.e-10) return 1;
  else if (normal[normal.dim()-1] < -1.e-10) return -1;
  return 0;
}


//
// Helper function to find the face opposite of f across c.
//
// Assumes this is well posed... e.g. that there is only one face of c that
// shares no nodes with f.   Returns -1 if this is not possible.
Entity_ID
findDownFace(const MeshCache& mesh, const Entity_ID c)
{
  auto cfaces = mesh.getCellFaces(c);
  // find the one that is down
  for (const auto& other_f : cfaces) {
    if (Impl::orientFace(mesh, other_f, c) == -1) return other_f;
  }
  return -1;
}

//
// Finds the cell in face cells that is not c, returning -1 if not possible.
//
Entity_ID
findOpposingCell(const MeshCache& mesh, const Entity_ID c, const Entity_ID f)
{
  auto fcells =  mesh.getFaceCells(f, Parallel_type::ALL);
  if (fcells.size() == 1) return -1;
  else if (fcells[0] == c) return fcells[1];
  else return fcells[0];
}

//
// Helper function for counting column size
//
std::size_t
countCellsInColumn(const MeshCache& mesh, Entity_ID f)
{
  std::size_t count = 0;
  Entity_ID c = mesh.getFaceCells(f, Parallel_type::ALL)[0]; // guaranteed size 1
  bool done = false;
  int i = 0;
  while (!done) {
    count++;
    Entity_ID f_opp = Impl::findDownFace(mesh, c);
    AMANZI_ASSERT(f_opp != f);
    Entity_ID c_opp = Impl::findOpposingCell(mesh, c, f_opp);
    if (c_opp < 0) {
      done = true;
    } else {
      f = f_opp;
      c = c_opp;
    }
  }
  return count;
}



} // namespace Impl
} // namespace AmanziMesh
} // namespace Amanzi


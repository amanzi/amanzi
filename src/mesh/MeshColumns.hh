/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
// Helper struct for resolving columnar structure of meshes.

#pragma once

#include "AmanziTypes.hh"
#include "MeshDefs.hh"

namespace Amanzi {
namespace AmanziMesh {

class MeshCache;

struct MeshColumns {
 public:
  MeshColumns()
    : num_columns_all(-1),
      num_columns_owned(-1) {}

  //
  // Constructor that guesses how to infer columnar structure.
  //
  void initialize(const MeshCache& mesh);

  //
  // Constructor that infers columnar structure from a mesh set.
  //
  void initialize(const MeshCache& mesh, const std::vector<std::string>& regions);

  //
  // Constructor that sets top faces directly
  //
  void initialize(const MeshCache& mesh, const Entity_ID_List& surface_faces);

  int num_columns_owned;
  int num_columns_all;

  // a vector of columns, each a list of cells in the column, ordered from top
  // down
  RaggedArray<Entity_ID> cells;

  // a vector of columns, each a list of faces in the column, ordered from top
  // down.  Note len(faces[i]) == len(cells[i])+1
  RaggedArray<Entity_ID> faces;

private:
  void buildColumn_(const MeshCache& mesh, const Entity_ID f, int col);

};

namespace Impl {

//
// Helper function to find the direction of f relative to c.
//
// Returns 1 if the up face, -1 if the down face, 0 if the side face.
//
int
orientFace(const MeshCache& mesh, const Entity_ID f, const Entity_ID c);

//
// Helper function to find the face of c pointing downward.
//
// Assumes this is well posed... e.g. that there is only one face of c that
// points down.   Returns -1 if none are found.
Entity_ID findDownFace(const MeshCache& mesh, const Entity_ID c);

//
// Finds the cell in face cells that is not c, returning -1 if not possible.
//
Entity_ID findOpposingCell(const MeshCache& mesh, const Entity_ID c, const Entity_ID f);

//
// Helper function for counting column size
//
std::size_t countCellsInColumn(const MeshCache& mesh, Entity_ID f);

}

} // namespace AmanziMesh
} // namespace Amanzi



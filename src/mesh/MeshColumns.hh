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
#include "MeshUtils.hh"

namespace Amanzi {
namespace AmanziMesh {

template<MemSpace_type>
class MeshCache;

struct MeshColumns {
 public:

  //
  // standard constructor
  //
  MeshColumns()
    : num_columns_all(-1),
      num_columns_owned(-1) {}

  //
  // Constructor that guesses how to infer columnar structure.
  //
  // Note these MUST be constructed on HOST caches.
  //
  void initialize(const MeshCache<MemSpace_type::HOST>& mesh);

  //
  // Constructor that infers columnar structure from a mesh set.
  //
  // Note these MUST be constructed on HOST caches.
  //
  void initialize(const MeshCache<MemSpace_type::HOST>& mesh, const std::vector<std::string>& regions);

  //
  // Constructor that sets top faces directly
  //
  // Note these MUST be constructed on HOST caches.
  //
  void initialize(const MeshCache<MemSpace_type::HOST>& mesh, const Entity_ID_List& surface_faces);

  int num_columns_owned;
  int num_columns_all;

  //
  // View all cells or faces
  //
  template<MemSpace_type MEM>
  KOKKOS_INLINE_FUNCTION
  decltype(auto) // cEntity_ID_View
  getCells(int col) {
    return cells_.getRow<MEM>(col);
  }

  template<MemSpace_type MEM>
  KOKKOS_INLINE_FUNCTION
  decltype(auto) // cEntity_ID_View
  getFaces(int col) {
    return faces_.getRow<MEM>(col);
  }

  //
  // Get one cell or face
  //
  template<MemSpace_type MEM>
  KOKKOS_INLINE_FUNCTION
  Entity_ID getCell(int col, int i) {
    return cells_.get<MEM>(col, i);
  }

  template<MemSpace_type MEM>
  KOKKOS_INLINE_FUNCTION
  Entity_ID getFace(int col, int i) {
    return faces_.get<MEM>(col, i);
  }

 private:
  void buildColumn_(const MeshCache<MemSpace_type::HOST>& mesh, const Entity_ID f, int col);

  public: 
  // a vector of columns, each a list of cells in the column, ordered from top
  // down
  RaggedArray_DualView<Entity_ID> cells_;

  // a vector of columns, each a list of faces in the column, ordered from top
  // down.  Note len(faces[i]) == len(cells[i])+1
  RaggedArray_DualView<Entity_ID> faces_;

};


namespace Impl {
//
// Helper function to find the direction of f relative to c.
//
// Returns 1 if the up face, -1 if the down face, 0 if the side face.
//
int
orientFace(const MeshCache<MemSpace_type::HOST>& mesh, const Entity_ID f, const Entity_ID c);

//
// Helper function to find the face of c pointing downward.
//
// Assumes this is well posed... e.g. that there is only one face of c that
// points down.   Returns -1 if none are found.
Entity_ID findDownFace(const MeshCache<MemSpace_type::HOST>& mesh, const Entity_ID c);

//
// Finds the cell in face cells that is not c, returning -1 if not possible.
//
Entity_ID findOpposingCell(const MeshCache<MemSpace_type::HOST>& mesh, const Entity_ID c, const Entity_ID f);

//
// Helper function for counting column size
//
std::size_t countCellsInColumn(const MeshCache<MemSpace_type::HOST>& mesh, Entity_ID f);

} // namespace Impl


} // namespace AmanziMesh
} // namespace Amanzi



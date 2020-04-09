/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

//
// Base mesh class for Amanzi
//
// Inline implementations
//

#include "Kokkos_Core.hpp"

// Downward Adjacencies
//---------------------
// void
// KOKKOS_INLINE_FUNCTION Mesh::cell_get_faces(const Entity_ID cellid,
// Kokkos::View<Entity_ID*>& faceids) const
//{
//  cell_get_faces_and_dirs(cellid, faceids, NULL, false);
//}


//
// Column information
//-----------------------
inline int
Mesh::num_columns(bool ghosted) const
{
  if (!columns_built_) {
    Errors::Message mesg("num_columns called before calling build_columns");
    Exceptions::amanzi_throw(mesg);
  }
  return ghosted ? column_cells_.size() - 1 :
                   num_owned_cols_; // number of vector of vectors
}


inline
const Entity_ID_List&
Mesh::cells_of_column(const int columnID) const
{
  if (!columns_built_) {
    Errors::Message mesg("cells_of_column called before calling build_columns");
    Exceptions::amanzi_throw(mesg);
  }
  return column_cells_[columnID];
}


inline
const Entity_ID_List&
Mesh::faces_of_column(const int columnID) const
{
  if (!columns_built_) {
    Errors::Message mesg("faces_of_columns called before calling build_columns");
    Exceptions::amanzi_throw(mesg);
  }
  return column_faces_[columnID];
}



inline int
Mesh::column_ID(const Entity_ID cellid) const
{
  if (!columns_built_) {
    Errors::Message mesg("column_ID called before calling build_columns");
    Exceptions::amanzi_throw(mesg);
  }
  return columnID_[cellid];
}


inline Entity_ID
Mesh::cell_get_cell_above(const Entity_ID cellid) const
{
  if (!columns_built_) {
    Errors::Message mesg(
      "cell_get_cell_above called before calling build_columns");
    Exceptions::amanzi_throw(mesg);
  }
  return cell_cellabove_[cellid];
}


inline Entity_ID
Mesh::cell_get_cell_below(const Entity_ID cellid) const
{
  if (!columns_built_) {
    Errors::Message mesg(
      "cell_get_cell_below called before calling build_columns");
    Exceptions::amanzi_throw(mesg);
  }
  return cell_cellbelow_[cellid];
}


inline Entity_ID
Mesh::node_get_node_above(const Entity_ID nodeid) const
{
  if (!columns_built_) {
    Errors::Message mesg(
      "node_get_node_above called before calling build_columns");
    Exceptions::amanzi_throw(mesg);
  }
  return node_nodeabove_[nodeid];
}


//
// Column information
//-----------------------
inline void
Mesh::get_set_entities(const std::string setname, const Entity_kind kind,
                       const Parallel_type ptype,
                       Entity_ID_List& entids) const
{
  std::vector<double> vofs;
  get_set_entities_and_vofs(setname, kind, ptype, entids, &vofs);
}

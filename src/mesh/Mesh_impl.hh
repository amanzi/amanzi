/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//
// Base mesh class for Amanzi
//
// Inline implementations
//


// Downward Adjacencies
//---------------------
inline
void
Mesh::cell_get_faces(const Entity_ID cellid, Entity_ID_List *faceids,
                          const bool ordered) const
{
  cell_get_faces_and_dirs(cellid, faceids, NULL, ordered);
}


//
// Column information
//-----------------------
inline
int
Mesh::num_columns(bool ghosted) const
{
  if (!columns_built_) build_columns_();
  return ghosted ? column_cells_.size() : num_owned_cols_; // number of vector of vectors
}


inline
const Entity_ID_List&
Mesh::cells_of_column(const int columnID_) const
{
  if (!columns_built_) build_columns_();
  return column_cells_[columnID_];
}


inline
const Entity_ID_List&
Mesh::faces_of_column(const int columnID_) const
{
  if (!columns_built_) build_columns_();
  return column_faces_[columnID_];
}


inline
int
Mesh::column_ID(const Entity_ID cellid) const
{
  if (!columns_built_) build_columns_();
  return columnID_[cellid];
}


inline
Entity_ID
Mesh::cell_get_cell_above(const Entity_ID cellid) const
{
  if (!columns_built_) build_columns_();
  return cell_cellabove_[cellid];
}


inline
Entity_ID
Mesh::cell_get_cell_below(const Entity_ID cellid) const
{
  if (!columns_built_) build_columns_();
  return cell_cellbelow_[cellid];
}


inline
Entity_ID
Mesh::node_get_node_above(const Entity_ID nodeid) const
{
  if (!columns_built_) build_columns_();
  return node_nodeabove_[nodeid];
}


//
// Column information
//-----------------------
inline
void
Mesh::get_set_entities(const std::string setname,
                       const Entity_kind kind,
                       const Parallel_type ptype,
                       Entity_ID_List *entids) const
{
  std::vector<double> vofs;
  get_set_entities_and_vofs(setname, kind, ptype, entids, &vofs);
}



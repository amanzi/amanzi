/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Rao Garimella, others
*/

//! <MISSING_ONELINE_DOCSTRING>

#include <algorithm>

#include "Teuchos_RCP.hpp"

#include "Mesh_simple.hh"
#include "GenerationSpec.hh"
#include "dbc.hh"
#include "errors.hh"


namespace Amanzi {
namespace AmanziMesh {

Mesh_simple::Mesh_simple(
  double x0, double y0, double z0, double x1, double y1, double z1, int nx,
  int ny, int nz, const Comm_ptr_type& comm,
  const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
  const Teuchos::RCP<const Teuchos::ParameterList>& plist,
  const bool request_faces, const bool request_edges)
  : Mesh(comm, gm, Teuchos::null, request_faces, request_edges),
    nx_(nx),
    ny_(ny),
    nz_(nz),
    x0_(x0),
    x1_(x1),
    y0_(y0),
    y1_(y1),
    z0_(z0),
    z1_(z1)
{
  // extract control parameters
  if (plist != Teuchos::null) {
    vo_ = Teuchos::rcp(new Amanzi::VerboseObject("Mesh:Simple", *plist));
  }

  Mesh::set_mesh_type(RECTANGULAR);
  Mesh::set_space_dimension(3);
  Mesh::set_manifold_dimension(3);
  if (gm != Teuchos::null) Mesh::set_geometric_model(gm);
  update();
  init_cache();
}


void
Mesh_simple::update()
{
  clear_internals_();
  update_internals_();
}


void
Mesh_simple::clear_internals_()
{
  coordinates_.resize(0);

  cell_to_face_.resize(0);
  cell_to_node_.resize(0);
  face_to_node_.resize(0);

  side_sets_.resize(0);
}


void
Mesh_simple::update_internals_()
{
  num_cells_ = nx_ * ny_ * nz_;
  num_nodes_ = (nx_ + 1) * (ny_ + 1) * (nz_ + 1);
  num_faces_ =
    (nx_ + 1) * ny_ * nz_ + nx_ * (ny_ + 1) * nz_ + nx_ * ny_ * (nz_ + 1);
  num_faces_bnd_ = 2 * nx_ * ny_ + 2 * nx_ * nz_ + 2 * ny_ * nz_;

  coordinates_.resize(3 * num_nodes_);

  double hx = (x1_ - x0_) / nx_;
  double hy = (y1_ - y0_) / ny_;
  double hz = (z1_ - z0_) / nz_;

  for (int iz = 0; iz <= nz_; iz++) {
    for (int iy = 0; iy <= ny_; iy++) {
      for (int ix = 0; ix <= nx_; ix++) {
        int istart = 3 * node_index_(ix, iy, iz);

        coordinates_[istart] = x0_ + ix * hx;
        coordinates_[istart + 1] = y0_ + iy * hy;
        coordinates_[istart + 2] = z0_ + iz * hz;
      }
    }
  }

  cell_to_face_.resize(6 * num_cells_);
  cell_to_face_dirs_.resize(6 * num_cells_);
  cell_to_node_.resize(8 * num_cells_);
  face_to_node_.resize(4 * num_faces_);
  node_to_face_.resize(13 * num_nodes_); // 1 extra for num faces
  node_to_cell_.resize(9 * num_nodes_);  // 1 extra for num cells
  face_to_cell_.resize(2 * num_faces_);
  face_to_cell_.assign(2 * num_faces_, -1);

  // loop over cells and initialize cell_to_node_
  for (int iz = 0; iz < nz_; iz++)
    for (int iy = 0; iy < ny_; iy++)
      for (int ix = 0; ix < nx_; ix++) {
        int jstart = 0;
        int istart = 8 * cell_index_(ix, iy, iz);
        int ncell = 0;

        cell_to_node_[istart] = node_index_(ix, iy, iz);
        cell_to_node_[istart + 1] = node_index_(ix + 1, iy, iz);
        cell_to_node_[istart + 2] = node_index_(ix + 1, iy + 1, iz);
        cell_to_node_[istart + 3] = node_index_(ix, iy + 1, iz);
        cell_to_node_[istart + 4] = node_index_(ix, iy, iz + 1);
        cell_to_node_[istart + 5] = node_index_(ix + 1, iy, iz + 1);
        cell_to_node_[istart + 6] = node_index_(ix + 1, iy + 1, iz + 1);
        cell_to_node_[istart + 7] = node_index_(ix, iy + 1, iz + 1);

        jstart = 9 * node_index_(ix, iy, iz);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart + 1 + ncell] = cell_index_(ix, iy, iz);
        (node_to_cell_[jstart])++;

        jstart = 9 * node_index_(ix + 1, iy, iz);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart + 1 + ncell] = cell_index_(ix, iy, iz);
        (node_to_cell_[jstart])++;

        jstart = 9 * node_index_(ix + 1, iy + 1, iz);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart + 1 + ncell] = cell_index_(ix, iy, iz);
        (node_to_cell_[jstart])++;

        jstart = 9 * node_index_(ix, iy + 1, iz);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart + 1 + ncell] = cell_index_(ix, iy, iz);
        (node_to_cell_[jstart])++;

        jstart = 9 * node_index_(ix, iy, iz + 1);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart + 1 + ncell] = cell_index_(ix, iy, iz);
        node_to_cell_[jstart]++;

        jstart = 9 * node_index_(ix + 1, iy, iz + 1); // 1 extra for num cells
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart + 1 + ncell] = cell_index_(ix, iy, iz);
        (node_to_cell_[jstart])++;

        jstart =
          9 * node_index_(ix + 1, iy + 1, iz + 1); // 1 extra for num cells
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart + 1 + ncell] = cell_index_(ix, iy, iz);
        (node_to_cell_[jstart])++;

        jstart = 9 * node_index_(ix, iy + 1, iz + 1); // 1 extra for num cells
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart + 1 + ncell] = cell_index_(ix, iy, iz);
        (node_to_cell_[jstart])++;
      }

  // loop over cells and initialize cell_to_face_
  for (int iz = 0; iz < nz_; iz++)
    for (int iy = 0; iy < ny_; iy++)
      for (int ix = 0; ix < nx_; ix++) {
        int istart = 6 * cell_index_(ix, iy, iz);
        int jstart = 0;

        cell_to_face_[istart] = xzface_index_(ix, iy, iz);
        cell_to_face_[istart + 1] = yzface_index_(ix + 1, iy, iz);
        cell_to_face_[istart + 2] = xzface_index_(ix, iy + 1, iz);
        cell_to_face_[istart + 3] = yzface_index_(ix, iy, iz);
        cell_to_face_[istart + 4] = xyface_index_(ix, iy, iz);
        cell_to_face_[istart + 5] = xyface_index_(ix, iy, iz + 1);

        cell_to_face_dirs_[istart] = 1;
        cell_to_face_dirs_[istart + 1] = 1;
        cell_to_face_dirs_[istart + 2] = -1;
        cell_to_face_dirs_[istart + 3] = -1;
        cell_to_face_dirs_[istart + 4] = -1;
        cell_to_face_dirs_[istart + 5] = 1;

        jstart = 2 * xzface_index_(ix, iy, iz);
        face_to_cell_[jstart + 1] = cell_index_(ix, iy, iz);

        jstart = 2 * yzface_index_(ix + 1, iy, iz);
        face_to_cell_[jstart + 1] = cell_index_(ix, iy, iz);

        jstart = 2 * xzface_index_(ix, iy + 1, iz);
        face_to_cell_[jstart] = cell_index_(ix, iy, iz);

        jstart = 2 * yzface_index_(ix, iy, iz);
        face_to_cell_[jstart] = cell_index_(ix, iy, iz);

        jstart = 2 * xyface_index_(ix, iy, iz);
        face_to_cell_[jstart] = cell_index_(ix, iy, iz);

        jstart = 2 * xyface_index_(ix, iy, iz + 1);
        face_to_cell_[jstart + 1] = cell_index_(ix, iy, iz);
      }

  // loop over faces and initialize face_to_node_
  // first we do the xy faces
  for (int iz = 0; iz <= nz_; iz++)
    for (int iy = 0; iy < ny_; iy++)
      for (int ix = 0; ix < nx_; ix++) {
        int istart = 4 * xyface_index_(ix, iy, iz);
        int jstart = 0;
        int nfaces = 0;

        face_to_node_[istart] = node_index_(ix, iy, iz);
        face_to_node_[istart + 1] = node_index_(ix + 1, iy, iz);
        face_to_node_[istart + 2] = node_index_(ix + 1, iy + 1, iz);
        face_to_node_[istart + 3] = node_index_(ix, iy + 1, iz);

        jstart = 13 * node_index_(ix, iy, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = 13 * node_index_(ix + 1, iy, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = 13 * node_index_(ix + 1, iy + 1, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = 13 * node_index_(ix, iy + 1, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;
      }

  // then we do the xz faces
  for (int iz = 0; iz < nz_; iz++)
    for (int iy = 0; iy <= ny_; iy++)
      for (int ix = 0; ix < nx_; ix++) {
        int istart = 4 * xzface_index_(ix, iy, iz);
        int jstart = 0;
        int nfaces = 0;

        face_to_node_[istart] = node_index_(ix, iy, iz);
        face_to_node_[istart + 1] = node_index_(ix + 1, iy, iz);
        face_to_node_[istart + 2] = node_index_(ix + 1, iy, iz + 1);
        face_to_node_[istart + 3] = node_index_(ix, iy, iz + 1);

        jstart = 13 * node_index_(ix, iy, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = 13 * node_index_(ix + 1, iy, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = 13 * node_index_(ix + 1, iy, iz + 1);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = 13 * node_index_(ix, iy, iz + 1);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;
      }

  // finally we do the yz faces
  for (int iz = 0; iz < nz_; iz++)
    for (int iy = 0; iy < ny_; iy++)
      for (int ix = 0; ix <= nx_; ix++) {
        int istart = 4 * yzface_index_(ix, iy, iz);
        int jstart = 0;
        int nfaces = 0;

        face_to_node_[istart] = node_index_(ix, iy, iz);
        face_to_node_[istart + 1] = node_index_(ix, iy + 1, iz);
        face_to_node_[istart + 2] = node_index_(ix, iy + 1, iz + 1);
        face_to_node_[istart + 3] = node_index_(ix, iy, iz + 1);

        jstart = 13 * node_index_(ix, iy, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = 13 * node_index_(ix, iy + 1, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = 13 * node_index_(ix, iy + 1, iz + 1);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = 13 * node_index_(ix, iy, iz + 1);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;
      }

  build_maps_();
}


void
Mesh_simple::build_maps_()
{
  std::vector<int> cells(num_cells_);
  for (int i = 0; i < num_cells_; i++) cells[i] = i;

  std::vector<int> nodes(num_nodes_);
  for (int i = 0; i < num_nodes_; i++) nodes[i] = i;

  std::vector<int> faces(num_faces_);
  for (int i = 0; i < num_faces_; i++) faces[i] = i;

  std::vector<int> faces_bnd(num_faces_bnd_);
  int i = 0;
  for (int ix = 0; ix < nx_; ix++)
    for (int iy = 0; iy < ny_; iy++) {
      faces_bnd[i++] = xyface_index_(ix, iy, 0);
      faces_bnd[i++] = xyface_index_(ix, iy, nz_);
    }
  for (int ix = 0; ix < nx_; ix++)
    for (int iz = 0; iz < nz_; iz++) {
      faces_bnd[i++] = xzface_index_(ix, 0, iz);
      faces_bnd[i++] = xzface_index_(ix, ny_, iz);
    }
  for (int iy = 0; iy < ny_; iy++)
    for (int iz = 0; iz < nz_; iz++) {
      faces_bnd[i++] = yzface_index_(0, iy, iz);
      faces_bnd[i++] = yzface_index_(nx_, iy, iz);
    }

  cell_map_ =
    Teuchos::rcp(new Map_type(-1, cells.data(), num_cells_, 0, get_comm()));
  face_map_ =
    Teuchos::rcp(new Map_type(-1, faces.data(), num_faces_, 0, get_comm()));
  node_map_ =
    Teuchos::rcp(new Map_type(-1, nodes.data(), num_nodes_, 0, get_comm()));
  extface_map_ = Teuchos::rcp(
    new Map_type(-1, faces_bnd.data(), num_faces_bnd_, 0, get_comm()));
}


Parallel_type
Mesh_simple::entity_get_ptype(const Entity_kind kind,
                              const Entity_ID entid) const
{
  return Parallel_type::OWNED; // Its a serial code
}


//--------------------------------------
// Get cell type
//--------------------------------------
AmanziMesh::Cell_type
Mesh_simple::cell_get_type(const AmanziMesh::Entity_ID cellid) const
{
  return HEX;
}


Entity_ID
Mesh_simple::GID(const AmanziMesh::Entity_ID lid,
                 const AmanziMesh::Entity_kind kind) const
{
  return lid; // Its a serial code
}


unsigned int
Mesh_simple::num_entities(AmanziMesh::Entity_kind kind,
                          AmanziMesh::Parallel_type ptype) const
{
  switch (kind) {
  case AmanziMesh::FACE:
    return (ptype != AmanziMesh::Parallel_type::GHOST) ? num_faces_ : 0;
    break;
  case AmanziMesh::NODE:
    return (ptype != AmanziMesh::Parallel_type::GHOST) ? num_nodes_ : 0;
    break;
  case AmanziMesh::CELL:
    return (ptype != AmanziMesh::Parallel_type::GHOST) ? num_cells_ : 0;
    break;
  default:
    throw std::exception();
    break;
  }
}


void
Mesh_simple::cell_get_faces_and_dirs_internal_(
  const AmanziMesh::Entity_ID cellid,
  Entity_ID_List& faceids,
  std::vector<int>& cfacedirs) const
{
  unsigned int offset = (unsigned int)6 * cellid;

  faceids.resize(6);
  cfacedirs.resize(6);

  for (int i = 0; i < 6; i++) {
    faceids[i] = *(cell_to_face_.begin() + offset);
    cfacedirs[i] = *(cell_to_face_dirs_.begin() + offset);
    offset++;
  }
}


void
Mesh_simple::cell_get_nodes(AmanziMesh::Entity_ID cell,
                            Entity_ID_List& nodeids) const
{
  unsigned int offset = (unsigned int)8 * cell;

  nodeids.resize(8); 

  for (int i = 0; i < 8; i++) {
    nodeids[i] = (*(cell_to_node_.begin() + offset));
    offset++;
  }
}


void
Mesh_simple::face_get_nodes(AmanziMesh::Entity_ID face,
                            Entity_ID_List& nodeids) const
{
  unsigned int offset = (unsigned int)4 * face;

  nodeids.resize(4); 

  for (int i = 0; i < 4; i++) {
    nodeids[i] = (*(face_to_node_.begin() + offset));
    offset++;
  }
}


//--------------------------------------
// Cooordinate Getters
//--------------------------------------
void
Mesh_simple::node_get_coordinates(const AmanziMesh::Entity_ID local_node_id,
                                  AmanziGeometry::Point* ncoords) const
{
  unsigned int offset = (unsigned int)3 * local_node_id;

  // ncoords->init(3);
  ncoords->set(3, &(coordinates_[offset]));
}


void
Mesh_simple::face_get_coordinates(
  AmanziMesh::Entity_ID local_face_id,
  std::vector<AmanziGeometry::Point>& fcoords) const
{
  Entity_ID_List node_indices(4);

  face_get_nodes(local_face_id, node_indices);

  fcoords.resize(node_indices.size());

  AmanziGeometry::Point xyz(3);
  for (int j = 0; j < node_indices.size(); ++j) {
    node_get_coordinates(node_indices[j], &xyz);
    fcoords[j] = xyz;
  }
}


void
Mesh_simple::cell_get_coordinates(
  AmanziMesh::Entity_ID local_cell_id,
  std::vector<AmanziGeometry::Point>& ccoords) const
{
  Entity_ID_List node_indices(8);

  cell_get_nodes(local_cell_id, node_indices);

  ccoords.resize(node_indices.size());

  AmanziGeometry::Point xyz(3);
  for (int j = 0; j < node_indices.size(); ++j) {
    node_get_coordinates(node_indices[j], &xyz);
    ccoords[j] = xyz;
  }
}


void
Mesh_simple::node_set_coordinates(const AmanziMesh::Entity_ID local_node_id,
                                  const double* ncoord)
{
  unsigned int offset = (unsigned int)3 * local_node_id;
  int spdim = Mesh::space_dimension();

  AMANZI_ASSERT(ncoord != NULL);

  std::vector<double>::iterator destination_begin =
    coordinates_.begin() + offset;
  for (int i = 0; i < spdim; i++) {
    *destination_begin = ncoord[i];
    destination_begin++;
  }
}


void
Mesh_simple::node_set_coordinates(const AmanziMesh::Entity_ID local_node_id,
                                  const AmanziGeometry::Point ncoord)
{
  unsigned int offset = (unsigned int)3 * local_node_id;

  std::vector<double>::iterator destination_begin =
    coordinates_.begin() + offset;
  // double coordarray[3] = {0.0,0.0,0.0}; // un-used
  int spdim = Mesh::space_dimension();
  // for (int i = 0; i < spdim; i++)
  //  coordarray[i] = ncoord[i];

  for (int i = 0; i < spdim; i++) {
    *destination_begin = ncoord[i];
    destination_begin++;
  }
}


void
Mesh_simple::node_get_cells(const AmanziMesh::Entity_ID nodeid,
                            const AmanziMesh::Parallel_type ptype,
                            Entity_ID_List& cellids) const
{
  unsigned int offset = (unsigned int)9 * nodeid;
  unsigned int ncells = node_to_cell_[offset];

  cellids.resize(ncells);

  for (int i = 0; i < ncells; i++) cellids[i] = (node_to_cell_[offset + i + 1]);
}


//--------------------------------------
// Faces of type 'ptype' connected to a node
//--------------------------------------
void
Mesh_simple::node_get_faces(const AmanziMesh::Entity_ID nodeid,
                            const AmanziMesh::Parallel_type ptype,
                            std::vector<AmanziMesh::Entity_ID>& faceids) const
{
  unsigned int offset = (unsigned int)13 * nodeid;
  unsigned int nfaces = node_to_face_[offset];

  faceids.resize(nfaces);

  for (int i = 0; i < nfaces; i++) faceids[i] = (node_to_face_[offset + i + 1]);
}


//--------------------------------------
// Get faces of ptype of a particular cell that are connected to the
// given node
//--------------------------------------
void
Mesh_simple::node_get_cell_faces(
  const AmanziMesh::Entity_ID nodeid, const AmanziMesh::Entity_ID cellid,
  const AmanziMesh::Parallel_type ptype,
  AmanziMesh::Entity_ID_List& faceids) const
{
  unsigned int offset = (unsigned int)6 * cellid;

  faceids.clear(); 

  for (int i = 0; i < 6; i++) {
    Entity_ID cellfaceid = cell_to_face_[offset + i];

    unsigned int offset2 = (unsigned int)4 * cellfaceid;

    Amanzi::AmanziMesh::Entity_ID_List fnodes;
    face_get_nodes(cellfaceid, fnodes);

    for (int j = 0; j < 4; j++) {
      if (face_to_node_[offset2 + j] == nodeid) {
        faceids.push_back(cellfaceid); 
        break;
      }
    }
  }
}


//--------------------------------------
// Cells connected to a face
//--------------------------------------
void
Mesh_simple::face_get_cells_internal_(
  const AmanziMesh::Entity_ID faceid, const AmanziMesh::Parallel_type ptype,
  AmanziMesh::Entity_ID_List& cellids) const
{
  unsigned int offset = (unsigned int)2 * faceid;
  cellids.resize(1); 

  if (face_to_cell_[offset] != -1) cellids[0] = face_to_cell_[offset];
  if (face_to_cell_[offset + 1] != -1) cellids[0] = face_to_cell_[offset + 1];
}


//-----------------------
// Same level adjacencies
//-----------------------

//--------------------------------------
// Face connected neighboring cells of given cell of a particular ptype
// (e.g. a hex has 6 face neighbors)

// The order in which the cellids are returned cannot be
// guaranteed in general except when ptype = ALL, in which case
// the cellids will correcpond to cells across the respective
// faces given by cell_get_faces
//--------------------------------------
void
Mesh_simple::cell_get_face_adj_cells(
  const AmanziMesh::Entity_ID cellid, const AmanziMesh::Parallel_type ptype,
  AmanziMesh::Entity_ID_List& fadj_cellids) const
{
  unsigned int offset = (unsigned int)6 * cellid;

  fadj_cellids.clear();
  for (int i = 0; i < 6; i++) {
    Entity_ID faceid = cell_to_face_[offset];

    unsigned int foffset = (unsigned int)2 * faceid;

    int adjcell0 = face_to_cell_[foffset];
    if (adjcell0 != -1 && adjcell0 != cellid) {
      fadj_cellids.push_back(adjcell0);
    } else {
      int adjcell1 = face_to_cell_[foffset + 1];
      if (adjcell1 != -1 && adjcell1 != cellid) {
        fadj_cellids.push_back(adjcell1);
      }
    }

    offset++;
  }
}


//--------------------------------------
// Node connected neighboring cells of given cell
// (a hex in a structured mesh has 26 node connected neighbors)
// The cells are returned in no particular order
//--------------------------------------
void
Mesh_simple::cell_get_node_adj_cells(
  const AmanziMesh::Entity_ID cellid, const AmanziMesh::Parallel_type ptype,
  AmanziMesh::Entity_ID_List& nadj_cellids) const
{
  unsigned int offset = (unsigned int)8 * cellid;

  nadj_cellids.resize(8); 

  for (int i = 0; i < 8; i++) {
    Entity_ID nodeid = cell_to_node_[offset + i];

    unsigned int offset2 = (unsigned int)9 * nodeid;
    unsigned int ncell = node_to_cell_[offset2];

    for (int j = 0; j < 8; j++) {
      Entity_ID nodecell = node_to_cell_[offset2 + j];

      unsigned int found = 0;
      unsigned int size = nadj_cellids.size();
      for (int k = 0; k < size; k++) {
        if (nadj_cellids[k] == nodecell) {
          found = 1;
          break;
        }
      }

      if (!found) {
        nadj_cellids.push_back(nodecell);
      }
    }
  }
}


Map_ptr_type
Mesh_simple::cell_map(bool include_ghost) const
{
  return cell_map_;
}


Map_ptr_type
Mesh_simple::face_map(bool include_ghost) const
{
  return face_map_;
}


Map_ptr_type
Mesh_simple::node_map(bool include_ghost) const
{
  return node_map_;
}


Map_ptr_type
Mesh_simple::exterior_face_map(bool include_ghost) const
{
  return extface_map_;
}


//--------------------------------------
// Epetra importer that will allow apps to import values from a Epetra
// vector defined on all owned faces into an Epetra vector defined
// only on exterior faces
//--------------------------------------
Import_ptr_type
Mesh_simple::exterior_face_importer(void) const
{
  Errors::Message mesg("exterior face importer is not implemented");
  amanzi_throw(mesg);
  return Teuchos::null;
}


void
Mesh_simple::get_set_entities_and_vofs(
  const std::string setname, const AmanziMesh::Entity_kind kind,
  const AmanziMesh::Parallel_type ptype,
  AmanziMesh::Entity_ID_List& setents,
  std::vector<double>* vofs) const
{
  // we ignore ptype since this is a serial implementation
  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm =
    Mesh::geometric_model();

  setents.clear(); 

  switch (kind) {
  case AmanziMesh::FACE: {
    Entity_ID_List ss;

    // Does this set exist?

    int nss = side_sets_.size(); // number of side sets
    bool found = false;
    for (int i = 0; i < nss; i++) {
      if ((side_set_regions_[i])->name() == setname) {
        found = true;
        ss = side_sets_[i];
      }
    }

    if (!found) { // create the side set from the region definition
      Teuchos::RCP<const AmanziGeometry::Region> rgn = gm->FindRegion(setname);

      if (rgn == Teuchos::null) {
        std::cerr << "Geometric model has no region named " << setname
                  << std::endl;
        std::cerr << "Cannot construct set by this name" << std::endl;
        throw std::exception();
      }

      if (rgn->manifold_dimension() != Mesh::manifold_dimension() - 1) {
        std::cerr << "Geometric model does not have a region named" << setname
                  << "with the appropriate dimension" << std::endl;
        std::cerr << "Cannot construct set by this name" << std::endl;
        throw std::exception();
      }

      if (rgn->type() == AmanziGeometry::BOX ||
          rgn->type() == AmanziGeometry::PLANE) {
        for (int ix = 0; ix < nx_; ix++)
          for (int iz = 0; iz < nz_; iz++) {
            int face;
            std::vector<AmanziGeometry::Point> fxyz;
            bool inbox;

            face = xzface_index_(ix, 0, iz);
            face_get_coordinates(face, fxyz);

            inbox = true;
            if (rgn->type() == AmanziGeometry::BOX) {
              AmanziGeometry::Point cenpnt(Mesh::space_dimension());
              for (int j = 0; j < 4; j++) cenpnt += fxyz[j];
              cenpnt /= 4.0;

              if (!rgn->inside(cenpnt)) inbox = false;
            } else {
              for (int j = 0; j < 4; j++) {
                if (!rgn->inside(fxyz[j])) {
                  inbox = false;
                  break;
                }
              }
            }

            if (inbox) ss.push_back(face);

            face = xzface_index_(ix, ny_, iz);
            face_get_coordinates(face, fxyz);

            inbox = true;
            if (rgn->type() == AmanziGeometry::BOX) {
              AmanziGeometry::Point cenpnt(Mesh::space_dimension());
              for (int j = 0; j < 4; j++) cenpnt += fxyz[j];
              cenpnt /= 4.0;

              if (!rgn->inside(cenpnt)) inbox = false;
            } else {
              for (int j = 0; j < 4; j++) {
                if (!rgn->inside(fxyz[j])) {
                  inbox = false;
                  break;
                }
              }
            }

            if (inbox) ss.push_back(face);
          }

        for (int ix = 0; ix < nx_; ix++)
          for (int iy = 0; iy < ny_; iy++) {
            int face;
            std::vector<AmanziGeometry::Point> fxyz;
            bool inbox;

            face = xyface_index_(ix, iy, 0);
            face_get_coordinates(face, fxyz);

            inbox = true;
            if (rgn->type() == AmanziGeometry::BOX) {
              AmanziGeometry::Point cenpnt(Mesh::space_dimension());
              for (int j = 0; j < 4; j++) cenpnt += fxyz[j];
              cenpnt /= 4.0;

              if (!rgn->inside(cenpnt)) inbox = false;
            } else {
              for (int j = 0; j < 4; j++) {
                if (!rgn->inside(fxyz[j])) {
                  inbox = false;
                  break;
                }
              }
            }

            if (inbox) ss.push_back(face);

            face = xyface_index_(ix, iy, nz_);
            face_get_coordinates(face, fxyz);

            inbox = true;
            if (rgn->type() == AmanziGeometry::BOX) {
              AmanziGeometry::Point cenpnt(Mesh::space_dimension());
              for (int j = 0; j < 4; j++) cenpnt += fxyz[j];
              cenpnt /= 4.0;

              if (!rgn->inside(cenpnt)) inbox = false;
            } else {
              for (int j = 0; j < 4; j++) {
                if (!rgn->inside(fxyz[j])) {
                  inbox = false;
                  break;
                }
              }
            }

            if (inbox) ss.push_back(face);
          }

        for (int iy = 0; iy < ny_; iy++)
          for (int iz = 0; iz < nz_; iz++) {
            int face;
            std::vector<AmanziGeometry::Point> fxyz;
            bool inbox;

            face = yzface_index_(0, iy, iz);
            face_get_coordinates(face, fxyz);

            inbox = true;
            if (rgn->type() == AmanziGeometry::BOX) {
              AmanziGeometry::Point cenpnt(Mesh::space_dimension());
              for (int j = 0; j < 4; j++) cenpnt += fxyz[j];
              cenpnt /= 4.0;

              if (!rgn->inside(cenpnt)) inbox = false;
            } else {
              for (int j = 0; j < 4; j++) {
                if (!rgn->inside(fxyz[j])) {
                  inbox = false;
                  break;
                }
              }
            }

            if (inbox) ss.push_back(face);

            face = yzface_index_(nx_, iy, iz);
            face_get_coordinates(face, fxyz);

            inbox = true;
            if (rgn->type() == AmanziGeometry::BOX) {
              AmanziGeometry::Point cenpnt(Mesh::space_dimension());
              for (int j = 0; j < 4; j++) cenpnt += fxyz[j];
              cenpnt /= 4.0;

              if (!rgn->inside(cenpnt)) inbox = false;
            } else {
              for (int j = 0; j < 4; j++) {
                if (!rgn->inside(fxyz[j])) {
                  inbox = false;
                  break;
                }
              }
            }

            if (inbox) ss.push_back(face);
          }

        side_sets_.push_back(ss);
        side_set_regions_.push_back(rgn);
      }

      else {
        std::cerr << "Region type not suitable/applicable for sidesets"
                  << std::endl;
        throw std::exception();
      }
    }

    setents = ss; 
    break;
  }

  case AmanziMesh::CELL: {
    Entity_ID_List cs; // cell set

    int ncs = element_blocks_.size();
    bool found = false;
    for (int i = 0; i < ncs; i++) {
      if ((element_block_regions_[i])->name() == setname) {
        found = true;
        cs = element_blocks_[i];
      }
    }

    if (!found) { // create the cell set from the region definition
      Teuchos::RCP<const AmanziGeometry::Region> rgn = gm->FindRegion(setname);

      if (rgn == Teuchos::null) {
        std::cerr << "Geometric model has no region named " << setname
                  << std::endl;
        std::cerr << "Cannot construct set by this name" << std::endl;
        throw std::exception();
      }

      if (rgn->manifold_dimension() > 0 &&
          rgn->manifold_dimension() != Mesh::manifold_dimension()) {
        std::cerr << "Geometric model does not have a region named" << setname
                  << "with the appropriate dimension" << std::endl;
        std::cerr << "Cannot construct set by this name" << std::endl;
        throw std::exception();
      }

      if (rgn->type() == AmanziGeometry::BOX ||
          rgn->type() == AmanziGeometry::COLORFUNCTION) {
        for (int ix = 0; ix < nx_; ix++)
          for (int iy = 0; iy < ny_; iy++)
            for (int iz = 0; iz < nz_; iz++) {
              int cell = cell_index_(ix, iy, iz);
              std::vector<AmanziGeometry::Point> cxyz;
              cell_get_coordinates(cell, cxyz);

              AmanziGeometry::Point cenpnt(Mesh::space_dimension());

              for (int j = 0; j < 8; j++) cenpnt += cxyz[j];
              cenpnt /= 8.0;

              if (rgn->inside(cenpnt)) cs.push_back(cell);
            }

        element_blocks_.push_back(cs);
        element_block_regions_.push_back(rgn);
      } else {
        std::cerr << "Region type not suitable/applicable for cellsets"
                  << std::endl;
        throw std::exception();
      }
    }

    setents = cs; 
    break;
  }

  case AmanziMesh::NODE: {
    Entity_ID_List ns;

    // Does this set exist?

    int nns = node_sets_.size(); // number of node sets
    bool found = false;
    for (int i = 0; i < nns; i++) {
      if ((node_set_regions_[i])->name() == setname) {
        found = true;
        ns = node_sets_[i];
      }
    }

    if (!found) { // create the side set from the region definition
      Teuchos::RCP<const AmanziGeometry::Region> rgn = gm->FindRegion(setname);

      if (rgn == Teuchos::null) {
        std::cerr << "Geometric model has no region named " << setname
                  << std::endl;
        std::cerr << "Cannot construct set by this name" << std::endl;
        throw std::exception();
      }

      bool done = false;
      for (int ix = 0; ix < nx_ + 1 && !done; ix++)
        for (int iy = 0; iy < ny_ + 1 && !done; iy++)
          for (int iz = 0; iz < nz_ + 1 && !done; iz++) {
            int node = node_index_(ix, iy, iz);
            AmanziGeometry::Point nxyz;
            node_get_coordinates(node, &nxyz);

            if (rgn->inside(nxyz)) {
              ns.push_back(node);

              if (rgn->type() == AmanziGeometry::POINT) done = true;
            }
          }

      node_sets_.push_back(ns);
      node_set_regions_.push_back(rgn);
    }

    setents = ns; 
    break;
  }

  default:
    // set type not recognized
    throw std::exception();
  }
}


//--------------------------------------
// Deform a mesh so that cell volumes conform as closely as possible
// to target volumes without dropping below the minimum volumes.  If
// move_vertical = true, nodes will be allowed to move only in the
// vertical direction (right now arbitrary node movement is not allowed)
//--------------------------------------
int
Mesh_simple::deform(const std::vector<double>& target_cell_volumes_in,
                    const std::vector<double>& min_cell_volumes_in,
                    const Entity_ID_List& fixed_nodes,
                    const bool move_vertical)
{
  Errors::Message mesg("Deformation not implemented for Mesh_simple");
  amanzi_throw(mesg);
  return 0;
}


//--------------------------------------
// Write mesh out to exodus file
//--------------------------------------
void
Mesh_simple::write_to_exodus_file(const std::string filename) const
{
  Errors::Message mesg("Not implemented");
  amanzi_throw(mesg);
}

} // namespace AmanziMesh
} // namespace Amanzi

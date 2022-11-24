/*
  Mesh

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, others
*/

#include <algorithm>

#include <Teuchos_RCP.hpp>

#include "dbc.hh"
#include "errors.hh"
#include "GenerationSpec.hh"
#include "RegionLogical.hh"

#include "Mesh_simple.hh"


namespace Amanzi {
namespace AmanziMesh {

//---------------------------------------------------------
// Constructor
//---------------------------------------------------------
Mesh_simple::Mesh_simple(double x0,
                         double y0,
                         double z0,
                         double x1,
                         double y1,
                         double z1,
                         int nx,
                         int ny,
                         int nz,
                         const Comm_ptr_type& comm,
                         const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
                         const Teuchos::RCP<const Teuchos::ParameterList>& plist,
                         const bool request_faces,
                         const bool request_edges)
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

  set_mesh_type(RECTANGULAR);
  set_space_dimension(3);
  set_manifold_dimension(3);
  if (gm != Teuchos::null) set_geometric_model(gm);

  CreateCache_();
  BuildMaps_();
}


//---------------------------------------------------------
// Update
//---------------------------------------------------------
void
Mesh_simple::CreateCache_()
{
  // clear old cache
  coordinates_.clear();

  cell_to_face_.clear();
  cell_to_edge_.clear();
  cell_to_node_.clear();
  face_to_node_.clear();
  edge_to_node_.clear();

  sets_.clear();

  // build new cache
  num_cells_ = nx_ * ny_ * nz_;
  num_nodes_ = (nx_ + 1) * (ny_ + 1) * (nz_ + 1);
  num_faces_ = (nx_ + 1) * ny_ * nz_ + nx_ * (ny_ + 1) * nz_ + nx_ * ny_ * (nz_ + 1);
  num_edges_ =
    nx_ * (ny_ + 1) * (nz_ + 1) + (nx_ + 1) * ny_ * (nz_ + 1) + (nx_ + 1) * (ny_ + 1) * nz_;
  num_faces_bnd_ = 2 * nx_ * ny_ + 2 * nx_ * nz_ + 2 * ny_ * nz_;

  // -- node coordinates
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

  // -- connectivity arrays
  cell_to_face_.resize(6 * num_cells_);
  cell_to_face_dirs_.resize(6 * num_cells_);
  face_to_cell_.assign(2 * num_faces_, -1);

  cell_to_node_.resize(8 * num_cells_);
  node_to_cell_.resize(9 * num_nodes_); // 1 extra for num cells

  face_to_node_.resize(4 * num_faces_);
  node_to_face_.resize(13 * num_nodes_); // 1 extra for num faces

  if (edges_requested_) {
    cell_to_edge_.resize(12 * num_edges_);
    face_to_edge_.resize(4 * num_faces_);
    face_to_edge_dirs_.resize(4 * num_faces_);

    edge_to_node_.resize(2 * num_edges_);
  }

  // loop over cells and initialize cell <-> node
  for (int iz = 0; iz < nz_; iz++) {
    for (int iy = 0; iy < ny_; iy++) {
      for (int ix = 0; ix < nx_; ix++) {
        for (int k = 0; k < 2; ++k) {
          int istart = 8 * cell_index_(ix, iy, iz) + 4 * k;

          cell_to_node_[istart] = node_index_(ix, iy, iz + k);
          cell_to_node_[istart + 1] = node_index_(ix + 1, iy, iz + k);
          cell_to_node_[istart + 2] = node_index_(ix + 1, iy + 1, iz + k);
          cell_to_node_[istart + 3] = node_index_(ix, iy + 1, iz + k);
        }

        for (int i = 0; i < 2; ++i) {
          for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
              int jstart = 9 * node_index_(ix + i, iy + j, iz + k);
              int ncell = node_to_cell_[jstart];
              node_to_cell_[jstart + 1 + ncell] = cell_index_(ix, iy, iz);
              (node_to_cell_[jstart])++;
            }
          }
        }
      }
    }
  }

  // loop over cells and initialize cell <-> face
  for (int iz = 0; iz < nz_; iz++) {
    for (int iy = 0; iy < ny_; iy++) {
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
    }
  }

  // loop over faces and initialize face <-> node
  // -- xy faces
  for (int iz = 0; iz <= nz_; iz++) {
    for (int iy = 0; iy < ny_; iy++) {
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
    }
  }

  // -- xz faces
  for (int iz = 0; iz < nz_; iz++) {
    for (int iy = 0; iy <= ny_; iy++) {
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
    }
  }

  // -- yz faces
  for (int iz = 0; iz < nz_; iz++) {
    for (int iy = 0; iy < ny_; iy++) {
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
    }
  }

  if (edges_requested_) {
    // loop over cells and initialize cell -> edges
    for (int iz = 0; iz < nz_; iz++) {
      for (int iy = 0; iy < ny_; iy++) {
        for (int ix = 0; ix < nx_; ix++) {
          int istart = 6 * cell_index_(ix, iy, iz);

          cell_to_edge_[istart] = xedge_index_(ix, iy, iz);
          cell_to_edge_[istart + 1] = xedge_index_(ix, iy + 1, iz);
          cell_to_edge_[istart + 2] = xedge_index_(ix, iy + 1, iz + 1);
          cell_to_edge_[istart + 3] = xedge_index_(ix, iy, iz + 1);

          cell_to_edge_[istart + 4] = yedge_index_(ix, iy, iz);
          cell_to_edge_[istart + 5] = yedge_index_(ix + 1, iy, iz);
          cell_to_edge_[istart + 6] = yedge_index_(ix + 1, iy, iz + 1);
          cell_to_edge_[istart + 7] = yedge_index_(ix, iy, iz + 1);

          cell_to_edge_[istart + 8] = zedge_index_(ix, iy, iz);
          cell_to_edge_[istart + 9] = zedge_index_(ix + 1, iy, iz);
          cell_to_edge_[istart + 10] = zedge_index_(ix + 1, iy + 1, iz);
          cell_to_edge_[istart + 11] = zedge_index_(ix, iy + 1, iz);
        }
      }
    }


    // loop over faces and initialize face -> edge
    // -- xy faces
    for (int iz = 0; iz <= nz_; iz++) {
      for (int iy = 0; iy < ny_; iy++) {
        for (int ix = 0; ix < nx_; ix++) {
          int istart = 4 * xyface_index_(ix, iy, iz);

          face_to_edge_[istart] = xedge_index_(ix, iy, iz);
          face_to_edge_[istart + 1] = yedge_index_(ix + 1, iy, iz);
          face_to_edge_[istart + 2] = xedge_index_(ix, iy + 1, iz);
          face_to_edge_[istart + 3] = yedge_index_(ix, iy, iz);

          face_to_edge_dirs_[istart] = 1;
          face_to_edge_dirs_[istart + 1] = 1;
          face_to_edge_dirs_[istart + 2] = -1;
          face_to_edge_dirs_[istart + 3] = -1;
        }
      }
    }

    // -- xz faces
    for (int iz = 0; iz < nz_; iz++) {
      for (int iy = 0; iy <= ny_; iy++) {
        for (int ix = 0; ix < nx_; ix++) {
          int istart = 4 * xzface_index_(ix, iy, iz);

          face_to_edge_[istart] = xedge_index_(ix, iy, iz);
          face_to_edge_[istart + 1] = zedge_index_(ix + 1, iy, iz);
          face_to_edge_[istart + 2] = xedge_index_(ix, iy, iz + 1);
          face_to_edge_[istart + 3] = zedge_index_(ix, iy, iz);

          face_to_edge_dirs_[istart] = 1;
          face_to_edge_dirs_[istart + 1] = 1;
          face_to_edge_dirs_[istart + 2] = -1;
          face_to_edge_dirs_[istart + 3] = -1;
        }
      }
    }

    // -- yz faces
    for (int iz = 0; iz < nz_; iz++) {
      for (int iy = 0; iy < ny_; iy++) {
        for (int ix = 0; ix <= nx_; ix++) {
          int istart = 4 * yzface_index_(ix, iy, iz);

          face_to_edge_[istart] = yedge_index_(ix, iy, iz);
          face_to_edge_[istart + 1] = zedge_index_(ix, iy + 1, iz);
          face_to_edge_[istart + 2] = yedge_index_(ix, iy, iz + 1);
          face_to_edge_[istart + 3] = zedge_index_(ix, iy, iz);

          face_to_edge_dirs_[istart] = 1;
          face_to_edge_dirs_[istart + 1] = 1;
          face_to_edge_dirs_[istart + 2] = -1;
          face_to_edge_dirs_[istart + 3] = -1;
        }
      }
    }

    // loop over edges and initialize edge -> nodes
    // -- x edges
    for (int iz = 0; iz <= nz_; iz++) {
      for (int iy = 0; iy <= ny_; iy++) {
        for (int ix = 0; ix < nx_; ix++) {
          int istart = 2 * xedge_index_(ix, iy, iz);

          edge_to_node_[istart] = node_index_(ix, iy, iz);
          edge_to_node_[istart + 1] = node_index_(ix + 1, iy, iz);
        }
      }
    }

    // -- y edges
    for (int iz = 0; iz <= nz_; iz++) {
      for (int iy = 0; iy < ny_; iy++) {
        for (int ix = 0; ix <= nx_; ix++) {
          int istart = 2 * yedge_index_(ix, iy, iz);

          edge_to_node_[istart] = node_index_(ix, iy, iz);
          edge_to_node_[istart + 1] = node_index_(ix, iy + 1, iz);
        }
      }
    }

    // -- z edges
    for (int iz = 0; iz < nz_; iz++) {
      for (int iy = 0; iy <= ny_; iy++) {
        for (int ix = 0; ix <= nx_; ix++) {
          int istart = 2 * zedge_index_(ix, iy, iz);

          edge_to_node_[istart] = node_index_(ix, iy, iz);
          edge_to_node_[istart + 1] = node_index_(ix, iy, iz + 1);
        }
      }
    }
  }
}


//---------------------------------------------------------
// Build Epetra maps
//---------------------------------------------------------
void
Mesh_simple::BuildMaps_()
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

  cell_map_ = Teuchos::rcp(new Epetra_Map(-1, num_cells_, &cells[0], 0, *get_comm()));
  face_map_ = Teuchos::rcp(new Epetra_Map(-1, num_faces_, &faces[0], 0, *get_comm()));
  node_map_ = Teuchos::rcp(new Epetra_Map(-1, num_nodes_, &nodes[0], 0, *get_comm()));
  extface_map_ = Teuchos::rcp(new Epetra_Map(-1, num_faces_bnd_, &faces_bnd[0], 0, *get_comm()));

  if (edges_requested_) {
    std::vector<int> edges(num_edges_);
    for (int n = 0; n < num_edges_; ++n) edges[n] = n;

    edge_map_ = Teuchos::rcp(new Epetra_Map(-1, num_edges_, &edges[0], 0, *get_comm()));
  }
}


//---------------------------------------------------------
// TBW
//---------------------------------------------------------
unsigned int
Mesh_simple::num_entities(AmanziMesh::Entity_kind kind, AmanziMesh::Parallel_type ptype) const
{
  switch (kind) {
  case FACE:
    return (ptype != AmanziMesh::Parallel_type::GHOST) ? num_faces_ : 0;
    break;
  case NODE:
    return (ptype != AmanziMesh::Parallel_type::GHOST) ? num_nodes_ : 0;
    break;
  case CELL:
    return (ptype != AmanziMesh::Parallel_type::GHOST) ? num_cells_ : 0;
    break;
  default:
    throw std::exception();
    break;
  }
}


//---------------------------------------------------------
// Connectivity: cell -> faces
//---------------------------------------------------------
void
Mesh_simple::cell_get_faces_and_dirs_internal_(const AmanziMesh::Entity_ID cellid,
                                               AmanziMesh::Entity_ID_List* faceids,
                                               std::vector<int>* cfacedirs,
                                               const bool ordered) const
{
  unsigned int offset = (unsigned int)6 * cellid;

  faceids->clear();
  auto it = cell_to_face_.begin() + offset;
  faceids->insert(faceids->begin(), it, it + 6);

  if (cfacedirs) {
    cfacedirs->clear();
    auto jt = cell_to_face_dirs_.begin() + offset;
    cfacedirs->insert(cfacedirs->begin(), jt, jt + 6);
  }
}


//---------------------------------------------------------
// Connectivity: cell -> edges
//---------------------------------------------------------
void
Mesh_simple::cell_get_edges_internal_(const Entity_ID cellid, Entity_ID_List* edgeids) const
{
  unsigned int offset = (unsigned int)12 * cellid;

  edgeids->clear();
  auto it = cell_to_edge_.begin() + offset;
  edgeids->insert(edgeids->begin(), it, it + 12);
}


//---------------------------------------------------------
// Connectivity: cell -> nodes
//---------------------------------------------------------
void
Mesh_simple::cell_get_nodes(AmanziMesh::Entity_ID cell, AmanziMesh::Entity_ID_List* nodeids) const
{
  unsigned int offset = (unsigned int)8 * cell;

  nodeids->clear();

  for (int i = 0; i < 8; i++) {
    nodeids->push_back(*(cell_to_node_.begin() + offset));
    offset++;
  }
}


//---------------------------------------------------------
// Connectivity: face -> nodes
//---------------------------------------------------------
void
Mesh_simple::face_get_nodes(AmanziMesh::Entity_ID face, AmanziMesh::Entity_ID_List* nodeids) const
{
  unsigned int offset = (unsigned int)4 * face;

  nodeids->clear();

  for (int i = 0; i < 4; i++) {
    nodeids->push_back(*(face_to_node_.begin() + offset));
    offset++;
  }
}


//---------------------------------------------------------
// Connectivity: face -> edges
//---------------------------------------------------------
void
Mesh_simple::face_get_edges_and_dirs_internal_(const Entity_ID faceid,
                                               Entity_ID_List* edgeids,
                                               std::vector<int>* fedgedirs,
                                               bool ordered) const
{
  unsigned int offset = (unsigned int)4 * faceid;

  edgeids->clear();
  auto it = face_to_edge_.begin() + offset;
  edgeids->insert(edgeids->begin(), it, it + 4);

  if (fedgedirs) {
    fedgedirs->clear();
    auto jt = face_to_edge_dirs_.begin();
    fedgedirs->insert(fedgedirs->begin(), jt, jt + 4);
  }
}


//---------------------------------------------------------
// Connectivity: edge -> nodes
//---------------------------------------------------------
void
Mesh_simple::edge_get_nodes(const Entity_ID edgeid, Entity_ID* nodeid0, Entity_ID* nodeid1) const
{
  unsigned int offset = (unsigned int)2 * edgeid;
  *nodeid0 = edge_to_node_[offset];
  *nodeid1 = edge_to_node_[offset + 1];
}


//---------------------------------------------------------
// Cooordinates of a node
//---------------------------------------------------------
void
Mesh_simple::node_get_coordinates(const AmanziMesh::Entity_ID local_node_id,
                                  AmanziGeometry::Point* ncoords) const
{
  unsigned int offset = (unsigned int)3 * local_node_id;

  // ncoords->init(3);
  ncoords->set(3, &(coordinates_[offset]));
}


//---------------------------------------------------------
// Modify cooordinates of a node: version 1
//---------------------------------------------------------
void
Mesh_simple::node_set_coordinates(const AmanziMesh::Entity_ID local_node_id, const double* ncoord)
{
  unsigned int offset = (unsigned int)3 * local_node_id;
  int spdim = Mesh::space_dimension();

  AMANZI_ASSERT(ncoord != NULL);

  std::vector<double>::iterator destination_begin = coordinates_.begin() + offset;
  for (int i = 0; i < spdim; i++) {
    *destination_begin = ncoord[i];
    destination_begin++;
  }
}


//---------------------------------------------------------
// Modify cooordinates of a node: version 2
//---------------------------------------------------------
void
Mesh_simple::node_set_coordinates(const AmanziMesh::Entity_ID local_node_id,
                                  const AmanziGeometry::Point ncoord)
{
  unsigned int offset = (unsigned int)3 * local_node_id;

  std::vector<double>::iterator destination_begin = coordinates_.begin() + offset;
  int spdim = Mesh::space_dimension();
  for (int i = 0; i < spdim; i++) {
    *destination_begin = ncoord[i];
    destination_begin++;
  }
}


//-------------------
// Upward adjacencies
//-------------------

//---------------------------------------------------------
// Connectivity: node -> cells
//---------------------------------------------------------
void
Mesh_simple::node_get_cells(const AmanziMesh::Entity_ID nodeid,
                            const AmanziMesh::Parallel_type ptype,
                            AmanziMesh::Entity_ID_List* cellids) const
{
  unsigned int offset = (unsigned int)9 * nodeid;
  unsigned int ncells = node_to_cell_[offset];

  cellids->clear();

  for (int i = 0; i < ncells; i++) cellids->push_back(node_to_cell_[offset + i + 1]);
}


//---------------------------------------------------------
// Faces of type 'ptype' connected to a node
//---------------------------------------------------------
void
Mesh_simple::node_get_faces(const AmanziMesh::Entity_ID nodeid,
                            const AmanziMesh::Parallel_type ptype,
                            AmanziMesh::Entity_ID_List* faceids) const
{
  unsigned int offset = (unsigned int)13 * nodeid;
  unsigned int nfaces = node_to_face_[offset];

  faceids->clear();

  for (int i = 0; i < nfaces; i++) faceids->push_back(node_to_face_[offset + i + 1]);
}


//---------------------------------------------------------
// Cells connected to a face
//---------------------------------------------------------
void
Mesh_simple::face_get_cells_internal_(const AmanziMesh::Entity_ID faceid,
                                      const AmanziMesh::Parallel_type ptype,
                                      AmanziMesh::Entity_ID_List* cellids) const
{
  unsigned int offset = (unsigned int)2 * faceid;

  cellids->clear();

  if (face_to_cell_[offset] != -1) cellids->push_back(face_to_cell_[offset]);
  if (face_to_cell_[offset + 1] != -1) cellids->push_back(face_to_cell_[offset + 1]);
}


//-----------------------
// Same level adjacencies
//-----------------------

//---------------------------------------------------------
// Face connected neighboring cells of given cell of a particular ptype
// (e.g. a hex has 6 face neighbors)

// The order in which the cellids are returned cannot be
// guaranteed in general except when ptype = ALL, in which case
// the cellids will correcpond to cells across the respective
// faces given by cell_get_faces
//---------------------------------------------------------
void
Mesh_simple::cell_get_face_adj_cells(const AmanziMesh::Entity_ID cellid,
                                     const AmanziMesh::Parallel_type ptype,
                                     AmanziMesh::Entity_ID_List* fadj_cellids) const
{
  unsigned int offset = (unsigned int)6 * cellid;

  fadj_cellids->clear();

  for (int i = 0; i < 6; i++) {
    Entity_ID faceid = cell_to_face_[offset];

    unsigned int foffset = (unsigned int)2 * faceid;

    int adjcell0 = face_to_cell_[foffset];
    if (adjcell0 != -1 && adjcell0 != cellid)
      fadj_cellids->push_back(adjcell0);
    else {
      int adjcell1 = face_to_cell_[foffset + 1];
      if (adjcell1 != -1 && adjcell1 != cellid) fadj_cellids->push_back(adjcell1);
    }

    offset++;
  }
}


//---------------------------------------------------------
// TBW
//---------------------------------------------------------
const Epetra_Map&
Mesh_simple::exterior_node_map(bool include_ghost) const
{
  Errors::Message mesg("Exterior node map is not implemented in this framework");
  Exceptions::amanzi_throw(mesg);
  throw(mesg);
}


//---------------------------------------------------------
// Epetra importer that will allow apps to import values from a Epetra
// vector defined on all owned faces into an Epetra vector defined
// only on exterior faces
//---------------------------------------------------------
const Epetra_Import&
Mesh_simple::exterior_face_importer(void) const
{
  Errors::Message mesg("exterior face importer is not implemented");
  amanzi_throw(mesg);
  throw(mesg); // this silences compiler warnings but is never called
}


//---------------------------------------------------------
// TBW
//---------------------------------------------------------
void
Mesh_simple::get_set_entities_and_vofs(const std::string& setname,
                                       const AmanziMesh::Entity_kind kind,
                                       const AmanziMesh::Parallel_type ptype,
                                       AmanziMesh::Entity_ID_List* setents,
                                       std::vector<double>* vofs) const
{
  // we ignore ptype since this is a serial implementation
  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm = Mesh::geometric_model();

  AMANZI_ASSERT(setents != NULL);

  std::string setname_internal = setname + std::to_string(kind);
  {
    auto it = sets_.find(setname_internal);
    if (it != sets_.end()) {
      *setents = it->second;
      return;
    }
  }

  // create the side set from the region definition
  Teuchos::RCP<const AmanziGeometry::Region> rgn = gm->FindRegion(setname);
  if (rgn == Teuchos::null) {
    std::cerr << "Geometric model has no region \"" << setname << "\"" << std::endl;
    throw std::exception();
  }

  setents->clear();

  switch (kind) {
  case FACE: {
    if (rgn->get_type() == AmanziGeometry::RegionType::BOX ||
        rgn->get_type() == AmanziGeometry::RegionType::PLANE) {
      // y = constant planes
      for (int iy = 0; iy <= ny_; ++iy) {
        for (int ix = 0; ix < nx_; ix++) {
          for (int iz = 0; iz < nz_; iz++) {
            std::vector<AmanziGeometry::Point> fxyz;

            int face = xzface_index_(ix, iy, iz);
            face_get_coordinates(face, &fxyz);

            bool inbox = true;
            if (rgn->get_type() == AmanziGeometry::RegionType::BOX) {
              auto cenpnt = geometric_center_(fxyz);
              if (!rgn->inside(cenpnt)) inbox = false;
            } else {
              inbox = all_inside_(fxyz, *rgn);
            }

            if (inbox) setents->push_back(face);
          }
        }
      }

      // z = constant planes
      for (int iz = 0; iz <= nz_; ++iz) {
        for (int ix = 0; ix < nx_; ix++) {
          for (int iy = 0; iy < ny_; iy++) {
            std::vector<AmanziGeometry::Point> fxyz;

            int face = xyface_index_(ix, iy, iz);
            face_get_coordinates(face, &fxyz);

            bool inbox = true;
            if (rgn->get_type() == AmanziGeometry::RegionType::BOX) {
              auto cenpnt = geometric_center_(fxyz);
              if (!rgn->inside(cenpnt)) inbox = false;
            } else {
              inbox = all_inside_(fxyz, *rgn);
            }

            if (inbox) setents->push_back(face);
          }
        }
      }

      // x = constant planes
      for (int ix = 0; ix <= nx_; ++ix) {
        for (int iy = 0; iy < ny_; iy++) {
          for (int iz = 0; iz < nz_; iz++) {
            std::vector<AmanziGeometry::Point> fxyz;

            int face = yzface_index_(ix, iy, iz);
            face_get_coordinates(face, &fxyz);

            bool inbox = true;
            if (rgn->get_type() == AmanziGeometry::RegionType::BOX) {
              auto cenpnt = geometric_center_(fxyz);
              if (!rgn->inside(cenpnt)) inbox = false;
            } else {
              inbox = all_inside_(fxyz, *rgn);
            }

            if (inbox) setents->push_back(face);
          }
        }
      }

      sets_[setname_internal] = *setents;
    }

    else if (rgn->get_type() != AmanziGeometry::RegionType::LOGICAL) {
      std::cerr << "Region \"" << rgn->get_name()
                << "\" type not applicable/supported for sidesets";
      throw std::exception();
    }

    break;
  }

  case CELL: {
    if (rgn->get_type() == AmanziGeometry::RegionType::BOX ||
        rgn->get_type() == AmanziGeometry::RegionType::COLORFUNCTION) {
      for (int ix = 0; ix < nx_; ix++) {
        for (int iy = 0; iy < ny_; iy++) {
          for (int iz = 0; iz < nz_; iz++) {
            int cell = cell_index_(ix, iy, iz);
            std::vector<AmanziGeometry::Point> cxyz;
            cell_get_coordinates(cell, &cxyz);

            auto cenpnt = geometric_center_(cxyz);
            if (rgn->inside(cenpnt)) setents->push_back(cell);
          }
        }
      }

      sets_[setname_internal] = *setents;
    } else {
      std::cerr << "Region \"" << rgn->get_name()
                << "\" type not applicable/supported for cellsets";
      throw std::exception();
    }

    break;
  }

  case NODE: {
    if (rgn->get_type() != AmanziGeometry::RegionType::LOGICAL) {
      bool done = false;
      for (int ix = 0; ix < nx_ + 1 && !done; ix++) {
        for (int iy = 0; iy < ny_ + 1 && !done; iy++) {
          for (int iz = 0; iz < nz_ + 1 && !done; iz++) {
            int node = node_index_(ix, iy, iz);
            AmanziGeometry::Point xyz;
            node_get_coordinates(node, &xyz);

            if (rgn->inside(xyz)) {
              setents->push_back(node);
              if (rgn->get_type() == AmanziGeometry::RegionType::POINT) done = true;
            }
          }
        }
      }

      sets_[setname_internal] = *setents;
    }
  } break;

  default:
    break;
  }

  if (rgn->get_type() == AmanziGeometry::RegionType::LOGICAL) {
    auto boolrgn = Teuchos::rcp_static_cast<const AmanziGeometry::RegionLogical>(rgn);
    const std::vector<std::string> rgn_names = boolrgn->get_component_regions();
    int nregs = rgn_names.size();

    std::vector<std::set<Entity_ID>> msets;

    for (int r = 0; r < nregs; r++) {
      auto rgn1 = gm->FindRegion(rgn_names[r]);

      // Did not find the rgn
      if (rgn1 == Teuchos::null) {
        std::stringstream ss;
        ss << "Geometric model has no region named " << rgn_names[r];
        Errors::Message msg(ss.str());
        Exceptions::amanzi_throw(msg);
      }

      std::string setname_internal1 = rgn1->get_name() + std::to_string(kind);

      if (sets_.find(setname_internal) == sets_.end()) {
        Entity_ID_List setents1;
        get_set_entities_and_vofs(rgn1->get_name(), kind, ptype, &setents1, vofs);
        sets_[setname_internal1] = setents1;
      }

      auto it = sets_.find(setname_internal1);
      msets.push_back(std::set<Entity_ID>(it->second.begin(), it->second.end()));
    }

    if (boolrgn->get_operation() == AmanziGeometry::BoolOpType::UNION) {
      for (int n = 1; n < msets.size(); ++n) {
        for (auto it = msets[n].begin(); it != msets[n].end(); ++it) msets[0].insert(*it);
      }
    }

    for (auto it = msets[0].begin(); it != msets[0].end(); ++it) setents->push_back(*it);
  }
}


//---------------------------------------------------------
// Deform a mesh so that cell volumes conform as closely as possible
// to target volumes without dropping below the minimum volumes.  If
// move_vertical = true, nodes will be allowed to move only in the
// vertical direction (right now arbitrary node movement is not allowed)
//---------------------------------------------------------
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


//---------------------------------------------------------
// Write mesh out to exodus file
//---------------------------------------------------------
void
Mesh_simple::write_to_exodus_file(const std::string filename) const
{
  Errors::Message mesg("Not implemented");
  amanzi_throw(mesg);
}


//---------------------------------------------------------
// Geometric center of a point cloud
//---------------------------------------------------------
AmanziGeometry::Point
Mesh_simple::geometric_center_(std::vector<AmanziGeometry::Point>& vxyz) const
{
  AmanziGeometry::Point gp(space_dimension());
  for (int j = 0; j < vxyz.size(); j++) gp += vxyz[j];
  gp /= vxyz.size();

  return gp;
}


//---------------------------------------------------------
// Check that all points are inside region
//---------------------------------------------------------
bool
Mesh_simple::all_inside_(std::vector<AmanziGeometry::Point>& vxyz,
                         const AmanziGeometry::Region& rgn) const
{
  for (int j = 0; j < vxyz.size(); j++) {
    if (!rgn.inside(vxyz[j])) return false;
  }
  return true;
}

} // namespace AmanziMesh
} // namespace Amanzi

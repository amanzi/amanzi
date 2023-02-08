/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Mesh

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
Mesh_simple::Mesh_simple(double x0, double y0, double z0,
                         double x1, double y1, double z1,
                         int nx, int ny, int nz,
                         const Comm_ptr_type& comm,
                         const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
                         const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : MeshFramework(comm, gm, plist),
    nx_(nx), ny_(ny), nz_(nz),
    x0_(x0), x1_(x1),
    y0_(y0), y1_(y1),
    z0_(z0), z1_(z1)
{
  setSpaceDimension(3);
  setManifoldDimension(3);
  edges_requested_ = plist_->get<bool>("request edges", false);
  CreateCache_();
}


//---------------------------------------------------------
// Update
//---------------------------------------------------------
void Mesh_simple::CreateCache_()
{

  // build new cache
  num_cells_ = nx_ * ny_ * nz_;
  num_nodes_ = (nx_ + 1) * (ny_ + 1) * (nz_ + 1);
  num_faces_ = (nx_ + 1) * ny_ * nz_ + nx_ * (ny_ + 1) * nz_ + nx_ * ny_ * (nz_ + 1);
  num_edges_ = nx_ * (ny_ + 1) * (nz_ + 1) + (nx_ + 1) * ny_ * (nz_ + 1) + (nx_ + 1) * (ny_ + 1) * nz_;

  // -- node coordinates
  Kokkos::resize(coordinates_, 3 * num_nodes_);

  double hx = (x1_ - x0_) / nx_;
  double hy = (y1_ - y0_) / ny_;
  double hz = (z1_ - z0_) / nz_;

  for (int iz = 0; iz <= nz_; iz++) {
    for (int iy = 0; iy <= ny_; iy++) {
      for (int ix = 0; ix <= nx_; ix++) {
        int istart = 3 * node_index_(ix, iy, iz);
        coordinates_[istart]     = x0_ + ix * hx;
        coordinates_[istart + 1] = y0_ + iy * hy;
        coordinates_[istart + 2] = z0_ + iz * hz;
      }
    }
  }

  // -- connectivity arrays
  Kokkos::resize(cell_to_face_,6 * num_cells_);
  Kokkos::resize(cell_to_face_dirs_, 6 * num_cells_);
  Kokkos::resize(face_to_cell_, 2 * num_faces_);
  initView(face_to_cell_, -1);  

  Kokkos::resize(face_to_node_,4 * num_faces_);
  Kokkos::resize(node_to_face_,13 * num_nodes_);  // 1 extra for num faces

  if (edges_requested_) {
    Kokkos::resize(face_to_edge_, 4 * num_faces_);
    Kokkos::resize(face_to_edge_dirs_, 4 * num_faces_);
    Kokkos::resize(edge_to_node_, 2 * num_edges_);
  }

  // loop over cells and initialize cell <-> face
  for (int iz = 0; iz < nz_; iz++) {
    for (int iy = 0; iy < ny_; iy++) {
      for (int ix = 0; ix < nx_; ix++) {
        int istart = 6 * cell_index_(ix,iy,iz);
        int jstart = 0;

        cell_to_face_[istart]     = xzface_index_(ix,  iy,  iz);
        cell_to_face_[istart + 1] = yzface_index_(ix+1,iy,  iz);
        cell_to_face_[istart + 2] = xzface_index_(ix,  iy+1,iz);
        cell_to_face_[istart + 3] = yzface_index_(ix,  iy,  iz);
        cell_to_face_[istart + 4] = xyface_index_(ix,  iy,  iz);
        cell_to_face_[istart + 5] = xyface_index_(ix,  iy,  iz+1);

        cell_to_face_dirs_[istart]     = 1;
        cell_to_face_dirs_[istart + 1] = 1;
        cell_to_face_dirs_[istart + 2] = -1;
        cell_to_face_dirs_[istart + 3] = -1;
        cell_to_face_dirs_[istart + 4] = -1;
        cell_to_face_dirs_[istart + 5] = 1;

        jstart = 2 * xzface_index_(ix, iy, iz);
        face_to_cell_[jstart + 1] = cell_index_(ix, iy, iz);

        jstart = 2 * yzface_index_(ix+1, iy, iz);
        face_to_cell_[jstart + 1] = cell_index_(ix, iy, iz);

        jstart = 2 * xzface_index_(ix, iy+1, iz);
        face_to_cell_[jstart] = cell_index_(ix, iy, iz);

        jstart = 2 * yzface_index_(ix, iy, iz);
        face_to_cell_[jstart] = cell_index_(ix, iy, iz);

        jstart = 2 * xyface_index_(ix, iy, iz);
        face_to_cell_[jstart] = cell_index_(ix, iy, iz);

        jstart = 2 * xyface_index_(ix, iy, iz+1);
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

        face_to_node_[istart]     = node_index_(ix,  iy,  iz);
        face_to_node_[istart + 1] = node_index_(ix+1,iy,  iz);
        face_to_node_[istart + 2] = node_index_(ix+1,iy+1,iz);
        face_to_node_[istart + 3] = node_index_(ix,  iy+1,iz);

        jstart = 13 * node_index_(ix,iy,iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = 13 * node_index_(ix+1, iy, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = 13 * node_index_(ix+1, iy+1, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = 13 * node_index_(ix, iy+1, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;
      }
    }
  }

  // -- xz faces
  for (int iz = 0; iz < nz_; iz++) {
    for (int iy = 0; iy <= ny_; iy++) {
      for (int ix=0; ix < nx_; ix++) {
        int istart = 4 * xzface_index_(ix, iy, iz);
        int jstart = 0;
        int nfaces = 0;

        face_to_node_[istart]     = node_index_(ix,  iy, iz);
        face_to_node_[istart + 1] = node_index_(ix+1,iy, iz);
        face_to_node_[istart + 2] = node_index_(ix+1,iy, iz+1);
        face_to_node_[istart + 3] = node_index_(ix,  iy, iz+1);

        jstart = 13 * node_index_(ix, iy, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = 13 * node_index_(ix+1, iy, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = 13 * node_index_(ix+1, iy, iz+1);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = 13 * node_index_(ix, iy, iz+1);
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

        face_to_node_[istart]     = node_index_(ix, iy,  iz);
        face_to_node_[istart + 1] = node_index_(ix, iy+1,iz);
        face_to_node_[istart + 2] = node_index_(ix, iy+1,iz+1);
        face_to_node_[istart + 3] = node_index_(ix, iy,  iz+1);

        jstart = 13 * node_index_(ix, iy, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = 13 * node_index_(ix, iy+1, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = 13 * node_index_(ix, iy+1, iz+1);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = 13 * node_index_(ix, iy, iz+1);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart + 1 + nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;
      }
    }
  }

  if (edges_requested_) { 
    // loop over faces and initialize face -> edge
    // -- xy faces
    for (int iz = 0; iz <= nz_; iz++) {
      for (int iy = 0; iy < ny_; iy++) {
        for (int ix = 0; ix < nx_; ix++) {
          int istart = 4 * xyface_index_(ix, iy, iz);

          face_to_edge_[istart]     = xedge_index_(ix,  iy,  iz);
          face_to_edge_[istart + 1] = yedge_index_(ix+1,iy,  iz);
          face_to_edge_[istart + 2] = xedge_index_(ix,  iy+1,iz);
          face_to_edge_[istart + 3] = yedge_index_(ix,  iy,  iz);

          face_to_edge_dirs_[istart]     = 1;
          face_to_edge_dirs_[istart + 1] = 1;
          face_to_edge_dirs_[istart + 2] = -1;
          face_to_edge_dirs_[istart + 3] = -1;
        }
      }
    }

    // -- xz faces
    for (int iz = 0; iz < nz_; iz++) {
      for (int iy = 0; iy <= ny_; iy++) {
        for (int ix=0; ix < nx_; ix++) {
          int istart = 4 * xzface_index_(ix, iy, iz);

          face_to_edge_[istart]     = xedge_index_(ix,  iy, iz);
          face_to_edge_[istart + 1] = zedge_index_(ix+1,iy, iz);
          face_to_edge_[istart + 2] = xedge_index_(ix,  iy, iz+1);
          face_to_edge_[istart + 3] = zedge_index_(ix,  iy, iz);

          face_to_edge_dirs_[istart]     = 1;
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

          face_to_edge_[istart]     = yedge_index_(ix, iy,  iz);
          face_to_edge_[istart + 1] = zedge_index_(ix, iy+1,iz);
          face_to_edge_[istart + 2] = yedge_index_(ix, iy,  iz+1);
          face_to_edge_[istart + 3] = zedge_index_(ix, iy,  iz);

          face_to_edge_dirs_[istart]     = 1;
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

          edge_to_node_[istart]     = node_index_(ix,  iy,  iz);
          edge_to_node_[istart + 1] = node_index_(ix+1,iy,  iz);
        }
      }
    }

    // -- y edges
    for (int iz = 0; iz <= nz_; iz++) {
      for (int iy = 0; iy < ny_; iy++) {
        for (int ix = 0; ix <= nx_; ix++) {
          int istart = 2 * yedge_index_(ix, iy, iz);

          edge_to_node_[istart]     = node_index_(ix, iy,   iz);
          edge_to_node_[istart + 1] = node_index_(ix, iy+1, iz);
        }
      }
    }

    // -- z edges
    for (int iz = 0; iz < nz_; iz++) {
      for (int iy = 0; iy <= ny_; iy++) {
        for (int ix = 0; ix <= nx_; ix++) {
          int istart = 2 * zedge_index_(ix, iy, iz);

          edge_to_node_[istart]     = node_index_(ix, iy, iz);
          edge_to_node_[istart + 1] = node_index_(ix, iy, iz+1);
        }
      }
    }
  }
}


//---------------------------------------------------------
// TBW
//---------------------------------------------------------
std::size_t Mesh_simple::getNumEntities(AmanziMesh::Entity_kind kind,
                                       AmanziMesh::Parallel_kind ptype) const
{
  switch (kind) {
    case Entity_kind::FACE:
      return (ptype != AmanziMesh::Parallel_kind::GHOST) ? num_faces_ : 0;
      break;
    case Entity_kind::NODE:
      return (ptype != AmanziMesh::Parallel_kind::GHOST) ? num_nodes_ : 0;
      break;
    case Entity_kind::CELL:
      return (ptype != AmanziMesh::Parallel_kind::GHOST) ? num_cells_ : 0;
      break;
    case Entity_kind::EDGE: 
      return (ptype != AmanziMesh::Parallel_kind::GHOST) ? num_edges_ : 0;
      break;
    default:
      throw std::exception();
      break;
  }
}


//---------------------------------------------------------
// Connectivity: cell -> faces
//---------------------------------------------------------
void Mesh_simple::getCellFacesAndDirs(const AmanziMesh::Entity_ID cellid,
        Entity_ID_View& faceids,
        Entity_Direction_View *cfacedirs) const
{
  unsigned int offset = (unsigned int) 6*cellid;

  Kokkos::resize(faceids, 6);
  for(int i = 0 ; i < 6 ; ++i) faceids[i] = cell_to_face_[offset+i];  

  if (cfacedirs) {
    Kokkos::resize(*cfacedirs, 6);
    for(int i = 0 ; i < 6 ; ++i) (*cfacedirs)[i] = cell_to_face_dirs_[offset+i];  
  }
}


//---------------------------------------------------------
// Connectivity: face -> nodes
//---------------------------------------------------------
void Mesh_simple::getFaceNodes(AmanziMesh::Entity_ID face,
        AmanziMesh::Entity_ID_View& nodeids) const
{
  unsigned int offset = (unsigned int) 4*face;
  Kokkos::resize(nodeids, 4); 
  for (int i = 0; i < 4; i++) {
    nodeids[i] = face_to_node_[offset+i];
  }
}


//---------------------------------------------------------
// Connectivity: face -> edges
//---------------------------------------------------------
void Mesh_simple::getFaceEdgesAndDirs(const Entity_ID faceid,
        Entity_ID_View& edgeids,
        Entity_Direction_View *fedgedirs) const
{
  unsigned int offset = (unsigned int) 4*faceid;

  Kokkos::resize(edgeids, 3); 
  for(int i = 0 ; i < 3 ; ++i) edgeids[i] = face_to_edge_[offset+i];  

  if (fedgedirs) {
    Kokkos::resize(*fedgedirs, 3); 
    for(int i = 0 ; i < 3 ; ++i) (*fedgedirs)[i] = face_to_edge_dirs_[offset+i];  
  }
}


//---------------------------------------------------------
// Connectivity: edge -> nodes
//---------------------------------------------------------
void Mesh_simple::getEdgeNodes(
  const Entity_ID edgeid, Entity_ID_View& nodes) const
{
  unsigned int offset = (unsigned int) 2*edgeid;
  Kokkos::resize(nodes, 2);
  nodes[0] = edge_to_node_[offset];
  nodes[1] = edge_to_node_[offset + 1];
}


//---------------------------------------------------------
// Cooordinates of a node
//---------------------------------------------------------
AmanziGeometry::Point
Mesh_simple::getNodeCoordinate(const AmanziMesh::Entity_ID local_node_id) const
{
  unsigned int offset = (unsigned int) 3*local_node_id;
  AmanziGeometry::Point ncoord;
  ncoord.set(3, &(coordinates_[offset]));
  return ncoord;
}


//---------------------------------------------------------
// Modify cooordinates of a node: version 2
//---------------------------------------------------------
void Mesh_simple::setNodeCoordinate(const AmanziMesh::Entity_ID local_node_id,
        const AmanziGeometry::Point& ncoord)
{
  unsigned int offset = (unsigned int) 3*local_node_id;

  int spdim = getSpaceDimension();
  for (int i = 0; i < spdim; i++) {
    coordinates_[offset+i] = ncoord[i]; 
  }
}


//---------------------------------------------------------
// Faces of type 'ptype' connected to a node
//---------------------------------------------------------
void Mesh_simple::getNodeFaces(const AmanziMesh::Entity_ID nodeid,
        const AmanziMesh::Parallel_kind ptype,
        AmanziMesh::Entity_ID_View& faceids) const
{
  unsigned int offset = (unsigned int) 13*nodeid;
  unsigned int nfaces = node_to_face_[offset];

  Kokkos::resize(faceids, nfaces); 

  for (int i = 0; i < nfaces; i++) 
    faceids[i] = node_to_face_[offset+i+1];
}


//---------------------------------------------------------
// Cells connected to a face
//---------------------------------------------------------
void Mesh_simple::getFaceCells(const AmanziMesh::Entity_ID faceid,
        const AmanziMesh::Parallel_kind ptype,
        AmanziMesh::Entity_ID_View& cellids) const
{
  unsigned int offset = (unsigned int) 2*faceid;

  Kokkos::resize(cellids, 2); 
  int cellids_ct = 0;
  if (face_to_cell_[offset] != -1)
    cellids[cellids_ct++] = face_to_cell_[offset];
  if (face_to_cell_[offset+1] != -1)
    cellids[cellids_ct++] = face_to_cell_[offset+1];
  Kokkos::resize(cellids, cellids_ct); 
}


}  // namespace AmanziMesh
}  // namespace Amanzi


/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//! A very simple structured 3D mesh, mostly for testing and serial applications.
/*

Note this implements the bare minimum functionality, and is pretty
inefficiently implemented.

*/

#ifndef AMANZI_MESH_SIMPLE_HH_
#define AMANZI_MESH_SIMPLE_HH_

#include <AmanziComm.hh>

#include <memory>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "MeshFramework.hh"
#include "Region.hh"

#include "GeometricModel.hh"
#include "VerboseObject.hh"

namespace Amanzi {
namespace AmanziMesh {

class Mesh_simple : public MeshFramework {
 public:
  // the request_faces and request_edges arguments have to be at the
  // end and not in the middle because if we omit them and specify a
  // pointer argument like gm or verbosity_obj, then there is implicit
  // conversion of the pointer to bool, thereby defeating the intent
  // of the call and making the pointer argument seem NULL. In C++11,
  // we could "delete" the illegal version of the call effectively
  // blocking the implicit conversion.
  Mesh_simple(double x0, double y0, double z0,
              double x1, double y1, double z1,
              int nx, int ny, int nz,
              const Comm_ptr_type& comm,
              const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm=Teuchos::null,
              const Teuchos::RCP<Teuchos::ParameterList>& plist=Teuchos::null);

  virtual ~Mesh_simple() = default;

  virtual bool hasEdges() const override { return edges_requested_; }

  // Number of entities of any kind (cell, face, node) and in a
  // particular category (OWNED, GHOST, ALL)
  virtual std::size_t getNumEntities(const Entity_kind kind,
                            const Parallel_type ptype) const override;


  // Node coordinates - 3 in 3D and 2 in 2D
  virtual AmanziGeometry::Point getNodeCoordinate(const Entity_ID nodeid) const override;
  virtual void setNodeCoordinate(const Entity_ID nodeid,
          const AmanziGeometry::Point& coord) override;

  virtual void getCellFacesAndDirs(
    const Entity_ID c,
    Entity_ID_List& faces,
    Entity_Direction_List * const dirs) const override;

  virtual void getFaceEdgesAndDirs(const Entity_ID f,
          Entity_ID_List& edges,
          Entity_Direction_List * const dirs=nullptr) const override;

  virtual void getFaceNodes(const Entity_ID f, Entity_ID_List& nodes) const override;

  virtual void getEdgeNodes(const Entity_ID e, Entity_ID_List& nodes) const override;


  virtual void getFaceCells(const Entity_ID f,
                          const Parallel_type ptype,
                          Entity_ID_List& cells) const override;

  virtual void getNodeFaces(const Entity_ID nodeid,
                            const Parallel_type ptype,
                            Entity_ID_List& faceids) const override;

 private:
  void CreateCache_();

private:
  bool edges_requested_;
  Double_List coordinates_;

  unsigned int node_index_(int i, int j, int k) const;
  unsigned int cell_index_(int i, int j, int k) const;

  unsigned int xyface_index_(int i, int j, int k) const;
  unsigned int yzface_index_(int i, int j, int k) const;
  unsigned int xzface_index_(int i, int j, int k) const;

  unsigned int xedge_index_(int i, int j, int k) const;
  unsigned int yedge_index_(int i, int j, int k) const;
  unsigned int zedge_index_(int i, int j, int k) const;

  int nx_, ny_, nz_;  // number of cells in the three coordinate directions
  double x0_, x1_, y0_, y1_, z0_, z1_;  // coordinates of lower left front and upper right back of brick

  int num_cells_, num_faces_, num_edges_, num_nodes_;

  // mesh connectivity arrays
  Entity_ID_List cell_to_face_;
  Entity_ID_List face_to_edge_;
  Entity_ID_List face_to_node_;
  Entity_ID_List edge_to_node_;

  Entity_ID_List node_to_face_;
  Entity_ID_List node_to_edge_;
  Entity_ID_List edge_to_face_;
  Entity_ID_List face_to_cell_;

  // orientation arrays
  Entity_Direction_View cell_to_face_dirs_;
  Entity_Direction_View face_to_edge_dirs_;
};


// -------------------------
// Template & inline members
// ------------------------
inline
unsigned int Mesh_simple::node_index_(int i, int j, int k) const {
  return i + j * (nx_ + 1) + k * (nx_ + 1) * (ny_ + 1);
}

inline
unsigned int Mesh_simple::cell_index_(int i, int j, int k) const {
  return i + j * nx_ + k * nx_ * ny_;
}

inline
unsigned int Mesh_simple::xyface_index_(int i, int j, int k) const {
  return i + j * nx_ + k * nx_ * ny_;
}

inline
unsigned int Mesh_simple::xzface_index_(int i, int j, int k) const {
  return i + j * nx_ + k * nx_ * (ny_ + 1) + xyface_index_(0, 0, nz_ + 1);
}

inline
unsigned int Mesh_simple::yzface_index_(int i, int j, int k) const {
  return i + j * (nx_ + 1) + k * (nx_ + 1) * ny_ + xzface_index_(0, 0, nz_);
}

inline
unsigned int Mesh_simple::xedge_index_(int i, int j, int k) const {
  return i + j * nx_ + k * nx_ * (ny_ + 1);
}

inline
unsigned int Mesh_simple::yedge_index_(int i, int j, int k) const {
  return i + j * (nx_ + 1) + k * (nx_ + 1) * ny_ + xedge_index_(0, 0, nz_ + 1);
}

inline
unsigned int Mesh_simple::zedge_index_(int i, int j, int k) const {
  return i + j * (nx_ + 1) + k * (nx_ + 1) * (ny_ + 1) + yedge_index_(0, 0, nz_ + 1);
}

}  // namespace AmanziMesh
}  // namespace Amanzi

#endif


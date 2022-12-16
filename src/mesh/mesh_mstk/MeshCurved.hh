/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Experimental mesh class for 2D meshes with curved faces.
*/

#ifndef AMANZI_MESH_CURVED_HH_
#define AMANZI_MESH_CURVED_HH_

#include <memory>

#include "Mesh_MSTK.hh"

namespace Amanzi {

namespace AmanziMesh {

class MeshCurved : public Mesh_MSTK {
 public:
  MeshCurved(double x0,
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
             const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : Mesh_MSTK(x0, y0, z0, x1, y1, z1, nx, ny, nz, comm, gm, plist, true, true),
      edge_ho_nodes_(nullptr),
      face_ho_nodes_(nullptr){};

  MeshCurved(double x0,
             double y0,
             double x1,
             double y1,
             int nx,
             int ny,
             const Comm_ptr_type& comm,
             const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
             const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : Mesh_MSTK(x0, y0, x1, y1, nx, ny, comm, gm, plist, true, false), face_ho_nodes_(nullptr){};

  MeshCurved(const std::string& filename,
             const Comm_ptr_type& comm,
             const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
             const Teuchos::RCP<Teuchos::ParameterList>& plist,
             bool request_faces = true,
             bool request_edges = false)
    : Mesh_MSTK(filename.c_str(), comm, gm, plist, request_faces, request_edges),
      face_ho_nodes_(nullptr){};

  ~MeshCurved(){};

  // new implementtion of some basis functions
  // -- volume/area of cell
  double cell_volume(const Entity_ID c, const bool recompute = false) const;
  double cell_volume_linear(const Entity_ID c) const { return Mesh::cell_volume(c); }

  // -- curvature information
  virtual void edge_get_ho_nodes(Entity_ID e, AmanziGeometry::Point_List* nodes) const override
  {
    if (edge_ho_nodes_ != nullptr) *nodes = (*edge_ho_nodes_)[e];
  }
  virtual void face_get_ho_nodes(Entity_ID f, AmanziGeometry::Point_List* nodes) const override
  {
    if (face_ho_nodes_ != nullptr) *nodes = (*face_ho_nodes_)[f];
  }

  // define a curved mesh using additional edge and face points
  void set_edge_ho_nodes(std::shared_ptr<std::vector<AmanziGeometry::Point_List>> edge_ho_nodes)
  {
    edge_ho_nodes_ = edge_ho_nodes;
  }
  void set_face_ho_nodes(std::shared_ptr<std::vector<AmanziGeometry::Point_List>> face_ho_nodes)
  {
    face_ho_nodes_ = face_ho_nodes;
  }

 private:
  std::shared_ptr<std::vector<AmanziGeometry::Point_List>> edge_ho_nodes_;
  std::shared_ptr<std::vector<AmanziGeometry::Point_List>> face_ho_nodes_;
};


inline double
MeshCurved::cell_volume(const Entity_ID c, const bool recompute) const
{
  if (space_dimension() == 3) return cell_volume_linear(c); // FIXME

  Entity_ID_List nodes, faces;
  AmanziGeometry::Point_List points;

  cell_get_nodes(c, &nodes);
  cell_get_faces(c, &faces);
  int nfaces = faces.size();

  double volume(0.0);
  AmanziGeometry::Point x0(2), x1(2), xf(2), tau(2), q2(2);

  for (int n = 0; n < nfaces; ++n) {
    int m = (n + 1) % nfaces;
    node_get_coordinates(nodes[n], &x0);
    node_get_coordinates(nodes[m], &x1);

    face_get_ho_nodes(faces[n], &points);
    int npoints = points.size();

    xf = (x1 + x0) / 2;
    tau = x1 - x0;
    volume += (xf[0] * tau[1] - xf[1] * tau[0]) / 2;

    if (npoints == 1) {
      q2 = points[0] - xf;
      volume += 2 * (q2[0] * tau[1] - q2[1] * tau[0]) / 3;
    }
  }

  return volume;
}

} // namespace AmanziMesh
} // namespace Amanzi
#endif

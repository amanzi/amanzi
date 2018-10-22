/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

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
  MeshCurved(double x0, double y0, double x1, double y1,
	     int nx, int ny, const Epetra_MpiComm *comm)
      : Mesh_MSTK(x0, y0, x1, y1, nx, ny, comm, 
                  Teuchos::null, Teuchos::null, true, false, PARTITIONER_DEFAULT),
        face_points_(NULL) {};
  ~MeshCurved() {};

  virtual 
  void face_get_curvature_points(Entity_ID f,
                                 const AmanziGeometry::Point_List* points) const override {
    points = &(*face_points_)[f];
  }

  // define a curved mesh using additional face points
  void set_face_points(std::shared_ptr<const std::vector<AmanziGeometry::Point_List> > face_points) {
    face_points_ = face_points;
  }

 private:
  std::shared_ptr<const std::vector<AmanziGeometry::Point_List> > face_points_;
};

}  // namespace AmanziMesh
}  // namesace Amanzi
#endif

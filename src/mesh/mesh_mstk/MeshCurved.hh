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
	     int nx, int ny, const Epetra_MpiComm *comm, Partitioner_type partitioner)
      : Mesh_MSTK(x0, y0, x1, y1, nx, ny, comm, 
                  Teuchos::null, Teuchos::null, true, false, partitioner),
        face_ho_nodes_(NULL) {};

  MeshCurved(const std::string& filename, const Epetra_MpiComm *comm, Partitioner_type partitioner)
      : Mesh_MSTK(filename.c_str(), comm, Teuchos::null, Teuchos::null, true, false, partitioner),
        face_ho_nodes_(NULL) {};

  ~MeshCurved() {};

  virtual 
  void face_get_ho_nodes(Entity_ID f,
                         AmanziGeometry::Point_List* nodes) const override {
    if (face_ho_nodes_ != NULL) *nodes = (*face_ho_nodes_)[f];
  }

  // define a curved mesh using additional face points
  void set_face_ho_nodes(std::shared_ptr<std::vector<AmanziGeometry::Point_List> > face_ho_nodes) {
    face_ho_nodes_ = face_ho_nodes;
  }

 private:
  std::shared_ptr<std::vector<AmanziGeometry::Point_List> > face_ho_nodes_;
};

}  // namespace AmanziMesh
}  // namesace Amanzi
#endif

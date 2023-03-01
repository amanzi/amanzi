/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, others
*/
/*

This is a semi-pseudo-1D, simlified mesh with a vertical column of prismatic
cells (the horizontal faces can be polygonal). Users of this class must note
several important assumptions and simplifications made in its implementation.

Preconditions:

1. It must be possible to call buildColumns() on the parent mesh.

2. All cells in the column are assumed to have the same prismatic topology
(horizontal faces can have polygonal topology, lateral faces are all quads).
Pinchouts are not allowed, but this is probably not enforced until you try to
deform.  Violating this may cause simply wrong information without errors.

Postconditions of this mesh:

1. The lateral faces of the cells are ignored - so each cell is
considered to have only two faces, a bottom face and a top face.

2. The normals of the two faces are pointing vertically down or up.

3. The base face of the column's nodal coordinates are enforced to be perfectly
horizontal.

4. The X-Y coordinates of nodes above the base face are enforced so that they
are perfectly stacked on top of the nodes of the base face.  The Z-coordinates
of the nodes of a face are adjusted so that they match the Z-coordinate of the
centroid of the original face.

5. The volumes of the resulting cells may or may not have the same volume as
the parent mesh.  They will have the same volume if and only if the parent mesh
is "terrain following."

*/

#ifndef AMANZI_MESH_COLUMN_HH_
#define AMANZI_MESH_COLUMN_HH_

#include <memory>
#include <vector>
#include <string>

#include "Teuchos_ParameterList.hpp"
#include "Epetra_Map.h"
#include "AmanziComm.hh"
#include "Epetra_SerialComm.h"

#include "dbc.hh"
#include "errors.hh"

#include "MeshFramework.hh"

namespace Amanzi {
namespace AmanziMesh {

struct MeshFrameworkColumnAlgorithms : public MeshFrameworkAlgorithms {
  // lumped things for more efficient calculation
  virtual std::pair<double, AmanziGeometry::Point>
  computeCellGeometry(const MeshFramework& mesh, const Entity_ID c) const override;

  // // replicated because of a lack of templated virtual functions
  // virtual std::pair<double, AmanziGeometry::Point>
  // computeCellGeometry(const Mesh& mesh, const Entity_ID c) const override;
};


class MeshFrameworkColumn : public MeshFramework {
 public:
  MeshFrameworkColumn(const Teuchos::RCP<MeshFramework>& col3D_mesh,
             const Teuchos::RCP<Teuchos::ParameterList>& plist);

  virtual ~MeshFrameworkColumn() = default;

  // reference for vis.
  virtual Teuchos::RCP<const MeshFramework> getParentMesh() const override {
    return col3D_mesh_->getParentMesh();
  }

  virtual bool isDeformable() const override { return true; }

  // Number of entities of any kind (cell, face, node) and in a
  // particular category (OWNED, GHOST, ALL)
  virtual
  std::size_t getNumEntities(const Entity_kind kind,
                             const Parallel_kind ptype) const override {
    std::size_t count;
    switch (kind) {
      case Entity_kind::FACE : {
        count = (ptype == Parallel_kind::GHOST) ? 0 : column_faces_.size();
        break;
      }
      case Entity_kind::BOUNDARY_FACE : {
        count = 2;
        break;
      }
      default : {
        count = col3D_mesh_->getNumEntities(kind, ptype);
        break;
      }
    }
    return count;
  }

  // Parent entity in the source mesh if mesh was derived from another mesh
  virtual Entity_ID getEntityParent(const Entity_kind kind, const Entity_ID lid) const override {
    Entity_ID ent = -1;
    switch (kind) {
      case Entity_kind::FACE:
        ent = col3D_mesh_->getEntityParent(kind, column_faces_(lid));
        break;
      default:
        ent = col3D_mesh_->getEntityParent(kind, lid);
        break;
    }
    return ent;
  }

  // Get cell type - UNKNOWN, TRI, QUAD, ... See MeshDefs.hh
  virtual
  Cell_kind getCellType(const Entity_ID lid) const override {
    return col3D_mesh_->getCellType(lid);
  }

  virtual AmanziGeometry::Point getNodeCoordinate(const Entity_ID node) const override {
    return col3D_mesh_->getNodeCoordinate(node);
  }

  // Get nodes of face
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal
  // In 2D, nfnodes is 2
  virtual
  void getFaceNodes(const Entity_ID faceid,
                    cEntity_ID_View& nodeids) const override {
    col3D_mesh_->getFaceNodes(column_faces_(faceid), nodeids);
  }

  // Upward adjacencies
  //-------------------

  // Cells of type 'ptype' connected to a node - The order of cells is
  // not guaranteed to be the same for corresponding nodes on
  // different processors
  virtual
  void getNodeCells(const Entity_ID nodeid,
                    const Parallel_kind ptype,
                    cEntity_ID_View& cellids) const override {
    col3D_mesh_->getNodeCells(nodeid, ptype, cellids);
  }


  // Faces of type 'ptype' connected to a node - The order of faces is
  // not guarnateed to be the same for corresponding nodes on
  // different processors
  virtual
  void getNodeFaces(const Entity_ID nodeid,
                    const Parallel_kind ptype,
                    cEntity_ID_View& faceids) const override {
    Errors::Message mesg("Not implemented");
    Exceptions::amanzi_throw(mesg);
  }

  //
  // Mesh modification
  //-------------------
  virtual
  void setNodeCoordinate(const Entity_ID nodeid,
                         const AmanziGeometry::Point& ncoord) override {
    col3D_mesh_->setNodeCoordinate(nodeid, ncoord);
  }

  // get faces and face dirs of a cell. This can be called by
  // cell_get_faces_and_dirs method of the base class and the data
  // cached or it can be called directly by the
  // cell_get_faces_and_dirs method of this class
  virtual
  void getCellFacesAndDirs(const Entity_ID cellid,
                           cEntity_ID_View& faceids,
                           cEntity_Direction_View * const face_dirs) const override
  {
    Entity_ID_View lfaceids("lfaceids",2); 
    Entity_Direction_View lface_dirs; 
    if (face_dirs) Kokkos::resize(lface_dirs,2);

    // NOTE: the face directions with respect to the cell may be at
    // odds with how it is in the parent mesh but within this mesh its
    // consistent - so we think everything will work as it should
    cEntity_ID_View faceids_extracted;
    cEntity_Direction_View face_dirs_extracted;
    col3D_mesh_->getCellFacesAndDirs(cellid, faceids_extracted,
            &face_dirs_extracted);
            
    int count = 0;
    for (int i=0; i!=faceids_extracted.size(); ++i) {
      if (face_in_column_[faceids_extracted[i]] >= 0) {
        lfaceids[count] = face_in_column_[faceids_extracted[i]];
        if (face_dirs) lface_dirs[count] = face_dirs_extracted[i];
        count++;
      }
    }
    faceids = lfaceids; 
    if(face_dirs) *face_dirs = lface_dirs; 
  }

  // Cells connected to a face - this function is implemented in each
  // mesh framework. The results are cached in the base class
  virtual
  void getFaceCells(const Entity_ID faceid,
                    const Parallel_kind ptype,
                    cEntity_ID_View& cellids) const override {
    col3D_mesh_->getFaceCells(column_faces_(faceid), ptype, cellids);
  }

 protected:
  void computeSpecialNodeCoordinates_();

 protected:
  Teuchos::RCP<MeshFramework> col3D_mesh_;
  int nfnodes_;
  int column_id_;
  Kokkos::MeshView<Entity_ID*, Kokkos::DefaultHostExecutionSpace> column_faces_;
  Entity_ID_View face_in_column_;
};


namespace MeshAlgorithms {

// helper function
template<class Mesh_type>
std::pair<double,AmanziGeometry::Point>
computeMeshColumnCellGeometry(const Mesh_type& mesh, const Entity_ID c)
{
  /* compute volume on the assumption that the top and bottom faces form
     a vertical columnar cell or in other words a polygonal prism */
  cEntity_ID_View cfaces;
  mesh.getCellFaces(c, cfaces);
  AMANZI_ASSERT(cfaces.size() == 2);
  double farea = mesh.getFaceArea(cfaces[0]);

  auto f0 = mesh.getFaceCentroid(cfaces[0]);
  auto f1 = mesh.getFaceCentroid(cfaces[1]);
  double height = std::abs(f0[2] - f1[2]);
  double volume = height * farea;
  auto cc = (f0 + f1)/2.;
  return std::make_pair(volume, cc);
}


} // namespace MeshAlgorithms
} // namespace AmanziMesh
} // namespace Amanzi

#endif

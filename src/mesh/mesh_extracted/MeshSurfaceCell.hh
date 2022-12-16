/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! This is a mesh for a single surface cell.

/*!

  This exists solely because we need "surface meshes" extracted from
  a MeshColumn.  This is really just 1 2D cell.  Really.
*/

#ifndef AMANZI_MESH_SURFACE_CELL_HH_
#define AMANZI_MESH_SURFACE_CELL_HH_

#include <vector>
#include <string>
#include <algorithm>

#include "Teuchos_ParameterList.hpp"
#include "Epetra_Map.h"
#include "AmanziComm.hh"
#include "Epetra_SerialComm.h"

#include "VerboseObject.hh"
#include "dbc.hh"
#include "errors.hh"

#include "Region.hh"
#include "Mesh.hh"

namespace Amanzi {
namespace AmanziMesh {

class MeshSurfaceCell : public MeshFramework {
 public:
  MeshSurfaceCell(const Teuchos::RCP<const MeshFramework>& col_mesh,
                  bool flatten=true);

  virtual ~MeshSurfaceCell() = default;

  // Parent entity in the source mesh if mesh was derived from another mesh
  virtual
  Entity_ID getEntityParent(const Entity_kind kind, const Entity_ID entid) const override {
    if (kind == Entity_kind::CELL) {
      AMANZI_ASSERT(entid == 0);
      return parent_face_;
    } else if (kind == Entity_kind::NODE) {
      return node_parents_[entid];
    } else {
      Errors::Message msg;
      msg << "MeshSurfaceCell: Cannot getEntityParent() for Entity kind " << to_string(kind);
      Exceptions::amanzi_throw(msg);
    }
    return -1;
  }

  virtual Teuchos::RCP<const MeshFramework> getParentMesh() const override {
    return parent_; }

// Get cell type - UNKNOWN, TRI, QUAD, ... See MeshDefs.hh
  virtual
  Cell_type getCellType(const Entity_ID lid) const override {
    return cell_type_;
  }

  //
  // General mesh information
  // -------------------------
  //
  // Number of entities of any kind (cell, face, node) and in a
  // particular category (OWNED, GHOST, ALL)
  virtual
  std::size_t getNumEntities(const Entity_kind kind,
                             const Parallel_type ptype) const override;


  // Get nodes of a cell
  virtual
  void getCellNodes(const Entity_ID cellid,
                    Entity_ID_List& nodeids) const override;

  virtual
  void getFaceNodes(const Entity_ID faceid,
                    Entity_ID_List& nodeids) const override;


  // Upward adjacencies
  //-------------------

  // Cells of type 'ptype' connected to a node - The order of cells is
  // not guaranteed to be the same for corresponding nodes on
  // different processors
  virtual
  void getNodeCells(const Entity_ID nodeid,
                    const Parallel_type ptype,
                    Entity_ID_List& cellids) const override;

  // Faces of type 'ptype' connected to a node - The order of faces is
  // not guarnateed to be the same for corresponding nodes on
  // different processors
  virtual
  void getNodeFaces(const Entity_ID nodeid,
                    const Parallel_type ptype,
                    Entity_ID_List& faceids) const override;

  // Node coordinates - 3 in 3D and 2 in 2D
  virtual
  AmanziGeometry::Point getNodeCoordinate(const Entity_ID nodeid) const override {
    return nodes_[nodeid];
  }

  //
  // Mesh modification
  //-------------------

  // Set coordinates of node
  virtual
  void setNodeCoordinate(const Entity_ID nodeid,
                         const AmanziGeometry::Point& ncoord) override {
    nodes_[nodeid] = ncoord;
  }


  // get faces and face dirs of a cell. This can be called by
  // cell_get_faces_and_dirs method of the base class and the data
  // cached or it can be called directly by the
  // cell_get_faces_and_dirs method of this class
  virtual
  void getCellFacesAndDirs(const Entity_ID cellid,
                           Entity_ID_List& faceids,
                           Entity_Direction_List * const dirs) const override;

  // Cells connected to a face - this function is implemented in each
  // mesh framework. The results are cached in the base class
  virtual
  void getFaceCells(const Entity_ID faceid,
                    const Parallel_type ptype,
                    Entity_ID_List& cellids) const override;

 protected:
  Teuchos::RCP<const MeshFramework> parent_;

  Point_List nodes_;
  Entity_ID_List node_parents_;
  std::map<Set_ID,bool> sets_;
  Entity_ID parent_face_;
  Cell_type cell_type_;

};


} // close namespace AmanziMesh
} // close namespace Amanzi


#endif /* _MESH_MAPS_H_ */

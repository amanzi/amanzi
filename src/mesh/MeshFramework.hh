/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, others
*/

//! The interface for meshes provided by external frameworks.
/*!

Developer note:

Most of this interface is here for testing.  Very little is likely to be used
in the final code, because largely the interface will be used to generate a
MeshCache object.  The MeshCache will then provide the full interface using
fast access, non-virtual methods.

A new Framework really must supply only a handful of methods, but may choose to
provide more, as long as they are consistent.

Note that the framework is split into two classes, MeshFramework and
MeshFrameworkAlgorithms, both of which must exist.  For many, the algorithms
will be the default class.  But some "special" frameworks may implement special
algorithms.  If they do so, the MeshCache object will need the algorithms even
if the Framework itself is deleted (hence the split).

*/

#pragma once

#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "GeometryDefs.hh"
#include "MeshDefs.hh"
#include "AmanziComm.hh"
#include "Mesh_Helpers.hh"

namespace Amanzi {

class VerboseObject;
namespace AmanziGeometry {
class Point;
class RegionLabeledSet;
class GeometricModel;
}

namespace AmanziMesh {

class MeshFramework;

//
// The framework class itself provides setters/getters/attributes, all
// topology, and coordinates.
//
class MeshFramework  {
 protected:
  MeshFramework(const Comm_ptr_type& comm,
                const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
                const Teuchos::RCP<Teuchos::ParameterList>& plist);

 public:
  virtual ~MeshFramework() = default;

  // ----------------------
  // Accessors and Mutators
  // ----------------------
  Comm_ptr_type getComm() const { return comm_; }
  void setComm(const Comm_ptr_type& comm) { comm_ = comm; }

  Teuchos::RCP<Teuchos::ParameterList> getParameterList() const { return plist_; }
  void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& plist) { plist_ = plist; }

  Teuchos::RCP<const VerboseObject> getVerboseObject() const { return vo_; }
  void setVerboseObject(const Teuchos::RCP<const VerboseObject>& vo) { vo_ = vo; }

  Teuchos::RCP<const AmanziGeometry::GeometricModel> getGeometricModel() const { return gm_; }
  void setGeometricModel(const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm) { gm_ = gm; }

  // space dimension describes the dimension of coordinates in space
  std::size_t getSpaceDimension() const { return space_dim_; }
  void setSpaceDimension(unsigned int dim) { space_dim_ = dim; }

  // manifold dimension describes the dimensionality of the corresponding R^n
  // manifold onto which this mesh can be projected.
  std::size_t getManifoldDimension() const { return manifold_dim_; }
  void setManifoldDimension(const unsigned int dim) { manifold_dim_ = dim; }

  // Some meshes are subsets of or derived from a parent mesh.
  // Usually this is null, but some meshes may provide it.
  virtual Teuchos::RCP<const MeshFramework> getParentMesh() const { return Teuchos::null; }

  // Some meshes have a corresponding mesh that is better for visualization.
  const MeshFramework& getVisMesh() const {
    if (vis_mesh_.get()) return *vis_mesh_;
    return *this;
  }
  void setVisMesh(const Teuchos::RCP<const MeshFramework>& vis_mesh) { vis_mesh_ = vis_mesh; }

  virtual void writeToExodusFile(const std::string& filename) const {}

  // Some meshes have edges
  //
  // DEVELOPER NOTE: frameworks that do not implement edges need not provide
  // any edge method -- defaults here all throw errors.
  virtual bool hasEdges() const { return false; }

  // DEVELOPER NOTE: frameworks that do not implement nodes DO need to provide
  // ALL node methods to have them throw errors.  The default here assumes
  // nodes exist.
  virtual bool hasNodes() const { return true; }

  // Some meshes may natively order in the ExodusII ordering
  virtual bool isOrdered() const { return false; }

  // Some meshes can be deformed.
  virtual bool isDeformable() const { return false; }

  // Some meshes are logical meshes and do not have coordinate info.
  virtual bool isLogical() const { return false; }

  // ----------------
  // Entity meta-data
  // ----------------
  virtual std::size_t getNumEntities(const Entity_kind kind, const Parallel_kind ptype) const = 0;

  // Parallel type of the entity.
  //
  // DEVELOPER NOTE: meshes which order entities by OWNED, GHOSTED need not
  // implement this method.
  virtual Parallel_kind getEntityPtype(const Entity_kind kind, const Entity_ID entid) const;

  // Global ID of any entity
  //
  // DEVELOPER NOTE: serial meshes need not provide this method -- the default
  // returns the LID.
  virtual Entity_GID getEntityGID(const Entity_kind kind, const Entity_ID lid) const;

  // DEVELOPER NOTE: serial meshes need not provide this method -- the default
  // returns a list of LIDs.
  virtual cEntity_GID_View getEntityGIDs(const Entity_kind kind, const Parallel_kind ptype) const;

  // corresponding entity in the parent mesh
  virtual Entity_ID getEntityParent(const Entity_kind kind, const Entity_ID entid) const;

  // Cell types: UNKNOWN, TRI, QUAD, etc. See MeshDefs.hh
  //
  // DEVELOPER NOTE: Default implementation guesses based on topology.
  virtual Cell_kind getCellType(const Entity_ID cellid) const;


  //---------------------
  // Geometry
  //---------------------
  // locations
  virtual AmanziGeometry::Point getNodeCoordinate(const Entity_ID node) const = 0;
  virtual void setNodeCoordinate(const Entity_ID nodeid, const AmanziGeometry::Point& ncoord);

  virtual cPoint_View getCellCoordinates(const Entity_ID c) const;
  virtual cPoint_View getFaceCoordinates(const Entity_ID f) const;
  virtual cPoint_View getEdgeCoordinates(const Entity_ID e) const;


  //---------------------
  // Downward adjacencies
  //---------------------
  // Get faces of a cell
  //
  // On a distributed mesh, this will return all the faces of the
  // cell, OWNED or GHOST. If the framework supports it, the faces will be
  // returned in a standard order according to Exodus II convention
  // for standard cells; in all other situations (not supported or
  // non-standard cells), the list of faces will be in arbitrary order
  //
  // EXTENSIONS: MSTK FRAMEWORK: by the way the parallel partitioning,
  // send-receive protocols and mesh query operators are designed, a side 
  // effect of this is that master and ghost entities will have the same
  // hierarchical topology.
  void getCellFaces(const Entity_ID c,
                    cEntity_ID_View& faces) const {
    getCellFacesAndDirs(c, faces, nullptr);
  }

  void getCellFaceDirs(const Entity_ID c,
                       cEntity_Direction_View& dirs) const {
    cEntity_ID_View faces;
    getCellFacesAndDirs(c, faces, &dirs);
  }

  // Get faces of a cell and directions in which the cell uses the face
  //
  // On a distributed mesh, this will return all the faces of the
  // cell, OWNED or GHOST. If the framework supports it, the faces will be
  // returned in a standard order according to Exodus II convention
  // for standard cells.
  //
  // In 3D, direction is 1 if face normal points out of cell
  // and -1 if face normal points into cell
  // In 2D, direction is 1 if face/edge is defined in the same
  // direction as the cell polygon, and -1 otherwise
  virtual void getCellFacesAndDirs(
    const Entity_ID c,
    cEntity_ID_View& faces,
    cEntity_Direction_View * const dirs) const = 0;

  virtual void getCellEdges(const Entity_ID c, cEntity_ID_View& edges) const;
  virtual void getCellNodes(const Entity_ID c, cEntity_ID_View& nodes) const;

  void getFaceEdges(const Entity_ID f,
                  cEntity_ID_View& edges) const {
    getFaceEdgesAndDirs(f, edges);
  }

  // Get edges of a face and directions in which the face uses the edges.
  //
  // In 3D, edge direction is 1 when it is oriented counter clockwise
  // with respect to the face natural normal.
  //
  // On a distributed mesh, this will return all the edges of the
  // face, OWNED or GHOST. If the framework supports it, the edges will be
  // returned in a ccw order around the face as it is naturally defined.
  //
  // IMPORTANT NOTE IN 2D: In meshes where the cells are two
  // dimensional, faces and edges are identical. For such cells, this
  // operator will return a single edge and a direction of 1. However,
  // this direction cannot be relied upon to compute, say, a contour
  // integral around the 2D cell.
  virtual void getFaceEdgesAndDirs(const Entity_ID f,
          cEntity_ID_View& edges,
          cEntity_Direction_View * const dirs=nullptr) const;

  // Get nodes of face
  //
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal.
  virtual void getFaceNodes(const Entity_ID f, cEntity_ID_View& nodes) const = 0;

  virtual void getEdgeNodes(const Entity_ID e, cEntity_ID_View& nodes) const;

  //-------------------
  // Upward adjacencies
  //-------------------
  // The cells are returned in no particular order. Also, the order of cells
  // is not guaranteed to be the same for corresponding faces on different
  // processors
  virtual void getFaceCells(const Entity_ID f,
                            const Parallel_kind ptype,
                            cEntity_ID_View& cells) const = 0;

  // Cells of a given Parallel_kind connected to an edge
  //
  // The order of cells is not guaranteed to be the same for corresponding
  // edges on different processors
  virtual void getEdgeCells(const Entity_ID edgeid,
                            const Parallel_kind ptype,
                            cEntity_ID_View& cellids) const;

  // Faces of type 'ptype' connected to an edge
  // NOTE: The order of faces is not guaranteed to be the same for
  // corresponding edges on different processors
  virtual void getEdgeFaces(const Entity_ID edgeid,
                            const Parallel_kind ptype,
                            cEntity_ID_View& faceids) const;

  // Cells of type 'ptype' connected to a node
  // NOTE: The order of cells is not guaranteed to be the same for
  // corresponding nodes on different processors
  virtual void getNodeCells(const Entity_ID nodeid,
                            const Parallel_kind ptype,
                            cEntity_ID_View& cellids) const;

  // Faces of type parallel 'ptype' connected to a node
  // NOTE: The order of faces is not guarnateed to be the same for
  // corresponding nodes on different processors
  virtual void getNodeFaces(const Entity_ID nodeid,
                            const Parallel_kind ptype,
                            cEntity_ID_View& faceids) const = 0;

  // Edges of type 'ptype' connected to a node
  //
  // The order of edges is not guaranteed to be the same for corresponding
  // node on different processors
  virtual void getNodeEdges(const Entity_ID nodeid,
                          const Parallel_kind ptype,
                          cEntity_ID_View& edgeids) const;

  //--------------------------------------------------------------
  // Mesh Sets for ICs, BCs, Material Properties and whatever else
  //--------------------------------------------------------------
  // Default implementation does not support framework sets.
  virtual bool isValidSetType(const AmanziGeometry::RegionType rtype, const Entity_kind kind) const {
    return false;
  }

  // Get list of entities of type kind from a framework set.
  //
  virtual void getSetEntities(const AmanziGeometry::RegionLabeledSet& region,
          const Entity_kind kind,
          const Parallel_kind ptype,
          cEntity_ID_View& entids) const;

  void hasEdgesOrThrow() const;

protected:
  void throwNotImplemented_(const std::string& fname) const;
  Cell_kind getCellType_(const Entity_ID c, const cEntity_ID_View& faces) const;

 protected:
  Comm_ptr_type comm_;
  Teuchos::RCP<Teuchos::ParameterList> plist_;
  Teuchos::RCP<const VerboseObject> vo_;
  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm_;
  Teuchos::RCP<const MeshFramework> vis_mesh_;

  std::size_t space_dim_;
  std::size_t manifold_dim_;
};




}  // namespace AmanziMesh
}  // namespace Amanzi


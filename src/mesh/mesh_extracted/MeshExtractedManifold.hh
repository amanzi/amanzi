/*
  Mesh Extracted

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov

  We assume that the parent mesh is a 3D mesh and the extracted mesh 
  lives on a 2D non-manifold. The extracted mesh has 3D geometry.

  Faces and edges are geometrycally identical in the extracted mesh
  and have same ids, but they parent ids are different. 
*/

#ifndef AMANZI_MESH_EXTRACTED_MANIFOLD_HH_
#define AMANZI_MESH_EXTRACTED_MANIFOLD_HH_

#include <memory>
#include <string>
#include <vector>

// TPLs
#include "Epetra_Map.h"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "AmanziComm.hh"
#include "dbc.hh"
#include "errors.hh"
#include "Mesh.hh"
#include "Region.hh"

namespace Amanzi {
namespace AmanziMesh {

class MeshExtractedManifold : public MeshFramework {
 public:
  // Construct a mesh by extracting a subset of entities from another
  // mesh. The subset may be specified by a list of entities. 
  MeshExtractedManifold(const Teuchos::RCP<Mesh>& parent_mesh,
                        const std::string& setname, 
                        const Entity_kind entity_kind,
                        const Comm_ptr_type& comm = Teuchos::null,
                        const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm = Teuchos::null,
                        const Teuchos::RCP<Teuchos::ParameterList>& plist = Teuchos::null,
                        bool request_faces = true,
                        bool request_edges = false,
                        bool flattened = false);
  ~MeshExtractedManifold() {};

  // initialization
  void InitParentMaps(const std::string& setname);
  void InitEpetraMaps();
  void InitExteriorEpetraMaps();

  //virtual Teuchos::RCP<const MeshFramework> getParentMesh() const override { return parent_mesh_; }

  // parallel type of entity - OWNED or GHOST
  virtual
  Parallel_type getEntityPtype(const Entity_kind kind, const Entity_ID id) const override {
    return (id < nents_owned_[kind]) ? Parallel_type::OWNED : Parallel_type::GHOST;
  }

  // parent entity if this mesh was extracted from another mesh
  virtual
  Entity_ID getEntityParent(const Entity_kind kind, const Entity_ID id) const override {
    return entid_to_parent_[kind][id];
  }

  // cell type - UNKNOWN, TRI, QUAD, ... See MeshDefs.hh
  virtual Cell_type getCellType(const Entity_ID c) const override;

  // general mesh information
  // -- number of entities of any kind (cell, face, node) and in a
  //    particular category (OWNED, GHOST, ALL)
  virtual std::size_t getNumEntities(const Entity_kind kind,
                                    const Parallel_type ptype) const override;

  // -- global ID of any entity
  virtual Entity_ID getEntityGID(const Entity_kind kind, const Entity_ID lid) const override {
    return ent_map_wghost_[kind]->GID(lid);
  }

  // downward adjacencies
  // -- cells
  virtual void getCellNodes(const Entity_ID c, Entity_ID_List& nodes) const override;

  // -- faces
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face are returned
  // In 3D, the nodes are returned in ccw order consistent with the face normal
  virtual void getFaceNodes(const Entity_ID f, Entity_ID_List& nodes) const override;

  // -- edges
  virtual void getEdgeNodes(const Entity_ID e, Entity_ID_List& nodes) const override {
    nodes.resize(2); 
    nodes[0] = e;
    nodes[1] = e; 
  }

  // upward adjacencies
  // -- cells of type 'ptype' connected to a node - The order of cells is not guaranteed
  //    to be the same for corresponding nodes on different processors
  virtual void getNodeCells(const Entity_ID n, const Parallel_type ptype,
                              Entity_ID_List& cells) const override;

  // -- faces of type 'ptype' connected to a node - The order of faces is not guaranteed 
  //    to be the same for corresponding nodes on different processors
  virtual void getNodeFaces(const Entity_ID n, const Parallel_type ptype,
                              Entity_ID_List& faces) const override {
    // parent_mesh_->node_get_edges() is not implemented, another algorithm is needed
    AMANZI_ASSERT(false);
  }

  // -- cells of type 'ptype' connected to an edge - The order of cells is not guaranteed
  //    to be the same for corresponding edges on different processors
  virtual void getEdgeCells(const Entity_ID e, const Parallel_type ptype,
                              Entity_ID_List& cells) const override;

  // same level adjacencies
  // -- face connected neighboring cells of given cell of a particular ptype
  //    (e.g. a hex has 6 face neighbors)
  // 
  // The order in which the cellids are returned cannot be guaranteed in general
  // except when ptype = ALL, in which case the cell ids will correspond to cells
  // across the respective faces given by cell_get_faces().
  virtual void getFaceCells(const Entity_ID c, const Parallel_type ptype,
                                       Entity_ID_List& cells) const override;

  // Mesh entity geometry
  // -- nodes
  virtual AmanziGeometry::Point getNodeCoordinate(const Entity_ID n) const override;

  // -- node modifications are blocked for the extracted mesh
  virtual void setNodeCoordinate(const Entity_ID n, const AmanziGeometry::Point& xyz) override {
    AMANZI_ASSERT(false);
  }
  //virtual void setNodeCoordinate(const Entity_ID n, const double *xyz) override {
  //  AMANZI_ASSERT(false);
  //}

  // -- edges

  // -- faces
  virtual Point_List getFaceCoordinates(const Entity_ID f) const override;

  // -- cells
  //    coordinates of nodes in cell in the standard order (Exodus II convention)
  //    STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES IN 3D
  //    For a general polyhedron this returns the node coordinates in arbitrary order
  virtual Point_List getCellCoordinates(const Entity_ID c) const override;

  // -- mesh modifications are blocked for the extracted mesh
  //virtual int deform(const std::vector<double>& target_cell_volumes,
  //                   const std::vector<double>& min_cell_volumes_in,
  //                   const Entity_ID_List& fixed_nodes,
  //                   const bool move_vertical) override {
  //  AMANZI_ASSERT(false);
  //  return 0;
  //}

  // Mesh Sets for ICs, BCs, Material Properties
  virtual
  bool isValidSetType(const AmanziGeometry::RegionType rtype, const Entity_kind kind) const override {
    return parent_mesh_->isValidSetType(rtype, kind);
  }

  // -- entities
  using MeshFramework::getSetEntities;

  decltype(auto) 
  getSetEntitiesAndVolumeFractions(const std::string& setname,
                                         const Entity_kind kind,
                                         const Parallel_type ptype) const;

  // miscellaneous functions
  //virtual void write_to_exodus_file(const std::string filename) const override {
  //  AMANZI_ASSERT(false);
  //}

  // low-level supporting functions
  // -- get faces of a cell and directions in which it is used - this function is
  //    implemented in each mesh framework. The results are cached in the base class
  virtual void getCellFacesAndDirs(const Entity_ID c,
                                                 Entity_ID_List& faces,
                                                 std::vector<int> *fdirs) const override;

  // -- cells connected to a face - this function is implemented in each mesh
  //    framework. The results are cached in the base class
  //virtual void getFaceCells(const Entity_ID f,
  //                                      const Parallel_type ptype,
  //                                      Entity_ID_List& cells) const override;

  // -- edges of a face - this function is implemented in each mesh
  //    framework. The results are cached in the base class
  virtual void getFaceEdgesAndDirs(const Entity_ID f,
                                                 Entity_ID_List& edges,
                                                 std::vector<int> *edirs) const override;

  // -- edges of a cell - this function is implemented in each mesh
  //    framework. The results are cached in the base class.
  virtual void getCellEdges(const Entity_ID c, Entity_ID_List& edges) const override;

 private:
  Entity_ID_List build_set_(const Teuchos::RCP<const AmanziGeometry::Region>& rgn,
                            const Entity_kind kind) const;

  Entity_ID_List build_set_cells_(const Teuchos::RCP<const AmanziGeometry::Region>& rgn, bool* missing) const;
  Entity_ID_List build_set_faces_(const Teuchos::RCP<const AmanziGeometry::Region>& rgn, bool* missing) const;
  Entity_ID_List build_set_nodes_(const Teuchos::RCP<const AmanziGeometry::Region>& rgn, bool* missing) const;

  Entity_ID_List build_from_parent_(const std::string& rgnname, const Entity_kind kind_d) const;

  void TryExtension_(const std::string& setname,
                     Entity_kind kind_p, Entity_kind kind_d, Entity_ID_List* setents) const;
  std::map<Entity_ID, int> EnforceOneLayerOfGhosts_(const std::string& setname, Entity_kind kind,
                                                    Entity_ID_List* setents) const;

  void PrintSets_() const;

 private: 
  Teuchos::RCP<const Mesh> parent_mesh_;

  // owned ids are enforced to be first in the child -> parent map
  mutable std::map<Entity_kind, Entity_ID> nents_owned_, nents_ghost_;
  mutable std::map<Entity_kind, Entity_ID_List> entid_to_parent_;
  mutable std::map<Entity_kind, std::map<Entity_ID, Entity_ID> > parent_to_entid_;  // reverse to previous map
  mutable std::map<Entity_kind, Teuchos::RCP<const Epetra_Map> > ent_map_owned_, ent_map_wghost_;
  mutable std::map<Entity_kind, Teuchos::RCP<const Epetra_Map> > ent_extmap_owned_, ent_extmap_wghost_;

  std::map<Set_ID, std::vector<int> > regions_;
  Teuchos::RCP<Epetra_Import> exterior_face_importer_;

  // sets
  mutable std::map<std::string, Entity_ID_List> sets_;
  mutable std::map<std::string, Entity_ID_List> parent_labeledsets_;

  // deformation
  mutable bool flattened_;
};

}  // namespace AmanziMesh
}  // namespace Amanzi

#endif


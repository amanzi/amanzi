/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  Mesh Extracted

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
  MeshExtractedManifold(
    const Teuchos::RCP<const Mesh>& parent_mesh,
    const std::string& setname,
    const Entity_kind entity_kind,
    const Comm_ptr_type& comm = Teuchos::null,
    const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm = Teuchos::null,
    const Teuchos::RCP<Teuchos::ParameterList>& plist = Teuchos::null,
    bool flattened = false);
  ~MeshExtractedManifold(){};

  // initialization
  void InitParentMaps(const std::string& setname);
  void InitEpetraMaps();

  virtual Teuchos::RCP<const MeshFramework> getParentMesh() const override
  {
    AMANZI_ASSERT(false);
    return Teuchos::null;
  }

  // Entity meta-data
  virtual Entity_ID getEntityGID(const Entity_kind kind, const Entity_ID lid) const override
  {
    return ent_map_wghost_[kind]->GID(lid);
  }

  // parent entity if this mesh was extracted from another mesh
  virtual Entity_ID getEntityParent(const Entity_kind kind, const Entity_ID id) const override
  {
    return entid_to_parent_[kind][id];
  }

  // general mesh information
  // -- number of entities of any kind (cell, face, node) and in a
  //    particular category (OWNED, GHOST, ALL)
  virtual std::size_t
  getNumEntities(const Entity_kind kind, const Parallel_kind ptype) const override;

  // -- faces
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face are returned
  // In 3D, the nodes are returned in ccw order consistent with the face normal
  virtual void getFaceNodes(const Entity_ID f, cEntity_ID_View& nodes) const override;

  // -- edges
  virtual void getEdgeNodes(const Entity_ID e, cEntity_ID_View& nodes) const override
  {
    Entity_ID_View lnodes("lnodes", 2);
    lnodes[0] = e;
    lnodes[1] = e;
    nodes = lnodes;
  }

  // -- faces of type 'ptype' connected to a node - The order of faces is not guaranteed
  //    to be the same for corresponding nodes on different processors
  virtual void getNodeFaces(const Entity_ID n, cEntity_ID_View& nfaces) const override
  {
    auto parent_nodes = entid_to_parent_.at(Entity_kind::NODE);
    auto parent_faces = entid_to_parent_.at(Entity_kind::FACE);
    auto my_parent_node = parent_nodes(n);
    auto my_parent_faces = parent_mesh_->getNodeFaces(my_parent_node);
    Entity_ID_View faces("faces", my_parent_faces.size());
    int i = 0;
    for (auto f : my_parent_faces) {
      for (auto pf : parent_faces)
        if (f == pf) faces(i++) = parent_to_entid_[Entity_kind::FACE].at(f);
    }
    Kokkos::resize(faces, i);
    nfaces = faces;
  }

  // -- cells of type 'ptype' connected to an edge - The order of cells is not guaranteed
  //    to be the same for corresponding edges on different processors
  virtual void getEdgeCells(const Entity_ID e, cEntity_ID_View& cells) const override;

  // same level adjacencies
  // -- face connected neighboring cells of given cell of a particular ptype
  //    (e.g. a hex has 6 face neighbors)
  //
  // The order in which the cellids are returned cannot be guaranteed in general
  // except when ptype = ALL, in which case the cell ids will correspond to cells
  // across the respective faces given by cell_get_faces().
  virtual void getFaceCells(const Entity_ID c, cEntity_ID_View& cells) const override;

  // Mesh entity geometry
  // -- nodes
  virtual AmanziGeometry::Point getNodeCoordinate(const Entity_ID n) const override;

  // Mesh Sets for ICs, BCs, Material Properties
  virtual bool
  isValidSetType(const AmanziGeometry::RegionType rtype, const Entity_kind kind) const override
  {
    return true;
  }

  // miscellaneous functions

  // low-level supporting functions
  // -- get faces of a cell and directions in which it is used - this function is
  //    implemented in each mesh framework. The results are cached in the base class
  virtual void getCellFacesAndDirs(const Entity_ID c,
                                   cEntity_ID_View& faces,
                                   cDirection_View* fdirs) const override;

  // -- edges of a face - this function is implemented in each mesh
  //    framework. The results are cached in the base class
  virtual void getFaceEdgesAndDirs(const Entity_ID f,
                                   cEntity_ID_View& edges,
                                   cDirection_View* edirs) const override;

  // -- edges of a cell - this function is implemented in each mesh
  //    framework. The results are cached in the base class.
  virtual void getCellEdges(const Entity_ID c, cEntity_ID_View& edges) const override;

 private:
  void TryExtension_(const std::string& setname,
                     Entity_kind kind_p,
                     Entity_kind kind_d,
                     Entity_ID_View* setents) const;
  template <class Entity_ID_View_Type>
  std::map<Entity_ID, int> EnforceOneLayerOfGhosts_(const std::string& setname,
                                                    Entity_kind kind,
                                                    Entity_ID_View_Type* setents) const;

 private:
  Teuchos::RCP<const Mesh> parent_mesh_;

  // owned ids are enforced to be first in the child -> parent map
  mutable std::map<Entity_kind, Entity_ID> nents_owned_, nents_ghost_;
  mutable std::map<Entity_kind, Entity_ID_View> entid_to_parent_;
  mutable std::map<Entity_kind, std::map<Entity_ID, Entity_ID>>
    parent_to_entid_; // reverse to previous map
  mutable std::map<Entity_kind, Teuchos::RCP<const Epetra_Map>> ent_map_wghost_;

  mutable bool flattened_;
};

} // namespace AmanziMesh
} // namespace Amanzi

#endif

/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  Mesh Extracted

*/

#include <set>
#include <utility>

// TPLs
#include "Epetra_IntVector.h"

// Amanzi
#include "dbc.hh"
#include "errors.hh"
#include "RegionLogical.hh"
#include "RegionPoint.hh"
#include "VerboseObject.hh"

// Amanzi::Mesh
#include "BlockMapUtils.hh"
#include "MeshExtractedManifold.hh"

#include "MeshExtractedDefs.hh"

namespace Amanzi {
namespace AmanziMesh {

/* ******************************************************************
* Light-weigthed constructor
****************************************************************** */
MeshExtractedManifold::MeshExtractedManifold(
    const Teuchos::RCP<const Mesh>& parent_mesh,
    const std::string& setname, 
    const Entity_kind entity_kind,
    const Comm_ptr_type& comm,
    const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
    const Teuchos::RCP<Teuchos::ParameterList>& plist,
    bool flattened)
  : MeshFramework(comm, gm, plist),
    parent_mesh_(parent_mesh),
    flattened_(flattened)
{  
  vo_ = Teuchos::rcp(new VerboseObject(comm_, "MeshExtractedManifold", *plist_));

  int d = parent_mesh_->getSpaceDimension();
  setSpaceDimension(d);
  setManifoldDimension(d - 1);
  if (flattened_) setSpaceDimension(d - 1);

  InitParentMaps(setname); 
}

/* ******************************************************************
* Number of OWNED, GHOST or ALL entities of different types
*
* Number of entities of any kind (cell, face, edge, node) and in a
* particular category (OWNED, GHOST, ALL)
****************************************************************** */
std::size_t MeshExtractedManifold::getNumEntities(const Entity_kind kind, 
                                                 const Parallel_kind ptype) const
{
  if (ptype == Parallel_kind::OWNED)
    return nents_owned_[kind];
  else if (ptype == Parallel_kind::ALL)
    return nents_owned_[kind] + nents_ghost_[kind];

  return nents_ghost_[kind]; 
}

/* ******************************************************************
* Connectivity list: cell -> faces. Base routine for getCellFaces().
****************************************************************** */
void MeshExtractedManifold::getCellFacesAndDirs(
    const Entity_ID c,
    Entity_ID_View& faces, Entity_Direction_View *fdirs) const
{
  int fp = entid_to_parent_[Entity_kind::CELL][c];
  parent_mesh_->getFaceEdgesAndDirs(fp, faces, fdirs);
  int nfaces = faces.size();
 
  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    faces[i] = parent_to_entid_[Entity_kind::FACE][f];
  }

  // algorithms on a non-manifold use multiple normals and special continuity
  // equations for fluxes, so that orientation does not play role.
  // Now if the mesh is flattened, the Mesh class algorithm uses 2D edge
  // orientation. In this case, the result is correct iff the 3D face normal
  // is exterior. 
  if (! flattened_ && fdirs) {
    for (int i = 0; i < nfaces; ++i) {
      (*fdirs)[i] = 1;
    }
  }
}


/* ******************************************************************
* Connectivity list: cell -> edges = cell -> faces
****************************************************************** */
void MeshExtractedManifold::getCellEdges(
    const Entity_ID c, Entity_ID_View& edges) const
{
  Entity_Direction_View edirs;

  int fp = entid_to_parent_[Entity_kind::CELL][c];
  parent_mesh_->getFaceEdgesAndDirs(fp, edges, &edirs);
  int nedges = edges.size();

  for (int i = 0; i < nedges; ++i) {
    int e = edges[i];
    edges[i] = parent_to_entid_[Entity_kind::FACE][e];
  }
}

/* ******************************************************************
* Connectivity list: face -> nodes
****************************************************************** */
void MeshExtractedManifold::getFaceNodes(const Entity_ID f, Entity_ID_View& nodes) const {
  Entity_ID_View v0;

  int ep = entid_to_parent_[Entity_kind::FACE][f];
  v0 = parent_mesh_->getEdgeNodes(ep);

  Kokkos::resize(nodes,2); 
  nodes[0] = parent_to_entid_[Entity_kind::NODE][v0[0]];
  nodes[1] = parent_to_entid_[Entity_kind::NODE][v0[1]];
}


/* ******************************************************************
* Connectivity list: face -> edges
****************************************************************** */
void MeshExtractedManifold::getFaceEdgesAndDirs(
    const Entity_ID f,
    Entity_ID_View& edges, Entity_Direction_View *edirs) const
{
  Entity_ID_View v0;

  int ep = entid_to_parent_[Entity_kind::FACE][f];
  v0 = parent_mesh_->getEdgeNodes(ep);

  Kokkos::resize(edges,2);
  edges[0] = parent_to_entid_[Entity_kind::NODE][v0[0]];
  edges[1] = parent_to_entid_[Entity_kind::NODE][v0[1]];
  Kokkos::resize(*edirs,2);
  (*edirs)[0] = 1; 
  (*edirs)[1] = 1; 
}

/* ******************************************************************
* Connectivity list: edge -> cells
****************************************************************** */
void MeshExtractedManifold::getEdgeCells(
   const Entity_ID e, const Parallel_kind ptype, Entity_ID_View& cells) const 
{
  Entity_ID_View faces;

  int ep = entid_to_parent_[Entity_kind::FACE][e];
  parent_mesh_->getEdgeFaces(ep, ptype, faces);
  int nfaces = faces.size();

  Kokkos::resize(cells, nfaces); 
  int cells_ct = 0; 
  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    auto it = parent_to_entid_[Entity_kind::CELL].find(f);
    if (it != parent_to_entid_[Entity_kind::CELL].end()) cells[cells_ct++] = it->second;
  }
  Kokkos::resize(cells,cells_ct); 
}


/* ******************************************************************
* Connectivity list: face -> cells
****************************************************************** */
void MeshExtractedManifold::getFaceCells(
    const Entity_ID f,
    const Parallel_kind ptype, Entity_ID_View& cells) const
{
  Entity_ID_View faces;

  int ep = entid_to_parent_[Entity_kind::FACE][f];
  parent_mesh_->getEdgeFaces(ep, ptype, faces);
  int nfaces = faces.size();

  Kokkos::resize(cells, nfaces); 
  int cells_ct = 0; 
  for (int i = 0; i < nfaces; ++i) {
    auto it = parent_to_entid_[Entity_kind::CELL].find(faces[i]);
    if (it != parent_to_entid_[Entity_kind::CELL].end()) cells[cells_ct++] = it->second;
  }
  Kokkos::resize(cells,cells_ct); 
}

/* ******************************************************************
* Position vector for a node
****************************************************************** */
AmanziGeometry::Point MeshExtractedManifold::getNodeCoordinate(
    const Entity_ID n) const
{
  AmanziGeometry::Point xyz; 
  auto np = entid_to_parent_[Entity_kind::NODE][n];
  xyz = parent_mesh_->getNodeCoordinate(np);

  if (flattened_) xyz.set(xyz[0], xyz[1]);
  return xyz; 
}

/* ******************************************************************
* Create internal maps for child->parent
****************************************************************** */
void MeshExtractedManifold::InitParentMaps(const std::string& setname)
{
  nents_owned_.clear();
  nents_ghost_.clear();
  entid_to_parent_.clear();
  parent_to_entid_.clear();

  std::vector<Entity_kind> kinds_extracted({Entity_kind::CELL, Entity_kind::FACE, Entity_kind::NODE});
  std::vector<Entity_kind> kinds_parent({Entity_kind::FACE, Entity_kind::EDGE, Entity_kind::NODE});

  for (int i = 0; i < 3; ++i) {
    auto kind_d = kinds_extracted[i];
    auto kind_p = kinds_parent[i];

    // build edge set from Exodus labeled face set
    Entity_ID_View setents;
    Double_View vofs;

    TryExtension_(setname, kind_p, kind_d, &setents);
    if (setents.size() == 0)
      std::tie(setents,vofs) = parent_mesh_->getSetEntitiesAndVolumeFractions(setname, kind_p, Parallel_kind::ALL);

    auto marked_ents = EnforceOneLayerOfGhosts_(setname, kind_p, &setents);

    // extract owned ids
    auto& ids_p = entid_to_parent_[kind_d];
    int ids_p_s = std::distance(marked_ents.begin(), marked_ents.end()); 
    Kokkos::resize(ids_p, ids_p_s); 
    ids_p_s = 0; 
    for (auto it = marked_ents.begin(); it != marked_ents.end(); ++it) {
      if (it->second == MASTER) ids_p[ids_p_s++] = it->first;
    }

    nents_owned_[kind_d] = ids_p_s;
    nents_ghost_[kind_d] = marked_ents.size() - ids_p_s;

    // extract ghost ids
    for (auto it = marked_ents.begin(); it != marked_ents.end(); ++it) {
      if (it->second == GHOST) ids_p[ids_p_s++] = it->first;
    }

    Kokkos::resize(ids_p, ids_p_s); 
    // create reverse ordered map
    auto& ids_d = parent_to_entid_[kind_d];
    ids_d.clear();
    for (int n = 0; n < ids_p.size(); ++n) {
      ids_d[ids_p[n]] = n;
    }
  }
}

/* ******************************************************************
* Exception due to limitations of the base mesh framework.
****************************************************************** */
void MeshExtractedManifold::TryExtension_(
    const std::string& setname,
    Entity_kind kind_p, Entity_kind kind_d, Entity_ID_View* setents) const
{
  // labeled set: extract edges

  const auto& gm = getGeometricModel();
  if (gm == Teuchos::null) return;

  auto rgn = gm->FindRegion(setname);
  if (rgn->get_type() != AmanziGeometry::RegionType::LABELEDSET) return;

  // populate list of edges
  auto [faceents,vofs] = parent_mesh_->getSetEntitiesAndVolumeFractions(setname, Entity_kind::FACE, Parallel_kind::ALL);
  auto marked_ents = EnforceOneLayerOfGhosts_(setname, Entity_kind::FACE, &faceents);

  Entity_ID_View edges, nodes;
  Entity_Direction_View dirs;
  std::set<Entity_ID> setents_tmp;

  for (auto it = marked_ents.begin(); it != marked_ents.end(); ++it) {
    int f = it->first;
    if (kind_p == Entity_kind::FACE) {
      setents_tmp.insert(f);
    }
    else if (kind_p == Entity_kind::EDGE) {
      parent_mesh_->getFaceEdgesAndDirs(f, edges, &dirs);
      for (int i = 0; i < edges.size(); ++i) {
        setents_tmp.insert(edges[i]);
      }
    }
    else if (kind_p == Entity_kind::NODE) {
      parent_mesh_->getFaceNodes(f, nodes);
      for (int i = 0; i < nodes.size(); ++i) {
        setents_tmp.insert(nodes[i]);
      }
    }
  }

  setToView(*setents, setents_tmp); 

}


/* ******************************************************************
* Limits the set of parent objects to only one layer of ghosts.
****************************************************************** */
std::map<Entity_ID, int> MeshExtractedManifold::EnforceOneLayerOfGhosts_(
    const std::string& setname, Entity_kind kind, Entity_ID_View* setents) const
{
  // base set is the set of master cells
  Entity_ID_View fullset;
  if (kind != Entity_kind::FACE) {
    Double_View vofs;
    std::tie(fullset,vofs) = parent_mesh_->getSetEntitiesAndVolumeFractions(setname, Entity_kind::FACE, Parallel_kind::ALL);
  } else { 
    fullset = *setents;
  }

  // initial set of entities is defined by master parent faces and is marked as 
  // potential master entities
  Entity_ID_View nodes, edges, faces;
  Entity_Direction_View dirs;
  int nfaces_owned = parent_mesh_->getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
  std::map<Entity_ID, int> nodeset0, nodeset, edgeset, faceset;

  for (int n = 0; n < fullset.size(); ++n) {
    int f = fullset[n];
    if (f < nfaces_owned) {
      parent_mesh_->getFaceNodes(f, nodes);
      for (int i = 0; i < nodes.size(); ++i) nodeset0[nodes[i]] = MASTER;
 
      if (kind == Entity_kind::EDGE) {
        parent_mesh_->getFaceEdgesAndDirs(f, edges, &dirs);
        for (int i = 0; i < edges.size(); ++i) edgeset[edges[i]] = MASTER;
      }
 
      faceset[f] = MASTER;
    } 
  }

  // ghosts entities are defined by neighboor faces of master faces. New
  // entities are marked as ghosts. Old entities are marked as undefined.
  nodeset = nodeset0;
  for (int n = 0; n < fullset.size(); ++n) {
    int f = fullset[n];
    if (f >= nfaces_owned) {
      bool found(false);
      parent_mesh_->getFaceNodes(f, nodes);
      for (int i = 0; i < nodes.size(); ++i) {
        if (nodeset0.find(nodes[i]) != nodeset0.end()) {
          found = true;
          break;
        }
      }
      if (found) {
        for (int i = 0; i < nodes.size(); ++i) {
          auto it = nodeset.find(nodes[i]);
          if (it == nodeset.end())
            nodeset[nodes[i]] = GHOST;
          else
            it->second |= GHOST;
        }

        if (kind == Entity_kind::EDGE) {
          parent_mesh_->getFaceEdgesAndDirs(f, edges, &dirs);
          for (int i = 0; i < edges.size(); ++i) {
            auto it = edgeset.find(edges[i]);
            if (it == edgeset.end())
              edgeset[edges[i]] = GHOST;
            else
              it->second |= GHOST;
          }
        }

        faceset[f] = GHOST;
      }
    }
  }

  // resolve master+ghost entities
  std::set<Entity_ID> auxset;
  for (int n = 0; n < fullset.size(); ++n) auxset.insert(fullset[n]);

  if (kind == Entity_kind::FACE) {
    return faceset;
  } else if (kind == Entity_kind::EDGE) {
    const auto& fmap = parent_mesh_->getMap(Entity_kind::FACE, true);

    int nowned = parent_mesh_->getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
    int gidmax = fmap.MaxAllGID();

    for (auto it = edgeset.begin(); it != edgeset.end(); ++it) {
      if (it->second == MASTER + GHOST) {
        parent_mesh_->getEdgeFaces(it->first, Parallel_kind::ALL, faces);
        int nfaces = faces.size();

        // compare maximum global ids for owned and all faces living on manifold
        Entity_ID gid_owned_min(gidmax + 1), gid_wghost_min(gidmax + 1);
        for (int n = 0; n < nfaces; ++n) {
          Entity_ID f = faces[n];
          if (auxset.find(f) != auxset.end()) {
            gid_wghost_min = std::min(gid_wghost_min, fmap.GID(f));
            if (f < nowned)
              gid_owned_min = std::min(gid_owned_min, fmap.GID(f));
          }
        }
        it->second = (gid_wghost_min == gid_owned_min) ? MASTER : GHOST;
      }
    }
    return edgeset;
  } else if (kind == Entity_kind::NODE) {
    const auto& fmap = parent_mesh_->getMap(Entity_kind::FACE, true);

    int nowned = parent_mesh_->getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
    int gidmax = fmap.MaxAllGID();

    for (auto it = nodeset.begin(); it != nodeset.end(); ++it) {
      if (it->second == MASTER + GHOST) {
        parent_mesh_->getNodeFaces(it->first, Parallel_kind::ALL, faces);
        int nfaces = faces.size();

        // compare maximum global ids for owned and all faces living on manifold
        Entity_ID gid_owned_min(gidmax + 1), gid_wghost_min(gidmax + 1);
        for (int n = 0; n < nfaces; ++n) {
          Entity_ID f = faces[n];
          if (auxset.find(f) != auxset.end()) {
            gid_wghost_min = std::min(gid_wghost_min, fmap.GID(f));
            if (f < nowned)
              gid_owned_min = std::min(gid_owned_min, fmap.GID(f));
          }
        }
        it->second = (gid_wghost_min == gid_owned_min) ? MASTER : GHOST;
      }
    }
    return nodeset;
  }
  return std::map<Entity_ID, int>();
}


/* ******************************************************************
* Global ID of any entity
****************************************************************** */
Entity_GID MeshExtractedManifold::getEntityGID(
    const Entity_kind kind, const Entity_ID lid) const
{
  Entity_kind kind_p;

  switch (kind) {
  case NODE:
    kind_p = Entity_kind::NODE;
    break;
  case EDGE:
    kind_p = Entity_kind::EDGE;
    break;
  case FACE:
    kind_p = Entity_kind::EDGE;
    break;
  case CELL:
    kind_p = Entity_kind::FACE;
    break;
  default:
    std::cerr << "Global ID requested for unknown entity type" << std::endl;
  }

  int id = entid_to_parent_[kind][lid];
  return parent_mesh_->getMap(kind_p, true).GID(id);
}


}  // namespace AmanziMesh
}  // namespace Amanzi


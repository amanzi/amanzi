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
  : MeshFramework(comm, gm, plist), parent_mesh_(parent_mesh), flattened_(flattened)
{
  vo_ = Teuchos::rcp(new VerboseObject(comm_, "MeshExtractedManifold", *plist_));

  int d = parent_mesh_->getSpaceDimension();
  setSpaceDimension(d);
  setManifoldDimension(d - 1);
  if (flattened_) setSpaceDimension(d - 1);

  InitParentMaps(setname);
  InitEpetraMaps();
}


/* ******************************************************************
* Epetra maps are structures specifying the global IDs of entities
* owned or used by this processor. This helps Epetra understand
* inter-partition dependencies of the data.
****************************************************************** */
void
MeshExtractedManifold::InitEpetraMaps()
{
  std::vector<Entity_kind> kinds_extracted({ CELL, FACE, NODE });
  std::vector<Entity_kind> kinds_parent({ FACE, EDGE, NODE });

  for (int i = 0; i < 3; ++i) {
    auto kind_d = kinds_extracted[i];
    auto kind_p = kinds_parent[i];

    // compute (discontinuous) owned global ids using the parent map
    Teuchos::RCP<const Epetra_BlockMap> parent_map =
      Teuchos::rcpFromRef(parent_mesh_->getMap(kind_p, false));
    Teuchos::RCP<const Epetra_BlockMap> parent_map_wghost =
      Teuchos::rcpFromRef(parent_mesh_->getMap(kind_p, true));

    int nents = nents_owned_[kind_d];
    int nents_wghost = nents_owned_[kind_d] + nents_ghost_[kind_d];
    auto gids = new int[nents_wghost];

    for (int n = 0; n < nents; ++n) {
      int id = entid_to_parent_[kind_d][n];
      gids[n] = parent_map_wghost->GID(id);
    }

    auto subset_map = Teuchos::rcp(new Epetra_Map(-1, nents, gids, 0, *comm_));

    // compute owned + ghost ids using the parent map and the minimum global id
    for (int n = 0; n < nents_wghost; ++n) {
      int id = entid_to_parent_[kind_d][n];
      gids[n] = parent_map_wghost->GID(id);
    }

    auto subset_map_wghost = Teuchos::rcp(new Epetra_Map(-1, nents_wghost, gids, 0, *comm_));
    delete[] gids;

    // create continuous maps
    auto mymesh = Teuchos::rcpFromRef(*this);
    auto tmp = createContiguousMaps(mymesh,
                                    std::make_pair(parent_map, parent_map_wghost),
                                    std::make_pair(subset_map, subset_map_wghost));

    ent_map_wghost_[kind_d] = tmp.second;
  }
}

/* ******************************************************************
* Number of OWNED, GHOST or ALL entities of different types
*
* Number of entities of any kind (cell, face, edge, node) and in a
* particular category (OWNED, GHOST, ALL)
****************************************************************** */
std::size_t
MeshExtractedManifold::getNumEntities(const Entity_kind kind, const Parallel_kind ptype) const
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
void
MeshExtractedManifold::getCellFacesAndDirs(const Entity_ID c,
                                           cEntity_ID_View& faces,
                                           cEntity_Direction_View* fdirs) const
{
  Entity_ID_View lfaces;
  int fp = entid_to_parent_[Entity_kind::CELL][c];
  parent_mesh_->getFaceEdgesAndDirs(fp, faces, fdirs);
  lfaces.fromConst(faces);
  int nfaces = lfaces.size();

  for (int i = 0; i < nfaces; ++i) {
    int f = lfaces[i];
    lfaces[i] = parent_to_entid_[Entity_kind::FACE][f];
  }

  // algorithms on a non-manifold use multiple normals and special continuity
  // equations for fluxes, so that orientation does not play role.
  // Now if the mesh is flattened, the Mesh class algorithm uses 2D edge
  // orientation. In this case, the result is correct iff the 3D face normal
  // is exterior.
  if (!flattened_ && fdirs) {
    Entity_Direction_View lfdirs;
    lfdirs.fromConst(*fdirs);
    for (int i = 0; i < nfaces; ++i) { (lfdirs)[i] = 1; }
    *fdirs = lfdirs;
  }
  faces = lfaces;
}


/* ******************************************************************
* Connectivity list: cell -> edges = cell -> faces
****************************************************************** */
void
MeshExtractedManifold::getCellEdges(const Entity_ID c, cEntity_ID_View& edges) const
{
  Entity_ID_View ledges;
  cEntity_Direction_View edirs;

  int fp = entid_to_parent_[Entity_kind::CELL][c];
  parent_mesh_->getFaceEdgesAndDirs(fp, edges, &edirs);
  ledges.fromConst(edges);
  int nedges = edges.size();

  for (int i = 0; i < nedges; ++i) {
    int e = ledges[i];
    ledges[i] = parent_to_entid_[Entity_kind::FACE][e];
  }
  edges = ledges;
}

/* ******************************************************************
* Connectivity list: face -> nodes
****************************************************************** */
void
MeshExtractedManifold::getFaceNodes(const Entity_ID f, cEntity_ID_View& nodes) const
{
  cEntity_ID_View v0;

  int ep = entid_to_parent_[Entity_kind::FACE][f];
  v0 = parent_mesh_->getEdgeNodes(ep);

  Entity_ID_View lnodes("lnodes", 2);
  lnodes[0] = parent_to_entid_[Entity_kind::NODE][v0[0]];
  lnodes[1] = parent_to_entid_[Entity_kind::NODE][v0[1]];
  nodes = lnodes;
}


/* ******************************************************************
* Connectivity list: face -> edges
****************************************************************** */
void
MeshExtractedManifold::getFaceEdgesAndDirs(const Entity_ID f,
                                           cEntity_ID_View& edges,
                                           cEntity_Direction_View* edirs) const
{
  cEntity_ID_View v0;

  int ep = entid_to_parent_[Entity_kind::FACE][f];
  v0 = parent_mesh_->getEdgeNodes(ep);

  Entity_ID_View ledges("ledges", 2);
  ledges[0] = parent_to_entid_[Entity_kind::NODE][v0[0]];
  ledges[1] = parent_to_entid_[Entity_kind::NODE][v0[1]];
  edges = ledges;
  Entity_ID_View ledirs("ledirs", 2);
  (ledirs)[0] = 1;
  (ledirs)[1] = 1;
  *edirs = ledirs;
}

/* ******************************************************************
* Connectivity list: edge -> cells
****************************************************************** */
void
MeshExtractedManifold::getEdgeCells(const Entity_ID e,
                                    const Parallel_kind ptype,
                                    cEntity_ID_View& cells) const
{
  cEntity_ID_View faces;

  int ep = entid_to_parent_[Entity_kind::FACE][e];
  parent_mesh_->getEdgeFaces(ep, ptype, faces);
  int nfaces = faces.size();

  Entity_ID_View lcells("lcells", nfaces);
  int cells_ct = 0;
  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    auto it = parent_to_entid_[Entity_kind::CELL].find(f);
    if (it != parent_to_entid_[Entity_kind::CELL].end()) lcells[cells_ct++] = it->second;
  }
  Kokkos::resize(lcells, cells_ct);
  cells = lcells;
}


/* ******************************************************************
* Connectivity list: face -> cells
****************************************************************** */
void
MeshExtractedManifold::getFaceCells(const Entity_ID f,
                                    const Parallel_kind ptype,
                                    cEntity_ID_View& cells) const
{
  cEntity_ID_View faces;

  int ep = entid_to_parent_[Entity_kind::FACE][f];
  parent_mesh_->getEdgeFaces(ep, ptype, faces);
  int nfaces = faces.size();

  Entity_ID_View lcells("lcells", nfaces);
  int cells_ct = 0;
  for (int i = 0; i < nfaces; ++i) {
    auto it = parent_to_entid_[Entity_kind::CELL].find(faces[i]);
    if (it != parent_to_entid_[Entity_kind::CELL].end()) lcells[cells_ct++] = it->second;
  }
  Kokkos::resize(lcells, cells_ct);
  cells = lcells;
}

/* ******************************************************************
* Position vector for a node
****************************************************************** */
AmanziGeometry::Point
MeshExtractedManifold::getNodeCoordinate(const Entity_ID n) const
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
void
MeshExtractedManifold::InitParentMaps(const std::string& setname)
{
  nents_owned_.clear();
  nents_ghost_.clear();
  entid_to_parent_.clear();
  parent_to_entid_.clear();

  std::vector<Entity_kind> kinds_extracted(
    { Entity_kind::CELL, Entity_kind::FACE, Entity_kind::NODE });
  std::vector<Entity_kind> kinds_parent(
    { Entity_kind::FACE, Entity_kind::EDGE, Entity_kind::NODE });

  for (int i = 0; i < 3; ++i) {
    auto kind_d = kinds_extracted[i];
    auto kind_p = kinds_parent[i];

    // build edge set from Exodus labeled face set
    Entity_ID_View setents;

    TryExtension_(setname, kind_p, kind_d, &setents);
    if (setents.size() == 0)
      std::tie(setents, std::ignore) =
        parent_mesh_->getSetEntitiesAndVolumeFractions(setname, kind_p, Parallel_kind::ALL);

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
    for (int n = 0; n < ids_p.size(); ++n) { ids_d[ids_p[n]] = n; }
  }
}

/* ******************************************************************
* Exception due to limitations of the base mesh framework.
****************************************************************** */
void
MeshExtractedManifold::TryExtension_(const std::string& setname,
                                     Entity_kind kind_p,
                                     Entity_kind kind_d,
                                     Entity_ID_View* setents) const
{
  // labeled set: extract edges

  const auto& gm = getGeometricModel();
  if (gm == Teuchos::null) return;

  auto rgn = gm->FindRegion(setname);
  if (rgn->get_type() != AmanziGeometry::RegionType::LABELEDSET) return;

  // populate list of edges
  auto [faceents, vofs] =
    parent_mesh_->getSetEntitiesAndVolumeFractions(setname, Entity_kind::FACE, Parallel_kind::ALL);
  auto marked_ents = EnforceOneLayerOfGhosts_(setname, Entity_kind::FACE, &faceents);

  cEntity_ID_View edges, nodes;
  cEntity_Direction_View dirs;
  std::set<Entity_ID> setents_tmp;

  for (auto it = marked_ents.begin(); it != marked_ents.end(); ++it) {
    int f = it->first;
    if (kind_p == Entity_kind::FACE) {
      setents_tmp.insert(f);
    } else if (kind_p == Entity_kind::EDGE) {
      parent_mesh_->getFaceEdgesAndDirs(f, edges, &dirs);
      for (int i = 0; i < edges.size(); ++i) { setents_tmp.insert(edges[i]); }
    } else if (kind_p == Entity_kind::NODE) {
      parent_mesh_->getFaceNodes(f, nodes);
      for (int i = 0; i < nodes.size(); ++i) { setents_tmp.insert(nodes[i]); }
    }
  }

  setToView(*setents, setents_tmp);
}


/* ******************************************************************
* Limits the set of parent objects to only one layer of ghosts.
****************************************************************** */
template <class Entity_ID_View_Type>
std::map<Entity_ID, int>
MeshExtractedManifold::EnforceOneLayerOfGhosts_(const std::string& setname,
                                                Entity_kind kind,
                                                Entity_ID_View_Type* setents) const
{
  // base set is the set of master cells
  Entity_ID_View fullset;
  if (kind != Entity_kind::FACE) {
    Double_View vofs;
    std::tie(fullset, vofs) = parent_mesh_->getSetEntitiesAndVolumeFractions(
      setname, Entity_kind::FACE, Parallel_kind::ALL);
  } else {
    fullset = *setents;
  }

  // initial set of entities is defined by master parent faces and is marked as
  // potential master entities
  cEntity_ID_View nodes, edges, faces;
  cEntity_Direction_View dirs;
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
            if (f < nowned) gid_owned_min = std::min(gid_owned_min, fmap.GID(f));
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
            if (f < nowned) gid_owned_min = std::min(gid_owned_min, fmap.GID(f));
          }
        }
        it->second = (gid_wghost_min == gid_owned_min) ? MASTER : GHOST;
      }
    }
    return nodeset;
  }
  return std::map<Entity_ID, int>();
}

} // namespace AmanziMesh
} // namespace Amanzi

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
    bool request_faces, bool request_edges,
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
  InitEpetraMaps(); 
  try {
    InitExteriorEpetraMaps();
  } catch(...) {
    if (vo_.get() && vo_->os_OK(Teuchos::VERB_HIGH)) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *(vo_->os()) << "Parent framework does not support buildng of exterior maps.\n";
    }
  }

  PrintSets_();
}


/* ******************************************************************
* Get cell type
****************************************************************** */
Cell_type MeshExtractedManifold::getCellType(const Entity_ID c) const
{
  Entity_ID_List faces; 
  getCellFaces(c,faces);
  int nfaces = faces.size();
  return (nfaces == 4) ? Amanzi::AmanziMesh::Cell_type::QUAD : 
    ((nfaces == 3) ? Amanzi::AmanziMesh::Cell_type::TRI : Amanzi::AmanziMesh::Cell_type::POLYGON);
}


/* ******************************************************************
* Number of OWNED, GHOST or ALL entities of different types
*
* Number of entities of any kind (cell, face, edge, node) and in a
* particular category (OWNED, GHOST, ALL)
****************************************************************** */
std::size_t MeshExtractedManifold::getNumEntities(const Entity_kind kind, 
                                                 const Parallel_type ptype) const
{
  if (ptype == Parallel_type::OWNED)
    return nents_owned_[kind];
  else if (ptype == Parallel_type::ALL)
    return nents_owned_[kind] + nents_ghost_[kind];

  return nents_ghost_[kind]; 
}

/* ******************************************************************
* Connectivity list: cell -> faces
****************************************************************** */
void MeshExtractedManifold::getCellFacesAndDirs(
    const Entity_ID c,
    Entity_ID_List& faces, std::vector<int> *fdirs) const
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
  if (! flattened_) {
    for (int i = 0; i < nfaces; ++i) {
      (*fdirs)[i] = 1;
    }
  }
}


/* ******************************************************************
* Connectivity list: cell -> edges = cell -> faces
****************************************************************** */
void MeshExtractedManifold::getCellEdges(
    const Entity_ID c, Entity_ID_List& edges) const
{
  std::vector<int> edirs;

  int fp = entid_to_parent_[Entity_kind::CELL][c];
  parent_mesh_->getFaceEdgesAndDirs(fp, edges, &edirs);
  int nedges = edges.size();

  for (int i = 0; i < nedges; ++i) {
    int e = edges[i];
    edges[i] = parent_to_entid_[Entity_kind::FACE][e];
  }
}


/* ******************************************************************
* Connectivity list: cell -> nodes
****************************************************************** */
void MeshExtractedManifold::getCellNodes(const Entity_ID c, Entity_ID_List& nodes) const
{
  int fp = entid_to_parent_[Entity_kind::CELL][c];
  parent_mesh_->getFaceNodes(fp, nodes);
  int nnodes = nodes.size();

  for (int i = 0; i < nnodes; ++i) {
    int v = nodes[i];
    nodes[i] = parent_to_entid_[Entity_kind::NODE][v];
  }
}


/* ******************************************************************
* Connectivity list: face -> nodes
****************************************************************** */
void MeshExtractedManifold::getFaceNodes(const Entity_ID f, Entity_ID_List& nodes) const {
  Entity_ID_List v0;

  int ep = entid_to_parent_[Entity_kind::FACE][f];
  v0 = parent_mesh_->getEdgeNodes(ep);

  nodes.clear();
  nodes.push_back(parent_to_entid_[Entity_kind::NODE][v0[0]]);
  nodes.push_back(parent_to_entid_[Entity_kind::NODE][v0[1]]);
}


/* ******************************************************************
* Connectivity list: face -> edges
****************************************************************** */
void MeshExtractedManifold::getFaceEdgesAndDirs(
    const Entity_ID f,
    Entity_ID_List& edges, std::vector<int> *edirs) const
{
  Entity_ID_List v0;

  int ep = entid_to_parent_[Entity_kind::FACE][f];
  v0 = parent_mesh_->getEdgeNodes(ep);

  edges.clear();
  edges.push_back(parent_to_entid_[Entity_kind::NODE][v0[0]]);
  edges.push_back(parent_to_entid_[Entity_kind::NODE][v0[1]]);
  edirs->resize(2, 1);
}


/* ******************************************************************
* Connectivity list: node -> cells
****************************************************************** */
void MeshExtractedManifold::getNodeCells(
    const Entity_ID n,
    const Parallel_type ptype, Entity_ID_List& cells) const 
{
  Entity_ID_List faces;

  int np = entid_to_parent_[Entity_kind::NODE][n];
  parent_mesh_->getNodeFaces(np, ptype, faces);
  int nfaces = faces.size();

  cells.clear();
  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    auto it = parent_to_entid_[Entity_kind::CELL].find(f);
    if (it != parent_to_entid_[Entity_kind::CELL].end()) cells.push_back(it->second);
  }
}


/* ******************************************************************
* Connectivity list: edge -> cells
****************************************************************** */
void MeshExtractedManifold::getEdgeCells(
   const Entity_ID e, const Parallel_type ptype, Entity_ID_List& cells) const 
{
  Entity_ID_List faces;

  int ep = entid_to_parent_[Entity_kind::FACE][e];
  parent_mesh_->getEdgeFaces(ep, ptype, faces);
  int nfaces = faces.size();

  cells.clear();
  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    auto it = parent_to_entid_[Entity_kind::CELL].find(f);
    if (it != parent_to_entid_[Entity_kind::CELL].end()) cells.push_back(it->second);
  }
}


/* ******************************************************************
* Connectivity list: face -> cells
****************************************************************** */
void MeshExtractedManifold::getFaceCells(
    const Entity_ID f,
    const Parallel_type ptype, Entity_ID_List& cells) const
{
  Entity_ID_List faces;

  int ep = entid_to_parent_[Entity_kind::FACE][f];
  parent_mesh_->getEdgeFaces(ep, ptype, faces);
  int nfaces = faces.size();

  cells.clear();
  for (int i = 0; i < nfaces; ++i) {
    auto it = parent_to_entid_[Entity_kind::CELL].find(faces[i]);
    if (it != parent_to_entid_[Entity_kind::CELL].end()) cells.push_back(it->second);
  }
}


/* ******************************************************************
* The node coordinates are returned in in arbitrary order
****************************************************************** */
Point_List MeshExtractedManifold::getCellCoordinates(
    const Entity_ID c) const
{
  Entity_ID_List nodes; 

  Point_List vxyz; 
  int fp = entid_to_parent_[Entity_kind::CELL][c];
  parent_mesh_->getFaceNodes(fp, nodes);
  int nnodes = nodes.size();

  AmanziGeometry::Point p(getSpaceDimension());

  vxyz.clear();
  for (int i = 0; i < nnodes; ++i) {
    p = parent_mesh_->getNodeCoordinate(nodes[i]);
    vxyz.push_back(p);
  }

  if (flattened_) {
    for (int i = 0; i < nnodes; ++i) {
      vxyz[i].set(vxyz[i][0], vxyz[i][1]);
    }
  }
  return vxyz; 
}


/* ******************************************************************
* Face coordinates use convention same as in face_to_nodes()
****************************************************************** */
Point_List MeshExtractedManifold::getFaceCoordinates(
    const Entity_ID f) const
{
  Point_List vxyz; 
  Entity_ID_List v0;

  int ep = entid_to_parent_[Entity_kind::FACE][f];
  v0 = parent_mesh_->getEdgeNodes(ep);

  AmanziGeometry::Point xyz(getSpaceDimension());

  vxyz.clear();
  xyz = parent_mesh_->getNodeCoordinate(v0[0]);
  vxyz.push_back(xyz);

  xyz = parent_mesh_->getNodeCoordinate(v0[1]);
  vxyz.push_back(xyz);

  if (flattened_) {
    for (int i = 0; i < 2; ++i) {
      vxyz[i].set(vxyz[i][0], vxyz[i][1]);
    }
  }
  return vxyz; 
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
* Get list of entities of type 'ptype' in set specified by setname
****************************************************************** */
decltype(auto) MeshExtractedManifold::getSetEntitiesAndVolumeFractions(
    const std::string& setname, 
    const Entity_kind kind, const Parallel_type ptype ) const
{

  std::vector<Entity_ID> setents; 
  std::vector<double> vofs; 
  
  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm = getGeometricModel();

  // is there an appropriate region by this name?
  Teuchos::RCP<const AmanziGeometry::Region> rgn;
  try {
    rgn = gm->FindRegion(setname);
  } catch (...) {
    // Catch
  }

  if (rgn == Teuchos::null) {
    std::stringstream ss;
    ss << "Geometric model has no region named \"" << setname <<"\", kind=" << static_cast<int>(kind) << "\n";
    Errors::Message msg(ss.str());
    Exceptions::amanzi_throw(msg);
  }

  std::string setname_internal = setname + std::to_string(static_cast<int>(kind));
  if (sets_.find(setname_internal) != sets_.end()) {
    setents = sets_[setname_internal];
  }

  // set does not exist - build it
  if (sets_.find(setname_internal) == sets_.end()) {
    sets_[setname_internal] = build_set_(rgn, kind);
    setents = sets_[setname_internal];
  }

  // extract entities only of specified parallel type
  int m(0);

  if (ptype == Parallel_type::ALL) {
    m = setents.size();
  } else if (ptype == Parallel_type::OWNED) {
    int nents = getNumEntities(kind, ptype); 
    for (int n = 0; n < setents.size(); ++n) {
      if ((setents)[n] < nents) (setents)[m++] = (setents)[n];
    }
  } else if (ptype == Parallel_type::GHOST) {
    int nents = getNumEntities(kind, Parallel_type::OWNED); 
    for (int n = 0; n < setents.size(); ++n) {
      if ((setents)[n] >= nents) (setents)[m++] = (setents)[n];
    }
  }
    
  setents.resize(m);
      
  // Check if there were no entities left on any processor after
  // extracting the appropriate category of entities
  int mglob; 
  getComm()->SumAll(&m, &mglob, 1);
  
  if (mglob == 0) {
    std::stringstream ss;
    ss << "Could not retrieve mesh entities of kind=" << static_cast<int>(kind) << " for set \"" << setname << "\"" << std::endl;
    Errors::Message msg(ss.str());
    Exceptions::amanzi_throw(msg);
  }

  return std::make_pair(setents, vofs);

}


/* ******************************************************************
* Create a set from similar set in mesh
****************************************************************** */
Entity_ID_List MeshExtractedManifold::build_set_(
    const Teuchos::RCP<const AmanziGeometry::Region>& rgn, const Entity_kind kind) const
{
  // modify rgn/set name by prefixing it with the type of entity requested
  std::string internal_name = rgn->get_name() + std::to_string(static_cast<int>(kind));

  // create entity set based on the region definition  
  Entity_ID_List mset;

  // special processing of regions
  if (rgn->get_type() == AmanziGeometry::RegionType::LABELEDSET) {
    if (parent_labeledsets_.find(internal_name) == parent_labeledsets_.end()) {
      Entity_kind kind_p(Entity_kind::NODE);
      if (kind == Entity_kind::CELL) kind_p = Entity_kind::FACE;
      if (kind == Entity_kind::FACE) kind_p = Entity_kind::EDGE;

      Entity_ID_List setents;
      TryExtension_(rgn->get_name(), kind_p, kind, &setents);
    }

    const auto& ids_p = parent_labeledsets_[internal_name];
    for (auto it = ids_p.begin(); it != ids_p.end(); ++it) {
      mset.push_back(parent_to_entid_[kind][*it]);
    }
    return mset;
  }

  // generic algorithm
  int nents_wghost = getNumEntities(kind, Parallel_type::ALL);
  if (rgn->get_type() == AmanziGeometry::RegionType::ALL)  {
    for (int n = 0; n < nents_wghost; ++n) {
      mset.push_back(n);
    }
    return mset;
  }

  bool missing(false);

  // check if this set exists in parent mesh
  if (flattened_ && rgn->get_type() != AmanziGeometry::RegionType::LOGICAL) {
    mset = build_from_parent_(rgn->get_name(), kind);
    if (mset.size() > 0) return mset;
  }

  if (kind == Entity_kind::CELL)
    mset = build_set_cells_(rgn, &missing);
  else if (kind == Entity_kind::FACE) 
    mset = build_set_faces_(rgn, &missing);
  else if (kind == Entity_kind::NODE) 
    mset = build_set_nodes_(rgn, &missing);

  // check if this set exists in parent mesh
  if (missing && !flattened_ && rgn->get_type() != AmanziGeometry::RegionType::LOGICAL) {
    mset = build_from_parent_(rgn->get_name(), kind);
    if (mset.size() == 0) missing = true;
  }

  if (missing) {
    if (vo_.get() && vo_->os_OK(Teuchos::VERB_HIGH)) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *(vo_->os()) << "Requested entities of kind=" << static_cast<int>(kind) 
        << " on region=\"" << rgn->get_name() << "\" of type " << rgn->get_type()  
        << " and dimension " << rgn->get_manifold_dimension() << ". Result is an empty set.\n";
    }
  }

  if (rgn->get_type() == AmanziGeometry::RegionType::LOGICAL) {
    Teuchos::RCP<const AmanziGeometry::GeometricModel> gm = getGeometricModel();

    auto boolrgn = Teuchos::rcp_static_cast<const AmanziGeometry::RegionLogical>(rgn);
    const std::vector<std::string> rgn_names = boolrgn->get_component_regions();
    int nregs = rgn_names.size();
    
    std::vector<std::set<Entity_ID> > msets;
    
    for (int r = 0; r < nregs; r++) {
      auto rgn1 = gm->FindRegion(rgn_names[r]);

      // Did not find the rgn
      if (rgn1 == Teuchos::null) {
        std::stringstream ss;
        ss << "Geometric model has no region named " << rgn_names[r];
        Errors::Message msg(ss.str());
        Exceptions::amanzi_throw(msg);
      }
        
      std::string setname_internal = rgn1->get_name() + std::to_string(static_cast<int>(kind));

      if (sets_.find(setname_internal) == sets_.end())
        sets_[setname_internal] = build_set_(rgn1, kind); 
     
      auto it = sets_.find(setname_internal);
      msets.push_back(std::set<Entity_ID>(it->second.begin(), it->second.end())); 
    }

    if (boolrgn->get_operation() == AmanziGeometry::BoolOpType::UNION) {
      for (int n = 1; n < msets.size(); ++n) {
        for (auto it = msets[n].begin(); it != msets[n].end(); ++it)
          msets[0].insert(*it);
      }
    }

    else if (boolrgn->get_operation() == AmanziGeometry::BoolOpType::SUBTRACT) {
      for (int n = 2; n < msets.size(); ++n) {
        for (auto it = msets[n].begin(); it != msets[n].end(); ++it)
          msets[0].insert(*it);
      }

      for (auto it = msets[1].begin(); it != msets[1].end(); ++it)
        msets[0].erase(*it);
    }

    else if (boolrgn->get_operation() == AmanziGeometry::BoolOpType::COMPLEMENT) {
      for (int n = 1; n < msets.size(); ++n) {
        for (auto it = msets[n].begin(); it != msets[n].end(); ++it)
          msets[0].insert(*it);
      } 

      auto tmp = msets[0];
      msets[0].clear();
      for (int n = 0; n < nents_wghost; ++n)
        if (tmp.find(n) == tmp.end()) msets[0].insert(n);
    } 

    else if (boolrgn->get_operation() == AmanziGeometry::BoolOpType::INTERSECT) {
      for (int n = 1; n < msets.size(); ++n) {
        for (auto it = msets[n].begin(); it != msets[n].end(); ++it)
          if (msets[0].find(*it) == msets[0].end()) msets[0].erase(*it);
      }
    } 

    for (auto it = msets[0].begin(); it != msets[0].end(); ++it)
      mset.push_back(*it);
  }

  return mset;
}


/* ******************************************************************
* Create a cell-set from similar set in mesh
****************************************************************** */
Entity_ID_List MeshExtractedManifold::build_set_cells_(
    const Teuchos::RCP<const AmanziGeometry::Region>& rgn, bool* missing) const
{
  int space_dim = getSpaceDimension();
  int manifold_dim = getManifoldDimension();

  Entity_ID_List mset;

  int ncells_wghost = getNumEntities(Entity_kind::CELL, Parallel_type::ALL);              
  *missing = false;

  if (rgn->get_type() == AmanziGeometry::RegionType::BOX ||
      rgn->get_type() == AmanziGeometry::RegionType::CYLINDER ||
      rgn->get_type() == AmanziGeometry::RegionType::COLORFUNCTION) {

    for (int c = 0; c < ncells_wghost; ++c)
      if (rgn->inside(getCellCentroid(c))) mset.push_back(c);
  }

  else if (rgn->get_type() == AmanziGeometry::RegionType::POINT) {
    auto rp = Teuchos::rcp_static_cast<const AmanziGeometry::RegionPoint>(rgn)->point();

    for (int c = 0; c < ncells_wghost; ++c)
      assert(false); 
      //if (point_in_cell(rp, c)) mset.push_back(c);
  }

  else if (((rgn->get_type() == AmanziGeometry::RegionType::PLANE) ||
            (rgn->get_type() == AmanziGeometry::RegionType::POLYGON)) && 
           manifold_dim == 2) {
    for (int c = 0; c < ncells_wghost; ++c) {
      std::vector<AmanziGeometry::Point> ccoords(space_dim);
      ccoords = getCellCoordinates(c);

      bool on_plane(true);
      for (int i = 0; i < ccoords.size(); ++i) {
        if (!rgn->inside(ccoords[i])) {
          on_plane = false;
          break;
        }
      }
      if (on_plane) mset.push_back(c);
    }
  }

  else if (rgn->get_type() != AmanziGeometry::RegionType::LOGICAL) {
    *missing = true;
  }

  return mset;
}


/* ******************************************************************
* Create a face-set from similar set in mesh
****************************************************************** */
Entity_ID_List MeshExtractedManifold::build_set_faces_(
    const Teuchos::RCP<const AmanziGeometry::Region>& rgn, bool* missing) const
{
  int space_dim = getSpaceDimension();
  int nfaces_wghost = getNumEntities(Entity_kind::FACE, Parallel_type::ALL);              

  Entity_ID_List mset;

  *missing = false;

  if (rgn->get_type() == AmanziGeometry::RegionType::BOX ||
      rgn->get_type() == AmanziGeometry::RegionType::CYLINDER) {
    for (int f = 0; f < nfaces_wghost; ++f) {
      if (rgn->inside(getFaceCentroid(f))) mset.push_back(f);
    }
  }

  else if (rgn->get_type() == AmanziGeometry::RegionType::PLANE ||
           rgn->get_type() == AmanziGeometry::RegionType::POLYGON) {
    for (int f = 0; f < nfaces_wghost; ++f) {
      std::vector<AmanziGeometry::Point> fcoords(space_dim);
      fcoords = getFaceCoordinates(f);
            
      bool on_plane(true);
      for (int i = 0; i < fcoords.size(); ++i) {
        if (!rgn->inside(fcoords[i])) {
          on_plane = false;
          break;
        }
      }
      if (on_plane) mset.push_back(f);
    }
  }

  else if (rgn->get_type() == AmanziGeometry::RegionType::BOUNDARY)  {
    const auto& fmap = ent_map_owned_[Entity_kind::FACE]; 
    const auto& map = ent_extmap_owned_[Entity_kind::FACE]; 

    int nfaces = map->NumMyElements(); 

    for (int f = 0; f < nfaces; ++f) {
      int lid = fmap->LID(map->GID(f));
      mset.push_back(lid);
    }
  }

  else if (rgn->get_type() != AmanziGeometry::RegionType::LOGICAL) {
    *missing = true;
  }

  return mset;
}


/* ******************************************************************
* Create a node-set from similar set in mesh
****************************************************************** */
Entity_ID_List MeshExtractedManifold::build_set_nodes_(
    const Teuchos::RCP<const AmanziGeometry::Region>& rgn, bool* missing) const
{
  int space_dim = getSpaceDimension();
  int nnodes_wghost = getNumEntities(Entity_kind::NODE, Parallel_type::ALL);

  Entity_ID_List mset;

  *missing = false;

  if (rgn->get_type() == AmanziGeometry::RegionType::BOX ||
      rgn->get_type() == AmanziGeometry::RegionType::PLANE ||
      rgn->get_type() == AmanziGeometry::RegionType::POLYGON ||
      rgn->get_type() == AmanziGeometry::RegionType::CYLINDER ||
      rgn->get_type() == AmanziGeometry::RegionType::POINT) {
    for (int v = 0; v < nnodes_wghost; ++v) {
      AmanziGeometry::Point xp(space_dim);
      xp = getNodeCoordinate(v);
                  
      if (rgn->inside(xp)) {
        mset.push_back(v);
        // Only one node per point rgn
        if (rgn->get_type() == AmanziGeometry::RegionType::POINT) break;      
      }
    }
  }

  else if (rgn->get_type() == AmanziGeometry::RegionType::BOUNDARY)  {
    const auto& vmap = ent_map_owned_[Entity_kind::NODE]; 
    const auto& map = ent_extmap_owned_[Entity_kind::NODE]; 

    int nnodes = ent_map_owned_[Entity_kind::NODE]->NumMyElements(); 

    for (int v = 0; v < nnodes; ++v) {
      int lid = vmap->LID(map->GID(v));
      mset.push_back(lid);
    }
  }

  else if (rgn->get_type() != AmanziGeometry::RegionType::LOGICAL) {
    *missing = true;
  }

  return mset;
}


/* ******************************************************************
* Create from a parent mesh
****************************************************************** */
Entity_ID_List MeshExtractedManifold::build_from_parent_(
  const std::string& rgnname, const Entity_kind kind_d) const
{
  Entity_ID_List mset, mset_tmp;
  Entity_kind kind_p;

  if (kind_d == Entity_kind::CELL) 
    kind_p = Entity_kind::FACE;
  else if (kind_d == Entity_kind::FACE) 
    kind_p = Entity_kind::EDGE;
  else if (kind_d == Entity_kind::NODE)
    kind_p = Entity_kind::NODE;
  else 
    return mset;

  try {
    mset_tmp = parent_mesh_->getSetEntities(rgnname, kind_p, Parallel_type::ALL);

    for (int i = 0; i < mset_tmp.size(); ++i) {
      int f = mset_tmp[i];
      auto it = parent_to_entid_[kind_d].find(f);
      if (it == parent_to_entid_[kind_d].end()) continue;
      mset.push_back(it->second);
    }
  } catch (...) {
    // pass
  }
  return mset;
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
    Entity_ID_List setents;
    std::vector<double> vofs;

    TryExtension_(setname, kind_p, kind_d, &setents);
    if (setents.size() == 0)
      std::tie(setents,vofs) = getSetEntitiesAndVolumeFractions(setname, kind_p, Parallel_type::ALL);

    auto marked_ents = EnforceOneLayerOfGhosts_(setname, kind_p, &setents);

    // extract owned ids
    auto& ids_p = entid_to_parent_[kind_d];
    ids_p.clear();
    for (auto it = marked_ents.begin(); it != marked_ents.end(); ++it) {
      if (it->second == MASTER) ids_p.push_back(it->first);
    }

    nents_owned_[kind_d] = ids_p.size();
    nents_ghost_[kind_d] = marked_ents.size() - ids_p.size();

    // extract ghost ids
    for (auto it = marked_ents.begin(); it != marked_ents.end(); ++it) {
      if (it->second == GHOST) ids_p.push_back(it->first);
    }

    // create reverse ordered map
    auto& ids_d = parent_to_entid_[kind_d];
    ids_d.clear();
    for (int n = 0; n < ids_p.size(); ++n) {
      ids_d[ids_p[n]] = n;
    }
  }
}


/* ******************************************************************
* Epetra maps are structures specifying the global IDs of entities
* owned or used by this processor. This helps Epetra understand
* inter-partition dependencies of the data.
****************************************************************** */
void MeshExtractedManifold::InitEpetraMaps()
{
  std::vector<Entity_kind> kinds_extracted({Entity_kind::CELL, Entity_kind::FACE, Entity_kind::NODE});
  std::vector<Entity_kind> kinds_parent({Entity_kind::FACE, Entity_kind::EDGE, Entity_kind::NODE});

  for (int i = 0; i < 3; ++i) {
    auto kind_d = kinds_extracted[i];
    auto kind_p = kinds_parent[i];

    // compute (discontinuous) owned global ids using the parent map 
    auto parent_map = Teuchos::rcpFromRef(parent_to_entid_[kind_p]);
    auto parent_map_wghost = Teuchos::rcpFromRef(parent_to_entid_[kind_p]);

    int nents = nents_owned_[kind_d];
    int nents_wghost = nents_owned_[kind_d] + nents_ghost_[kind_d];
    auto gids = new int[nents_wghost];

    for (int n = 0; n < nents; ++n) {
      int id = entid_to_parent_[kind_d][n];
      gids[n] = parent_map_wghost->at(id);
    }

    auto subset_map = Teuchos::rcp(new Epetra_Map(-1, nents, gids, 0, *comm_));

    // compute owned + ghost ids using the parent map and the minimum global id
    for (int n = 0; n < nents_wghost; ++n) {
      int id = entid_to_parent_[kind_d][n];
      gids[n] = parent_map_wghost->at(id);
    }

    auto subset_map_wghost = Teuchos::rcp(new Epetra_Map(-1, nents_wghost, gids, 0, *comm_));
    delete [] gids;

    // creare continuous maps
    auto mymesh = Teuchos::rcpFromRef(*this);
    Mesh mcache(mymesh,plist_); 
    assert(false); 
    // Which method to use here? 
    //auto tmp = createContinuousMaps(Teuchos::rcpFromRef(mcache),
    //                                std::make_pair(parent_map, parent_map_wghost),
    //                                std::make_pair(subset_map, subset_map_wghost));

    //ent_map_owned_[kind_d] = tmp.first;
    //ent_map_wghost_[kind_d] = tmp.second;
  }
}


/* ******************************************************************
* Exterior Epetra maps are cannot be always extracted from a
* parent mesh, so we build them explicitly.
****************************************************************** */
void MeshExtractedManifold::InitExteriorEpetraMaps()
{
  int nents = nents_owned_[Entity_kind::FACE];
  Entity_ID_List cells, nodes;

  const auto& fmap_owned = ent_map_owned_[Entity_kind::FACE];
  const auto& fmap_wghost = ent_map_wghost_[Entity_kind::FACE];

  const auto& vmap_owned = ent_map_owned_[Entity_kind::NODE];
  const auto& vmap_wghost = ent_map_wghost_[Entity_kind::NODE];
  
  std::vector<int> gids;
  Epetra_IntVector counts(*fmap_wghost.get()), flags(*vmap_wghost.get());

  // collect information about boundary faces and nodes
  for (int n = 0; n < nents; ++n) {
    getFaceCells(n, Parallel_type::ALL, cells);
    int ncells = cells.size();
    if (ncells == 1) {
      gids.push_back(fmap_owned->GID(n));

      getFaceNodes(n, nodes);
      for (int i = 0; i < nodes.size(); ++i) flags[nodes[i]] = 1;
    }
    counts[n] = ncells; 
  }

  // process faces
  ent_extmap_owned_[Entity_kind::FACE] = Teuchos::rcp(new Epetra_Map(-1, gids.size(), gids.data(), 0, *comm_));
  exterior_face_importer_ = Teuchos::rcp(new Epetra_Import(*ent_extmap_owned_[Entity_kind::FACE], *ent_map_owned_[Entity_kind::FACE]));

#ifdef HAVE_MPI
  {
    auto importer = Teuchos::rcp(new Epetra_Import(*fmap_wghost.get(), *fmap_owned.get()));
    int* vdata;
    counts.ExtractView(&vdata);
    Epetra_IntVector tmp(View, *fmap_owned.get(), vdata);
    counts.Import(tmp, *importer, Insert);
  }
#endif

  int nents_wghost = nents + nents_ghost_[Entity_kind::FACE];
  for (int n = nents; n < nents_wghost; ++n) {
    if (counts[n] == 1) gids.push_back(fmap_wghost->GID(n));
  }

  ent_extmap_wghost_[Entity_kind::FACE] = Teuchos::rcp(new Epetra_Map(-1, gids.size(), gids.data(), 0, *comm_));

  // process nodes
  gids.clear();

  nents = nents_owned_[Entity_kind::NODE];
  for (int n = 0; n < nents; ++n) {
    if (flags[n] == 1) gids.push_back(vmap_owned->GID(n));
  }

  ent_extmap_owned_[Entity_kind::NODE] = Teuchos::rcp(new Epetra_Map(-1, gids.size(), gids.data(), 0, *comm_));

#ifdef HAVE_MPI
  {
    auto importer = Teuchos::rcp(new Epetra_Import(*vmap_wghost.get(), *vmap_owned.get()));
    int* vdata;
    flags.ExtractView(&vdata);
    Epetra_IntVector tmp(View, *vmap_owned.get(), vdata);
    flags.Import(tmp, *importer, Insert);
  }
#endif

  nents_wghost = nents + nents_ghost_[Entity_kind::NODE];
  for (int n = nents; n < nents_wghost; ++n) {
    if (flags[n] == 1) gids.push_back(vmap_wghost->GID(n));
  }

  ent_extmap_wghost_[Entity_kind::NODE] = Teuchos::rcp(new Epetra_Map(-1, gids.size(), gids.data(), 0, *comm_));
}


/* ******************************************************************
* Exception due to limitations of the base mesh framework.
****************************************************************** */
void MeshExtractedManifold::TryExtension_(
    const std::string& setname,
    Entity_kind kind_p, Entity_kind kind_d, Entity_ID_List* setents) const
{
  // labeled set: extract edges
  setents->clear();

  const auto& gm = getGeometricModel();
  if (gm == Teuchos::null) return;

  auto rgn = gm->FindRegion(setname);
  if (rgn->get_type() != AmanziGeometry::RegionType::LABELEDSET) return;

  // populate list of edges
  std::vector<Entity_ID> faceents;
  std::vector<double> vofs;
  std::tie(faceents,vofs) = getSetEntitiesAndVolumeFractions(setname, Entity_kind::FACE, Parallel_type::ALL);
  auto marked_ents = EnforceOneLayerOfGhosts_(setname, Entity_kind::FACE, &faceents);

  Entity_ID_List edges, nodes;
  std::vector<int> dirs;
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

  for (auto it = setents_tmp.begin(); it != setents_tmp.end(); ++it)
    setents->push_back(*it);

  std::string setname_internal = setname + std::to_string(static_cast<int>(kind_d));
  parent_labeledsets_[setname_internal] = *setents;
}


/* ******************************************************************
* Limits the set of parent objects to only one layer of ghosts.
****************************************************************** */
std::map<Entity_ID, int> MeshExtractedManifold::EnforceOneLayerOfGhosts_(
    const std::string& setname, Entity_kind kind, Entity_ID_List* setents) const
{
  // base set is the set of master cells
  Entity_ID_List fullset;
  if (kind != Entity_kind::FACE) {
    std::vector<double> vofs;
    std::tie(fullset,vofs) = getSetEntitiesAndVolumeFractions(setname, Entity_kind::FACE, Parallel_type::ALL);
  } else { 
    fullset = *setents;
  }

  // initial set of entities is defined by master parent faces and is marked as 
  // potential master entities
  Entity_ID_List nodes, edges, faces;
  std::vector<int> dirs;
  int nfaces_owned = parent_mesh_->getNumEntities(Entity_kind::FACE, Parallel_type::OWNED);
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
    const auto& fmap = parent_to_entid_[Entity_kind::FACE];

    int nowned = parent_mesh_->getNumEntities(Entity_kind::FACE, Parallel_type::OWNED);
    int gidmax = std::max_element(fmap.begin(), fmap.end())->second; 

    for (auto it = edgeset.begin(); it != edgeset.end(); ++it) {
      if (it->second == MASTER + GHOST) {
        parent_mesh_->getEdgeFaces(it->first, Parallel_type::ALL, faces);
        int nfaces = faces.size();

        // compare maximum global ids for owned and all faces living on manifold
        Entity_ID gid_owned_min(gidmax + 1), gid_wghost_min(gidmax + 1);
        for (int n = 0; n < nfaces; ++n) {
          Entity_ID f = faces[n];
          if (auxset.find(f) != auxset.end()) {
            gid_wghost_min = std::min(gid_wghost_min, fmap.at(f));
            if (f < nowned)
              gid_owned_min = std::min(gid_owned_min, fmap.at(f));
          }
        }
        it->second = (gid_wghost_min == gid_owned_min) ? MASTER : GHOST;
      }
    }
    return edgeset;
  } else if (kind == Entity_kind::NODE) {
    const auto& fmap = parent_to_entid_[Entity_kind::FACE];

    int nowned = parent_mesh_->getNumEntities(Entity_kind::FACE, Parallel_type::OWNED);
    int gidmax = std::max_element(fmap.begin(), fmap.end())->second; 

    for (auto it = nodeset.begin(); it != nodeset.end(); ++it) {
      if (it->second == MASTER + GHOST) {
        parent_mesh_->getNodeFaces(it->first, Parallel_type::ALL, faces);
        int nfaces = faces.size();

        // compare maximum global ids for owned and all faces living on manifold
        Entity_ID gid_owned_min(gidmax + 1), gid_wghost_min(gidmax + 1);
        for (int n = 0; n < nfaces; ++n) {
          Entity_ID f = faces[n];
          if (auxset.find(f) != auxset.end()) {
            gid_wghost_min = std::min(gid_wghost_min, fmap.at(f));
            if (f < nowned)
              gid_owned_min = std::min(gid_owned_min, fmap.at(f));
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
* Statistics
****************************************************************** */
void MeshExtractedManifold::PrintSets_() const
{
  if (vo_.get() && vo_->os_OK(Teuchos::VERB_LOW)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *(vo_->os()) << "sets: ";
    for (auto it : sets_) {
      int i1, i0(it.second.size());
      getComm()->SumAll(&i0, &i1, 1);
      *(vo_->os()) << "\"" << it.first << "\" (" << i1 << ") ";
    }
    *(vo_->os())  << "\n";

    *(vo_->os()) << "parent labeledsets: ";
    for (auto it : parent_labeledsets_) {
      int i1, i0(it.second.size());
      getComm()->SumAll(&i0, &i1, 1);
      *(vo_->os()) << "\"" << it.first << "\" (" << i1 << ") ";
    }
    *(vo_->os())  << "\n";
  }
}

}  // namespace AmanziMesh
}  // namespace Amanzi


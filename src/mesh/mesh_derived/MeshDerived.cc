/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

#include <utility>

// Amanzi
#include "dbc.hh"
#include "errors.hh"
#include "RegionLogical.hh"
#include "RegionPoint.hh"
#include "VerboseObject.hh"

// Amanzi::Mesh
#include "BlockMapUtils.hh"
#include "MeshDerived.hh"

namespace Amanzi {
namespace AmanziMesh {

/* ******************************************************************
* Light-weigthed constructor
****************************************************************** */
MeshDerived::MeshDerived(
    const Teuchos::RCP<const Mesh>& parent_mesh,
    const std::string& setname, 
    const Entity_kind entity_kind,
    const Comm_ptr_type& comm,
    const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
    const Teuchos::RCP<const Teuchos::ParameterList>& plist,
    bool request_faces, bool request_edges)
  : Mesh(comm, gm, plist, request_faces, request_edges),
    parent_mesh_(parent_mesh)
{
  set_space_dimension(parent_mesh_->space_dimension());
  set_manifold_dimension(parent_mesh_->space_dimension() - 1);

  InitParentMaps(setname); 
  InitEpetraMaps(); 
}


/* ******************************************************************
* Get cell type
****************************************************************** */
Cell_type MeshDerived::cell_get_type(const Entity_ID c) const
{
  int nfaces = cell_face_ids_[c].size();
  return (nfaces == 4) ? QUAD : ((nfaces == 3) ? TRI : POLYGON);
}


/* ******************************************************************
* Number of OWNED, GHOST or ALL entities of different types
*
* Number of entities of any kind (cell, face, edge, node) and in a
* particular category (OWNED, GHOST, ALL)
****************************************************************** */
unsigned int MeshDerived::num_entities(const Entity_kind kind, 
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
void MeshDerived::cell_get_faces_and_dirs_internal_(
    const Entity_ID c,
    Entity_ID_List *faces, std::vector<int> *fdirs, const bool ordered) const
{
  int fp = entid_to_parent_[CELL][c];
  parent_mesh_->face_get_edges_and_dirs(fp, faces, fdirs);
  int nfaces = faces->size();

  for (int i = 0; i < nfaces; ++i) {
    int f = (*faces)[i];
    (*faces)[i] = parent_to_entid_[CELL][f];
  }
}


/* ******************************************************************
* Connectivity list: face -> nodes
****************************************************************** */
void MeshDerived::face_get_nodes(const Entity_ID f, Entity_ID_List *nodes) const {
  Entity_ID v0, v1;

  int ep = entid_to_parent_[FACE][f];
  parent_mesh_->edge_get_nodes(ep, &v0, &v1);

  nodes->clear();
  nodes->push_back(parent_to_entid_[NODE][v0]);
  nodes->push_back(parent_to_entid_[NODE][v1]);
}


/* ******************************************************************
* Connectivity list: face -> edges
****************************************************************** */
void MeshDerived::face_get_edges_and_dirs_internal_(
    const Entity_ID f,
    Entity_ID_List *edges, std::vector<int> *edirs, const bool ordered) const
{
  Entity_ID v0, v1;

  int ep = entid_to_parent_[FACE][f];
  parent_mesh_->edge_get_nodes(ep, &v0, &v1);

  edges->clear();
  edges->push_back(parent_to_entid_[NODE][v0]);
  edges->push_back(parent_to_entid_[NODE][v1]);
  edirs->resize(2, 1);
}


/* ******************************************************************
* Connectivity list: face -> cells
****************************************************************** */
void MeshDerived::face_get_cells_internal_(
    const Entity_ID f,
    const Parallel_type ptype, Entity_ID_List *cells) const
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  int ep = entid_to_parent_[FACE][f];
  // parent_mesh_->edge_get_faces_and_dirs(ep, ptype, &faces, &dirs);
  AMANZI_ASSERT(false);
  int nfaces = faces.size();

  cells->clear();
  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    auto it = parent_to_entid_[CELL].find(f);
    if (it != parent_to_entid_[CELL].end()) cells->push_back(it->second);
  }
}


/* ******************************************************************
* The node coordinates are returned in in arbitrary order
****************************************************************** */
void MeshDerived::cell_get_coordinates(const Entity_ID c,
                                       std::vector<AmanziGeometry::Point> *vxyz) const
{
  Entity_ID_List nodes; 
  parent_mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  AmanziGeometry::Point p(space_dimension());

  vxyz->clear();
  for (int i = 0; i < nnodes; ++i) {
    parent_mesh_->node_get_coordinates(nodes[i], &p);
    vxyz->push_back(p);
  }
}


/* ******************************************************************
* Face coordinates use convention same as in face_to_nodes()
****************************************************************** */
void MeshDerived::face_get_coordinates(const Entity_ID f,
                                       std::vector<AmanziGeometry::Point>* vxyz) const
{
  Entity_ID_List nodes; 
  parent_mesh_->face_get_nodes(f, &nodes);
  int nnodes = nodes.size();

  AmanziGeometry::Point p(space_dimension());

  vxyz->clear();
  for (int i = 0; i < nnodes; ++i) {
    parent_mesh_->node_get_coordinates(nodes[i], &p);
    vxyz->push_back(p);
  }
}


/* ******************************************************************
* Get list of entities of type 'ptype' in set specified by setname
****************************************************************** */
void MeshDerived::get_set_entities_and_vofs(const std::string setname, 
                                            const Entity_kind kind, 
                                            const Parallel_type ptype, 
                                            std::vector<Entity_ID> *setents,
                                            std::vector<double> *vofs) const
{
  assert(setents != nullptr);

  setents->clear();
  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm = geometric_model();

  // is there an appropriate region by this name?
  Teuchos::RCP<const AmanziGeometry::Region> rgn;
  try {
    rgn = gm->FindRegion(setname);
  } catch (...) {
    valid_set_name(setname, kind);
  }

  if (rgn == Teuchos::null) {
    std::stringstream ss;
    ss << "Geometric model has no region named \"" << setname <<"\", kind=" << kind << "\n";
    Errors::Message msg(ss.str());
    Exceptions::amanzi_throw(msg);
  }

  std::string setname_internal = setname + std::to_string(kind);
  if (sets_.find(setname_internal) != sets_.end()) {
    Entity_ID_List block;
    std::vector<double> vofs;
    parent_mesh_->get_set_entities_and_vofs(setname, AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL, &block, &vofs);

    setents->clear();
    for (int c = 0; c < num_entities(CELL, Parallel_type::ALL); ++c) {
      int f = entid_to_parent_[CELL][c];
      if (std::find(block.begin(), block.end(), f) != block.end()) setents->push_back(c);
    }
    sets_[setname_internal] = *setents;
  }

  // all attempts to find the set failed so it must not exist - build it
  if (sets_.find(setname_internal) == sets_.end()) {
    *setents = sets_[setname_internal] = build_set_(rgn, kind);
  }

  // reset and count to get the real number
  int m(0);

  if (ptype == Parallel_type::ALL) {
    m = setents->size();
  } else if (ptype == Parallel_type::OWNED) {
    int nents = num_entities(kind, ptype); 
    for (int n = 0; n < setents->size(); ++n) {
      if ((*setents)[n] < nents) (*setents)[m++] = (*setents)[n];
    }
  } else if (ptype == Parallel_type::GHOST) {
    int nents = num_entities(kind, Parallel_type::OWNED); 
    for (int n = 0; n < setents->size(); ++n) {
      if ((*setents)[n] >= nents) (*setents)[m++] = (*setents)[n];
    }
  }
    
  setents->resize(m);
      
  // Check if there were no entities left on any processor after
  // extracting the appropriate category of entities
  int mglob; 
  get_comm()->SumAll(&m, &mglob, 1);
  
  if (mglob == 0) {
    std::stringstream ss;
    ss << "Could not retrieve any mesh entities of type " << kind << " for set " << setname << std::endl;
    Errors::Message msg(ss.str());
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
* Create a set from similar set in mesh
****************************************************************** */
Entity_ID_List MeshDerived::build_set_(
    const Teuchos::RCP<const AmanziGeometry::Region>& rgn, const Entity_kind kind) const
{
  int manifold_dim = manifold_dimension();
  int space_dim_ = space_dimension();
  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm = geometric_model();

  // modify rgn/set name by prefixing it with the type of entity requested
  std::string internal_name = rgn->name() + std::to_string(kind);

  // create entity set based on the region definition  
  int ncells_wghost = num_entities(CELL, Parallel_type::ALL);              
  int nfaces_wghost = num_entities(FACE, Parallel_type::ALL);              
  int nnodes_wghost = num_entities(NODE, Parallel_type::ALL);
  Entity_ID_List mset;

  switch (kind) {      
  // create a set of cells
  case CELL:
    if (rgn->type() == AmanziGeometry::BOX ||
        rgn->type() == AmanziGeometry::CYLINDER ||
        rgn->type() == AmanziGeometry::COLORFUNCTION) {

      for (int c = 0; c < ncells_wghost; ++c)
        if (rgn->inside(cell_centroid(c))) mset.push_back(c);
    }

    else if (rgn->type() == AmanziGeometry::ALL)  {
      for (int c = 0; c < ncells_wghost; ++c)
        mset.push_back(c);
    }

    else if (rgn->type() == AmanziGeometry::POINT) {
      auto rp = Teuchos::rcp_static_cast<const AmanziGeometry::RegionPoint>(rgn)->point();

      for (int c = 0; c < ncells_wghost; ++c)
        if (point_in_cell(rp, c)) mset.push_back(c);
    }

    else if (((rgn->type() == AmanziGeometry::PLANE) ||
              (rgn->type() == AmanziGeometry::POLYGON)) && 
             manifold_dim == 2) {
      for (int c = 0; c < ncells_wghost; ++c) {
        std::vector<AmanziGeometry::Point> ccoords(space_dim_);
        cell_get_coordinates(c, &ccoords);

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

    else if (rgn->type() == AmanziGeometry::LOGICAL) {
      // will process later in this subroutine
    }

    else {
      if (vo_.get() && vo_->os_OK(Teuchos::VERB_HIGH)) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *(vo_->os()) << "Requested CELLS on rgn " << rgn->name() 
            << " of type " << rgn->type()  
            << " and dimension " << rgn->manifold_dimension() << ".\n"
            << "This request will result in an empty set";
      }
    }

    break;

  // create a set of faces
  case FACE:
    if (rgn->type() == AmanziGeometry::BOX ||
        rgn->type() == AmanziGeometry::CYLINDER) {
      for (int f = 0; f < nfaces_wghost; ++f) {
        if (rgn->inside(face_centroid(f))) mset.push_back(f);
      }
    }

    else if (rgn->type() == AmanziGeometry::ALL)  {
      for (int f = 0; f < nfaces_wghost; ++f)
        mset.push_back(f);
    }

    else if (rgn->type() == AmanziGeometry::PLANE ||
             rgn->type() == AmanziGeometry::POLYGON) {
      for (int f = 0; f < nfaces_wghost; ++f) {
        std::vector<AmanziGeometry::Point> fcoords(space_dim_);
        face_get_coordinates(f, &fcoords);
            
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

    else if (rgn->type() == AmanziGeometry::LOGICAL) {
      // Will handle it later in the routine
    }

    else if (rgn->type() == AmanziGeometry::BOUNDARY)  {
      const Epetra_Map& fmap = face_map(true); 
      const Epetra_Map& map = exterior_face_map(true); 

      int nfaces = map.NumMyElements(); 

      for (int f = 0; f < nfaces; ++f) {
        int lid = fmap.LID(map.GID(f));
        mset.push_back(lid);
      }
    }

    else {
      if (vo_.get() && vo_->os_OK(Teuchos::VERB_HIGH)) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *(vo_->os()) << "Requested FACES on rgn " << rgn->name()
            << " of type " << rgn->type() << " and dimension "
            << rgn->manifold_dimension() << ".\n" 
            << "This request will result in an empty set\n";
      }
    }

    break;

  // create a set of nodes
  case NODE:
    if (rgn->type() == AmanziGeometry::BOX ||
        rgn->type() == AmanziGeometry::PLANE ||
        rgn->type() == AmanziGeometry::POLYGON ||
        rgn->type() == AmanziGeometry::CYLINDER ||
        rgn->type() == AmanziGeometry::POINT) {
      for (int v = 0; v < nnodes_wghost; ++v) {
        AmanziGeometry::Point xp(space_dim_);
        node_get_coordinates(v, &xp);
                  
        if (rgn->inside(xp)) {
          mset.push_back(v);

          // Only one node per point rgn
          if (rgn->type() == AmanziGeometry::POINT) break;      
        }
      }
    }

    else if (rgn->type() == AmanziGeometry::ALL)  {
      for (int v = 0; v < nnodes_wghost; ++v)
        mset.push_back(v);
    }

    else if (rgn->type() == AmanziGeometry::LOGICAL) {
      // We will handle it later in the routine
    }

    else if (rgn->type() == AmanziGeometry::BOUNDARY)  {
      const Epetra_Map& vmap = node_map(true); 
      const Epetra_Map& map = exterior_node_map(true); 

      int nnodes = map.NumMyElements(); 

      for (int v = 0; v < nnodes; ++v) {
        int lid = vmap.LID(map.GID(v));
        mset.push_back(lid);
      }
    }

    else {
      if (vo_.get() && vo_->os_OK(Teuchos::VERB_HIGH)) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *(vo_->os()) << "Requested POINTS on rgn " << rgn->name() 
            << " of type " << rgn->type() << " and dimension " 
            << rgn->manifold_dimension() << ".\n" 
            << "This request will result in an empty set\n";
      }
    }
      
    break;

  default:
    break;
  }


  if (rgn->type() == AmanziGeometry::LOGICAL) {
    auto boolrgn = Teuchos::rcp_static_cast<const AmanziGeometry::RegionLogical>(rgn);
    const std::vector<std::string> rgn_names = boolrgn->component_regions();
    int nregs = rgn_names.size();
    
    std::vector<Entity_ID_List> msets;
    std::vector<Teuchos::RCP<const AmanziGeometry::Region> > rgns;
    
    for (int r = 0; r < nregs; r++) {
      auto rgn1 = gm->FindRegion(rgn_names[r]);

      // Did not find the rgn
      if (rgn1 == Teuchos::null) {
        std::stringstream ss;
        ss << "Geometric model has no region named " << rgn_names[r];
        Errors::Message msg(ss.str());
        Exceptions::amanzi_throw(msg);
      }
        
      rgns.push_back(rgn1);
      std::string setname_internal = rgn1->name() + std::to_string(kind);
      if (sets_.find(setname_internal) != sets_.end())
        msets.push_back(build_set_(rgn1, kind)); 
      else
        msets.push_back(sets_[setname_internal]);
    }

    // Check the entity types of the sets are consistent with the
    // entity type of the requested set
    /*
    for (int ms = 0; ms < msets.size(); ms++) {
      if (MSet_EntDim(msets[ms]) != enttype) {
        Errors::Message msg("Amanzi cannot operate on sets of different entity types");
        Exceptions::amanzi_throw(msg);               
      }
    }
    */
    
    if (boolrgn->operation() == AmanziGeometry::UNION) {
      for (int ms = 0; ms < msets.size(); ms++) {
        for (auto it = msets[ms].begin(); it != msets[ms].end(); ++it) {
          mset.push_back(*it);
        }
      }
    }
    else if (boolrgn->operation() == AmanziGeometry::SUBTRACT ||
             boolrgn->operation() == AmanziGeometry::COMPLEMENT ||
             boolrgn->operation() == AmanziGeometry::INTERSECT) {
      Errors::Message msg("Missing code for logical operation");
      Exceptions::amanzi_throw(msg);
    }
  }

  return mset;
}


/* ******************************************************************
* Create internal maps for child->parent
****************************************************************** */
void MeshDerived::InitParentMaps(const std::string& setname)
{
  nents_owned_.clear();
  nents_ghost_.clear();
  entid_to_parent_.clear();
  parent_to_entid_.clear();

  std::vector<Entity_kind> kinds_derived({CELL, FACE, NODE});
  std::vector<Entity_kind> kinds_parent({FACE, EDGE, NODE});

  for (int i = 0; i < 3; ++i) {
    auto kind_d = kinds_derived[i];
    auto kind_p = kinds_parent[i];

    Entity_ID_List setents;
    std::vector<double> vofs;
    parent_mesh_->get_set_entities_and_vofs(setname, kind_p, Parallel_type::ALL, &setents, &vofs);

    // extract owned ids
    int nowned_p = parent_mesh_->num_entities(kind_p, Parallel_type::OWNED);

    auto& ids_p = entid_to_parent_[kind_d];
    for (int n = 0; n < setents.size(); ++n) {
      if (setents[n] < nowned_p) 
        ids_p.push_back(setents[n]);
    }

    nents_owned_[kind_d] = ids_p.size();
    nents_ghost_[kind_d] = setents.size() - ids_p.size();

    // extract ghost ids
    for (int n = 0; n < setents.size(); ++n) {
      if (setents[n] >= nowned_p)
        ids_p.push_back(setents[n]);
    }

    // create reverse ordered map
    auto& ids_d = parent_to_entid_[kind_d];
    for (int n = 0; n < setents.size(); ++n) {
      ids_d[setents[n]] = n;
    }
  }
}


/* ******************************************************************
* Epetra maps are structures specifying the global IDs of entities 
* owned or used by this processor. This helps Epetra understand 
* inter-partition dependencies of the data.
****************************************************************** */
void MeshDerived::InitEpetraMaps()
{
  std::vector<Entity_kind> kinds_derived({CELL, FACE, NODE});
  std::vector<Entity_kind> kinds_parent({FACE, EDGE, NODE});

  for (int i = 0; i < 3; ++i) {
    auto kind_d = kinds_derived[i];
    auto kind_p = kinds_parent[i];

    std::vector<Entity_ID> setents;
    std::vector<double> vofs;
    get_set_entities_and_vofs("ALL", kind_d, Parallel_type::OWNED, &setents, &vofs);

    // compute owned global ids using the parent map and the minimum global id
    Teuchos::RCP<const Epetra_BlockMap> parent_map = Teuchos::rcpFromRef(parent_mesh_->map(kind_p, false));

    int nents = nents_owned_[kind_d];
    int nents_wghost = nents_owned_[kind_d] + nents_ghost_[kind_d];
    auto gids = new int[nents_wghost];

    for (int n = 0; n < nents; ++n) {
      int id = entid_to_parent_[kind_d][n];
      gids[n] = parent_map->GID(id);
    }

    auto subset_map = Teuchos::rcp(new Epetra_Map(-1, nents, gids, 0, *comm_));

    // compute owned + ghost ids using the parent map and the minimum global id
    Teuchos::RCP<const Epetra_BlockMap> parent_map_wghost = Teuchos::rcpFromRef(parent_mesh_->map(kind_p, true));

    for (int n = 0; n < nents_wghost; ++n) {
      int id = entid_to_parent_[kind_d][n];
      gids[n] = parent_map_wghost->GID(id);
    }

    auto subset_map_wghost = Teuchos::rcp(new Epetra_Map(-1, nents_wghost, gids, 0, *comm_));
    delete [] gids;

    // creare continuous maps
    auto mymesh = Teuchos::rcpFromRef(*this);
    auto tmp = AmanziMesh::CreateContinuousMaps(mymesh,
                                    std::make_pair(parent_map, parent_map_wghost),
                                    std::make_pair(subset_map, subset_map_wghost));

    ent_map_owned_[kind_d] = tmp.first;
    ent_map_wghost_[kind_d] = tmp.second;
  }
}

}  // namespace AmanziMesh
}  // namespace Amanzi


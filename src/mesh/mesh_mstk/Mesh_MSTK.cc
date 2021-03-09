/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, others
*/

//! Implementation of the Mesh interface leveraging MSTK.

#include "dbc.hh"
#include "errors.hh"

#include "VerboseObject.hh"
#include "Point.hh"
#include "GeometricModel.hh"

#include "RegionLogical.hh"
#include "RegionPoint.hh"
#include "RegionLabeledSet.hh"

#include "Mesh_MSTK.hh"

using namespace std;

namespace Amanzi {

namespace AmanziMesh {

//--------------------------------------
// Constructor - load up mesh from file
//--------------------------------------
Mesh_MSTK::Mesh_MSTK(const std::string& filename,
                     const Comm_ptr_type& comm,
                     const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
                     const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : MeshFramework(comm, gm, plist),
    cells_initialized_(false),
    faces_initialized_(false),
    edges_initialized_(false)
{
  read_plist_();

  if (vo_.get() && vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *(vo_->os()) << "Construct mesh from file" << std::endl;
  }

  // Pre-processing (init, MPI queries etc)
  int space_dim = 3;
  pre_create_steps_(space_dim);
  init_mesh_from_file_(filename);

  int cell_dim = MESH_Num_Regions(mesh_) ? 3 : 2;
  int max;
  comm->MaxAll(&cell_dim,&max,1);

  if (max != cell_dim) {
    Errors::Message mesg("cell dimension on this processor is different from max cell dimension across all processors");
    Exceptions::amanzi_throw(mesg);
  }
  set_manifold_dimension(cell_dim);

  if (cell_dim == 2 && space_dim == 3) {
    // Check if this is a completely planar mesh
    // in which case one can label the space dimension as 2
    //
    // cannot use getNodeCoordinate() yet!
    MVertex_ptr mv = nullptr, mv0 = MESH_Vertex(mesh_,0);
    double vxyz[3], z0;
    MV_Coords(mv0,vxyz);
    z0 = vxyz[2];

    bool planar = true;
    int idx = 0;
    while ((mv = MESH_Next_Vertex(mesh_,&idx))) {
      MV_Coords(mv,vxyz);
      if (z0 != vxyz[2]) {
        planar = false;
        break;
      }
    }

    if (planar) space_dim = 2;
    comm->MaxAll(&space_dim,&max,1);
    space_dim = max;
    set_space_dimension(space_dim);
  }

  // Verify mesh and geometric model compatibility
  if (gm != Teuchos::null && gm->dimension() != get_space_dimension()) {
    Errors::Message msg("Geometric model and mesh have different dimensions.");
    Exceptions::amanzi_throw(msg);
  }

  // Do all the processing required for setting up the mesh for Amanzi
  post_create_steps_();
}


//--------------------------------------
// Construct a 3D regular hexahedral mesh internally
//--------------------------------------
Mesh_MSTK::Mesh_MSTK(const double x0, const double y0, const double z0,
                     const double x1, const double y1, const double z1,
                     const unsigned int nx, const unsigned int ny,
                     const unsigned int nz,
                     const Comm_ptr_type& comm,
                     const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
                     const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : MeshFramework(comm, gm, plist),
    edges_requested_(false),
    cells_initialized_(false),
    faces_initialized_(false),
    edges_initialized_(false)
{
  read_plist_();

  int ok = 1;
  int space_dimension = 3;
  pre_create_steps_(space_dimension);

  if (vo_.get() && vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << "Construct mesh from low/hi coords - 3D" << std::endl;
  }

  if (serial_run) {
    // Load serial mesh
    mesh_ = MESH_New(F1);
    ok &= generate_regular_mesh_(mesh_,x0,y0,z0,x1,y1,z1,nx,ny,nz);
    set_manifold_dimension(3);
    myprocid = 0;

  } else {
    Mesh_ptr globalmesh = nullptr;
    int topo_dim=3; // What is the topological dimension of the mesh
    int ring = 1; // One layer of ghost cells in parallel meshes
    int with_attr = 1;  // update of attributes in parallel meshes
    int del_inmesh = 1; // delete input mesh as soon as possible
    int method = static_cast<int>(partitioner_);

    if (myprocid == 0) {
      globalmesh = MESH_New(F1);
      ok &= generate_regular_mesh_(globalmesh,x0,y0,z0,x1,y1,z1,nx,ny,nz);
      topo_dim = (MESH_Num_Regions(globalmesh) == 0) ? 2 : 3;

    } else {
      globalmesh = nullptr;
    }

    ok = ok & MSTK_Mesh_Distribute(globalmesh,&mesh_,&topo_dim,ring,with_attr,
                                   method,del_inmesh,mpicomm_);
    set_manifold_dimension(topo_dim);
  }

  if (!ok) {
    Errors::Message mesg;
    mesg << "Failed to generate mesh on processor " << myprocid;
    Exceptions::amanzi_throw(mesg);
  }

  // Do all the processing required for setting up the mesh for Amanzi
  post_create_steps_();
}


//--------------------------------------
// Construct a 2D regular quadrilateral mesh internally
//--------------------------------------
Mesh_MSTK::Mesh_MSTK(const double x0, const double y0,
                     const double x1, const double y1,
                     const int nx, const int ny,
                     const Comm_ptr_type& comm,
                     const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
                     const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : MeshFramework(comm, gm, plist),
    edges_requested_(false),
    cells_initialized_(false),
    faces_initialized_(false),
    edges_initialized_(false)
{
  read_plist_();

  int ok = 1;
  int space_dim = 2;
  pre_create_steps_(space_dim);

  if (vo_.get() && vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << "Construct mesh from low/hi coords - 2D" << std::endl;
  }

  int topo_dim = space_dim; // What is the topological dimension of the mesh
  set_manifold_dimension(topo_dim);

  if (serial_run) {
    // Load serial mesh
    mesh_ = MESH_New(F1);
    ok &= generate_regular_mesh_(mesh_,x0,y0,x1,y1,nx,ny);
    myprocid = 0;

  } else {
    Mesh_ptr globalmesh = nullptr;
    int ring = 1; // One layer of ghost cells in parallel meshes
    int with_attr = 1;  // update of attributes in parallel meshes
    int del_inmesh = 1; // delete input mesh at the earliest
    int method = static_cast<int>(partitioner_);

    if (myprocid == 0) {
      globalmesh = MESH_New(F1);
      ok &= generate_regular_mesh_(globalmesh,x0,y0,x1,y1,nx,ny);
      topo_dim = (MESH_Num_Regions(globalmesh) == 0) ? 2 : 3;

    } else {
      globalmesh = nullptr;
    }

    ok &= MSTK_Mesh_Distribute(globalmesh,&mesh_,&topo_dim,ring,with_attr,
                                   method,del_inmesh,mpicomm_);
  }

  if (!ok) {
    Errors::Message mesg;
    mesg << "Failed to generate mesh on processor " << myprocid;
    Exceptions::amanzi_throw(mesg);
  }

  // Do all the processing required for setting up the mesh for Amanzi
  post_create_steps_();
}


//---------------------------------------------------------
// Extract MSTK entities from an ID list and make a new MSTK mesh
//---------------------------------------------------------
Mesh_MSTK::Mesh_MSTK(const Teuchos::RCP<const MeshFramework>& parent_mesh,
                     const Entity_ID_List& entity_ids,
                     const Entity_kind entity_kind,
                     const bool flatten,
                     const Comm_ptr_type& comm,
                     const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
                     const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : MeshFramework(comm,
                  gm == Teuchos::null ? parent_mesh->get_geometric_model() : gm,
                  plist == Teuchos::null ? Teuchos::rcp(
                    new Teuchos::ParameterList(*parent_mesh->get_parameter_list())) : plist),
    edges_requested_(false),
    cells_initialized_(false),
    faces_initialized_(false),
    edges_initialized_(false)
{
  read_plist_();
  auto parent_mesh_as_mstk = Teuchos::rcp_dynamic_cast<const Mesh_MSTK>(parent_mesh);
  if (!parent_mesh_as_mstk.get()) {
    Errors::Message mesg("Cannot extract an MSTK mesh from a non-MSTK mesh.");
    Exceptions::amanzi_throw(mesg);
  }
  set_parent(parent_mesh_as_mstk);
  auto parent_mesh_mstk = parent_mesh_as_mstk->mesh_;

  // store pointers to the MESH_XXXFromID functions so that they can
  // be called without a switch statement
  static MEntity_ptr (*MEntFromID[4])(Mesh_ptr,int) =
    {MESH_VertexFromID, MESH_EdgeFromID, MESH_FaceFromID, MESH_RegionFromID};
  MType entity_dim = parent_mesh_as_mstk->entity_kind_to_mtype_(entity_kind);

  // Also make sure that the mesh object can do fast queries on local IDs
  MESH_Enable_LocalIDSearch(parent_mesh_mstk);

  // collect MSTK entities and extract
  int nent = entity_ids.size();
  List_ptr src_ents = List_New(nent);
  for (int i = 0; i < nent; ++i) {
    MEntity_ptr ent = MEntFromID[entity_dim](parent_mesh_mstk,entity_ids[i]+1);
    List_Add(src_ents,ent);
  }
  extract_mstk_mesh_(src_ents, entity_dim, flatten, faces_requested_, edges_requested_);
  List_Delete(src_ents);
}


//---------------------------------------------------------
// Destructor with cleanup
//---------------------------------------------------------
Mesh_MSTK::~Mesh_MSTK() {
  if (faceflip_) delete [] faceflip_;
  if (edgeflip_) delete [] edgeflip_;

  if (owned_verts_) MSet_Delete(owned_verts_);
  if (ghost_verts_) MSet_Delete(ghost_verts_);
  if (owned_edges_) MSet_Delete(owned_edges_);
  if (ghost_edges_) MSet_Delete(ghost_edges_);
  if (owned_faces_) MSet_Delete(owned_faces_);
  if (ghost_faces_) MSet_Delete(ghost_faces_);
  if (owned_cells_) MSet_Delete(owned_cells_);
  if (ghost_cells_) MSet_Delete(ghost_cells_);

  if (entities_deleted_) {
    if (deleted_vertices_) List_Delete(deleted_vertices_);
    if (deleted_edges_) List_Delete(deleted_edges_);
    if (deleted_faces_) List_Delete(deleted_faces_);
    if (deleted_regions_) List_Delete(deleted_regions_);
  }

  MAttrib_Delete(celltype_att_);
  if (vparentatt_) MAttrib_Delete(vparentatt_);
  if (eparentatt_) MAttrib_Delete(eparentatt_);
  if (fparentatt_) MAttrib_Delete(fparentatt_);
  if (rparentatt_) MAttrib_Delete(rparentatt_);

  MESH_Delete(mesh_);
}

//
// parse the parameterlist
//
void Mesh_MSTK::read_plist_()
{
  // extract optional control parameters, but first specify defaults
  partitioner_ = createPartitionerType(plist_->get<std::string>("partitioner", "metis"));
  contiguous_gids_ = plist_->get<bool>("contiguous global ids", true);
  edges_requested_ = plist_->get<bool>("request edges", false);
  faces_requested_ = plist_->get<bool>("request faces", true);
}


//---------------------------------------------------------
// Translate a setname into a special string with decorations
// indicating which type of entity is in that set
//---------------------------------------------------------
std::string
Mesh_MSTK::internal_name_of_set_(
  const AmanziGeometry::RegionLabeledSet& rgn,
  const Entity_kind entity_kind) const
{
  std::string label = rgn.label();
  std::string internal_name;
  if (entity_kind == Entity_kind::CELL)
    internal_name = "matset_" + label;
  else if (entity_kind == Entity_kind::FACE)
    internal_name = "sideset_" + label;
  else if (entity_kind == Entity_kind::EDGE)
    internal_name = "edgeset_not_supported";
  else if (entity_kind == Entity_kind::NODE)
    internal_name = "nodeset_" + label;
  return internal_name;
}


//---------------------------------------------------------
// Get an alternate name (elemset_N instead of matset_N) for sets of type
// Labeled Set and entity kind Cell. For everything else return regular name
//---------------------------------------------------------
std::string
Mesh_MSTK::other_internal_name_of_set_(
  const AmanziGeometry::RegionLabeledSet& rgn,
  const Entity_kind entity_kind) const
{
  std::string label = rgn.label();
  std::string internal_name = "elemset_" + label;
  return internal_name;
}


//---------------------------------------------------------
// Number of OWNED, GHOST or ALL entities of different types
//
// Number of entities of any kind (cell, face, edge, node) and in a
// particular category (OWNED, GHOST, ALL)
//---------------------------------------------------------
std::size_t
Mesh_MSTK::getNumEntities(const Entity_kind kind, const Parallel_type ptype) const
{
  switch (kind) {
  case Entity_kind::NODE:
    switch (ptype) {
    case Parallel_type::OWNED:
      return MSet_Num_Entries(owned_verts_);
      break;
    case Parallel_type::GHOST:
      return !serial_run ? MSet_Num_Entries(ghost_verts_) : 0;
      break;
    case Parallel_type::ALL:
      return MESH_Num_Vertices(mesh_);
      break;
    default:
      return 0;
    }
    break;

  case Entity_kind::EDGE:
    AMANZI_ASSERT(edges_initialized_);
    switch (ptype) {
    case Parallel_type::OWNED:
      return MSet_Num_Entries(owned_edges_);
      break;
    case Parallel_type::GHOST:
      return !serial_run ? MSet_Num_Entries(ghost_edges_) : 0;
      break;
    case Parallel_type::ALL:
      return MESH_Num_Edges(mesh_);
      break;
    default:
      return 0;
    }
    break;

  case Entity_kind::FACE:
    AMANZI_ASSERT(faces_initialized_);
    switch (ptype) {
    case Parallel_type::OWNED:
      return MSet_Num_Entries(owned_faces_);
      break;
    case Parallel_type::GHOST:
      return !serial_run ? MSet_Num_Entries(ghost_faces_) : 0;
      break;
    case Parallel_type::ALL:
      return (get_manifold_dimension() == 2 ? MESH_Num_Edges(mesh_) : MESH_Num_Faces(mesh_));
      break;
    default:
      return 0;
    }
    break;

  case Entity_kind::CELL:
    switch (ptype) {
    case Parallel_type::OWNED:
      return MSet_Num_Entries(owned_cells_);
      break;
    case Parallel_type::GHOST:
      return !serial_run ? MSet_Num_Entries(ghost_cells_) : 0;
      break;
    case Parallel_type::ALL:
      return ((get_manifold_dimension() == 2) ? MESH_Num_Faces(mesh_) : MESH_Num_Regions(mesh_));
      break;
    default:
      return 0;
    }
    break;
  default:
    std::cerr << "Count requested for unknown entity type" << std::endl;
  }
  return 0;
}


//---------------------------------------------------------
// Get cell type
//---------------------------------------------------------
Cell_type Mesh_MSTK::getCellType(const Entity_ID cellid) const
{
  MEntity_ptr cell = cell_id_to_handle_[cellid];
  int ival;
  MEnt_Get_AttVal(cell,celltype_att_,&ival,NULL,NULL);
  return (Cell_type) ival;
}


//---------------------------------------------------------
// Get faces of a cell and directions in which the cell uses the face

// The Amanzi coding guidelines regarding function arguments is purposely
// violated here to allow for a default input argument

// On a distributed mesh, this will return all the faces of the
// cell, OWNED or GHOST. If ordered = true, the faces will be
// returned in a standard order according to Exodus II convention
// for standard cells; in all other situations (ordered = false or
// non-standard cells), the list of faces will be in arbitrary order

// In 3D, direction is 1 if face normal points out of cell
// and -1 if face normal points into cell
// In 2D, direction is 1 if face/edge is defined in the same
// direction as the cell polygon, and -1 otherwise
//---------------------------------------------------------
void Mesh_MSTK::getCellFacesAndDirs_ordered_(const Entity_ID cellid,
                                                Entity_ID_List& faceids,
                                                Entity_Direction_List * const face_dirs) const
{
  if (get_manifold_dimension() == 3) {
    Cell_type celltype = getCellType(cellid);
    if (celltype != Cell_type::UNKNOWN) {
      int lid, nf;
      MEntity_ptr cell = cell_id_to_handle_[cellid];

      List_ptr rfaces = MR_Faces((MRegion_ptr)cell);
      nf = List_Num_Entries(rfaces);

      faceids.resize(nf);
      if (face_dirs) face_dirs->resize(nf);

      /* base face */
      MFace_ptr face0 = nullptr;
      int fdir0 = 0;

      if (celltype == Cell_type::TET || celltype == Cell_type::HEX) {
        face0 = List_Entry(rfaces,0);
        fdir0 = MR_FaceDir_i((MRegion_ptr)cell,0);

      } else if (celltype == Cell_type::PRISM) { /* Find the first triangular face */
        for (int i = 0; i < 5; ++i) {
          MFace_ptr face = List_Entry(rfaces,i);
          if (MF_Num_Edges(face) == 3) {
            face0 = face;
            fdir0 = MR_FaceDir_i((MRegion_ptr)cell,i);
            break;
          }
        }

      } else if (celltype == Cell_type::PYRAMID) { /* Find the quad face */
        for (int i = 0; i < 5; ++i) {
          MFace_ptr face = List_Entry(rfaces,i);
          if (MF_Num_Edges(face) == 4) {
            face0 = face;
            fdir0 = MR_FaceDir_i((MRegion_ptr)cell,i);
            break;
          }
        }
      }

      /* Markers for faces to avoid searching */
      int mkid = MSTK_GetMarker();
      MEnt_Mark(face0,mkid);

      /* Add all lateral faces first (faces adjacent to the base face) */
      List_ptr fedges0 = MF_Edges(face0,!fdir0,0);
      int idx = 0;
      MEdge_ptr fe;
      nf = 0;
      while ((fe = List_Next_Entry(fedges0,&idx))) {
        /* Is there an unprocessed face in this region that is
           adjacent to this edge */
        int idx2 = 0;
        MFace_ptr fadj = nullptr;
        int i = 0;
        while ((fadj = List_Next_Entry(rfaces,&idx2))) {
          if (fadj != face0 && !MEnt_IsMarked(fadj,mkid)) {
            if (MF_UsesEntity(fadj,fe,MEDGE)) {
              lid = MEnt_ID(fadj);
              faceids[nf] = lid-1;

              if (face_dirs) {
                int fdir = (MR_FaceDir_i((MRegion_ptr)cell,i) == 1) ? 1 : -1;
                if (faceflip_[lid-1]) fdir *= -1;
                (*face_dirs)[nf] = fdir;
              }

              MEnt_Mark(fadj,mkid);
              nf++;
            }
          }
          ++i;
        }
      }
      List_Delete(fedges0);

      /* Add the base face */
      lid = MEnt_ID(face0);
      faceids[nf] = lid-1;

      if (face_dirs) {
        fdir0 = fdir0 ? 1 : -1;
        if (faceflip_[lid-1]) fdir0 *= -1;
        (*face_dirs)[nf] = fdir0;
      }
      nf++;

      /* If there is a last remaining face, it is the top face */
      MFace_ptr fopp;
      idx = 0;
      int i = 0;
      while ((fopp = List_Next_Entry(rfaces,&idx))) {
        if (fopp != face0 && !MEnt_IsMarked(fopp,mkid)) {
          lid = MEnt_ID(fopp);
          faceids[nf] = lid-1;

          if (face_dirs) {
            int fdir = (MR_FaceDir_i((MRegion_ptr)cell,i) == 1) ? 1 : -1;
            if (faceflip_[lid-1]) fdir *= -1;
            (*face_dirs)[nf] = fdir;
          }
          nf++;
          break;
        }
        ++i;
      }

      List_Unmark(rfaces,mkid);
      MSTK_FreeMarker(mkid);
      List_Delete(rfaces);

    } else {
      getCellFacesAndDirs_unordered_(cellid,faceids,face_dirs);
    }
  } else {
    getCellFacesAndDirs_unordered_(cellid,faceids,face_dirs);
  }
}


void Mesh_MSTK::getCellFacesAndDirs_unordered_(const Entity_ID cellid,
        Entity_ID_List& faceids,
        Entity_Direction_List * const face_dirs) const
{
  MEntity_ptr cell = cell_id_to_handle_[cellid];

  if (get_manifold_dimension() == 3) {
    int nrf;
    List_ptr rfaces;
    rfaces = MR_Faces((MRegion_ptr)cell);
    nrf = List_Num_Entries(rfaces);
    faceids.resize(nrf);

    Entity_ID_List::iterator itf = faceids.begin();
    for (int i = 0; i < nrf; ++i) {
      MFace_ptr face = List_Entry(rfaces,i);
      int lid = MEnt_ID(face);
      *itf = lid-1;  // assign to next spot by dereferencing iterator
      ++itf;
    }

    List_Delete(rfaces);
    if (face_dirs) {
      face_dirs->resize(nrf);
      std::vector<int>::iterator itd = face_dirs->begin();
      for (int i = 0; i < nrf; ++i) {
        int lid = faceids[i];
        int fdir = 2*MR_FaceDir_i((MRegion_ptr)cell,i) - 1;
        fdir = faceflip_[lid] ? -fdir : fdir;
        *itd = fdir;  // assign to next spot by dereferencing iterator
        ++itd;
      }
    }

  } else {  // get_manifold_dimension() = 2; surface or 2D mesh
    int nfe;
    List_ptr fedges;
    fedges = MF_Edges((MFace_ptr)cell,1,0);
    nfe = List_Num_Entries(fedges);

    faceids.resize(nfe);

    Entity_ID_List::iterator itf = faceids.begin();
    for (int i = 0; i < nfe; ++i) {
      MEdge_ptr edge = List_Entry(fedges,i);
      int lid = MEnt_ID(edge);
      *itf = lid-1;  // assign to next spot by dereferencing iterator
      ++itf;
    }

    List_Delete(fedges);
    if (face_dirs) {
      face_dirs->resize(nfe);
      std::vector<int>::iterator itd = face_dirs->begin();
      for (int i = 0; i < nfe; ++i) {
        int lid = faceids[i];
        int fdir = 2*MF_EdgeDir_i((MFace_ptr)cell,i) - 1;
        fdir = faceflip_[lid] ? -fdir : fdir;
        *itd = fdir;  // assign to next spot by dereferencing iterator
        itd++;
      }
    }
  }
}


void Mesh_MSTK::getCellFacesAndDirs(const Entity_ID cellid,
        Entity_ID_List& faceids,
        Entity_Direction_List * const face_dirs) const
{
  AMANZI_ASSERT(faces_initialized_);
  if (cells_initialized_) {
    getCellFacesAndDirs_ordered_(cellid, faceids, face_dirs);
  } else {
    getCellFacesAndDirs_unordered_(cellid, faceids, face_dirs);
  }
}


void Mesh_MSTK::getCellEdges(const Entity_ID cellid,
        Entity_ID_List& edgeids) const
{
  AMANZI_ASSERT(edges_initialized_);
  MEntity_ptr cell;
  cell = cell_id_to_handle_[cellid];

  if (get_manifold_dimension() == 3) {
    int nre;
    List_ptr redges;

    redges = MR_Edges((MRegion_ptr)cell);
    nre = List_Num_Entries(redges);
    edgeids.resize(nre);

    Entity_ID_List::iterator ite = edgeids.begin();
    for (int i = 0; i < nre; ++i) {
      MEdge_ptr edge = List_Entry(redges,i);
      int lid = MEnt_ID(edge);
      *ite = lid-1;  // assign to next spot by dereferencing iterator
      ++ite;
    }

    List_Delete(redges);

  } else {  // get_manifold_dimension() = 2; surface or 2D mesh
    int nfe;

    List_ptr fedges;
    fedges = MF_Edges((MFace_ptr)cell,1,0);
    nfe = List_Num_Entries(fedges);

    edgeids.resize(nfe);

    Entity_ID_List::iterator ite = edgeids.begin();
    for (int i = 0; i < nfe; ++i) {
      MEdge_ptr edge = List_Entry(fedges,i);
      int lid = MEnt_ID(edge);
      *ite = lid-1;  // assign to next spot by dereferencing iterator
      ++ite;
    }

    List_Delete(fedges);
  }
}


//---------------------------------------------------------
// Get nodes of cell
// On a distributed mesh, all nodes (OWNED or GHOST) of the cell
// are returned
// Nodes are returned in a standard order (Exodus II convention)
// STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES in 3D
// For a general polyhedron this will return the nodes in
// arbitrary order
// In 2D, the nodes of the polygon will be returned in ccw order
// consistent with the face normal
//---------------------------------------------------------
void Mesh_MSTK::getCellNodes(const Entity_ID cellid,
                               Entity_ID_List& nodeids) const
{
  int nn, lid;
  MEntity_ptr cell = cell_id_to_handle_[cellid];

  if (get_manifold_dimension() == 3) {                    // Volume mesh
    List_ptr rverts = MR_Vertices(cell);

    nn = List_Num_Entries(rverts);
    nodeids.resize(nn);
    Entity_ID_List::iterator it = nodeids.begin();
    for (int i = 0; i < nn; ++i) {
      lid = MEnt_ID(List_Entry(rverts,i));
      *it = lid-1;  // assign to next spot by dereferencing iterator
      ++it;
    }
    List_Delete(rverts);

  } else {                                 // Surface mesh
    List_ptr fverts = MF_Vertices(cell,1,0);
    nn = List_Num_Entries(fverts);
    nodeids.resize(nn);
    Entity_ID_List::iterator it = nodeids.begin();
    for (int i = 0; i < nn; ++i) {
      lid = MEnt_ID(List_Entry(fverts,i));
      *it = lid-1;  // assign to next spot by dereferencing iterator
      it++;
    }
    List_Delete(fverts);
  }
}


void Mesh_MSTK::getFaceEdgesAndDirs(const Entity_ID faceid,
                                                  Entity_ID_List& edgeids,
                                                  Entity_Direction_List * const edge_dirs) const
{
  AMANZI_ASSERT(faces_initialized_);
  AMANZI_ASSERT(edges_initialized_);

  MEntity_ptr face = face_id_to_handle_[faceid];
  if (get_manifold_dimension() == 3) {
    int nfe;
    List_ptr fedges;

    fedges = MF_Edges((MFace_ptr)face,1,0);
    nfe = List_Num_Entries(fedges);
    edgeids.resize(nfe);

    Entity_ID_List::iterator ite = edgeids.begin();
    for (int i = 0; i < nfe; ++i) {
      MEdge_ptr edge = List_Entry(fedges,i);
      int lid = MEnt_ID(edge);
      *ite = lid-1;  // assign to next spot by dereferencing iterator
      ++ite;
    }

    List_Delete(fedges);

    if (edge_dirs) {
      edge_dirs->resize(nfe);

      std::vector<int>::iterator itd = edge_dirs->begin();
      for (int i = 0; i < nfe; ++i) {
        int lid = edgeids[i];
        int edir = 2*MF_EdgeDir_i((MFace_ptr)face,i) - 1;
        edir = edgeflip_[lid] ? -edir : edir;
        *itd = edir;  // assign to next spot by dereferencing iterator
        ++itd;
      }
    }

  } else {  // get_manifold_dimension() = 2; surface or 2D mesh
    // face is same dimension as edge; just return the edge with a
    // direction of 1
    MEdge_ptr edge = (MEdge_ptr) face;

    edgeids.resize(1);
    edgeids[0] = MEnt_ID(edge)-1;

    if (edge_dirs) {
      edge_dirs->resize(1);
      (*edge_dirs)[0] = 1;
    }
  }
}


//---------------------------------------------------------
// Get nodes of face
// On a distributed mesh, all nodes (OWNED or GHOST) of the face
// are returned
// In 3D, the nodes of the face are returned in ccw order consistent
// with the face normal
// In 2D, nfnodes is 2
//---------------------------------------------------------
void Mesh_MSTK::getFaceNodes(const Entity_ID faceid,
                           Entity_ID_List& nodeids) const
{
  AMANZI_ASSERT(faces_initialized_);
  int nn, lid;
  MEntity_ptr genface = face_id_to_handle_[faceid];

  if (get_manifold_dimension() == 3) {  // Volume mesh
    int dir = !faceflip_[faceid];

    List_ptr fverts = MF_Vertices(genface,dir,0);
    AMANZI_ASSERT(fverts != nullptr);

    nn = List_Num_Entries(fverts);
    nodeids.resize(nn);
    Entity_ID_List::iterator it = nodeids.begin();

    for (int i = 0; i < nn; ++i) {
      lid = MEnt_ID(List_Entry(fverts,i));
      *it = lid-1;  // assign to next spot by dereferencing iterator
      ++it;
    }

    List_Delete(fverts);

  } else {                // Surface mesh or 2D mesh
    nodeids.resize(2);
    if (faceflip_[faceid]) {
      nodeids[0] = MEnt_ID(ME_Vertex(genface,1))-1;
      nodeids[1] = MEnt_ID(ME_Vertex(genface,0))-1;
    } else {
      nodeids[0] = MEnt_ID(ME_Vertex(genface,0))-1;
      nodeids[1] = MEnt_ID(ME_Vertex(genface,1))-1;
    }
  }
}


//---------------------------------------------------------
// Get nodes of an edge
//---------------------------------------------------------
void Mesh_MSTK::getEdgeNodes(const Entity_ID edgeid,
                           Entity_ID_List& nodes) const
{
  AMANZI_ASSERT(edges_initialized_);
  MEdge_ptr edge = (MEdge_ptr) edge_id_to_handle_[edgeid];
  nodes.resize(2);
  if (edgeflip_[edgeid]) {
    nodes[0] = MEnt_ID(ME_Vertex(edge,1))-1;
    nodes[1] = MEnt_ID(ME_Vertex(edge,0))-1;
  } else {
    nodes[0] = MEnt_ID(ME_Vertex(edge,0))-1;
    nodes[1] = MEnt_ID(ME_Vertex(edge,1))-1;
  }
}


//---------------------------------------------------------
// Cells of type 'ptype' connected to a node. This routine uses
// push_back on or near the partition boundary since we cannot tell at
// the outset how many entries will be put into the list
//---------------------------------------------------------
void Mesh_MSTK::getNodeCells(const Entity_ID nodeid,
                           const Parallel_type ptype,
                           Entity_ID_List& cellids) const
{
  int idx, lid, nc;
  List_ptr cell_list;
  MEntity_ptr ment;

  MVertex_ptr mv = (MVertex_ptr) vtx_id_to_handle_[nodeid];
  // mesh vertex on a processor boundary may be connected to owned
  // and ghost cells. So depending on the requested cell type, we
  // may have to omit some entries
  if (get_manifold_dimension() == 3)
    cell_list = MV_Regions(mv);
  else
    cell_list = MV_Faces(mv);

  nc = List_Num_Entries(cell_list);
  cellids.resize(nc); // resize to maximum size possible
  Entity_ID_List::iterator it = cellids.begin();

  int n = 0;
  idx = 0;
  while ((ment = List_Next_Entry(cell_list,&idx))) {
    if (MEnt_PType(ment) == PGHOST) {
      if (ptype == Parallel_type::GHOST || ptype == Parallel_type::ALL) {
        lid = MEnt_ID(ment);
        *it = lid-1;  // assign to next spot by dereferencing iterator
        ++it;
        ++n;
      }
    }
    else {
      if (ptype == Parallel_type::OWNED || ptype == Parallel_type::ALL) {
        lid = MEnt_ID(ment);
        *it = lid-1;  // assign to next spot by dereferencing iterator
        ++it;
        ++n;
      }
    }
  }
  cellids.resize(n); // resize to the actual number of cells being returned
  List_Delete(cell_list);
}


//---------------------------------------------------------
// Faces of type 'ptype' connected to a node. This routine uses
// push_back on or near the partition boundary since we cannot tell at
// the outset how many entries will be put into the list
//---------------------------------------------------------
void Mesh_MSTK::getNodeFaces(const Entity_ID nodeid,
                           const Parallel_type ptype,
                           Entity_ID_List& faceids) const
{
  int idx, lid, n;
  List_ptr face_list;
  MEntity_ptr ment;

  AMANZI_ASSERT(faces_initialized_);
  MVertex_ptr mv = (MVertex_ptr) vtx_id_to_handle_[nodeid];

  if (get_manifold_dimension() == 3)
    face_list = MV_Faces(mv);
  else
    face_list = MV_Edges(mv);

  int nf = List_Num_Entries(face_list);
  faceids.resize(nf); // resize to the maximum
  Entity_ID_List::iterator it = faceids.begin();
  idx = 0; n = 0;
  while ((ment = List_Next_Entry(face_list,&idx))) {
    if (MEnt_PType(ment) == PGHOST) {
      if (ptype == Parallel_type::GHOST || ptype == Parallel_type::ALL) {
        lid = MEnt_ID(ment);
        *it = lid-1;  // assign to next spot by dereferencing iterator
        ++it;
        ++n;
      }
    }
    else {
      if (ptype == Parallel_type::OWNED || ptype == Parallel_type::ALL) {
        lid = MEnt_ID(ment);
        *it = lid-1;  // assign to next spot by dereferencing iterator
        ++it;
        ++n;
      }
    }
  }
  faceids.resize(n); // resize to the actual number of faces being returned

  List_Delete(face_list);
}


//---------------------------------------------------------
// Edges of type 'ptype' connected to a node.
//---------------------------------------------------------
void Mesh_MSTK::getNodeEdges(const Entity_ID nodeid,
                           const Parallel_type ptype,
                           Entity_ID_List& edgeids) const
{
  int idx, lid, nc;
  List_ptr edge_list;
  MEntity_ptr ment;

  MVertex_ptr mv = (MVertex_ptr) vtx_id_to_handle_[nodeid];

  // mesh vertex on a processor boundary may be connected to owned
  // and ghost cells. So depending on the requested cell type, we
  // may have to omit some entries
  edge_list = MV_Edges(mv);
  nc = List_Num_Entries(edge_list);

  edgeids.resize(nc); // resize to maximum size possible
  Entity_ID_List::iterator it = edgeids.begin();

  int n = 0;
  idx = 0;
  while ((ment = List_Next_Entry(edge_list, &idx))) {
    if (MEnt_PType(ment) == PGHOST) {
      if (ptype == Parallel_type::GHOST || ptype == Parallel_type::ALL) {
        lid = MEnt_ID(ment);
        *it = lid-1;  // assign to next spot by dereferencing iterator
        ++it;
        ++n;
      }
    }
    else {
      if (ptype == Parallel_type::OWNED || ptype == Parallel_type::ALL) {
        lid = MEnt_ID(ment);
        *it = lid-1;  // assign to next spot by dereferencing iterator
        ++it;
        ++n;
      }
    }
  }
  edgeids.resize(n); // resize to the actual number of cells being returned
  List_Delete(edge_list);
}


//---------------------------------------------------------
// Faces of type 'ptype' connected to an edge.
//---------------------------------------------------------
void Mesh_MSTK::getEdgeFaces(const Entity_ID edgeid,
                           const Parallel_type ptype,
                           Entity_ID_List& faceids) const
{
  int idx, lid, nc;
  List_ptr face_list;
  MEntity_ptr ment;

  AMANZI_ASSERT(get_manifold_dimension() == 3);

  MEdge_ptr me = (MEdge_ptr) edge_id_to_handle_[edgeid];
  face_list = ME_Faces(me);

  nc = List_Num_Entries(face_list);
  faceids.resize(nc); // resize to maximum size possible
  Entity_ID_List::iterator it = faceids.begin();

  int n = 0;
  idx = 0;
  while ((ment = List_Next_Entry(face_list,&idx))) {
    if (MEnt_PType(ment) == PGHOST) {
      if (ptype == Parallel_type::GHOST || ptype == Parallel_type::ALL) {
        lid = MEnt_ID(ment);
        *it = lid-1;  // assign to next spot by dereferencing iterator
        ++it;
        ++n;
      }

    } else {
      if (ptype == Parallel_type::OWNED || ptype == Parallel_type::ALL) {
        lid = MEnt_ID(ment);
        *it = lid-1;  // assign to next spot by dereferencing iterator
        ++it;
        ++n;
      }
    }
  }
  faceids.resize(n); // resize to the actual number of cells being returned
  List_Delete(face_list);
}


//---------------------------------------------------------
// Cells of type 'ptype' connected to an edge. This routine uses
// push_back on or near the partition boundary since we cannot tell at
// the outset how many entries will be put into the list
//---------------------------------------------------------
void Mesh_MSTK::getEdgeCells(const Entity_ID edgeid,
                           const Parallel_type ptype,
                           Entity_ID_List& cellids) const
{
  MEdge_ptr me = (MEdge_ptr) edge_id_to_handle_[edgeid];

  // mesh edge on a processor boundary may be connected to owned
  // and ghost cells. So depending on the requested cell type, we
  // may have to omit some entries
  List_ptr cell_list;
  if (get_manifold_dimension() == 3)
    cell_list = ME_Regions(me);
  else
    cell_list = ME_Faces(me);

  int nc = List_Num_Entries(cell_list);
  cellids.resize(nc); // resize to maximum size possible

  int n = 0;
  int idx = 0;
  MEntity_ptr ment;
  while ((ment = List_Next_Entry(cell_list,&idx))) {
    if (MEnt_PType(ment) == PGHOST) {
      if (ptype == Parallel_type::GHOST || ptype == Parallel_type::ALL) {
        int lid = MEnt_ID(ment);
        cellids[n] = lid-1;
        ++n;
      }

    } else {
      if (ptype == Parallel_type::OWNED || ptype == Parallel_type::ALL) {
        int lid = MEnt_ID(ment);
        cellids[n] = lid-1;
        ++n;
      }
    }
  }
  cellids.resize(n); // resize to the actual number of cells being returned
  List_Delete(cell_list);
}


//---------------------------------------------------------
// Cells connected to a face
//---------------------------------------------------------
void Mesh_MSTK::getFaceCells(const Entity_ID faceid,
        const Parallel_type ptype,
        Entity_ID_List& cellids) const
{
  AMANZI_ASSERT(faces_initialized_);
  cellids.clear();

  if (get_manifold_dimension() == 3) {
    MFace_ptr mf = (MFace_ptr) face_id_to_handle_[faceid];

    List_ptr fregs = MF_Regions(mf);
    MRegion_ptr mr;
    if (ptype == Parallel_type::ALL) {
      int idx = 0;
      while ((mr = List_Next_Entry(fregs,&idx)))
        cellids.push_back(MR_ID(mr)-1);

    } else {
      int idx = 0;
      while ((mr = List_Next_Entry(fregs,&idx))) {
        if (MEnt_PType(mr) == PGHOST) {
          if (ptype == Parallel_type::GHOST)
            cellids.push_back(MR_ID(mr)-1);
        } else if (ptype == Parallel_type::OWNED) {
            cellids.push_back(MR_ID(mr)-1);
        }
      }
    }
    List_Delete(fregs);

  } else {
    MEdge_ptr me = (MEdge_ptr) face_id_to_handle_[faceid];

    List_ptr efaces = ME_Faces(me);
    MFace_ptr mf;
    if (ptype == Parallel_type::ALL) {
      int idx = 0;
      while ((mf = List_Next_Entry(efaces,&idx)))
        cellids.push_back(MF_ID(mf)-1);
    }
    else {
      int idx = 0;
      while ((mf = List_Next_Entry(efaces,&idx))) {
        if (MEnt_PType(mf) == PGHOST) {
          if (ptype == Parallel_type::GHOST)
            cellids.push_back(MF_ID(mf)-1);
        } else if (ptype == Parallel_type::OWNED) {
          cellids.push_back(MF_ID(mf)-1);
        }
      }
    }
    List_Delete(efaces);
  }
}


//---------------------------------------------------------
// Node coordinates - 3 in 3D and 2 in 2D
//---------------------------------------------------------
AmanziGeometry::Point
Mesh_MSTK::getNodeCoordinate(const Entity_ID nodeid) const
{
  MEntity_ptr vtx = vtx_id_to_handle_[nodeid];
  double coords[3];
  MV_Coords(vtx,coords);
  if (get_space_dimension() == 2) {
    return AmanziGeometry::Point(coords[0], coords[1]);
  } else {
    return AmanziGeometry::Point(coords[0], coords[1], coords[2]);
  }
}


//---------------------------------------------------------
// Modify a node's coordinates
//---------------------------------------------------------
void Mesh_MSTK::setNodeCoordinate(const AmanziMesh::Entity_ID nodeid,
                                     const AmanziGeometry::Point& coords)
{
  MVertex_ptr v = vtx_id_to_handle_[nodeid];
  double coordarray[3] = {0.0,0.0,0.0};
  for (int i = 0; i < get_space_dimension(); i++)
    coordarray[i] = coords[i];
  MV_Set_Coords(v,(double *)coordarray);
}


//---------------------------------------------------------
// Parent entity in the source mesh if mesh was derived from another mesh
//---------------------------------------------------------
Entity_ID
Mesh_MSTK::getEntityParent(const Entity_kind kind, const Entity_ID entid) const
{
  int ival;
  double rval;
  void *pval = nullptr;
  MEntity_ptr ment = nullptr;
  MAttrib_ptr att = nullptr;

  switch(kind) {
  case Entity_kind::CELL:
    att = (get_manifold_dimension() == 3) ? rparentatt_ : fparentatt_;
    ment = (MEntity_ptr) cell_id_to_handle_[entid];
    break;
  case Entity_kind::FACE:
    att = (get_manifold_dimension() == 3) ? fparentatt_ : eparentatt_;
    ment = (MEntity_ptr) face_id_to_handle_[entid];
    break;
  case Entity_kind::EDGE:
    att = eparentatt_;
    ment = (MEntity_ptr) edge_id_to_handle_[entid];
    break;
  case Entity_kind::NODE:
    if (!vparentatt_) return 0;
    att = vparentatt_;
    ment = (MEntity_ptr) vtx_id_to_handle_[entid];
    break;
  default:
    {}
  }
  if (!att) return 0;
  MEnt_Get_AttVal(ment,att,&ival,&rval,&pval);
  if (pval)
    return MEnt_ID((MEntity_ptr)pval)-1;
  else
    return 0;
}


//---------------------------------------------------------
// Global ID of any entity
//---------------------------------------------------------
Entity_GID Mesh_MSTK::getEntityGID(const Entity_kind kind, const Entity_ID lid) const
{
  MEntity_ptr ent;
  switch (kind) {
  case Entity_kind::NODE:
    ent = vtx_id_to_handle_[lid];
    break;
  case Entity_kind::EDGE:
    ent = edge_id_to_handle_[lid];
    break;
  case Entity_kind::FACE:
    ent = face_id_to_handle_[lid];
    break;
  case Entity_kind::CELL:
    ent = cell_id_to_handle_[lid];
    break;
  default:
    std::cerr << "Global ID requested for unknown entity type" << std::endl;
  }

  if (serial_run) return MEnt_ID(ent)-1;
  else return MEnt_GlobalID(ent)-1;
}


//---------------------------------------------------------
// Procedure to perform all the post-mesh creation steps in a constructor
//---------------------------------------------------------
void Mesh_MSTK::post_create_steps_()
{
  // Initialize data structures for various entities - vertices/nodes
  // and cells are always initialized; edges and faces only if
  // requested
  init_nodes_();
  edgeflip_ = nullptr;
  if (edges_requested_) init_edges_();
  if (faces_requested_) init_faces_();
  init_cells_();
}


//---------------------------------------------------------
// Some initializations
//---------------------------------------------------------
void Mesh_MSTK::clear_internals_()
{
  faceflip_ = nullptr;
  mesh_ = nullptr;
  owned_verts_ = ghost_verts_ = nullptr;
  owned_edges_ = ghost_edges_ = nullptr;
  owned_faces_ = ghost_faces_ = nullptr;
  owned_cells_ = ghost_cells_ = nullptr;
  celltype_att_ = nullptr;
  rparentatt_ = fparentatt_ = eparentatt_ = vparentatt_ = nullptr;

  edges_initialized_ = false;
  faces_initialized_ = false;
  owned_verts_ = ghost_verts_ = nullptr;
  owned_edges_ = ghost_edges_ = nullptr;
  owned_faces_ = ghost_faces_ = nullptr;
  owned_cells_ = ghost_cells_ = nullptr;
  deleted_vertices_ = deleted_edges_ = deleted_faces_ = deleted_regions_ = nullptr;
  entities_deleted_ = false;
}


//---------------------------------------------------------
// initialize vertex info
//---------------------------------------------------------
void Mesh_MSTK::init_nodes_()
{
  // create owned and not owned vertex lists
  init_pvert_lists_();
  // create maps from IDs to handles
  init_vertex_id2handle_maps_();
}


//---------------------------------------------------------
// Initialize edge info
//---------------------------------------------------------
void Mesh_MSTK::init_edges_()
{
  edges_initialized_ = true;
  // Create owned and not owned lists
  init_pedge_lists_();
  // Create maps from IDs to handles
  init_edge_id2handle_maps_();
  // Initialize boolean flag indicating whether slave edges are reversed in
  // direction from the master and must be flipped
  init_pedge_dirs_();
}


//---------------------------------------------------------
// Initialize face info
//---------------------------------------------------------
void Mesh_MSTK::init_faces_()
{
  faces_initialized_ = true;
  // Create owned and not owned lists
  init_pface_lists_();
  // Create maps from IDs to handles
  init_face_id2handle_maps_();
  // Initialize boolean flag indicating whether slave faces are reversed in
  // direction from the master and must be flipped
  init_pface_dirs_();
}


//---------------------------------------------------------
// Initialize cell info
//---------------------------------------------------------
void Mesh_MSTK::init_cells_()
{
  label_celltype_();
  // create owned and not owned cell lists
  init_pcell_lists_();
  // create maps from IDs to handles
  init_cell_id2handle_maps_();
  cells_initialized_ = true;
}


//---------------------------------------------------------
// ID to handle/pointer map for vertices
//---------------------------------------------------------
void Mesh_MSTK::init_vertex_id2handle_maps_()
{
  // If the mesh is dynamic, then this code has to be revisited
  // Amanzi has IDs starting from 0, MSTK has IDs starting from 1
  int nv = MESH_Num_Vertices(mesh_);
  vtx_id_to_handle_.resize(nv);

  int idx = 0;
  int lid = 1;
  MVertex_ptr vtx;
  while ((vtx = MSet_Next_Entry(owned_verts_,&idx))) {
    MEnt_Set_ID(vtx,lid);
    vtx_id_to_handle_[lid-1] = vtx;
    lid++;
  }

  idx = 0;
  while ((vtx = MSet_Next_Entry(ghost_verts_,&idx))) {
    MEnt_Set_ID(vtx,lid);
    vtx_id_to_handle_[lid-1] = vtx;
    lid++;
  }
}


//---------------------------------------------------------
// ID to handle/pointer map for edges
//---------------------------------------------------------
void Mesh_MSTK::init_edge_id2handle_maps_()
{
  // If the mesh is dynamic, then this code has to be revisited
  // Amanzi has IDs starting from 0, MSTK has IDs starting from 1
  int ne = MESH_Num_Edges(mesh_);
  edge_id_to_handle_.resize(ne);

  int idx = 0;
  int lid = 1;
  MEdge_ptr edge;
  while ((edge = MSet_Next_Entry(owned_edges_,&idx))) {
    MEnt_Set_ID(edge,lid);
    edge_id_to_handle_[lid-1] = edge;
    lid++;
  }

  idx = 0;
  while ((edge = MSet_Next_Entry(ghost_edges_,&idx))) {
    MEnt_Set_ID(edge,lid);
    edge_id_to_handle_[lid-1] = edge;
    lid++;
  }
}


//---------------------------------------------------------
// ID to handle/pointer map for faces
//---------------------------------------------------------
void Mesh_MSTK::init_face_id2handle_maps_()
{
  // If the mesh is dynamic, then this code has to be revisited
  // Amanzi has IDs starting from 0, MSTK has IDs starting from 1
  int nf = (get_manifold_dimension() == 2) ? MESH_Num_Edges(mesh_) : MESH_Num_Faces(mesh_);
  face_id_to_handle_.resize(nf);

  int idx = 0;
  int lid = 1;
  MEntity_ptr genface;
  while ((genface = MSet_Next_Entry(owned_faces_,&idx))) {
    MEnt_Set_ID(genface,lid);
    face_id_to_handle_[lid-1] = genface;
    lid++;
  }

  idx = 0;
  while ((genface = MSet_Next_Entry(ghost_faces_,&idx))) {
    MEnt_Set_ID(genface,lid);
    face_id_to_handle_[lid-1] = genface;
    lid++;
  }
}


//---------------------------------------------------------
// ID to handle/pointer map for cells
//---------------------------------------------------------
void Mesh_MSTK::init_cell_id2handle_maps_()
{
  // If the mesh is dynamic, then this code has to be revisited
  // Amanzi has IDs starting from 0, MSTK has IDs starting from 1
  int nc = (get_manifold_dimension() == 2) ? MESH_Num_Faces(mesh_) : MESH_Num_Regions(mesh_);
  cell_id_to_handle_.resize(nc);

  int idx = 0;
  int lid = 1;
  MEntity_ptr gencell;  // Mesh region in 3D, face in 2D
  while ((gencell = MSet_Next_Entry(owned_cells_,&idx))) {
    MEnt_Set_ID(gencell,lid);
    cell_id_to_handle_[lid-1] = gencell;
    lid++;
  }

  idx = 0;
  while ((gencell = MSet_Next_Entry(ghost_cells_,&idx))) {
    MEnt_Set_ID(gencell,lid);
    cell_id_to_handle_[lid-1] = gencell;
    lid++;
  }
}


//---------------------------------------------------------
// create lists of owned and not owned vertices
//---------------------------------------------------------
void Mesh_MSTK::init_pvert_lists_()
{
  // Get all vertices on this processor
  ghost_verts_ = MSet_New(mesh_,"ghost_verts_",MVERTEX);
  owned_verts_ = MSet_New(mesh_,"owned_verts_",MVERTEX);

  int idx = 0;
  MVertex_ptr vtx;
  while ((vtx = MESH_Next_Vertex(mesh_,&idx))) {
    if (MV_PType(vtx) == PGHOST)
      MSet_Add(ghost_verts_,vtx);
    else
      MSet_Add(owned_verts_,vtx);
  }
}


//---------------------------------------------------------
// create lists of owned and not owned edges
//---------------------------------------------------------
void Mesh_MSTK::init_pedge_lists_()
{
  // Get all vertices on this processor
  ghost_edges_ = MSet_New(mesh_,"ghost_edges_",MEDGE);
  owned_edges_ = MSet_New(mesh_,"owned_edges_",MEDGE);

  int idx = 0;
  MEdge_ptr edge;
  while ((edge = MESH_Next_Edge(mesh_,&idx))) {
    if (ME_PType(edge) == PGHOST)
      MSet_Add(ghost_edges_,edge);
    else
      MSet_Add(owned_edges_,edge);
  }
} // Mesh_MSTK::init_pedge_lists_


void Mesh_MSTK::init_pedge_dirs_()
{
  int ne = MESH_Num_Edges(mesh_);

  if (serial_run) {
    edgeflip_ = new bool[ne];
    for (int i = 0; i < ne; ++i) edgeflip_[i] = false;

  } else {
    // Do some additional processing to see if ghost edges and their masters
    // are oriented the same way; if not, turn on flag to flip the directions
    // when returning to the application code
    MAttrib_ptr attev0 = MAttrib_New(mesh_,"TMP_EV0_ATT",INT,MEDGE);
    MAttrib_ptr attev1 = MAttrib_New(mesh_,"TMP_EV1_ATT",INT,MEDGE);

    int idx = 0;
    MEdge_ptr edge;
    while ((edge = MESH_Next_Edge(mesh_,&idx))) {
      if (ME_PType(edge) != PINTERIOR) {
        MVertex_ptr vertex0 = ME_Vertex(edge,0);
        MVertex_ptr vertex1 = ME_Vertex(edge,1);

        MEnt_Set_AttVal(edge,attev0,MEnt_GlobalID(vertex0),0.0,NULL);
        MEnt_Set_AttVal(edge,attev1,MEnt_GlobalID(vertex1),0.0,NULL);
      }
    }

    MESH_UpdateAttributes(mesh_,mpicomm_);

    edgeflip_ = new bool[ne];
    for (int i = 0; i < ne; ++i) edgeflip_[i] = false;

    double rval;
    void *pval;
    idx = 0;
    while ((edge = MSet_Next_Entry(ghost_edges_,&idx))) {
      int remote_vertexid0, remote_vertexid1;

      MEnt_Get_AttVal(edge,attev0,&remote_vertexid0,&rval,&pval);
      MEnt_Get_AttVal(edge,attev1,&remote_vertexid1,&rval,&pval);

      int local_vertexid0 = MEnt_GlobalID(ME_Vertex(edge,0));
      int local_vertexid1 = MEnt_GlobalID(ME_Vertex(edge,1));

      if (remote_vertexid1 == local_vertexid0 ||
          remote_vertexid0 == local_vertexid1) {
        int lid = MEnt_ID(edge);
        edgeflip_[lid-1] = true;

      } else { // Sanity Check
        if (remote_vertexid1 != local_vertexid1 &&
            remote_vertexid0 != local_vertexid0) {
          std::stringstream mesg_stream;
          mesg_stream << "Edge vertices mismatch between master and ghost (processor " << myprocid << ")";
          Errors::Message mesg(mesg_stream.str());
          Exceptions::amanzi_throw(mesg);
        }
      }
    }
  }
}


//---------------------------------------------------------
// Create lists of owned and not owned faces
//---------------------------------------------------------
void Mesh_MSTK::init_pface_lists_()
{
  // Get all faces on this processor
  if (get_manifold_dimension() == 3) {
    ghost_faces_ = MSet_New(mesh_,"ghost_faces_",MFACE);
    owned_faces_ = MSet_New(mesh_,"owned_faces_",MFACE);

    int idx = 0;
    MFace_ptr face;
    while ((face = MESH_Next_Face(mesh_,&idx))) {
      if (MF_PType(face) == PGHOST)
        MSet_Add(ghost_faces_,face);
      else
        MSet_Add(owned_faces_,face);
    }

  } else if (get_manifold_dimension() == 2) {
    ghost_faces_ = MSet_New(mesh_,"ghost_faces_",MFACE);
    owned_faces_ = MSet_New(mesh_,"owned_faces_",MFACE);

    int idx = 0;
    MEdge_ptr edge;
    while ((edge = MESH_Next_Edge(mesh_,&idx))) {
      if (ME_PType(edge) == PGHOST)
        MSet_Add(ghost_faces_,edge);
      else
        MSet_Add(owned_faces_,edge);
    }
  }
  return;
}

// Detect whether ghost faces are in opposite direction of owned faces
// on processor boundaries
void Mesh_MSTK::init_pface_dirs_()
{
  int nf = (get_manifold_dimension() == 2) ? MESH_Num_Edges(mesh_) : MESH_Num_Faces(mesh_);

  if (serial_run) {
    faceflip_ = new bool[nf];
    for (int i = 0; i < nf; ++i) faceflip_[i] = false;

  } else {
    if (get_manifold_dimension() == 3)
      init_pface_dirs_3_();
    else if (get_manifold_dimension() == 2)
      init_pface_dirs_2_();
  }
}


// Detect whether ghost faces are in opposite direction of owned faces
// on processor boundaries - Version for solid meshes
void Mesh_MSTK::init_pface_dirs_3_()
{
  int nf = MESH_Num_Faces(mesh_);

  // Do some additional processing to see if ghost faces and their masters
  // are oriented the same way; if not, turn on flag to flip the directions
  // when returning to the application code

  // attributes to store
  MAttrib_ptr attfc0 = MAttrib_New(mesh_,"TMP_FC0_ATT",INT,MFACE);
  MAttrib_ptr attfc1 = MAttrib_New(mesh_,"TMP_FC1_ATT",INT,MFACE);

  int idx = 0;
  MFace_ptr face;
  while ((face = MESH_Next_Face(mesh_,&idx))) {
    if (MF_PType(face) != PINTERIOR) {
      MRegion_ptr region0 = MF_Region(face,0);
      if (region0)
        MEnt_Set_AttVal(face,attfc0,MEnt_GlobalID(region0),0.0,NULL);

      MRegion_ptr region1 = MF_Region(face,1);
      if (region1)
        MEnt_Set_AttVal(face,attfc1,MEnt_GlobalID(region1),0.0,NULL);
    }
  }

  MESH_UpdateAttributes(mesh_,mpicomm_);

  faceflip_ = new bool[nf];
  for (int i = 0; i < nf; ++i) faceflip_[i] = false;

  idx = 0;
  while ((face = MSet_Next_Entry(ghost_faces_,&idx))) {
    int remote_regid0, remote_regid1;
    double rval;
    void *pval;
    MEnt_Get_AttVal(face,attfc0,&remote_regid0,&rval,&pval);
    MEnt_Get_AttVal(face,attfc1,&remote_regid1,&rval,&pval);

    MRegion_ptr region0 = MF_Region(face,0);
    int local_regid0 = region0 ? MEnt_GlobalID(region0) : 0;
    MRegion_ptr region1 = MF_Region(face,1);
    int local_regid1 = region1 ? MEnt_GlobalID(region1) : 0;

    if (remote_regid1 == local_regid0 ||
        remote_regid0 == local_regid1) {
      int lid = MEnt_ID(face);
      faceflip_[lid-1] = true;

    } else { // Sanity Check
      if (remote_regid1 != local_regid1 &&
          remote_regid0 != local_regid0) {
        Errors::Message mesg;
        mesg << "Face cells mismatch between master and ghost (processor " << myprocid << ")";
        Exceptions::amanzi_throw(mesg);
      }
    }
  }
  MAttrib_Delete(attfc0);
  MAttrib_Delete(attfc1);
}


// Detect whether ghost faces are in opposite direction of owned faces
// on processor boundaries - Version for surface meshes
void Mesh_MSTK::init_pface_dirs_2_()
{
  int ne = MESH_Num_Edges(mesh_);

  // Do some additional processing to see if ghost faces and their masters
  // are oriented the same way; if not, turn on flag to flip the directions
  // when returning to the application code
  MAttrib_ptr attev0 = MAttrib_New(mesh_,"TMP_EV0_ATT",INT,MEDGE);
  MAttrib_ptr attev1 = MAttrib_New(mesh_,"TMP_EV1_ATT",INT,MEDGE);

  int idx = 0;
  MEdge_ptr edge;
  while ((edge = MESH_Next_Edge(mesh_,&idx))) {
    if (ME_PType(edge) != PINTERIOR) {
      MVertex_ptr ev0 = ME_Vertex(edge, 0);
      MVertex_ptr ev1 = ME_Vertex(edge, 1);

      MEnt_Set_AttVal(edge,attev0,MEnt_GlobalID(ev0),0.0,NULL);
      MEnt_Set_AttVal(edge,attev1,MEnt_GlobalID(ev1),0.0,NULL);
    }
  }

  MESH_UpdateAttributes(mesh_,mpicomm_);

  faceflip_ = new bool[ne];
  for (int i = 0; i < ne; ++i) faceflip_[i] = false;

  idx = 0;
  while ((edge = MSet_Next_Entry(ghost_faces_,&idx))) {
    int remote_evgid0, remote_evgid1;
    double rval;
    void *pval;
    MEnt_Get_AttVal(edge,attev0,&remote_evgid0,&rval,&pval);
    MEnt_Get_AttVal(edge,attev1,&remote_evgid1,&rval,&pval);

    MVertex_ptr ev0 = ME_Vertex(edge, 0);
    MVertex_ptr ev1 = ME_Vertex(edge, 1);
    int local_evgid0 = MV_GlobalID(ev0);
    int local_evgid1 = MV_GlobalID(ev1);

    if (remote_evgid1 == local_evgid0 ||
        remote_evgid0 == local_evgid1) {
      int lid = MEnt_ID(edge);
      faceflip_[lid-1] = true;
    }
  }
  MAttrib_Delete(attev0);
  MAttrib_Delete(attev1);
}


//---------------------------------------------------------
// create lists of owned and not owned cells
//---------------------------------------------------------
void Mesh_MSTK::init_pcell_lists_()
{
  int idx = 0;
  if (get_manifold_dimension() == 3) {
    MRegion_ptr region;
    owned_cells_ = MSet_New(mesh_,"owned_cells_",MREGION);
    ghost_cells_ = MSet_New(mesh_,"ghost_cells_",MREGION);

    idx = 0;
    while ((region = MESH_Next_Region(mesh_,&idx))) {
      if (MR_PType(region) == PGHOST)
        MSet_Add(ghost_cells_,region);
      else
        MSet_Add(owned_cells_,region);
    }

  } else if (get_manifold_dimension() == 2) {
    MFace_ptr face;
    owned_cells_ = MSet_New(mesh_,"owned_cells_",MFACE);
    ghost_cells_ = MSet_New(mesh_,"ghost_cells_",MFACE);

    idx = 0;
    while ((face = MESH_Next_Face(mesh_,&idx))) {
      if (MF_PType(face) == PGHOST)
        MSet_Add(ghost_cells_,face);
      else
        MSet_Add(owned_cells_,face);
    }

  } else {
    Errors::Message mesg("Implemented only for 2D and 3D");
    Exceptions::amanzi_throw(mesg);
  }
  return;
}

void Mesh_MSTK::label_celltype_()
{
  if (get_manifold_dimension() == 2) {
    celltype_att_ = MAttrib_New(mesh_,"Cell_type",INT,MFACE);
    int idx = 0;
    MFace_ptr face;
    while ((face = MESH_Next_Face(mesh_,&idx))) {
      Entity_ID_List edges;
      List_ptr fedges = MF_Edges(face, 0, 0);
      int idx2 = 0;
      MEdge_ptr edge;
      while ((edge = List_Next_Entry(fedges,&idx2)))
        edges.push_back(MEnt_ID(edge)-1);

      Cell_type ctype = MeshFramework::getCellType_(MEnt_ID(face), edges);
      MEnt_Set_AttVal(face, celltype_att_, (int) ctype, 0.0, NULL);
    }

  } else if (get_manifold_dimension() == 3) {
    celltype_att_ = MAttrib_New(mesh_,"Cell_type",INT,MREGION);
    int idx = 0;
    MRegion_ptr region;
    while ((region = MESH_Next_Region(mesh_,&idx))) {
      List_ptr rfaces = MR_Faces(region);
      Entity_ID_List faces;
      int idx2 = 0;
      MFace_ptr face;
      while ((face = List_Next_Entry(rfaces,&idx2)))
        faces.push_back(MEnt_ID(face)-1);

      Cell_type ctype = MeshFramework::getCellType_(MEnt_ID(region), faces);
      MEnt_Set_AttVal(region, celltype_att_, (int) ctype, 0.0, NULL);
    }
  }
}



void
Mesh_MSTK::getSetEntities(const AmanziGeometry::RegionLabeledSet& region,
                            const Entity_kind kind,
                            const Parallel_type ptype,
                            Entity_ID_List& entids) const
{
  if (kind != createEntityKind(region.entity_str())) {
    Errors::Message msg;
    msg << "Inconsistent request of labeled set for region \"" << region.name()
        << "\" of kind \"" << region.entity_str() << "\" requested as type \""
        << to_string(kind) << "\"";
    Exceptions::amanzi_throw(msg);
  }
  std::string internal_name = internal_name_of_set_(region, kind);
  auto mset = MESH_MSetByName(mesh_, internal_name.c_str());

  // check for the wierd case of both Material and CellSet
  if (kind == Entity_kind::CELL) {
    std::string other_internal_name = other_internal_name_of_set_(region, kind);
    auto mset2 = MESH_MSetByName(mesh_, other_internal_name.c_str());

    if (mset) {
      if (mset2) {
        Errors::Message msg;
        msg << "Exodus II file has both element block and element set with ID "
            << region.label() << " - Amanzi cannot handle this case.";
        Exceptions::amanzi_throw(msg);
      }

    } else {
      if (mset2) {
        mset = mset2;
      } else {
        // from the old code, it is unclear whether this is right.  It seems
        // possible that both are empty because no entities on this partition.
        //
        // But why wouldn't this have errored previously then? We need a
        // parallel test with a mesh such that only one cell is labeled with a
        // given material ID to test. --etc
        Errors::Message msg;
        msg << "Exodus II file has no labeled cell set with ID " << region.label();
        Exceptions::amanzi_throw(msg);
      }
    }
  }

  if (mset) {
    // Its possible some sets won't exist on some partitions
    int entdim = MSet_EntDim(mset);
    if (get_manifold_dimension() == 3) {
      if ((region.entity_str() == "CELL" && entdim != MREGION) ||
          (region.entity_str() == "FACE" && entdim != MFACE) ||
          (region.entity_str() == "NODE" && entdim != MVERTEX)) {
        Errors::Message mesg("Mismatch of entity type in labeled set region and mesh set");
        Exceptions::amanzi_throw(mesg);
      }

    } else if (get_manifold_dimension() == 2) {
      if ((region.entity_str() == "CELL" && entdim != MFACE) ||
          (region.entity_str() == "FACE" && entdim != MEDGE) ||
          (region.entity_str() == "NODE" && entdim != MVERTEX)) {
        Errors::Message msg("Mismatch of entity type in labeled set region and mesh set");
        Exceptions::amanzi_throw(msg);
      }
    }

    if (entities_deleted_) {
      int idx = 0;
      MEntity_ptr ent;
      while ((ent = MSet_Next_Entry(mset,&idx))) {
        if (MEnt_Dim(ent) == MDELETED)
          MSet_Rem(mset, ent);
      }
    }

    int nent_loc = MSet_Num_Entries(mset);
    entids.resize(nent_loc);

    if (nent_loc) {
      nent_loc = 0; // reset and count to get the real number
      int idx = 0;
      MEntity_ptr ment;
      switch (ptype) {
        case Parallel_type::OWNED:
          idx = 0;
          while ((ment = MSet_Next_Entry(mset,&idx))) {
            if (MEnt_PType(ment) != PGHOST) {
              entids[nent_loc] = MEnt_ID(ment)-1;
              ++nent_loc;
            }
          }
          break;
        case Parallel_type::GHOST:
          idx = 0;
          while ((ment = MSet_Next_Entry(mset,&idx))) {
            if (MEnt_PType(ment) == PGHOST) {
              entids[nent_loc] = MEnt_ID(ment)-1;
              ++nent_loc;
            }
          }
          break;
        case Parallel_type::ALL:
          idx = 0;
          while ((ment = MSet_Next_Entry(mset,&idx))) {
            entids[nent_loc] = MEnt_ID(ment)-1;
            ++nent_loc;
          }
          break;
        default:
        {}
      }
      entids.resize(nent_loc);
    }

  } else {
    entids.resize(0);
  }
}


void Mesh_MSTK::collapse_degen_edges_()
{
  const int topoflag=0; // Don't worry about violation of model classification
  int idx, evgid0, evgid1;
  MVertex_ptr vertex, ev0, ev1, vkeep, vdel;
  MEdge_ptr edge;
  MFace_ptr face;
  List_ptr deleted_ents_all = List_New(10);
  List_ptr merged_entity_pairs_all = List_New(10);
  double len2;
  std::vector<int> merged_ents_info;

  idx = 0;
  while ((edge = MESH_Next_Edge(mesh_,&idx))) {

    len2 = ME_Len(edge);

    if (len2 <= 1.0e-15) {
      /* Degenerate edge  - must collapse */
      entities_deleted_ = true;

      /* Collapse, choosing the vertex to be deleted and vertex to
         be kept consistently. If topological constraints permit,
         collapse the vertex with the higher global ID to the vertex
         with the lower global ID. If they do not, reverse the
         order. Since global IDs and topological constraints are the
         same for master and slave edges and their nodes, we will not
         have conflict between processors */
      ev0 = ME_Vertex(edge,0); evgid0 = MEnt_GlobalID(ev0);
      ev1 = ME_Vertex(edge,1); evgid1 = MEnt_GlobalID(ev1);

      if (evgid0 < evgid1) {
        vkeep = ev0;
        vdel = ev1;
      }
      else {
        vkeep = ev1;
        vdel = ev0;
      }

      List_ptr deleted_ents = NULL, merged_entity_pairs = NULL;
      vkeep = ME_Collapse(edge, vkeep, topoflag, &deleted_ents,
                          &merged_entity_pairs);

      if (!vkeep) {
        vkeep = vdel;
        vdel = (vkeep == ev0) ? ev1 : ev0;

        vkeep = ME_Collapse(edge, vkeep, topoflag, &deleted_ents,
                            &merged_entity_pairs);
      }

      if (!vkeep) {
        Errors::Message mesg("Could not collapse degenerate edge. Expect computational issues with connected elements");
        Exceptions::amanzi_throw(mesg);
      }
      else {
        if (deleted_ents) {
          List_Cat(deleted_ents_all, deleted_ents);
          List_Delete(deleted_ents);
        }
        if (merged_entity_pairs) {
          List_Cat(merged_entity_pairs_all, merged_entity_pairs);
          List_Delete(merged_entity_pairs);
        }
      }
    }
  }
  int nmerged = List_Num_Entries(merged_entity_pairs_all)/2;
  for (int j = 0; j < nmerged; j++) {
    MEntity_ptr delent = List_Entry(merged_entity_pairs_all, 2*j);
    MEntity_ptr keepent = List_Entry(merged_entity_pairs_all, 2*j+1);
    merged_ents_info.push_back(static_cast<int>(MEnt_Dim(keepent)));
    merged_ents_info.push_back(MEnt_GlobalID(delent));
    merged_ents_info.push_back(MEnt_GlobalID(keepent));
  }

  int *nmerged_proc = new int[numprocs];
  int *nmerged_proc_x3 = new int[numprocs];
  MPI_Allgather(&nmerged, 1, MPI_INT, nmerged_proc, 1, MPI_INT, mpicomm_);

  int *offset = new int[numprocs];
  int nmerged_global = 0;
  for (int p = 0; p < numprocs; p++) {
    offset[p] = 3*nmerged_global;
    nmerged_global += nmerged_proc[p];
    nmerged_proc_x3[p] = 3*nmerged_proc[p];
  }

  // We probably can make this more efficient by using point-to-point
  // communication

  int *merged_ents_info_global = new int[3*nmerged_global];
  MPI_Allgatherv(&(merged_ents_info[0]), 3*nmerged, MPI_INT,
                 merged_ents_info_global, nmerged_proc_x3, offset,
                 MPI_INT, mpicomm_);

  idx = 0;
  while ((vertex = MESH_Next_Vertex(mesh_, &idx))) {
    if (MV_PType(vertex) == PGHOST) {
      int vgid = MV_GlobalID(vertex);
      for (int i = 0; i < nmerged_global; i++) {
        if (merged_ents_info_global[3*i] == MVERTEX &&
            merged_ents_info_global[3*i+1] == vgid) {
          // Found vertex that got deleted and replaced by another vtx
          // on a different proc
          MV_Set_GlobalID(vertex, merged_ents_info_global[3*i+2]);
          break;
        }
      }
    }
  }

  idx = 0;
  while ((edge = MESH_Next_Edge(mesh_, &idx))) {
    if (ME_PType(edge) == PGHOST) {
      int egid = ME_GlobalID(edge);
      for (int i = 0; i < nmerged_global; i++) {
        if (merged_ents_info_global[3*i] == MEDGE &&
            merged_ents_info_global[3*i+1] == egid) {
          // Found edge that got deleted and replaced by another edge
          // on a different proc
          ME_Set_GlobalID(edge, merged_ents_info_global[3*i+2]);
          break;
        }
      }
    }
  }

  idx = 0;
  while ((face = MESH_Next_Face(mesh_, &idx))) {
    if (MF_PType(face) == PGHOST) {
      int fgid = MF_GlobalID(face);
      for (int i = 0; i < nmerged_global; i++) {
        if (merged_ents_info_global[3*i] == MFACE &&
            merged_ents_info_global[3*i+1] == fgid) {
          // Found face that got deleted and replaced by another face
          // on a different proc
          MF_Set_GlobalID(face, merged_ents_info_global[3*i+2]);
          break;
        }
      }
    }
  }

  delete [] nmerged_proc;
  delete [] nmerged_proc_x3;
  delete [] merged_ents_info_global;
  delete [] offset;

  // Go through all mesh sets and replace any merged entities
  MEntity_ptr delent;
  int nsets = MESH_Num_MSets(mesh_);
  idx = 0;
  while ((delent = List_Next_Entry(merged_entity_pairs_all, &idx))) {
    MEntity_ptr keepent = List_Next_Entry(merged_entity_pairs_all, &idx);
    int entdim = MEnt_Dim(keepent);

    for (int j = 0; j < nsets; j++) {
      MSet_ptr mset = MESH_MSet(mesh_, j);
      if (MSet_EntDim(mset) != entdim) continue;

      int iloc = MSet_Locate(mset, delent);
      if (iloc != -1)  // found deleted entity in set; replace it with keepent
        MSet_Replacei(mset, iloc, keepent);
    }
  }

  // Go through all mesh sets and remove entities that were deleted
  // (not merged)
  idx = 0;
  while ((delent = List_Next_Entry(deleted_ents_all, &idx))) {
    int entdim = MEnt_OrigDim(delent);

    for (int j = 0; j < nsets; j++) {
      MSet_ptr mset = MESH_MSet(mesh_, j);
      if (MSet_EntDim(mset) != entdim) continue;

      int iloc = MSet_Locate(mset, delent);
      if (iloc != -1)  // found deleted entity in set; replace it with keepent
        MSet_Remi(mset, iloc);
    }
  }

  // ME_Collapse only marked these entities as DELETED but now
  // delete them for good
  idx = 0;
  while ((delent = List_Next_Entry(deleted_ents_all, &idx)))
    MEnt_Delete(delent, 0);

  List_Delete(deleted_ents_all);
  List_Delete(merged_entity_pairs_all);

  // Now renumber global IDs to make them contiguous
  //  if (entities_deleted) {
  //   std::cerr << "Entities deleted in collapse_degen_edges ..." << "\n";
  //   MESH_Renumber_GlobalIDs(mesh, MALLTYPE, 0, NULL, mpicomm_);
  //  }

#ifdef DEBUG
  if (MESH_Num_Regions(mesh_) > 0) {  // 3D mesh
    idx = 0;
    while ((face = MESH_Next_Face(mesh_, &idx))) {
      List_ptr fregs = MF_Regions(face);
      if (fregs)
        List_Delete(fregs);
      else {
        std::cerr << "Dangling mesh face with no connected cells AFTER COLLAPSE\n on P" << myprocid << "\n";
        MF_Print(face,3);
      }
    }
  }
#endif
}


//
// Construction function
//
void Mesh_MSTK::init_mesh_from_file_(const std::string& filename)
{
  int ok = 0;
  mesh_ = MESH_New(F1);

  if (filename.find(".exo") != std::string::npos) {  // Exodus file
    // Read the mesh on processor 0
    ok = MESH_ImportFromExodusII(mesh_, filename.c_str(), NULL, mpicomm_);

    // Collapse any degenerate edges in the mesh
    collapse_degen_edges_();  // Assumes its operating on member var 'mesh'

    // Renumber local IDs to be contiguous
    MESH_Renumber(mesh_, 0, MALLTYPE);

    if (numprocs > 1) {
      // Distribute the mesh to all the processors
      int topo_dim = MESH_Num_Regions(mesh_) ? 3 : 2;
      int num_ghost_layers = 1;
      int with_attr = 1;  // Redistribute any attributes and sets
      int method = static_cast<int>(partitioner_);
      int del_inmesh = 1;  // Delete input mesh (on P0) after distribution

      Mesh_ptr globalmesh = mesh_;
      mesh_ = MESH_New(F1);

      ok &= MSTK_Mesh_Distribute(globalmesh, &mesh_, &topo_dim,
                                 num_ghost_layers, with_attr, method,
                                 del_inmesh, mpicomm_);
      if (contiguous_gids_) {
        ok &= MESH_Renumber_GlobalIDs(mesh_, MALLTYPE, 0, NULL, mpicomm_);
      }
    }

  } else if (filename.find(".par") != std::string::npos) {  // Nemesis file
    // Read the individual partitions on each processor
    ok = MESH_ImportFromNemesisI(mesh_, filename.c_str(), NULL, mpicomm_);

    // Collapse any degenerate edges in the mesh
    collapse_degen_edges_();

    // Renumber local IDs to be contiguous
    MESH_Renumber(mesh_, 0, MALLTYPE);

    // Weave the meshes together to form interprocessor connections
    int num_ghost_layers = 1;
    int input_type = 1;  // We are given partitioned meshes with a
                         // unique global ID on each mesh vertex
    int topo_dim = MESH_Num_Regions(mesh_) ? 3 : 2;
    ok &= MSTK_Weave_DistributedMeshes(mesh_, topo_dim, num_ghost_layers,
                                       input_type, mpicomm_);

    if (contiguous_gids_) {
      ok &= MESH_Renumber_GlobalIDs(mesh_, MALLTYPE, 0, NULL, mpicomm_);
    }

  } else {
    Errors::Message msg;
    msg << "Cannot identify file type from extension of input file " <<
      filename << " on processor " << myprocid;
    Exceptions::amanzi_throw(msg);
  }

  if (!ok) {
    Errors::Message msg;
    msg << "Failed to load " << filename << " on processor " <<
      myprocid;
    Exceptions::amanzi_throw(msg);
  }
}


int Mesh_MSTK::generate_regular_mesh_(Mesh_ptr mesh, double x0, double y0,
                                     double z0, double x1, double y1,
                                     double z1, int nx, int ny, int nz)
{
/*

  Index directions for classification templates

  k   j
  |  /
  | /
  |/__ i


  Model vertex, edge and face enumeration for classification templates


         MODEL                   MODEL                  MODEL
         VERTICES                EDGES                  FACES

     7 ______ 8          ___7___           ______
      /|          /|          /|          /|         /|      2   /|
     / |         / |       12/ |8      11/ |             / |  4      / |
   5/______/6 |        /___3___/  |6           /______/  |
    |  |        |  |        |  |        |  |            |  |        | 5|
    |  |____|_|        |  |___5_|_|            |6 |_1___|_|
    |  /3       |  /4      4|  /        |  /            |  /        |  /
    | /         | /         | /9       2| /10           | /      3  | /
    |/_____|/          |/_____|/              |/_____|/
   1             2                1

                                                    Front  - Face1
                                                    Back   - Face2
                                                    Bottom - Face3
                                                    Top    - Face4
                                                    Left   - Face6
                                                    Right  - Face5

  Classification of mesh regions onto multiple material regions is not done
  here since the "geometric model" could have overlapping regions. Instead
  mesh sets are created as necessary based on point location in regions.

*/

  int i, j, k, ii, jj, kk, gid, gdim;
  double xyz[3], dx, dy, dz;
  MVertex_ptr ***verts, mv, rverts[8], fverts[4], everts[2];
  MEdge_ptr me;
  MFace_ptr mf;
  MRegion_ptr mr;
  int vgid_tmpl[3][3][3] = {{{1,4,5},{9,6,12},{3,8,7}},{{1,1,3},{3,1,4},{5,2,7}},{{2,2,6},{10,5,11},{4,6,8}}};
  int vgdim_tmpl[3][3][3]= {{{0,1,0},{1,2,1}, {0,1,0}},{{1,2,1},{2,3,2},{1,2,1}},{{0,1,0},{1,2,1},{0,1,0}}};
  int egdim_tmpl[3][3] = {{1,2,1},{2,3,2},{1,2,1}};
  int egid_tmpl2[3][3] = {{4,6,8},{1,1,2},{2,5,6}};  /* Y direction edges (iterating over i,k) */
  int egid_tmpl1[3][3] = {{9,6,12},{3,1,4},{10,5,11}}; /* Z direction edges (iterating over i,j)*/
  int egid_tmpl0[3][3] = {{1,1,3},{3,1,4},{5,2,7}}; /* X direction edges (iterating over j,k) */
  int fgdim_tmpl[3] = {2,3,2};
  int fgid_tmpl0[3] = {6,1,5};
  int fgid_tmpl1[3] = {1,1,2};
  int fgid_tmpl2[3] = {3,1,4};

  dx = (x1-x0)/nx;
  dy = (y1-y0)/ny;
  dz = (z1-z0)/nz;

  verts = (MVertex_ptr ***) malloc((nx+1)*sizeof(MVertex_ptr **));
  for (j = 0; j < nx+1; ++j) {
    verts[j] = (MVertex_ptr **) malloc((ny+1)*sizeof(MVertex_ptr *));
    for (k = 0; k < ny+1; ++k)
      verts[j][k] = (MVertex_ptr *) malloc((nz+1)*sizeof(MVertex_ptr));
  }

  for (k = 0; k < nz+1; ++k) {
    xyz[2] = (k == nz) ? z1 : z0 + k*dz;
    kk =  (k%nz) ? 1 : (k ? 2 : 0);

    for (j = 0; j < ny+1; ++j) {
      xyz[1] = (j == ny) ? y1 : y0 + j*dy;
      jj = (j%ny) ? 1 : (j ? 2 : 0);

      for (i = 0; i < nx+1; ++i) {
        xyz[0] = (i == nx) ? x1 : x0 + i*dx;
        ii = (i%nx) ? 1 : (i ? 2 : 0);

        mv = MV_New(mesh);
        MV_Set_Coords(mv,xyz);
        verts[i][j][k] = mv;

        gdim  = vgdim_tmpl[ii][jj][kk];
        MV_Set_GEntDim(mv,gdim);

        gid = vgid_tmpl[ii][jj][kk];
        MV_Set_GEntID(mv,gid);
      }
    }
  }


  /* Create the edges explicitly to get the classification right */
  for (i = 0; i < nx+1; ++i) {
    for (j = 0; j < ny+1; ++j) {
      for (k = 0; k < nz; ++k) {
        me = ME_New(mesh);

        everts[0] = verts[i][j][k];
        everts[1] = verts[i][j][k+1];
        ME_Set_Vertex(me,0,everts[0]);
        ME_Set_Vertex(me,1,everts[1]);

        ii = (i%nx) ? 1 : (i ? 2 : 0);
        jj = (j%ny) ? 1 : (j ? 2 : 0);
        gdim = egdim_tmpl[ii][jj];
        gid = egid_tmpl2[ii][jj];

        ME_Set_GEntDim(me,gdim);
        ME_Set_GEntID(me,gid);
      }
    }
  }

  for (i = 0; i < nx+1; ++i) {
    for (k = 0; k < nz+1; ++k) {
      for (j = 0; j < ny; ++j) {
        me = ME_New(mesh);

        everts[0] = verts[i][j][k];
        everts[1] = verts[i][j+1][k];
        ME_Set_Vertex(me,0,everts[0]);
        ME_Set_Vertex(me,1,everts[1]);

        ii = (i%nx) ? 1 : (i ? 2 : 0);
        kk = (k%nz) ? 1 : (k ? 2 : 0);
        gdim = egdim_tmpl[ii][kk];
        gid = egid_tmpl1[ii][kk];

        ME_Set_GEntDim(me,gdim);
        ME_Set_GEntID(me,gid);
      }
    }
  }

  for (j = 0; j < ny+1; ++j) {
    for (k = 0; k < nz+1; ++k) {
      for (i = 0; i < nx; ++i) {
        me = ME_New(mesh);

        everts[0] = verts[i][j][k];
        everts[1] = verts[i+1][j][k];
        ME_Set_Vertex(me,0,everts[0]);
        ME_Set_Vertex(me,1,everts[1]);

        jj = (j%ny) ? 1 : (j ? 2 : 0);
        kk = (k%nz) ? 1 : (k ? 2 : 0);
        gdim = egdim_tmpl[jj][kk];
        gid = egid_tmpl0[jj][kk];

        ME_Set_GEntDim(me,gdim);
        ME_Set_GEntID(me,gid);
      }
    }
  }


  /* Create the faces explicitly to get the classification right */
  for (i = 0; i < nx+1; ++i) {
    for (j = 0; j < ny; ++j) {
      for (k = 0; k < nz; ++k) {
        mf = MF_New(mesh);

        fverts[0] = verts[i][j][k];
        fverts[1] = verts[i][j+1][k];
        fverts[2] = verts[i][j+1][k+1];
        fverts[3] = verts[i][j][k+1];
        MF_Set_Vertices(mf,4,fverts);

        ii = (i%nx) ? 1 : (i ? 2 : 0);
        gdim = fgdim_tmpl[ii];
        gid = fgid_tmpl0[ii];

        MF_Set_GEntDim(mf,gdim);
        MF_Set_GEntID(mf,gid);
      }
    }
  }

  for (j = 0; j < ny+1; ++j) {
    for (i = 0; i < nx; ++i) {
      for (k = 0; k < nz; ++k) {
        mf = MF_New(mesh);

        fverts[0] = verts[i][j][k];
        fverts[1] = verts[i+1][j][k];
        fverts[2] = verts[i+1][j][k+1];
        fverts[3] = verts[i][j][k+1];
        MF_Set_Vertices(mf,4,fverts);

        jj = (j%ny) ? 1 : (j ? 2 : 0);
        gdim = fgdim_tmpl[jj];
        gid = fgid_tmpl1[jj];

        MF_Set_GEntDim(mf,gdim);
        MF_Set_GEntID(mf,gid);
      }
    }
  }

  for (k = 0; k < nz+1; ++k) {
    for (i = 0; i < nx; ++i) {
      for (j = 0; j < ny; ++j) {
        mf = MF_New(mesh);

        fverts[0] = verts[i][j][k];
        fverts[1] = verts[i+1][j][k];
        fverts[2] = verts[i+1][j+1][k];
        fverts[3] = verts[i][j+1][k];
        MF_Set_Vertices(mf,4,fverts);

        kk = (k%nz) ? 1 : (k ? 2 : 0);
        gdim = fgdim_tmpl[kk];
        gid = fgid_tmpl2[kk];

        MF_Set_GEntDim(mf,gdim);
        MF_Set_GEntID(mf,gid);
      }
    }
  }


  /* Not the most efficient way but the easiest to code */

  for (i = 0; i < nx; ++i) {
    for (j = 0; j < ny; ++j) {
      for (k = 0; k < nz; ++k) {
        mr = MR_New(mesh);
        MR_Set_GEntID(mr,1);

        rverts[0] = verts[i][j][k];       rverts[1] = verts[i+1][j][k];
        rverts[2] = verts[i+1][j+1][k];   rverts[3] = verts[i][j+1][k];
        rverts[4] = verts[i][j][k+1];     rverts[5] = verts[i+1][j][k+1];
        rverts[6] = verts[i+1][j+1][k+1]; rverts[7] = verts[i][j+1][k+1];

        MR_Set_Vertices(mr, 8, rverts, 6, NULL);
      }
    }
  }

  for (i = 0; i < nx+1; ++i) {
    for (j = 0; j < ny+1; ++j)
      free(verts[i][j]);
    free(verts[i]);
  }
  free(verts);

  return 1;
}


int Mesh_MSTK::generate_regular_mesh_(Mesh_ptr mesh, double x0, double y0,
                                     double x1, double y1, int nx, int ny)
{
  int i, j, dir[4];
  double xyz[3], dx, dy;
  MVertex_ptr **verts, v0, v1, mv;
  MEdge_ptr fedges[4], me;
  MFace_ptr mf;

  dx = (x1-x0)/nx;
  dy = (y1-y0)/ny;

  verts = (MVertex_ptr **) malloc((nx+1)*sizeof(MVertex_ptr *));
  for (i = 0; i < nx+1; ++i)
    verts[i] = (MVertex_ptr *) malloc((ny+1)*sizeof(MVertex_ptr));

  xyz[2] = 0.0;
  for (j = 0; j < ny+1; ++j) {
    xyz[1] = (j == ny) ? y1 : y0 + j*dy;

    for (i = 0; i < nx+1; ++i) {
      xyz[0] = (i == nx) ? x1 : x0 + i*dx;

      mv = MV_New(mesh);
      MV_Set_Coords(mv,xyz);

      if (i == 0) {
        if (j == 0) {
          MV_Set_GEntDim(mv,0);
          MV_Set_GEntID(mv,1);
        }
        else if (j == ny) {
          MV_Set_GEntDim(mv,0);
          MV_Set_GEntID(mv,4);
        }
        else {
          MV_Set_GEntDim(mv,1);
          MV_Set_GEntID(mv,4);
        }
      }
      else if (i == nx) {
        if (j == 0) {
          MV_Set_GEntDim(mv,0);
          MV_Set_GEntID(mv,2);
        }
        else if (j == ny) {
          MV_Set_GEntDim(mv,0);
          MV_Set_GEntID(mv,3);
        }
        else {
          MV_Set_GEntDim(mv,1);
          MV_Set_GEntID(mv,2);
        }
      }
      else {
        if (j == 0) {
          MV_Set_GEntDim(mv,1);
          MV_Set_GEntID(mv,1);
        }
        else if (j == ny) {
          MV_Set_GEntDim(mv,1);
          MV_Set_GEntID(mv,3);
        }
        else {
          MV_Set_GEntDim(mv,2);
          MV_Set_GEntID(mv,1);
        }
      }

      verts[i][j] = mv;
    }
  }


  for (i = 0; i < nx; ++i) {
    for (j = 0; j < ny; ++j) {
      mf = MF_New(mesh);

      /* edge 0 */
      v0 = verts[i][j];
      v1 = verts[i+1][j];
      fedges[0] = MVs_CommonEdge(v0,v1);
      if (fedges[0])
        dir[0] = (ME_Vertex(fedges[0],0) == v0) ? 1 : 0;
      else {
        me = ME_New(mesh);

        ME_Set_Vertex(me,0,v0);
        ME_Set_Vertex(me,1,v1);

        if (j == 0) {
          ME_Set_GEntDim(me,1);
          ME_Set_GEntID(me,1);
        }
        else {
          ME_Set_GEntDim(me,2);
          ME_Set_GEntID(me,1);
        }

        fedges[0] = me;
        dir[0] = 1;
      }


      /* edge 1 */
      v0 = verts[i+1][j];
      v1 = verts[i+1][j+1];
      fedges[1] = MVs_CommonEdge(v0,v1);
      if (fedges[1])
        dir[1] = (ME_Vertex(fedges[1],0) == v0) ? 1 : 0;
      else {
        me = ME_New(mesh);

        ME_Set_Vertex(me,0,v0);
        ME_Set_Vertex(me,1,v1);

        if (i+1 == nx) {
          ME_Set_GEntDim(me,1);
          ME_Set_GEntID(me,2);
        }
        else {
          ME_Set_GEntDim(me,2);
          ME_Set_GEntID(me,1);
        }

        fedges[1] = me;
        dir[1] = 1;
      }


      /* edge 2 */
      v0 = verts[i+1][j+1];
      v1 = verts[i][j+1];
      fedges[2] = MVs_CommonEdge(v0,v1);
      if (fedges[2])
        dir[2] = (ME_Vertex(fedges[2],0) == v0) ? 1 : 0;
      else {
        me = ME_New(mesh);

        ME_Set_Vertex(me,0,v0);
        ME_Set_Vertex(me,1,v1);

        if (j+1 == nx) {
          ME_Set_GEntDim(me,1);
          ME_Set_GEntID(me,3);
        }
        else {
          ME_Set_GEntDim(me,2);
          ME_Set_GEntID(me,1);
        }

        fedges[2] = me;
        dir[2] = 1;
      }


      /* edge 3 */
      v0 = verts[i][j+1];
      v1 = verts[i][j];
      fedges[3] = MVs_CommonEdge(v0,v1);
      if (fedges[3])
        dir[3] = (ME_Vertex(fedges[3],0) == v0) ? 1 : 0;
      else {
        me = ME_New(mesh);

        ME_Set_Vertex(me,0,v0);
        ME_Set_Vertex(me,1,v1);

        if (i == 0) {
          ME_Set_GEntDim(me,1);
          ME_Set_GEntID(me,4);
        }
        else {
          ME_Set_GEntDim(me,2);
          ME_Set_GEntID(me,1);
        }

        fedges[3] = me;
        dir[3] = 1;
      }


      MF_Set_Edges(mf,4,fedges,dir);

      MF_Set_GEntDim(mf,2);
      MF_Set_GEntID(mf,1);
    }
  }

  for (i = 0; i < nx+1; ++i)
    free(verts[i]);
  free(verts);

  return 1;
}


void Mesh_MSTK::pre_create_steps_(const int space_dimension)
{
  clear_internals_();
  MSTK_Init();
  set_space_dimension(space_dimension);

  auto mpicomm = Teuchos::rcp_dynamic_cast<const MpiComm_type>(get_comm());
  if (!mpicomm.get()) {
    mpicomm_ = MPI_COMM_SELF; // this should never be!
    serial_run = false;
    myprocid = 0;
    numprocs = 1;
  } else {
    mpicomm_ = mpicomm->GetMpiComm();
    myprocid = comm_->MyPID();
    numprocs = comm_->NumProc();
    serial_run = (numprocs == 1);
  }
}


void Mesh_MSTK::inherit_labeled_sets_(MAttrib_ptr copyatt,
                                     List_ptr src_entities)
{
  int idx, idx2, diffdim;
  MSet_ptr mset;

  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm = get_geometric_model();

  if (gm == Teuchos::null) {
    std::cerr << "Need region definitions to initialize sets" << std::endl;
    return;
  }

  Mesh_ptr parent_mstk_mesh = parent_mesh_->mesh_;

  // Difference in cell dimension of this mesh and its parent
  // Labeled set entity dimensions will be similarly dialed down

  diffdim = parent_mesh_->get_manifold_dimension() - get_manifold_dimension();
  if (diffdim > 1) {
    Errors::Message mesg("Dimension of mesh and its parent differ by more than 1");
    Exceptions::amanzi_throw(mesg);
  }

  unsigned int ngr = gm->size();

  for (int i = 0; i < ngr; ++i) {
    Teuchos::RCP<const AmanziGeometry::Region> rgn = gm->FindRegion(i);

    if (rgn->type() == AmanziGeometry::LABELEDSET) {

      // Get the set from the parent mesh

      Teuchos::RCP<const AmanziGeometry::RegionLabeledSet> lsrgn =
          Teuchos::rcp_static_cast<const AmanziGeometry::RegionLabeledSet>(rgn);

      std::string internal_name;
      std::string label = lsrgn->label();

      if (lsrgn->entity_str() == "CELL")
        internal_name = internal_name_of_set_(*lsrgn,Entity_kind::CELL);
      else if (lsrgn->entity_str() == "FACE")
        internal_name = internal_name_of_set_(*lsrgn,Entity_kind::FACE);
      else if (lsrgn->entity_str() == "NODE")
        internal_name = internal_name_of_set_(*lsrgn,Entity_kind::NODE);


      MSet_ptr mset_parent = MESH_MSetByName(parent_mstk_mesh,
                                             internal_name.c_str());
      if (!mset_parent)
        continue;

      // Also, if this is a lower dimensional mesh (like a surface
      // mesh created from a solid mesh) and the set contains entities
      // from which it was created (like a face set) then don't
      // inherit this set - otherwise we will get odd things like
      // internal edges in the surface mesh being labeled as "face
      // sets"

      if (diffdim > 0) {
        int found = 0;
        idx = 0;
        MEntity_ptr ent;
        while ((ent = List_Next_Entry(src_entities, &idx))) {
          if (MSet_Contains(mset_parent, ent)) {
            found = 1;
            break;
          }
        }
        if (found) continue;
      }

      // Create the set in this mesh

      MType subentdim;
      MType entdim = MSet_EntDim(mset_parent);
      if (entdim == MVERTEX)
        subentdim = MVERTEX;
      else
        subentdim = (MType) (entdim-diffdim);

      mset = MSet_New(mesh_,internal_name.c_str(),subentdim);


      // Populate the set

      int mkid = MSTK_GetMarker();

      MEntity_ptr ent;
      idx = 0;
      while ((ent = MSet_Next_Entry(mset_parent,&idx))) {
        if (MEnt_Dim(ent) == MDELETED)
          continue;
        MEntity_ptr copyent;
        int ival;
        double rval;

        if (subentdim == entdim) {
          MEnt_Get_AttVal(ent,copyatt,&ival,&rval,&copyent);
          if (!copyent) continue;

          MSet_Add(mset,copyent);
        }
        else {
          if (entdim == MREGION) {
            MFace_ptr rf;
            List_ptr rfaces = MR_Faces((MRegion_ptr)ent);
            idx2 = 0;
            while ((rf = List_Next_Entry(rfaces,&idx2))) {
              MEnt_Get_AttVal(rf,copyatt,&ival,&rval,&copyent);
              if (!copyent) continue;

              if (!MEnt_IsMarked(copyent,mkid)) {
                MSet_Add(mset,copyent);
                MEnt_Mark(copyent,mkid);
              }
            }
            List_Delete(rfaces);
          }
          else if (entdim == MFACE) {
            MEdge_ptr fe;
            List_ptr fedges = MF_Edges((MFace_ptr)ent,1,0);
            idx2 = 0;
            while ((fe = List_Next_Entry(fedges,&idx2))) {
              MEnt_Get_AttVal(fe,copyatt,&ival,&rval,&copyent);
              if (!copyent) continue;

              if (!MEnt_IsMarked(copyent,mkid)) {
                MSet_Add(mset,copyent);
                MEnt_Mark(copyent,mkid);
              }
            }
            List_Delete(fedges);
          }
        }

      }

      MSet_Unmark(mset,mkid);
      MSTK_FreeMarker(mkid);

    }
  }
}


//---------------------------------------------------------
// Extract a list of MSTK entities and make a new MSTK mesh
// For private use of Mesh_MSTK class only
//---------------------------------------------------------
void Mesh_MSTK::extract_mstk_mesh_(List_ptr src_entities,
                                  const MType entity_dim,
                                  const bool flatten,
                                  const bool request_faces,
                                  const bool request_edges)
{
  int ival = 0, idx;
  double rval = 0., xyz[3];
  void *pval;

  AMANZI_ASSERT(parent_mesh_.get());
  Mesh_ptr parent_mesh_mstk = parent_mesh_->mesh_;

  // Make sure Global ID searches are enabled
  MESH_Enable_GlobalIDSearch(parent_mesh_mstk);

  if (flatten) {
    if (entity_dim == MREGION || entity_dim == MVERTEX) {
      Errors::Message mesg("Flattening or extruding allowed only for sets of FACEs in volume mesh or CELLs in surface meshes");
      Exceptions::amanzi_throw(mesg);
    }
  }

  if (entity_dim == MEDGE) {
    Errors::Message mesg("Requested mesh constructor produces 1D mesh which is not supported by Amanzi");
    Exceptions::amanzi_throw(mesg);
  }

  // Pre-processing (init, MPI queries etc)
  if (flatten)
    pre_create_steps_(get_parent()->get_space_dimension()-1);
  else
    pre_create_steps_(get_parent()->get_space_dimension());

  // What is the cell dimension of new mesh
  switch (entity_dim) {
  case MREGION:
    set_manifold_dimension(3); // extract regions/cells from mesh
    break;
  case MFACE:
    set_manifold_dimension(2); // extract faces from mesh
    break;
  case MEDGE: {
    Errors::Message mesg("Edge list passed into extract mesh. Cannot extract a wire or point mesh");
    Exceptions::amanzi_throw(mesg);
    break;
  }
  case MVERTEX: {
    Errors::Message mesg("Vertex list passed into extract mesh. Cannot extract a point mesh");
    Exceptions::amanzi_throw(mesg);
    break;
  }
  default: {
    Errors::Message mesg1("Unrecognized Entity_kind");
    Exceptions::amanzi_throw(mesg1);
  }
  }

  // Create new mesh in MSTK
  mesh_ = MESH_New(MESH_RepType(parent_mesh_mstk));

  // Have to do some additional work for extruding an extracted mesh
  // Extrusion applicable only in the case of entdim = MFACE/MEDGE
  MAttrib_ptr copyatt = MAttrib_New(parent_mesh_mstk,"copyatt",POINTER,MALLTYPE);
  vparentatt_ = MAttrib_New(mesh_,"vparentatt_",POINTER,MVERTEX);
  eparentatt_ = MAttrib_New(mesh_,"eparentatt_",POINTER,MEDGE);
  fparentatt_ = MAttrib_New(mesh_,"fparentatt_",POINTER,MFACE);
  rparentatt_ = MAttrib_New(mesh_,"rparentatt_",POINTER,MREGION);

  switch (entity_dim) {
    case MREGION: {  // Extracting a subset of a solid mesh
      idx = 0;
      MRegion_ptr mr;
      while ((mr = (MRegion_ptr) List_Next_Entry(src_entities,&idx))) {

        List_ptr rfaces = MR_Faces(mr);
        int nrf = List_Num_Entries(rfaces);
        MFace_ptr rfaces_new[MAXPF3];
        int rfdirs_new[MAXPF3];
        for (int i = 0; i < nrf; ++i) {
          MFace_ptr mf = List_Entry(rfaces,i);

          MEnt_Get_AttVal(mf,copyatt,&ival,&rval,&pval);
          if (pval) {
            rfaces_new[i] = pval;
            rfdirs_new[i] = MR_FaceDir_i(mr,i);
          }
          else {

            List_ptr fverts = MF_Vertices(mf,1,0);
            int nfv = List_Num_Entries(fverts);
            MVertex_ptr fverts_new[MAXPV2];
            for (int j = 0; j < nfv; ++j) {
              MVertex_ptr mv = List_Entry(fverts,j);
              MEnt_Get_AttVal(mv,copyatt,&ival,&rval,&pval);
              if (pval)
                fverts_new[j] = pval;
              else {
                fverts_new[j] = MV_New(mesh_);
                MV_Coords(mv,xyz);
                MV_Set_Coords(fverts_new[j],xyz);
                MV_Set_GEntDim(fverts_new[j],MV_GEntDim(mv));
                MV_Set_GEntID(fverts_new[j],MV_GEntID(mv));
                MEnt_Set_AttVal(mv,copyatt,ival,rval,fverts_new[j]);
                MEnt_Set_AttVal(fverts_new[j],vparentatt_,0,0.0,mv);
              }
            }
            List_Delete(fverts);

            rfaces_new[i] = MF_New(mesh_);
            MF_Set_Vertices(rfaces_new[i],nfv,fverts_new);
            MF_Set_GEntDim(rfaces_new[i],MF_GEntDim(mf));
            MF_Set_GEntID(rfaces_new[i],MF_GEntID(mf));
            rfdirs_new[i] = MR_FaceDir_i(mr,i);

            MEnt_Set_AttVal(mf,copyatt,ival,rval,rfaces_new[i]);
            MEnt_Set_AttVal(rfaces_new[i],fparentatt_,0,0.0,mf);
          }
        }
        List_Delete(rfaces);

        MRegion_ptr mr_new = MR_New(mesh_);
        MR_Set_Faces(mr_new,nrf,rfaces_new,rfdirs_new);
        MR_Set_GEntID(mr_new,MR_GEntID(mr));

        MEnt_Set_AttVal(mr,copyatt,ival,rval,mr_new);
        MEnt_Set_AttVal(mr_new,rparentatt_,0,0.0,mr);
      }

      break;
    }
    case MFACE: {  // Extracting a surface from a solid mesh or subset of
      //           // a surface mesh
      idx = 0;
      MFace_ptr mf = nullptr;
      while ((mf = (MFace_ptr) List_Next_Entry(src_entities,&idx))) {
        List_ptr fedges = MF_Edges(mf,1,0);
        int nfe = List_Num_Entries(fedges);
        int fedirs[MAXPV2];
        MEdge_ptr fedges_new[MAXPV2];
        for (int j = 0; j < nfe; ++j) {
          MEdge_ptr me = List_Entry(fedges,j);
          MEnt_Get_AttVal(me,copyatt,&ival,&rval,&pval);
          if (pval)
            fedges_new[j] = pval;
          else {
            fedges_new[j] = ME_New(mesh_);

            for (int k = 0; k < 2; ++k) {
              MVertex_ptr mv = ME_Vertex(me,k);
              MVertex_ptr mv_new = nullptr;
              MEnt_Get_AttVal(mv,copyatt,&ival,&rval,&pval);
              if (pval)
                mv_new = pval;
              else {
                MV_Coords(mv,xyz);
                if (flatten) xyz[2] = 0.0;
                mv_new = MV_New(mesh_);
                MV_Set_Coords(mv_new,xyz);
                MV_Set_GEntDim(mv_new,MV_GEntDim(mv));
                MV_Set_GEntID(mv_new,MV_GEntID(mv));
                MEnt_Set_AttVal(mv,copyatt,ival,rval,mv_new);
                MEnt_Set_AttVal(mv_new,vparentatt_,0,0.0,mv);
              }

              ME_Set_Vertex(fedges_new[j],k,mv_new);
              ME_Set_GEntDim(fedges_new[j],ME_GEntDim(me));
              ME_Set_GEntID(fedges_new[j],ME_GEntID(me));
              MEnt_Set_AttVal(me,copyatt,ival,rval,fedges_new[j]);
              MEnt_Set_AttVal(fedges_new[j],eparentatt_,0,0.0,me);
            }
          }
          fedirs[j] = MF_EdgeDir_i(mf,j);
        }
        List_Delete(fedges);

        MFace_ptr mf_new = MF_New(mesh_);
        MF_Set_Edges(mf_new,nfe,fedges_new,fedirs);
        MF_Set_GEntDim(mf_new,2);  // This has to be surface mesh
        if (MF_GEntDim(mf) == 2)
          MF_Set_GEntID(mf_new,MF_GEntID(mf));

        MEnt_Set_AttVal(mf,copyatt,ival,rval,mf_new);
        MEnt_Set_AttVal(mf_new,fparentatt_,0,0.0,mf);
      }

      break;
    }
    case MEDGE: {  // Extracting a wire mesh from a solid or surface mesh
      idx = 0;
      MEdge_ptr me = nullptr;
      while ((me = (MEdge_ptr) List_Next_Entry(src_entities,&idx))) {

        MEdge_ptr me_new = ME_New(mesh_);

        for (int j = 0; j < 2; ++j)  {
          MVertex_ptr mv = ME_Vertex(me,j);

          MVertex_ptr mv_new = nullptr;

          MEnt_Get_AttVal(mv,copyatt,&ival,&rval,&pval);
          if (pval)
            mv_new = pval;
          else {
            MV_Coords(mv,xyz);
            if (flatten) {
              xyz[1] = 0.0;
              xyz[2] = 0.0;
            }
            mv_new = MV_New(mesh_);
            MV_Set_Coords(mv_new,xyz);
            MV_Set_GEntDim(mv_new,MV_GEntDim(mv));
            MV_Set_GEntID(mv_new,MV_GEntID(mv));

            MEnt_Set_AttVal(mv,copyatt,ival,rval,mv_new);
            MEnt_Set_AttVal(mv_new,vparentatt_,0,0.0,mv);
          }

          ME_Set_Vertex(me_new,j,mv_new);
        }

        if (ME_GEntDim(me) == 1)
          ME_Set_GEntDim(me_new, 1);
        MEnt_Set_AttVal(me,copyatt,ival,rval,me_new);
        MEnt_Set_AttVal(me_new,eparentatt_,0,0.0,me);
      }
      break;
    }
    case MVERTEX: {

      idx = 0;
      MVertex_ptr mv = nullptr;
      while ((mv = (MVertex_ptr) List_Next_Entry(src_entities,&idx))) {
        MVertex_ptr mv_new = MV_New(mesh_);
        MV_Set_Coords(mv_new,xyz);
        if (flatten) xyz[2] = 0.0;
        MV_Set_GEntDim(mv_new,MV_GEntDim(mv));
        MV_Set_GEntID(mv_new,MV_GEntID(mv));

        MEnt_Set_AttVal(mv,copyatt,ival,rval,mv_new);
        MEnt_Set_AttVal(mv_new,vparentatt_,0,0.0,mv);
      }
      break;
    }
    default: {
      Errors::Message mesg("Unknown entity type");
      Exceptions::amanzi_throw(mesg);
    }
  }

  if (!serial_run) {
    // Have to assign global IDs and build ghost entities
    int num_ghost_layers = 1;
    int input_type = 0; /* No parallel info is given */
    int status = MSTK_Weave_DistributedMeshes(mesh_, get_manifold_dimension(),
                                              num_ghost_layers, input_type, mpicomm_);

    // Now we have to build parent information for global entities
    MAttrib_ptr vparentgid_att = MAttrib_New(mesh_,"vparent_gid",INT,MVERTEX);
    MAttrib_ptr eparentgid_att = MAttrib_New(mesh_,"eparent_gid",INT,MEDGE);
    MAttrib_ptr fparentgid_att = MAttrib_New(mesh_,"fparent_gid",INT,MFACE);
    MAttrib_ptr rparentgid_att = MAttrib_New(mesh_,"rparent_gid",INT,MREGION);

    // Attach parent global ID info to entities used by other processors
    idx = 0;
    MVertex_ptr mv = nullptr;
    while ((mv = (MVertex_ptr) MESH_Next_Vertex(mesh_,&idx)))
      if (MV_PType(mv) == POVERLAP) {
        MEnt_Get_AttVal(mv,vparentatt_,&ival,&rval,&pval);
        MEnt_Set_AttVal(mv,vparentgid_att,MV_GlobalID((MVertex_ptr)pval),0.0,
                        NULL);
      }
    idx = 0;
    MEdge_ptr me = nullptr;
    while ((me = (MEdge_ptr) MESH_Next_Edge(mesh_,&idx)))
      if (ME_PType(me) == POVERLAP) {
        MEnt_Get_AttVal(me,eparentatt_,&ival,&rval,&pval);
        MEnt_Set_AttVal(me,eparentgid_att,ME_GlobalID((MEdge_ptr)pval),0.0,
                        NULL);
      }
    idx = 0;
    MFace_ptr mf = nullptr;
    while ((mf = (MFace_ptr) MESH_Next_Face(mesh_,&idx)))
      if (MF_PType(mf) == POVERLAP) {
        MEnt_Get_AttVal(mf,fparentatt_,&ival,&rval,&pval);
        MEnt_Set_AttVal(mf,fparentgid_att,MF_GlobalID((MFace_ptr)pval),0.0,
                        NULL);
      }
    idx = 0;
    MRegion_ptr mr = nullptr;
    while ((mr = (MRegion_ptr) MESH_Next_Region(mesh_,&idx)))
      if (MR_PType(mr) == POVERLAP) {
        MEnt_Get_AttVal(mr,rparentatt_,&ival,&rval,&pval);
        MEnt_Set_AttVal(mr,rparentgid_att,MR_GlobalID((MRegion_ptr)pval),0.0,
                        NULL);
      }

    // Update attributes on ghost entities - this will ensure that
    // ghost entities have their parent global ID information
    status &= MESH_UpdateAttributes(mesh_,mpicomm_);

    // Now reverse engineer the parents of ghost entities from the global IDs
    idx = 0;
    while ((mv = (MVertex_ptr) MESH_Next_GhostVertex(mesh_,&idx))) {
      MEnt_Get_AttVal(mv,vparentgid_att,&ival,&rval,&pval);
      MVertex_ptr mv_parent = MESH_VertexFromGlobalID(parent_mesh_mstk,ival);
      if (!mv_parent) {
        Errors::Message
          mesg("Cannot find ghost vertex with given global ID");
        Exceptions::amanzi_throw(mesg);
      }
      MEnt_Set_AttVal(mv,vparentatt_,0,0.0,mv_parent);
    }
    idx = 0;
    while ((me = (MEdge_ptr) MESH_Next_GhostEdge(mesh_,&idx))) {
      MEnt_Get_AttVal(me,eparentgid_att,&ival,&rval,&pval);
      MEdge_ptr me_parent = MESH_EdgeFromGlobalID(parent_mesh_mstk,ival);
      if (!me_parent) {
        Errors::Message
          mesg("Cannot find ghost edge with given global ID");
        Exceptions::amanzi_throw(mesg);
      }
      MEnt_Set_AttVal(me,eparentatt_,0,0.0,me_parent);
    }
    idx = 0;
    while ((mf = (MFace_ptr) MESH_Next_GhostFace(mesh_,&idx))) {
      MEnt_Get_AttVal(mf,fparentgid_att,&ival,&rval,&pval);
      MFace_ptr mf_parent = MESH_FaceFromGlobalID(parent_mesh_mstk,ival);
      if (!mf_parent) {
        Errors::Message
          mesg("Cannot find ghost face with given global ID");
        Exceptions::amanzi_throw(mesg);
      }
      MEnt_Set_AttVal(mf,fparentatt_,0,0.0,mf_parent);
    }
    idx = 0;
    while ((mr = (MRegion_ptr) MESH_Next_GhostRegion(mesh_,&idx))) {
      MEnt_Get_AttVal(mr,rparentgid_att,&ival,&rval,&pval);
      MRegion_ptr mr_parent = MESH_RegionFromGlobalID(parent_mesh_mstk,ival);
      if (!mr_parent) {
        Errors::Message
          mesg("Cannot find ghost region with given global ID");
        Exceptions::amanzi_throw(mesg);
      }
      MEnt_Set_AttVal(mr,rparentatt_,0,0.0,mr_parent);
    }

    MAttrib_Delete(vparentgid_att);
    MAttrib_Delete(eparentgid_att);
    MAttrib_Delete(fparentgid_att);
    MAttrib_Delete(rparentgid_att);
  }

  // We have to do an extra step to build new labeled sets based on
  // labeled sets of the base mesh
  inherit_labeled_sets_(copyatt, src_entities);

  // Do all the processing required for setting up the mesh for Amanzi
  post_create_steps_();

  // Clean up
  switch (entity_dim) {
  case MREGION: {
    MRegion_ptr mr = nullptr;
    idx = 0;
    while ((mr = (MRegion_ptr) List_Next_Entry(src_entities,&idx))) {

      List_ptr rfaces = MR_Faces(mr);
      int nrf = List_Num_Entries(rfaces);

      for (int i = 0; i < nrf; ++i) {
        MFace_ptr mf = List_Entry(rfaces,i);
        MEnt_Rem_AttVal(mf,copyatt);

        List_ptr fverts = MF_Vertices(mf,1,0);
        int nfv = List_Num_Entries(fverts);

        for (int j = 0; j < nfv; ++j) {
          MVertex_ptr mv = List_Entry(fverts,j);
          MEnt_Rem_AttVal(mv,copyatt);
        }
        List_Delete(fverts);

        MEnt_Rem_AttVal(mf,copyatt);
      }
      List_Delete(rfaces);

      MEnt_Rem_AttVal(mr,copyatt);
    }
    break;
  }
    case MFACE: {
    MFace_ptr mf = nullptr;
    idx = 0;
    while ((mf = (MFace_ptr) List_Next_Entry(src_entities,&idx))) {

      List_ptr fedges = MF_Edges(mf,1,0);
      int nfe = List_Num_Entries(fedges);
      for (int j = 0; j < nfe; ++j) {
        MEdge_ptr me = List_Entry(fedges,j);
        MEnt_Rem_AttVal(me,copyatt);
        MVertex_ptr mv = ME_Vertex(me,MF_EdgeDir_i(mf,j));
        MEnt_Rem_AttVal(mv,copyatt);
      }
      List_Delete(fedges);

      MEnt_Rem_AttVal(mf,copyatt);
    }

    break;
  }
    case MEDGE: {
    MEdge_ptr me;
    idx = 0;
    while ((me = (MEdge_ptr) List_Next_Entry(src_entities,&idx))) {
      for (int j = 0; j < 2; ++j)  {
        MVertex_ptr mv = ME_Vertex(me,j);
        MEnt_Rem_AttVal(mv,copyatt);
      }
      MEnt_Rem_AttVal(me,copyatt);
    }

    break;
  }
  case MVERTEX: {
    MVertex_ptr mv = nullptr;
    idx = 0;
    while ((mv = (MVertex_ptr) List_Next_Entry(src_entities,&idx)))
      MEnt_Rem_AttVal(mv,copyatt);

    break;
  }
  default: {
    Errors::Message mesg("Unknown entity type");
    Exceptions::amanzi_throw(mesg);
  }
  }

  MAttrib_Delete(copyatt);
}


//---------------------------------------------------------
// Write mesh out to exodus file
//---------------------------------------------------------
void
Mesh_MSTK::writeToExodusFile(const std::string& filename) const {
  MESH_ExportToExodusII(mesh_, filename.c_str(), -1, NULL, NULL, mpicomm_);
}


// Run MSTK's internal checks - meant for debugging only
// Returns true if everything is ok, false otherwise
bool
Mesh_MSTK::run_internal_mstk_checks() const {
  return MESH_CheckTopo(mesh_) && MESH_Parallel_Check(mesh_, mpicomm_);
}


}  // namespace AmanziMesh
}  // namespace Amanzi


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
  comm->MaxAll(&cell_dim, &max, 1);

  if (max != cell_dim) {
    Errors::Message mesg("cell dimension on this processor is different from max cell dimension "
                         "across all processors");
    Exceptions::amanzi_throw(mesg);
  }
  setManifoldDimension(cell_dim);

  if (cell_dim == 2 && space_dim == 3) {
    // Check if this is a completely planar mesh
    // in which case one can label the space dimension as 2
    //
    // cannot use getNodeCoordinate() yet!
    MVertex_ptr mv = nullptr, mv0 = MESH_Vertex(mesh_, 0);
    double vxyz[3], z0;
    MV_Coords(mv0, vxyz);
    z0 = vxyz[2];

    bool planar = true;
    int idx = 0;
    while ((mv = MESH_Next_Vertex(mesh_, &idx))) {
      MV_Coords(mv, vxyz);
      if (z0 != vxyz[2]) {
        planar = false;
        break;
      }
    }

    if (planar) space_dim = 2;
    comm->MaxAll(&space_dim, &max, 1);
    space_dim = max;
    setSpaceDimension(space_dim);
  }

  // Verify mesh and geometric model compatibility
  if (gm != Teuchos::null && gm->dimension() != getSpaceDimension()) {
    Errors::Message msg("Geometric model and mesh have different dimensions.");
    Exceptions::amanzi_throw(msg);
  }

  // Do all the processing required for setting up the mesh for Amanzi
  post_create_steps_();
}


//--------------------------------------
// Construct a 3D regular hexahedral mesh internally
//--------------------------------------
Mesh_MSTK::Mesh_MSTK(const double x0,
                     const double y0,
                     const double z0,
                     const double x1,
                     const double y1,
                     const double z1,
                     const unsigned int nx,
                     const unsigned int ny,
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
    ok &= generate_regular_mesh_(mesh_, x0, y0, z0, x1, y1, z1, nx, ny, nz);
    setManifoldDimension(3);
    myprocid = 0;

  } else {
    Mesh_ptr globalmesh = nullptr;
    int topo_dim = 3;   // What is the topological dimension of the mesh
    int ring = 1;       // One layer of ghost cells in parallel meshes
    int with_attr = 1;  // update of attributes in parallel meshes
    int del_inmesh = 1; // delete input mesh as soon as possible
    int method = static_cast<int>(partitioner_);

    if (myprocid == 0) {
      globalmesh = MESH_New(F1);
      ok &= generate_regular_mesh_(globalmesh, x0, y0, z0, x1, y1, z1, nx, ny, nz);
      topo_dim = (MESH_Num_Regions(globalmesh) == 0) ? 2 : 3;

    } else {
      globalmesh = nullptr;
    }

    ok = ok & MSTK_Mesh_Distribute(
                globalmesh, &mesh_, &topo_dim, ring, with_attr, method, del_inmesh, mpicomm_);
    setManifoldDimension(topo_dim);
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
Mesh_MSTK::Mesh_MSTK(const double x0,
                     const double y0,
                     const double x1,
                     const double y1,
                     const int nx,
                     const int ny,
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
  setManifoldDimension(topo_dim);

  if (serial_run) {
    // Load serial mesh
    mesh_ = MESH_New(F1);
    ok &= generate_regular_mesh_(mesh_, x0, y0, x1, y1, nx, ny);
    myprocid = 0;

  } else {
    Mesh_ptr globalmesh = nullptr;
    int ring = 1;       // One layer of ghost cells in parallel meshes
    int with_attr = 1;  // update of attributes in parallel meshes
    int del_inmesh = 1; // delete input mesh at the earliest
    int method = static_cast<int>(partitioner_);

    if (myprocid == 0) {
      globalmesh = MESH_New(F1);
      ok &= generate_regular_mesh_(globalmesh, x0, y0, x1, y1, nx, ny);
      topo_dim = (MESH_Num_Regions(globalmesh) == 0) ? 2 : 3;

    } else {
      globalmesh = nullptr;
    }

    ok &= MSTK_Mesh_Distribute(
      globalmesh, &mesh_, &topo_dim, ring, with_attr, method, del_inmesh, mpicomm_);
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
                     const Entity_ID_View& entity_ids,
                     const Entity_kind entity_kind,
                     const bool flatten,
                     const Comm_ptr_type& comm,
                     const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
                     const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : MeshFramework(comm,
                  gm == Teuchos::null ? parent_mesh->getGeometricModel() : gm,
                  plist == Teuchos::null ?
                    Teuchos::rcp(new Teuchos::ParameterList(*parent_mesh->getParameterList())) :
                    plist),
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
  parent_mesh_ = parent_mesh_as_mstk;
  auto parent_mesh_mstk = parent_mesh_as_mstk->mesh_;

  // store pointers to the MESH_XXXFromID functions so that they can
  // be called without a switch statement
  static MEntity_ptr (*MEntFromID[4])(
    Mesh_ptr, int) = { MESH_VertexFromID, MESH_EdgeFromID, MESH_FaceFromID, MESH_RegionFromID };
  MType entity_dim = parent_mesh_as_mstk->entity_kind_to_mtype_(entity_kind);

  // Also make sure that the mesh object can do fast queries on local IDs
  MESH_Enable_LocalIDSearch(parent_mesh_mstk);

  // collect MSTK entities and extract
  int nent = entity_ids.size();
  List_ptr src_ents = List_New(nent);
  for (int i = 0; i < nent; ++i) {
    MEntity_ptr ent = MEntFromID[entity_dim](parent_mesh_mstk, entity_ids[i] + 1);
    List_Add(src_ents, ent);
  }
  extract_mstk_mesh_(src_ents, entity_dim, flatten, faces_requested_, edges_requested_);
  List_Delete(src_ents);
}


//---------------------------------------------------------
// Destructor with cleanup
//---------------------------------------------------------
Mesh_MSTK::~Mesh_MSTK()
{
  if (faceflip_) delete[] faceflip_;
  if (edgeflip_) delete[] edgeflip_;

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
void
Mesh_MSTK::read_plist_()
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
Mesh_MSTK::internal_name_of_set_(const AmanziGeometry::RegionLabeledSet& rgn,
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
Mesh_MSTK::other_internal_name_of_set_(const AmanziGeometry::RegionLabeledSet& rgn,
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
Mesh_MSTK::getNumEntities(const Entity_kind kind, const Parallel_kind ptype) const
{
  switch (kind) {
  case Entity_kind::NODE:
    switch (ptype) {
    case Parallel_kind::OWNED:
      return MSet_Num_Entries(owned_verts_);
      break;
    case Parallel_kind::GHOST:
      return !serial_run ? MSet_Num_Entries(ghost_verts_) : 0;
      break;
    case Parallel_kind::ALL:
      return MESH_Num_Vertices(mesh_);
      break;
    default:
      return 0;
    }
    break;

  case Entity_kind::EDGE:
    AMANZI_ASSERT(edges_initialized_);
    switch (ptype) {
    case Parallel_kind::OWNED:
      return MSet_Num_Entries(owned_edges_);
      break;
    case Parallel_kind::GHOST:
      return !serial_run ? MSet_Num_Entries(ghost_edges_) : 0;
      break;
    case Parallel_kind::ALL:
      return MESH_Num_Edges(mesh_);
      break;
    default:
      return 0;
    }
    break;

  case Entity_kind::FACE:
    AMANZI_ASSERT(faces_initialized_);
    switch (ptype) {
    case Parallel_kind::OWNED:
      return MSet_Num_Entries(owned_faces_);
      break;
    case Parallel_kind::GHOST:
      return !serial_run ? MSet_Num_Entries(ghost_faces_) : 0;
      break;
    case Parallel_kind::ALL:
      return (getManifoldDimension() == 2 ? MESH_Num_Edges(mesh_) : MESH_Num_Faces(mesh_));
      break;
    default:
      return 0;
    }
    break;

  case Entity_kind::CELL:
    switch (ptype) {
    case Parallel_kind::OWNED:
      return MSet_Num_Entries(owned_cells_);
      break;
    case Parallel_kind::GHOST:
      return !serial_run ? MSet_Num_Entries(ghost_cells_) : 0;
      break;
    case Parallel_kind::ALL:
      return ((getManifoldDimension() == 2) ? MESH_Num_Faces(mesh_) : MESH_Num_Regions(mesh_));
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
Cell_kind
Mesh_MSTK::getCellType(const Entity_ID cellid) const
{
  MEntity_ptr cell = cell_id_to_handle_[cellid];
  int ival;
  MEnt_Get_AttVal(cell, celltype_att_, &ival, NULL, NULL);
  return (Cell_kind)ival;
}


//---------------------------------------------------------
// Get faces of a cell and directions in which the cell uses the face

// The Amanzi coding guidelines regarding function arguments is purposely
// violated here to allow for a default input argument

// On a distributed mesh, this will return all the faces of the
// cell, OWNED or GHOST. If cells_initialized_ = true, the faces will be
// returned in a standard order according to Exodus II convention
// for standard cells; in all other situations (cells_initialized_ = false or
// non-standard cells), the list of faces will be in arbitrary order

// In 3D, direction is 1 if face normal points out of cell
// and -1 if face normal points into cell
// In 2D, direction is 1 if face/edge is defined in the same
// direction as the cell polygon, and -1 otherwise
//---------------------------------------------------------
void
Mesh_MSTK::getCellFacesAndDirs_ordered_(
  const Entity_ID cellid,
  View_type<const Entity_ID, MemSpace_kind::HOST>& faceids,
  View_type<const Direction_type, MemSpace_kind::HOST>* face_dirs) const
{
  if (getManifoldDimension() == 3) {
    Cell_kind celltype = getCellType(cellid);
    if (celltype == Cell_kind::TET || celltype == Cell_kind::PRISM ||
        celltype == Cell_kind::PYRAMID || celltype == Cell_kind::HEX) {
      Entity_ID_View lfaceids;
      Entity_Direction_View lface_dirs;
      int lid, nf;
      MEntity_ptr cell = cell_id_to_handle_[cellid];

      List_ptr rfaces = MR_Faces((MRegion_ptr)cell);
      nf = List_Num_Entries(rfaces);

      Kokkos::resize(lfaceids, nf);
      if (face_dirs) Kokkos::resize(lface_dirs, nf);

      /* base face */
      MFace_ptr face0 = nullptr;
      int fdir0 = 0;

      if (celltype == Cell_kind::TET || celltype == Cell_kind::HEX) {
        face0 = List_Entry(rfaces, 0);
        fdir0 = MR_FaceDir_i((MRegion_ptr)cell, 0);

      } else if (celltype == Cell_kind::PRISM) { /* Find the first triangular face */
        for (int i = 0; i < 5; ++i) {
          MFace_ptr face = List_Entry(rfaces, i);
          if (MF_Num_Edges(face) == 3) {
            face0 = face;
            fdir0 = MR_FaceDir_i((MRegion_ptr)cell, i);
            break;
          }
        }

      } else if (celltype == Cell_kind::PYRAMID) { /* Find the quad face */
        for (int i = 0; i < 5; ++i) {
          MFace_ptr face = List_Entry(rfaces, i);
          if (MF_Num_Edges(face) == 4) {
            face0 = face;
            fdir0 = MR_FaceDir_i((MRegion_ptr)cell, i);
            break;
          }
        }
      }

      /* Markers for faces to avoid searching */
      int mkid = MSTK_GetMarker();
      MEnt_Mark(face0, mkid);

      /* Add all lateral faces first (faces adjacent to the base face) */
      List_ptr fedges0 = MF_Edges(face0, !fdir0, 0);
      int idx = 0;
      MEdge_ptr fe;
      nf = 0;
      while ((fe = List_Next_Entry(fedges0, &idx))) {
        /* Is there an unprocessed face in this region that is
           adjacent to this edge */
        int idx2 = 0;
        MFace_ptr fadj = nullptr;
        int i = 0;
        while ((fadj = List_Next_Entry(rfaces, &idx2))) {
          if (fadj != face0 && !MEnt_IsMarked(fadj, mkid)) {
            if (MF_UsesEntity(fadj, fe, MEDGE)) {
              lid = MEnt_ID(fadj);
              lfaceids[nf] = lid - 1;

              if (face_dirs) {
                int fdir = (MR_FaceDir_i((MRegion_ptr)cell, i) == 1) ? 1 : -1;
                if (faceflip_[lid - 1]) fdir *= -1;
                (lface_dirs)[nf] = fdir;
              }

              MEnt_Mark(fadj, mkid);
              nf++;
            }
          }
          ++i;
        }
      }
      List_Delete(fedges0);

      /* Add the base face */
      lid = MEnt_ID(face0);
      lfaceids[nf] = lid - 1;

      if (face_dirs) {
        fdir0 = fdir0 ? 1 : -1;
        if (faceflip_[lid - 1]) fdir0 *= -1;
        (lface_dirs)[nf] = fdir0;
      }
      nf++;

      /* If there is a last remaining face, it is the top face */
      MFace_ptr fopp;
      idx = 0;
      int i = 0;
      while ((fopp = List_Next_Entry(rfaces, &idx))) {
        if (fopp != face0 && !MEnt_IsMarked(fopp, mkid)) {
          lid = MEnt_ID(fopp);
          lfaceids[nf] = lid - 1;

          if (face_dirs) {
            int fdir = (MR_FaceDir_i((MRegion_ptr)cell, i) == 1) ? 1 : -1;
            if (faceflip_[lid - 1]) fdir *= -1;
            (lface_dirs)[nf] = fdir;
          }
          nf++;
          break;
        }
        ++i;
      }

      List_Unmark(rfaces, mkid);
      MSTK_FreeMarker(mkid);
      List_Delete(rfaces);

      faceids = lfaceids;
      if (face_dirs) *face_dirs = lface_dirs;

    } else {
      getCellFacesAndDirs_unordered_(cellid, faceids, face_dirs);
    }
  } else {
    getCellFacesAndDirs_unordered_(cellid, faceids, face_dirs);
  }
}


void
Mesh_MSTK::getCellFacesAndDirs_unordered_(
  const Entity_ID cellid,
  View_type<const Entity_ID, MemSpace_kind::HOST>& faceids,
  View_type<const Direction_type, MemSpace_kind::HOST>* face_dirs) const
{
  Entity_ID_View lfaceids;
  Entity_Direction_View lface_dirs;
  MEntity_ptr cell = cell_id_to_handle_[cellid];

  if (getManifoldDimension() == 3) {
    int nrf;
    List_ptr rfaces;
    rfaces = MR_Faces((MRegion_ptr)cell);
    nrf = List_Num_Entries(rfaces);
    Kokkos::resize(lfaceids, nrf);

    for (int i = 0; i < nrf; ++i) {
      MFace_ptr face = List_Entry(rfaces, i);
      int lid = MEnt_ID(face);
      lfaceids[i] = lid - 1;
    }

    List_Delete(rfaces);

    /* Reserved for next major MSTK release
    int rfaceids[MAXPF3];

    MR_FaceIDs((MRegion_ptr)cell,&nrf,rfaceids);
    faceids->resize(nrf);
    Entity_ID_List::iterator it = faceids->begin();
    for (int i = 0; i < nrf; ++i) {
      *it = rfaceids[i]-1;
      ++it;
    }
    */

    if (face_dirs) {
      Kokkos::resize(lface_dirs, nrf);
      for (int i = 0; i < nrf; ++i) {
        int lid = lfaceids[i];
        int fdir = 2 * MR_FaceDir_i((MRegion_ptr)cell, i) - 1;
        fdir = faceflip_[lid] ? -fdir : fdir;
        (lface_dirs)[i] = fdir; // assign to next spot by dereferencing iterator
      }
    }

  } else { // getManifoldDimension() = 2; surface or 2D mesh
    int nfe;
    List_ptr fedges;
    fedges = MF_Edges((MFace_ptr)cell, 1, 0);
    nfe = List_Num_Entries(fedges);

    Kokkos::resize(lfaceids, nfe);
    for (int i = 0; i < nfe; ++i) {
      MEdge_ptr edge = List_Entry(fedges, i);
      int lid = MEnt_ID(edge);
      lfaceids[i] = lid - 1; // assign to next spot by dereferencing iterator
    }

    List_Delete(fedges);

    /* Reserved for next major MSTK release

    int fedgeids[MAXPV2];
    MF_EdgeIDs((MFace_ptr)cell,1,0,&nfe,fedgeids);

    faceids->resize(nfe);
    Entity_ID_List::iterator itf = faceids->begin();
    for (int i = 0; i < nfe; ++i) {
      *itf = fedgeids[i]-1;
      ++itf;
    }
    */

    if (face_dirs) {
      Kokkos::resize(lface_dirs, nfe);
      for (int i = 0; i < nfe; ++i) {
        int lid = lfaceids[i];
        int fdir = 2 * MF_EdgeDir_i((MFace_ptr)cell, i) - 1;
        fdir = faceflip_[lid] ? -fdir : fdir;
        (lface_dirs)[i] = fdir; // assign to next spot by dereferencing iterator
      }
    }
  }
  faceids = lfaceids;
  if (face_dirs) *face_dirs = lface_dirs;
}


void
Mesh_MSTK::getCellFacesAndDirs(
  const Entity_ID cellid,
  View_type<const Entity_ID, MemSpace_kind::HOST>& faceids,
  View_type<const Direction_type, MemSpace_kind::HOST>* const face_dirs) const
{
  AMANZI_ASSERT(faces_initialized_);
  if (cells_initialized_) {
    getCellFacesAndDirs_ordered_(cellid, faceids, face_dirs);
  } else {
    getCellFacesAndDirs_unordered_(cellid, faceids, face_dirs);
  }
}


void
Mesh_MSTK::getCellEdges(const Entity_ID cellid,
                        View_type<const Entity_ID, MemSpace_kind::HOST>& edgeids) const
{
  Entity_ID_View ledgeids;
  AMANZI_ASSERT(edges_initialized_);
  MEntity_ptr cell;
  cell = cell_id_to_handle_[cellid];

  if (getManifoldDimension() == 3) {
    int nre;
    List_ptr redges;

    redges = MR_Edges((MRegion_ptr)cell);
    nre = List_Num_Entries(redges);
    Kokkos::resize(ledgeids, nre);

    for (int i = 0; i < nre; ++i) {
      MEdge_ptr edge = List_Entry(redges, i);
      int lid = MEnt_ID(edge);
      ledgeids[i] = lid - 1; // assign to next spot by dereferencing iterator
    }

    List_Delete(redges);

    /* Reserved for next major MSTK release
    int redgeids[MAXPF3];

    MR_EdgeIDs((MRegion_ptr)cell,&nre,redgeids);
    edgeids->resize(nre);
    Entity_ID_List::iterator it = edgeids->begin();
    for (int i = 0; i < nre; ++i) {
      *it = redgeids[i]-1;
      ++it;
    }
    */

  } else { // manifold_dimension() = 2; surface or 2D mesh
    int nfe;

    List_ptr fedges;
    fedges = MF_Edges((MFace_ptr)cell, 1, 0);
    nfe = List_Num_Entries(fedges);

    Kokkos::resize(ledgeids, nfe);

    for (int i = 0; i < nfe; ++i) {
      MEdge_ptr edge = List_Entry(fedges, i);
      int lid = MEnt_ID(edge);
      ledgeids[i] = lid - 1; // assign to next spot by dereferencing iterator
    }

    List_Delete(fedges);

    /* Reserved for next major MSTK release

    int fedgeids[MAXPV2];
    MF_EdgeIDs((MFace_ptr)cell,1,0,&nfe,fedgeids);

    edgeids->resize(nfe);
    Entity_ID_List::iterator ite = edgeids->begin();
    for (int i = 0; i < nfe; ++i) {
      *ite = fedgeids[i]-1;
      ++ite;
    }
    */
  }
  edgeids = ledgeids;
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
void
Mesh_MSTK::getCellNodes(const Entity_ID cellid,
                        View_type<const Entity_ID, MemSpace_kind::HOST>& nodeids) const
{
  Entity_ID_View lnodeids;
  int nn, lid;
  MEntity_ptr cell = cell_id_to_handle_[cellid];

  if (getManifoldDimension() == 3) { // Volume mesh
    List_ptr rverts = MR_Vertices(cell);

    nn = List_Num_Entries(rverts);
    Kokkos::resize(lnodeids, nn);
    for (int i = 0; i < nn; ++i) {
      lid = MEnt_ID(List_Entry(rverts, i));
      lnodeids[i] = lid - 1; // assign to next spot by dereferencing iterator
    }
    List_Delete(rverts);

  } else { // Surface mesh
    List_ptr fverts = MF_Vertices(cell, 1, 0);
    nn = List_Num_Entries(fverts);
    Kokkos::resize(lnodeids, nn);
    for (int i = 0; i < nn; ++i) {
      lid = MEnt_ID(List_Entry(fverts, i));
      lnodeids[i] = lid - 1; // assign to next spot by dereferencing iterator
    }
    List_Delete(fverts);
  }
  nodeids = lnodeids;
}


void
Mesh_MSTK::getFaceEdgesAndDirs(
  const Entity_ID faceid,
  View_type<const Entity_ID, MemSpace_kind::HOST>& edgeids,
  View_type<const Direction_type, MemSpace_kind::HOST>* edge_dirs) const
{
  AMANZI_ASSERT(faces_initialized_);
  AMANZI_ASSERT(edges_initialized_);

  Entity_ID_View ledgeids;
  Entity_Direction_View ledge_dirs;

  MEntity_ptr face = face_id_to_handle_[faceid];
  if (getManifoldDimension() == 3) {
    int nfe;
    List_ptr fedges;

    fedges = MF_Edges((MFace_ptr)face, 1, 0);
    nfe = List_Num_Entries(fedges);
    Kokkos::resize(ledgeids, nfe);

    for (int i = 0; i < nfe; ++i) {
      MEdge_ptr edge = List_Entry(fedges, i);
      int lid = MEnt_ID(edge);
      ledgeids[i] = lid - 1; // assign to next spot by dereferencing iterator
    }

    List_Delete(fedges);

    /* Reserved for next major MSTK release
    int fedgeids[MAXPF3];

    MF_EdgeIDs((MFace_ptr)face,&nfe,fedgeids);
    fedgeids->resize(nfe);
    Entity_ID_List::iterator it = fedgeids->begin();
    for (int i = 0; i < nfe; ++i) {
      *it = fedgeids[i]-1;
      ++it;
    }
    */

    if (edge_dirs) {
      Kokkos::resize(ledge_dirs, nfe);

      for (int i = 0; i < nfe; ++i) {
        int lid = ledgeids[i];
        int edir = 2 * MF_EdgeDir_i((MFace_ptr)face, i) - 1;
        edir = edgeflip_[lid] ? -edir : edir;
        (ledge_dirs)[i] = edir; // assign to next spot by dereferencing iterator
      }
    }

  } else { // getManifoldDimension() = 2; surface or 2D mesh
    // face is same dimension as edge; just return the edge with a
    // direction of 1
    MEdge_ptr edge = (MEdge_ptr)face;

    Kokkos::resize(ledgeids, 1);
    ledgeids[0] = MEnt_ID(edge) - 1;

    if (edge_dirs) {
      Kokkos::resize(ledge_dirs, 1);
      (ledge_dirs)[0] = 1;
    }
  }
  edgeids = ledgeids;
  if (edge_dirs) *edge_dirs = ledge_dirs;
}


//---------------------------------------------------------
// Get nodes of face
// On a distributed mesh, all nodes (OWNED or GHOST) of the face
// are returned
// In 3D, the nodes of the face are returned in ccw order consistent
// with the face normal
// In 2D, nfnodes is 2
//---------------------------------------------------------
void
Mesh_MSTK::getFaceNodes(const Entity_ID faceid,
                        View_type<const Entity_ID, MemSpace_kind::HOST>& nodeids) const
{
  AMANZI_ASSERT(faces_initialized_);
  Entity_ID_View lnodeids;
  int nn, lid;
  MEntity_ptr genface = face_id_to_handle_[faceid];

  if (getManifoldDimension() == 3) { // Volume mesh
    int dir = !faceflip_[faceid];

    List_ptr fverts = MF_Vertices(genface, dir, 0);
    AMANZI_ASSERT(fverts != nullptr);

    nn = List_Num_Entries(fverts);
    Kokkos::resize(lnodeids, nn);

    for (int i = 0; i < nn; ++i) {
      lid = MEnt_ID(List_Entry(fverts, i));
      lnodeids[i] = lid - 1; // assign to next spot by dereferencing iterator
    }

    List_Delete(fverts);

  } else { // Surface mesh or 2D mesh
    Kokkos::resize(lnodeids, 2);
    if (faceflip_[faceid]) {
      lnodeids[0] = MEnt_ID(ME_Vertex(genface, 1)) - 1;
      lnodeids[1] = MEnt_ID(ME_Vertex(genface, 0)) - 1;
    } else {
      lnodeids[0] = MEnt_ID(ME_Vertex(genface, 0)) - 1;
      lnodeids[1] = MEnt_ID(ME_Vertex(genface, 1)) - 1;
    }
  }
  nodeids = lnodeids;
}


//---------------------------------------------------------
// Get nodes of an edge
//---------------------------------------------------------
void
Mesh_MSTK::getEdgeNodes(const Entity_ID edgeid,
                        View_type<const Entity_ID, MemSpace_kind::HOST>& nodes) const
{
  Entity_ID_View lnodes;
  AMANZI_ASSERT(edges_initialized_);
  MEdge_ptr edge = (MEdge_ptr)edge_id_to_handle_[edgeid];
  Kokkos::resize(lnodes, 2);
  if (edgeflip_[edgeid]) {
    lnodes[0] = MEnt_ID(ME_Vertex(edge, 1)) - 1;
    lnodes[1] = MEnt_ID(ME_Vertex(edge, 0)) - 1;
  } else {
    lnodes[0] = MEnt_ID(ME_Vertex(edge, 0)) - 1;
    lnodes[1] = MEnt_ID(ME_Vertex(edge, 1)) - 1;
  }
  nodes = lnodes;
}


//---------------------------------------------------------
// Cells of type 'ptype' connected to a node. This routine uses
// push_back on or near the partition boundary since we cannot tell at
// the outset how many entries will be put into the list
//---------------------------------------------------------
void
Mesh_MSTK::getNodeCells(const Entity_ID nodeid,
                        const Parallel_kind ptype,
                        View_type<const Entity_ID, MemSpace_kind::HOST>& cellids) const
{
  Entity_ID_View lcellids;
  int idx, lid, nc;
  List_ptr cell_list;
  MEntity_ptr ment;

  MVertex_ptr mv = (MVertex_ptr)vtx_id_to_handle_[nodeid];

  /* Reserved for next major release of MSTK
  if (MV_PType(mv) == PINTERIOR && ptype != Parallel_type::GHOST) {

    if (manifold_dimension() == 3) {
      int nvr, regionids[200];
      MV_RegionIDs(mv,&nvr,regionids);
      AMANZI_ASSERT(nvr < 200);
      cellids->resize(nvr);
      Entity_ID_List::iterator it = cellids->begin();
      for (int i = 0; i < nvr; ++i) {
        *it = regionids[i]-1;  // assign to next spot by dereferencing iterator
        ++it;
      }
    }
    else {
      int nvf, faceids[200];
      MV_FaceIDs(mv,&nvf,faceids);
      AMANZI_ASSERT(nvf < 200);
      cellids->resize(nvf);
      Entity_ID_List::iterator it = cellids->begin();
      for (int i = 0; i < nvf; ++i) {
        *it = faceids[i]-1;  // assign to next spot by dereferencing iterator
        ++it;
      }
    }

  }
  else {
  */
  // mesh vertex on a processor boundary may be connected to owned
  // and ghost cells. So depending on the requested cell type, we
  // may have to omit some entries
  if (getManifoldDimension() == 3)
    cell_list = MV_Regions(mv);
  else
    cell_list = MV_Faces(mv);

  nc = List_Num_Entries(cell_list);
  Kokkos::resize(lcellids, nc); // resize to maximum size possible

  int n = 0;
  idx = 0;
  while ((ment = List_Next_Entry(cell_list, &idx))) {
    if (MEnt_PType(ment) == PGHOST) {
      if (ptype == Parallel_kind::GHOST || ptype == Parallel_kind::ALL) {
        lid = MEnt_ID(ment);
        lcellids[n++] = lid - 1;
      }
    } else {
      if (ptype == Parallel_kind::OWNED || ptype == Parallel_kind::ALL) {
        lid = MEnt_ID(ment);
        lcellids[n++] = lid - 1;
      }
    }
  }
  Kokkos::resize(lcellids, n); // resize to the actual number of cells being returned
  List_Delete(cell_list);
  cellids = lcellids;
}


//---------------------------------------------------------
// Faces of type 'ptype' connected to a node. This routine uses
// push_back on or near the partition boundary since we cannot tell at
// the outset how many entries will be put into the list
//---------------------------------------------------------
void
Mesh_MSTK::getNodeFaces(const Entity_ID nodeid,
                        const Parallel_kind ptype,
                        View_type<const Entity_ID, MemSpace_kind::HOST>& faceids) const
{
  Entity_ID_View lfaceids;
  int idx, lid, n;
  List_ptr face_list;
  MEntity_ptr ment;

  AMANZI_ASSERT(faces_initialized_);
  MVertex_ptr mv = (MVertex_ptr)vtx_id_to_handle_[nodeid];

  /* Reserved for next major release of MSTK
  if (MV_PType(mv) == PINTERIOR && ptype != Parallel_type::GHOST) {
    if (manifold_dimension() == 3) {
      int nvf, vfaceids[200];

      MV_FaceIDs(mv,&nvf,vfaceids);
      AMANZI_ASSERT(nvf < 200);

      faceids->resize(nvf);
      Entity_ID_List::iterator it = faceids->begin();
      for (int i = 0; i < nvf; ++i) {
        *it = vfaceids[i]-1;  // assign to next spot by dereferencing iterator
        ++it;
      }
    }
    else if (manifold_dimension() == 2) {
      int nve, vedgeids[200];

      MV_EdgeIDs(mv,&nve,vedgeids);
      AMANZI_ASSERT(nve < 200);

      faceids->resize(nve);
      Entity_ID_List::iterator it = faceids->begin();
      for (int i = 0; i < nve; ++i) {
        *it = vedgeids[i]-1;  // assign to next spot by dereferencing iterator
        ++it;
      }
    }
  }
  else {
  */

  if (getManifoldDimension() == 3)
    face_list = MV_Faces(mv);
  else
    face_list = MV_Edges(mv);

  int nf = List_Num_Entries(face_list);
  Kokkos::resize(lfaceids, nf); // resize to the maximum
  idx = 0;
  n = 0;
  while ((ment = List_Next_Entry(face_list, &idx))) {
    if (MEnt_PType(ment) == PGHOST) {
      if (ptype == Parallel_kind::GHOST || ptype == Parallel_kind::ALL) {
        lid = MEnt_ID(ment);
        lfaceids[n++] = lid - 1;
      }
    } else {
      if (ptype == Parallel_kind::OWNED || ptype == Parallel_kind::ALL) {
        lid = MEnt_ID(ment);
        lfaceids[n++] = lid - 1;
      }
    }
  }
  Kokkos::resize(lfaceids, n); // resize to the actual number of faces being returned
  faceids = lfaceids;
  List_Delete(face_list);
}


//---------------------------------------------------------
// Edges of type 'ptype' connected to a node.
//---------------------------------------------------------
void
Mesh_MSTK::getNodeEdges(const Entity_ID nodeid,
                        const Parallel_kind ptype,
                        View_type<const Entity_ID, MemSpace_kind::HOST>& edgeids) const
{
  Entity_ID_View ledgeids;
  int idx, lid, nc;
  List_ptr edge_list;
  MEntity_ptr ment;

  MVertex_ptr mv = (MVertex_ptr)vtx_id_to_handle_[nodeid];

  // mesh vertex on a processor boundary may be connected to owned
  // and ghost cells. So depending on the requested cell type, we
  // may have to omit some entries
  edge_list = MV_Edges(mv);
  nc = List_Num_Entries(edge_list);

  Kokkos::resize(ledgeids, nc); // resize to maximum size possible

  int n = 0;
  idx = 0;
  while ((ment = List_Next_Entry(edge_list, &idx))) {
    if (MEnt_PType(ment) == PGHOST) {
      if (ptype == Parallel_kind::GHOST || ptype == Parallel_kind::ALL) {
        lid = MEnt_ID(ment);
        ledgeids[n++] = lid - 1;
      }
    } else {
      if (ptype == Parallel_kind::OWNED || ptype == Parallel_kind::ALL) {
        lid = MEnt_ID(ment);
        ledgeids[n++] = lid - 1;
      }
    }
  }
  Kokkos::resize(ledgeids, n); // resize to the actual number of cells being returned
  edgeids = ledgeids;
  List_Delete(edge_list);
}


//---------------------------------------------------------
// Faces of type 'ptype' connected to an edge.
//---------------------------------------------------------
void
Mesh_MSTK::getEdgeFaces(const Entity_ID edgeid,
                        const Parallel_kind ptype,
                        View_type<const Entity_ID, MemSpace_kind::HOST>& faceids) const
{
  Entity_ID_View lfaceids;
  int idx, lid, nc;
  List_ptr face_list;
  MEntity_ptr ment;

  //AMANZI_ASSERT(getManifoldDimension() == 3);

  MEdge_ptr me = (MEdge_ptr)edge_id_to_handle_[edgeid];
  face_list = ME_Faces(me);

  nc = List_Num_Entries(face_list);
  Kokkos::resize(lfaceids, nc); // resize to maximum size possible

  int n = 0;
  idx = 0;
  while ((ment = List_Next_Entry(face_list, &idx))) {
    if (MEnt_PType(ment) == PGHOST) {
      if (ptype == Parallel_kind::GHOST || ptype == Parallel_kind::ALL) {
        lid = MEnt_ID(ment);
        lfaceids[n++] = lid - 1;
      }

    } else {
      if (ptype == Parallel_kind::OWNED || ptype == Parallel_kind::ALL) {
        lid = MEnt_ID(ment);
        lfaceids[n++] = lid - 1;
      }
    }
  }
  Kokkos::resize(lfaceids, n); // resize to the actual number of cells being returned
  faceids = lfaceids;
  List_Delete(face_list);
}


//---------------------------------------------------------
// Cells of type 'ptype' connected to an edge. This routine uses
// push_back on or near the partition boundary since we cannot tell at
// the outset how many entries will be put into the list
//---------------------------------------------------------
void
Mesh_MSTK::getEdgeCells(const Entity_ID edgeid,
                        const Parallel_kind ptype,
                        View_type<const Entity_ID, MemSpace_kind::HOST>& cellids) const
{
  Entity_ID_View lcellids;
  MEdge_ptr me = (MEdge_ptr)edge_id_to_handle_[edgeid];

  // mesh edge on a processor boundary may be connected to owned
  // and ghost cells. So depending on the requested cell type, we
  // may have to omit some entries
  List_ptr cell_list;
  if (getManifoldDimension() == 3)
    cell_list = ME_Regions(me);
  else
    cell_list = ME_Faces(me);

  int nc = List_Num_Entries(cell_list);
  Kokkos::resize(lcellids, nc); // resize to maximum size possible

  int n = 0;
  int idx = 0;
  MEntity_ptr ment;
  while ((ment = List_Next_Entry(cell_list, &idx))) {
    if (MEnt_PType(ment) == PGHOST) {
      if (ptype == Parallel_kind::GHOST || ptype == Parallel_kind::ALL) {
        int lid = MEnt_ID(ment);
        lcellids[n++] = lid - 1;
      }

    } else {
      if (ptype == Parallel_kind::OWNED || ptype == Parallel_kind::ALL) {
        int lid = MEnt_ID(ment);
        lcellids[n++] = lid - 1;
      }
    }
  }
  Kokkos::resize(lcellids, n); // resize to the actual number of cells being returned
  cellids = lcellids;
  List_Delete(cell_list);
}


//---------------------------------------------------------
// Cells connected to a face
//---------------------------------------------------------
void
Mesh_MSTK::getFaceCells(const Entity_ID faceid,
                        const Parallel_kind ptype,
                        View_type<const Entity_ID, MemSpace_kind::HOST>& cellids) const
{
  AMANZI_ASSERT(faces_initialized_);
  Entity_ID_List vcellids;

  if (getManifoldDimension() == 3) {
    MFace_ptr mf = (MFace_ptr)face_id_to_handle_[faceid];

    List_ptr fregs = MF_Regions(mf);
    MRegion_ptr mr;
    if (ptype == Parallel_kind::ALL) {
      int idx = 0;
      while ((mr = List_Next_Entry(fregs, &idx))) vcellids.push_back(MR_ID(mr) - 1);

    } else {
      int idx = 0;
      while ((mr = List_Next_Entry(fregs, &idx))) {
        if (MEnt_PType(mr) == PGHOST) {
          if (ptype == Parallel_kind::GHOST) vcellids.push_back(MR_ID(mr) - 1);
        } else if (ptype == Parallel_kind::OWNED) {
          vcellids.push_back(MR_ID(mr) - 1);
        }
      }
    }
    List_Delete(fregs);

  } else {
    MEdge_ptr me = (MEdge_ptr)face_id_to_handle_[faceid];

    List_ptr efaces = ME_Faces(me);
    MFace_ptr mf;
    if (ptype == Parallel_kind::ALL) {
      int idx = 0;
      while ((mf = List_Next_Entry(efaces, &idx))) vcellids.push_back(MF_ID(mf) - 1);
    } else {
      int idx = 0;
      while ((mf = List_Next_Entry(efaces, &idx))) {
        if (MEnt_PType(mf) == PGHOST) {
          if (ptype == Parallel_kind::GHOST) vcellids.push_back(MF_ID(mf) - 1);
        } else if (ptype == Parallel_kind::OWNED) {
          vcellids.push_back(MF_ID(mf) - 1);
        }
      }
    }
    List_Delete(efaces);
  }
  vectorToConstView(cellids, vcellids);
}


//---------------------------------------------------------
// Node coordinates - 3 in 3D and 2 in 2D
//---------------------------------------------------------
AmanziGeometry::Point
Mesh_MSTK::getNodeCoordinate(const Entity_ID nodeid) const
{
  MEntity_ptr vtx = vtx_id_to_handle_[nodeid];
  double coords[3];
  MV_Coords(vtx, coords);
  if (getSpaceDimension() == 2) {
    return AmanziGeometry::Point(coords[0], coords[1]);
  } else {
    return AmanziGeometry::Point(coords[0], coords[1], coords[2]);
  }
}


//---------------------------------------------------------
// Modify a node's coordinates
//---------------------------------------------------------
void
Mesh_MSTK::setNodeCoordinate(const AmanziMesh::Entity_ID nodeid,
                             const AmanziGeometry::Point& coords)
{
  MVertex_ptr v = vtx_id_to_handle_[nodeid];
  double coordarray[3] = { 0.0, 0.0, 0.0 };
  for (int i = 0; i < getSpaceDimension(); i++) coordarray[i] = coords[i];
  MV_Set_Coords(v, (double*)coordarray);
}


//---------------------------------------------------------
// Private routine creating mesh sets for GM regions
//---------------------------------------------------------
MSet_ptr
Mesh_MSTK::build_set(const Teuchos::RCP<const AmanziGeometry::Region>& region,
                     const Entity_kind kind) const
{
  int celldim = manifold_dimension();
  int space_dim = space_dimension();
  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm = geometric_model();

  // Modify region/set name by prefixing it with the type of entity requested

  std::string internal_name = internal_name_of_set(region, kind);

  // Create entity set based on the region defintion

  MSet_ptr mset;
  MType enttype;
  switch (kind) {
  case CELL: // cellsets

    enttype = (celldim == 3) ? MREGION : MFACE;
    mset = MSet_New(mesh_, internal_name.c_str(), enttype);

    if (region->get_type() == AmanziGeometry::RegionType::BOX ||
        region->get_type() == AmanziGeometry::RegionType::CYLINDER ||
        region->get_type() == AmanziGeometry::RegionType::COLORFUNCTION) {
      int ncell = num_entities(CELL, Parallel_type::ALL);

      for (int icell = 0; icell < ncell; icell++)
        if (region->inside(cell_centroid(icell))) MSet_Add(mset, cell_id_to_handle[icell]);

    } else if (region->get_type() == AmanziGeometry::RegionType::ALL) {
      int ncell = num_entities(CELL, Parallel_type::ALL);

      for (int icell = 0; icell < ncell; icell++) MSet_Add(mset, cell_id_to_handle[icell]);

    } else if (region->get_type() == AmanziGeometry::RegionType::ENUMERATED) {
      auto rgn = Teuchos::rcp_static_cast<const AmanziGeometry::RegionEnumerated>(region);
      int ncell = num_entities(CELL, Parallel_type::ALL);

      for (int icell = 0; icell < ncell; icell++) {
        Entity_ID gid = MEnt_GlobalID(cell_id_to_handle[icell]);
        for (const auto& jset : rgn->entities()) {
          if (jset == gid) {
            MSet_Add(mset, cell_id_to_handle[icell]);
            break;
          }
        }
      }

    } else if (region->get_type() == AmanziGeometry::RegionType::POINT) {
      AmanziGeometry::Point vpnt(space_dim);
      AmanziGeometry::Point rgnpnt(space_dim);

      rgnpnt = Teuchos::rcp_static_cast<const AmanziGeometry::RegionPoint>(region)->point();

      int nnode = num_entities(NODE, Parallel_type::ALL);
      double mindist2 = 1.e+16;
      int minnode = -1;

      int inode;
      for (inode = 0; inode < nnode; inode++) {
        node_get_coordinates(inode, &vpnt);
        double dist2 = (vpnt - rgnpnt) * (vpnt - rgnpnt);

        if (dist2 < mindist2) {
          mindist2 = dist2;
          minnode = inode;
          if (mindist2 <= 1.0e-14) break;
        }
      }

      Entity_ID_List cells;
      node_get_cells(minnode, Parallel_type::ALL, &cells);

      int ncells = cells.size();
      for (int ic = 0; ic < ncells; ic++) {
        Entity_ID icell = cells[ic];

        // Check if point is contained in cell
        if (point_in_cell(rgnpnt, icell)) MSet_Add(mset, cell_id_to_handle[icell]);
      }

      // finally check all cells, typical for anisotropic meshes
      if (MSet_Num_Entries(mset) == 0) {
        int ncells_wghost = num_entities(CELL, Parallel_type::ALL);
        for (int c = 0; c < ncells_wghost; ++c)
          if (point_in_cell(rgnpnt, c)) MSet_Add(mset, cell_id_to_handle[c]);
      }

    } else if ((region->get_type() == AmanziGeometry::RegionType::PLANE) ||
               (region->get_type() == AmanziGeometry::RegionType::POLYGON)) {
      if (celldim == 2) {
        Entity_ID_List nodes;
        AmanziGeometry::Point xyz(space_dim);

        int ncells = num_entities(CELL, Parallel_type::ALL);
        for (int ic = 0; ic < ncells; ic++) {
          cell_get_nodes(ic, &nodes);

          bool on_plane = true;
          for (int j = 0; j < nodes.size(); ++j) {
            node_get_coordinates(nodes[j], &xyz);
            if (!region->inside(xyz)) {
              on_plane = false;
              break;
            }
          }

          if (on_plane) MSet_Add(mset, cell_id_to_handle[ic]);
        }
      }

    } else if (region->get_type() == AmanziGeometry::RegionType::LOGICAL) {
      // will process later in this subroutine
    } else if (region->get_type() == AmanziGeometry::RegionType::LABELEDSET) {
      if (parent_mesh_.get() != NULL) {
        auto lsrgn = Teuchos::rcp_static_cast<const AmanziGeometry::RegionLabeledSet>(region);
        std::string label = lsrgn->label();
        std::string entity_type = lsrgn->entity_str();

        if (kind == CELL && entity_type != "FACE") {
          if (vo_.get() && vo_->os_OK(Teuchos::VERB_MEDIUM)) {
            *(vo_->os()) << "Found labeled set region \"" << region->get_name()
                         << "\" but it contains entities of type " << entity_type
                         << ", not the requested type.\n";
          }
        } else {
          MSet_ptr mset2;

          // Build set on a fly.
          // -- first request the set on the parent to make sure it was constructed in MSTK in all cases.
          AmanziMesh::Entity_ID_List parent_ids;
          parent_mesh_->get_set_entities(region->get_name(), FACE, Parallel_type::ALL, &parent_ids);

          int ival;
          double rval;
          void* pval = nullptr;
          MAttrib_ptr att = (manifold_dimension() == 3) ? rparentatt : fparentatt;

          std::string internal_parent_name = internal_name_of_set(region, FACE);
          mset2 = MESH_MSetByName(parent_mesh_->mesh_, internal_parent_name.c_str());

          for (int c = 0; c < num_entities(CELL, Parallel_type::ALL); ++c) {
            auto ment = cell_id_to_handle[c];
            MEnt_Get_AttVal(ment, att, &ival, &rval, &pval);
            if (MSet_Locate(mset2, pval) >= 0) MSet_Add(mset, ment);
          }

          // Due to the parallel partitioning its possible that this
          // set is not on this processor
          if (!mset) {
            if (comm_->NumProc() == 1) {
              Errors::Message msg;
              msg << "Could not find labeled set \"" << label
                  << "\" in mesh file to initialize mesh set \"" << region->get_name()
                  << "\". Verify mesh file.";
              amanzi_throw(msg);
            }
          }
        }
      } else {
        // Just retrieve and return the set

        auto lsrgn = Teuchos::rcp_static_cast<const AmanziGeometry::RegionLabeledSet>(region);
        std::string label = lsrgn->label();
        std::string entity_type = lsrgn->entity_str();

        if (entity_type != "CELL") {
          Errors::Message mesg;
          mesg << "Entity type of labeled set region \"" << region->get_name()
               << "\" and build_set request (cell) do not match";
          Exceptions::amanzi_throw(mesg);
        }

        mset = MESH_MSetByName(mesh_, internal_name.c_str());

        std::string other_internal_name = other_internal_name_of_set(region, kind);
        MSet_ptr mset2 = MESH_MSetByName(mesh_, other_internal_name.c_str());

        if (mset) {
          if (mset2) {
            std::stringstream mesg_stream;
            mesg_stream << "Exodus II file has element block and element set with the same ID "
                        << label << " - Amanzi cannot handle this case.";
            Errors::Message mesg(mesg_stream.str());
            Exceptions::amanzi_throw(mesg);
          }
        } else {
          if (mset2)
            mset = mset2;
          else {
            std::stringstream mesg_stream;
            mesg_stream << "Exodus II file has no labeled cell set with ID " << label;
            Errors::Message mesg(mesg_stream.str());
            Exceptions::amanzi_throw(mesg);
          }
        }
      }
    } else {
      if (vo_.get() && vo_->os_OK(Teuchos::VERB_HIGH)) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *(vo_->os()) << "Requested CELLS on region " << region->get_name() << " of type "
                     << region->get_type() << " and dimension " << region->get_manifold_dimension()
                     << ".\n"
                     << "This request will result in an empty set";
      }
    }

    break;

  case FACE: // sidesets

    //
    // Commented out so that we can ask for a face set in a 3D box
    //
    //          if (region->get_manifold_dimension() != celldim-1)
    //            {
    //              std::cerr << "No region of dimension " << celldim-1 << " defined in geometric model" << std::endl;
    //              std::cerr << "Cannot construct cell set from region " << setname << std::endl;
    //            }

    enttype = (celldim == 3) ? MFACE : MEDGE;
    mset = MSet_New(mesh_, internal_name.c_str(), enttype);

    if (region->get_type() == AmanziGeometry::RegionType::BOX ||
        region->get_type() == AmanziGeometry::RegionType::CYLINDER) {
      int nface = num_entities(FACE, Parallel_type::ALL);

      if (nface > 0) {
        if (!kdtree_faces_initialized_) {
          face_centroid(0);
          kdtree_faces_.Init(&face_centroids_);
          kdtree_faces_initialized_ = true;
        }

        auto box = Teuchos::rcp_dynamic_cast<const AmanziGeometry::RegionBox>(region);
        AmanziGeometry::Point query = (box->point0() + box->point1()) / 2;
        double radius = AmanziGeometry::norm(box->point0() - query);
        double radius_sqr = std::pow(radius + AmanziGeometry::TOL, 2);

        std::vector<double> dist_sqr;
        auto idx = kdtree_faces_.SearchInSphere(query, dist_sqr, radius_sqr);

        for (int i = 0; i < idx.size(); ++i) {
          int iface = idx[i];
          if (region->inside(face_centroid(iface))) MSet_Add(mset, face_id_to_handle[iface]);
        }
      }
    } else if (region->get_type() == AmanziGeometry::RegionType::ALL) {
      int nface = num_entities(FACE, Parallel_type::ALL);

      for (int iface = 0; iface < nface; iface++) { MSet_Add(mset, face_id_to_handle[iface]); }
    } else if (region->get_type() == AmanziGeometry::RegionType::PLANE ||
               region->get_type() == AmanziGeometry::RegionType::POLYGON) {
      int nface = num_entities(FACE, Parallel_type::ALL);

      for (int iface = 0; iface < nface; iface++) {
        std::vector<AmanziGeometry::Point> fcoords(space_dim);

        face_get_coordinates(iface, &fcoords);

        bool on_plane = true;
        for (int j = 0; j < fcoords.size(); ++j) {
          if (!region->inside(fcoords[j])) {
            on_plane = false;
            break;
          }
        }

        if (on_plane) MSet_Add(mset, face_id_to_handle[iface]);
      }

    } else if (region->get_type() == AmanziGeometry::RegionType::LABELEDSET) {
      // Just retrieve and return the set

      Teuchos::RCP<const AmanziGeometry::RegionLabeledSet> lsrgn =
        Teuchos::rcp_static_cast<const AmanziGeometry::RegionLabeledSet>(region);
      std::string label = lsrgn->label();
      std::string entity_type = lsrgn->entity_str();

      if (entity_type != "FACE") {
        Errors::Message mesg;
        mesg << "Entity type of labeled set region \"" << region->get_name()
             << "\" and build_set request (face) do not match";
        Exceptions::amanzi_throw(mesg);
      }

      mset = MESH_MSetByName(mesh_, internal_name.c_str());
    } else if (region->get_type() == AmanziGeometry::RegionType::LOGICAL) {
      // Will handle it later in the routine
    } else if (region->get_type() == AmanziGeometry::RegionType::BOUNDARY) {
      const Epetra_Map& fmap = face_map(true);
      const Epetra_Map& map = exterior_face_map(true);

      int nface = map.NumMyElements();

      for (int iface = 0; iface < nface; iface++) {
        int lid = fmap.LID(map.GID(iface));
        MSet_Add(mset, face_id_to_handle[lid]);
      }
    } else {
      Teuchos::RCP<const VerboseObject> verbobj = verbosity_obj();
      if (verbobj.get() && verbobj->os_OK(Teuchos::VERB_HIGH)) {
        Teuchos::OSTab tab = verbobj->getOSTab();
        *(verbobj->os()) << "Requested FACES on region " << region->get_name() << " of type "
                         << region->get_type() << " and dimension "
                         << region->get_manifold_dimension() << ".\n"
                         << "This request will result in an empty set\n";
      }
    }
    break;

  case EDGE: // Edgesets

    enttype = MEDGE;
    mset = MSet_New(mesh_, internal_name.c_str(), enttype);

    if (region->get_type() == AmanziGeometry::RegionType::BOX ||
        region->get_type() == AmanziGeometry::RegionType::PLANE ||
        region->get_type() == AmanziGeometry::RegionType::POLYGON) {
      int nedge = num_entities(EDGE, Parallel_type::ALL);

      for (int iedge = 0; iedge < nedge; iedge++) {
        const auto& epnt = edge_centroid(iedge);

        if (region->inside(epnt)) { MSet_Add(mset, edge_id_to_handle[iedge]); }
      }
    } else if (region->get_type() == AmanziGeometry::RegionType::LOGICAL) {
      // We will handle it later in the routine
    } else {
      Teuchos::RCP<const VerboseObject> verbobj = verbosity_obj();
      if (verbobj.get() && verbobj->os_OK(Teuchos::VERB_HIGH)) {
        Teuchos::OSTab tab = verbobj->getOSTab();
        *(verbobj->os()) << "Requested EDGEs on region " << region->get_name() << " of type "
                         << region->get_type() << " and dimension "
                         << region->get_manifold_dimension() << ".\n"
                         << "This request will result in an empty set\n";
      }
    }

    break;

  case NODE: // Nodesets

    enttype = MVERTEX;
    mset = MSet_New(mesh_, internal_name.c_str(), enttype);

    if (region->get_type() == AmanziGeometry::RegionType::BOX ||
        region->get_type() == AmanziGeometry::RegionType::PLANE ||
        region->get_type() == AmanziGeometry::RegionType::POLYGON ||
        region->get_type() == AmanziGeometry::RegionType::CYLINDER ||
        region->get_type() == AmanziGeometry::RegionType::POINT) {
      int nnode = num_entities(NODE, Parallel_type::ALL);

      for (int inode = 0; inode < nnode; inode++) {
        AmanziGeometry::Point vpnt(space_dim);
        node_get_coordinates(inode, &vpnt);

        if (region->inside(vpnt)) {
          MSet_Add(mset, vtx_id_to_handle[inode]);

          // Only one node per point region
          if (region->get_type() == AmanziGeometry::RegionType::POINT) break;
        }
      }
    } else if (region->get_type() == AmanziGeometry::RegionType::ALL) {
      int nnode = num_entities(NODE, Parallel_type::ALL);

      for (int inode = 0; inode < nnode; inode++) MSet_Add(mset, vtx_id_to_handle[inode]);

    } else if (region->get_type() == AmanziGeometry::RegionType::LABELEDSET) {
      // Just retrieve and return the set

      Teuchos::RCP<const AmanziGeometry::RegionLabeledSet> lsrgn =
        Teuchos::rcp_static_cast<const AmanziGeometry::RegionLabeledSet>(region);
      std::string label = lsrgn->label();
      std::string entity_type = lsrgn->entity_str();

      if (entity_type != "FACE") {
        Errors::Message mesg;
        mesg << "Entity type of labeled set region \"" << region->get_name()
             << "\" and build_set request (face) do not match";
        Exceptions::amanzi_throw(mesg);
      }

      mset = MESH_MSetByName(mesh_, internal_name.c_str());
    } else if (region->get_type() == AmanziGeometry::RegionType::LOGICAL) {
      // We will handle it later in the routine
    } else if (region->get_type() == AmanziGeometry::RegionType::BOUNDARY) {
      const Epetra_Map& vmap = node_map(true);
      const Epetra_Map& map = exterior_node_map(true);

      int nnode = map.NumMyElements();

      for (int inode = 0; inode < nnode; inode++) {
        int lid = vmap.LID(map.GID(inode));
        MSet_Add(mset, vtx_id_to_handle[lid]);
      }
    } else {
      Teuchos::RCP<const VerboseObject> verbobj = verbosity_obj();
      if (verbobj.get() && verbobj->os_OK(Teuchos::VERB_HIGH)) {
        Teuchos::OSTab tab = verbobj->getOSTab();
        *(verbobj->os()) << "Requested POINTS on region " << region->get_name() << " of type "
                         << region->get_type() << " and dimension "
                         << region->get_manifold_dimension() << ".\n"
                         << "This request will result in an empty set\n";
      }
    }

    break;

  default:
    break;
  }


  if (region->get_type() == AmanziGeometry::RegionType::LOGICAL) {
    Teuchos::RCP<const AmanziGeometry::RegionLogical> boolregion =
      Teuchos::rcp_static_cast<const AmanziGeometry::RegionLogical>(region);
    const std::vector<std::string> region_names = boolregion->get_component_regions();
    int nreg = region_names.size();

    std::vector<MSet_ptr> msets;
    std::vector<Teuchos::RCP<const AmanziGeometry::Region>> regions;

    for (int r = 0; r < nreg; r++) {
      Teuchos::RCP<const AmanziGeometry::Region> rgn1 = gm->FindRegion(region_names[r]);
      regions.push_back(rgn1);

      // Did not find the region
      if (rgn1 == Teuchos::null) {
        std::stringstream mesg_stream;
        mesg_stream << "Geometric model has no region named " << region_names[r];
        Errors::Message mesg(mesg_stream.str());
        Exceptions::amanzi_throw(mesg);
      }

      internal_name = internal_name_of_set(rgn1, kind);
      MSet_ptr mset1 = MESH_MSetByName(mesh_, internal_name.c_str());
      if (!mset1) mset1 = build_set(rgn1, kind); // Recursive call

      msets.push_back(mset1);
    }

    // Check the entity types of the sets are consistent with the
    // entity type of the requested set

    for (int ms = 0; ms < msets.size(); ms++)
      if (MSet_EntDim(msets[ms]) != enttype) {
        // Validate the dimensionality of the object
        if (MSet_EntDim(msets[ms]) == space_dim) {
          //MSet_ptr testerr = construct_logical(region,gm,kind,enttype);
          continue;
        }

        Errors::Message mesg("Amanzi cannot operate on sets of different entity types");
        Exceptions::amanzi_throw(mesg);
      }

    int mkid = MSTK_GetMarker();

    if (boolregion->get_operation() == AmanziGeometry::BoolOpType::COMPLEMENT) {
      for (int ms = 0; ms < msets.size(); ms++) MSet_Mark(msets[ms], mkid);

      int idx = 0;
      switch (enttype) {
      case MREGION:
        MRegion_ptr mr;
        while ((mr = MESH_Next_Region(mesh_, &idx)))
          if (!MEnt_IsMarked(mr, mkid)) MSet_Add(mset, mr);
        break;
      case MFACE:
        MFace_ptr mf;
        while ((mf = MESH_Next_Face(mesh_, &idx)))
          if (!MEnt_IsMarked(mf, mkid)) MSet_Add(mset, mf);
        break;
      case MEDGE:
        MEdge_ptr me;
        while ((me = MESH_Next_Edge(mesh_, &idx)))
          if (!MEnt_IsMarked(me, mkid)) MSet_Add(mset, me);
      case MVERTEX:
        MVertex_ptr mv;
        while ((mv = MESH_Next_Vertex(mesh_, &idx)))
          if (!MEnt_IsMarked(mv, mkid)) MSet_Add(mset, mv);
        break;
      default:
        break;
      }

      for (int ms = 0; ms < msets.size(); ms++) MSet_Unmark(msets[ms], mkid);

    } else if (boolregion->get_operation() == AmanziGeometry::BoolOpType::UNION) {
      for (int ms = 0; ms < msets.size(); ms++) {
        MEntity_ptr ment;
        int idx = 0;
        while ((ment = MSet_Next_Entry(msets[ms], &idx)))
          if (!MEnt_IsMarked(ment, mkid)) {
            MSet_Add(mset, ment);
            MEnt_Mark(ment, mkid);
          }
      }
      MSet_Unmark(mset, mkid);

    } else if (boolregion->get_operation() == AmanziGeometry::BoolOpType::SUBTRACT) {
      /* Mark entities in all sets except the first */

      for (int ms = 1; ms < msets.size(); ms++) MSet_Mark(msets[ms], mkid);

      /* Look for entities in the first set but not in
         any of the other sets */
      MEntity_ptr ment;
      int idx = 0;
      while ((ment = MSet_Next_Entry(msets[0], &idx)))
        if (!MEnt_IsMarked(ment, mkid)) {
          MSet_Add(mset, ment);
          MEnt_Mark(ment, mkid);
        }

      for (int ms = 1; ms < msets.size(); ms++) MSet_Unmark(msets[ms], mkid);

    } else if (boolregion->get_operation() == AmanziGeometry::BoolOpType::INTERSECT) {
      /* Can't do this using markers alone - need attributes */

      MAttrib_ptr matt = MAttrib_New(mesh_, "XSECTATT", INT, MALLTYPE);

      for (int ms = 0; ms < msets.size(); ms++) {
        MEntity_ptr ment;
        int idx = 0;
        while ((ment = MSet_Next_Entry(msets[ms], &idx))) {
          int ival;
          double rval;
          void* pval;
          MEnt_Get_AttVal(ment, matt, &ival, &rval, &pval);
          ival++;
          MEnt_Set_AttVal(ment, matt, ival, rval, pval);
        }
      }

      for (int ms = 0; ms < msets.size(); ms++) {
        MEntity_ptr ment;
        int idx = 0;
        while ((ment = MSet_Next_Entry(msets[ms], &idx))) {
          int ival;
          double rval;
          void* pval;
          MEnt_Get_AttVal(ment, matt, &ival, &rval, &pval);
          if ((ival == msets.size()) && !MEnt_IsMarked(ment, mkid)) {
            /* entity is in all sets */
            MSet_Add(mset, ment);
            MEnt_Mark(ment, mkid);
          }
        }
      }

      MSet_Unmark(mset, mkid);

      for (int ms = 0; ms < msets.size(); ms++) {
        MEntity_ptr ment;
        int idx = 0;
        while ((ment = MSet_Next_Entry(msets[ms], &idx))) MEnt_Rem_AttVal(ment, matt);
      }
      MAttrib_Delete(matt);
    }

    MSTK_FreeMarker(mkid);

    for (int ms = 0; ms < msets.size(); ms++) {
      MSet_Unmark(msets[ms], mkid);
      if (regions[ms]->get_lifecycle() == AmanziGeometry::LifeCycleType::TEMPORARY)
        MSet_Delete(msets[ms]);
    }
  }
  return mset;
}


//---------------------------------------------------------
// Get list of entities of type 'category' in set specified by setname
//---------------------------------------------------------
void
Mesh_MSTK::get_set_entities_and_vofs(const std::string& setname,
                                     const Entity_kind kind,
                                     const Parallel_type ptype,
                                     std::vector<Entity_ID>* setents,
                                     std::vector<double>* vofs) const
{
  int idx;
  MSet_ptr mset1(NULL);
  MEntity_ptr ment;
  int celldim = manifold_dimension();
  Teuchos::RCP<const VerboseObject> verbobj = verbosity_obj();

  assert(setents != NULL);

  setents->clear();

  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm = geometric_model();

  // Is there an appropriate region by this name?

  Teuchos::RCP<const AmanziGeometry::Region> rgn;
  try {
    rgn = gm->FindRegion(setname);
  } catch (...) {
    valid_set_name(setname, kind);
  }

  // Did not find the region

  if (rgn == Teuchos::null) {
    std::stringstream mesg_stream;
    mesg_stream << "Geometric model has no region named \"" << setname << "\"\n";
    Errors::Message mesg(mesg_stream.str());
    Exceptions::amanzi_throw(mesg);
  }


  std::string internal_name = internal_name_of_set(rgn, kind);

  // If region is of type labeled set, a mesh set should have been
  // initialized from the input file. If region requires volume
  // fractions or is a segment, use base class capabilities to
  // build a mesh set. Otherwise, if a mesh set exists, search
  // the database for it. This is a two step procedure, which shows
  // probably defficiency of the internal naming convention (KL).
  // Finally, build a new mesh set for the region.

  if (rgn->get_type() == AmanziGeometry::RegionType::LABELEDSET && parent_mesh_.get() == NULL) {
    auto lsrgn = Teuchos::rcp_static_cast<const AmanziGeometry::RegionLabeledSet>(rgn);
    std::string label = lsrgn->label();
    std::string entity_type = lsrgn->entity_str();

    if ((kind == CELL && entity_type != "CELL") || (kind == FACE && entity_type != "FACE") ||
        (kind == NODE && entity_type != "NODE")) {
      if (verbobj.get() && verbobj->os_OK(Teuchos::VERB_MEDIUM)) {
        *(verbobj->os()) << "Found labeled set region \"" << setname
                         << "\" but it contains entities of type " << entity_type
                         << ", not the requested type\n";
      }
    } else {
      mset1 = MESH_MSetByName(mesh_, internal_name.c_str());

      if (!mset1 && kind == CELL) {
        // Since both element blocks and cell sets are referenced
        // with the region type 'Labeled Set' and Entity kind 'Cell'
        // we have to account for both possibilities. NOTE: THIS
        // MEANS THAT IF AN ELEMENT BLOCK AND ELEMENT SET HAVE THE
        // SAME ID, ONLY THE ELEMENT BLOCK WILL GET PICKED UP - WE
        // CHECKED FOR THIS IN BUILD SET

        std::string internal_name2 = other_internal_name_of_set(rgn, kind);
        mset1 = MESH_MSetByName(mesh_, internal_name2.c_str());
      }

      // Due to the parallel partitioning its possible that this
      // set is not on this processor

      if (!mset1) {
        if (comm_->NumProc() == 1) {
          Errors::Message msg;
          msg << "Could not find labeled set \"" << label
              << "\" in mesh file to initialize mesh set \"" << setname << "\". Verify mesh file.";
          amanzi_throw(msg);
        }
      }
    }
  } else if ((rgn->get_type() == AmanziGeometry::RegionType::BOX_VOF) ||
             (rgn->get_type() == AmanziGeometry::RegionType::LINE_SEGMENT)) {
    // Call routine from the base class and exit.
    Mesh::get_set_entities_box_vofs_(rgn, kind, ptype, setents, vofs);
    return;
  } else {
    // Modify region/set name by prefixing it with the type of
    // entity requested

    mset1 = MESH_MSetByName(mesh_, internal_name.c_str());

    // Make sure we retrieved a mesh set with the right kind of entities

    MType entdim;

    switch (kind) {
    case CELL:
      if (celldim == 3)
        entdim = MREGION;
      else if (celldim == 2)
        entdim = MFACE;
      break;

    case FACE:
      if (celldim == 3)
        entdim = MFACE;
      else if (celldim == 2)
        entdim = MEDGE;
      break;

    case NODE:
      entdim = MVERTEX;
      break;

    default:
      entdim = MUNKNOWNTYPE;
    }

    // If not, can we find a mesh set with the right name and right
    // kind of entities

    char setname1[256];

    if (mset1 && MSet_EntDim(mset1) != entdim) {
      idx = 0;
      while ((mset1 = MESH_Next_MSet(mesh_, &idx))) {
        MSet_Name(mset1, setname1);

        if (MSet_EntDim(mset1) == entdim && strcmp(setname1, internal_name.c_str()) == 0) break;
      }
    }
  }

  // All attempts to find the set failed so it must not exist - build it

  if (mset1 == NULL) { mset1 = build_set(rgn, kind); }

  // Check if no processor got any mesh entities

  int nent_loc = (mset1 == NULL) ? 0 : MSet_Num_Entries(mset1);

  setents->resize(nent_loc);
  Entity_ID_List::iterator it = setents->begin();

  if (nent_loc) {
    nent_loc = 0; // reset and count to get the real number

    switch (ptype) {
    case Parallel_type::OWNED:
      idx = 0;
      while ((ment = MSet_Next_Entry(mset1, &idx))) {
        if (MEnt_PType(ment) != PGHOST) {
          *it = MEnt_ID(ment) - 1; // assign to next spot by dereferencing iterator
          ++it;
          ++nent_loc;
        }
      }
      break;
    case Parallel_type::GHOST:
      idx = 0;
      while ((ment = MSet_Next_Entry(mset1, &idx))) {
        if (MEnt_PType(ment) == PGHOST) {
          *it = MEnt_ID(ment) - 1; // assign to next spot by dereferencing iterator
          ++it;
          ++nent_loc;
        }
      }
      break;
    case Parallel_type::ALL:
      idx = 0;
      while ((ment = MSet_Next_Entry(mset1, &idx))) {
        *it = MEnt_ID(ment) - 1; // assign to next spot by dereferencing iterator
        ++it;
        ++nent_loc;
      }
      break;
    default: {
    }
    }

    setents->resize(nent_loc);
  }
}


//---------------------------------------------------------
// Parent entity in the source mesh if mesh was derived from another mesh
//---------------------------------------------------------
Entity_ID
Mesh_MSTK::getEntityParent(const Entity_kind kind, const Entity_ID entid) const
{
  int ival;
  double rval;
  void* pval = nullptr;
  MEntity_ptr ment = nullptr;
  MAttrib_ptr att = nullptr;

  switch (kind) {
  case Entity_kind::CELL:
    att = (getManifoldDimension() == 3) ? rparentatt_ : fparentatt_;
    ment = (MEntity_ptr)cell_id_to_handle_[entid];
    break;
  case Entity_kind::FACE:
    att = (getManifoldDimension() == 3) ? fparentatt_ : eparentatt_;
    ment = (MEntity_ptr)face_id_to_handle_[entid];
    break;
  case Entity_kind::EDGE:
    att = eparentatt_;
    ment = (MEntity_ptr)edge_id_to_handle_[entid];
    break;
  case Entity_kind::NODE:
    if (!vparentatt_) return 0;
    att = vparentatt_;
    ment = (MEntity_ptr)vtx_id_to_handle_[entid];
    break;
  default: {
  }
  }
  if (!att) return 0;
  MEnt_Get_AttVal(ment, att, &ival, &rval, &pval);
  if (pval)
    return MEnt_ID((MEntity_ptr)pval) - 1;
  else
    return 0;
}


//---------------------------------------------------------
// Global ID of any entity
//---------------------------------------------------------
Entity_GID
Mesh_MSTK::getEntityGID(const Entity_kind kind, const Entity_ID lid) const
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

  if (serial_run)
    return MEnt_ID(ent) - 1;
  else
    return MEnt_GlobalID(ent) - 1;
}


//---------------------------------------------------------
// Procedure to perform all the post-mesh creation steps in a constructor
//---------------------------------------------------------
void
Mesh_MSTK::post_create_steps_()
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
void
Mesh_MSTK::clear_internals_()
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
void
Mesh_MSTK::init_nodes_()
{
  // create owned and not owned vertex lists
  init_pvert_lists_();
  // create maps from IDs to handles
  init_vertex_id2handle_maps_();
}


//---------------------------------------------------------
// Initialize edge info
//---------------------------------------------------------
void
Mesh_MSTK::init_edges_()
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
void
Mesh_MSTK::init_faces_()
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
void
Mesh_MSTK::init_cells_()
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
void
Mesh_MSTK::init_vertex_id2handle_maps_()
{
  // If the mesh is dynamic, then this code has to be revisited
  // Amanzi has IDs starting from 0, MSTK has IDs starting from 1
  int nv = MESH_Num_Vertices(mesh_);
  vtx_id_to_handle_.resize(nv);

  int idx = 0;
  int lid = 1;
  MVertex_ptr vtx;
  while ((vtx = MSet_Next_Entry(owned_verts_, &idx))) {
    MEnt_Set_ID(vtx, lid);
    vtx_id_to_handle_[lid - 1] = vtx;
    lid++;
  }

  idx = 0;
  while ((vtx = MSet_Next_Entry(ghost_verts_, &idx))) {
    MEnt_Set_ID(vtx, lid);
    vtx_id_to_handle_[lid - 1] = vtx;
    lid++;
  }
}


//---------------------------------------------------------
// ID to handle/pointer map for edges
//---------------------------------------------------------
void
Mesh_MSTK::init_edge_id2handle_maps_()
{
  // If the mesh is dynamic, then this code has to be revisited
  // Amanzi has IDs starting from 0, MSTK has IDs starting from 1
  int ne = MESH_Num_Edges(mesh_);
  edge_id_to_handle_.resize(ne);

  int idx = 0;
  int lid = 1;
  MEdge_ptr edge;
  while ((edge = MSet_Next_Entry(owned_edges_, &idx))) {
    MEnt_Set_ID(edge, lid);
    edge_id_to_handle_[lid - 1] = edge;
    lid++;
  }

  idx = 0;
  while ((edge = MSet_Next_Entry(ghost_edges_, &idx))) {
    MEnt_Set_ID(edge, lid);
    edge_id_to_handle_[lid - 1] = edge;
    lid++;
  }
}


//---------------------------------------------------------
// ID to handle/pointer map for faces
//---------------------------------------------------------
void
Mesh_MSTK::init_face_id2handle_maps_()
{
  // If the mesh is dynamic, then this code has to be revisited
  // Amanzi has IDs starting from 0, MSTK has IDs starting from 1
  int nf = (getManifoldDimension() == 2) ? MESH_Num_Edges(mesh_) : MESH_Num_Faces(mesh_);
  face_id_to_handle_.resize(nf);

  int idx = 0;
  int lid = 1;
  MEntity_ptr genface;
  while ((genface = MSet_Next_Entry(owned_faces_, &idx))) {
    MEnt_Set_ID(genface, lid);
    face_id_to_handle_[lid - 1] = genface;
    lid++;
  }

  idx = 0;
  while ((genface = MSet_Next_Entry(ghost_faces_, &idx))) {
    MEnt_Set_ID(genface, lid);
    face_id_to_handle_[lid - 1] = genface;
    lid++;
  }
}


//---------------------------------------------------------
// ID to handle/pointer map for cells
//---------------------------------------------------------
void
Mesh_MSTK::init_cell_id2handle_maps_()
{
  // If the mesh is dynamic, then this code has to be revisited
  // Amanzi has IDs starting from 0, MSTK has IDs starting from 1
  int nc = (getManifoldDimension() == 2) ? MESH_Num_Faces(mesh_) : MESH_Num_Regions(mesh_);
  cell_id_to_handle_.resize(nc);

  int idx = 0;
  int lid = 1;
  MEntity_ptr gencell; // Mesh region in 3D, face in 2D
  while ((gencell = MSet_Next_Entry(owned_cells_, &idx))) {
    MEnt_Set_ID(gencell, lid);
    cell_id_to_handle_[lid - 1] = gencell;
    lid++;
  }

  idx = 0;
  while ((gencell = MSet_Next_Entry(ghost_cells_, &idx))) {
    MEnt_Set_ID(gencell, lid);
    cell_id_to_handle_[lid - 1] = gencell;
    lid++;
  }
}


//---------------------------------------------------------
// create lists of owned and not owned vertices
//---------------------------------------------------------
void
Mesh_MSTK::init_pvert_lists_()
{
  // Get all vertices on this processor
  ghost_verts_ = MSet_New(mesh_, "ghost_verts_", MVERTEX);
  owned_verts_ = MSet_New(mesh_, "owned_verts_", MVERTEX);

  int idx = 0;
  MVertex_ptr vtx;
  while ((vtx = MESH_Next_Vertex(mesh_, &idx))) {
    if (MV_PType(vtx) == PGHOST)
      MSet_Add(ghost_verts_, vtx);
    else
      MSet_Add(owned_verts_, vtx);
  }
}


//---------------------------------------------------------
// create lists of owned and not owned edges
//---------------------------------------------------------
void
Mesh_MSTK::init_pedge_lists_()
{
  // Get all vertices on this processor
  ghost_edges_ = MSet_New(mesh_, "ghost_edges_", MEDGE);
  owned_edges_ = MSet_New(mesh_, "owned_edges_", MEDGE);

  int idx = 0;
  MEdge_ptr edge;
  while ((edge = MESH_Next_Edge(mesh_, &idx))) {
    if (ME_PType(edge) == PGHOST)
      MSet_Add(ghost_edges_, edge);
    else
      MSet_Add(owned_edges_, edge);
  }
} // Mesh_MSTK::init_pedge_lists_


void
Mesh_MSTK::init_pedge_dirs_()
{
  int ne = MESH_Num_Edges(mesh_);

  if (serial_run) {
    edgeflip_ = new bool[ne];
    for (int i = 0; i < ne; ++i) edgeflip_[i] = false;

  } else {
    // Do some additional processing to see if ghost edges and their masters
    // are oriented the same way; if not, turn on flag to flip the directions
    // when returning to the application code
    MAttrib_ptr attev0 = MAttrib_New(mesh_, "TMP_EV0_ATT", INT, MEDGE);
    MAttrib_ptr attev1 = MAttrib_New(mesh_, "TMP_EV1_ATT", INT, MEDGE);

    int idx = 0;
    MEdge_ptr edge;
    while ((edge = MESH_Next_Edge(mesh_, &idx))) {
      if (ME_PType(edge) != PINTERIOR) {
        MVertex_ptr vertex0 = ME_Vertex(edge, 0);
        MVertex_ptr vertex1 = ME_Vertex(edge, 1);

        MEnt_Set_AttVal(edge, attev0, MEnt_GlobalID(vertex0), 0.0, NULL);
        MEnt_Set_AttVal(edge, attev1, MEnt_GlobalID(vertex1), 0.0, NULL);
      }
    }

    MESH_UpdateAttributes(mesh_, mpicomm_);

    edgeflip_ = new bool[ne];
    for (int i = 0; i < ne; ++i) edgeflip_[i] = false;

    double rval;
    void* pval;
    idx = 0;
    while ((edge = MSet_Next_Entry(ghost_edges_, &idx))) {
      int remote_vertexid0, remote_vertexid1;

      MEnt_Get_AttVal(edge, attev0, &remote_vertexid0, &rval, &pval);
      MEnt_Get_AttVal(edge, attev1, &remote_vertexid1, &rval, &pval);

      int local_vertexid0 = MEnt_GlobalID(ME_Vertex(edge, 0));
      int local_vertexid1 = MEnt_GlobalID(ME_Vertex(edge, 1));

      if (remote_vertexid1 == local_vertexid0 || remote_vertexid0 == local_vertexid1) {
        int lid = MEnt_ID(edge);
        edgeflip_[lid - 1] = true;

      } else { // Sanity Check
        if (remote_vertexid1 != local_vertexid1 && remote_vertexid0 != local_vertexid0) {
          std::stringstream mesg_stream;
          mesg_stream << "Edge vertices mismatch between master and ghost (processor " << myprocid
                      << ")";
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
void
Mesh_MSTK::init_pface_lists_()
{
  // Get all faces on this processor
  if (getManifoldDimension() == 3) {
    ghost_faces_ = MSet_New(mesh_, "ghost_faces_", MFACE);
    owned_faces_ = MSet_New(mesh_, "owned_faces_", MFACE);

    int idx = 0;
    MFace_ptr face;
    while ((face = MESH_Next_Face(mesh_, &idx))) {
      if (MF_PType(face) == PGHOST)
        MSet_Add(ghost_faces_, face);
      else
        MSet_Add(owned_faces_, face);
    }

  } else if (getManifoldDimension() == 2) {
    ghost_faces_ = MSet_New(mesh_, "ghost_faces_", MFACE);
    owned_faces_ = MSet_New(mesh_, "owned_faces_", MFACE);

    int idx = 0;
    MEdge_ptr edge;
    while ((edge = MESH_Next_Edge(mesh_, &idx))) {
      if (ME_PType(edge) == PGHOST)
        MSet_Add(ghost_faces_, edge);
      else
        MSet_Add(owned_faces_, edge);
    }
  }
  return;
}

// Detect whether ghost faces are in opposite direction of owned faces
// on processor boundaries
void
Mesh_MSTK::init_pface_dirs_()
{
  int nf = (getManifoldDimension() == 2) ? MESH_Num_Edges(mesh_) : MESH_Num_Faces(mesh_);

  if (serial_run) {
    faceflip_ = new bool[nf];
    for (int i = 0; i < nf; ++i) faceflip_[i] = false;

  } else {
    if (getManifoldDimension() == 3)
      init_pface_dirs_3_();
    else if (getManifoldDimension() == 2)
      init_pface_dirs_2_();
  }
}


// Detect whether ghost faces are in opposite direction of owned faces
// on processor boundaries - Version for solid meshes
void
Mesh_MSTK::init_pface_dirs_3_()
{
  int nf = MESH_Num_Faces(mesh_);

  // Do some additional processing to see if ghost faces and their masters
  // are oriented the same way; if not, turn on flag to flip the directions
  // when returning to the application code

  // attributes to store
  MAttrib_ptr attfc0 = MAttrib_New(mesh_, "TMP_FC0_ATT", INT, MFACE);
  MAttrib_ptr attfc1 = MAttrib_New(mesh_, "TMP_FC1_ATT", INT, MFACE);

  int idx = 0;
  MFace_ptr face;
  while ((face = MESH_Next_Face(mesh_, &idx))) {
    if (MF_PType(face) != PINTERIOR) {
      MRegion_ptr region0 = MF_Region(face, 0);
      if (region0) MEnt_Set_AttVal(face, attfc0, MEnt_GlobalID(region0), 0.0, NULL);

      MRegion_ptr region1 = MF_Region(face, 1);
      if (region1) MEnt_Set_AttVal(face, attfc1, MEnt_GlobalID(region1), 0.0, NULL);
    }
  }

  MESH_UpdateAttributes(mesh_, mpicomm_);

  faceflip_ = new bool[nf];
  for (int i = 0; i < nf; ++i) faceflip_[i] = false;

  idx = 0;
  while ((face = MSet_Next_Entry(ghost_faces_, &idx))) {
    int remote_regid0, remote_regid1;
    double rval;
    void* pval;
    MEnt_Get_AttVal(face, attfc0, &remote_regid0, &rval, &pval);
    MEnt_Get_AttVal(face, attfc1, &remote_regid1, &rval, &pval);

    MRegion_ptr region0 = MF_Region(face, 0);
    int local_regid0 = region0 ? MEnt_GlobalID(region0) : 0;
    MRegion_ptr region1 = MF_Region(face, 1);
    int local_regid1 = region1 ? MEnt_GlobalID(region1) : 0;

    if (remote_regid1 == local_regid0 || remote_regid0 == local_regid1) {
      int lid = MEnt_ID(face);
      faceflip_[lid - 1] = true;

    } else { // Sanity Check
      if (remote_regid1 != local_regid1 && remote_regid0 != local_regid0) {
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
void
Mesh_MSTK::init_pface_dirs_2_()
{
  int ne = MESH_Num_Edges(mesh_);

  // Do some additional processing to see if ghost faces and their masters
  // are oriented the same way; if not, turn on flag to flip the directions
  // when returning to the application code
  MAttrib_ptr attev0 = MAttrib_New(mesh_, "TMP_EV0_ATT", INT, MEDGE);
  MAttrib_ptr attev1 = MAttrib_New(mesh_, "TMP_EV1_ATT", INT, MEDGE);

  int idx = 0;
  MEdge_ptr edge;
  while ((edge = MESH_Next_Edge(mesh_, &idx))) {
    if (ME_PType(edge) != PINTERIOR) {
      MVertex_ptr ev0 = ME_Vertex(edge, 0);
      MVertex_ptr ev1 = ME_Vertex(edge, 1);

      MEnt_Set_AttVal(edge, attev0, MEnt_GlobalID(ev0), 0.0, NULL);
      MEnt_Set_AttVal(edge, attev1, MEnt_GlobalID(ev1), 0.0, NULL);
    }
  }

  MESH_UpdateAttributes(mesh_, mpicomm_);

  faceflip_ = new bool[ne];
  for (int i = 0; i < ne; ++i) faceflip_[i] = false;

  idx = 0;
  while ((edge = MSet_Next_Entry(ghost_faces_, &idx))) {
    int remote_evgid0, remote_evgid1;
    double rval;
    void* pval;
    MEnt_Get_AttVal(edge, attev0, &remote_evgid0, &rval, &pval);
    MEnt_Get_AttVal(edge, attev1, &remote_evgid1, &rval, &pval);

    MVertex_ptr ev0 = ME_Vertex(edge, 0);
    MVertex_ptr ev1 = ME_Vertex(edge, 1);
    int local_evgid0 = MV_GlobalID(ev0);
    int local_evgid1 = MV_GlobalID(ev1);

    if (remote_evgid1 == local_evgid0 || remote_evgid0 == local_evgid1) {
      int lid = MEnt_ID(edge);
      faceflip_[lid - 1] = true;
    }
  }
  MAttrib_Delete(attev0);
  MAttrib_Delete(attev1);
}


//---------------------------------------------------------
// create lists of owned and not owned cells
//---------------------------------------------------------
void
Mesh_MSTK::init_pcell_lists_()
{
  int idx = 0;
  if (getManifoldDimension() == 3) {
    MRegion_ptr region;
    owned_cells_ = MSet_New(mesh_, "owned_cells_", MREGION);
    ghost_cells_ = MSet_New(mesh_, "ghost_cells_", MREGION);

    idx = 0;
    while ((region = MESH_Next_Region(mesh_, &idx))) {
      if (MR_PType(region) == PGHOST)
        MSet_Add(ghost_cells_, region);
      else
        MSet_Add(owned_cells_, region);
    }

  } else if (getManifoldDimension() == 2) {
    MFace_ptr face;
    owned_cells_ = MSet_New(mesh_, "owned_cells_", MFACE);
    ghost_cells_ = MSet_New(mesh_, "ghost_cells_", MFACE);

    idx = 0;
    while ((face = MESH_Next_Face(mesh_, &idx))) {
      if (MF_PType(face) == PGHOST)
        MSet_Add(ghost_cells_, face);
      else
        MSet_Add(owned_cells_, face);
    }

  } else {
    Errors::Message mesg("Implemented only for 2D and 3D");
    Exceptions::amanzi_throw(mesg);
  }
  return;
}

void
Mesh_MSTK::label_celltype_()
{
  if (getManifoldDimension() == 2) {
    celltype_att_ = MAttrib_New(mesh_, "Cell_kind", INT, MFACE);
    int idx = 0;
    MFace_ptr face;
    while ((face = MESH_Next_Face(mesh_, &idx))) {
      Entity_ID_List vedges;
      List_ptr fedges = MF_Edges(face, 0, 0);
      int idx2 = 0;
      MEdge_ptr edge;
      while ((edge = List_Next_Entry(fedges, &idx2))) vedges.push_back(MEnt_ID(edge) - 1);
      Entity_ID_View edges;
      vectorToView(edges, vedges);
      Cell_kind ctype = MeshFramework::getCellType_(MEnt_ID(face), edges);
      MEnt_Set_AttVal(face, celltype_att_, (int)ctype, 0.0, NULL);
    }

  } else if (getManifoldDimension() == 3) {
    celltype_att_ = MAttrib_New(mesh_, "Cell_kind", INT, MREGION);
    int idx = 0;
    MRegion_ptr region;
    while ((region = MESH_Next_Region(mesh_, &idx))) {
      List_ptr rfaces = MR_Faces(region);
      Entity_ID_List vfaces;
      int idx2 = 0;
      MFace_ptr face;
      while ((face = List_Next_Entry(rfaces, &idx2))) vfaces.push_back(MEnt_ID(face) - 1);
      Entity_ID_View faces;
      vectorToView(faces, vfaces);

      Cell_kind ctype = MeshFramework::getCellType_(MEnt_ID(region), faces);
      MEnt_Set_AttVal(region, celltype_att_, (int)ctype, 0.0, NULL);
    }
  }
}


void
Mesh_MSTK::getSetEntities(const AmanziGeometry::RegionLabeledSet& region,
                          const Entity_kind kind,
                          const Parallel_kind ptype,
                          View_type<const Entity_ID, MemSpace_kind::HOST>& entids) const
{
  Entity_ID_View lentids;
  if (kind != createEntityKind(region.entity_str())) {
    Errors::Message msg;
    msg << "Inconsistent request of labeled set for region \"" << region.get_name()
        << "\" of kind \"" << region.entity_str() << "\" requested as type \"" << to_string(kind)
        << "\"";
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
        msg << "Exodus II file has both element block and element set with ID " << region.label()
            << " - Amanzi cannot handle this case.";
        Exceptions::amanzi_throw(msg);
      }

    } else {
      if (mset2) { mset = mset2; }
    }
  }

  if (mset) {
    // Its possible some sets won't exist on some partitions
    int entdim = MSet_EntDim(mset);
    if (getManifoldDimension() == 3) {
      if ((region.entity_str() == "CELL" && entdim != MREGION) ||
          (region.entity_str() == "FACE" && entdim != MFACE) ||
          (region.entity_str() == "NODE" && entdim != MVERTEX)) {
        Errors::Message mesg("Mismatch of entity type in labeled set region and mesh set");
        Exceptions::amanzi_throw(mesg);
      }
    } else if (getManifoldDimension() == 2) {
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
      while ((ent = MSet_Next_Entry(mset, &idx))) {
        if (MEnt_Dim(ent) == MDELETED) MSet_Rem(mset, ent);
      }
    }

    int nent_loc = MSet_Num_Entries(mset);
    Kokkos::resize(lentids, nent_loc);

    if (nent_loc) {
      nent_loc = 0; // reset and count to get the real number
      int idx = 0;
      MEntity_ptr ment;
      switch (ptype) {
      case Parallel_kind::OWNED:
        idx = 0;
        while ((ment = MSet_Next_Entry(mset, &idx))) {
          if (MEnt_PType(ment) != PGHOST) {
            lentids[nent_loc] = MEnt_ID(ment) - 1;
            ++nent_loc;
          }
        }
        break;
      case Parallel_kind::GHOST:
        idx = 0;
        while ((ment = MSet_Next_Entry(mset, &idx))) {
          if (MEnt_PType(ment) == PGHOST) {
            lentids[nent_loc] = MEnt_ID(ment) - 1;
            ++nent_loc;
          }
        }
        break;
      case Parallel_kind::ALL:
        idx = 0;
        while ((ment = MSet_Next_Entry(mset, &idx))) {
          lentids[nent_loc] = MEnt_ID(ment) - 1;
          ++nent_loc;
        }
        break;
      default: {
      }
      }
      Kokkos::resize(lentids, nent_loc);
    }
  }
  entids = lentids;
}


void
Mesh_MSTK::collapse_degen_edges_()
{
  const int topoflag = 0; // Don't worry about violation of model classification
  int idx, evgid0, evgid1;
  MVertex_ptr vertex, ev0, ev1, vkeep, vdel;
  MEdge_ptr edge;
  MFace_ptr face;
  List_ptr deleted_ents_all = List_New(10);
  List_ptr merged_entity_pairs_all = List_New(10);
  double len2;
  std::vector<int> merged_ents_info;

  idx = 0;
  while ((edge = MESH_Next_Edge(mesh_, &idx))) {
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
      ev0 = ME_Vertex(edge, 0);
      evgid0 = MEnt_GlobalID(ev0);
      ev1 = ME_Vertex(edge, 1);
      evgid1 = MEnt_GlobalID(ev1);

      if (evgid0 < evgid1) {
        vkeep = ev0;
        vdel = ev1;
      } else {
        vkeep = ev1;
        vdel = ev0;
      }

      List_ptr deleted_ents = NULL, merged_entity_pairs = NULL;
      vkeep = ME_Collapse(edge, vkeep, topoflag, &deleted_ents, &merged_entity_pairs);

      if (!vkeep) {
        vkeep = vdel;
        vdel = (vkeep == ev0) ? ev1 : ev0;

        vkeep = ME_Collapse(edge, vkeep, topoflag, &deleted_ents, &merged_entity_pairs);
      }

      if (!vkeep) {
        Errors::Message mesg("Could not collapse degenerate edge. Expect computational issues with "
                             "connected elements");
        Exceptions::amanzi_throw(mesg);
      } else {
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
  int nmerged = List_Num_Entries(merged_entity_pairs_all) / 2;
  for (int j = 0; j < nmerged; j++) {
    MEntity_ptr delent = List_Entry(merged_entity_pairs_all, 2 * j);
    MEntity_ptr keepent = List_Entry(merged_entity_pairs_all, 2 * j + 1);
    merged_ents_info.push_back(static_cast<int>(MEnt_Dim(keepent)));
    merged_ents_info.push_back(MEnt_GlobalID(delent));
    merged_ents_info.push_back(MEnt_GlobalID(keepent));
  }

  int* nmerged_proc = new int[numprocs];
  int* nmerged_proc_x3 = new int[numprocs];
  MPI_Allgather(&nmerged, 1, MPI_INT, nmerged_proc, 1, MPI_INT, mpicomm_);

  int* offset = new int[numprocs];
  int nmerged_global = 0;
  for (int p = 0; p < numprocs; p++) {
    offset[p] = 3 * nmerged_global;
    nmerged_global += nmerged_proc[p];
    nmerged_proc_x3[p] = 3 * nmerged_proc[p];
  }

  // We probably can make this more efficient by using point-to-point
  // communication

  int* merged_ents_info_global = new int[3 * nmerged_global];
  MPI_Allgatherv(&(merged_ents_info[0]),
                 3 * nmerged,
                 MPI_INT,
                 merged_ents_info_global,
                 nmerged_proc_x3,
                 offset,
                 MPI_INT,
                 mpicomm_);

  idx = 0;
  while ((vertex = MESH_Next_Vertex(mesh_, &idx))) {
    if (MV_PType(vertex) == PGHOST) {
      int vgid = MV_GlobalID(vertex);
      for (int i = 0; i < nmerged_global; i++) {
        if (merged_ents_info_global[3 * i] == MVERTEX &&
            merged_ents_info_global[3 * i + 1] == vgid) {
          // Found vertex that got deleted and replaced by another vtx
          // on a different proc
          MV_Set_GlobalID(vertex, merged_ents_info_global[3 * i + 2]);
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
        if (merged_ents_info_global[3 * i] == MEDGE && merged_ents_info_global[3 * i + 1] == egid) {
          // Found edge that got deleted and replaced by another edge
          // on a different proc
          ME_Set_GlobalID(edge, merged_ents_info_global[3 * i + 2]);
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
        if (merged_ents_info_global[3 * i] == MFACE && merged_ents_info_global[3 * i + 1] == fgid) {
          // Found face that got deleted and replaced by another face
          // on a different proc
          MF_Set_GlobalID(face, merged_ents_info_global[3 * i + 2]);
          break;
        }
      }
    }
  }

  delete[] nmerged_proc;
  delete[] nmerged_proc_x3;
  delete[] merged_ents_info_global;
  delete[] offset;

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
      if (iloc != -1) // found deleted entity in set; replace it with keepent
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
      if (iloc != -1) // found deleted entity in set; replace it with keepent
        MSet_Remi(mset, iloc);
    }
  }

  // ME_Collapse only marked these entities as DELETED but now
  // delete them for good
  idx = 0;
  while ((delent = List_Next_Entry(deleted_ents_all, &idx))) MEnt_Delete(delent, 0);

  List_Delete(deleted_ents_all);
  List_Delete(merged_entity_pairs_all);

  // Now renumber global IDs to make them contiguous
  //  if (entities_deleted) {
  //   std::cerr << "Entities deleted in collapse_degen_edges ..." << "\n";
  //   MESH_Renumber_GlobalIDs(mesh, MALLTYPE, 0, NULL, mpicomm_);
  //  }

#ifdef DEBUG
  if (MESH_Num_Regions(mesh_) > 0) { // 3D mesh
    idx = 0;
    while ((face = MESH_Next_Face(mesh_, &idx))) {
      List_ptr fregs = MF_Regions(face);
      if (fregs)
        List_Delete(fregs);
      else {
        std::cerr << "Dangling mesh face with no connected cells AFTER COLLAPSE\n on P" << myprocid
                  << "\n";
        MF_Print(face, 3);
      }
    }
  }
#endif
}


//
// Construction function
//
void
Mesh_MSTK::init_mesh_from_file_(const std::string& filename)
{
  int ok = 0;
  mesh_ = MESH_New(F1);

  if (filename.find(".exo") != std::string::npos) { // Exodus file
    // Read the mesh on processor 0
    ok = MESH_ImportFromExodusII(mesh_, filename.c_str(), NULL, mpicomm_);

    // Collapse any degenerate edges in the mesh
    collapse_degen_edges_(); // Assumes its operating on member var 'mesh'

    // Renumber local IDs to be contiguous
    MESH_Renumber(mesh_, 0, MALLTYPE);

    if (numprocs > 1) {
      // Distribute the mesh to all the processors
      int topo_dim = MESH_Num_Regions(mesh_) ? 3 : 2;
      int num_ghost_layers = 1;
      int with_attr = 1; // Redistribute any attributes and sets
      int method = static_cast<int>(partitioner_);
      int del_inmesh = 1; // Delete input mesh (on P0) after distribution

      Mesh_ptr globalmesh = mesh_;
      mesh_ = MESH_New(F1);

      ok &= MSTK_Mesh_Distribute(
        globalmesh, &mesh_, &topo_dim, num_ghost_layers, with_attr, method, del_inmesh, mpicomm_);
      if (contiguous_gids_) { ok &= MESH_Renumber_GlobalIDs(mesh_, MALLTYPE, 0, NULL, mpicomm_); }
    }

  } else if (filename.find(".par") != std::string::npos) { // Nemesis file
    // Read the individual partitions on each processor
    ok = MESH_ImportFromNemesisI(mesh_, filename.c_str(), NULL, mpicomm_);

    // Collapse any degenerate edges in the mesh
    collapse_degen_edges_();

    // Renumber local IDs to be contiguous
    MESH_Renumber(mesh_, 0, MALLTYPE);

    // Weave the meshes together to form interprocessor connections
    int num_ghost_layers = 1;
    int input_type = 1; // We are given partitioned meshes with a
                        // unique global ID on each mesh vertex
    int topo_dim = MESH_Num_Regions(mesh_) ? 3 : 2;
    ok &= MSTK_Weave_DistributedMeshes(mesh_, topo_dim, num_ghost_layers, input_type, mpicomm_);

    if (contiguous_gids_) { ok &= MESH_Renumber_GlobalIDs(mesh_, MALLTYPE, 0, NULL, mpicomm_); }

  } else {
    Errors::Message msg;
    msg << "Cannot identify file type from extension of input file " << filename << " on processor "
        << myprocid;
    Exceptions::amanzi_throw(msg);
  }

  if (!ok) {
    Errors::Message msg;
    msg << "Failed to load " << filename << " on processor " << myprocid;
    Exceptions::amanzi_throw(msg);
  }
}


int
Mesh_MSTK::generate_regular_mesh_(Mesh_ptr mesh,
                                  double x0,
                                  double y0,
                                  double z0,
                                  double x1,
                                  double y1,
                                  double z1,
                                  int nx,
                                  int ny,
                                  int nz)
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
  int vgid_tmpl[3][3][3] = { { { 1, 4, 5 }, { 9, 6, 12 }, { 3, 8, 7 } },
                             { { 1, 1, 3 }, { 3, 1, 4 }, { 5, 2, 7 } },
                             { { 2, 2, 6 }, { 10, 5, 11 }, { 4, 6, 8 } } };
  int vgdim_tmpl[3][3][3] = { { { 0, 1, 0 }, { 1, 2, 1 }, { 0, 1, 0 } },
                              { { 1, 2, 1 }, { 2, 3, 2 }, { 1, 2, 1 } },
                              { { 0, 1, 0 }, { 1, 2, 1 }, { 0, 1, 0 } } };
  int egdim_tmpl[3][3] = { { 1, 2, 1 }, { 2, 3, 2 }, { 1, 2, 1 } };
  int egid_tmpl2[3][3] = { { 4, 6, 8 },
                           { 1, 1, 2 },
                           { 2, 5, 6 } }; /* Y direction edges (iterating over i,k) */
  int egid_tmpl1[3][3] = { { 9, 6, 12 },
                           { 3, 1, 4 },
                           { 10, 5, 11 } }; /* Z direction edges (iterating over i,j)*/
  int egid_tmpl0[3][3] = { { 1, 1, 3 },
                           { 3, 1, 4 },
                           { 5, 2, 7 } }; /* X direction edges (iterating over j,k) */
  int fgdim_tmpl[3] = { 2, 3, 2 };
  int fgid_tmpl0[3] = { 6, 1, 5 };
  int fgid_tmpl1[3] = { 1, 1, 2 };
  int fgid_tmpl2[3] = { 3, 1, 4 };

  dx = (x1 - x0) / nx;
  dy = (y1 - y0) / ny;
  dz = (z1 - z0) / nz;

  verts = (MVertex_ptr***)malloc((nx + 1) * sizeof(MVertex_ptr**));
  for (j = 0; j < nx + 1; ++j) {
    verts[j] = (MVertex_ptr**)malloc((ny + 1) * sizeof(MVertex_ptr*));
    for (k = 0; k < ny + 1; ++k) verts[j][k] = (MVertex_ptr*)malloc((nz + 1) * sizeof(MVertex_ptr));
  }

  for (k = 0; k < nz + 1; ++k) {
    xyz[2] = (k == nz) ? z1 : z0 + k * dz;
    kk = (k % nz) ? 1 : (k ? 2 : 0);

    for (j = 0; j < ny + 1; ++j) {
      xyz[1] = (j == ny) ? y1 : y0 + j * dy;
      jj = (j % ny) ? 1 : (j ? 2 : 0);

      for (i = 0; i < nx + 1; ++i) {
        xyz[0] = (i == nx) ? x1 : x0 + i * dx;
        ii = (i % nx) ? 1 : (i ? 2 : 0);

        mv = MV_New(mesh);
        MV_Set_Coords(mv, xyz);
        verts[i][j][k] = mv;

        gdim = vgdim_tmpl[ii][jj][kk];
        MV_Set_GEntDim(mv, gdim);

        gid = vgid_tmpl[ii][jj][kk];
        MV_Set_GEntID(mv, gid);
      }
    }
  }


  /* Create the edges explicitly to get the classification right */
  for (i = 0; i < nx + 1; ++i) {
    for (j = 0; j < ny + 1; ++j) {
      for (k = 0; k < nz; ++k) {
        me = ME_New(mesh);

        everts[0] = verts[i][j][k];
        everts[1] = verts[i][j][k + 1];
        ME_Set_Vertex(me, 0, everts[0]);
        ME_Set_Vertex(me, 1, everts[1]);

        ii = (i % nx) ? 1 : (i ? 2 : 0);
        jj = (j % ny) ? 1 : (j ? 2 : 0);
        gdim = egdim_tmpl[ii][jj];
        gid = egid_tmpl2[ii][jj];

        ME_Set_GEntDim(me, gdim);
        ME_Set_GEntID(me, gid);
      }
    }
  }

  for (i = 0; i < nx + 1; ++i) {
    for (k = 0; k < nz + 1; ++k) {
      for (j = 0; j < ny; ++j) {
        me = ME_New(mesh);

        everts[0] = verts[i][j][k];
        everts[1] = verts[i][j + 1][k];
        ME_Set_Vertex(me, 0, everts[0]);
        ME_Set_Vertex(me, 1, everts[1]);

        ii = (i % nx) ? 1 : (i ? 2 : 0);
        kk = (k % nz) ? 1 : (k ? 2 : 0);
        gdim = egdim_tmpl[ii][kk];
        gid = egid_tmpl1[ii][kk];

        ME_Set_GEntDim(me, gdim);
        ME_Set_GEntID(me, gid);
      }
    }
  }

  for (j = 0; j < ny + 1; ++j) {
    for (k = 0; k < nz + 1; ++k) {
      for (i = 0; i < nx; ++i) {
        me = ME_New(mesh);

        everts[0] = verts[i][j][k];
        everts[1] = verts[i + 1][j][k];
        ME_Set_Vertex(me, 0, everts[0]);
        ME_Set_Vertex(me, 1, everts[1]);

        jj = (j % ny) ? 1 : (j ? 2 : 0);
        kk = (k % nz) ? 1 : (k ? 2 : 0);
        gdim = egdim_tmpl[jj][kk];
        gid = egid_tmpl0[jj][kk];

        ME_Set_GEntDim(me, gdim);
        ME_Set_GEntID(me, gid);
      }
    }
  }


  /* Create the faces explicitly to get the classification right */
  for (i = 0; i < nx + 1; ++i) {
    for (j = 0; j < ny; ++j) {
      for (k = 0; k < nz; ++k) {
        mf = MF_New(mesh);

        fverts[0] = verts[i][j][k];
        fverts[1] = verts[i][j + 1][k];
        fverts[2] = verts[i][j + 1][k + 1];
        fverts[3] = verts[i][j][k + 1];
        MF_Set_Vertices(mf, 4, fverts);

        ii = (i % nx) ? 1 : (i ? 2 : 0);
        gdim = fgdim_tmpl[ii];
        gid = fgid_tmpl0[ii];

        MF_Set_GEntDim(mf, gdim);
        MF_Set_GEntID(mf, gid);
      }
    }
  }

  for (j = 0; j < ny + 1; ++j) {
    for (i = 0; i < nx; ++i) {
      for (k = 0; k < nz; ++k) {
        mf = MF_New(mesh);

        fverts[0] = verts[i][j][k];
        fverts[1] = verts[i + 1][j][k];
        fverts[2] = verts[i + 1][j][k + 1];
        fverts[3] = verts[i][j][k + 1];
        MF_Set_Vertices(mf, 4, fverts);

        jj = (j % ny) ? 1 : (j ? 2 : 0);
        gdim = fgdim_tmpl[jj];
        gid = fgid_tmpl1[jj];

        MF_Set_GEntDim(mf, gdim);
        MF_Set_GEntID(mf, gid);
      }
    }
  }

  for (k = 0; k < nz + 1; ++k) {
    for (i = 0; i < nx; ++i) {
      for (j = 0; j < ny; ++j) {
        mf = MF_New(mesh);

        fverts[0] = verts[i][j][k];
        fverts[1] = verts[i + 1][j][k];
        fverts[2] = verts[i + 1][j + 1][k];
        fverts[3] = verts[i][j + 1][k];
        MF_Set_Vertices(mf, 4, fverts);

        kk = (k % nz) ? 1 : (k ? 2 : 0);
        gdim = fgdim_tmpl[kk];
        gid = fgid_tmpl2[kk];

        MF_Set_GEntDim(mf, gdim);
        MF_Set_GEntID(mf, gid);
      }
    }
  }


  /* Not the most efficient way but the easiest to code */

  for (i = 0; i < nx; ++i) {
    for (j = 0; j < ny; ++j) {
      for (k = 0; k < nz; ++k) {
        mr = MR_New(mesh);
        MR_Set_GEntID(mr, 1);

        rverts[0] = verts[i][j][k];
        rverts[1] = verts[i + 1][j][k];
        rverts[2] = verts[i + 1][j + 1][k];
        rverts[3] = verts[i][j + 1][k];
        rverts[4] = verts[i][j][k + 1];
        rverts[5] = verts[i + 1][j][k + 1];
        rverts[6] = verts[i + 1][j + 1][k + 1];
        rverts[7] = verts[i][j + 1][k + 1];

        MR_Set_Vertices(mr, 8, rverts, 6, NULL);
      }
    }
  }

  for (i = 0; i < nx + 1; ++i) {
    for (j = 0; j < ny + 1; ++j) free(verts[i][j]);
    free(verts[i]);
  }
  free(verts);

  return 1;
}


int
Mesh_MSTK::generate_regular_mesh_(Mesh_ptr mesh,
                                  double x0,
                                  double y0,
                                  double x1,
                                  double y1,
                                  int nx,
                                  int ny)
{
  int i, j, dir[4];
  double xyz[3], dx, dy;
  MVertex_ptr **verts, v0, v1, mv;
  MEdge_ptr fedges[4], me;
  MFace_ptr mf;

  dx = (x1 - x0) / nx;
  dy = (y1 - y0) / ny;

  verts = (MVertex_ptr**)malloc((nx + 1) * sizeof(MVertex_ptr*));
  for (i = 0; i < nx + 1; ++i) verts[i] = (MVertex_ptr*)malloc((ny + 1) * sizeof(MVertex_ptr));

  xyz[2] = 0.0;
  for (j = 0; j < ny + 1; ++j) {
    xyz[1] = (j == ny) ? y1 : y0 + j * dy;

    for (i = 0; i < nx + 1; ++i) {
      xyz[0] = (i == nx) ? x1 : x0 + i * dx;

      mv = MV_New(mesh);
      MV_Set_Coords(mv, xyz);

      if (i == 0) {
        if (j == 0) {
          MV_Set_GEntDim(mv, 0);
          MV_Set_GEntID(mv, 1);
        } else if (j == ny) {
          MV_Set_GEntDim(mv, 0);
          MV_Set_GEntID(mv, 4);
        } else {
          MV_Set_GEntDim(mv, 1);
          MV_Set_GEntID(mv, 4);
        }
      } else if (i == nx) {
        if (j == 0) {
          MV_Set_GEntDim(mv, 0);
          MV_Set_GEntID(mv, 2);
        } else if (j == ny) {
          MV_Set_GEntDim(mv, 0);
          MV_Set_GEntID(mv, 3);
        } else {
          MV_Set_GEntDim(mv, 1);
          MV_Set_GEntID(mv, 2);
        }
      } else {
        if (j == 0) {
          MV_Set_GEntDim(mv, 1);
          MV_Set_GEntID(mv, 1);
        } else if (j == ny) {
          MV_Set_GEntDim(mv, 1);
          MV_Set_GEntID(mv, 3);
        } else {
          MV_Set_GEntDim(mv, 2);
          MV_Set_GEntID(mv, 1);
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
      v1 = verts[i + 1][j];
      fedges[0] = MVs_CommonEdge(v0, v1);
      if (fedges[0])
        dir[0] = (ME_Vertex(fedges[0], 0) == v0) ? 1 : 0;
      else {
        me = ME_New(mesh);

        ME_Set_Vertex(me, 0, v0);
        ME_Set_Vertex(me, 1, v1);

        if (j == 0) {
          ME_Set_GEntDim(me, 1);
          ME_Set_GEntID(me, 1);
        } else {
          ME_Set_GEntDim(me, 2);
          ME_Set_GEntID(me, 1);
        }

        fedges[0] = me;
        dir[0] = 1;
      }


      /* edge 1 */
      v0 = verts[i + 1][j];
      v1 = verts[i + 1][j + 1];
      fedges[1] = MVs_CommonEdge(v0, v1);
      if (fedges[1])
        dir[1] = (ME_Vertex(fedges[1], 0) == v0) ? 1 : 0;
      else {
        me = ME_New(mesh);

        ME_Set_Vertex(me, 0, v0);
        ME_Set_Vertex(me, 1, v1);

        if (i + 1 == nx) {
          ME_Set_GEntDim(me, 1);
          ME_Set_GEntID(me, 2);
        } else {
          ME_Set_GEntDim(me, 2);
          ME_Set_GEntID(me, 1);
        }

        fedges[1] = me;
        dir[1] = 1;
      }


      /* edge 2 */
      v0 = verts[i + 1][j + 1];
      v1 = verts[i][j + 1];
      fedges[2] = MVs_CommonEdge(v0, v1);
      if (fedges[2])
        dir[2] = (ME_Vertex(fedges[2], 0) == v0) ? 1 : 0;
      else {
        me = ME_New(mesh);

        ME_Set_Vertex(me, 0, v0);
        ME_Set_Vertex(me, 1, v1);

        if (j + 1 == nx) {
          ME_Set_GEntDim(me, 1);
          ME_Set_GEntID(me, 3);
        } else {
          ME_Set_GEntDim(me, 2);
          ME_Set_GEntID(me, 1);
        }

        fedges[2] = me;
        dir[2] = 1;
      }


      /* edge 3 */
      v0 = verts[i][j + 1];
      v1 = verts[i][j];
      fedges[3] = MVs_CommonEdge(v0, v1);
      if (fedges[3])
        dir[3] = (ME_Vertex(fedges[3], 0) == v0) ? 1 : 0;
      else {
        me = ME_New(mesh);

        ME_Set_Vertex(me, 0, v0);
        ME_Set_Vertex(me, 1, v1);

        if (i == 0) {
          ME_Set_GEntDim(me, 1);
          ME_Set_GEntID(me, 4);
        } else {
          ME_Set_GEntDim(me, 2);
          ME_Set_GEntID(me, 1);
        }

        fedges[3] = me;
        dir[3] = 1;
      }


      MF_Set_Edges(mf, 4, fedges, dir);

      MF_Set_GEntDim(mf, 2);
      MF_Set_GEntID(mf, 1);
    }
  }

  for (i = 0; i < nx + 1; ++i) free(verts[i]);
  free(verts);

  return 1;
}


void
Mesh_MSTK::pre_create_steps_(const int space_dimension)
{
  clear_internals_();
  MSTK_Init();
  setSpaceDimension(space_dimension);

  auto mpicomm = Teuchos::rcp_dynamic_cast<const MpiComm_type>(getComm());
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


void
Mesh_MSTK::inherit_labeled_sets_(MAttrib_ptr copyatt, List_ptr src_entities)
{
  int idx, idx2, diffdim;
  MSet_ptr mset;

  auto gm = getGeometricModel();
  if (gm == Teuchos::null) return;

  Mesh_ptr parent_mstk_mesh = parent_mesh_->mesh_;

  // Difference in cell dimension of this mesh and its parent
  // Labeled set entity dimensions will be similarly dialed down
  diffdim = parent_mesh_->getManifoldDimension() - getManifoldDimension();
  if (diffdim > 1) {
    Errors::Message mesg("Dimension of mesh and its parent differ by more than 1");
    Exceptions::amanzi_throw(mesg);
  }

  unsigned int ngr = gm->size();
  for (int i = 0; i < ngr; ++i) {
    auto rgn = gm->FindRegion(i);
    if (rgn->get_type() == AmanziGeometry::RegionType::LABELEDSET) {
      // Get the set from the parent mesh
      Teuchos::RCP<const AmanziGeometry::RegionLabeledSet> lsrgn =
        Teuchos::rcp_static_cast<const AmanziGeometry::RegionLabeledSet>(rgn);

      std::string internal_name;
      std::string label = lsrgn->label();

      if (lsrgn->entity_str() == "CELL")
        internal_name = internal_name_of_set_(*lsrgn, Entity_kind::CELL);
      else if (lsrgn->entity_str() == "FACE")
        internal_name = internal_name_of_set_(*lsrgn, Entity_kind::FACE);
      else if (lsrgn->entity_str() == "NODE")
        internal_name = internal_name_of_set_(*lsrgn, Entity_kind::NODE);

      MSet_ptr msetParentMesh = MESH_MSetByName(parent_mstk_mesh, internal_name.c_str());
      if (!msetParentMesh) continue;

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
          if (MSet_Contains(msetParentMesh, ent)) {
            found = 1;
            break;
          }
        }
        if (found) continue;
      }

      // Create the set in this mesh
      MType subentdim;
      MType entdim = MSet_EntDim(msetParentMesh);
      if (entdim == MVERTEX)
        subentdim = MVERTEX;
      else
        subentdim = (MType)(entdim - diffdim);
      mset = MSet_New(mesh_, internal_name.c_str(), subentdim);

      // Populate the set
      int mkid = MSTK_GetMarker();

      MEntity_ptr ent;
      idx = 0;
      while ((ent = MSet_Next_Entry(msetParentMesh, &idx))) {
        if (MEnt_Dim(ent) == MDELETED) continue;
        MEntity_ptr copyent;
        int ival;
        double rval;
        if (subentdim == entdim) {
          MEnt_Get_AttVal(ent, copyatt, &ival, &rval, &copyent);
          if (!copyent) continue;
          MSet_Add(mset, copyent);

        } else {
          if (entdim == MREGION) {
            MFace_ptr rf;
            List_ptr rfaces = MR_Faces((MRegion_ptr)ent);
            idx2 = 0;
            while ((rf = List_Next_Entry(rfaces, &idx2))) {
              MEnt_Get_AttVal(rf, copyatt, &ival, &rval, &copyent);
              if (!copyent) continue;
              if (!MEnt_IsMarked(copyent, mkid)) {
                MSet_Add(mset, copyent);
                MEnt_Mark(copyent, mkid);
              }
            }
            List_Delete(rfaces);

          } else if (entdim == MFACE) {
            MEdge_ptr fe;
            List_ptr fedges = MF_Edges((MFace_ptr)ent, 1, 0);
            idx2 = 0;
            while ((fe = List_Next_Entry(fedges, &idx2))) {
              MEnt_Get_AttVal(fe, copyatt, &ival, &rval, &copyent);
              if (!copyent) continue;

              if (!MEnt_IsMarked(copyent, mkid)) {
                MSet_Add(mset, copyent);
                MEnt_Mark(copyent, mkid);
              }
            }
            List_Delete(fedges);
          }
        }
      }

      MSet_Unmark(mset, mkid);
      MSTK_FreeMarker(mkid);
    }
  }
}


//---------------------------------------------------------
// Extract a list of MSTK entities and make a new MSTK mesh
// For private use of Mesh_MSTK class only
//---------------------------------------------------------
void
Mesh_MSTK::extract_mstk_mesh_(List_ptr src_entities,
                              const MType entity_dim,
                              const bool flatten,
                              const bool request_faces,
                              const bool request_edges)
{
  int ival = 0, idx;
  double rval = 0., xyz[3];
  void* pval;

  AMANZI_ASSERT(parent_mesh_.get());
  Mesh_ptr parent_mesh_mstk = parent_mesh_->mesh_;

  // Make sure Global ID searches are enabled
  MESH_Enable_GlobalIDSearch(parent_mesh_mstk);

  if (flatten) {
    if (entity_dim == MREGION || entity_dim == MVERTEX) {
      Errors::Message mesg("Flattening or extruding allowed only for sets of FACEs in volume mesh "
                           "or CELLs in surface meshes");
      Exceptions::amanzi_throw(mesg);
    }
  }

  if (entity_dim == MEDGE) {
    Errors::Message mesg(
      "Requested mesh constructor produces 1D mesh which is not supported by Amanzi");
    Exceptions::amanzi_throw(mesg);
  }

  // Pre-processing (init, MPI queries etc)
  if (flatten)
    pre_create_steps_(getParentMesh()->getSpaceDimension() - 1);
  else
    pre_create_steps_(getParentMesh()->getSpaceDimension());

  // What is the cell dimension of new mesh
  switch (entity_dim) {
  case MREGION:
    setManifoldDimension(3); // extract regions/cells from mesh
    break;
  case MFACE:
    setManifoldDimension(2); // extract faces from mesh
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
  MAttrib_ptr copyatt = MAttrib_New(parent_mesh_mstk, "copyatt", POINTER, MALLTYPE);
  vparentatt_ = MAttrib_New(mesh_, "vparentatt_", POINTER, MVERTEX);
  eparentatt_ = MAttrib_New(mesh_, "eparentatt_", POINTER, MEDGE);
  fparentatt_ = MAttrib_New(mesh_, "fparentatt_", POINTER, MFACE);
  rparentatt_ = MAttrib_New(mesh_, "rparentatt_", POINTER, MREGION);

  switch (entity_dim) {
  case MREGION: { // Extracting a subset of a solid mesh
    idx = 0;
    MRegion_ptr mr;
    while ((mr = (MRegion_ptr)List_Next_Entry(src_entities, &idx))) {
      List_ptr rfaces = MR_Faces(mr);
      int nrf = List_Num_Entries(rfaces);
      MFace_ptr rfaces_new[MAXPF3];
      int rfdirs_new[MAXPF3];
      for (int i = 0; i < nrf; ++i) {
        MFace_ptr mf = List_Entry(rfaces, i);

        MEnt_Get_AttVal(mf, copyatt, &ival, &rval, &pval);
        if (pval) {
          rfaces_new[i] = pval;
          rfdirs_new[i] = MR_FaceDir_i(mr, i);
        } else {
          List_ptr fverts = MF_Vertices(mf, 1, 0);
          int nfv = List_Num_Entries(fverts);
          MVertex_ptr fverts_new[MAXPV2];
          for (int j = 0; j < nfv; ++j) {
            MVertex_ptr mv = List_Entry(fverts, j);
            MEnt_Get_AttVal(mv, copyatt, &ival, &rval, &pval);
            if (pval)
              fverts_new[j] = pval;
            else {
              fverts_new[j] = MV_New(mesh_);
              MV_Coords(mv, xyz);
              MV_Set_Coords(fverts_new[j], xyz);
              MV_Set_GEntDim(fverts_new[j], MV_GEntDim(mv));
              MV_Set_GEntID(fverts_new[j], MV_GEntID(mv));
              MEnt_Set_AttVal(mv, copyatt, ival, rval, fverts_new[j]);
              MEnt_Set_AttVal(fverts_new[j], vparentatt_, 0, 0.0, mv);
            }
          }
          List_Delete(fverts);

          rfaces_new[i] = MF_New(mesh_);
          MF_Set_Vertices(rfaces_new[i], nfv, fverts_new);
          MF_Set_GEntDim(rfaces_new[i], MF_GEntDim(mf));
          MF_Set_GEntID(rfaces_new[i], MF_GEntID(mf));
          rfdirs_new[i] = MR_FaceDir_i(mr, i);

          MEnt_Set_AttVal(mf, copyatt, ival, rval, rfaces_new[i]);
          MEnt_Set_AttVal(rfaces_new[i], fparentatt_, 0, 0.0, mf);
        }
      }
      List_Delete(rfaces);

      MRegion_ptr mr_new = MR_New(mesh_);
      MR_Set_Faces(mr_new, nrf, rfaces_new, rfdirs_new);
      MR_Set_GEntID(mr_new, MR_GEntID(mr));

      MEnt_Set_AttVal(mr, copyatt, ival, rval, mr_new);
      MEnt_Set_AttVal(mr_new, rparentatt_, 0, 0.0, mr);
    }

    break;
  }
  case MFACE: { // Extracting a surface from a solid mesh or subset of
    //           // a surface mesh
    idx = 0;

    MFace_ptr mf = nullptr;
    while ((mf = (MFace_ptr)List_Next_Entry(src_entities, &idx))) {
      List_ptr fedges = MF_Edges(mf, 1, 0);
      int nfe = List_Num_Entries(fedges);
      int fedirs[MAXPV2];
      MEdge_ptr fedges_new[MAXPV2];
      for (int j = 0; j < nfe; ++j) {
        MEdge_ptr me = List_Entry(fedges, j);
        MEnt_Get_AttVal(me, copyatt, &ival, &rval, &pval);
        if (pval)
          fedges_new[j] = pval;
        else {
          fedges_new[j] = ME_New(mesh_);

          for (int k = 0; k < 2; ++k) {
            MVertex_ptr mv = ME_Vertex(me, k);
            MVertex_ptr mv_new = nullptr;
            MEnt_Get_AttVal(mv, copyatt, &ival, &rval, &pval);
            if (pval) {
              mv_new = pval;
            } else {
              MV_Coords(mv, xyz);
              if (flatten) xyz[2] = 0.0;
              mv_new = MV_New(mesh_);
              MV_Set_Coords(mv_new, xyz);
              MV_Set_GEntDim(mv_new, MV_GEntDim(mv));
              MV_Set_GEntID(mv_new, MV_GEntID(mv));
              MEnt_Set_AttVal(mv, copyatt, ival, rval, mv_new);
              MEnt_Set_AttVal(mv_new, vparentatt_, 0, 0.0, mv);
            }

            ME_Set_Vertex(fedges_new[j], k, mv_new);
          }
          ME_Set_GEntDim(fedges_new[j], ME_GEntDim(me));
          ME_Set_GEntID(fedges_new[j], ME_GEntID(me));
          MEnt_Set_AttVal(me, copyatt, ival, rval, fedges_new[j]);
          MEnt_Set_AttVal(fedges_new[j], eparentatt_, 0, 0.0, me);
        }
        fedirs[j] = MF_EdgeDir_i(mf, j);
      }
      List_Delete(fedges);

      MFace_ptr mf_new = MF_New(mesh_);
      MF_Set_Edges(mf_new, nfe, fedges_new, fedirs);
      MF_Set_GEntDim(mf_new, 2); // This has to be surface mesh
      if (MF_GEntDim(mf) == 2) MF_Set_GEntID(mf_new, MF_GEntID(mf));

      MEnt_Set_AttVal(mf, copyatt, ival, rval, mf_new);
      MEnt_Set_AttVal(mf_new, fparentatt_, 0, 0.0, mf);
    }

    break;
  }
  case MEDGE: { // Extracting a wire mesh from a solid or surface mesh
    idx = 0;
    MEdge_ptr me = nullptr;
    while ((me = (MEdge_ptr)List_Next_Entry(src_entities, &idx))) {
      MEdge_ptr me_new = ME_New(mesh_);

      for (int j = 0; j < 2; ++j) {
        MVertex_ptr mv = ME_Vertex(me, j);

        MVertex_ptr mv_new = nullptr;

        MEnt_Get_AttVal(mv, copyatt, &ival, &rval, &pval);
        if (pval)
          mv_new = pval;
        else {
          MV_Coords(mv, xyz);
          if (flatten) {
            xyz[1] = 0.0;
            xyz[2] = 0.0;
          }
          mv_new = MV_New(mesh_);
          MV_Set_Coords(mv_new, xyz);
          MV_Set_GEntDim(mv_new, MV_GEntDim(mv));
          MV_Set_GEntID(mv_new, MV_GEntID(mv));

          MEnt_Set_AttVal(mv, copyatt, ival, rval, mv_new);
          MEnt_Set_AttVal(mv_new, vparentatt_, 0, 0.0, mv);
        }

        ME_Set_Vertex(me_new, j, mv_new);
      }

      if (ME_GEntDim(me) == 1) ME_Set_GEntDim(me_new, 1);
      MEnt_Set_AttVal(me, copyatt, ival, rval, me_new);
      MEnt_Set_AttVal(me_new, eparentatt_, 0, 0.0, me);
    }
    break;
  }
  case MVERTEX: {
    idx = 0;
    MVertex_ptr mv = nullptr;
    while ((mv = (MVertex_ptr)List_Next_Entry(src_entities, &idx))) {
      MVertex_ptr mv_new = MV_New(mesh_);
      MV_Set_Coords(mv_new, xyz);
      if (flatten) xyz[2] = 0.0;
      MV_Set_GEntDim(mv_new, MV_GEntDim(mv));
      MV_Set_GEntID(mv_new, MV_GEntID(mv));

      MEnt_Set_AttVal(mv, copyatt, ival, rval, mv_new);
      MEnt_Set_AttVal(mv_new, vparentatt_, 0, 0.0, mv);
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
    int status = MSTK_Weave_DistributedMeshes(
      mesh_, getManifoldDimension(), num_ghost_layers, input_type, mpicomm_);

    // Now we have to build parent information for global entities
    MAttrib_ptr vparentgid_att = MAttrib_New(mesh_, "vparent_gid", INT, MVERTEX);
    MAttrib_ptr eparentgid_att = MAttrib_New(mesh_, "eparent_gid", INT, MEDGE);
    MAttrib_ptr fparentgid_att = MAttrib_New(mesh_, "fparent_gid", INT, MFACE);
    MAttrib_ptr rparentgid_att = MAttrib_New(mesh_, "rparent_gid", INT, MREGION);

    // Attach parent global ID info to entities used by other processors
    int size = getComm()->NumProc();
    int rank = getComm()->MyPID();

    idx = 0;
    MVertex_ptr mv = nullptr;

    while ((mv = (MVertex_ptr)MESH_Next_Vertex(mesh_, &idx))) {
      auto ptype = MV_PType(mv);
      if (ptype == POVERLAP) {
        MEnt_Get_AttVal(mv, vparentatt_, &ival, &rval, &pval);
        MEnt_Set_AttVal(mv, vparentgid_att, MV_GlobalID((MVertex_ptr)pval), 0.0, NULL);
      }
    }

    MEdge_ptr me = nullptr;
    if (entity_dim !=
        MREGION) { // edge parents not set on 3D extraction -- maybe they should be --etc
      idx = 0;
      while ((me = (MEdge_ptr)MESH_Next_Edge(mesh_, &idx)))
        if (ME_PType(me) == POVERLAP) {
          MEnt_Get_AttVal(me, eparentatt_, &ival, &rval, &pval);
          MEnt_Set_AttVal(me, eparentgid_att, ME_GlobalID((MEdge_ptr)pval), 0.0, NULL);
        }
    }

    idx = 0;
    MFace_ptr mf = nullptr;
    while ((mf = (MFace_ptr)MESH_Next_Face(mesh_, &idx))) {
      auto ptype = MF_PType(mf);
      if (ptype == POVERLAP) {
        MEnt_Get_AttVal(mf, fparentatt_, &ival, &rval, &pval);
        MEnt_Set_AttVal(mf, fparentgid_att, MF_GlobalID((MFace_ptr)pval), 0.0, NULL);
      }
    }

    idx = 0;
    MRegion_ptr mr = nullptr;
    while ((mr = (MRegion_ptr)MESH_Next_Region(mesh_, &idx)))
      if (MR_PType(mr) == POVERLAP) {
        MEnt_Get_AttVal(mr, rparentatt_, &ival, &rval, &pval);
        MEnt_Set_AttVal(mr, rparentgid_att, MR_GlobalID((MRegion_ptr)pval), 0.0, NULL);
      }

    // Update attributes on ghost entities - this will ensure that
    // ghost entities have their parent global ID information
    status &= MESH_UpdateAttributes(mesh_, mpicomm_);

    // Now reverse engineer the parents of ghost entities from the global IDs
    idx = 0;
    while ((mv = (MVertex_ptr)MESH_Next_GhostVertex(mesh_, &idx))) {
      MEnt_Get_AttVal(mv, vparentgid_att, &ival, &rval, &pval);
      MVertex_ptr mv_parent = MESH_VertexFromGlobalID(parent_mesh_mstk, ival);
      if (!mv_parent) {
        Errors::Message mesg("Cannot find ghost vertex with given global ID");
        Exceptions::amanzi_throw(mesg);
      }
      MEnt_Set_AttVal(mv, vparentatt_, 0, 0.0, mv_parent);
    }
    if (entity_dim !=
        MREGION) { // edge parents not set on 3D extraction -- maybe they should be --etc
      idx = 0;
      while ((me = (MEdge_ptr)MESH_Next_GhostEdge(mesh_, &idx))) {
        MEnt_Get_AttVal(me, eparentgid_att, &ival, &rval, &pval);
        MEdge_ptr me_parent = MESH_EdgeFromGlobalID(parent_mesh_mstk, ival);
        if (!me_parent) {
          Errors::Message mesg("Cannot find ghost edge with given global ID");
          Exceptions::amanzi_throw(mesg);
        }
        MEnt_Set_AttVal(me, eparentatt_, 0, 0.0, me_parent);
      }
    }
    idx = 0;
    while ((mf = (MFace_ptr)MESH_Next_GhostFace(mesh_, &idx))) {
      MEnt_Get_AttVal(mf, fparentgid_att, &ival, &rval, &pval);
      MFace_ptr mf_parent = MESH_FaceFromGlobalID(parent_mesh_mstk, ival);
      if (!mf_parent) {
        Errors::Message mesg("Cannot find ghost face with given global ID");
        Exceptions::amanzi_throw(mesg);
      }
      MEnt_Set_AttVal(mf, fparentatt_, 0, 0.0, mf_parent);
    }
    idx = 0;
    while ((mr = (MRegion_ptr)MESH_Next_GhostRegion(mesh_, &idx))) {
      MEnt_Get_AttVal(mr, rparentgid_att, &ival, &rval, &pval);
      MRegion_ptr mr_parent = MESH_RegionFromGlobalID(parent_mesh_mstk, ival);
      if (!mr_parent) {
        Errors::Message mesg("Cannot find ghost region with given global ID");
        Exceptions::amanzi_throw(mesg);
      }
      MEnt_Set_AttVal(mr, rparentatt_, 0, 0.0, mr_parent);
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
    while ((mr = (MRegion_ptr)List_Next_Entry(src_entities, &idx))) {
      List_ptr rfaces = MR_Faces(mr);
      int nrf = List_Num_Entries(rfaces);

      for (int i = 0; i < nrf; ++i) {
        MFace_ptr mf = List_Entry(rfaces, i);
        MEnt_Rem_AttVal(mf, copyatt);

        List_ptr fverts = MF_Vertices(mf, 1, 0);
        int nfv = List_Num_Entries(fverts);

        for (int j = 0; j < nfv; ++j) {
          MVertex_ptr mv = List_Entry(fverts, j);
          MEnt_Rem_AttVal(mv, copyatt);
        }
        List_Delete(fverts);

        MEnt_Rem_AttVal(mf, copyatt);
      }
      List_Delete(rfaces);

      MEnt_Rem_AttVal(mr, copyatt);
    }
    break;
  }
  case MFACE: {
    MFace_ptr mf = nullptr;
    idx = 0;
    while ((mf = (MFace_ptr)List_Next_Entry(src_entities, &idx))) {
      List_ptr fedges = MF_Edges(mf, 1, 0);
      int nfe = List_Num_Entries(fedges);
      for (int j = 0; j < nfe; ++j) {
        MEdge_ptr me = List_Entry(fedges, j);
        MEnt_Rem_AttVal(me, copyatt);
        MVertex_ptr mv = ME_Vertex(me, MF_EdgeDir_i(mf, j));
        MEnt_Rem_AttVal(mv, copyatt);
      }
      List_Delete(fedges);

      MEnt_Rem_AttVal(mf, copyatt);
    }

    break;
  }
  case MEDGE: {
    MEdge_ptr me;
    idx = 0;
    while ((me = (MEdge_ptr)List_Next_Entry(src_entities, &idx))) {
      for (int j = 0; j < 2; ++j) {
        MVertex_ptr mv = ME_Vertex(me, j);
        MEnt_Rem_AttVal(mv, copyatt);
      }
      MEnt_Rem_AttVal(me, copyatt);
    }

    break;
  }
  case MVERTEX: {
    MVertex_ptr mv = nullptr;
    idx = 0;
    while ((mv = (MVertex_ptr)List_Next_Entry(src_entities, &idx))) MEnt_Rem_AttVal(mv, copyatt);

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
Mesh_MSTK::writeToExodusFile(const std::string& filename) const
{
  MESH_ExportToExodusII(mesh_, filename.c_str(), -1, NULL, NULL, mpicomm_);
}


// Run MSTK's internal checks - meant for debugging only
// Returns true if everything is ok, false otherwise
bool
Mesh_MSTK::run_internal_mstk_checks() const
{
  return MESH_CheckTopo(mesh_) && MESH_Parallel_Check(mesh_, mpicomm_);
}


} // namespace AmanziMesh
} // namespace Amanzi

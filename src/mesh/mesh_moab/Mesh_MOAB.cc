/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, Konstantin Lipnikov, others
*/

//! Implementation of the Mesh interface leveraging MOAB.
#include "dbc.hh"
#include "errors.hh"

#include "RegionLabeledSet.hh"
#include "RegionPoint.hh"
#include "RegionLogical.hh"
#include "Mesh_MOAB.hh"

using namespace moab;

namespace Amanzi {
namespace AmanziMesh {

//--------------------------------------------------------------------
// Constructor - load up mesh from file
//--------------------------------------------------------------------
Mesh_MOAB::Mesh_MOAB(const std::string& filename,
                     const Comm_ptr_type& comm,
                     const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
                     const Teuchos::RCP<const Teuchos::ParameterList>& plist,
                     const bool request_faces,
                     const bool request_edges)
    : Mesh(comm, gm, plist, request_faces, request_edges),
      extface_map_w_ghosts_(NULL), extface_map_wo_ghosts_(NULL)
{
  int result, rank;
  
  clear_internals_();

  // Core MOAB object
  mbcore_ = std::unique_ptr<moab::Core>(new moab::Core());

  if (comm_.get()) {
    auto mpi_comm = Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm_);
    if (mpi_comm.get()) {    
      // MOAB's parallel communicator
      int mbcomm_id;
      mbcomm_ = std::unique_ptr<ParallelComm>(new ParallelComm(&*mbcore_, mpi_comm->GetMpiComm(), &mbcomm_id));

      if (!mbcomm_.get()) {
        Errors::Message message("Failed to initialize MOAB communicator");
        Exceptions::amanzi_throw(message);
      }
    }
  }

  if (!mbcomm_.get() || mbcomm_->size() == 1) 
    serial_run = true;
  else
    serial_run = false;

  if (!serial_run) {
    // Load partitioned mesh - serial read of mesh with partition
    // info, deletion of non-local entities, resolution of
    // interprocessor connections. If we need ghosts we have to add
    // the option "PARALLEL_GHOSTS=A.B.C.D" where A is usually the
    // topological dimension of the mesh cells, B is the bridge
    // dimension or the dimension of entities across which we want
    // ghosts (0 for vertex connected ghost cells) and C indicates
    // the number of layers of ghost cells we want, D indicates if
    // we want edges/faces bounding the ghost cells as well (1 for
    // edges, 2 for faces, 3 for faces and edges)

    // In the specification for the Ghosts we made the assumption 
    // that we are dealing with 3D meshes only

    result = mbcore_->load_file(
        filename.c_str(), NULL,
        "PARALLEL=READ_DELETE;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION=PARALLEL_PARTITION;PARALLEL_GHOSTS=3.0.1.2",
        NULL, NULL, 0);
      
    rank = mbcomm_->rank();
  }
  else {
    // Load serial mesh
    result = mbcore_->load_file(filename.c_str(), 0, 0, NULL, NULL, 0);
    rank = 0;
  }

  if (result != MB_SUCCESS) {
    std::cerr << "FAILED" << std::endl;
    std::cerr << "Failed to load " << filename << " on processor " << rank << std::endl;
    std::cerr << "MOAB error code " << result << std::endl;
    AMANZI_ASSERT(result == MB_SUCCESS);
  }
      
      
  // Dimension of space, mesh cells, faces etc
  int ndim;  
  result = mbcore_->get_dimension(ndim);
  space_dim_ = ndim;
    

  // Highest topological dimension
  int nent;
  result = mbcore_->get_number_entities_by_dimension(0, 3, nent, false);
  ErrorCheck_(result, "Problem getting number of entities of dim 3");

  if (nent) {
    celldim = 3;
    facedim = 2;
  }
  else {
    result = mbcore_->get_number_entities_by_dimension(0, 2, nent, false);
    ErrorCheck_(result, "Problem getting number of entities of dim 2");

    if (nent) {
      celldim = 2;
      facedim = 1;
    }
    else {
      std::cerr << "Flow code works only on 2D and 3D meshes" << std::endl;
      AMANZI_ASSERT(nent > 0);
    }
  }

  setManifoldDimension(celldim);

  // redefine space dimension  // FIXME
  mbcore_->set_dimension(celldim);
  space_dim_ = celldim;
      
  { // Keep together and in this order 
    init_pvert_lists();
    init_pcell_lists(); // cells MUST be initialized before faces
    init_pface_lists();
    
    // Create maps from local IDs to MOAB entity handles (must be after
    // the various init_p*_list calls)
    init_id_handle_maps();
  }
    
  init_global_ids();

  init_pface_dirs();

  // Create Epetra_maps
  init_cell_map();
  init_face_map();
  init_node_map();
  
  // Initialize some info about the global number of sets, global set
  // IDs and set types
  if (geometric_model() != Teuchos::null) init_set_info();
}


//--------------------------------------------------------------------
// Clean up 
//--------------------------------------------------------------------
Mesh_MOAB::~Mesh_MOAB() {
  delete cell_map_wo_ghosts_;
  delete cell_map_w_ghosts_;
  delete face_map_wo_ghosts_;
  delete face_map_w_ghosts_;
  delete node_map_wo_ghosts_;
  delete node_map_w_ghosts_;
  if (extface_map_wo_ghosts_) delete extface_map_wo_ghosts_;
  if (extface_map_w_ghosts_) delete extface_map_w_ghosts_;
  delete [] setids_;
  delete [] setdims_;
  delete [] faceflip;
}


//--------------------------------------------------------------------
// Some initializations
//--------------------------------------------------------------------
void Mesh_MOAB::clear_internals_() 
{ 
  mbcore_ = nullptr;
  mbcomm_ = nullptr;

  AllVerts.clear();
  OwnedVerts.clear();
  NotOwnedVerts.clear();
  AllFaces.clear();
  OwnedFaces.clear();
  NotOwnedFaces.clear();
  AllCells.clear();
  OwnedCells.clear();
  GhostCells.clear();
    
  lid_tag = 0;
  gid_tag = 0;
  cstag = 0;
  sstag = 0;
  nstag = 0;

  celldim = -1;
  facedim = -1;

  faceflip = NULL;

  cell_map_w_ghosts_ = cell_map_wo_ghosts_ = NULL;
  face_map_w_ghosts_ = face_map_wo_ghosts_ = NULL;
  node_map_w_ghosts_ = node_map_wo_ghosts_ = NULL;

  nsets = 0;
  setids_ = setdims_ = NULL;
}


//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
void Mesh_MOAB::init_id_handle_maps()
{
  int i, nv, nf, nc;
  int result;

  // Assign local IDs to entities
  // -- nodes
  int tagval = 0;
  result = mbcore_->tag_get_handle("LOCAL_ID", 1,
                                  MB_TYPE_INTEGER, lid_tag,
                                  MB_TAG_CREAT|MB_TAG_DENSE, &tagval);
  ErrorCheck_(result, "Problem getting tag handle for LOCAL_ID");
      
  nv = AllVerts.size();
  node_id_to_handle.reserve(nv);

  i = 0;
  for (auto it = OwnedVerts.begin(); it != OwnedVerts.end(); ++it) {
    moab::EntityHandle handle = *it;
    result = mbcore_->tag_set_data(lid_tag, &handle, 1, &i);
    ErrorCheck_(result, "Problem getting local ID for vertex");
    node_id_to_handle[i++] = handle;
  }    

  for (auto it = NotOwnedVerts.begin(); it != NotOwnedVerts.end(); ++it) {
    moab::EntityHandle handle = *it;
    result = mbcore_->tag_set_data(lid_tag, &handle, 1, &i);
    ErrorCheck_(result, "Problem getting local ID for vertex");
    node_id_to_handle[i++] = handle;
  }
    
  // -- faces
  nf = AllFaces.size();

  face_id_to_handle.reserve(nf);

  i = 0;
  for (auto it = OwnedFaces.begin(); it != OwnedFaces.end(); ++it) {
    moab::EntityHandle face = *it;
    result = mbcore_->tag_set_data(lid_tag, &face, 1, &i);
    ErrorCheck_(result, "Problem getting local ID for face");
    face_id_to_handle[i++] = face;
  }

  for (auto it = NotOwnedFaces.begin(); it != NotOwnedFaces.end(); ++it) {
    moab::EntityHandle face = *it;
    result = mbcore_->tag_set_data(lid_tag, &face, 1, &i);
    ErrorCheck_(result, "Problem getting local ID for face");
    face_id_to_handle[i++] = face;
  }
    
  // -- cells
  nc = AllCells.size();

  cell_id_to_handle.reserve(nc);

  i = 0;
  for (auto it = OwnedCells.begin(); it != OwnedCells.end(); ++it) {
    moab::EntityHandle handle = *it;
    result = mbcore_->tag_set_data(lid_tag, &handle, 1, &i);
    ErrorCheck_(result, "Problem getting local ID for cell");
    cell_id_to_handle[i++] = handle;
  }

  for (auto it = GhostCells.begin(); it != GhostCells.end(); ++it) {
    moab::EntityHandle handle = *it;
    result = mbcore_->tag_set_data(lid_tag, &handle, 1, &i);
    ErrorCheck_(result, "Problem getting local ID for cell");
    cell_id_to_handle[i++] = handle;
  }
}


//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
void Mesh_MOAB::init_global_ids() {
  int result;
    
  if (mbcomm_) {
    // Ask Parallel Communicator to assign global IDs to entities
    bool largest_dim_only=false;
    int start_id = 0;
    int largest_dim=celldim;
    result = mbcomm_->assign_global_ids(0, largest_dim, start_id, largest_dim_only);
    ErrorCheck_(result, "Problem assigning global IDs");
    
    // Exchange global IDs across all processors
    result = mbcore_->tag_get_handle("GLOBAL_ID", gid_tag);
    ErrorCheck_(result, "Could not get tag handle for GLOBAL_ID data");

    mbcomm_->exchange_tags(gid_tag, AllVerts);
    mbcomm_->exchange_tags(gid_tag, AllFaces);
    mbcomm_->exchange_tags(gid_tag, AllCells);
  }
  else {
    // Case without MPI - we assign global IDs ourselves
    int tagval = 0;
    result = mbcore_->tag_get_handle("GLOBAL_ID", 1,
                                    MB_TYPE_INTEGER, gid_tag,
                                    MB_TAG_CREAT|MB_TAG_DENSE, &tagval);
    ErrorCheck_(result, "Problem getting tag handle for GLOBAL_ID");
      
    int nent = AllVerts.size();
    int *gids = new int[nent];
    for (int i = 0; i < nent; i++) gids[i] = i;

    result = mbcore_->tag_set_data(gid_tag, AllVerts, gids);
    ErrorCheck_(result, "Problem setting global IDs for vertices");

    delete [] gids;

    nent = AllFaces.size();
    gids = new int[nent];
    for (int i = 0; i < nent; i++) gids[i] = i;

    result = mbcore_->tag_set_data(gid_tag, AllFaces, gids);
    ErrorCheck_(result, "Problem setting global IDs for faces");

    delete [] gids;

    nent = AllCells.size();
    gids = new int[nent];
    for (int i = 0; i < nent; i++) gids[i] = i;

    result = mbcore_->tag_set_data(gid_tag, AllCells, gids);
    ErrorCheck_(result, "Problem setting global IDs for cells");

    delete [] gids;
  }
}


//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
void Mesh_MOAB::init_pvert_lists() {
  int result;

  // Get all vertices on this processor 
  result = mbcore_->get_entities_by_dimension(0, 0, AllVerts, false);
  ErrorCheck_(result, "Could not get vertices");

  // Get not owned vertices
  result = mbcomm_->get_pstatus_entities(0, PSTATUS_NOT_OWNED, NotOwnedVerts);
  ErrorCheck_(result, "Could not get NotOwned vertices");

  // Subtract from all vertices on processor to get owned vertices only
  OwnedVerts = AllVerts;  // I think we DO want a data copy here
  OwnedVerts -= NotOwnedVerts;
}


//--------------------------------------------------------------------
// init_pface_lists is more complicated than init_pvert_lists and
// init_pcell_lists because of the way MOAB is setting up shared
// entities and ghost entities. When we ask MOAB to resolve shared
// entities, then MOAB sets up faces on interprocessor boundaries and
// assigns each of them to some processor. Therefore, the pstatus tags
// on these faces are correctly set. On the other hand when we ask for
// ghost cells, MOAB does not automatically create ghost faces. Also,
// when we go through ghost cells and create their faces, MOAB does
// not tag them as ghost faces, tagging them as owned faces
// instead. So we have to process them specially.
//--------------------------------------------------------------------
void Mesh_MOAB::init_pface_lists() {
  int result;

  // Make MOAB create the missing 'faces' (faces in 3D, edges in
  // 2D). We do this by looping over the cells and asking for their
  // faces with the create_if_missing=true option
  for (auto it = AllCells.begin(); it != AllCells.end(); it++) {
    moab::EntityHandle cell = *it;
    moab::Range cfaces;
      
    result = mbcore_->get_adjacencies(&cell, 1, facedim, true, cfaces, Core::UNION);
    ErrorCheck_(result, "Could not get faces of cell");
  }

  // Get all "faces" (edges in 2D, faces in 3D) on this processor
  result = mbcore_->get_entities_by_dimension(0, facedim, AllFaces, false);
  ErrorCheck_(result, "Could not get 'faces'");

  // Get not owned faces 
  result = mbcomm_->get_pstatus_entities(facedim, PSTATUS_NOT_OWNED, NotOwnedFaces);
  ErrorCheck_(result, "Could not get NotOwned 'faces'");

  // Subtract from all faces on processor to get owned faces only
  OwnedFaces = AllFaces;  // I think we DO want a data copy here
  OwnedFaces -= NotOwnedFaces;
}


//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
void Mesh_MOAB::init_pface_dirs()
{
  int result, zero(0);
  int face_lid, face_gid, cell_gid;
  int sidenum, offset, facedir;
  moab::EntityHandle cell, face;

  // Do some additional processing to see if ghost faces and their masters
  // are oriented the same way; if not, turn on flag to flip them

  // In this code, we increment local values of global IDs by 1 so
  // that we can distinguish between the lowest gid and no data
  moab::Tag tmp_fc0_tag, tmp_fc1_tag;
  result = mbcore_->tag_get_handle("TMP_FC0_TAG", 1,
                                  MB_TYPE_INTEGER, tmp_fc0_tag,
                                  MB_TAG_CREAT | MB_TAG_DENSE, &zero);
  ErrorCheck_(result, "Problem getting new tag handle");

  result = mbcore_->tag_get_handle("TMP_FC1_TAG", 1,
                                  MB_TYPE_INTEGER, tmp_fc1_tag,
                                  MB_TAG_CREAT | MB_TAG_DENSE, &zero);
  ErrorCheck_(result, "Problem getting new tag handle");

  for (auto it = OwnedFaces.begin(); it != OwnedFaces.end(); it++) {
    moab::Range fcells;
    face = *it;

    result = mbcore_->get_adjacencies(&face, 1, celldim, false, fcells, Core::UNION);
    ErrorCheck_(result, "Could not get cells of face");

    result = mbcore_->tag_set_data(tmp_fc0_tag, &face, 1, &zero);
    ErrorCheck_(result, "Problem setting tag data");

    result = mbcore_->tag_set_data(tmp_fc1_tag, &face, 1, &zero);
    ErrorCheck_(result, "Problem setting tag data");

    for (auto jt = fcells.begin(); jt != fcells.end(); ++jt) {
      cell = *jt;
      result = mbcore_->side_number(cell, face, sidenum, facedir, offset);
      ErrorCheck_(result, "Could not get face dir w.r.t. cell");

      result = mbcore_->tag_get_data(gid_tag, &cell, 1, &cell_gid);
      ErrorCheck_(result, "Problem getting tag data");

      cell_gid += 1; 

      if (facedir == 1) {
        result = mbcore_->tag_set_data(tmp_fc0_tag, &face, 1, &cell_gid);
        ErrorCheck_(result, "Problem setting tag data");
      }
      else {
        result = mbcore_->tag_set_data(tmp_fc1_tag, &face, 1, &cell_gid);
        ErrorCheck_(result, "Problem setting tag data");
      }
    }
  }

  result = mbcomm_->exchange_tags(tmp_fc0_tag, AllFaces);
  ErrorCheck_(result, "Could not get exchange tag data successfully");

  result = mbcomm_->exchange_tags(tmp_fc1_tag, AllFaces);
  ErrorCheck_(result, "Could not get exchange tag data successfully");

  faceflip = new bool[AllFaces.size()];
  for (int i = 0; i < AllFaces.size(); i++) faceflip[i] = false;

  for (auto it = NotOwnedFaces.begin(); it != NotOwnedFaces.end(); it++) {
    moab::Range fcells;
    int ghost_cell0_gid = 0, ghost_cell1_gid = 0;
    int master_cell0_gid = 0, master_cell1_gid = 0;

    face = *it;

    result = mbcore_->tag_get_data(tmp_fc0_tag, &face, 1, &master_cell0_gid);
    ErrorCheck_(result, "Could not get face tag data");

    result = mbcore_->tag_get_data(tmp_fc1_tag, &face, 1, &master_cell1_gid);
    ErrorCheck_(result, "Could not get face tag data");

    result = mbcore_->get_adjacencies(&face, 1, celldim, false, fcells, Core::UNION);
    ErrorCheck_(result, "Could not get cells of face");

    for (auto jt = fcells.begin(); jt != fcells.end(); ++jt) {
      cell = *jt;

      result = mbcore_->side_number(cell, face, sidenum, facedir, offset);
      ErrorCheck_(result, "Could not get face dir w.r.t. cell");

      if (facedir == 1) {
        result = mbcore_->tag_get_data(gid_tag, &cell, 1, &ghost_cell0_gid);
        ErrorCheck_(result, "Problem getting tag data");
        ghost_cell0_gid += 1;
      }
      else {
        result = mbcore_->tag_get_data(gid_tag, &cell, 1, &ghost_cell1_gid);
        ErrorCheck_(result, "Problem getting tag data");
        ghost_cell1_gid += 1;
      }
    }

    if (ghost_cell0_gid == master_cell1_gid || 
        ghost_cell1_gid == master_cell0_gid) {

      // Both cells don't have to match because a ghost face may 
      // not have the cell on the other side
      result = mbcore_->tag_get_data(lid_tag, &face, 1, &face_lid);
      ErrorCheck_(result, "Could not get face tag data");
      faceflip[face_lid] = true;

    }
    else { // Sanity check
      if (ghost_cell0_gid != master_cell0_gid &&
          ghost_cell1_gid != master_cell1_gid) {

        // Problem if there is no match at all
        result = mbcore_->tag_get_data(gid_tag, &face, 1, &face_gid);
        ErrorCheck_(result, "Problem getting tag data");

        std::cout << "Face cells mismatch between master and ghost (processor " 
                  << mbcomm_->rank() << ")" << std::endl;
        std::cout << " Face " << face_gid << std::endl;
        std::cout << "Master cells " << master_cell0_gid << " " << master_cell1_gid << std::endl;
        std::cout << "Ghost cells " << ghost_cell0_gid << " " << ghost_cell1_gid << std::endl;
      }
    }
  }
}


//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
void Mesh_MOAB::init_pcell_lists() {
  int result;

  // Get all cells (faces in 2D, regions in 3D) on this processor
  result = mbcore_->get_entities_by_dimension(0, celldim, AllCells, false);
  ErrorCheck_(result, "Could not get cells");
    
  // Get not owned cells (which is the same as ghost cells)
  result = mbcomm_->get_pstatus_entities(celldim, PSTATUS_GHOST, GhostCells);
  ErrorCheck_(result, "Could not get ghost cells");
    
  // Subtract from all cells on processor to get owned cells only
  OwnedCells = AllCells;  // I think we DO want a data copy here
  OwnedCells -= GhostCells;
}


//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
void Mesh_MOAB::init_set_info()
{
  int result;
  moab::Tag tag;

  // Get element block, sideset and nodeset tags
  result = mbcore_->tag_get_handle(MATERIAL_SET_TAG_NAME, cstag);
  ErrorCheck_(result, "Could not get tag for material sets");

  result = mbcore_->tag_get_handle(NEUMANN_SET_TAG_NAME, 1, MB_TYPE_INTEGER, sstag);
  ErrorCheck_(result, "Could not get tag for side sets");

  result = mbcore_->tag_get_handle(DIRICHLET_SET_TAG_NAME, nstag);
  ErrorCheck_(result, "Could not get tag for node sets");

  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm = Mesh::geometric_model();

  if (gm == Teuchos::null) { 
    Errors::Message mesg("Need region definitions to initialize sets");
    amanzi_throw(mesg);
  }
    
  unsigned int ngr = gm->size();

  for (int i = 0; i < ngr; i++) {
    Teuchos::RCP<const AmanziGeometry::Region> rgn = gm->FindRegion(i);

    if (rgn->type() == AmanziGeometry::LABELEDSET) {
      Teuchos::RCP<const AmanziGeometry::RegionLabeledSet> lsrgn =
          Teuchos::rcp_static_cast<const AmanziGeometry::RegionLabeledSet>(rgn);

      std::string internal_name;
      std::string label = lsrgn->label();
      std::string entity_type_str = lsrgn->entity_str();

      if (entity_type_str == "CELL")
        internal_name = internal_name_of_set(rgn, CELL);
      else if (entity_type_str == "FACE")
        internal_name = internal_name_of_set(rgn, FACE);
      else if (entity_type_str == "NODE")
        internal_name = internal_name_of_set(rgn, NODE);

      result = mbcore_->tag_get_handle(internal_name.c_str(), 1, MB_TYPE_INTEGER,
                                      tag, MB_TAG_CREAT|MB_TAG_SPARSE);
      ErrorCheck_(result, "Problem getting labeled set");
    }
    //    else { /* General region - we have to account for all kinds of
    //              entities being queried in a set defined by this 
    //              region */
    //      Entity_kind int_to_kind[3] = {NODE, FACE, CELL};
    //
    //      for (int k = 0; k < 3; k++) {
    //        Entity_kind kind = int_to_kind[k];
    //
    //      std::string internal_name = internal_name_of_set(rgn, kind);
    //
    //  result = mbcore_->tag_get_handle(internal_name.c_str(), 1, MB_TYPE_INTEGER,
    //                                  tag, MB_TAG_CREAT|MB_TAG_SPARSE);
    //    if (result != MB_SUCCESS) {
    //      std::cerr << "Could not create tag with name " << rgn->name() << std::endl;
    //      AMANZI_ASSERT(result != MB_SUCCESS);
    //    }
    //  }
    //}
  }
}


//--------------------------------------
// Number of OWNED, GHOST or ALL entities of different types
//--------------------------------------
unsigned int Mesh_MOAB::num_entities(Entity_kind kind, 
                                     Parallel_kind ptype) const
{
  switch (kind) {
  case NODE:
    switch (ptype) {
    case Parallel_kind::OWNED:
      return !serial_run ? OwnedVerts.size() : AllVerts.size();
      break;
    case Parallel_kind::GHOST:
      return !serial_run ? NotOwnedVerts.size() : 0;
      break;
    case Parallel_kind::ALL:
      return AllVerts.size();
      break;
    default:
      return 0;
    }
    break;

  case FACE:
    switch (ptype) {
    case Parallel_kind::OWNED:
      return !serial_run ? OwnedFaces.size() : AllFaces.size();
      break;
    case Parallel_kind::GHOST:
      return !serial_run ? NotOwnedFaces.size() : 0;
      break;
    case Parallel_kind::ALL:
      return AllFaces.size();
      break;
    default:
      return 0;
    }
    break;

  case CELL:
    switch (ptype) {
    case Parallel_kind::OWNED:
      return !serial_run ? OwnedCells.size() : AllCells.size();
      break;
    case Parallel_kind::GHOST:
      return !serial_run ? GhostCells.size() : 0;
      break;
    case Parallel_kind::ALL:
      return AllCells.size();
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


//--------------------------------------------------------------------
// Get faces of a cell and directions in which the cell uses the face 
//
// On a distributed mesh, this will return all the faces of the
// cell, OWNED or GHOST. If ordered = true, the faces will be
// returned in a standard order according to Exodus II convention
// for standard cells; in all other situations (ordered = false or
// non-standard cells), the list of faces will be in arbitrary order
//
// In 3D, direction is 1 if face normal points out of cell
// and -1 if face normal points into cell
// In 2D, direction is 1 if face/edge is defined in the same
// direction as the cell polygon, and -1 otherwise
//--------------------------------------------------------------------
void Mesh_MOAB::cell_get_faces_and_dirs_internal_(const Entity_ID cellid,
                                                  Entity_ID_View *faceids,
                                                  std::vector<int> *face_dirs,
                                                  const bool ordered) const
{
  moab::EntityHandle cell;
  moab::Range cell_faces;
  std::vector<moab::EntityHandle> cell_nodes, face_nodes;
  int *cell_faceids, *cell_facedirs;
  int nf, result;
  int cfstd[6][4] = {{0,1,5,4},    // Expected cell-face-node pattern
                     {1,2,6,5},
                     {2,3,7,6},
                     {0,4,7,3},
                     {0,3,2,1},
                     {4,5,6,7}};

  cell = cell_id_to_handle[cellid];

  result = mbcore_->get_adjacencies(&cell, 1, facedim, true, cell_faces, 
                                   moab::Interface::INTERSECT);
  ErrorCheck_(result, "Problem getting faces of cell");

  nf = cell_faces.size();
  faceids->resize(nf);
  if (face_dirs) face_dirs->resize(nf);

  cell_faceids = new int[nf];                        
  if (face_dirs) cell_facedirs = new int[nf];

  // Have to re-sort the faces according a specific template for hexes
  if (ordered && nf == 6) { // Hex
    moab::EntityHandle *ordfaces, face;

    ordfaces = new moab::EntityHandle[6];

    result = mbcore_->get_connectivity(&cell, 1, cell_nodes);
    ErrorCheck_(result, "Problem getting nodes of cell");

    for (int i = 0; i < nf; i++) {
      // Search for a face that has all the expected nodes
      bool found = false;
      int j;
      for (j = 0; j < nf; j++) {
        face = cell_faces[j];
        result = mbcore_->get_connectivity(&face, 1, face_nodes);
        ErrorCheck_(result, "Problem getting nodes of face");

        // Check if this face has all the expected nodes
        bool all_present = true;

        for (int k = 0; k < 4; k++) {
          Entity_ID node = cell_nodes[cfstd[i][k]];

          if (face_nodes[0] != node && face_nodes[1] != node &&
              face_nodes[2] != node && face_nodes[3] != node) {
            all_present = false;
            break;
          }
        }

        if (all_present) {
          found = true;
          break;
        }
      }

      AMANZI_ASSERT(found);

      if (found)
        ordfaces[i] = face;
    }

    result = mbcore_->tag_get_data(lid_tag, ordfaces, 6, cell_faceids);
    ErrorCheck_(result, "Problem getting tag data");

    if (face_dirs) {
      for (int i = 0; i < nf; i++) {
        face = ordfaces[i];
        int sidenum, offset;
        
        result = mbcore_->side_number(cell, face, sidenum, cell_facedirs[i], offset);
        ErrorCheck_(result, "Could not find face dir in cell");
        
        // If this is a ghost face and the master has the opposite direction
        // we are supposed to flip it
        if (faceflip[cell_faceids[i]]) cell_facedirs[i] *= -1;
      }
    }

    delete [] ordfaces;
  }
  else {
    result = mbcore_->tag_get_data(lid_tag, cell_faces, cell_faceids);
    ErrorCheck_(result, "Problem getting tag data");

    if (face_dirs) {
      for (int i = 0; i < nf; i++) {
        moab::EntityHandle face = cell_faces[i];
        int sidenum, offset;
        
        result = mbcore_->side_number(cell, face, sidenum, cell_facedirs[i], offset);
        ErrorCheck_(result, "Could not find face dir in cell");
        
        // If this is a ghost face and the master has the opposite direction
        // we are supposed to flip it
        if (faceflip[cell_faceids[i]]) cell_facedirs[i] *= -1;
      }
    }
  }

  Entity_ID_View::iterator itf = faceids->begin();
  for (int i = 0; i < nf; i++) {
    *itf = cell_faceids[i];
    ++itf;
  }
  if (face_dirs) {
    std::vector<int>::iterator itd = face_dirs->begin();
    for (int i = 0; i < nf; i++) {
      *itd = cell_facedirs[i];
      ++itd;
    }
  }

  delete [] cell_faceids;
  if (face_dirs) delete [] cell_facedirs;
}


//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
void Mesh_MOAB::cell_get_nodes(Entity_ID cellid, Entity_ID_View *cnodes) const
{
  moab::EntityHandle cell;
  std::vector<moab::EntityHandle> cell_nodes;
  int *cell_nodeids;
  int nn, result;

  cell = cell_id_to_handle[cellid];
      
  result = mbcore_->get_connectivity(&cell, 1, cell_nodes);
  ErrorCheck_(result, "Problem getting nodes of cell");

  nn = cell_nodes.size();
  cell_nodeids = new int[nn];

  cnodes->resize(nn);
    
  for (int i = 0; i < nn; i++) {
    result = mbcore_->tag_get_data(lid_tag, &(cell_nodes[i]), 1, &(cell_nodeids[i]));
    ErrorCheck_(result, "Problem getting tag data");
  }

  Entity_ID_View::iterator itn = cnodes->begin();
  for (int i = 0; i < nn; i++) {
    *itn = cell_nodeids[i];
    ++itn;
  }

  delete [] cell_nodeids;
}


//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
void Mesh_MOAB::face_get_nodes(Entity_ID faceid, Entity_ID_View *fnodes) const 
{
  moab::EntityHandle face;
  std::vector<moab::EntityHandle> face_nodes;
  int *face_nodeids;
  int nn, result;

  face = face_id_to_handle[faceid];

  result = mbcore_->get_connectivity(&face, 1, face_nodes, true);
  ErrorCheck_(result, "Problem getting nodes of face");

  nn = face_nodes.size();

  face_nodeids = new int[nn];
  if (faceflip[faceid]) {
    for (int i = nn-1; i >= 0; i--) {
      result = mbcore_->tag_get_data(lid_tag, &(face_nodes[i]), 1, &(face_nodeids[nn-i-1]));
      ErrorCheck_(result, "Problem getting tag data");
    }
  }
  else {
    for (int i = 0; i < nn; i++) {
      result = mbcore_->tag_get_data(lid_tag, &(face_nodes[i]), 1, &(face_nodeids[i]));
      ErrorCheck_(result, "Problem getting tag data");
    }
  }

  fnodes->resize(nn);
  Entity_ID_View::iterator itn = fnodes->begin();
  for (int i = 0; i < nn; i++) {
    *itn = face_nodeids[i];
    ++itn;
  }

  delete [] face_nodeids;
}
  

//--------------------------------------------------------------------
// Copy node coordinates to Point
//--------------------------------------------------------------------
void Mesh_MOAB::node_get_coordinates(Entity_ID node_id, AmanziGeometry::Point *ncoord) const 
{
  moab::EntityHandle node;
  double coords[3];

  node = node_id_to_handle[node_id];

  int result = mbcore_->get_coords(&node, 1, coords);
  ErrorCheck_(result, "Problem getting node coordinates");

  ncoord->set(space_dim_, coords);
}


//--------------------------------------------------------------------
// Modify node coordinates
//--------------------------------------------------------------------
void Mesh_MOAB::node_set_coordinates(const AmanziMesh::Entity_ID nodeid, 
                                     const double *coords)
{
  moab::EntityHandle v = node_id_to_handle[nodeid];

  int result = mbcore_->set_coords(&v, 1, coords);
  ErrorCheck_(result, "Problem setting node coordinates");
}


//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
void Mesh_MOAB::node_set_coordinates(const AmanziMesh::Entity_ID nodeid, 
                                     const AmanziGeometry::Point coords)
{
  moab::EntityHandle v = node_id_to_handle[nodeid];

  double coordarray[3] = {0.0, 0.0, 0.0};

  for (int i = 0; i < space_dim_; i++)
    coordarray[i] = coords[i];

  int result = mbcore_->set_coords(&v, 1, coordarray);
  ErrorCheck_(result, "Problem setting node coordinates");
}


//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
void Mesh_MOAB::cell_get_coordinates(
    Entity_ID cellid, Point_List *ccoords) const
{
  moab::EntityHandle cell;
  std::vector<moab::EntityHandle> cell_nodes;
  double *coords;
  int nn, result;

  cell = cell_id_to_handle[cellid];

  ccoords->clear();
      
  result = mbcore_->get_connectivity(&cell, 1, cell_nodes);
  ErrorCheck_(result, "Problem getting nodes of a cell");

  nn = cell_nodes.size();

  coords = new double[3];
  
  ccoords->resize(nn);
  Point_List::iterator it = ccoords->begin();

  for (int i = 0; i < nn; i++) {
    result = mbcore_->get_coords(&(cell_nodes[i]), 1, coords);
    ErrorCheck_(result, "Problem getting coordinates of a node");

    it->set(space_dim_, coords);
    ++it;
  }

  delete [] coords;
}


//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
void Mesh_MOAB::face_get_coordinates(Entity_ID faceid,
                                     Point_List *fcoords) const
{
  moab::EntityHandle face;
  std::vector<moab::EntityHandle> face_nodes;
  double *coords;
  int nn, result;

  face = face_id_to_handle[faceid];

  fcoords->clear();

  result = mbcore_->get_connectivity(&face, 1, face_nodes, true);
  ErrorCheck_(result, "Problem getting nodes of face");

  nn = face_nodes.size();

  coords = new double[3];
    
  fcoords->resize(nn);
  Point_List::iterator it = fcoords->begin();

  if (faceflip[faceid]) {
    for (int i = nn-1; i >=0; i--) {
      result = mbcore_->get_coords(&(face_nodes[i]), 1, coords);
      ErrorCheck_(result, "Problem getting coordinates of node");

      it->set(space_dim_, coords);
      ++it;
    }
  }
  else {
    for (int i = 0; i < nn; i++) {
      result = mbcore_->get_coords(&(face_nodes[i]), 1, coords);
      ErrorCheck_(result, "Problem getting tag data");

      it->set(space_dim_, coords);
      ++it;
    }
  }

  delete [] coords;
}
  

//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
moab::Tag Mesh_MOAB::build_set(
    const Teuchos::RCP<const AmanziGeometry::Region>& region,
    const Entity_kind kind) const
{
  int celldim_ = Mesh::manifold_dimension();
  int space_dim = Mesh::space_dimension();
  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm = Mesh::geometric_model();
  int one = 1;
  moab::Tag tag(0);

  // Modify region/set name by prefixing it with the type of entity requested
  std::string internal_name = internal_name_of_set(region, kind);

  // Create entity set based on the region defintion  
  switch (kind) {      
  case CELL:    // cellsets      
    if (region->type() == AmanziGeometry::BOX ||
        region->type() == AmanziGeometry::COLORFUNCTION) {

      mbcore_->tag_get_handle(internal_name.c_str(), 1, MB_TYPE_INTEGER, tag,
                             MB_TAG_CREAT|MB_TAG_SPARSE);
      
      int ncell = num_entities(CELL, Parallel_kind::ALL);              

      for (int icell = 0; icell < ncell; icell++) {
        if (region->inside(cell_centroid(icell))) {
          mbcore_->tag_set_data(tag, &(cell_id_to_handle[icell]), 1, &one);
        }
      }
    }
    else if (region->type() == AmanziGeometry::POINT) {
      AmanziGeometry::Point vpnt(space_dim);
      AmanziGeometry::Point rgnpnt(space_dim);

      mbcore_->tag_get_handle(internal_name.c_str(), 1, MB_TYPE_INTEGER, tag,
                             MB_TAG_CREAT|MB_TAG_SPARSE);
      
      rgnpnt = Teuchos::rcp_static_cast<const AmanziGeometry::RegionPoint>(region)->point();
        
      int nnode = num_entities(NODE, Parallel_kind::ALL);
      double mindist2 = 1.e+16;
      int minnode = -1;
        
      int inode;
      for (inode = 0; inode < nnode; inode++) {
                  
        node_get_coordinates(inode, &vpnt);                  
        double dist2 = (vpnt - rgnpnt) * (vpnt - rgnpnt);
 
        if (dist2 < mindist2) {
          mindist2 = dist2;
          minnode = inode;
          if (mindist2 <= 1.0e-14)
            break;
        }
      }

      Entity_ID_View cells, cells1;

      node_get_cells(minnode, Parallel_kind::ALL, &cells);
      int ncells = cells.size();

      for (int ic = 0; ic < ncells; ic++) {
        Entity_ID icell = cells[ic];
        
        // Check if point is contained in cell            
        if (point_in_cell(rgnpnt, icell)) {
          mbcore_->tag_set_data(tag, &(cell_id_to_handle[icell]), 1, &one);
        }
      }
    }

    else if (region->type() == AmanziGeometry::PLANE) {

      mbcore_->tag_get_handle(internal_name.c_str(), 1, MB_TYPE_INTEGER, tag,
                             MB_TAG_CREAT|MB_TAG_SPARSE);
      
      if (celldim_ == 2) {
        int ncells = num_entities(CELL, Parallel_kind::ALL);              

        for (int ic = 0; ic < ncells; ic++) {
          Point_List ccoords(space_dim);

          cell_get_coordinates(ic, &ccoords);

          bool on_plane = true;
          for (int j = 0; j < ccoords.size(); j++) {
            if (!region->inside(ccoords[j])) {
              on_plane = false;
              break;
            }
          }
                  
          if (on_plane)
            mbcore_->tag_set_data(tag, &(cell_id_to_handle[ic]), 1, &one);
        }
      }

    }
    else if (region->type() == AmanziGeometry::LOGICAL) {
      // will process later in this subroutine
    }
    else if (region->type() == AmanziGeometry::LABELEDSET) {
      // Nothing to do
      tag = cstag;
    }
    else {
      Errors::Message mesg;
      mesg << "Region \"" << region->name() << "\" type not applicable/supported for cell sets";
      amanzi_throw(mesg);
    }
      
    break;
      
  case FACE:  // sidesets
    if (region->type() == AmanziGeometry::BOX)  {

      mbcore_->tag_get_handle(internal_name.c_str(), 1, MB_TYPE_INTEGER, tag,
                             MB_TAG_CREAT|MB_TAG_SPARSE);

      int nface = num_entities(FACE, Parallel_kind::ALL);
        
      for (int iface = 0; iface < nface; iface++) {
        if (region->inside(face_centroid(iface)))
          mbcore_->tag_set_data(tag, &(face_id_to_handle[iface]), 1, &one);

      }
    }
    else if (region->type() == AmanziGeometry::PLANE ||
             region->type() == AmanziGeometry::POLYGON) {

      mbcore_->tag_get_handle(internal_name.c_str(), 1, MB_TYPE_INTEGER, tag,
                             MB_TAG_CREAT | MB_TAG_SPARSE);

      int nface = num_entities(FACE, Parallel_kind::ALL);
              
      for (int iface = 0; iface < nface; iface++) {
        Point_List fcoords(space_dim);
            
        face_get_coordinates(iface, &fcoords);
            
        bool on_plane = true;
        for (int j = 0; j < fcoords.size(); j++) {
          if (!region->inside(fcoords[j])) {
            on_plane = false;
            break;
          }
        }
                  
        if (on_plane)
          mbcore_->tag_set_data(tag, &(face_id_to_handle[iface]), 1, &one);
      }
    }
    else if (region->type() == AmanziGeometry::LABELEDSET) {
      // Nothing to do
      tag = sstag;
    }
    else if (region->type() == AmanziGeometry::LOGICAL) {
      // Will handle it later in the routine
    }
    else {
      Errors::Message mesg;
      mesg << "Region \"" << region->name() << "\" type not applicable/supported for face sets";
      amanzi_throw(mesg);
    }
    break;

  case NODE: // Nodesets
    if (region->type() == AmanziGeometry::BOX ||
        region->type() == AmanziGeometry::PLANE ||
        region->type() == AmanziGeometry::POLYGON ||
        region->type() == AmanziGeometry::POINT) {

      mbcore_->tag_get_handle(internal_name.c_str(), 1, MB_TYPE_INTEGER, tag,
                             MB_TAG_CREAT|MB_TAG_SPARSE);

      int nnode = num_entities(NODE, Parallel_kind::ALL);

      for (int inode = 0; inode < nnode; inode++) {

        AmanziGeometry::Point vpnt(space_dim);
        node_get_coordinates(inode, &vpnt);
                  
        if (region->inside(vpnt)) {
          mbcore_->tag_set_data(tag, &(node_id_to_handle[inode]), 1, &one);

          // Only one node per point region
          if (region->type() == AmanziGeometry::POINT)
            break;      
        }
      }
    }
    else if (region->type() == AmanziGeometry::LABELEDSET) {
      // Just retrieve and return the set
      tag = nstag;
    }
    else if (region->type() == AmanziGeometry::LOGICAL) {
      // We will handle it later in the routine
    }
    else {
      Errors::Message mesg("Region type not applicable/supported for node sets");
      amanzi_throw(mesg);
    }
      
    break;

  default:
    Errors::Message mesg("Unknown geometric entity");
    amanzi_throw(mesg);
    break;
  }

  if (region->type() == AmanziGeometry::LOGICAL) {
    std::string new_internal_name;

    Teuchos::RCP<const AmanziGeometry::RegionLogical> boolregion =
        Teuchos::rcp_static_cast<const AmanziGeometry::RegionLogical>(region);
    std::vector<std::string> region_names = boolregion->component_regions();
    int nreg = region_names.size();
    
    std::vector<moab::Tag> tags;
    std::vector<Teuchos::RCP<const AmanziGeometry::Region> > regions;
    moab::Tag tag1;
    moab::Range entset;
    
    for (int r = 0; r < nreg; r++) {
      Teuchos::RCP<const AmanziGeometry::Region> rgn1 = gm->FindRegion(region_names[r]);
      regions.push_back(rgn1);

      // Did not find the region
      if (rgn1 == Teuchos::null) {
        std::stringstream mesg_stream;
        mesg_stream << "Geometric model has no region named " << 
          region_names[r];
        Errors::Message mesg(mesg_stream.str());
        amanzi_throw(mesg);
      }
        
      internal_name = internal_name_of_set(rgn1, kind);
      mbcore_->tag_get_handle(internal_name.c_str(), 1, MB_TYPE_INTEGER,
                             tag1, MB_TAG_SPARSE);
      if (!tag1)        
        tag1 = build_set(rgn1, kind);  // Recursive call

      tags.push_back(tag1);
    }

    // Check the entity types of the sets are consistent with the
    // entity type of the requested set
    if (boolregion->operation() == AmanziGeometry::COMPLEMENT) {
      int *values[1] = {&one};
      moab::Range entset1, entset2;

      switch (kind) {
      case CELL:
        for (int i = 0; i < tags.size(); i++) {
          mbcore_->get_entities_by_type_and_tag(0, MBHEX, &(tags[i]),
                                               (void **)values, 1, entset2);
          
          entset1.merge(entset2);
        }
        entset = AllCells;
        entset -= entset1;
        break;
      case FACE:
        for (int i = 0; i < tags.size(); i++) {
          mbcore_->get_entities_by_type_and_tag(0, MBQUAD, &(tags[i]),
                                               (void **)values, 1, entset2);
          
          entset1.merge(entset2);
        }
        entset = AllFaces;
        entset -= entset1;
        break;
      case NODE:
        for (int i = 0; i < tags.size(); i++) {
          mbcore_->get_entities_by_type_and_tag(0, MBVERTEX, &(tags[i]),
                                               (void **)values, 1, entset2);
          
          entset1.merge(entset2);
        }
        entset = AllVerts;
        entset -= entset1;
        break;

      default:
        Errors::Message mesg("Unknown geometric entity");
        amanzi_throw(mesg);
        break;
      }

      for (int r = 0; r < nreg; r++)
        new_internal_name = new_internal_name + "+" + region_names[r];
      new_internal_name = "NOT_" + new_internal_name;
    }

    else if (boolregion->operation() == AmanziGeometry::UNION) {
      moab::Range entset1;
      switch (kind) {
      case CELL:
        for (int i = 0; i < tags.size(); i++) {
          int *values[1] = {&one};
          mbcore_->get_entities_by_type_and_tag(0, MBHEX, &(tags[i]),
                                               (void **)values, 1, entset1);

          entset.merge(entset1);
        }
        break;
      case FACE:
        for (int i = 0; i < tags.size(); i++) {
          int *values[1] = {&one};
          mbcore_->get_entities_by_type_and_tag(0, MBQUAD, &(tags[i]),
                                               (void **)values, 1, entset1);

          entset.merge(entset1);
        }
        break;
      case NODE:
        for (int i = 0; i < tags.size(); i++) {
          int *values[1] = {&one};
          mbcore_->get_entities_by_type_and_tag(0, MBVERTEX, &(tags[i]),
                                               (void **)values, 1, entset1);

          entset.merge(entset1);
        }
        break;

      default:
        Errors::Message mesg("Unknown geometric entity");
        amanzi_throw(mesg);
        break;
      }
      
      for (int r = 0; r < nreg; r++)
        new_internal_name = new_internal_name + "+" + region_names[r];
    }

    else if (boolregion->operation() == AmanziGeometry::SUBTRACT) {
      int *values[1] = { &one };
      moab::Range entset1, entset2;

      switch (kind) {
      case CELL:
        mbcore_->get_entities_by_type_and_tag(0, MBHEX, &(tags[0]),
                                             (void **)values, 1, entset);
        for (int i = 1; i < tags.size(); i++) {
          mbcore_->get_entities_by_type_and_tag(0, MBHEX, &(tags[i]),
                                               (void **)values, 1, entset1);
          
          entset.merge(entset1);
        }
        break;
      case FACE:
        mbcore_->get_entities_by_type_and_tag(0, MBQUAD, &(tags[0]),
                                             (void **)values, 1, entset);
        for (int i = 1; i < tags.size(); i++) {
          mbcore_->get_entities_by_type_and_tag(0, MBQUAD, &(tags[i]),
                                               (void **)values, 1, entset1);
          
          entset.merge(entset1);
        }
        break;
      case NODE:
        mbcore_->get_entities_by_type_and_tag(0, MBVERTEX, &(tags[0]),
                                             (void **)values, 1, entset);
        for (int i = 1; i < tags.size(); i++) {
          mbcore_->get_entities_by_type_and_tag(0, MBVERTEX, &(tags[i]),
                                               (void **)values, 1, entset1);
          
          entset.merge(entset1);
        }
        break;

      default:
        Errors::Message mesg("Unknown geometric entity");
        amanzi_throw(mesg);
        break;
      }
      
      new_internal_name = region_names[0];
      for (int r = 0; r < nreg; r++)
        new_internal_name = new_internal_name + "-" + region_names[r];
    }
    else if (boolregion->operation() == AmanziGeometry::INTERSECT) {
      Errors::Message mesg("INTERSECT region not implemented in MOAB");
      amanzi_throw(mesg);
    }

    mbcore_->tag_get_handle(new_internal_name.c_str(), 1, MB_TYPE_INTEGER,
                           tag, MB_TAG_CREAT | MB_TAG_SPARSE);

    for (auto it = entset.begin(); it != entset.end(); ++it) {
      moab::EntityHandle ent = *it;
      mbcore_->tag_set_data(tag, &ent, 1, &one);
    }
  }

  return tag;
}


//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
void Mesh_MOAB::get_set_entities_and_vofs(const std::string& setname, 
                                          const Entity_kind kind, 
                                          const Parallel_kind ptype,
                                          Entity_ID_View *setents,
                                          Double_List *vofs) const
{
  int lid, one=1;
  int space_dim = Mesh::space_dimension();

  AMANZI_ASSERT(setents != NULL);
  
  setents->clear();

  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm = Mesh::geometric_model();

  // Is there an appropriate region by this name?
  Teuchos::RCP<const AmanziGeometry::Region> rgn = gm->FindRegion(setname);

  moab::Range mset1, sets;

  // Did not find the region
  if (rgn == Teuchos::null) {
    std::stringstream mesg_stream;
    mesg_stream << "Geometric model has no region named " << setname;
    Errors::Message mesg(mesg_stream.str());
    amanzi_throw(mesg);
  }

  // If region is of type labeled set and a mesh set should have been
  // initialized from the input file
  if (rgn->type() == AmanziGeometry::LABELEDSET) {
    Teuchos::RCP<const AmanziGeometry::RegionLabeledSet> lsrgn =
        Teuchos::rcp_static_cast<const AmanziGeometry::RegionLabeledSet>(rgn);
    std::string label = lsrgn->label();
    std::stringstream labelstream(label);
    int labelint;
    labelstream >> labelint;
    std::string entity_type = lsrgn->entity_str();

    if ((kind == CELL && entity_type != "CELL") ||
        (kind == FACE && entity_type != "FACE") ||
        (kind == NODE && entity_type != "NODE")) {
      std::stringstream mesg_stream;
      mesg_stream << "Found labeled set region named " << setname 
                  << " but it contains entities of type " << entity_type << ", not the requested type";
      Errors::Message mesg(mesg_stream.str());
      amanzi_throw(mesg);
    } 

    Tag tag;
    if (kind == CELL) tag = cstag;
    else if (kind == FACE) tag = sstag;
    else if (kind == NODE) tag = nstag;

    mbcore_->get_entities_by_type_and_tag(0, MBENTITYSET, &tag, NULL, 1, sets);
    for (auto it = sets.begin(); it != sets.end(); ++it) {
      int set_id;
      mbcore_->tag_get_data(tag, &*it, 1, &set_id);
      if (labelint == set_id) {
        mbcore_->get_entities_by_handle(*it, mset1);
        break;
      }
    }
  }
  else {
    // Modify region/set name by prefixing it with the type of entity requested
    moab::Tag tag(0);
    int* values[1] = { &one };

    std::string internal_name = internal_name_of_set(rgn, kind);
    mbcore_->tag_get_handle(internal_name.c_str(), 1, MB_TYPE_INTEGER, tag, MB_TAG_SPARSE);
    if (!tag)
      tag = build_set(rgn, kind);

    if (kind == NODE) {
      mbcore_->get_entities_by_type_and_tag(0, MBVERTEX, &tag, (void **)values, 1, mset1);
    }
    else if (space_dim == 3 && kind == CELL) {
      mbcore_->get_entities_by_type_and_tag(0, MBHEX, &tag, (void **)values, 1, mset1);
    } else if (space_dim == 3 && kind == FACE) {
      mbcore_->get_entities_by_type_and_tag(0, MBQUAD, &tag, (void **)values, 1, mset1);
    }
    else if (space_dim == 2 && kind == CELL) {
      mbcore_->get_entities_by_type_and_tag(0, MBQUAD, &tag, (void **)values, 1, mset1);
    } else if (space_dim == 2 && kind == FACE) {
      mbcore_->get_entities_by_type_and_tag(0, MBEDGE, &tag, (void **)values, 1, mset1);
    }
  }

  // Check if no processor got any mesh entities
  int nent_loc = mset1.size();

#ifdef DEBUG
  int nent_glob;

  getComm()->SumAll(&nent_loc, &nent_glob, 1);
  if (nent_glob == 0) {
    std::stringstream mesg_stream;
    mesg_stream << "Could not retrieve any mesh entities for set " << setname << std::endl;
    Errors::Message mesg(mesg_stream.str());
    Exceptions::amanzi_throw(mesg);
  }
#endif
  
  setents->resize(nent_loc);
  if (nent_loc) {
    unsigned char pstatus;
    nent_loc = 0; // reset and count to get the real number

    switch (ptype) {
    case Parallel_kind::OWNED:
      for (auto it = mset1.begin(); it != mset1.end(); ++it) {
        moab::EntityHandle ent = *it;

        mbcomm_->get_pstatus(ent, pstatus);
        if ((pstatus & PSTATUS_NOT_OWNED) == 0) {
          mbcore_->tag_get_data(lid_tag, &ent, 1, &lid);
          (*setents)[nent_loc++] = lid;
        }
      }
      break;
    case Parallel_kind::GHOST:
      for (auto it = mset1.begin(); it != mset1.end(); ++it) {
        moab::EntityHandle ent = *it;

        mbcomm_->get_pstatus(ent, pstatus);
        if ((pstatus & PSTATUS_NOT_OWNED) == 1) {
          mbcore_->tag_get_data(lid_tag, &ent, 1, &lid);
          (*setents)[nent_loc++] = lid;
        }
      }
      break;
    case Parallel_kind::ALL:
      for (auto it = mset1.begin(); it != mset1.end(); ++it) {
        moab::EntityHandle ent = *it;

        mbcore_->tag_get_data(lid_tag, &ent, 1, &lid);
        (*setents)[nent_loc++] = lid;
      }
      break;

    default:
      Errors::Message mesg("Unknown geometric entity");
      amanzi_throw(mesg);
      break;
    }
    
    setents->resize(nent_loc);
  }
      
  // Check if there were no entities left on any processor after
  // extracting the appropriate category of entities 
#ifdef DEBUG
  getComm()->SumAll(&nent_loc, &nent_glob, 1);
  
  if (nent_glob == 0) {
    std::stringstream mesg_stream;
    mesg_stream << "Could not retrieve any mesh entities of type " << setkind 
                << " for set " << setname << std::endl;
    Errors::Message mesg(mesg_stream.str());
    Exceptions::amanzi_throw(mesg);
  }
#endif
}


//-------------------
// Upward adjacencies
//-------------------

//--------------------------------------------------------------------
// Cells of type 'ptype' connected to a node
//--------------------------------------------------------------------
void Mesh_MOAB::node_get_cells(const Entity_ID nodeid, 
                               const Parallel_kind ptype,
                               Entity_ID_View *cellids) const 
{
  unsigned char pstatus;
  std::vector<EntityHandle> ids;

  mbcore_->get_adjacencies(&(node_id_to_handle[nodeid]), 1, celldim, true, ids);
  int lid, nids = ids.size();

  cellids->clear();
  for (int i = 0; i < nids; ++i) {
    mbcore_->tag_get_data(lid_tag, &(ids[i]), 1, &lid);

    if (ptype == Parallel_kind::ALL) {
      cellids->push_back(lid);
    } if (ptype == Parallel_kind::OWNED) {
      mbcomm_->get_pstatus(ids[i], pstatus);
      if ((pstatus & PSTATUS_NOT_OWNED) == 0) {
        cellids->push_back(lid);
      }
    } else if (ptype == Parallel_kind::GHOST) {
      mbcomm_->get_pstatus(ids[i], pstatus);
      if ((pstatus & PSTATUS_NOT_OWNED) == 1) {
        cellids->push_back(lid);
      }
    }
  }
}
    

//--------------------------------------------------------------------
// Faces of type 'ptype' connected to a node
//--------------------------------------------------------------------
void Mesh_MOAB::node_get_faces(const Entity_ID nodeid, 
                               const Parallel_kind ptype,
                               Entity_ID_View *faceids) const
{
  throw std::exception();
}
    

//--------------------------------------------------------------------
// Cells connected to a face
//--------------------------------------------------------------------
void Mesh_MOAB::face_get_cells_internal_(const Entity_ID faceid, 
                                         const Parallel_kind ptype,
                                         Entity_ID_View *cellids) const
{
  int result;
  moab::EntityHandle face = face_id_to_handle[faceid];
  moab::Range fcells;
  
  result = mbcore_->get_adjacencies(&face, 1, celldim, true, fcells, Core::UNION);
  if (result != MB_SUCCESS) {
    std::cerr << "Could not get cells of face" << faceid << std::endl;
    AMANZI_ASSERT(result == MB_SUCCESS);
  }
  
  int nc = fcells.size();
  int fcellids[2];

  result = mbcore_->tag_get_data(lid_tag, fcells, (void *)fcellids);
  ErrorCheck_(result, "Problem getting id tag data");

  cellids->resize(2);
  Entity_ID_View::iterator it = cellids->begin();

  unsigned char pstatus;
  int n = 0;
  switch (ptype) {
  case Parallel_kind::ALL:
    for (int i = 0; i < nc; i++) {
      *it = fcellids[i];
      ++it;
      ++n;
    }
    break;
  case Parallel_kind::OWNED:
    for (int i = 0; i < nc; i++) {
      result = mbcomm_->get_pstatus(fcells[i], pstatus);
      if ((pstatus & PSTATUS_NOT_OWNED) == 0) {
        *it = fcellids[i];
        ++it;
        ++n;
      }
    }
    break;
  case Parallel_kind::GHOST:
    for (int i = 0; i < nc; i++) {
      result = mbcomm_->get_pstatus(fcells[i], pstatus);
      if ((pstatus & PSTATUS_NOT_OWNED) == 1) {
        *it = fcellids[i];
        ++it;
        ++n;
      }
    }
    break;
  default:
    Errors::Message mesg("Unknown geometric entity");
    amanzi_throw(mesg);
    break;
  }

  cellids->resize(n);
}
    

//-----------------------
// Same level adjacencies
//-----------------------

//--------------------------------------------------------------------
// Face connected neighboring cells of given cell of a particular ptype
// (e.g. a hex has 6 face neighbors)

// The order in which the cellids are returned cannot be
// guaranteed in general except when ptype = ALL, in which case
// the cellids will correcpond to cells across the respective
// faces given by cell_get_faces
//--------------------------------------------------------------------
void Mesh_MOAB::cell_get_face_adj_cells(const Entity_ID cellid,
                                        const Parallel_kind ptype,
                                        Entity_ID_View *fadj_cellids) const
{
  throw std::exception();
}


//--------------------------------------------------------------------
// Epetra map for cells - basically a structure specifying the
// global IDs of cells owned or used by this processor
//--------------------------------------------------------------------
void Mesh_MOAB::init_cell_map()
{
  int *cell_gids;
  int ncell, result;

  if (!serial_run) {
    // For parallel runs create map without and with ghost cells included
    // Also, put in owned cells before the ghost cells
    cell_gids = new int[OwnedCells.size() + GhostCells.size()];
    
    result = mbcore_->tag_get_data(gid_tag, OwnedCells, cell_gids);
    ErrorCheck_(result, "Problem getting tag data");
    ncell = OwnedCells.size();
    
    cell_map_wo_ghosts_ = new Epetra_Map(-1, ncell, cell_gids, 0, *getComm());
    
    result = mbcore_->tag_get_data(gid_tag, GhostCells, &(cell_gids[ncell]));
    ErrorCheck_(result, "Problem getting tag data");
    
    ncell += GhostCells.size();

    cell_map_w_ghosts_ = new Epetra_Map(-1, ncell, cell_gids, 0, *getComm());
  }
  else {
    cell_gids = new int[AllCells.size()];

    result = mbcore_->tag_get_data(gid_tag, AllCells, cell_gids);
    ErrorCheck_(result, "Problem getting tag data");

    ncell = AllCells.size();
    cell_map_wo_ghosts_ = new Epetra_Map(-1, ncell, cell_gids, 0, *getComm());
  }

  delete [] cell_gids;
}


//--------------------------------------------------------------------
// Epetra map for faces - basically a structure specifying the
// global IDs of cells owned or used by this processor
//--------------------------------------------------------------------
void Mesh_MOAB::init_face_map()
{
  int *face_gids, *extface_gids;
  int nface, result;

  if (!serial_run) {
    // For parallel runs create map without and with ghost cells included
    // Also, put in owned cells before the ghost cells
    int nowned = OwnedFaces.size();
    int nnotowned = NotOwnedFaces.size();
    face_gids = new int[nowned + nnotowned];
    extface_gids = new int[nowned + nnotowned];
    
    result = mbcore_->tag_get_data(gid_tag, OwnedFaces, face_gids);
    ErrorCheck_(result, "Problem getting tag data");

    nface = nowned;
    face_map_wo_ghosts_ = new Epetra_Map(-1, nface, face_gids, 0, *getComm());

    result = mbcore_->tag_get_data(gid_tag, NotOwnedFaces, &(face_gids[nface]));
    ErrorCheck_(result, "Problem getting tag data");
    
    nface += nnotowned;
    face_map_w_ghosts_ = new Epetra_Map(-1, nface, face_gids, 0, *getComm());

    // domain boundary faces: owned
    int n_extface_owned, n_extface = 0;
    for (int f = 0; f < nface; ++f) {
      moab::EntityHandle face = face_id_to_handle[f];
      moab::Range fcells;
  
      result = mbcore_->get_adjacencies(&face, 1, celldim, true, fcells, Core::UNION);
      ErrorCheck_(result, "Could not get cells of face");
  
      if (fcells.size() == 1)
        extface_gids[n_extface++] = face_gids[f];
      
      if (f == nowned - 1) n_extface_owned = n_extface;
    }

    extface_map_wo_ghosts_ = new Epetra_Map(-1, n_extface_owned, extface_gids, 0, *getComm());

    // domain boundary faces: owned + ghost. We filter out incorrect ghost faces
    int n_extface_notowned = n_extface - n_extface_owned;
    std::vector<int> pr_id(n_extface_notowned), lc_id(n_extface_notowned);
    extface_map_wo_ghosts_->RemoteIDList(n_extface_notowned, &(extface_gids[n_extface_owned]), pr_id.data(), lc_id.data());

    n_extface = n_extface_owned;
    for (int k = 0; k < n_extface_notowned; k++) {
      if (pr_id[k] >= 0)
        extface_gids[n_extface++] = extface_gids[n_extface_owned + k];
    }

    extface_map_w_ghosts_ = new Epetra_Map(-1, n_extface, extface_gids, 0, *getComm());
  }
  else {
    // all faces
    nface = AllFaces.size();
    face_gids = new int[nface];

    result = mbcore_->tag_get_data(gid_tag, AllFaces, face_gids);
    ErrorCheck_(result, "Problem getting tag data");

    face_map_wo_ghosts_ = new Epetra_Map(-1, nface, face_gids, 0, *getComm());

    // boundary faces
    int n_extface = 0;
    extface_gids = new int[nface];

    Entity_ID_View cellids;
    for (int f = 0; f < nface; ++f) {
      face_get_cells_internal_(f, Parallel_kind::ALL, &cellids);
      if (cellids.size() == 1)
        extface_gids[n_extface++] = face_gids[f];
    }

    extface_map_wo_ghosts_ = new Epetra_Map(-1, n_extface, extface_gids, 0, *getComm());
  }

  delete [] face_gids;
  delete [] extface_gids;
}


//--------------------------------------------------------------------
// Epetra map for nodes - basically a structure specifying the
// global IDs of cells owned or used by this processor
//--------------------------------------------------------------------
void Mesh_MOAB::init_node_map()
{
  int *vert_gids;
  int nvert, result;

  if (!serial_run) {
    // For parallel runs create map without and with ghost verts included
    // Also, put in owned cells before the ghost verts
    vert_gids = new int[OwnedVerts.size()+NotOwnedVerts.size()];
    
    result = mbcore_->tag_get_data(gid_tag, OwnedVerts, vert_gids);
    ErrorCheck_(result, "Problem getting tag data");

    nvert = OwnedVerts.size();
    
    node_map_wo_ghosts_ = new Epetra_Map(-1, nvert, vert_gids, 0, *getComm());
    
    result = mbcore_->tag_get_data(gid_tag, NotOwnedVerts, &(vert_gids[nvert]));
    ErrorCheck_(result, "Problem getting tag data");
    
    nvert += NotOwnedVerts.size();

    node_map_w_ghosts_ = new Epetra_Map(-1, nvert, vert_gids, 0, *getComm());
  }
  else {
    vert_gids = new int[AllVerts.size()];

    result = mbcore_->tag_get_data(gid_tag, AllVerts, vert_gids);
    ErrorCheck_(result, "Problem getting tag data");

    nvert = AllVerts.size();
    node_map_wo_ghosts_ = new Epetra_Map(-1, nvert, vert_gids, 0, *getComm());
  }

  delete [] vert_gids;
}


//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
Entity_ID Mesh_MOAB::GID(Entity_ID lid, Entity_kind kind) const
{
  moab::EntityHandle ent;
  Entity_ID gid;

  switch (kind) {
  case NODE:
    ent = node_id_to_handle[lid];
    break;

  case FACE:
    ent = face_id_to_handle[lid];
    break;

  case CELL:
    ent = cell_id_to_handle[lid];
    break;

  default:
    std::cerr << "Global ID requested for unknown entity type" << std::endl;
  }

  int result = mbcore_->tag_get_data(gid_tag, &ent, 1, &gid);
  ErrorCheck_(result, "Problem getting tag data");

  return gid;
}


//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
const Epetra_Map& Mesh_MOAB::cell_map(bool include_ghost) const
{
  if (serial_run)
    return *cell_map_wo_ghosts_;
  else
    return (include_ghost ? *cell_map_w_ghosts_ : *cell_map_wo_ghosts_);
}


//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
const Epetra_Map& Mesh_MOAB::face_map(bool include_ghost) const
{
  if (serial_run)
    return *face_map_wo_ghosts_;
  else
    return (include_ghost ? *face_map_w_ghosts_ : *face_map_wo_ghosts_);
}


//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
const Epetra_Map& Mesh_MOAB::node_map(bool include_ghost) const
{
  if (serial_run)
    return *node_map_wo_ghosts_;
  else
    return (include_ghost ? *node_map_w_ghosts_ : *node_map_wo_ghosts_);
}


//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
const Epetra_Map& Mesh_MOAB::exterior_face_map(bool include_ghost) const {
  if (serial_run)
    return *extface_map_wo_ghosts_;
  else
    return (include_ghost ? *extface_map_w_ghosts_ : *extface_map_wo_ghosts_);
}


//--------------------------------------------------------------------
// Epetra importer that will allow apps to import values from a Epetra
// vector defined on all owned faces into an Epetra vector defined
// only on exterior faces
//--------------------------------------------------------------------
const Epetra_Import& Mesh_MOAB::exterior_face_importer(void) const
{
  Errors::Message mesg("Exterior face importer is not implemented");
  amanzi_throw(mesg);
}


//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
const Epetra_Map& Mesh_MOAB::exterior_node_map(bool include_ghost) const {
  Errors::Message mesg("Exterior node map is not implemented");
  Exceptions::amanzi_throw(mesg);
  throw(mesg);
}


//--------------------------------------------------------------------
// Get parallel type of enetity
//--------------------------------------------------------------------
Parallel_kind Mesh_MOAB::entity_get_ptype(const Entity_kind kind, 
                                          const Entity_ID entid) const
{
  moab::EntityHandle ent(0);
  unsigned char pstatus;

  switch (kind) {
  case NODE:
    ent = node_id_to_handle[entid];
    break;
      
  case FACE:
    ent = face_id_to_handle[entid];
    break;
      
  case CELL:
    ent = cell_id_to_handle[entid];
    break;

  default:
    AMANZI_ASSERT("Global ID requested for unknown entity type");
  }

  mbcomm_->get_pstatus(ent, pstatus);
  return ((pstatus & PSTATUS_NOT_OWNED) == 1) ? Parallel_kind::GHOST : Parallel_kind::OWNED;
}
  
  
//--------------------------------------------------------------------
// Get cell type
//--------------------------------------------------------------------
Cell_kind Mesh_MOAB::cell_get_type(const Entity_ID cellid) const
{
  if (space_dim_ == 2) return QUAD;
  return HEX;
}
    

//--------------------------------------------------------------------
// TBW
//--------------------------------------------------------------------
std::string 
Mesh_MOAB::internal_name_of_set(const Teuchos::RCP<const AmanziGeometry::Region>& r,
                                const Entity_kind entity_kind) const
{
  std::string internal_name;
  
  if (r->type() == AmanziGeometry::LABELEDSET) {
    
    Teuchos::RCP<const AmanziGeometry::RegionLabeledSet> lsrgn =
        Teuchos::rcp_static_cast<const AmanziGeometry::RegionLabeledSet>(r);
    std::string label = lsrgn->label();

    if (entity_kind == CELL)
      internal_name = "matset_" + label;
    else if (entity_kind == FACE)
      internal_name = "sideset_" + label;
    else if (entity_kind == NODE)
      internal_name = "nodeset_" + label;      
  }
  else {
    if (entity_kind == CELL)
      internal_name = "CELLSET_" + r->name();
    else if (entity_kind == FACE)
      internal_name = "FACESET_" + r->name();
    else if (entity_kind == NODE)
      internal_name = "NODESET_" + r->name();
  }

  return internal_name;
}


//--------------------------------------------------------------------
// Deform a mesh so that cell volumes conform as closely as possible
// to target volumes without dropping below the minimum volumes.  If
// move_vertical = true, nodes will be allowed to move only in the
// vertical direction (right now arbitrary node movement is not allowed)
//--------------------------------------------------------------------
int Mesh_MOAB::deform(const Double_List& target_cell_volumes_in, 
                      const Double_List& min_cell_volumes_in, 
                      const Entity_ID_View& fixed_nodes,
                      const bool move_vertical) 
{
  Errors::Message mesg("Deformation not implemented for Mesh_MOAB");
  amanzi_throw(mesg);
  return 1;
}


//--------------------------------------------------------------------
// Miscellaneous
//--------------------------------------------------------------------
void Mesh_MOAB::write_to_exodus_file(const std::string filename) const {
  throw std::exception();
}

}  // namespace AmanziMesh
}  // namespace Amanzi

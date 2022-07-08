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
#include "RegionEnumerated.hh"

#include "Mesh_MSTK.hh"

using namespace std;

namespace Amanzi {

namespace AmanziMesh {

char kind_to_string[4][256] = {"NODE","EDGE","FACE","CELL"};

void Mesh_MSTK::init_mesh_from_file_(const std::string& filename,
                                     const Partitioner_type partitioner) {

  int ok = 0;

  mesh_ = MESH_New(F1);

  if (filename.find(".exo") != std::string::npos) {  // Exodus file

    // Read the mesh on processor 0
    ok = MESH_ImportFromExodusII(mesh_, filename.c_str(), NULL, mpicomm_);

    // Collapse any degenerate edges in the mesh
    collapse_degen_edges();  // Assumes its operating on member var 'mesh'

    // Renumber local IDs to be contiguous
    MESH_Renumber(mesh_, 0, MALLTYPE);

    if (numprocs > 1) {
      // Distribute the mesh to all the processors
      int topo_dim = MESH_Num_Regions(mesh_) ? 3 : 2;
      int num_ghost_layers = 1;
      int with_attr = 1;  // Redistribute any attributes and sets
      int method = static_cast<int>(partitioner);
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
    collapse_degen_edges();

    // Renumber local IDs to be contiguous
    MESH_Renumber(mesh_, 0, MALLTYPE);

    // Weave the meshes together to form interprocessor connections
    int num_ghost_layers = 1;
    int input_type = 1;  // We are given partitioned meshes with a
    //                   // unique global ID on each mesh vertex
    int topo_dim = MESH_Num_Regions(mesh_) ? 3 : 2;

    ok &= MSTK_Weave_DistributedMeshes(mesh_, topo_dim, num_ghost_layers,
                                       input_type, mpicomm_);

    // Global IDs are discontinuous due to elimination of degeneracies
    // but we cannot renumber global IDs until interprocessor
    // connectivity has been established

    if (contiguous_gids_) {
      ok &= MESH_Renumber_GlobalIDs(mesh_, MALLTYPE, 0, NULL, mpicomm_);
    }

  } else {
    std::stringstream mesg_stream;
    mesg_stream << "Cannot identify file type from extension of input file " <<
        filename << " on processor " << myprocid << std::endl;
    Errors::Message mesg(mesg_stream.str());
    Exceptions::amanzi_throw(mesg);
  }

  if (!ok) {
    std::stringstream mesg_stream;
    mesg_stream << "Failed to load " << filename << " on processor " <<
        myprocid << std::endl;
    Errors::Message mesg(mesg_stream.str());
    Exceptions::amanzi_throw(mesg);
  }
}

//--------------------------------------
// Constructor - load up mesh from file
//--------------------------------------
Mesh_MSTK::Mesh_MSTK(const std::string& filename,
                     const Comm_ptr_type& comm,
                     const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
                     const Teuchos::RCP<const Teuchos::ParameterList>& plist,
                     const bool request_faces,
                     const bool request_edges) :
    Mesh(comm, gm, plist, request_faces, request_edges),
    faces_initialized(false), edges_initialized(false),
    extface_map_w_ghosts_(nullptr), extface_map_wo_ghosts_(nullptr),
    extnode_map_w_ghosts_(nullptr), extnode_map_wo_ghosts_(nullptr),
    owned_to_extface_importer_(nullptr),
    meshxyz(nullptr), 
    target_cell_volumes_(nullptr), min_cell_volumes_(nullptr)
{
  // extract optional control parameters, but first specify defaults
  contiguous_gids_ = true;
  Partitioner_type partitioner = PARTITIONER_DEFAULT;
    
  if (plist != Teuchos::null) {
    if (plist->isSublist("unstructured")) {
      const auto tmp = Teuchos::sublist(plist, "unstructured");
      if (tmp->isSublist("expert")) {
        const auto expert_list = Teuchos::sublist(tmp, "expert");

        // -- partitioner
        if (expert_list->isParameter("partitioner")) {
          std::string partitioner_str = expert_list->get<std::string>("partitioner");
          if (partitioner_str == "METIS" || partitioner_str == "metis")
            partitioner = Partitioner_type::METIS;
          else if (partitioner_str == "ZOLTAN_GRAPH" || partitioner_str == "zoltan_graph")
            partitioner = Partitioner_type::ZOLTAN_GRAPH;
          else if (partitioner_str == "ZOLTAN_RCB" || partitioner_str == "zoltan_rcb")
            partitioner = Partitioner_type::ZOLTAN_RCB;
        }
        if (expert_list->isParameter("contiguous global ids")) {
          contiguous_gids_ = expert_list->get<bool>("contiguous global ids");
        }
        if (expert_list->isParameter("request edges")) {
          edges_requested_ = expert_list->get<bool>("request edges");
        }
      }
    }
  }

  // Assume three dimensional problem if constructor called without 
  // the space_dimension parameter
  if (vo_.get() && vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *(vo_->os()) << "Construct mesh from file" << std::endl;
  }

  // Pre-processing (init, MPI queries etc)
  int space_dim = 3;
  pre_create_steps_(space_dim);

  init_mesh_from_file_(filename, partitioner);

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
    
    if (planar)
      space_dim = 2;

    comm->MaxAll(&space_dim,&max,1);

    space_dim = max;
    set_space_dimension(space_dim);      
  }

  // Verify mesh and geometric model compatibility

  if (gm != Teuchos::null && gm->dimension() != space_dimension()) {
    Exceptions::amanzi_throw(Errors::Message("Geometric model and mesh have different dimensions."));
  }

  // Do all the processing required for setting up the mesh for Amanzi 
  
  post_create_steps_(request_faces, edges_requested_);
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
                     const Teuchos::RCP<const Teuchos::ParameterList>& plist,
                     const bool request_faces,
                     const bool request_edges) :
    Mesh(comm, gm, plist, request_faces,request_edges), 
    faces_initialized(false), edges_initialized(false),
    extface_map_w_ghosts_(nullptr), extface_map_wo_ghosts_(nullptr),
    extnode_map_w_ghosts_(nullptr), extnode_map_wo_ghosts_(nullptr),
    owned_to_extface_importer_(nullptr),
    meshxyz(nullptr), 
    target_cell_volumes_(nullptr), min_cell_volumes_(nullptr)
{
  // extract optional control parameters, but first specify defaults
  contiguous_gids_ = true;
  Partitioner_type partitioner = PARTITIONER_DEFAULT;
    
  if (plist != Teuchos::null) {
    if (plist->isSublist("unstructured")) {
      const auto tmp = Teuchos::sublist(plist, "unstructured");
      if (tmp->isSublist("expert")) {
        const auto expert_list = Teuchos::sublist(tmp, "expert");

        // -- partitioner
        if (expert_list->isParameter("partitioner")) {
          std::string partitioner_str = expert_list->get<std::string>("partitioner");
          if (partitioner_str == "METIS" || partitioner_str == "metis")
            partitioner = Partitioner_type::METIS;
          else if (partitioner_str == "ZOLTAN_GRAPH" || partitioner_str == "zoltan_graph")
            partitioner = Partitioner_type::ZOLTAN_GRAPH;
          else if (partitioner_str == "ZOLTAN_RCB" || partitioner_str == "zoltan_rcb")
            partitioner = Partitioner_type::ZOLTAN_RCB;
        }
        if (expert_list->isParameter("contiguous global ids")) {
          contiguous_gids_ = expert_list->get<bool>("contiguous global ids");
        }
        if (expert_list->isParameter("request edges")) {
          edges_requested_ = expert_list->get<bool>("request edges");
        }
      }
    }
  }

  int ok;
  int space_dimension = 3;
  pre_create_steps_(space_dimension);

  set_mesh_type(RECTANGULAR);   // Discretizations can use this info if they want

  if (serial_run) {
    // Load serial mesh
    mesh_ = MESH_New(F1);
    ok = generate_regular_mesh(mesh_,x0,y0,z0,x1,y1,z1,nx,ny,nz);

    set_manifold_dimension(3);

    myprocid = 0;
  } else {
    Mesh_ptr globalmesh = nullptr;
    int topo_dim=3; // What is the topological dimension of the mesh
    int ring = 1; // One layer of ghost cells in parallel meshes
    int with_attr = 1;  // update of attributes in parallel meshes
    int del_inmesh = 1; // delete input mesh as soon as possible
    int method = static_cast<int>(partitioner);

    
    if (myprocid == 0) {
      globalmesh = MESH_New(F1);

      ok = generate_regular_mesh(globalmesh,x0,y0,z0,x1,y1,z1,nx,ny,nz);
      
      topo_dim = (MESH_Num_Regions(globalmesh) == 0) ? 2 : 3;
    }
    else {
      globalmesh = nullptr;
      ok = 1;
    }

#ifdef MSTK_2_21rc1_OR_NEWER
    ok = ok & MSTK_Mesh_Distribute(globalmesh,&mesh_,&topo_dim,ring,with_attr,
                                   method,del_inmesh,mpicomm_);
#else
    ok = ok & MSTK_Mesh_Distribute(globalmesh,&mesh_,&topo_dim,ring,with_attr,
                                   method,mpicomm_);
    if (myprocid == 0)
      MESH_Delete(globalmesh);
#endif

    set_manifold_dimension(topo_dim);
  }

  if (!ok) {
    std::stringstream mesg_stream;
    mesg_stream << "Failed to generate mesh on processor " << myprocid;
    Errors::Message mesg(mesg_stream.str());
    Exceptions::amanzi_throw(mesg);
    AMANZI_ASSERT(ok);
  }

  // Do all the processing required for setting up the mesh for Amanzi 
  
  post_create_steps_(request_faces, edges_requested_);
}


//--------------------------------------
// Construct a 2D regular quadrilateral mesh internally
//--------------------------------------
Mesh_MSTK::Mesh_MSTK(const double x0, const double y0,
                     const double x1, const double y1,
                     const int nx, const int ny, 
                     const Comm_ptr_type& comm,
                     const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
                     const Teuchos::RCP<const Teuchos::ParameterList>& plist,
                     const bool request_faces,
                     const bool request_edges) :
    Mesh(comm, gm, plist, request_faces, request_edges), 
    faces_initialized(false), edges_initialized(false),
    extface_map_w_ghosts_(nullptr), extface_map_wo_ghosts_(nullptr),
    extnode_map_w_ghosts_(nullptr), extnode_map_wo_ghosts_(nullptr),
    owned_to_extface_importer_(nullptr),
    meshxyz(nullptr), 
    target_cell_volumes_(nullptr), min_cell_volumes_(nullptr)
{
  // extract optional control parameters, but first specify defaults
  contiguous_gids_ = true;
  Partitioner_type partitioner = PARTITIONER_DEFAULT;
    
  if (plist != Teuchos::null) {
    if (plist->isSublist("unstructured")) {
      const auto tmp = Teuchos::sublist(plist, "unstructured");
      if (tmp->isSublist("expert")) {
        const auto expert_list = Teuchos::sublist(tmp, "expert");

        // -- partitioner
        if (expert_list->isParameter("partitioner")) {
          std::string partitioner_str = expert_list->get<std::string>("partitioner");
          if (partitioner_str == "METIS" || partitioner_str == "metis")
            partitioner = Partitioner_type::METIS;
          else if (partitioner_str == "ZOLTAN_GRAPH" || partitioner_str == "zoltan_graph")
            partitioner = Partitioner_type::ZOLTAN_GRAPH;
          else if (partitioner_str == "ZOLTAN_RCB" || partitioner_str == "zoltan_rcb")
            partitioner = Partitioner_type::ZOLTAN_RCB;
        }
        if (expert_list->isParameter("contiguous global ids")) {
          contiguous_gids_ = expert_list->get<bool>("contiguous global ids");
        }
      }
    }
  }

  int ok;
  int space_dim = 2;
  pre_create_steps_(space_dim);

  if (vo_.get() && vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << "Construct mesh from low/hi coords - 2D" << std::endl;
  }

  set_mesh_type(RECTANGULAR);   // Discretizations can use this info if they want

  int topo_dim=space_dim; // What is the topological dimension of the mesh
  set_manifold_dimension(topo_dim);

  if (serial_run) {

    // Load serial mesh

    mesh_ = MESH_New(F1);
    ok = generate_regular_mesh(mesh_,x0,y0,x1,y1,nx,ny);

    myprocid = 0;
  }
  else {
    Mesh_ptr globalmesh = nullptr;
    int ring = 1; // One layer of ghost cells in parallel meshes
    int with_attr = 1;  // update of attributes in parallel meshes
    int del_inmesh = 1; // delete input mesh at the earliest
    int method = static_cast<int>(partitioner);

    if (myprocid == 0) {
      globalmesh = MESH_New(F1);

      ok = generate_regular_mesh(globalmesh,x0,y0,x1,y1,nx,ny);
      
      topo_dim = (MESH_Num_Regions(globalmesh) == 0) ? 2 : 3;
    }
    else {
      globalmesh = nullptr;
      ok = 1;
    }

#ifdef MSTK_2_21rc1_OR_NEWER
    ok = ok & MSTK_Mesh_Distribute(globalmesh,&mesh_,&topo_dim,ring,with_attr,
                                   method,del_inmesh,mpicomm_);
#else
    ok = ok & MSTK_Mesh_Distribute(globalmesh,&mesh_,&topo_dim,ring,with_attr,
                                   method,mpicomm_);
    if (myprocid == 0)
      MESH_Delete(globalmesh);
#endif

  }

  if (!ok) {
    std::stringstream mesg_stream;
    mesg_stream << "Failed to generate mesh on processor " << myprocid;
    Errors::Message mesg(mesg_stream.str());
    Exceptions::amanzi_throw(mesg);
    AMANZI_ASSERT(ok);
  }

  // Do all the processing required for setting up the mesh for Amanzi 
  post_create_steps_(request_faces, request_edges);
}


//---------------------------------------------------------
// Extract MSTK entities from an ID list and make a new MSTK mesh
//---------------------------------------------------------
Mesh_MSTK::Mesh_MSTK(const Teuchos::RCP<const Mesh>& parent_mesh,
                     const Entity_ID_List& entity_ids, 
                     const Entity_kind entity_kind,
                     const bool flatten,
                     const Comm_ptr_type& comm,
                     const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
                     const Teuchos::RCP<const Teuchos::ParameterList>& plist,
                     const bool request_faces,
                     const bool request_edges) :
    Mesh(comm, gm == Teuchos::null ? parent_mesh->geometric_model() : gm,
         plist == Teuchos::null ? parent_mesh->parameter_list() : plist,
         request_faces, request_edges),
    extface_map_w_ghosts_(nullptr), extface_map_wo_ghosts_(nullptr),
    extnode_map_w_ghosts_(nullptr), extnode_map_wo_ghosts_(nullptr),
    owned_to_extface_importer_(nullptr),
    parent_mesh_(Teuchos::rcp_dynamic_cast<const Mesh_MSTK>(parent_mesh))
{  
  if (!parent_mesh_.get()) {
    Errors::Message mesg("Cannot extract an MSTK mesh from a non-MSTK mesh.");
    Exceptions::amanzi_throw(mesg);
  }
    
  auto parent_mesh_mstk = parent_mesh_->mesh_;

  // store pointers to the MESH_XXXFromID functions so that they can
  // be called without a switch statement 
  static MEntity_ptr (*MEntFromID[4])(Mesh_ptr,int) =
    {MESH_VertexFromID, MESH_EdgeFromID, MESH_FaceFromID, MESH_RegionFromID};

  MType entity_dim = parent_mesh_->entity_kind_to_mtype(entity_kind);

  // Also make sure that the mesh object can do fast queries on local IDs
  //
  // Commented out for now to avoid another backwards incompatible
  // change and if the mesh has no modifications it is fast. Also, its
  // a one time setup cost but it would be good to enable it sometime
  //
  // MESH_Enable_LocalIDSearch(parent_mesh_mstk);

  int nent = entity_ids.size();
  List_ptr src_ents = List_New(nent);
  for (int i = 0; i < nent; ++i) {
    MEntity_ptr ent = MEntFromID[entity_dim](parent_mesh_mstk,entity_ids[i]+1);
    List_Add(src_ents,ent);
  }
  
  extract_mstk_mesh(src_ents, entity_dim, flatten, request_faces, request_edges);

  List_Delete(src_ents);
}


//---------------------------------------------------------
// Translate a setname into a special string with decorations
// indicating which type of entity is in that set
//---------------------------------------------------------
std::string 
Mesh_MSTK::internal_name_of_set(const Teuchos::RCP<const AmanziGeometry::Region>& r,
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
    else if (entity_kind == EDGE)
      internal_name = "edgeset_not_supported";
    else if (entity_kind == NODE)
      internal_name = "nodeset_" + label;      
  }
  else {
    if (entity_kind == CELL)
      internal_name = "CELLSET_" + r->name();
    else if (entity_kind == FACE)
      internal_name = "FACESET_" + r->name();
    else if (entity_kind == EDGE)
      internal_name = "EDGESET_" + r->name();
    else if (entity_kind == NODE)
      internal_name = "NODESET_" + r->name();
  }

  return internal_name;
}


//---------------------------------------------------------
// Get an alternate name (elemset_N instead of matset_N) for sets of type 
// Labeled Set and entity kind Cell. For everything else return regular name 
//---------------------------------------------------------
std::string 
Mesh_MSTK::other_internal_name_of_set(const Teuchos::RCP<const AmanziGeometry::Region>& r,
                                      const Entity_kind entity_kind) const
{

  std::string internal_name;
  
  if (r->type() == AmanziGeometry::LABELEDSET && entity_kind == CELL) {
    
    Teuchos::RCP<const AmanziGeometry::RegionLabeledSet> lsrgn =
        Teuchos::rcp_static_cast<const AmanziGeometry::RegionLabeledSet>(r);
    std::string label = lsrgn->label();

    internal_name = "elemset_" + label;
    return internal_name;
  }
  else
    return internal_name_of_set(r,entity_kind);
}


//---------------------------------------------------------
// Extract a list of MSTK entities and make a new MSTK mesh
// For private use of Mesh_MSTK class only
//---------------------------------------------------------
void Mesh_MSTK::extract_mstk_mesh(List_ptr src_entities, 
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
    pre_create_steps_(parent_mesh_->space_dimension()-1);
  else
    pre_create_steps_(parent_mesh_->space_dimension());

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
  vparentatt = MAttrib_New(mesh_,"vparentatt",POINTER,MVERTEX);
  eparentatt = MAttrib_New(mesh_,"eparentatt",POINTER,MEDGE);
  fparentatt = MAttrib_New(mesh_,"fparentatt",POINTER,MFACE);
  rparentatt = MAttrib_New(mesh_,"rparentatt",POINTER,MREGION);
  
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
                MEnt_Set_AttVal(fverts_new[j],vparentatt,0,0.0,mv);
              }
            }
            List_Delete(fverts);
            
            rfaces_new[i] = MF_New(mesh_);
            MF_Set_Vertices(rfaces_new[i],nfv,fverts_new);
            MF_Set_GEntDim(rfaces_new[i],MF_GEntDim(mf));
            MF_Set_GEntID(rfaces_new[i],MF_GEntID(mf));
            rfdirs_new[i] = MR_FaceDir_i(mr,i);

            MEnt_Set_AttVal(mf,copyatt,ival,rval,rfaces_new[i]);
            MEnt_Set_AttVal(rfaces_new[i],fparentatt,0,0.0,mf);
          }
        }
        List_Delete(rfaces);

        MRegion_ptr mr_new = MR_New(mesh_);
        MR_Set_Faces(mr_new,nrf,rfaces_new,rfdirs_new);
        MR_Set_GEntID(mr_new,MR_GEntID(mr));

        MEnt_Set_AttVal(mr,copyatt,ival,rval,mr_new);
        MEnt_Set_AttVal(mr_new,rparentatt,0,0.0,mr);
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
                MEnt_Set_AttVal(mv_new,vparentatt,0,0.0,mv);
              }

              ME_Set_Vertex(fedges_new[j],k,mv_new);
              ME_Set_GEntDim(fedges_new[j],ME_GEntDim(me));
              ME_Set_GEntID(fedges_new[j],ME_GEntID(me));
              MEnt_Set_AttVal(me,copyatt,ival,rval,fedges_new[j]);
              MEnt_Set_AttVal(fedges_new[j],eparentatt,0,0.0,me);
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
        MEnt_Set_AttVal(mf_new,fparentatt,0,0.0,mf);
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
            MEnt_Set_AttVal(mv_new,vparentatt,0,0.0,mv);
          }

          ME_Set_Vertex(me_new,j,mv_new);
        }

        if (ME_GEntDim(me) == 1)
          ME_Set_GEntDim(me_new, 1);        
        MEnt_Set_AttVal(me,copyatt,ival,rval,me_new);
        MEnt_Set_AttVal(me_new,eparentatt,0,0.0,me);
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
        MEnt_Set_AttVal(mv_new,vparentatt,0,0.0,mv);
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
    int status = MSTK_Weave_DistributedMeshes(mesh_, manifold_dimension(),
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
        MEnt_Get_AttVal(mv,vparentatt,&ival,&rval,&pval);
        MEnt_Set_AttVal(mv,vparentgid_att,MV_GlobalID((MVertex_ptr)pval),0.0,
                        NULL);
      }
    MEdge_ptr me = nullptr;
    if (entity_dim != MREGION) { // edge parents not set on 3D extraction -- maybe they should be --etc
      idx = 0;
      while ((me = (MEdge_ptr) MESH_Next_Edge(mesh_,&idx)))
        if (ME_PType(me) == POVERLAP) {
          MEnt_Get_AttVal(me,eparentatt,&ival,&rval,&pval);
          MEnt_Set_AttVal(me,eparentgid_att,ME_GlobalID((MEdge_ptr)pval),0.0,
                          NULL);
        }
    }
    idx = 0;
    MFace_ptr mf = nullptr;
    while ((mf = (MFace_ptr) MESH_Next_Face(mesh_,&idx)))
      if (MF_PType(mf) == POVERLAP) {
        MEnt_Get_AttVal(mf,fparentatt,&ival,&rval,&pval);
        MEnt_Set_AttVal(mf,fparentgid_att,MF_GlobalID((MFace_ptr)pval),0.0,
                        NULL);
      }
    idx = 0;
    MRegion_ptr mr = nullptr;
    while ((mr = (MRegion_ptr) MESH_Next_Region(mesh_,&idx)))
      if (MR_PType(mr) == POVERLAP) {
        MEnt_Get_AttVal(mr,rparentatt,&ival,&rval,&pval);
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
      MEnt_Set_AttVal(mv,vparentatt,0,0.0,mv_parent);
    }
    if (entity_dim != MREGION) {
      idx = 0;
      while ((me = (MEdge_ptr) MESH_Next_GhostEdge(mesh_,&idx))) {
        MEnt_Get_AttVal(me,eparentgid_att,&ival,&rval,&pval);
        MEdge_ptr me_parent = MESH_EdgeFromGlobalID(parent_mesh_mstk,ival);
        if (!me_parent) {
          Errors::Message 
            mesg("Cannot find ghost edge with given global ID");
          Exceptions::amanzi_throw(mesg);
        }
        MEnt_Set_AttVal(me,eparentatt,0,0.0,me_parent);
      }
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
      MEnt_Set_AttVal(mf,fparentatt,0,0.0,mf_parent);
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
      MEnt_Set_AttVal(mr,rparentatt,0,0.0,mr_parent);
    }

    MAttrib_Delete(vparentgid_att);
    MAttrib_Delete(eparentgid_att);
    MAttrib_Delete(fparentgid_att);
    MAttrib_Delete(rparentgid_att);
  }


  // We have to do an extra step to build new labeled sets based on
  // labeled sets of the base mesh

  inherit_labeled_sets(copyatt, src_entities);


  // Do all the processing required for setting up the mesh for Amanzi 
  
  post_create_steps_(request_faces, request_edges);


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
// Destructor with cleanup
//---------------------------------------------------------
Mesh_MSTK::~Mesh_MSTK() {
  delete cell_map_wo_ghosts_;
  delete cell_map_w_ghosts_;
  if (face_map_wo_ghosts_) delete face_map_wo_ghosts_;
  if (face_map_w_ghosts_) delete face_map_w_ghosts_;
  if (edge_map_wo_ghosts_) delete edge_map_wo_ghosts_;
  if (edge_map_w_ghosts_) delete edge_map_w_ghosts_;
  delete node_map_wo_ghosts_;
  delete node_map_w_ghosts_;
  if (extface_map_wo_ghosts_) delete extface_map_wo_ghosts_;
  if (extface_map_w_ghosts_) delete extface_map_w_ghosts_;
  if (owned_to_extface_importer_) delete owned_to_extface_importer_;
  if (extnode_map_wo_ghosts_) delete extnode_map_wo_ghosts_;
  if (extnode_map_w_ghosts_) delete extnode_map_w_ghosts_;
  delete [] faceflip;
  if (edgeflip) delete [] edgeflip;

  if (OwnedVerts) MSet_Delete(OwnedVerts);
  if (NotOwnedVerts) MSet_Delete(NotOwnedVerts);
  if (OwnedEdges) MSet_Delete(OwnedEdges);
  if (NotOwnedEdges) MSet_Delete(NotOwnedEdges);
  if (OwnedFaces) MSet_Delete(OwnedFaces);
  if (NotOwnedFaces) MSet_Delete(NotOwnedFaces);
  if (OwnedCells) MSet_Delete(OwnedCells);
  if (GhostCells) MSet_Delete(GhostCells);

  if (entities_deleted) {
    if (deleted_vertices) List_Delete(deleted_vertices);
    if (deleted_edges) List_Delete(deleted_edges);
    if (deleted_faces) List_Delete(deleted_faces);
    if (deleted_regions) List_Delete(deleted_regions);
  }

  MAttrib_Delete(celltype_att);
  if (vparentatt) MAttrib_Delete(vparentatt);
  if (eparentatt) MAttrib_Delete(eparentatt);
  if (fparentatt) MAttrib_Delete(fparentatt);
  if (rparentatt) MAttrib_Delete(rparentatt);

  MESH_Delete(mesh_);
}


//---------------------------------------------------------
// Number of OWNED, GHOST or ALL entities of different types
//
// Number of entities of any kind (cell, face, edge, node) and in a
// particular category (OWNED, GHOST, ALL)
//---------------------------------------------------------
unsigned int Mesh_MSTK::num_entities(const Entity_kind kind, 
                                     const Parallel_type ptype) const
{
  switch (kind) {
  case NODE:

    switch (ptype) {
    case Parallel_type::OWNED:
      return MSet_Num_Entries(OwnedVerts);
      break;
    case Parallel_type::GHOST:
      return !serial_run ? MSet_Num_Entries(NotOwnedVerts) : 0;
      break;
    case Parallel_type::ALL:
      return MESH_Num_Vertices(mesh_);
      break;
    default:
      return 0;
    }
    break;

  case EDGE:

    AMANZI_ASSERT(edges_initialized);

    switch (ptype) {
    case Parallel_type::OWNED:
      return MSet_Num_Entries(OwnedEdges);
      break;
    case Parallel_type::GHOST:
      return !serial_run ? MSet_Num_Entries(NotOwnedEdges) : 0;
      break;
    case Parallel_type::ALL:
      return MESH_Num_Edges(mesh_);
      break;
    default:
      return 0;
    }
    break;


  case FACE:

    AMANZI_ASSERT(faces_initialized);

    switch (ptype) {
    case Parallel_type::OWNED:
      return MSet_Num_Entries(OwnedFaces);
      break;
    case Parallel_type::GHOST:
      return !serial_run ? MSet_Num_Entries(NotOwnedFaces) : 0;
      break;
    case Parallel_type::ALL:
      return (manifold_dimension() == 2 ? MESH_Num_Edges(mesh_) : MESH_Num_Faces(mesh_));
      break;
    default:
      return 0;
    }
    break;


  case CELL:

    switch (ptype) {
    case Parallel_type::OWNED:
      return MSet_Num_Entries(OwnedCells);
      break;
    case Parallel_type::GHOST:
      return !serial_run ? MSet_Num_Entries(GhostCells) : 0;
      break;
    case Parallel_type::ALL:
      return ((manifold_dimension() == 2) ? MESH_Num_Faces(mesh_) : MESH_Num_Regions(mesh_));
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
Cell_type Mesh_MSTK::cell_get_type(const Entity_ID cellid) const
{
  MEntity_ptr cell = nullptr;
  int ival;
  Cell_type celltype;
  
  cell = cell_id_to_handle[cellid];
  
  MEnt_Get_AttVal(cell,celltype_att,&ival,NULL,NULL);
  celltype = (Cell_type) ival;

  return celltype;
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
void Mesh_MSTK::cell_get_faces_and_dirs_ordered(const Entity_ID cellid,
                                                Entity_ID_List *faceids,
                                                std::vector<int> *face_dirs) const 
{

  MEntity_ptr cell = nullptr;

  if (manifold_dimension() == 3) {

    int celltype = cell_get_type(cellid);
      
    if (celltype >= TET && celltype <= HEX) {
      int lid, nf;
      
      cell = cell_id_to_handle[cellid];

      List_ptr rfaces = MR_Faces((MRegion_ptr)cell);   
      nf = List_Num_Entries(rfaces);
      
      faceids->resize(nf);
      if (face_dirs)
        face_dirs->resize(nf);
      
      /* base face */
      
      MFace_ptr face0 = nullptr;
      int fdir0 = 0;
      
      if (celltype == TET || celltype == HEX) {
        face0 = List_Entry(rfaces,0);
        fdir0 = MR_FaceDir_i((MRegion_ptr)cell,0);
      }
      else if (celltype == PRISM) { /* Find the first triangular face */
        for (int i = 0; i < 5; ++i) {
          MFace_ptr face = List_Entry(rfaces,i);
          if (MF_Num_Edges(face) == 3) {
            face0 = face;
            fdir0 = MR_FaceDir_i((MRegion_ptr)cell,i);
            break;
          }
        }
      }
      else if (celltype == PYRAMID) { /* Find the quad face */
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
              (*faceids)[nf] = lid-1;
              
              if (face_dirs) {
                int fdir = (MR_FaceDir_i((MRegion_ptr)cell,i) == 1) ? 1 : -1;
                if (faceflip[lid-1]) fdir *= -1;
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
      (*faceids)[nf] = lid-1;
      
      if (face_dirs) {
        fdir0 = fdir0 ? 1 : -1;
        if (faceflip[lid-1]) fdir0 *= -1;
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
          (*faceids)[nf] = lid-1;
          
          if (face_dirs) {
            int fdir = (MR_FaceDir_i((MRegion_ptr)cell,i) == 1) ? 1 : -1;
            if (faceflip[lid-1]) fdir *= -1;
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
    }
    else 
      cell_get_faces_and_dirs_unordered(cellid,faceids,face_dirs);
  }
  else
    cell_get_faces_and_dirs_unordered(cellid,faceids,face_dirs);

}


void Mesh_MSTK::cell_get_faces_and_dirs_unordered(const Entity_ID cellid,
                                                  Entity_ID_List *faceids,
                                                  std::vector<int> *face_dirs) const
{
  MEntity_ptr cell;

  AMANZI_ASSERT(faceids != nullptr);

  cell = cell_id_to_handle[cellid];

  if (manifold_dimension() == 3) {
    int nrf;
    List_ptr rfaces;

    rfaces = MR_Faces((MRegion_ptr)cell);
    nrf = List_Num_Entries(rfaces);
    faceids->resize(nrf);

    Entity_ID_List::iterator itf = faceids->begin();
    for (int i = 0; i < nrf; ++i) {
      MFace_ptr face = List_Entry(rfaces,i);
      int lid = MEnt_ID(face);
      *itf = lid-1;  // assign to next spot by dereferencing iterator
      ++itf;
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
      face_dirs->resize(nrf);

      std::vector<int>::iterator itd = face_dirs->begin();
      for (int i = 0; i < nrf; ++i) {
        int lid = (*faceids)[i];
        int fdir = 2*MR_FaceDir_i((MRegion_ptr)cell,i) - 1;
        fdir = faceflip[lid] ? -fdir : fdir;
        *itd = fdir;  // assign to next spot by dereferencing iterator
        ++itd;
      }
    }
    
  }
  else {  // manifold_dimension() = 2; surface or 2D mesh
    int nfe;

    List_ptr fedges;
    fedges = MF_Edges((MFace_ptr)cell,1,0);
    nfe = List_Num_Entries(fedges);

    faceids->resize(nfe);

    Entity_ID_List::iterator itf = faceids->begin();
    for (int i = 0; i < nfe; ++i) {
      MEdge_ptr edge = List_Entry(fedges,i);
      int lid = MEnt_ID(edge);
      *itf = lid-1;  // assign to next spot by dereferencing iterator
      ++itf;
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
      face_dirs->resize(nfe);
      std::vector<int>::iterator itd = face_dirs->begin();
      for (int i = 0; i < nfe; ++i) { 
        int lid = (*faceids)[i];
        int fdir = 2*MF_EdgeDir_i((MFace_ptr)cell,i) - 1;
        fdir = faceflip[lid] ? -fdir : fdir;
        *itd = fdir;  // assign to next spot by dereferencing iterator
        itd++;
      }
    }

  }
}


void Mesh_MSTK::cell_get_faces_and_dirs_internal_(const Entity_ID cellid,
                                                  Entity_ID_List *faceids,
                                                  std::vector<int> *face_dirs,
                                                  const bool ordered) const 
{
  AMANZI_ASSERT(faces_initialized);

  if (ordered)
    cell_get_faces_and_dirs_ordered(cellid, faceids, face_dirs);
  else
    cell_get_faces_and_dirs_unordered(cellid, faceids, face_dirs);
}


void Mesh_MSTK::cell_get_edges_internal_(const Entity_ID cellid,
                                         Entity_ID_List *edgeids) const 
{
  AMANZI_ASSERT(edges_initialized);

  MEntity_ptr cell;

  AMANZI_ASSERT(edgeids != nullptr);

  cell = cell_id_to_handle[cellid];

  if (manifold_dimension() == 3) {
    int nre;
    List_ptr redges;

    redges = MR_Edges((MRegion_ptr)cell);
    nre = List_Num_Entries(redges);
    edgeids->resize(nre);

    Entity_ID_List::iterator ite = edgeids->begin();
    for (int i = 0; i < nre; ++i) {
      MEdge_ptr edge = List_Entry(redges,i);
      int lid = MEnt_ID(edge);
      *ite = lid-1;  // assign to next spot by dereferencing iterator
      ++ite;
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
    
  }
  else {  // manifold_dimension() = 2; surface or 2D mesh
    int nfe;

    List_ptr fedges;
    fedges = MF_Edges((MFace_ptr)cell,1,0);
    nfe = List_Num_Entries(fedges);

    edgeids->resize(nfe);

    Entity_ID_List::iterator ite = edgeids->begin();
    for (int i = 0; i < nfe; ++i) {
      MEdge_ptr edge = List_Entry(fedges,i);
      int lid = MEnt_ID(edge);
      *ite = lid-1;  // assign to next spot by dereferencing iterator
      ++ite;
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
void Mesh_MSTK::cell_get_nodes(const Entity_ID cellid, 
                               std::vector<Entity_ID> *nodeids) const
{
  MEntity_ptr cell;
  int nn, lid;

  AMANZI_ASSERT(nodeids != nullptr);

  cell = cell_id_to_handle[cellid];
      
  /* Reserved for next major MSTK release
  int vertids[MAXPV3];

  if (manifold_dimension() == 3)            // Volume mesh
    MR_VertexIDs(cell,&nn,vertids);
  else                                  // Surface mesh
    MF_VertexIDs(cell,1,0,&nn,vertids);

  nodeids->resize(nn);
  Entity_ID_List::iterator it = nodeids->begin();
  for (int i = 0; i < nn; ++i) {
    *it = vertids[i]-1;
    ++it;
  }
  */
    
  if (manifold_dimension() == 3) {                    // Volume mesh
    List_ptr rverts = MR_Vertices(cell);
 
    nn = List_Num_Entries(rverts);
    nodeids->resize(nn);
    Entity_ID_List::iterator it = nodeids->begin();
    for (int i = 0; i < nn; ++i) {
      lid = MEnt_ID(List_Entry(rverts,i));
      *it = lid-1;  // assign to next spot by dereferencing iterator
      ++it;
    }
    
    List_Delete(rverts);
  }
  else {                                 // Surface mesh
    List_ptr fverts = MF_Vertices(cell,1,0);

    nn = List_Num_Entries(fverts);
    nodeids->resize(nn);
    Entity_ID_List::iterator it = nodeids->begin();
    for (int i = 0; i < nn; ++i) {
      lid = MEnt_ID(List_Entry(fverts,i));
      *it = lid-1;  // assign to next spot by dereferencing iterator
      it++;
    }
    
    List_Delete(fverts);
  }
} 


void Mesh_MSTK::face_get_edges_and_dirs_internal_(const Entity_ID faceid,
                                                  Entity_ID_List *edgeids,
                                                  std::vector<int> *edge_dirs,
                                                  bool ordered) const
{
  AMANZI_ASSERT(edgeids != nullptr);

  AMANZI_ASSERT(faces_initialized);
  AMANZI_ASSERT(edges_initialized);

  MEntity_ptr face;

  face = face_id_to_handle[faceid];

  if (manifold_dimension() == 3) {
    int nfe;
    List_ptr fedges;

    fedges = MF_Edges((MFace_ptr)face,1,0);
    nfe = List_Num_Entries(fedges);
    edgeids->resize(nfe);

    Entity_ID_List::iterator ite = edgeids->begin();
    for (int i = 0; i < nfe; ++i) {
      MEdge_ptr edge = List_Entry(fedges,i);
      int lid = MEnt_ID(edge);
      *ite = lid-1;  // assign to next spot by dereferencing iterator
      ++ite;
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
      edge_dirs->resize(nfe);

      std::vector<int>::iterator itd = edge_dirs->begin();
      for (int i = 0; i < nfe; ++i) {
        int lid = (*edgeids)[i];
        int edir = 2*MF_EdgeDir_i((MFace_ptr)face,i) - 1;
        edir = edgeflip[lid] ? -edir : edir;
        *itd = edir;  // assign to next spot by dereferencing iterator
        ++itd;
      }
    }
    
  }
  else {  // manifold_dimension() = 2; surface or 2D mesh

    // face is same dimension as edge; just return the edge with a
    // direction of 1

    MEdge_ptr edge = (MEdge_ptr) face;

    edgeids->resize(1);
    (*edgeids)[0] = MEnt_ID(edge)-1;

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
void Mesh_MSTK::face_get_nodes(const Entity_ID faceid, 
                               std::vector<Entity_ID> *nodeids) const
{
  MEntity_ptr genface;
  int nn, lid;

  AMANZI_ASSERT(faces_initialized);

  AMANZI_ASSERT(nodeids != nullptr);

  genface = face_id_to_handle[faceid];
  
  if (manifold_dimension() == 3) {  // Volume mesh
    int dir = !faceflip[faceid];

    /* Reserved for next major MSTK release
    int vertids[MAXPV2];

    MF_VertexIDs((MFace_ptr) genface,dir,0,&nn,vertids);

    nodeids->resize(nn);
    for (int i = 0; i < nn; ++i) 
      (*nodeids)[i] = vertids[i]-1;
    */

    List_ptr fverts = MF_Vertices(genface,dir,0);
    AMANZI_ASSERT(fverts != nullptr);

    nn = List_Num_Entries(fverts);
    nodeids->resize(nn);
    Entity_ID_List::iterator it = nodeids->begin();

    for (int i = 0; i < nn; ++i) {
      lid = MEnt_ID(List_Entry(fverts,i));
      *it = lid-1;  // assign to next spot by dereferencing iterator
      ++it;
    }

    List_Delete(fverts);
    
  }
  else {                // Surface mesh or 2D mesh
    nodeids->resize(2);

    if (faceflip[faceid]) {
      (*nodeids)[0] = MEnt_ID(ME_Vertex(genface,1))-1;
      (*nodeids)[1] = MEnt_ID(ME_Vertex(genface,0))-1;
    }
    else {
      (*nodeids)[0] = MEnt_ID(ME_Vertex(genface,0))-1;
      (*nodeids)[1] = MEnt_ID(ME_Vertex(genface,1))-1;
    }
  }
}
  

//---------------------------------------------------------
// Get nodes of an edge
//---------------------------------------------------------
void Mesh_MSTK::edge_get_nodes(const Entity_ID edgeid,
                               Entity_ID *nodeid0, Entity_ID *nodeid1) const
{
  AMANZI_ASSERT(edges_initialized);

  MEdge_ptr edge = (MEdge_ptr) edge_id_to_handle[edgeid];

  if (edgeflip[edgeid]) {
    *nodeid0 = MEnt_ID(ME_Vertex(edge,1))-1;
    *nodeid1 = MEnt_ID(ME_Vertex(edge,0))-1;
  }
  else {
    *nodeid0 = MEnt_ID(ME_Vertex(edge,0))-1;
    *nodeid1 = MEnt_ID(ME_Vertex(edge,1))-1;
  }
}


//---------------------------------------------------------
// Cells of type 'ptype' connected to a node. This routine uses
// push_back on or near the partition boundary since we cannot tell at
// the outset how many entries will be put into the list
//---------------------------------------------------------
void Mesh_MSTK::node_get_cells(const Entity_ID nodeid, 
                               const Parallel_type ptype,
                               std::vector<Entity_ID> *cellids) const
{
  int idx, lid, nc;
  List_ptr cell_list;
  MEntity_ptr ment;

  AMANZI_ASSERT (cellids != nullptr);

  MVertex_ptr mv = (MVertex_ptr) vtx_id_to_handle[nodeid];
  
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

    if (manifold_dimension() == 3)
      cell_list = MV_Regions(mv);
    else
      cell_list = MV_Faces(mv);

    nc = List_Num_Entries(cell_list);
    cellids->resize(nc); // resize to maximum size possible
    Entity_ID_List::iterator it = cellids->begin();

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
    cellids->resize(n); // resize to the actual number of cells being returned

    List_Delete(cell_list);

    /*
  }
    */
}


//---------------------------------------------------------
// Faces of type 'ptype' connected to a node. This routine uses
// push_back on or near the partition boundary since we cannot tell at
// the outset how many entries will be put into the list
//---------------------------------------------------------
void Mesh_MSTK::node_get_faces(const Entity_ID nodeid, 
                               const Parallel_type ptype,
                               std::vector<Entity_ID> *faceids) const
{
  int idx, lid, n;
  List_ptr face_list;
  MEntity_ptr ment;

  AMANZI_ASSERT(faces_initialized);
  AMANZI_ASSERT(faceids != nullptr);

  MVertex_ptr mv = (MVertex_ptr) vtx_id_to_handle[nodeid];

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
    
    if (manifold_dimension() == 3)
      face_list = MV_Faces(mv);
    else
      face_list = MV_Edges(mv);
    
    int nf = List_Num_Entries(face_list);
    faceids->resize(nf); // resize to the maximum
    Entity_ID_List::iterator it = faceids->begin();
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
    faceids->resize(n); // resize to the actual number of faces being returned
    
    List_Delete(face_list);

    /*
  }
    */
}


//---------------------------------------------------------
// Edges of type 'ptype' connected to a node.
//---------------------------------------------------------
void Mesh_MSTK::node_get_edges(const Entity_ID nodeid, 
                               const Parallel_type ptype,
                               std::vector<Entity_ID> *edgeids) const
{
  int idx, lid, nc;
  List_ptr edge_list;
  MEntity_ptr ment;

  AMANZI_ASSERT (edgeids != nullptr);

  MVertex_ptr mv = (MVertex_ptr) vtx_id_to_handle[nodeid];
  
  // mesh vertex on a processor boundary may be connected to owned
  // and ghost cells. So depending on the requested cell type, we
  // may have to omit some entries

  edge_list = MV_Edges(mv);
  nc = List_Num_Entries(edge_list);

  edgeids->resize(nc); // resize to maximum size possible
  Entity_ID_List::iterator it = edgeids->begin();

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
  edgeids->resize(n); // resize to the actual number of cells being returned

  List_Delete(edge_list);
}


//---------------------------------------------------------
// Faces of type 'ptype' connected to an edge.
//---------------------------------------------------------
void Mesh_MSTK::edge_get_faces(const Entity_ID edgeid, 
                               const Parallel_type ptype,
                               std::vector<Entity_ID> *faceids) const
{
  int idx, lid, nc;
  List_ptr face_list;
  MEntity_ptr ment;

  AMANZI_ASSERT(faceids != nullptr && manifold_dimension() == 3);

  MEdge_ptr me = (MEdge_ptr) edge_id_to_handle[edgeid];
  face_list = ME_Faces(me);

  nc = List_Num_Entries(face_list);
  faceids->resize(nc); // resize to maximum size possible
  Entity_ID_List::iterator it = faceids->begin();

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
  faceids->resize(n); // resize to the actual number of cells being returned

  List_Delete(face_list);
}


//---------------------------------------------------------
// Cells of type 'ptype' connected to an edge. This routine uses
// push_back on or near the partition boundary since we cannot tell at
// the outset how many entries will be put into the list
//---------------------------------------------------------
void Mesh_MSTK::edge_get_cells(const Entity_ID edgeid, 
                               const Parallel_type ptype,
                               std::vector<Entity_ID> *cellids) const
{
  int idx, lid, nc;
  List_ptr cell_list;
  MEntity_ptr ment;

  AMANZI_ASSERT (cellids != nullptr);

  MEdge_ptr me = (MEdge_ptr) edge_id_to_handle[edgeid];
  
  // mesh edge on a processor boundary may be connected to owned
  // and ghost cells. So depending on the requested cell type, we
  // may have to omit some entries

  if (manifold_dimension() == 3)
    cell_list = ME_Regions(me);
  else
    cell_list = ME_Faces(me);

  nc = List_Num_Entries(cell_list);
  cellids->resize(nc); // resize to maximum size possible
  Entity_ID_List::iterator it = cellids->begin();

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
  cellids->resize(n); // resize to the actual number of cells being returned

  List_Delete(cell_list);
}


//---------------------------------------------------------
// Cells connected to a face
//---------------------------------------------------------
void Mesh_MSTK::face_get_cells_internal_(const Entity_ID faceid, 
                                         const Parallel_type ptype,
                                         std::vector<Entity_ID> *cellids) const
{
  AMANZI_ASSERT(faces_initialized);
  AMANZI_ASSERT(cellids != nullptr);
  cellids->clear();

  if (manifold_dimension() == 3) {
    MFace_ptr mf = (MFace_ptr) face_id_to_handle[faceid];
   
    List_ptr fregs = MF_Regions(mf);
    MRegion_ptr mr;
    if (ptype == Parallel_type::ALL) {      
      int idx = 0;
      while ((mr = List_Next_Entry(fregs,&idx)))
        cellids->push_back(MR_ID(mr)-1);
    }
    else {
      int idx = 0;
      while ((mr = List_Next_Entry(fregs,&idx))) {
        if (MEnt_PType(mr) == PGHOST) {
          if (ptype == Parallel_type::GHOST)
            cellids->push_back(MR_ID(mr)-1);
        }
        else
          if (ptype == Parallel_type::OWNED)
            cellids->push_back(MR_ID(mr)-1);
      }
    }
    List_Delete(fregs);

  }
  else {
    MEdge_ptr me = (MEdge_ptr) face_id_to_handle[faceid];

    List_ptr efaces = ME_Faces(me);
    MFace_ptr mf;
    if (ptype == Parallel_type::ALL) {
      int idx = 0;
      while ((mf = List_Next_Entry(efaces,&idx)))
        cellids->push_back(MF_ID(mf)-1);
    }
    else {
      int idx = 0;
      while ((mf = List_Next_Entry(efaces,&idx))) {
        if (MEnt_PType(mf) == PGHOST) {
          if (ptype == Parallel_type::GHOST)
            cellids->push_back(MF_ID(mf)-1);
        }
        else
          if (ptype == Parallel_type::OWNED)
            cellids->push_back(MF_ID(mf)-1);
      }
    }
    List_Delete(efaces);

  }
}
    

//-----------------------
// Same level adjacencies
//-----------------------

//---------------------------------------------------------
// Face connected neighboring cells of given cell. This routine uses
// push_back since we cannot tell at the outset how many entries will
// be put into the list
//---------------------------------------------------------
void Mesh_MSTK::cell_get_face_adj_cells(const Entity_ID cellid,
                                        const Parallel_type ptype,
                                        std::vector<Entity_ID> *fadj_cellids) const
{
  int lid;

  AMANZI_ASSERT(faces_initialized);

  AMANZI_ASSERT(fadj_cellids != nullptr);

  fadj_cellids->clear();

  if (manifold_dimension() == 3) {

    MRegion_ptr mr = (MRegion_ptr) cell_id_to_handle[cellid];

    List_ptr rfaces = MR_Faces(mr);
    int idx = 0;
    MFace_ptr mf;
    while ((mf = List_Next_Entry(rfaces,&idx))) {
      List_ptr fregs = MF_Regions(mf);
      int idx2 = 0;
      MRegion_ptr mr2;
      while ((mr2 = List_Next_Entry(fregs,&idx2))) {
        if (mr2 != mr) {
          if (MEnt_PType(mr2) == PGHOST) {
            if (ptype == Parallel_type::GHOST || ptype == Parallel_type::ALL) {
              lid = MEnt_ID(mr2);
              fadj_cellids->push_back(lid-1);
            }
          }
          else {
            if (ptype == Parallel_type::OWNED || ptype == Parallel_type::ALL) {
              lid = MEnt_ID(mr2);
              fadj_cellids->push_back(lid-1);
            }
          }
        }
      }
      List_Delete(fregs);
    }

    List_Delete(rfaces);

  }
  else if (manifold_dimension() == 2) {

    MFace_ptr mf = (MFace_ptr) cell_id_to_handle[cellid];

    List_ptr fedges = MF_Edges(mf,1,0);
    int idx = 0;
    MEdge_ptr me;
    while ((me = List_Next_Entry(fedges,&idx))) {
      List_ptr efaces = ME_Faces(me);
      int idx2 = 0;
      MFace_ptr mf2;
      while ((mf2 = List_Next_Entry(efaces,&idx2))) {
        if (mf2 != mf) {
          if (MEnt_PType(mf2) == PGHOST) {
            if (ptype == Parallel_type::GHOST || ptype == Parallel_type::ALL) {
              lid = MEnt_ID(mf2);
              fadj_cellids->push_back(lid-1);
            }
          }
          else {
            if (ptype == Parallel_type::OWNED || ptype == Parallel_type::ALL) {
              lid = MEnt_ID(mf2);
              fadj_cellids->push_back(lid-1);
            }
          }
        }
      }
      List_Delete(efaces);
    }

    List_Delete(fedges);
  }
}


//---------------------------------------------------------
// Node connected neighboring cells of given cell.  This routine uses
// push_back since we cannot tell at the outset how many entries will
// be put into the list
//---------------------------------------------------------
void Mesh_MSTK::cell_get_node_adj_cells(const Entity_ID cellid,
                                        const Parallel_type ptype,
                                        std::vector<Entity_ID> *nadj_cellids) const
{
  int lid, mkid;
  List_ptr cell_list;

  AMANZI_ASSERT(nadj_cellids != nullptr);

  nadj_cellids->clear();

  mkid = MSTK_GetMarker();

  cell_list = List_New(0);
  if (manifold_dimension() == 3) {

    MRegion_ptr mr = (MRegion_ptr) cell_id_to_handle[cellid];

    List_ptr rvertices = MR_Vertices(mr);
    int idx = 0;
    MVertex_ptr mv;
    while ((mv = List_Next_Entry(rvertices,&idx))) {
      List_ptr vregs = MV_Regions(mv);
      int idx2 = 0;
      MRegion_ptr mr2;
      while ((mr2 = List_Next_Entry(vregs,&idx2))) {
        if (mr2 != mr && !MEnt_IsMarked(mr2,mkid)) {
          MEnt_Mark(mr2,mkid);
          List_Add(cell_list,mr2);
          if (MEnt_PType(mr2) == PGHOST) {
            if (ptype == Parallel_type::GHOST || ptype == Parallel_type::ALL) {
              lid = MEnt_ID(mr2);
              nadj_cellids->push_back(lid-1);
            }
          }
          else {
            if (ptype == Parallel_type::OWNED || ptype == Parallel_type::ALL) {
              lid = MEnt_ID(mr2);
              nadj_cellids->push_back(lid-1);
            }
          }
        }
      }
      List_Delete(vregs);
    }

    List_Delete(rvertices);
    
  }
  else if (manifold_dimension() == 2) {

    MFace_ptr mf = (MFace_ptr) cell_id_to_handle[cellid];

    List_ptr fverts = MF_Vertices(mf,1,0);
    int idx = 0;
    MVertex_ptr mv;
    while ((mv = List_Next_Entry(fverts,&idx))) {
      List_ptr vfaces = MV_Faces(mv);
      int idx2 = 0;
      MFace_ptr mf2;
      while ((mf2 = List_Next_Entry(vfaces,&idx2))) {
        if (mf2 != mf && !MEnt_IsMarked(mf2,mkid)) {
          MEnt_Mark(mf2,mkid);
          List_Add(cell_list,mf2);
          if (MEnt_PType(mf2) == PGHOST) {
            if (ptype == Parallel_type::GHOST || ptype == Parallel_type::ALL) {
              lid = MEnt_ID(mf2);
              nadj_cellids->push_back(lid-1);
            }
          }
          else {
            if (ptype == Parallel_type::OWNED || ptype == Parallel_type::ALL) {
              lid = MEnt_ID(mf2);
              nadj_cellids->push_back(lid-1);
            }
          }
        }
      }
      List_Delete(vfaces);
    }

    List_Delete(fverts);

  }

  List_Unmark(cell_list,mkid);
  List_Delete(cell_list);
  MSTK_FreeMarker(mkid);
}


//---------------------------------------------------------
// Node coordinates - 3 in 3D and 2 in 2D
//---------------------------------------------------------
void Mesh_MSTK::node_get_coordinates(const Entity_ID nodeid, 
                                     AmanziGeometry::Point *ncoords) const
{    
  MEntity_ptr vtx;
  double coords[3];
  int spdim = space_dimension();
  
  AMANZI_ASSERT(ncoords != nullptr);

  vtx = vtx_id_to_handle[nodeid];

  MV_Coords(vtx,coords);
  ncoords->set(spdim,coords);
}


//---------------------------------------------------------
// Coordinates of cells in standard order (Exodus II convention)
// STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES IN 3D
// For a general polyhedron this will return the node coordinates in
// arbitrary order
// Number of nodes is vector size divided by number of spatial dimensions
//---------------------------------------------------------
void Mesh_MSTK::cell_get_coordinates(const Entity_ID cellid, 
                                     std::vector<AmanziGeometry::Point> *ccoords) const
{    
  MEntity_ptr cell;
  double coords[3];
  int nn;
  int spdim = space_dimension(), celldim = manifold_dimension();

  AMANZI_ASSERT(ccoords != NULL);

  cell = cell_id_to_handle[cellid];
      
  if (celldim == 3) {
    List_ptr rverts = MR_Vertices(cell);
    
    nn = List_Num_Entries(rverts);

    ccoords->resize(nn);
    std::vector<AmanziGeometry::Point>::iterator it = ccoords->begin();

    for (int i = 0; i < nn; ++i) {
      MV_Coords(List_Entry(rverts,i),coords);
      it->set(spdim,coords); // same as (*it).set()
      ++it; 
    }    

    List_Delete(rverts);
  }
  else if (celldim == 2) {
    List_ptr fverts = MF_Vertices(cell,1,0);

    nn = List_Num_Entries(fverts);

    ccoords->resize(nn);
    std::vector<AmanziGeometry::Point>::iterator it = ccoords->begin();

    for (int i = 0; i < nn; ++i) {
      MV_Coords(List_Entry(fverts,i),coords);
      it->set(spdim,coords); // same as (*it).set()
      ++it;
    }

    List_Delete(fverts);
  }
}


//---------------------------------------------------------
// Face coordinates - conventions same as face_to_nodes call 
// Number of nodes is the vector size divided by number of spatial dimensions
//---------------------------------------------------------
void Mesh_MSTK::face_get_coordinates(const Entity_ID faceid, 
                                     std::vector<AmanziGeometry::Point> *fcoords) const
{
  MEntity_ptr genface;
  double coords[3];
  int spdim = space_dimension(), celldim = manifold_dimension();

  AMANZI_ASSERT(faces_initialized);
  AMANZI_ASSERT(fcoords != NULL);

  genface = face_id_to_handle[faceid];

  if (celldim == 3) {
    int dir = !faceflip[faceid];

    List_ptr fverts = MF_Vertices((MFace_ptr) genface,dir,0);

    int nn = List_Num_Entries(fverts);
    fcoords->resize(nn);
    std::vector<AmanziGeometry::Point>::iterator it = fcoords->begin();
        
    for (int i = 0; i < nn; ++i) {
      MV_Coords(List_Entry(fverts,i),coords);
      it->set(spdim,coords); // same as (*it).set()
      ++it;
    }

    List_Delete(fverts);
  }
  else { // Planar mesh or Surface mesh embedded in 3D
    MVertex_ptr ev[2];
    if (!faceflip[faceid]) {
      ev[0] = ME_Vertex((MEdge_ptr)genface,0);
      ev[1] = ME_Vertex((MEdge_ptr)genface,1);
    }
    else {
      ev[1] = ME_Vertex((MEdge_ptr)genface,0);
      ev[0] = ME_Vertex((MEdge_ptr)genface,1);
    }

    fcoords->resize(2);

    MV_Coords(ev[0],coords);
    ((*fcoords)[0]).set(spdim,coords);
        
    MV_Coords(ev[1],coords);
    ((*fcoords)[1]).set(spdim,coords);
  }
}
  

//---------------------------------------------------------
// Modify a node's coordinates
//---------------------------------------------------------
void Mesh_MSTK::node_set_coordinates(const AmanziMesh::Entity_ID nodeid, 
                                     const double *coords)
{
  MVertex_ptr v = vtx_id_to_handle[nodeid];
  MV_Set_Coords(v,(double *)coords);
}


void Mesh_MSTK::node_set_coordinates(const AmanziMesh::Entity_ID nodeid, 
                                     const AmanziGeometry::Point coords)
{
  MVertex_ptr v = vtx_id_to_handle[nodeid];

  double coordarray[3] = {0.0,0.0,0.0};
  for (int i = 0; i < space_dimension(); i++)
    coordarray[i] = coords[i];

  MV_Set_Coords(v,(double *)coordarray);
}


//---------------------------------------------------------
// Private routine creating mesh sets for GM regions
//---------------------------------------------------------
MSet_ptr Mesh_MSTK::build_set(const Teuchos::RCP<const AmanziGeometry::Region>& region,
                              const Entity_kind kind) const
{
  int celldim = manifold_dimension();
  int space_dim = space_dimension();
  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm = geometric_model();

  // Modify region/set name by prefixing it with the type of entity requested

  std::string internal_name = internal_name_of_set(region,kind);

  // Create entity set based on the region defintion  

  MSet_ptr mset;
  MType enttype;
  switch (kind) {      
  case CELL:    // cellsets      

    enttype = (celldim == 3) ? MREGION : MFACE;
    mset = MSet_New(mesh_,internal_name.c_str(),enttype);
      
    if (region->type() == AmanziGeometry::BOX ||
        region->type() == AmanziGeometry::CYLINDER ||
        region->type() == AmanziGeometry::COLORFUNCTION) {

      int ncell = num_entities(CELL, Parallel_type::ALL);              

      for (int icell = 0; icell < ncell; icell++)
        if (region->inside(cell_centroid(icell)))
          MSet_Add(mset,cell_id_to_handle[icell]);

    }
    else if (region->type() == AmanziGeometry::ALL)  {

      int ncell = num_entities(CELL, Parallel_type::ALL);              

      for (int icell = 0; icell < ncell; icell++)
        MSet_Add(mset,cell_id_to_handle[icell]);

    }
    else if (region->type() == AmanziGeometry::ENUMERATED)  {
      auto rgn = Teuchos::rcp_static_cast<const AmanziGeometry::RegionEnumerated>(region);
      int ncell = num_entities(CELL, Parallel_type::ALL);

      for (int icell = 0; icell < ncell; icell++) {
        Entity_ID gid = MEnt_GlobalID(cell_id_to_handle[icell]);
        for (const auto& jset : rgn->entities()) {
          if (jset == gid) {
            MSet_Add(mset,cell_id_to_handle[icell]);
            break;
          }
        }
      }

    }
    else if (region->type() == AmanziGeometry::POINT) {
      AmanziGeometry::Point vpnt(space_dim);
      AmanziGeometry::Point rgnpnt(space_dim);

      rgnpnt = Teuchos::rcp_static_cast<const AmanziGeometry::RegionPoint>(region)->point();
        
      int nnode = num_entities(NODE, Parallel_type::ALL);
      double mindist2 = 1.e+16;
      int minnode = -1;
        
      int inode;
      for (inode = 0; inode < nnode; inode++) {
        node_get_coordinates(inode, &vpnt);                  
        double dist2 = (vpnt-rgnpnt)*(vpnt-rgnpnt);
 
        if (dist2 < mindist2) {
          mindist2 = dist2;
          minnode = inode;
          if (mindist2 <= 1.0e-14)
            break;
        }
      }

      Entity_ID_List cells;
      node_get_cells(minnode, Parallel_type::ALL, &cells);
      
      int ncells = cells.size();
      for (int ic = 0; ic < ncells; ic++) {
        Entity_ID icell = cells[ic];
        
        // Check if point is contained in cell            
        if (point_in_cell(rgnpnt, icell))
          MSet_Add(mset, cell_id_to_handle[icell]);
      }

      // finally check all cells, typical for anisotropic meshes
      if (MSet_Num_Entries(mset) == 0) {
        int ncells_wghost = num_entities(CELL, Parallel_type::ALL);
        for (int c = 0; c < ncells_wghost; ++c)
          if (point_in_cell(rgnpnt, c)) 
            MSet_Add(mset, cell_id_to_handle[c]);
      }

    }
    else if ((region->type() == AmanziGeometry::PLANE)||
             (region->type() == AmanziGeometry::POLYGON)) {

      if (celldim == 2) {

        int ncells = num_entities(CELL, Parallel_type::ALL);              
        for (int ic = 0; ic < ncells; ic++) {

          std::vector<AmanziGeometry::Point> ccoords(space_dim);

          cell_get_coordinates(ic, &ccoords);

          bool on_plane = true;
          for (int j = 0; j < ccoords.size(); ++j) {
            if (!region->inside(ccoords[j])) {
              on_plane = false;
              break;
            }
          }
                  
          if (on_plane)
            MSet_Add(mset,cell_id_to_handle[ic]);
        }
      }

    }
    else if (region->type() == AmanziGeometry::LOGICAL) {
      // will process later in this subroutine
    }
    else if (region->type() == AmanziGeometry::LABELEDSET) {
      if (parent_mesh_.get() != NULL) {
        auto lsrgn = Teuchos::rcp_static_cast<const AmanziGeometry::RegionLabeledSet>(region);
        std::string label = lsrgn->label();
        std::string entity_type = lsrgn->entity_str();

        if (kind == CELL && entity_type != "FACE") {
          if (vo_.get() && vo_->os_OK(Teuchos::VERB_MEDIUM)) {
            *(vo_->os()) << "Found labeled set region \"" << region->name() 
                             << "\" but it contains entities of type " << entity_type 
                             << ", not the requested type.\n";
          }
        } 
        else {
          MSet_ptr mset2;

          // Build set on a fly.
          // -- first request the set on the parent to make sure it was constructed in MSTK in all cases.
          AmanziMesh::Entity_ID_List parent_ids;
          parent_mesh_->get_set_entities(region->name(), FACE, Parallel_type::ALL, &parent_ids);

          int ival;
          double rval;
          void *pval = nullptr;
          MAttrib_ptr att = (manifold_dimension() == 3) ? rparentatt : fparentatt;

          std::string internal_parent_name = internal_name_of_set(region,FACE);
          mset2 = MESH_MSetByName(parent_mesh_->mesh_, internal_parent_name.c_str());

          for (int c = 0; c < num_entities(CELL, Parallel_type::ALL); ++c) {
            auto ment = cell_id_to_handle[c];
            MEnt_Get_AttVal(ment,att,&ival,&rval,&pval);
            if (MSet_Locate(mset2, pval) >= 0) MSet_Add(mset, ment);
          }

          // Due to the parallel partitioning its possible that this
          // set is not on this processor
          if (!mset) {
            if (comm_->NumProc() == 1) {
              Errors::Message msg;
              msg << "Could not find labeled set \"" << label 
                  << "\" in mesh file to initialize mesh set \"" << region->name()  
                  << "\". Verify mesh file.";
              amanzi_throw(msg);
            }
          }
        }
      }
      else {
        // Just retrieve and return the set

        auto lsrgn = Teuchos::rcp_static_cast<const AmanziGeometry::RegionLabeledSet>(region);
        std::string label = lsrgn->label();
        std::string entity_type = lsrgn->entity_str();

        if (entity_type != "CELL") {
          Errors::Message mesg;
          mesg << "Entity type of labeled set region \"" << region->name() 
               << "\" and build_set request (cell) do not match";
          Exceptions::amanzi_throw(mesg);
        }

        mset = MESH_MSetByName(mesh_,internal_name.c_str());

        std::string other_internal_name = other_internal_name_of_set(region,kind);
        MSet_ptr mset2 = MESH_MSetByName(mesh_,other_internal_name.c_str());

        if (mset) {
          if (mset2) {
            std::stringstream mesg_stream;
            mesg_stream << "Exodus II file has element block and element set with the same ID " << label << " - Amanzi cannot handle this case.";
            Errors::Message mesg(mesg_stream.str());
            Exceptions::amanzi_throw(mesg);
          }
        } 
        else {
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
    }
    else {
      if (vo_.get() && vo_->os_OK(Teuchos::VERB_HIGH)) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *(vo_->os()) << "Requested CELLS on region " << region->name() 
            << " of type " << region->type()  
            << " and dimension " << region->manifold_dimension() << ".\n"
            << "This request will result in an empty set";
      }
    }
      
    break;
      
  case FACE:  // sidesets

    //
    // Commented out so that we can ask for a face set in a 3D box
    //
    //          if (region->manifold_dimension() != celldim-1) 
    //            {
    //              std::cerr << "No region of dimension " << celldim-1 << " defined in geometric model" << std::endl;
    //              std::cerr << "Cannot construct cell set from region " << setname << std::endl;
    //            }

    enttype = (celldim == 3) ? MFACE : MEDGE;
    mset = MSet_New(mesh_,internal_name.c_str(),enttype);

    if (region->type() == AmanziGeometry::BOX ||
        region->type() == AmanziGeometry::CYLINDER) {
      int nface = num_entities(FACE, Parallel_type::ALL);

      if (nface > 0) { 
        if (! kdtree_faces_initialized_) {
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
          if (region->inside(face_centroid(iface)))
            MSet_Add(mset, face_id_to_handle[iface]);
        }
      }
    }
    else if (region->type() == AmanziGeometry::ALL)  {

      int nface = num_entities(FACE, Parallel_type::ALL);
        
      for (int iface = 0; iface < nface; iface++) {
        MSet_Add(mset, face_id_to_handle[iface]);
      }
    }
    else if (region->type() == AmanziGeometry::PLANE ||
             region->type() == AmanziGeometry::POLYGON) {

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
                  
        if (on_plane)
          MSet_Add(mset,face_id_to_handle[iface]);
      }

    }
    else if (region->type() == AmanziGeometry::LABELEDSET) {
      // Just retrieve and return the set

      Teuchos::RCP<const AmanziGeometry::RegionLabeledSet> lsrgn =
          Teuchos::rcp_static_cast<const AmanziGeometry::RegionLabeledSet>(region);
      std::string label = lsrgn->label();
      std::string entity_type = lsrgn->entity_str();

      if (entity_type != "FACE") {
        Errors::Message mesg;
        mesg << "Entity type of labeled set region \"" << region->name() 
             << "\" and build_set request (face) do not match";
        Exceptions::amanzi_throw(mesg);
      }

      mset = MESH_MSetByName(mesh_,internal_name.c_str());
    }
    else if (region->type() == AmanziGeometry::LOGICAL) {
      // Will handle it later in the routine
    }
    else if (region->type() == AmanziGeometry::BOUNDARY)  {

      const Epetra_Map& fmap = face_map(true); 
      const Epetra_Map& map = exterior_face_map(true); 

      int nface = map.NumMyElements(); 

      for (int iface = 0; iface < nface; iface++) {
        int lid = fmap.LID(map.GID(iface));
        MSet_Add(mset,face_id_to_handle[lid]);
      }
    }
    else {
      Teuchos::RCP<const VerboseObject> verbobj = verbosity_obj();
      if (verbobj.get() && verbobj->os_OK(Teuchos::VERB_HIGH)) {
        Teuchos::OSTab tab = verbobj->getOSTab();
        *(verbobj->os()) << "Requested FACES on region " << region->name()
            << " of type " << region->type() << " and dimension "
            << region->manifold_dimension() << ".\n" 
            << "This request will result in an empty set\n";
      }
    }
    break;

  case EDGE: // Edgesets

    enttype = MEDGE;
    mset = MSet_New(mesh_,internal_name.c_str(),enttype);

    if (region->type() == AmanziGeometry::BOX ||
        region->type() == AmanziGeometry::PLANE ||
        region->type() == AmanziGeometry::POLYGON) {

      int nedge = num_entities(EDGE, Parallel_type::ALL);

      for (int iedge = 0; iedge < nedge; iedge++) {
        const auto& epnt = edge_centroid(iedge);
                  
        if (region->inside(epnt)) {
          MSet_Add(mset,edge_id_to_handle[iedge]);
        }
      }
    }
    else if (region->type() == AmanziGeometry::LOGICAL) {
      // We will handle it later in the routine
    }
    else {
      Teuchos::RCP<const VerboseObject> verbobj = verbosity_obj();
      if (verbobj.get() && verbobj->os_OK(Teuchos::VERB_HIGH)) {
        Teuchos::OSTab tab = verbobj->getOSTab();
        *(verbobj->os()) << "Requested EDGEs on region " << region->name() 
            << " of type " << region->type() << " and dimension " 
            << region->manifold_dimension() << ".\n" 
            << "This request will result in an empty set\n";
      }
    }
      
    break;

  case NODE: // Nodesets

    enttype = MVERTEX;
    mset = MSet_New(mesh_,internal_name.c_str(),enttype);

    if (region->type() == AmanziGeometry::BOX ||
        region->type() == AmanziGeometry::PLANE ||
        region->type() == AmanziGeometry::POLYGON ||
        region->type() == AmanziGeometry::CYLINDER ||
        region->type() == AmanziGeometry::POINT) {

      int nnode = num_entities(NODE, Parallel_type::ALL);

      for (int inode = 0; inode < nnode; inode++) {

        AmanziGeometry::Point vpnt(space_dim);
        node_get_coordinates(inode, &vpnt);
                  
        if (region->inside(vpnt)) {
          MSet_Add(mset,vtx_id_to_handle[inode]);

          // Only one node per point region
          if (region->type() == AmanziGeometry::POINT)
            break;      
        }
      }
    }
    else if (region->type() == AmanziGeometry::ALL)  {

      int nnode = num_entities(NODE, Parallel_type::ALL);

      for (int inode = 0; inode < nnode; inode++)
        MSet_Add(mset,vtx_id_to_handle[inode]);

    }
    else if (region->type() == AmanziGeometry::LABELEDSET) {
      // Just retrieve and return the set

      Teuchos::RCP<const AmanziGeometry::RegionLabeledSet> lsrgn =
          Teuchos::rcp_static_cast<const AmanziGeometry::RegionLabeledSet>(region);
      std::string label = lsrgn->label();
      std::string entity_type = lsrgn->entity_str();

      if (entity_type != "FACE") {
        Errors::Message mesg;
        mesg << "Entity type of labeled set region \"" << region->name() 
             << "\" and build_set request (face) do not match";
        Exceptions::amanzi_throw(mesg);
      }

      mset = MESH_MSetByName(mesh_,internal_name.c_str());
    }
    else if (region->type() == AmanziGeometry::LOGICAL) {
      // We will handle it later in the routine
    }
    else if (region->type() == AmanziGeometry::BOUNDARY)  {

      const Epetra_Map& vmap = node_map(true); 
      const Epetra_Map& map = exterior_node_map(true); 

      int nnode = map.NumMyElements(); 

      for (int inode = 0; inode < nnode; inode++) {
        int lid = vmap.LID(map.GID(inode));
        MSet_Add(mset,vtx_id_to_handle[lid]);
      }
    }
    else {
      Teuchos::RCP<const VerboseObject> verbobj = verbosity_obj();
      if (verbobj.get() && verbobj->os_OK(Teuchos::VERB_HIGH)) {
        Teuchos::OSTab tab = verbobj->getOSTab();
        *(verbobj->os()) << "Requested POINTS on region " << region->name() 
            << " of type " << region->type() << " and dimension " 
            << region->manifold_dimension() << ".\n" 
            << "This request will result in an empty set\n";
      }
    }
      
    break;

  default:
    break;
  }


  if (region->type() == AmanziGeometry::LOGICAL) {
    Teuchos::RCP<const AmanziGeometry::RegionLogical> boolregion =
        Teuchos::rcp_static_cast<const AmanziGeometry::RegionLogical>(region);
    const std::vector<std::string> region_names = boolregion->component_regions();
    int nreg = region_names.size();
    
    std::vector<MSet_ptr> msets;
    std::vector<Teuchos::RCP<const AmanziGeometry::Region> > regions;
    
    for (int r = 0; r < nreg; r++) {
      Teuchos::RCP<const AmanziGeometry::Region> rgn1 = gm->FindRegion(region_names[r]);
      regions.push_back(rgn1);

      // Did not find the region
      if (rgn1 == Teuchos::null) {
        std::stringstream mesg_stream;
        mesg_stream << "Geometric model has no region named " << 
          region_names[r];
        Errors::Message mesg(mesg_stream.str());
        Exceptions::amanzi_throw(mesg);
      }
        
      internal_name = internal_name_of_set(rgn1,kind);
      MSet_ptr mset1 = MESH_MSetByName(mesh_,internal_name.c_str());
      if (!mset1)        
        mset1 = build_set(rgn1,kind);  // Recursive call

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

        Errors::Message 
          mesg("Amanzi cannot operate on sets of different entity types");
        Exceptions::amanzi_throw(mesg);               
      }
    
    int mkid = MSTK_GetMarker();
      
    if (boolregion->operation() == AmanziGeometry::COMPLEMENT) {
      
      for (int ms = 0; ms < msets.size(); ms++)
        MSet_Mark(msets[ms],mkid);
      
      int idx = 0;
      switch (enttype) {
      case MREGION:
        MRegion_ptr mr;
        while ((mr = MESH_Next_Region(mesh_,&idx))) 
          if (!MEnt_IsMarked(mr,mkid))
            MSet_Add(mset,mr);
        break;        
      case MFACE: 
        MFace_ptr mf;
        while ((mf = MESH_Next_Face(mesh_,&idx))) 
          if (!MEnt_IsMarked(mf,mkid))
            MSet_Add(mset,mf);
        break;
      case MEDGE:
        MEdge_ptr me;
        while ((me = MESH_Next_Edge(mesh_,&idx))) 
          if (!MEnt_IsMarked(me,mkid))
            MSet_Add(mset,me);
      case MVERTEX:
        MVertex_ptr mv;
        while ((mv = MESH_Next_Vertex(mesh_,&idx))) 
          if (!MEnt_IsMarked(mv,mkid))
            MSet_Add(mset,mv);
        break;
      default:
        break;
      }
      
      for (int ms = 0; ms < msets.size(); ms++)
        MSet_Unmark(msets[ms],mkid);

    }
    else if (boolregion->operation() == AmanziGeometry::UNION) {

      for (int ms = 0; ms < msets.size(); ms++) {
        MEntity_ptr ment;
        int idx = 0;
        while ((ment = MSet_Next_Entry(msets[ms],&idx))) 
          if (!MEnt_IsMarked(ment,mkid)) {
            MSet_Add(mset,ment);
            MEnt_Mark(ment,mkid);
          }
      }
      MSet_Unmark(mset,mkid);
      
    }
    else if (boolregion->operation() == AmanziGeometry::SUBTRACT) {

      /* Mark entities in all sets except the first */
      
      for (int ms = 1; ms < msets.size(); ms++)
        MSet_Mark(msets[ms],mkid);
      
      /* Look for entities in the first set but not in 
         any of the other sets */
      MEntity_ptr ment;
      int idx = 0;
      while ((ment = MSet_Next_Entry(msets[0],&idx))) 
        if (!MEnt_IsMarked(ment,mkid)) {
          MSet_Add(mset,ment);
          MEnt_Mark(ment,mkid);
        }
      
      for (int ms = 1; ms < msets.size(); ms++)
        MSet_Unmark(msets[ms],mkid);

    }
    else if (boolregion->operation() == AmanziGeometry::INTERSECT) {

      /* Can't do this using markers alone - need attributes */
      
      MAttrib_ptr matt = MAttrib_New(mesh_,"XSECTATT",INT,MALLTYPE);
      
      for (int ms = 0; ms < msets.size(); ms++) {
        MEntity_ptr ment;
        int idx = 0;
        while ((ment = MSet_Next_Entry(msets[ms],&idx))) {
          int ival;
          double rval;
          void *pval;
          MEnt_Get_AttVal(ment,matt,&ival,&rval,&pval);
          ival++;
          MEnt_Set_AttVal(ment,matt,ival,rval,pval);
        }
      }
      
      for (int ms = 0; ms < msets.size(); ms++) {
        MEntity_ptr ment;
        int idx = 0;
        while ((ment = MSet_Next_Entry(msets[ms],&idx))) {
          int ival;
          double rval;
          void *pval;
          MEnt_Get_AttVal(ment,matt,&ival,&rval,&pval);
          if ((ival == msets.size()) && !MEnt_IsMarked(ment,mkid)) {
            /* entity is in all sets */
            MSet_Add(mset,ment);
            MEnt_Mark(ment,mkid);
          }
        }
      }
      
      MSet_Unmark(mset,mkid);
      
      for (int ms = 0; ms < msets.size(); ms++) {
        MEntity_ptr ment;
        int idx = 0;
        while ((ment = MSet_Next_Entry(msets[ms],&idx)))
          MEnt_Rem_AttVal(ment,matt);
      }
      MAttrib_Delete(matt);
    }

    MSTK_FreeMarker(mkid);

    for (int ms = 0; ms < msets.size(); ms++) {
      MSet_Unmark(msets[ms],mkid);
      if (regions[ms]->lifecycle() == AmanziGeometry::TEMPORARY)
        MSet_Delete(msets[ms]);
    }
  }
  return mset;
}


//---------------------------------------------------------
// Get list of entities of type 'category' in set specified by setname
//---------------------------------------------------------
void Mesh_MSTK::get_set_entities_and_vofs(const std::string& setname, 
                                          const Entity_kind kind, 
                                          const Parallel_type ptype, 
                                          std::vector<Entity_ID> *setents,
                                          std::vector<double> *vofs) const
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
    mesg_stream << "Geometric model has no region named \"" << setname <<"\"\n";
    Errors::Message mesg(mesg_stream.str());
    Exceptions::amanzi_throw(mesg);
  }


  std::string internal_name = internal_name_of_set(rgn,kind);

  // If region is of type labeled set, a mesh set should have been
  // initialized from the input file. If region requires volume 
  // fractions or is a segment, use base class capabilities to 
  // build a mesh set. Otherwise, if a mesh set exists, search 
  // the database for it. This is a two step procedure, which shows
  // probably defficiency of the internal naming convention (KL).
  // Finally, build a new mesh set for the region.
  
  if (rgn->type() == AmanziGeometry::LABELEDSET && parent_mesh_.get() == NULL) {
    auto lsrgn = Teuchos::rcp_static_cast<const AmanziGeometry::RegionLabeledSet>(rgn);
    std::string label = lsrgn->label();
    std::string entity_type = lsrgn->entity_str();

    if ((kind == CELL && entity_type != "CELL") ||
        (kind == FACE && entity_type != "FACE") ||
        (kind == NODE && entity_type != "NODE")) {
      if (verbobj.get() && verbobj->os_OK(Teuchos::VERB_MEDIUM)) {
        *(verbobj->os()) << "Found labeled set region \"" << setname 
                         << "\" but it contains entities of type " << entity_type 
                         << ", not the requested type\n";
      }
    } 
    else {
      mset1 = MESH_MSetByName(mesh_, internal_name.c_str());
      
      if (!mset1 && kind == CELL) {
        // Since both element blocks and cell sets are referenced
        // with the region type 'Labeled Set' and Entity kind 'Cell'
        // we have to account for both possibilities. NOTE: THIS
        // MEANS THAT IF AN ELEMENT BLOCK AND ELEMENT SET HAVE THE
        // SAME ID, ONLY THE ELEMENT BLOCK WILL GET PICKED UP - WE
        // CHECKED FOR THIS IN BUILD SET

        std::string internal_name2 = other_internal_name_of_set(rgn,kind);
        mset1 = MESH_MSetByName(mesh_, internal_name2.c_str());
      }

      // Due to the parallel partitioning its possible that this
      // set is not on this processor
      
      if (!mset1) {
        if (comm_->NumProc() == 1) {
          Errors::Message msg;
          msg << "Could not find labeled set \"" << label 
              << "\" in mesh file to initialize mesh set \"" << setname 
              << "\". Verify mesh file.";
          amanzi_throw(msg);
        }
      }
    }
  }
  else if ((rgn->type() == AmanziGeometry::BOX_VOF) || 
           (rgn->type() == AmanziGeometry::LINE_SEGMENT)) {
    // Call routine from the base class and exit.
    Mesh::get_set_entities_box_vofs_(rgn, kind, ptype, setents, vofs);
    return;
  }
  else {
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
        MSet_Name(mset1,setname1);
              
        if (MSet_EntDim(mset1) == entdim &&
            strcmp(setname1,internal_name.c_str()) == 0)
          break;
      }
    }
  }

  // All attempts to find the set failed so it must not exist - build it

  if (mset1 == NULL) {
    mset1 = build_set(rgn, kind);
  }

  // Check if no processor got any mesh entities

  int nent_loc = (mset1 == NULL) ? 0 : MSet_Num_Entries(mset1);

  setents->resize(nent_loc);
  Entity_ID_List::iterator it = setents->begin();

  if (nent_loc) {
    nent_loc = 0; // reset and count to get the real number

    switch (ptype) {
    case Parallel_type::OWNED:
      idx = 0;
      while ((ment = MSet_Next_Entry(mset1,&idx))) {
        if (MEnt_PType(ment) != PGHOST) {
          *it = MEnt_ID(ment)-1;  // assign to next spot by dereferencing iterator
          ++it;
          ++nent_loc;
        }
      }
      break;
    case Parallel_type::GHOST:
      idx = 0;
      while ((ment = MSet_Next_Entry(mset1,&idx))) {
        if (MEnt_PType(ment) == PGHOST) {
          *it = MEnt_ID(ment)-1;  // assign to next spot by dereferencing iterator
          ++it;
          ++nent_loc;
        }
      }
      break;
    case Parallel_type::ALL:
      idx = 0;
      while ((ment = MSet_Next_Entry(mset1,&idx))) {
        *it = MEnt_ID(ment)-1;  // assign to next spot by dereferencing iterator
        ++it;
        ++nent_loc;
      }
      break;
    default:
      {}
    }
    
    setents->resize(nent_loc);
  }
}


//---------------------------------------------------------
// Parent entity in the source mesh if mesh was derived from another mesh
//---------------------------------------------------------
Entity_ID Mesh_MSTK::entity_get_parent(const Entity_kind kind, const Entity_ID entid) const
{
  int ival;
  double rval;
  void *pval = nullptr;
  MEntity_ptr ment = nullptr;
  MAttrib_ptr att = nullptr;

  switch(kind) {
  case CELL:
    att = (manifold_dimension() == 3) ? rparentatt : fparentatt;
    ment = (MEntity_ptr) cell_id_to_handle[entid];
    break;
  case FACE:
    att = (manifold_dimension() == 3) ? fparentatt : eparentatt;
    ment = (MEntity_ptr) face_id_to_handle[entid];
    break;
  case EDGE:
    att = eparentatt;
    ment = (MEntity_ptr) edge_id_to_handle[entid];
    break;
  case NODE:
    if (!vparentatt) return 0;
    att = vparentatt;
    ment = (MEntity_ptr) vtx_id_to_handle[entid];
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
// Epetra map for cells - basically a structure specifying the global
// IDs of cells owned or used by this processor. This helps Epetra
// understand inter-partition dependencies of the data.
//
// Amanzi/Epetra want global IDs to start at 0
//---------------------------------------------------------
void Mesh_MSTK::init_cell_map()
{
  int *cell_gids;
  int ncell, idx, i;
  MEntity_ptr ment;

  if (!serial_run) {

    // For parallel runs create map without and with ghost cells included
    // Also, put in owned cells before the ghost cells

    int nowned = MSet_Num_Entries(OwnedCells);
    int nnotowned = MSet_Num_Entries(GhostCells);

    cell_gids = new int[nowned+nnotowned];
    
    ncell = nowned;
   
    idx = 0; i = 0;
    while ((ment = MSet_Next_Entry(OwnedCells,&idx)))
      cell_gids[i++] = MEnt_GlobalID(ment)-1;

    cell_map_wo_ghosts_ = new Epetra_Map(-1,ncell,cell_gids,0,*comm_);
    

    ncell += nnotowned;

    idx = 0; 
    while ((ment = MSet_Next_Entry(GhostCells,&idx)))
      cell_gids[i++] = MEnt_GlobalID(ment)-1;
    
    cell_map_w_ghosts_ = new Epetra_Map(-1,ncell,cell_gids,0,*comm_);

  }
  else {    
    ncell = MSet_Num_Entries(OwnedCells);
    cell_gids = new int[ncell];

    idx = 0; i = 0;
    while ((ment = MSet_Next_Entry(OwnedCells,&idx)))      
      cell_gids[i++] = MEnt_ID(ment)-1;

    cell_map_wo_ghosts_ = new Epetra_Map(-1,ncell,cell_gids,0,*comm_);
  }

  delete [] cell_gids;
}


//---------------------------------------------------------
// Epetra map for faces - basically a structure specifying the global
// IDs of faces owned or used by this processor. This helps Epetra
// understand inter-partition dependencies of the data.
//
// Amanzi/Epetra want global IDs to start at 0
//---------------------------------------------------------
void Mesh_MSTK::init_face_map()
{
  int *face_gids = nullptr, *extface_gids = nullptr;
  int nface = -1, n_extface = -1, idx, i = -1, j = -1;
  MEntity_ptr ment = nullptr;

  if (!serial_run) {

    // For parallel runs create map without and with ghost cells included
    // Also, put in owned cells before the ghost cells
    // Additionally, create a map of exterior faces only

    int nowned = MSet_Num_Entries(OwnedFaces);
    int nnotowned = MSet_Num_Entries(NotOwnedFaces);

    face_gids = new int[nowned+nnotowned];
    extface_gids = new int[nowned+nnotowned]; // Exterior faces
    
    idx = 0; i = 0; j = 0;
    while ((ment = MSet_Next_Entry(OwnedFaces,&idx))) {
      int gid = MEnt_GlobalID(ment);

      face_gids[i++] = gid-1;
     
      if (manifold_dimension() == 3) {
        List_ptr fregs = MF_Regions((MFace_ptr) ment);
        if (List_Num_Entries(fregs) == 1)
          extface_gids[j++] = gid-1;
        if (fregs)
          List_Delete(fregs);
      }
      else if (manifold_dimension() == 2) {
        List_ptr efaces = ME_Faces((MEdge_ptr) ment);
        if (List_Num_Entries(efaces) == 1)
          extface_gids[j++] = gid-1;
        if (efaces)
          List_Delete(efaces);
      }
    }
    n_extface = j;
    nface = nowned;
    
    face_map_wo_ghosts_ = new Epetra_Map(-1,nface,face_gids,0,*comm_);
    extface_map_wo_ghosts_ = new Epetra_Map(-1,n_extface,extface_gids,0,*comm_);

   
    idx = 0;
    while ((ment = MSet_Next_Entry(NotOwnedFaces,&idx))) {
      int gid = MEnt_GlobalID(ment);
      face_gids[i++] = gid-1;
    }
    nface += nnotowned;

    face_map_w_ghosts_ = new Epetra_Map(-1,nface,face_gids,0,*comm_);

    // Build a list of global IDs of ghost faces on a boundary which may be
    // wither exterior or on-processor boundary
    
    idx = 0;
    int nnotowned_bnd = 0;
    std::vector<int> gl_id(nnotowned), pr_id(nnotowned), lc_id(nnotowned);

    while ((ment = MSet_Next_Entry(NotOwnedFaces,&idx))) {
      int gid = MEnt_GlobalID(ment);
      if (manifold_dimension() == 3) {
        List_ptr fregs = MF_Regions((MFace_ptr) ment);
        if (List_Num_Entries(fregs) == 1) {
          gl_id[nnotowned_bnd++] = gid-1;
        }                
        if (fregs)
          List_Delete(fregs);
      }
      else if (manifold_dimension() == 2) {
        List_ptr efaces = ME_Faces((MEdge_ptr) ment);
        if (List_Num_Entries(efaces) == 1)
          gl_id[nnotowned_bnd++] = gid-1;
        if (efaces)
          List_Delete(efaces);
      }
    }

    // Get the local IDs of  (lc_id) copies of owned boundary faces on remote processors (pr_id).
    // In effect we are checking if a ghost face that claims to be on the boundary is in the
    // owned boundary face list on another processor (pr_id >= 0)
    
    extface_map_wo_ghosts_->RemoteIDList(nnotowned, gl_id.data(), pr_id.data(), lc_id.data());

    int n_extface_w_ghosts = extface_map_wo_ghosts_->NumMyElements();

    // Add to maping only external faces (which belong to local mapping on other processors
    for (int k=0; k < nnotowned_bnd; k++) {
      if (pr_id[k] >= 0) {
        n_extface_w_ghosts++;
      }
    }

    std::vector<int> global_id_ghosted(n_extface_w_ghosts);
    for (int k=0; k<n_extface; k++)  {
      global_id_ghosted[k] = extface_gids[k];  
    }

    // Add to maping only external faces (which belong to local mapping on other processors
    int l=0;
    for (int k=0; k < nnotowned_bnd; k++) {
      if (pr_id[k] >= 0) {
        global_id_ghosted[n_extface + l] = gl_id[k];
        l++;
      }
    }
    
    extface_map_w_ghosts_ = new Epetra_Map(-1, n_extface_w_ghosts, global_id_ghosted.data(), 0, *comm_);
        
  }
  else {

    if (manifold_dimension() == 3) {

      nface = MESH_Num_Faces(mesh_);
      face_gids = new int[nface];
      extface_gids = new int[nface];
      
      idx = 0; i = 0; j = 0;
      while ((ment = MESH_Next_Face(mesh_,&idx))) {
        int gid = MEnt_ID(ment);
        face_gids[i++] = gid-1;
        
        List_ptr fregs = MF_Regions((MFace_ptr) ment);
        if (List_Num_Entries(fregs) == 1)
          extface_gids[j++] = gid-1;
        if (fregs)
          List_Delete(fregs);
      }
      
    }
    else if (manifold_dimension() == 2) {
      
      nface = MESH_Num_Edges(mesh_);
      face_gids = new int[nface];
      extface_gids = new int[nface];
      
      idx = 0; i = 0; j = 0;
      while ((ment = MESH_Next_Edge(mesh_,&idx))) {
        int gid = MEnt_ID(ment);
        face_gids[i++] = gid-1;
        
        List_ptr efaces = ME_Faces((MEdge_ptr) ment);
        if (List_Num_Entries(efaces) == 1)
          extface_gids[j++] = gid-1;
        if (efaces)
          List_Delete(efaces);
      }
    }
    n_extface = j;

    face_map_wo_ghosts_ = new Epetra_Map(-1,nface,face_gids,0,*comm_);
    extface_map_wo_ghosts_ = new Epetra_Map(-1,n_extface,extface_gids,0,*comm_);
  }

  owned_to_extface_importer_ = new Epetra_Import(*extface_map_wo_ghosts_,*face_map_wo_ghosts_);

  delete [] face_gids;
  delete [] extface_gids;
}


//---------------------------------------------------------
// Epetra map for edges - basically a structure specifying the global
// IDs of edges owned or used by this processor. This helps Epetra
// understand inter-partition dependencies of the data.
//
// Amanzi/Epetra want global IDs to start at 0
//---------------------------------------------------------
void Mesh_MSTK::init_edge_map()
{
  int *edge_gids;
  int nedge, idx, i;
  MEntity_ptr ment;

  if (!serial_run) {

    // For parallel runs create map without and with not-owned edges included
    // Also, put in owned edges before the not-owned edges

    int nowned = MSet_Num_Entries(OwnedEdges);
    int nnotowned = MSet_Num_Entries(NotOwnedEdges);

    edge_gids = new int[nowned+nnotowned];
    
    idx = 0; i = 0;
    while ((ment = MSet_Next_Entry(OwnedEdges,&idx))) {
      int gid = MEnt_GlobalID(ment);

      edge_gids[i++] = gid-1;
    }
    nedge = nowned;
    
    edge_map_wo_ghosts_ = new Epetra_Map(-1,nedge,edge_gids,0,*comm_);

    idx = 0;
    while ((ment = MSet_Next_Entry(NotOwnedEdges,&idx))) 
      edge_gids[i++] = MEnt_GlobalID(ment)-1;

    nedge += nnotowned;

    edge_map_w_ghosts_ = new Epetra_Map(-1,nedge,edge_gids,0,*comm_);

  }
  else {

    nedge = MESH_Num_Edges(mesh_);
    edge_gids = new int[nedge];
      
    idx = 0; i = 0;
    while ((ment = MESH_Next_Edge(mesh_,&idx))) {
      int gid = MEnt_ID(ment);
      edge_gids[i++] = gid-1;
    }

    edge_map_wo_ghosts_ = new Epetra_Map(-1,nedge,edge_gids,0,*comm_);
  }

  delete [] edge_gids;
}


//---------------------------------------------------------
// Epetra map for nodes - basically a structure specifying the global
// IDs of nodes owned or used by this processor. This helps Epetra
// understand inter-partition dependencies of the data.
//
// Amanzi/Epetra want global IDs to start at 0
//---------------------------------------------------------
void Mesh_MSTK::init_node_map()
{
  int *node_gids, *extnode_gids;
  int nnode, n_extnode, idx, i, j;
  MEntity_ptr ment;

  if (!serial_run) {

    // For parallel runs create map without and with ghost nodes included
    // Also, put in owned nodes before the ghost nodes

    int nowned = MSet_Num_Entries(OwnedVerts);
    int nnotowned = MSet_Num_Entries(NotOwnedVerts);

    node_gids = new int[nowned+nnotowned];
    extnode_gids = new int[nowned+nnotowned]; // Exterior nodes
    
    idx = 0; i = 0; j = 0;
    while ((ment = MSet_Next_Entry(OwnedVerts,&idx))) {
      int gid = MEnt_GlobalID(ment);
      node_gids[i++] = gid-1;

      if (is_boundary_node_(ment))
        extnode_gids[j++] = gid-1;
    }
    n_extnode = j;
    nnode = nowned;
    
    node_map_wo_ghosts_ = new Epetra_Map(-1,nnode,node_gids,0,*comm_);
    extnode_map_wo_ghosts_ = new Epetra_Map(-1,n_extnode,extnode_gids,0,*comm_);
    

    idx = 0;
    while ((ment = MSet_Next_Entry(NotOwnedVerts,&idx))) {
      node_gids[i++] = MEnt_GlobalID(ment)-1;
    }
    nnode += nnotowned;

    node_map_w_ghosts_ = new Epetra_Map(-1,nnode,node_gids,0,*comm_);


    // Build a list of global IDs of ghost nodes on the boundary which may be 
    // either exterior or processor boundary
    
    idx = 0;
    int nnotowned_bnd = 0;
    std::vector<int> gl_id(nnotowned), pr_id(nnotowned), lc_id(nnotowned);

    while ((ment = MSet_Next_Entry(NotOwnedVerts,&idx))) {
      int gid = MEnt_GlobalID(ment);
      if (is_boundary_node_(ment))
        gl_id[nnotowned_bnd++] = gid-1;
    }

    // Get the local IDs of (lc_id) copies of owned boundary nodes on remote processors (pr_id).
    // In effect we are checking if a ghost node that claims to be on the boundary is in the
    // owned boundary node list on another processor (pr_id >= 0)
    
    extnode_map_wo_ghosts_->RemoteIDList(nnotowned, gl_id.data(), pr_id.data(), lc_id.data());

    int n_extnode_w_ghosts = extnode_map_wo_ghosts_->NumMyElements();

    // Add to maping only external faces (which belong to local mapping on other processors
    for (int k=0; k < nnotowned_bnd; k++) {
      if (pr_id[k] >= 0) {
        n_extnode_w_ghosts++;
      }
    }

    std::vector<int> global_id_ghosted(n_extnode_w_ghosts);
    for (int k=0; k<n_extnode; k++)  {
      global_id_ghosted[k] = extnode_gids[k];  
    }

    //Add to maping only external faces (which belong to local mapping on other processors
    int l=0;
    for (int k=0; k < nnotowned_bnd; k++) {
      if (pr_id[k] >= 0) {
        global_id_ghosted[n_extnode + l] = gl_id[k];
        l++;
      }
    }

    extnode_map_w_ghosts_ = new Epetra_Map(-1, n_extnode_w_ghosts, global_id_ghosted.data(), 0, *comm_);

  } else {
    nnode = MSet_Num_Entries(OwnedVerts);

    node_gids = new int[nnode];
    extnode_gids = new int[nnode];

    idx = 0; i = 0; j = 0;
    while ((ment = MSet_Next_Entry(OwnedVerts,&idx))) {
      int gid = MEnt_ID(ment);
      node_gids[i++] = gid-1;

      if (is_boundary_node_(ment))
        extnode_gids[j++] = gid-1;
    }
    n_extnode = j;

    node_map_wo_ghosts_ = new Epetra_Map(-1,nnode,node_gids,0,*comm_);
    extnode_map_wo_ghosts_ = new Epetra_Map(-1,n_extnode,extnode_gids,0,*comm_);
  }

  delete [] node_gids;
  delete [] extnode_gids;
}


//---------------------------------------------------------
// Global ID of any entity
//---------------------------------------------------------
Entity_ID Mesh_MSTK::GID(const Entity_ID lid, const Entity_kind kind) const
{
  MEntity_ptr ent;

  switch (kind) {
  case NODE:
    ent = vtx_id_to_handle[lid];
    break;

  case EDGE:
    ent = edge_id_to_handle[lid];
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

  if (serial_run)
    return MEnt_ID(ent)-1;
  else
    return MEnt_GlobalID(ent)-1;
}


//---------------------------------------------------------
// Procedure to perform all the post-mesh creation steps in a constructor
//---------------------------------------------------------
void Mesh_MSTK::post_create_steps_(const bool request_faces, 
                                   const bool request_edges)
{
  label_celltype();

  // Initialize data structures for various entities - vertices/nodes
  // and cells are always initialized; edges and faces only if
  // requested

  init_nodes();

  edgeflip = NULL; 
  if (request_edges) init_edges();
  if (request_faces) init_faces();
  init_cells();

  if (geometric_model() != Teuchos::null)
    init_set_info();

}


//---------------------------------------------------------
// Some initializations
//---------------------------------------------------------
void Mesh_MSTK::clear_internals_() 
{ 

  faceflip = nullptr;

  cell_map_w_ghosts_ = cell_map_wo_ghosts_ = nullptr;
  edge_map_w_ghosts_ = edge_map_wo_ghosts_ = nullptr;
  face_map_w_ghosts_ = face_map_wo_ghosts_ = nullptr;
  node_map_w_ghosts_ = node_map_wo_ghosts_ = nullptr;

  mesh_ = nullptr;

  OwnedVerts = NotOwnedVerts = nullptr;
  OwnedEdges = NotOwnedEdges = nullptr;
  OwnedFaces = NotOwnedFaces = nullptr;
  OwnedCells = GhostCells = nullptr;

  celltype_att = nullptr;
  rparentatt = fparentatt = eparentatt = vparentatt = nullptr;
}


//---------------------------------------------------------
// initialize vertex info
//---------------------------------------------------------
void Mesh_MSTK::init_nodes()
{
  // create owned and not owned vertex lists

  init_pvert_lists();

  // create maps from IDs to handles

  init_vertex_id2handle_maps();

  // Create Epetra_maps indicating global IDs of owned and not owned nodes 
  
  init_node_map();
}


//---------------------------------------------------------
// Initialize edge info
//---------------------------------------------------------
void Mesh_MSTK::init_edges()
{
  edges_initialized = true;

  // Create owned and not owned lists

  init_pedge_lists();

  // Create maps from IDs to handles 

  init_edge_id2handle_maps();

  // Initialize boolean flag indicating whether slave edges are reversed in
  // direction from the master and must be flipped

  init_pedge_dirs();

  // Create epetra map containing global IDs of owned and not owned edges

  init_edge_map();
}


//---------------------------------------------------------
// Initialize face info
//---------------------------------------------------------
void Mesh_MSTK::init_faces() {

  faces_initialized = true;

  // Create owned and not owned lists

  init_pface_lists();

  // Create maps from IDs to handles 

  init_face_id2handle_maps();

  // Initialize boolean flag indicating whether slave faces are reversed in
  // direction from the master and must be flipped

  init_pface_dirs();

  // Create epetra map containing global IDs of owned and not owned faces

  init_face_map();
}


//---------------------------------------------------------
// Initialize cell info
//---------------------------------------------------------
void Mesh_MSTK::init_cells()
{
  // create owned and not owned cell lists

  init_pcell_lists(); 

  // create maps from IDs to handles

  init_cell_id2handle_maps();

  // Create Epetra_maps indicating global IDs of owned and not owned cells

  init_cell_map();
}


//---------------------------------------------------------
// ID to handle/pointer map for vertices
//---------------------------------------------------------
void Mesh_MSTK::init_vertex_id2handle_maps()
{
  int lid, nv, idx;
  MVertex_ptr vtx;

  // If the mesh is dynamic, then this code has to be revisited
  
  // Amanzi has IDs starting from 0, MSTK has IDs starting from 1

  nv = MESH_Num_Vertices(mesh_);

  vtx_id_to_handle.resize(nv);

  idx = 0; lid = 1;
  while ((vtx = MSet_Next_Entry(OwnedVerts,&idx))) {
    MEnt_Set_ID(vtx,lid);
    vtx_id_to_handle[lid-1] = vtx;
    lid++;
  }
    
  idx = 0;
  while ((vtx = MSet_Next_Entry(NotOwnedVerts,&idx))) {
    MEnt_Set_ID(vtx,lid);
    vtx_id_to_handle[lid-1] = vtx;
    lid++;
  }

}


//---------------------------------------------------------
// ID to handle/pointer map for edges
//---------------------------------------------------------
void Mesh_MSTK::init_edge_id2handle_maps()
{
  int lid, ne, idx;
  MEdge_ptr edge;

  // If the mesh is dynamic, then this code has to be revisited
  
  // Amanzi has IDs starting from 0, MSTK has IDs starting from 1

  ne = MESH_Num_Edges(mesh_);

  edge_id_to_handle.resize(ne);

  idx = 0; lid = 1;
  while ((edge = MSet_Next_Entry(OwnedEdges,&idx))) {
    MEnt_Set_ID(edge,lid);
    edge_id_to_handle[lid-1] = edge;
    lid++;
  }
    
  idx = 0;
  while ((edge = MSet_Next_Entry(NotOwnedEdges,&idx))) {
    MEnt_Set_ID(edge,lid);
    edge_id_to_handle[lid-1] = edge;
    lid++;
  }
}


//---------------------------------------------------------
// ID to handle/pointer map for faces
//---------------------------------------------------------
void Mesh_MSTK::init_face_id2handle_maps()
{
  int lid, nf, idx;
  MEntity_ptr genface;  // Mesh face in 3D, edge in 2D

  // If the mesh is dynamic, then this code has to be revisited
  
  // Amanzi has IDs starting from 0, MSTK has IDs starting from 1

  nf = (manifold_dimension() == 2) ? MESH_Num_Edges(mesh_) : MESH_Num_Faces(mesh_);

  face_id_to_handle.resize(nf);

  idx = 0; lid = 1;
  while ((genface = MSet_Next_Entry(OwnedFaces,&idx))) {
    MEnt_Set_ID(genface,lid);
    face_id_to_handle[lid-1] = genface;
    lid++;
  }
  
  idx = 0;
  while ((genface = MSet_Next_Entry(NotOwnedFaces,&idx))) {
    MEnt_Set_ID(genface,lid);
    face_id_to_handle[lid-1] = genface;
    lid++;
  }
}


//---------------------------------------------------------
// ID to handle/pointer map for cells
//---------------------------------------------------------
void Mesh_MSTK::init_cell_id2handle_maps()
{
  int lid, nc, idx;
  MEntity_ptr gencell;  // Mesh region in 3D, face in 2D

  // If the mesh is dynamic, then this code has to be revisited
  
  // Amanzi has IDs starting from 0, MSTK has IDs starting from 1

  nc = (manifold_dimension() == 2) ? MESH_Num_Faces(mesh_) : MESH_Num_Regions(mesh_);

  cell_id_to_handle.resize(nc);

  idx = 0; lid = 1;
  while ((gencell = MSet_Next_Entry(OwnedCells,&idx))) {
    MEnt_Set_ID(gencell,lid);
    cell_id_to_handle[lid-1] = gencell;
    lid++;
  }
    
  idx = 0;
  while ((gencell = MSet_Next_Entry(GhostCells,&idx))) {
    MEnt_Set_ID(gencell,lid);
    cell_id_to_handle[lid-1] = gencell;
    lid++;
  }
}


//---------------------------------------------------------
// create lists of owned and not owned vertices
//---------------------------------------------------------
void Mesh_MSTK::init_pvert_lists()
{
  int idx = 0;
  MVertex_ptr vtx;

  // Get all vertices on this processor 

  NotOwnedVerts = MSet_New(mesh_,"NotOwnedVerts",MVERTEX);
  OwnedVerts = MSet_New(mesh_,"OwnedVerts",MVERTEX);

  idx = 0;
  while ((vtx = MESH_Next_Vertex(mesh_,&idx))) {
    if (MV_PType(vtx) == PGHOST)
      MSet_Add(NotOwnedVerts,vtx);
    else
      MSet_Add(OwnedVerts,vtx);
  }
}


//---------------------------------------------------------
// create lists of owned and not owned edges
//---------------------------------------------------------
void Mesh_MSTK::init_pedge_lists()
{
  int idx = 0;
  MEdge_ptr edge;

  // Get all vertices on this processor 

  NotOwnedEdges = MSet_New(mesh_,"NotOwnedEdges",MEDGE);
  OwnedEdges = MSet_New(mesh_,"OwnedEdges",MEDGE);

  idx = 0;
  while ((edge = MESH_Next_Edge(mesh_,&idx))) {
    if (ME_PType(edge) == PGHOST)
      MSet_Add(NotOwnedEdges,edge);
    else
      MSet_Add(OwnedEdges,edge);
  }

} // Mesh_MSTK::init_pedge_lists


void Mesh_MSTK::init_pedge_dirs() {
  MEdge_ptr edge;
  MAttrib_ptr attev0, attev1;
  int idx;

  int ne = MESH_Num_Edges(mesh_);

  if (serial_run) {
    edgeflip = new bool[ne];
    for (int i = 0; i < ne; ++i) edgeflip[i] = false;
  }
  else {
    // Do some additional processing to see if ghost edges and their masters
    // are oriented the same way; if not, turn on flag to flip the directions
    // when returning to the application code

    attev0 = MAttrib_New(mesh_,"TMP_EV0_ATT",INT,MEDGE);
    attev1 = MAttrib_New(mesh_,"TMP_EV1_ATT",INT,MEDGE);
  

    idx = 0;
    while ((edge = MESH_Next_Edge(mesh_,&idx))) {
      if (ME_PType(edge) != PINTERIOR) {
        MVertex_ptr vertex0 = ME_Vertex(edge,0);
        MVertex_ptr vertex1 = ME_Vertex(edge,1);

        MEnt_Set_AttVal(edge,attev0,MEnt_GlobalID(vertex0),0.0,NULL);
        MEnt_Set_AttVal(edge,attev1,MEnt_GlobalID(vertex1),0.0,NULL);
      }
    }  


    MESH_UpdateAttributes(mesh_,mpicomm_);
    
    
    edgeflip = new bool[ne];
    for (int i = 0; i < ne; ++i) edgeflip[i] = false;
    
    double rval;
    void *pval;

    idx = 0;
    while ((edge = MSet_Next_Entry(NotOwnedEdges,&idx))) {
      int remote_vertexid0, remote_vertexid1;

      MEnt_Get_AttVal(edge,attev0,&remote_vertexid0,&rval,&pval);
      MEnt_Get_AttVal(edge,attev1,&remote_vertexid1,&rval,&pval);
      
      int local_vertexid0 = MEnt_GlobalID(ME_Vertex(edge,0));
      int local_vertexid1 = MEnt_GlobalID(ME_Vertex(edge,1));
      
      if (remote_vertexid1 == local_vertexid0 || 
          remote_vertexid0 == local_vertexid1) {
        int lid = MEnt_ID(edge);
        edgeflip[lid-1] = true;
      }
      else { // Sanity Check

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
void Mesh_MSTK::init_pface_lists()
{
  int idx = 0;

  // Get all faces on this processor 

  if (manifold_dimension() == 3) {

    MFace_ptr face;

    NotOwnedFaces = MSet_New(mesh_,"NotOwnedFaces",MFACE);
    OwnedFaces = MSet_New(mesh_,"OwnedFaces",MFACE);

    idx = 0;
    while ((face = MESH_Next_Face(mesh_,&idx))) {
      if (MF_PType(face) == PGHOST)
        MSet_Add(NotOwnedFaces,face);
      else
        MSet_Add(OwnedFaces,face);
    }
  }
  else if (manifold_dimension() == 2) {

    MEdge_ptr edge;

    NotOwnedFaces = MSet_New(mesh_,"NotOwnedFaces",MFACE);
    OwnedFaces = MSet_New(mesh_,"OwnedFaces",MFACE);

    idx = 0;
    while ((edge = MESH_Next_Edge(mesh_,&idx))) {
      if (ME_PType(edge) == PGHOST)
        MSet_Add(NotOwnedFaces,edge);
      else
        MSet_Add(OwnedFaces,edge);
    }
  }
  else {
    std::cerr << "Not implemented for face dimension" << std::endl;
  }

  return;
}

// Detect whether ghost faces are in opposite direction of owned faces
// on processor boundaries

void Mesh_MSTK::init_pface_dirs()
{
  int nf = (manifold_dimension() == 2) ? MESH_Num_Edges(mesh_) : MESH_Num_Faces(mesh_);

  if (serial_run) {
    faceflip = new bool[nf];
    for (int i = 0; i < nf; ++i) faceflip[i] = false;
  }
  else {
    if (manifold_dimension() == 3)
      init_pface_dirs_3();
    else if (manifold_dimension() == 2)
      init_pface_dirs_2();
  }
}


// Detect whether ghost faces are in opposite direction of owned faces
// on processor boundaries - Version for solid meshes

void Mesh_MSTK::init_pface_dirs_3() {
  MRegion_ptr region0, region1;
  MFace_ptr face;
  MAttrib_ptr attfc0, attfc1;
  int idx;
  int local_regid0, local_regid1;
  int remote_regid0, remote_regid1;
  
  int nf = MESH_Num_Faces(mesh_);

  // Do some additional processing to see if ghost faces and their masters
  // are oriented the same way; if not, turn on flag to flip the directions
  // when returning to the application code
  
  // attributes to store 
  attfc0 = MAttrib_New(mesh_,"TMP_FC0_ATT",INT,MFACE);
  attfc1 = MAttrib_New(mesh_,"TMP_FC1_ATT",INT,MFACE);

  idx = 0;
  while ((face = MESH_Next_Face(mesh_,&idx))) {
    if (MF_PType(face) != PINTERIOR) {
      region0 = MF_Region(face,0);
      if (region0)
        MEnt_Set_AttVal(face,attfc0,MEnt_GlobalID(region0),0.0,NULL);
      
      region1 = MF_Region(face,1);
      if (region1)
        MEnt_Set_AttVal(face,attfc1,MEnt_GlobalID(region1),0.0,NULL);
    }
  }    
  
  MESH_UpdateAttributes(mesh_,mpicomm_);


  faceflip = new bool[nf];
  for (int i = 0; i < nf; ++i) faceflip[i] = false;
  
  double rval;
  void *pval;
  
  idx = 0;
  while ((face = MSet_Next_Entry(NotOwnedFaces,&idx))) {
    
    MEnt_Get_AttVal(face,attfc0,&remote_regid0,&rval,&pval);
    MEnt_Get_AttVal(face,attfc1,&remote_regid1,&rval,&pval);
    
    region0 = MF_Region(face,0);
    local_regid0 = region0 ? MEnt_GlobalID(region0) : 0;
    region1 = MF_Region(face,1);
    local_regid1 = region1 ? MEnt_GlobalID(region1) : 0;
    
    if (remote_regid1 == local_regid0 || 
        remote_regid0 == local_regid1) {
      int lid = MEnt_ID(face);
      faceflip[lid-1] = true;
    }
    else { // Sanity Check
      
      if (remote_regid1 != local_regid1 &&
          remote_regid0 != local_regid0) {
        
        std::stringstream mesg_stream;
        mesg_stream << "Face cells mismatch between master and ghost (processor " << myprocid << ")";
        Errors::Message mesg(mesg_stream.str());
        Exceptions::amanzi_throw(mesg);
      }
    }
  }
  MAttrib_Delete(attfc0);
  MAttrib_Delete(attfc1);
}


// Detect whether ghost faces are in opposite direction of owned faces
// on processor boundaries - Version for surface meshes

void Mesh_MSTK::init_pface_dirs_2() {
  MEdge_ptr edge;
  MAttrib_ptr attev0, attev1;
  int idx;

  int ne = MESH_Num_Edges(mesh_);

  // Do some additional processing to see if ghost faces and their masters
  // are oriented the same way; if not, turn on flag to flip the directions
  // when returning to the application code
  
  attev0 = MAttrib_New(mesh_,"TMP_EV0_ATT",INT,MEDGE);
  attev1 = MAttrib_New(mesh_,"TMP_EV1_ATT",INT,MEDGE);

  idx = 0;
  while ((edge = MESH_Next_Edge(mesh_,&idx))) {
    if (ME_PType(edge) != PINTERIOR) {
      MVertex_ptr ev0 = ME_Vertex(edge, 0);
      MVertex_ptr ev1 = ME_Vertex(edge, 1);
      
      MEnt_Set_AttVal(edge,attev0,MEnt_GlobalID(ev0),0.0,NULL);
      MEnt_Set_AttVal(edge,attev1,MEnt_GlobalID(ev1),0.0,NULL);
    }
  }

  MESH_UpdateAttributes(mesh_,mpicomm_);

  faceflip = new bool[ne];
  for (int i = 0; i < ne; ++i) faceflip[i] = false;
    
  double rval;
  void *pval;
  
  idx = 0;
  while ((edge = MSet_Next_Entry(NotOwnedFaces,&idx))) {
    int remote_evgid0, remote_evgid1;
    MEnt_Get_AttVal(edge,attev0,&remote_evgid0,&rval,&pval);
    MEnt_Get_AttVal(edge,attev1,&remote_evgid1,&rval,&pval);
    
    MVertex_ptr ev0 = ME_Vertex(edge, 0);
    MVertex_ptr ev1 = ME_Vertex(edge, 1);
    int local_evgid0 = MV_GlobalID(ev0);
    int local_evgid1 = MV_GlobalID(ev1);
    
    if (remote_evgid1 == local_evgid0 || 
        remote_evgid0 == local_evgid1) {
      int lid = MEnt_ID(edge);
      faceflip[lid-1] = true;
    }
  }
  MAttrib_Delete(attev0);
  MAttrib_Delete(attev1);
}


//---------------------------------------------------------
// create lists of owned and not owned cells
//---------------------------------------------------------
void Mesh_MSTK::init_pcell_lists()
{
  int idx = 0;

  if (manifold_dimension() == 3) {
    MRegion_ptr region;

    OwnedCells = MSet_New(mesh_,"OwnedCells",MREGION);
    GhostCells = MSet_New(mesh_,"GhostCells",MREGION);

    idx = 0;
    while ((region = MESH_Next_Region(mesh_,&idx))) {
      if (MR_PType(region) == PGHOST)
        MSet_Add(GhostCells,region);
      else
        MSet_Add(OwnedCells,region);
    }
  }
  else if (manifold_dimension() == 2) {
    MFace_ptr face;

    OwnedCells = MSet_New(mesh_,"OwnedCells",MFACE);
    GhostCells = MSet_New(mesh_,"GhostCells",MFACE);

    idx = 0;
    while ((face = MESH_Next_Face(mesh_,&idx))) {
      if (MF_PType(face) == PGHOST)
        MSet_Add(GhostCells,face);
      else
        MSet_Add(OwnedCells,face);
    }
  }
  else {
    Errors::Message mesg("Implemented only for 2D and 3D");
    Exceptions::amanzi_throw(mesg);
  }

  return;
}


void Mesh_MSTK::init_set_info()
{
  MSet_ptr mset;
  
  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm = geometric_model();

  if (gm == Teuchos::null) { 
    Errors::Message mesg("Need region definitions to initialize sets");
    Exceptions::amanzi_throw(mesg);
  }
    

  unsigned int ngr = gm->size();

  for (int i = 0; i < ngr; ++i) {
    Teuchos::RCP<const AmanziGeometry::Region> rgn = gm->FindRegion(i);

    MType entdim;
    if (rgn->type() == AmanziGeometry::LABELEDSET) {

      Teuchos::RCP<const AmanziGeometry::RegionLabeledSet> lsrgn =
          Teuchos::rcp_static_cast<const AmanziGeometry::RegionLabeledSet>(rgn);

      std::string internal_name;
      std::string label = lsrgn->label();
      std::string entity_type_str = lsrgn->entity_str();

      if (entity_type_str == "CELL")
        internal_name = internal_name_of_set(rgn,CELL);
      else if (entity_type_str == "FACE")
        internal_name = internal_name_of_set(rgn,FACE);
      else if (entity_type_str == "NODE")
        internal_name = internal_name_of_set(rgn,NODE);

      mset = MESH_MSetByName(mesh_,internal_name.c_str());
   
      if (!mset) {
        continue;  // Its possible some sets won't exist on some partitions

        //        Errors::Message mesg("Missing labeled set \"" + label + "\" or error in input");
        //        Exceptions::amanzi_throw(mesg);
      }

      entdim = MSet_EntDim(mset);
      if (manifold_dimension() == 3) {

        if ((entity_type_str == "CELL" && entdim != MREGION) ||
            (entity_type_str == "FACE" && entdim != MFACE) ||
            (entity_type_str == "NODE" && entdim != MVERTEX)) {
          Errors::Message mesg("Mismatch of entity type in labeled set region and mesh set");
          Exceptions::amanzi_throw(mesg);
        }
      }
      else if (manifold_dimension() == 2) {
        if ((entity_type_str == "CELL" && entdim != MFACE) ||
            (entity_type_str == "FACE" && entdim != MEDGE) ||
            (entity_type_str == "NODE" && entdim != MVERTEX)) {
          std::cerr << "Mismatch of entity type in labeled set region and mesh set" << std::endl;
          throw std::exception();
        }
      }

      if (mset) {
        if (entities_deleted) {
          int idx = 0;
          MEntity_ptr ent;
          while ((ent = MSet_Next_Entry(mset,&idx))) {
            if (MEnt_Dim(ent) == MDELETED)
              MSet_Rem(mset, ent);
          }
        }
      }
    }
    else { /* General region - we have to account for all kinds of
              entities being queried in a set defined by this 
              region */
      Entity_kind int_to_kind[4] = {NODE,EDGE,FACE,CELL};

      for (int k = 0; k < 4; ++k) {
        Entity_kind kind = int_to_kind[k];
        
        std::string internal_name = internal_name_of_set(rgn,kind);

        mset = MESH_MSetByName(mesh_,internal_name.c_str());

        if (mset) {
          if (entities_deleted) {
            int idx = 0;
            MEntity_ptr ent;
            while ((ent = MSet_Next_Entry(mset,&idx))) {
              if (MEnt_Dim(ent) == MDELETED)
                MSet_Rem(mset, ent);
            }
          }
        }
      }
    }
  }
}


void Mesh_MSTK::collapse_degen_edges()
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
#ifdef MSTK_3_00_OR_NEWER
      /* Degenerate edge  - must collapse */

      entities_deleted = true;

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
#else
      Errors::Message mesg("Mesh contains a degenerate edge. MSTK version 3.0.0 or later is required to support collapsing of degenerate edges!");
      Exceptions::amanzi_throw(mesg);
#endif
    }
  }
#ifdef MSTK_3_00_OR_NEWER
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

#endif

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


Cell_type Mesh_MSTK::MFace_Celltype(MFace_ptr face)
{
  int nfv = MF_Num_Vertices(face);

  switch (nfv) {
  case 3:
    return TRI;
    break;
  case 4:
    return QUAD;
    break;
  default:
    return POLYGON;
  }
}


Cell_type Mesh_MSTK::MRegion_Celltype(MRegion_ptr region)
{
  List_ptr rverts, rfaces;
  MFace_ptr face;
  int nrv, nrf, idx2, nquads;

  rverts = MR_Vertices(region);
  nrv = List_Num_Entries(rverts);
  List_Delete(rverts);
  
  nrf = MR_Num_Faces(region);
  
  switch (nrf) {
  case 4:
    if (nrv == 4)
      return TET;
    else
      return POLYHED;
    break;
  case 5:
    
    nquads = 0;
    rfaces = MR_Faces(region);          
    idx2 = 0;
    while ((face = List_Next_Entry(rfaces,&idx2)))
      if (MF_Num_Vertices(face) == 4)
        nquads++;
    List_Delete(rfaces);
    
    switch (nquads) {
    case 1:
      return PYRAMID;
      break;
    case 3:
      return PRISM;
      break;
    default:
      return POLYHED;
    }

    break;

  case 6:

    nquads = 0;
    rfaces = MR_Faces(region);          
    idx2 = 0;
    while ((face = List_Next_Entry(rfaces,&idx2)))
      if (MF_Num_Vertices(face) == 4)
        nquads++;
    List_Delete(rfaces);
    
    if (nquads == 6) 
      return HEX;
    else 
      return POLYHED;

    break;

  default:
    return POLYHED;
  }
}


void Mesh_MSTK::label_celltype()
{
  Cell_type ctype;
  int idx;
  MFace_ptr face;
  MRegion_ptr region;

  if (manifold_dimension() == 2) 
    celltype_att = MAttrib_New(mesh_,"Cell_type",INT,MFACE);
  else
    celltype_att = MAttrib_New(mesh_,"Cell_type",INT,MREGION);

  if (manifold_dimension() == 2) {

    idx = 0;
    while ((face = MESH_Next_Face(mesh_,&idx))) {
      ctype = MFace_Celltype(face);
      MEnt_Set_AttVal(face,celltype_att,ctype,0.0,NULL);
    }
      
  }
  else if (manifold_dimension() == 3) {

    idx = 0;
    while ((region = MESH_Next_Region(mesh_,&idx))) {
      ctype = MRegion_Celltype(region);
      MEnt_Set_AttVal(region,celltype_att,ctype,0.0,NULL);
    }
  }
}


int Mesh_MSTK::generate_regular_mesh(Mesh_ptr mesh, double x0, double y0, 
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


int Mesh_MSTK::generate_regular_mesh(Mesh_ptr mesh, double x0, double y0, 
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

  edges_initialized = false;
  faces_initialized = false;
  OwnedVerts = NotOwnedVerts = NULL;
  OwnedEdges = NotOwnedEdges = NULL;
  OwnedFaces = NotOwnedFaces = NULL;
  OwnedCells = GhostCells = NULL;
  node_map_w_ghosts_ = node_map_wo_ghosts_ = NULL;
  edge_map_w_ghosts_ = edge_map_wo_ghosts_ = NULL;
  face_map_w_ghosts_ = face_map_wo_ghosts_ = NULL;
  cell_map_w_ghosts_ = cell_map_wo_ghosts_ = NULL;
  deleted_vertices = deleted_edges = deleted_faces = deleted_regions = NULL;
  entities_deleted = false;
}


void Mesh_MSTK::inherit_labeled_sets(MAttrib_ptr copyatt,
                                     List_ptr src_entities)
{
  int idx, idx2, diffdim;
  MSet_ptr mset;
  
  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm = geometric_model();

  if (gm == Teuchos::null) { 
    std::cerr << "Need region definitions to initialize sets" << std::endl;
    return;
  }

  Mesh_ptr parent_mstk_mesh = parent_mesh_->mesh_;

  // Difference in cell dimension of this mesh and its parent
  // Labeled set entity dimensions will be similarly dialed down
  
  diffdim = parent_mesh_->manifold_dimension() - manifold_dimension();
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
        internal_name = internal_name_of_set(rgn,CELL);
      else if (lsrgn->entity_str() == "FACE")
        internal_name = internal_name_of_set(rgn,FACE);
      else if (lsrgn->entity_str() == "NODE")
        internal_name = internal_name_of_set(rgn,NODE);


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
// Write mesh out to exodus file
//---------------------------------------------------------
void
Mesh_MSTK::write_to_exodus_file(const std::string filename) const {
  MESH_ExportToExodusII(mesh_,filename.c_str(),-1,NULL,NULL,mpicomm_);
}


// Run MSTK's internal checks - meant for debugging only
// Returns true if everything is ok, false otherwise

bool
Mesh_MSTK::run_internal_mstk_checks() const {
  return MESH_CheckTopo(mesh_) && MESH_Parallel_Check(mesh_, mpicomm_);
}


//---------------------------------------------------------
// Check if node is a boundary (physical or on-processor) node.
//---------------------------------------------------------
bool
Mesh_MSTK::is_boundary_node_(const MEntity_ptr ment) const
{
  if (manifold_dimension() == 3) {
    List_ptr vfaces = MV_Faces((MVertex_ptr) ment);
    if (vfaces) {
      int nvfaces = List_Num_Entries(vfaces);
      for (int k = 0; k < nvfaces; ++k) {
        List_ptr fregs = MF_Regions((MFace_ptr) List_Entry(vfaces, k));
        if (fregs) {
          int nfregs = List_Num_Entries(fregs);
          List_Delete(fregs);
          if (nfregs == 1) {
            List_Delete(vfaces);
            return true;
          }
        }
      }
      List_Delete(vfaces);
    }
  }

  else if (manifold_dimension() == 2) {
    List_ptr vedges = MV_Edges((MVertex_ptr) ment);
    if (vedges) {
      int nvedges = List_Num_Entries(vedges);
      for (int k = 0; k < nvedges; ++k) {
        List_ptr efaces = ME_Faces((MEdge_ptr) List_Entry(vedges, k));
        if (efaces) {
          int nefaces = List_Num_Entries(efaces);
          List_Delete(efaces);
          if (nefaces == 1) {
            List_Delete(vedges);
            return true;
          }
        }
      }
      List_Delete(vedges);
    }
  }

  return false;
}

}  // namespace AmanziMesh
}  // namespace Amanzi


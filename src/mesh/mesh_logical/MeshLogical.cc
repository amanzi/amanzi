#include "EnumeratedSetRegion.hh"
#include "MeshLogical.hh"

namespace Amanzi {
namespace AmanziMesh {

//
// MeshLogical Constructor
//  Includes gravity.
//
//  - cell_volume             : length ncells array of cell volumes
//  - face_cell_list          : length nfaces array of length 2 arrays
//                              defining the topology
//  - face_cell_lengths       : length of the cell-to-face connection
//  - face_area_normals       : length nfaces array of normals of the
//                              face, points from cell 1 to 2 in
//                              face_cell_list topology, magnitude
//                              is area
//  - cell_centroids          : (optional, for plotting) length ncell
//                              array of centroids
//
// Breaks standards following the rest of the mesh infrastructure.
//
MeshLogical::MeshLogical(const Epetra_MpiComm *incomm,
			 const std::vector<double>& cell_volumes_,
			 const std::vector<Entity_ID_List>& face_cell_ids_,
			 const std::vector<std::vector<double> >& face_cell_lengths_,
			 const std::vector<AmanziGeometry::Point>& face_area_normals_,
			 const std::vector<AmanziGeometry::Point>* cell_centroids_,
			 const VerboseObject *verbosity_obj)
  : Mesh(verbosity_obj) {

  ASSERT(face_cell_ids_.size() == face_cell_lengths_.size());
  ASSERT(face_cell_ids_.size() == face_area_normals_.size());
  
  set_comm(incomm);
  cell_volumes = cell_volumes_;
  face_cell_ids = face_cell_ids_;
  
  // normal1 is negative normal0
  face_normal0 = face_area_normals_;
  face_normal1 = face_area_normals_;
  for (int f=0; f!=face_area_normals_.size(); ++f) {
    if (cell_centroids_) {
      if (face_cell_ids[f].size() == 2) {
	if (face_normal0[f] * ((*cell_centroids_)[face_cell_ids[f][1]] -
			       (*cell_centroids_)[face_cell_ids[f][0]]) > 0.) {
	  // normal is outward from cell 0 to cell 1
	  face_normal1[f] = -face_normal1[f];
	} else {
	  // normal is outward from cell 1 to cell 0
	  face_normal0[f] = -face_normal0[f];
	}
      } else {
	// pass, boundary face and we must assume normal given was correct
      }
    } else {
      // no centroid info, doesn't matter what we choose
      face_normal1[f] = -face_normal1[f];
    }
  }

  // optional centroids
  if (cell_centroids_) {
    cell_centroids = *cell_centroids_;

    // face centroids constructed from normals, lengths, cell
    // centroids using assumption of perpendicular bisector
    face_centroids.resize(face_cell_ids.size());
    for (int f=0; f!=face_cell_ids.size(); ++f) {
      if (face_cell_ids[f].size() == 2) {
	face_centroids[f] = (cell_centroids[face_cell_ids[f][0]] +
			     cell_centroids[face_cell_ids[f][1]]) / 2.0;
      } else {
	ASSERT(face_cell_ids[f].size() == 1);
	face_centroids[f] = cell_centroids[face_cell_ids[f][0]]
	  + (face_cell_lengths_[f][0] / AmanziGeometry::norm(face_normal0[f]))
	  * face_normal0[f];
      }
    }
  }

  // populate cell, extra face info
  face_cell_ptype.resize(face_cell_ids.size());
  face_areas.resize(face_cell_ids.size());
  
  cell_face_ids.resize(cell_volumes.size());
  cell_face_dirs.resize(cell_volumes.size());
  cell_face_bisectors.resize(cell_volumes.size());
  int f_id=0;
  for (std::vector<Entity_ID_List>::const_iterator f=face_cell_ids.begin();
       f!=face_cell_ids.end(); ++f) {
    face_cell_ptype[f_id].push_back(OWNED);
    face_cell_ptype[f_id].push_back(f->size() == 2 ?
				    OWNED : PTYPE_UNKNOWN);
    face_areas[f_id] = AmanziGeometry::norm(face_normal0[f_id]);
    
    cell_face_ids[(*f)[0]].push_back(f_id);
    cell_face_dirs[(*f)[0]].push_back(1);

    AmanziGeometry::Point unit_normal(face_normal0[f_id]);
    unit_normal /= AmanziGeometry::norm(unit_normal);
    unit_normal *= face_cell_lengths_[f_id][0];
    cell_face_bisectors[(*f)[0]].push_back(unit_normal);

    if (f->size() > 1 && (*f)[1] >= 0) {
      cell_face_ids[(*f)[1]].push_back(f_id);
      cell_face_dirs[(*f)[1]].push_back(-1);

      AmanziGeometry::Point unit_normal(face_normal1[f_id]);
      unit_normal /= face_areas[f_id];
      cell_face_bisectors[(*f)[0]].push_back(unit_normal * face_cell_lengths_[f_id][1]);
    }      
    
    f_id++;
  }

  // toggle flags
  cell_geometry_precomputed = true;
  face_geometry_precomputed = true;
  cell2face_info_cached = true;
  faces_requested = true;
  face2cell_info_cached = true;

  // build epetra maps
  init_maps();
}  

  
// build maps
void
MeshLogical::init_maps() {
  // cell map
  maps_[CELL] = Teuchos::rcp(new Epetra_Map(-1, cell_face_ids.size(), 0, *comm));

  // face map
  Teuchos::RCP<Epetra_Map> face_map =
    Teuchos::rcp(new Epetra_Map(-1, face_cell_ids.size(), 0, *comm));
  maps_[FACE] = face_map;

  // exterior face map
  std::vector<int> extface_ids;
  int f_id = 0;
  for (std::vector<Entity_ID_List>::iterator f=face_cell_ids.begin();
       f!=face_cell_ids.end(); ++f) {
    if (f->size() == 1) {
      extface_ids.push_back(face_map->GID(f_id));
    }
    f_id++;
  }
  maps_[BOUNDARY_FACE] = Teuchos::rcp(new Epetra_Map(-1, extface_ids.size(),
						     &extface_ids[0], 0, *comm));

  exterior_face_importer_ = Teuchos::rcp(new Epetra_Import(*maps_[BOUNDARY_FACE],
							   *face_map));  


  num_entities_[CELL] = maps_[CELL]->NumMyElements();
  num_entities_[FACE] = maps_[FACE]->NumMyElements();
  num_entities_[BOUNDARY_FACE] = maps_[BOUNDARY_FACE]->NumMyElements();
  num_entities_[NODE] = 0;
}


// testing purposes -- checks if the caches match
bool
MeshLogical::operator==(const MeshLogical& other) {
  double _eps = 1.e-10;
  
  if (&other == this) return true;
  if (cell_face_ids != other.cell_face_ids) return false;
  if (face_cell_ids != other.face_cell_ids) return false;

  if (cell_volumes.size() != other.cell_volumes.size()) return false;
  for (size_t i=0; i!=cell_volumes.size(); ++i) {
    if (std::abs(cell_volumes[i] - other.cell_volumes[i]) > _eps) return false;
  }

  if (face_normal0.size() != other.face_normal0.size()) return false;
  for (size_t i=0; i!=face_normal0.size(); ++i) {
    if (AmanziGeometry::norm(face_normal0[i] - other.face_normal0[i]) > _eps) return false;
  }

  if (face_normal1.size() != other.face_normal1.size()) return false;  
  for (size_t i=0; i!=face_normal1.size(); ++i) {
    if (AmanziGeometry::norm(face_normal1[i] - other.face_normal1[i]) > _eps) return false;
  }

  if (cell_centroids.size() != other.cell_centroids.size()) return false;  
  for (size_t i=0; i!=cell_centroids.size(); ++i) {
    if (AmanziGeometry::norm(cell_centroids[i] - other.cell_centroids[i]) > _eps) return false;
  }

  if (face_centroids.size() != other.face_centroids.size()) return false;  
  for (size_t i=0; i!=face_centroids.size(); ++i) {
    if (AmanziGeometry::norm(face_centroids[i] - other.face_centroids[i]) > _eps) return false;
  }
  return true;
}

  
// Get parallel type of entity - OWNED, GHOST, USED (See MeshDefs.hh)
Parallel_type
MeshLogical::entity_get_ptype(const Entity_kind kind,
			      const Entity_ID entid) const {
  return OWNED;
}

// Parent entity in the source mesh if mesh was derived from another mesh
Entity_ID
MeshLogical::entity_get_parent(const Entity_kind kind, const Entity_ID entid) const {
  return -1;
}

// Get cell type - UNKNOWN, TRI, QUAD, POLYGON, TET, PRISM, PYRAMID, HEX, POLYHED 
// See MeshDefs.hh
Cell_type
MeshLogical::cell_get_type(const Entity_ID cellid) const {
  return CELLTYPE_UNKNOWN;
}


//
// General mesh information
// -------------------------
//
// Number of entities of any kind (cell, face, node) and in a
// particular category (OWNED, GHOST, USED)
unsigned int
MeshLogical::num_entities (const Entity_kind kind,
			   const Parallel_type ptype) const {
  return num_entities_.at(kind);
}

// Global ID of any entity
Entity_ID
MeshLogical::GID(const Entity_ID lid, const Entity_kind kind) const {
  return maps_.at(kind)->GID(lid);
}


//
// Mesh Entity Adjacencies
//-------------------------
// Downward Adjacencies
//---------------------
// Get nodes of a cell
void
MeshLogical::cell_get_nodes (const Entity_ID cellid, 
			     Entity_ID_List *nodeids) const {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}


// Get the bisectors, i.e. vectors from cell centroid to face centroids.
void
MeshLogical::cell_get_faces_and_bisectors (const Entity_ID cellid,
			   Entity_ID_List *faceids,
			   std::vector<AmanziGeometry::Point> *bisectors,
			   const bool ordered) const {
  if (faceids) *faceids = cell_face_ids[cellid];
  if (bisectors) *bisectors = cell_face_bisectors[cellid];
}
  

// Get nodes of face
// On a distributed mesh, all nodes (OWNED or GHOST) of the face
// are returned
// In 3D, the nodes of the face are returned in ccw order consistent
// with the face normal
// In 2D, nfnodes is 2
void
MeshLogical::face_get_nodes (const Entity_ID faceid,
			     Entity_ID_List *nodeids) const {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}


// Get nodes of edge
void
MeshLogical::edge_get_nodes (const Entity_ID edgeid, 
			     Entity_ID *nodeid0, Entity_ID *nodeid1) const {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}

// Upward adjacencies
//-------------------

// Cells of type 'ptype' connected to a node - The order of cells is
// not guaranteed to be the same for corresponding nodes on
// different processors
void
MeshLogical::node_get_cells (const Entity_ID nodeid,
			     const Parallel_type ptype,
			     Entity_ID_List *cellids) const {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}  

// Faces of type 'ptype' connected to a node - The order of faces is
// not guarnateed to be the same for corresponding nodes on
// different processors
void
MeshLogical::node_get_faces (const Entity_ID nodeid,
			     const Parallel_type ptype,
			     Entity_ID_List *faceids) const {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}  

// Get faces of ptype of a particular cell that are connected to the
// given node - The order of faces is not guarnateed to be the same
// for corresponding nodes on different processors
void
MeshLogical::node_get_cell_faces (const Entity_ID nodeid,
				  const Entity_ID cellid,
				  const Parallel_type ptype,
				  Entity_ID_List *faceids) const {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}  

// Same level adjacencies
//-----------------------

// Face connected neighboring cells of given cell of a particular ptype
// (e.g. a hex has 6 face neighbors)

// The order in which the cellids are returned cannot be
// guaranteed in general except when ptype = USED, in which case
// the cellids will correcpond to cells across the respective
// faces given by cell_get_faces
void
MeshLogical::cell_get_face_adj_cells(const Entity_ID cellid,
				     const Parallel_type ptype,
				     Entity_ID_List *fadj_cellids) const {
  Errors::Message mesg("Not implemented."); // this could be implemented, just not sure we need it yet. --etc
  Exceptions::amanzi_throw(mesg);
}  

// Node connected neighboring cells of given cell
// (a hex in a structured mesh has 26 node connected neighbors)
// The cells are returned in no particular order
void
MeshLogical::cell_get_node_adj_cells(const Entity_ID cellid,
				     const Parallel_type ptype,
				     Entity_ID_List *nadj_cellids) const {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}

//
// Mesh entity geometry
//--------------
//
// Node coordinates - 3 in 3D and 2 in 2D
void
MeshLogical::node_get_coordinates (const Entity_ID nodeid,
				   AmanziGeometry::Point *ncoord) const {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}  


// Face coordinates - conventions same as face_to_nodes call
// Number of nodes is the vector size divided by number of spatial dimensions
void
MeshLogical::face_get_coordinates (const Entity_ID faceid,
				   std::vector<AmanziGeometry::Point> *fcoords) const {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}

// Coordinates of cells in standard order (Exodus II convention)
// STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES IN 3D
// For a general polyhedron this will return the node coordinates in
// arbitrary order
// Number of nodes is vector size divided by number of spatial dimensions
void
MeshLogical::cell_get_coordinates (const Entity_ID cellid,
				   std::vector<AmanziGeometry::Point> *ccoords) const {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}

//
// Mesh modification
//-------------------

// Set coordinates of node
void
MeshLogical::node_set_coordinates (const Entity_ID nodeid,
				   const AmanziGeometry::Point ncoord) {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}  

void
MeshLogical::node_set_coordinates (const Entity_ID nodeid,
				   const double *ncoord) {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}


// Deformation not supported.
int
MeshLogical::deform (const std::vector<double>& target_cell_volumes_in,
		     const std::vector<double>& min_cell_volumes_in,
		     const Entity_ID_List& fixed_nodes,
		     const bool move_vertical) {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
  return -1;
}  

//
// Epetra maps
//------------
const Epetra_Map&
MeshLogical::cell_map (const bool include_ghost) const {
  return *maps_.at(CELL);
}

const Epetra_Map&
MeshLogical::face_map (const bool include_ghost) const {
  return *maps_.at(FACE);
}  

const Epetra_Map&
MeshLogical::node_map (const bool include_ghost) const {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
  return *maps_.at(NODE);
}


const Epetra_Map&
MeshLogical::exterior_face_map (void) const {
  return *maps_.at(BOUNDARY_FACE);
}

// Epetra importer that will allow apps to import values from a
// Epetra vector defined on all owned faces into an Epetra vector
// defined only on exterior faces
const Epetra_Import&
MeshLogical::exterior_face_importer (void) const {
  return *exterior_face_importer_;
}


//
// Mesh Sets for ICs, BCs, Material Properties and whatever else
//--------------------------------------------------------------
//
// Get number of entities of type 'category' in set
unsigned int
MeshLogical::get_set_size (const Set_ID setid,
			   const Entity_kind kind,
			   const Parallel_type ptype) const {
  Entity_ID_List ents;
  get_set_entities(setid, kind, ptype, &ents);
  return ents.size();
}

unsigned int
MeshLogical::get_set_size (const std::string setname,
			   const Entity_kind kind,
			   const Parallel_type ptype) const {
  return get_set_size(geometric_model_->FindRegion(setname)->id(),kind,ptype);
}

unsigned int
MeshLogical::get_set_size (const char *setname,
			   const Entity_kind kind,
			   const Parallel_type ptype) const {
  std::string name(setname);
  return get_set_size(name,kind,ptype);
}

// Get list of entities of type 'category' in set
void
MeshLogical::get_set_entities (const Set_ID setid,
			       const Entity_kind kind,
			       const Parallel_type ptype,
			       Entity_ID_List *entids) const {
  AmanziGeometry::RegionPtr rgn = geometric_model_->FindRegion(setid);

  if (rgn->name() == "All" || rgn->name() == "all" || rgn->name() == "ALL") {
    int nent = num_entities(kind, ptype);
    entids->resize(num_entities(kind, ptype));
    for (int i=0; i!=nent; ++i) {
      (*entids)[i] = i;
    }
    return;
  }

  if (rgn->type() == AmanziGeometry::ENUMERATEDSET) {
    AmanziGeometry::EnumeratedSetRegionPtr esrgn =
      dynamic_cast<AmanziGeometry::EnumeratedSetRegionPtr>(rgn);

    if ((esrgn->entity_str() == "CELL" && kind == CELL) ||
	(esrgn->entity_str() == "FACE" && kind == FACE)) {
      *entids = esrgn->entities();
    }
  }
  return;
}

void
MeshLogical::get_set_entities (const std::string setname,
			       const Entity_kind kind,
			       const Parallel_type ptype,
			       Entity_ID_List *entids) const {
  get_set_entities(geometric_model_->FindRegion(setname)->id(), kind, ptype, entids);
  return;
}
  
  
void
MeshLogical::get_set_entities (const char *setname,
			       const Entity_kind kind,
			       const Parallel_type ptype,
			       Entity_ID_List *entids) const {
  std::string name(setname);
  get_set_entities(name, kind, ptype, entids);
  return;
}

// Miscellaneous functions
void
MeshLogical::write_to_exodus_file(const std::string filename) const {
  // don't know how to do this! FIXME! --etc
}


// Geometry
int
MeshLogical::compute_cell_geometry(const Entity_ID cellid, 
				   double *volume, 
				   AmanziGeometry::Point *centroid) const {
  // this is a placeholder, these cannot be recomputed
  if (volume) *volume = cell_volumes[cellid];

  if (centroid) {
    if (cell_centroids.size() > 0) {
      *centroid = cell_centroids[cellid];
    } else {
      *centroid = AmanziGeometry::Point();
    }
  }
  return 1;
}


int
MeshLogical::compute_face_geometry(const Entity_ID faceid, 
				   double *area, 
				   AmanziGeometry::Point *centroid, 
				   AmanziGeometry::Point *normal0,
				   AmanziGeometry::Point *normal1) const {
  // this is a placeholder, these cannot be recomputed
  if (area) *area = face_areas[faceid];
  if (centroid) *centroid = AmanziGeometry::Point();
  if (normal0) *normal0 = face_normal0[faceid];
  if (normal1) *normal1 = face_normal1[faceid];
  return 1;
}
  
  
// get faces of a cell and directions in which it is used - this function
// is implemented in each mesh framework. The results are cached in 
// the base class
void
MeshLogical::cell_get_faces_and_dirs_internal (const Entity_ID cellid,
					       Entity_ID_List *faceids,
					       std::vector<int> *face_dirs,
					       const bool ordered) const {
  Errors::Message mesg("DEVELOPER ERROR: cell_get_faces_and_dirs_internal() should not be called");
  Exceptions::amanzi_throw(mesg);
}

// Cells connected to a face - this function is implemented in each
// mesh framework. The results are cached in the base class
void
MeshLogical::face_get_cells_internal (const Entity_ID faceid,
				      const Parallel_type ptype,
				      Entity_ID_List *cellids) const {
  Errors::Message mesg("DEVELOPER ERROR: face_get_cells_internal() should not be called");
  Exceptions::amanzi_throw(mesg);
}

// edges of a face - this function is implemented in each mesh
// framework. The results are cached in the base class
void
MeshLogical::face_get_edges_and_dirs_internal (const Entity_ID faceid,
					       Entity_ID_List *edgeids,
					       std::vector<int> *edge_dirs,
					       const bool ordered) const {
  Errors::Message mesg("No edges in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}
  

// edges of a cell - this function is implemented in each mesh
// framework. The results are cached in the base class. 
void
MeshLogical::cell_get_edges_internal (const Entity_ID cellid,
				      Entity_ID_List *edgeids) const {
  Errors::Message mesg("No edges in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}
  

// edges and directions of a 2D cell - this function is implemented
// in each mesh framework. The results are cached in the base class.
void
MeshLogical::cell_2D_get_edges_and_dirs_internal (const Entity_ID cellid,
						  Entity_ID_List *edgeids,
						  std::vector<int> *edge_dirs) const {
  Errors::Message mesg("No edges in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}
  

int
MeshLogical::build_columns() const {
  Errors::Message mesg("No columns are buildable in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}


// Cache connectivity info.
void
MeshLogical::cache_cell2face_info() const {
  Errors::Message mesg("DEVELOPER ERROR: cache should be created in finalize()");
  Exceptions::amanzi_throw(mesg);
}

void
MeshLogical::cache_face2cell_info() const {
  Errors::Message mesg("DEVELOPER ERROR: cache should be created in finalize()");
  Exceptions::amanzi_throw(mesg);
}
  
  
int
MeshLogical::compute_cell_geometric_quantities() const {
  Errors::Message mesg("DEVELOPER ERROR: cache should be created in finalize()");
  Exceptions::amanzi_throw(mesg);
}

int
MeshLogical::compute_face_geometric_quantities() const {
  Errors::Message mesg("DEVELOPER ERROR: cache should be created in finalize()");
  Exceptions::amanzi_throw(mesg);
}
  
} // close namespace AmanziMesh
} // close namespace Amanzi

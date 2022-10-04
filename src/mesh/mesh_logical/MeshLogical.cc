/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

/*!
 Logical mesh that can be modified and constructed on the fly.

 Logical mesh is a topologically defined mesh with no real coordinate
 geometry.  By definition it is perfectly parallel with no ghost entities,
 as it is intended to be used along with a normal mesh as a subgrid model.
 As it is not a geomtric mesh, it cannot work with all (many) spatial
 discretizations -- currently only Finite Volume.

 In particular:
  1. nodes do not exist
*/

#include <set>
#include "RegionEnumerated.hh"
#include "MeshLogical.hh"

namespace Amanzi {
namespace AmanziMesh {

//
// Topological constructor of a MeshLogical splits topology from geometry.
//
MeshLogical::MeshLogical(const Comm_ptr_type& comm,
                         const std::vector<Entity_ID_List>& face_cell_ids,
                         const std::vector<std::vector<int> >& face_cell_dirs,
                         const Teuchos::RCP<const Teuchos::ParameterList>& plist)
    : Mesh(comm, Teuchos::null, plist, true, false)
{
  logical_ = true;
  if (face_cell_ids.size() != face_cell_dirs.size()) {
    Errors::Message mesg("MeshLogical topology created with bad data");
    Exceptions::amanzi_throw(mesg);
  }

  set_space_dimension(3);
  set_manifold_dimension(1);

  // face-cell connectivity topology
  face_cell_ids_ = face_cell_ids;

  // Count number of cells referenced, and check that the number of cells
  // referenced is equal to the largest id referenced (+1 for 0)
  int c_max = -1;
  std::set<int> cells;
  for (auto f : face_cell_ids_) {
    for (auto c : f) {
      cells.insert(c);
      c_max = std::max(c, c_max);
    }
  }
  int num_cells = cells.size();
  if (num_cells != c_max + 1) {
    Errors::Message mesg("MeshLogical expects cells to be 0-indexed and contiguous.");
    Exceptions::amanzi_throw(mesg);
  }    

  // geometric info that must get set later
  cell_volumes_.resize(num_cells, -1.);
  face_areas_.resize(face_cell_ids_.size(), -1.);
  
  face_normals_.resize(face_cell_ids_.size());
  face_cell_ptype_.resize(face_cell_ids_.size());

  // populate cell, extra face info
  cell_face_ids_.resize(num_cells);
  cell_face_dirs_.resize(num_cells);
  cell_face_bisectors_.resize(num_cells);

  int f_id=0;
  for (auto & f : face_cell_ids_) {
    cell_face_ids_[f[0]].push_back(f_id);
    cell_face_dirs_[f[0]].push_back(face_cell_dirs[f_id][0]);
    cell_face_bisectors_[f[0]].emplace_back(AmanziGeometry::Point());

    face_cell_ptype_[f_id].push_back(Parallel_type::OWNED);
    face_normals_[f_id].resize(2);

    // 1s complement the cell to make it consistent with the dir.  While this
    // replicates information, it gives a lookup by both face AND cell of the
    // face dir wrt to the cell with no extra storage.  That said, it means
    // that indices MUST be signed ints, which limits the size of the mesh
    // that can fit 32 bit ints...
    if (face_cell_dirs[f_id][0] < 0) f[0] = ~f[0];

    if (f.size() > 1 && f[1] >= 0) {
      cell_face_ids_[f[1]].push_back(f_id);
      cell_face_dirs_[f[1]].push_back(face_cell_dirs[f_id][1]);
      cell_face_bisectors_[f[1]].emplace_back(AmanziGeometry::Point());
      
      face_cell_ptype_[f_id].push_back(Parallel_type::OWNED);
      if (face_cell_dirs[f_id][1] < 0) f[1] = ~f[1];

    } else {
      face_cell_ptype_[f_id].push_back(Parallel_type::PTYPE_UNKNOWN);
    }      
    f_id++;
  }

  // toggle flags
  cell_geometry_precomputed_ = false;
  face_geometry_precomputed_ = false;
  cell2face_info_cached_ = true;
  faces_requested_ = true;
  face2cell_info_cached_ = true;

  // build epetra maps
  init_maps();
}
  

//
// MeshLogical constructor includes geometry.
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
MeshLogical::MeshLogical(const Comm_ptr_type& comm,
                         const std::vector<Entity_ID_List>& face_cell_ids,
                         const std::vector<std::vector<int> >& face_cell_dirs,
                         const std::vector<double>& cell_volumes,
                         const std::vector<double>& face_areas,
                         const std::vector<std::vector<AmanziGeometry::Point> >& face_cell_bisectors,
                         const std::vector<AmanziGeometry::Point>* cell_centroids,
                         const Teuchos::RCP<const Teuchos::ParameterList>& plist)
  : MeshLogical(comm, face_cell_ids, face_cell_dirs, plist)
{
  // set geometry
  set_logical_geometry(&cell_volumes, &face_areas, &face_cell_bisectors, cell_centroids);
}


void MeshLogical::get_logical_geometry(std::vector<double>* const cell_volumes,
        std::vector<double>* const face_areas,
        std::vector<std::vector<AmanziGeometry::Point> >* const face_cell_bisectors,
        std::vector<AmanziGeometry::Point>* const cell_centroids) const
{
  if (cell_volumes) *cell_volumes = cell_volumes_;
  if (face_areas) *face_areas = face_areas_;
  if (cell_centroids) *cell_centroids = cell_centroids_;

  // not sure what would be useful for user here... wait til we have a user and fix it!
  if (face_cell_bisectors) {
    face_cell_bisectors->resize(num_entities(FACE, Parallel_type::ALL));
  }
}


void MeshLogical::set_logical_geometry(std::vector<double> const* const cell_volumes,
        std::vector<double> const* const face_areas,
        std::vector<std::vector<AmanziGeometry::Point> > const* const face_cell_bisectors,
        std::vector<AmanziGeometry::Point> const* const cell_centroids)
{
  auto n_cells = num_entities(CELL, Parallel_type::OWNED);
  auto n_faces = num_entities(FACE, Parallel_type::OWNED);
  
  if (cell_volumes && n_cells != cell_volumes->size()) {
    Errors::Message mesg("MeshLogical::set_logical_geometry() called with bad data");
    Exceptions::amanzi_throw(mesg);
  }
  if (face_areas && n_faces != face_areas->size()) {
    Errors::Message mesg("MeshLogical::set_logical_geometry() called with bad data");
    Exceptions::amanzi_throw(mesg);
  }
  if (face_cell_bisectors && n_faces != face_cell_bisectors->size()) {
    Errors::Message mesg("MeshLogical::set_logical_geometry() called with bad data");
    Exceptions::amanzi_throw(mesg);
  }
  if (cell_centroids && n_cells != cell_centroids->size()) {
    Errors::Message mesg("MeshLogical::set_logical_geometry() called with bad data");
    Exceptions::amanzi_throw(mesg);
  }


  if (cell_volumes) cell_volumes_ = *cell_volumes;
  if (face_areas) face_areas_ = *face_areas;

  if (face_cell_bisectors) {
    cell_face_bisectors_.clear();
    cell_face_bisectors_.resize(n_cells);

    for (Entity_ID f=0; f!=n_faces; ++f) {
      Entity_ID_List f_cells;
      face_get_cells(f, Parallel_type::ALL, &f_cells);
      AMANZI_ASSERT((*face_cell_bisectors)[f].size() == f_cells.size());
      
      for (int c_index=0; c_index!=f_cells.size(); ++c_index) {
        Entity_ID c = f_cells[c_index];
        cell_face_bisectors_[c].push_back((*face_cell_bisectors)[f][c_index]);

        face_normals_[f][c_index] = face_areas_[f] / AmanziGeometry::norm((*face_cell_bisectors)[f][c_index])
                                    * (*face_cell_bisectors)[f][c_index];
      }
    }
  }

  // do this after bisectors
  if (cell_centroids) {
    cell_centroids_ = *cell_centroids;

    // face centroids constructed from normals, lengths, cell
    // centroids using assumption of perpendicular bisector
    face_centroids_.resize(face_cell_ids_.size());
    for (int f=0; f!=face_cell_ids_.size(); ++f) {
      AmanziMesh::Entity_ID_List f_cells;
      face_get_cells(f, AmanziMesh::Parallel_type::ALL, &f_cells);

      auto c0 = f_cells[0];
      int f_index = 0;
      for (auto ff : cell_face_ids_[c0]) {
        if (ff == f) break;
        f_index++;
      }
      AMANZI_ASSERT(f_index < cell_face_ids_[c0].size());
      auto face_centroid_left = cell_centroids_[c0] + cell_face_bisectors_[c0][f_index];

      if (f_cells.size() == 2) {
        auto c1 = f_cells[1];
        int f_index2 = 0;
        for (auto ff : cell_face_ids_[c1]) {
          if (ff == f) break;
          f_index2++;
        }
        AMANZI_ASSERT(f_index2 < cell_face_ids_[c1].size());
        auto face_centroid_right = cell_centroids_[c1] + cell_face_bisectors_[c1][f_index2];
        face_centroids_[f] = (face_centroid_right + face_centroid_left) / 2.;
      } else {
        face_centroids_[f] = face_centroid_left;
      }
    }
  }
  
  cell_geometry_precomputed_ = true;
  face_geometry_precomputed_ = true;
}




// build maps
void
MeshLogical::init_maps() {
  // cell map
  maps_[CELL] = Teuchos::rcp(new Epetra_Map(-1, cell_face_ids_.size(), 0, *comm_));

  // face map
  Teuchos::RCP<Epetra_Map> face_map =
      Teuchos::rcp(new Epetra_Map(-1, face_cell_ids_.size(), 0, *comm_));
  maps_[FACE] = face_map;

  // exterior face map
  std::vector<int> extface_ids;
  int f_id = 0;
  for (std::vector<Entity_ID_List>::iterator f=face_cell_ids_.begin();
       f!=face_cell_ids_.end(); ++f) {
    if (f->size() == 1) {
      extface_ids.push_back(face_map->GID(f_id));
    }
    f_id++;
  }
  maps_[BOUNDARY_FACE] = Teuchos::rcp(new Epetra_Map(-1, extface_ids.size(),
          &extface_ids[0], 0, *comm_));

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
  if (cell_face_ids_ != other.cell_face_ids_) return false;
  if (face_cell_ids_ != other.face_cell_ids_) return false;

  if (cell_volumes_.size() != other.cell_volumes_.size()) return false;
  for (size_t i=0; i!=cell_volumes_.size(); ++i) {
    if (std::abs(cell_volumes_[i] - other.cell_volumes_[i]) > _eps) return false;
  }

  if (face_normals_.size() != other.face_normals_.size()) return false;
  for (size_t i=0; i!=face_normals_.size(); ++i) {
    if (AmanziGeometry::norm(face_normals_[i][0] - other.face_normals_[i][0]) > _eps) return false;
  }

  for (size_t i=0; i!=face_normals_.size(); ++i) {
    if (AmanziGeometry::norm(face_normals_[i][1] - other.face_normals_[i][1]) > _eps) return false;
  }

  if (cell_centroids_.size() != other.cell_centroids_.size()) return false;
  for (size_t i=0; i!=cell_centroids_.size(); ++i) {
    if (AmanziGeometry::norm(cell_centroids_[i] - other.cell_centroids_[i]) > _eps) return false;
  }

  if (face_centroids_.size() != other.face_centroids_.size()) return false;
  for (size_t i=0; i!=face_centroids_.size(); ++i) {
    if (AmanziGeometry::norm(face_centroids_[i] - other.face_centroids_[i]) > _eps) return false;
  }

  if (cell_face_bisectors_.size() != other.cell_face_bisectors_.size()) return false;
  for (size_t i=0; i!=cell_face_bisectors_.size(); ++i) {
    if (cell_face_bisectors_[i].size() != other.cell_face_bisectors_[i].size()) return false;
    for (size_t j=0; j!=cell_face_bisectors_[i].size(); ++j)
      if (AmanziGeometry::norm(cell_face_bisectors_[i][j] - other.cell_face_bisectors_[i][j]) > _eps) return false;
  }

  return true;
}


// Get parallel type of entity - OWNED, GHOST, ALL (See MeshDefs.hh)
Parallel_type
MeshLogical::entity_get_ptype(const Entity_kind kind,
                              const Entity_ID entid) const {
  return Parallel_type::OWNED;
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
// particular category (OWNED, Parallel_type::GHOST, Parallel_type::ALL)
unsigned int
MeshLogical::num_entities(const Entity_kind kind,
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
MeshLogical::cell_get_nodes(const Entity_ID cellid,
                            Entity_ID_List *nodeids) const {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}


// Get the bisectors, i.e. vectors from cell centroid to face centroids.
void
MeshLogical::cell_get_faces_and_bisectors(
    const Entity_ID cellid,
    Entity_ID_List *faceids,
    std::vector<AmanziGeometry::Point> *bisectors,
    const bool ordered) const {
  if (faceids) *faceids = cell_face_ids_[cellid];
  if (bisectors) *bisectors = cell_face_bisectors_[cellid];
}


// Get nodes of face
// On a distributed mesh, all nodes (OWNED or GHOST) of the face
// are returned
// In 3D, the nodes of the face are returned in ccw order consistent
// with the face normal
// In 2D, nfnodes is 2
void
MeshLogical::face_get_nodes(const Entity_ID faceid,
                            Entity_ID_List *nodeids) const {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}


// Get nodes of edge
void
MeshLogical::edge_get_nodes(const Entity_ID edgeid,
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
MeshLogical::node_get_cells(const Entity_ID nodeid,
                            const Parallel_type ptype,
                            Entity_ID_List *cellids) const {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}


// Faces of type 'ptype' connected to a node - The order of faces is
// not guarnateed to be the same for corresponding nodes on
// different processors
void
MeshLogical::node_get_faces(const Entity_ID nodeid,
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
// guaranteed in general except when ptype = ALL, in which case
// the cellids will correcpond to cells across the respective
// faces given by cell_get_faces
void
MeshLogical::cell_get_face_adj_cells(const Entity_ID cellid,
                                     const Parallel_type ptype,
                                     Entity_ID_List *fadj_cellids) const {
  fadj_cellids->clear();
  Entity_ID_List faces;
  cell_get_faces_and_bisectors(cellid, &faces, NULL, false);
  fadj_cellids->reserve(faces.size());
  for (auto f : faces) {
    Entity_ID_List cells;
    face_get_cells(f, ptype, &cells);
    if (cells[0] == cellid) {
      if (cells.size() == 2) {
        fadj_cellids->emplace_back(cells[1]);
      }
    } else if (cells.size() == 2 && cells[1] == cellid) {
      fadj_cellids->emplace_back(cells[0]);
    } else {
      Errors::Message mesg("Topological issue in MeshLogical.");
      Exceptions::amanzi_throw(mesg);
    }
  }
}


//
// Mesh entity geometry
//--------------

// Node coordinates - 3 in 3D and 2 in 2D
void
MeshLogical::node_get_coordinates(const Entity_ID nodeid,
                                  AmanziGeometry::Point *ncoord) const {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}


//
// Mesh modification
//-------------------

// Set coordinates of node
void
MeshLogical::node_set_coordinates(const Entity_ID nodeid,
                                  const AmanziGeometry::Point ncoord) {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}


void
MeshLogical::node_set_coordinates(const Entity_ID nodeid,
                                  const double *ncoord) {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}


// Deformation not supported.
int
MeshLogical::deform(const std::vector<double>& target_cell_volumes_in,
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
MeshLogical::cell_map(bool include_ghost) const {
  return *maps_.at(CELL);
}

const Epetra_Map&
MeshLogical::face_map(bool include_ghost) const {
  return *maps_.at(FACE);
}

const Epetra_Map&
MeshLogical::node_map(bool include_ghost) const {
  Errors::Message mesg("No nodes in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
  return *maps_.at(NODE);
}


const Epetra_Map&
MeshLogical::exterior_face_map(bool include_ghost) const {
  return *maps_.at(BOUNDARY_FACE);
}


const Epetra_Map&
MeshLogical::exterior_node_map(bool include_ghost) const {
  Errors::Message mesg("No exterior node map in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
  throw(mesg);
}


// Epetra importer that will allow apps to import values from a
// Epetra vector defined on all owned faces into an Epetra vector
// defined only on exterior faces
const Epetra_Import&
MeshLogical::exterior_face_importer(void) const {
  return *exterior_face_importer_;
}


//
// Mesh Sets for ICs, BCs, Material Properties and whatever else
//--------------------------------------------------------------

// Get list of entities of type 'category' in set
void
MeshLogical::get_set_entities(const std::string& setname,
                              const Entity_kind kind,
                              const Parallel_type ptype,
                              Entity_ID_List *entids) const {

  Teuchos::RCP<const AmanziGeometry::Region> rgn = geometric_model_->FindRegion(setname);

  if (rgn->get_type() == AmanziGeometry::RegionType::ALL) {
    int nent = num_entities(kind, ptype);
    entids->resize(num_entities(kind, ptype));
    for (int i=0; i!=nent; ++i) {
      (*entids)[i] = i;
    }
    return;

  } else if (rgn->get_type() == AmanziGeometry::RegionType::ENUMERATED) {
    Teuchos::RCP<const AmanziGeometry::RegionEnumerated> esrgn =
        Teuchos::rcp_static_cast<const AmanziGeometry::RegionEnumerated>(rgn);

    if ((esrgn->entity_str() == "CELL" && kind == CELL) ||
        (esrgn->entity_str() == "FACE" && kind == FACE)) {
      *entids = esrgn->entities();
    }
  } else {
    Errors::Message mesg("MeshLogical currently only supports regions of type \"region: all\" or \"region: enumerated\"");
    Exceptions::amanzi_throw(mesg);
  }
  return;
}


void
MeshLogical::get_set_entities_and_vofs(const std::string& setname,
                                       const Entity_kind kind,
                                       const Parallel_type ptype,
                                       Entity_ID_List *entids,
                                       std::vector<double> *vofs) const {
  get_set_entities(setname, kind, ptype, entids);
  return;
}


// Miscellaneous functions
void
MeshLogical::write_to_exodus_file(const std::string filename) const {
  // don't know how to do this! FIXME! --etc
}


// Geometry
int
MeshLogical::compute_cell_geometry_(const Entity_ID cellid,
        double *volume,
        AmanziGeometry::Point *centroid) const {
  // this is a placeholder, these cannot be recomputed
  if (volume) *volume = cell_volumes_[cellid];

  if (centroid) {
    if (cell_centroids_.size() > 0) {
      *centroid = cell_centroids_[cellid];
    } else {
      *centroid = AmanziGeometry::Point();
    }
  }
  return 1;
}


int
MeshLogical::compute_face_geometry_(
    const Entity_ID faceid,
    double *area,
    AmanziGeometry::Point *centroid,
    std::vector<AmanziGeometry::Point> *normals) const {
  // this is a placeholder, these cannot be recomputed
  if (area) *area = face_areas_[faceid];
  if (centroid) *centroid = AmanziGeometry::Point();
  if (normals) *normals = face_normals_[faceid];
  return 1;
}


// get faces of a cell and directions in which it is used - this function
// is implemented in each mesh framework. The results are cached in
// the base class
void
MeshLogical::cell_get_faces_and_dirs_internal_(
    const Entity_ID cellid,
    Entity_ID_List *faceids,
    std::vector<int> *face_dirs,
    const bool ordered) const {
  Errors::Message mesg("DEVELOPER ERROR: cell_get_faces_and_dirs_internal_() should not be called");
  Exceptions::amanzi_throw(mesg);
}


// Cells connected to a face - this function is implemented in each
// mesh framework. The results are cached in the base class
void
MeshLogical::face_get_cells_internal_(const Entity_ID faceid,
                                      const Parallel_type ptype,
                                      Entity_ID_List *cellids) const {
  Errors::Message mesg("DEVELOPER ERROR: face_get_cells_internal_() should not be called");
  Exceptions::amanzi_throw(mesg);
}


// edges of a face - this function is implemented in each mesh
// framework. The results are cached in the base class
void
MeshLogical::face_get_edges_and_dirs_internal_(
    const Entity_ID faceid,
    Entity_ID_List *edgeids,
    std::vector<int> *edge_dirs,
    const bool ordered) const {
  Errors::Message mesg("No edges in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}


// edges of a cell - this function is implemented in each mesh
// framework. The results are cached in the base class.
void
MeshLogical::cell_get_edges_internal_(const Entity_ID cellid,
                                      Entity_ID_List *edgeids) const {
  Errors::Message mesg("No edges in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}


int
MeshLogical::build_columns_() const {
  Errors::Message mesg("No columns are buildable in MeshLogical.");
  Exceptions::amanzi_throw(mesg);
  return -1;
}


//
// Note this works on Mesh, but is less useful for a general mesh
// --------------------------------------------------------------------------------
bool viewMeshLogical(const Mesh& m, std::ostream& os) {
  if (m.get_comm()->NumProc() != 1) {
    return true;
  }

  os << "cell_centroids, volumes =" << std::endl;
  for (int c=0; c!=m.num_entities(CELL, Parallel_type::OWNED); ++c) {
    os << m.cell_centroid(c) << " " << m.cell_volume(c) << std::endl;
  }
  os << "face_connections, areas =" << std::endl;
  for (int f=0; f!=m.num_entities(FACE, Parallel_type::OWNED); ++f) {
    AmanziMesh::Entity_ID_List fcells;
    m.face_get_cells(f, Parallel_type::ALL, &fcells);
    for (auto c : fcells) os << c << " ";
    os << m.face_area(f) << std::endl;
  }

  os << "cell_sets =" << std::endl;
  for (auto& r : *m.geometric_model()) {
    if (r->get_type() == AmanziGeometry::RegionType::ENUMERATED) {
      AmanziMesh::Entity_ID_List set;
      m.get_set_entities(r->get_name(), AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &set);
      os << r->get_name() << " ";
      for (auto e : set) os << e << " ";
      os << std::endl;
    }
  }

  return false;
  
}
  


}  // namespace AmanziMesh
}  // namespace Amanzi

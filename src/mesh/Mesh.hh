#ifndef _AMANZI_MESH_H_
#define _AMANZI_MESH_H_

#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>

#include <memory>

#include "MeshDefs.hh"
#include "Cell_topology.hh"
#include "Point.hh"
#include "GeometricModel.hh"
#include "Region.hh"


#include <map>

namespace Amanzi
{

namespace AmanziMesh
{

class Mesh
{

 private:

  unsigned int celldim, spacedim;
  mutable bool geometry_precomputed;
  mutable std::vector<double> cell_volumes, face_areas;
  mutable std::vector<AmanziGeometry::Point> cell_centroids,
    face_centroids, face_normals;
  AmanziGeometry::GeometricModelPtr geometric_model_;

  const Epetra_MpiComm *comm; // temporary until we get an amanzi communicator


  int compute_geometric_quantities() const;
  int compute_cell_geometry(const Entity_ID cellid, 
                            double *volume, 
                            AmanziGeometry::Point *centroid) const;
  int compute_face_geometry(const Entity_ID faceid, 
                            double *area, 
                            AmanziGeometry::Point *centroid, 
                            AmanziGeometry::Point *normal) const;


 protected:

 public:

  // constructor

  Mesh()
    : spacedim(3), celldim(3), geometry_precomputed(false), comm(NULL),
      geometric_model_(NULL)
  {
  }

  // destructor

  virtual ~Mesh() {}

  inline
  void set_space_dimension(const unsigned int dim) {
    spacedim = dim;
  }

  inline
  unsigned int space_dimension() const
  {
    return spacedim;
  }

  inline
  void set_cell_dimension(const unsigned int dim) {
    celldim = dim;   // 3 is solid mesh, 2 is surface mesh
  }

  inline
  unsigned int cell_dimension() const
  {
    return celldim;
  }

  inline
  void set_geometric_model(const AmanziGeometry::GeometricModelPtr &gm) {
    geometric_model_ = gm;
  }

  inline
  AmanziGeometry::GeometricModelPtr geometric_model() const
  {
    return geometric_model_;
  }



  // Get parallel type of entity

  virtual
  Parallel_type entity_get_ptype(const Entity_kind kind,
                                 const Entity_ID entid) const = 0;




  // Get cell type

  virtual
  Cell_type cell_get_type(const Entity_ID cellid) const = 0;




  //
  // General mesh information
  // -------------------------
  //

  // Number of entities of any kind (cell, face, node) and in a
  // particular category (OWNED, GHOST, USED)

  virtual
  unsigned int num_entities (const Entity_kind kind,
                             const Parallel_type ptype) const = 0;


  // Global ID of any entity

  virtual
  unsigned int GID(const Entity_ID lid, const Entity_kind kind) const = 0;



  //
  // Mesh Entity Adjacencies
  //-------------------------


  // Downward Adjacencies
  //---------------------

  // Get faces of a cell.

  // On a distributed mesh, this will return all the faces of the
  // cell, OWNED or GHOST. The faces will be returned in a standard
  // order according to Exodus II convention.

  virtual
  void cell_get_faces (const Entity_ID cellid,
                       Entity_ID_List *faceids) const = 0;


  // Get directions in which a cell uses face
  // In 3D, direction is 1 if face normal points out of cell
  // and -1 if face normal points into cell
  // In 2D, direction is 1 if face/edge is defined in the same
  // direction as the cell polygon, and -1 otherwise

  virtual
  void cell_get_face_dirs (const Entity_ID cellid,
                           std::vector<int> *face_dirs) const = 0;


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

  virtual
  void cell_get_faces_and_dirs (const Entity_ID cellid,
                                Entity_ID_List *faceids,
                                std::vector<int> *face_dirs,
				const bool ordered=false) const = 0;



  // Get nodes of cell
  // On a distributed mesh, all nodes (OWNED or GHOST) of the cell
  // are returned
  // Nodes are returned in a standard order (Exodus II convention)
  // STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES in 3D
  // For a general polyhedron this will return the nodes in
  // arbitrary order
  // In 2D, the nodes of the polygon will be returned in ccw order
  // consistent with the face normal

  virtual
  void cell_get_nodes (const Entity_ID cellid,
                       Entity_ID_List *nodeids) const = 0;


  // Get nodes of face
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal
  // In 2D, nfnodes is 2

  virtual
  void face_get_nodes (const Entity_ID faceid,
                       Entity_ID_List *nodeids) const = 0;



  // Upward adjacencies
  //-------------------

  // Cells of type 'ptype' connected to a node

  virtual
  void node_get_cells (const Entity_ID nodeid,
                       const Parallel_type ptype,
                       Entity_ID_List *cellids) const = 0;

  // Faces of type 'ptype' connected to a node

  virtual
  void node_get_faces (const Entity_ID nodeid,
                       const Parallel_type ptype,
                       Entity_ID_List *faceids) const = 0;

  // Get faces of ptype of a particular cell that are connected to the
  // given node

  virtual
  void node_get_cell_faces (const Entity_ID nodeid,
                            const Entity_ID cellid,
                            const Parallel_type ptype,
                            Entity_ID_List *faceids) const = 0;

  // Cells connected to a face

  virtual
  void face_get_cells (const Entity_ID faceid,
                       const Parallel_type ptype,
                       Entity_ID_List *cellids) const = 0;



  // Same level adjacencies
  //-----------------------

  // Face connected neighboring cells of given cell of a particular ptype
  // (e.g. a hex has 6 face neighbors)

  // The order in which the cellids are returned cannot be
  // guaranteed in general except when ptype = USED, in which case
  // the cellids will correcpond to cells across the respective
  // faces given by cell_get_faces

  virtual
  void cell_get_face_adj_cells(const Entity_ID cellid,
                               const Parallel_type ptype,
                               Entity_ID_List *fadj_cellids) const = 0;

  // Node connected neighboring cells of given cell
  // (a hex in a structured mesh has 26 node connected neighbors)
  // The cells are returned in no particular order

  virtual
  void cell_get_node_adj_cells(const Entity_ID cellid,
                               const Parallel_type ptype,
                               Entity_ID_List *nadj_cellids) const = 0;



  //
  // Mesh Topology for viz
  //----------------------
  //
  // We need a special function because certain types of degenerate
  // hexes will not be recognized as any standard element type (hex,
  // pyramid, prism or tet). The original topology of this element
  // without any collapsed nodes will be returned by this call.


  // Original cell type

  virtual
  Cell_type cell_get_type_4viz(const Entity_ID cellid) const = 0;


  // See cell_get_nodes for details on node ordering

  virtual
  void cell_get_nodes_4viz (const Entity_ID cellid,
                            Entity_ID_List *nodeids) const = 0;



  //
  // Mesh entity geometry
  //--------------
  //

  // Node coordinates - 3 in 3D and 2 in 2D

  virtual
  void node_get_coordinates (const Entity_ID nodeid,
                             AmanziGeometry::Point *ncoord) const = 0;


  // Face coordinates - conventions same as face_to_nodes call
  // Number of nodes is the vector size divided by number of spatial dimensions

  virtual
  void face_get_coordinates (const Entity_ID faceid,
                             std::vector<AmanziGeometry::Point> *fcoords) const = 0;

  // Coordinates of cells in standard order (Exodus II convention)
  // STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES IN 3D
  // For a general polyhedron this will return the node coordinates in
  // arbitrary order
  // Number of nodes is vector size divided by number of spatial dimensions

  virtual
  void cell_get_coordinates (const Entity_ID cellid,
                             std::vector<AmanziGeometry::Point> *ccoords) const = 0;


  // Mesh entity geometry
  //--------------
  //


  // Volume/Area of cell

  double cell_volume (const Entity_ID cellid, const bool recompute=false) const;

  // Area/length of face

  double face_area(const Entity_ID faceid, const bool recompute=false) const;

  // Centroid of cell

  AmanziGeometry::Point cell_centroid (const Entity_ID cellid, const bool recompute=false) const;

  // Centroid of face

  AmanziGeometry::Point face_centroid (const Entity_ID faceid, const bool recompute=false) const;

  // Normal to face
  // The vector is normalized and then weighted by the area of the face

  AmanziGeometry::Point face_normal (const Entity_ID faceid, const bool recompute=false) const;

  // Point in cell

  bool point_in_cell (const AmanziGeometry::Point &p, 
                      const Entity_ID cellid) const;


  //
  // Mesh modification
  //-------------------

  // Set coordinates of node

  virtual
  void node_set_coordinates (const Entity_ID nodeid,
                             const AmanziGeometry::Point ncoord) = 0;


  virtual
  void node_set_coordinates (const Entity_ID nodeid,
                             const double *ncoord) = 0;



  // Deform the mesh by moving given nodes to given coordinates
  // If the flag keep_valid is true, then the nodes are moved
  // only as much as possible without making the mesh invalid
  // The final positions of the nodes is returned in final_positions


  int deform (const Entity_ID_List nodeids,
              const AmanziGeometry::Point_List new_positions,
              const bool keep_valid,
              AmanziGeometry::Point_List *final_positions);

  //
  // Epetra maps
  //------------


  virtual
  const Epetra_Map& cell_epetra_map (const bool include_ghost) const = 0;

  virtual
  const Epetra_Map& face_epetra_map (const bool include_ghost) const = 0;

  virtual
  const Epetra_Map& node_epetra_map (const bool include_ghost) const = 0;




  //
  // Mesh Sets for ICs, BCs, Material Properties and whatever else
  //--------------------------------------------------------------
  //

  // Number of sets containing entities of type 'kind' in mesh
  // 
  // DEPRECATED due to ambiguity in determining what types of sets
  // some regions are supposed to create (a planar region can 
  // result in sidesets or nodesets

  unsigned int num_sets(const Entity_kind kind) const;


  // Ids of sets containing entities of 'kind'
  // 
  // DEPRECATED due to ambiguity in determining what types of sets
  // some regions are supposed to create (a planar region can 
  // result in sidesets or nodesets

  void get_set_ids (const Entity_kind kind,
                    Set_ID_List *setids) const;




  // Is this is a valid ID of a set containing entities of 'kind'

  bool valid_set_id (const Set_ID setid,
                     const Entity_kind kind) const;

  // Is this is a valid ID of a set containing entities of 'kind'

  bool valid_set_name (const std::string setname,
                       const Entity_kind kind) const;


  // Get set ID from set name - returns 0 if no match is found
  
  unsigned int set_id_from_name(const std::string setname) const;


  // Get set name from set ID - returns 0 if no match is found
  
  std::string set_name_from_id(const int setid) const;


  // Get number of entities of type 'category' in set

  virtual
  unsigned int get_set_size (const Set_ID setid,
                             const Entity_kind kind,
                             const Parallel_type ptype) const = 0;

  virtual
  unsigned int get_set_size (const Set_Name setname,
                             const Entity_kind kind,
                             const Parallel_type ptype) const = 0;

  virtual
  unsigned int get_set_size (const char *setname,
                             const Entity_kind kind,
                             const Parallel_type ptype) const = 0;


  // Get list of entities of type 'category' in set

  virtual
  void get_set_entities (const Set_ID setid,
                         const Entity_kind kind,
                         const Parallel_type ptype,
                         Entity_ID_List *entids) const = 0;

  virtual
  void get_set_entities (const Set_Name setname,
                         const Entity_kind kind,
                         const Parallel_type ptype,
                         Entity_ID_List *entids) const = 0;

  virtual
  void get_set_entities (const char *setname,
                         const Entity_kind kind,
                         const Parallel_type ptype,
                         Entity_ID_List *entids) const = 0;



  // communicator access
  // temporary until we set up an amanzi communicator

  // Changing this to return Epetra_MpiComm * because the 
  // stk::ParalllelMachine cannot be initialized with anything
  // other than an MpiComm type - so we can't do serial builds anyway 

  inline
  const Epetra_MpiComm* get_comm() const {
    return comm;
  }

  inline
  void set_comm(const Epetra_MpiComm *incomm) {
    comm = incomm;
  }


  // Temporary routines for backward compatibility

  void cell_to_faces (unsigned int cell,
                      std::vector<unsigned int>::iterator begin,
                      std::vector<unsigned int>::iterator end) const;

  void cell_to_faces (unsigned int cell,
                      unsigned int* begin, unsigned int *end) const;


  void cell_to_face_dirs (unsigned int cell,
                          std::vector<int>::iterator begin,
                          std::vector<int>::iterator end) const;
  void cell_to_face_dirs (unsigned int cell,
                          int * begin, int * end) const;



  void cell_to_nodes (unsigned int cell,
                      std::vector<unsigned int>::iterator begin,
                      std::vector<unsigned int>::iterator end) const;
  void cell_to_nodes (unsigned int cell,
                      unsigned int * begin, unsigned int * end) const;




  void face_to_nodes (unsigned int face,
                      std::vector<unsigned int>::iterator begin,
                      std::vector<unsigned int>::iterator end) const;
  void face_to_nodes (unsigned int face,
                      unsigned int * begin, unsigned int * end) const;



  void node_to_coordinates (unsigned int node,
                            std::vector<double>::iterator begin,
                            std::vector<double>::iterator end) const;
  void node_to_coordinates (unsigned int node,
                            double * begin,
                            double * end) const;

  void face_to_coordinates (unsigned int face,
                            std::vector<double>::iterator begin,
                            std::vector<double>::iterator end) const;
  void face_to_coordinates (unsigned int face,
                            double * begin,
                            double * end) const;

  void cell_to_coordinates (unsigned int cell,
                            std::vector<double>::iterator begin,
                            std::vector<double>::iterator end) const;
  void cell_to_coordinates (unsigned int cell,
                            double * begin,
                            double * end) const;


  const Epetra_Map& cell_map (bool include_ghost) const
  {
    return cell_epetra_map (include_ghost);
  };

  const Epetra_Map& face_map (bool include_ghost) const
  {
    return face_epetra_map (include_ghost);
  };

  const Epetra_Map& node_map (bool include_ghost) const
  {
    return node_epetra_map (include_ghost);
  };


  unsigned int count_entities (Entity_kind kind,
                               Parallel_type ptype) const;


  void get_set (unsigned int set_id, Entity_kind kind,
                Parallel_type ptype,
                std::vector<unsigned int>::iterator begin,
                std::vector<unsigned int>::iterator end) const;
  void get_set (unsigned int set_id, Entity_kind kind,
                Parallel_type ptype,
                unsigned int * begin,
                unsigned int * end) const;

  // Id numbers
  // DEPRECATED - DO NOT USE
  void get_set_ids (Entity_kind kind,
                    std::vector<unsigned int>::iterator begin,
                    std::vector<unsigned int>::iterator end) const;
  void get_set_ids (Entity_kind kind,
                    unsigned int * begin,
                    unsigned int * end) const;

}; // End class Mesh


} // end namespace AmanziMesh

} // end namespace Amanzi




#endif /* _AMANZI_MESH_H_ */

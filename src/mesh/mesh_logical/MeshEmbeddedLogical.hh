/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_EMBEDDED_LOGICAL_MESH_H_
#define AMANZI_EMBEDDED_LOGICAL_MESH_H_

#include <Epetra_Map.h>
#include <AmanziComm.hh>
#include <Epetra_SerialComm.h>

#include <memory>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Mesh.hh"
#include "Region.hh"

#include "VerboseObject.hh"
#include "dbc.hh"
#include "errors.hh"

namespace Amanzi {

namespace AmanziMesh {

// Logical mesh that can be modified and constructed on the fly.
//
// Logical mesh is a topologically defined mesh with no real
// coordinate geometry.  By definition it is perfectly parallel with
// no ghost entities, as it is intended to be used along with a normal
// mesh as a subgrid model.  As it is not a geomtric mesh, it cannot
// work with all (many) discretizations -- currently only Finite
// Volume.
//
// In particular:
//  1. nodes do not exist

class MeshEmbeddedLogical : public Mesh {
 public:
  //
  // MeshEmbeddedLogical Constructor
  //
  // Combines two meshes, each representing their own domain, into a
  // single mesh via set of specified connections.  The result is a
  // logical mesh in the sense that it provides a limited interface.
  //
  //  - face_cell_list          : length nfaces array of length 2 arrays
  //                              defining the topology
  //  - face_cell_lengths       : length of the cell-to-face connection
  //  - face_area_normals       : length nfaces array of normals of the
  //                              face, points from cell 1 to 2 in
  //                              face_cell_list topology, magnitude
  //                              is area
  MeshEmbeddedLogical(const Comm_ptr_type& comm,
                      Teuchos::RCP<Mesh> bg_mesh,
                      Teuchos::RCP<Mesh> log_mesh,
                      const std::vector<std::vector<int> >& face_cell_list,
                      const std::vector<std::vector<double> >& face_cell_lengths,
                      const std::vector<AmanziGeometry::Point>& face_area_normals,
                      const Teuchos::RCP<const Teuchos::ParameterList>& plist=Teuchos::null);

  // Get parallel type of entity - OWNED, GHOST, ALL (See MeshDefs.hh)
  virtual
  Parallel_type entity_get_ptype(const Entity_kind kind,
                                 const Entity_ID entid) const override;

  // Parent entity in the source mesh if mesh was derived from another mesh
  virtual
  Entity_ID entity_get_parent(const Entity_kind kind, const Entity_ID entid) const override;

  // Get cell type - UNKNOWN, TRI, QUAD, POLYGON, TET, PRISM, PYRAMID, HEX, POLYHED
  // See MeshDefs.hh
  virtual
  Cell_type cell_get_type(const Entity_ID cellid) const override;

  //
  // General mesh information
  // -------------------------
  //
  // Number of entities of any kind (cell, face, node) and in a
  // particular category (OWNED, GHOST, ALL)
  virtual
  unsigned int num_entities(const Entity_kind kind,
                            const Parallel_type ptype) const override;


  // Global ID of any entity
  virtual
  Entity_ID GID(const Entity_ID lid, const Entity_kind kind) const override;



  //
  // Mesh Entity Adjacencies
  //-------------------------


  // Downward Adjacencies
  //---------------------
  // Get nodes of a cell
  virtual
  void cell_get_nodes(const Entity_ID cellid,
                      Entity_ID_List *nodeids) const override;

  // Get the bisectors, i.e. vectors from cell centroid to face centroids.
  virtual
  void cell_get_faces_and_bisectors(const Entity_ID cellid,
          Entity_ID_List *faceids,
          std::vector<AmanziGeometry::Point> *bisectors,
          const bool ordered=false) const override;

  // Get nodes of face
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal
  // In 2D, nfnodes is 2
  virtual
  void face_get_nodes(const Entity_ID faceid,
                      Entity_ID_List *nodeids) const override;


  // Get nodes of edge
  virtual
  void edge_get_nodes(const Entity_ID edgeid,
                      Entity_ID *nodeid0, Entity_ID *nodeid1) const override;

  // Upward adjacencies
  //-------------------

  // Cells of type 'ptype' connected to a node - The order of cells is
  // not guaranteed to be the same for corresponding nodes on
  // different processors
  virtual
  void node_get_cells(const Entity_ID nodeid,
                      const Parallel_type ptype,
                      Entity_ID_List *cellids) const override;


  // Faces of type 'ptype' connected to a node - The order of faces is
  // not guarnateed to be the same for corresponding nodes on
  // different processors
  virtual
  void node_get_faces(const Entity_ID nodeid,
                      const Parallel_type ptype,
                      Entity_ID_List *faceids) const override;

  // Cells of type 'ptype' connected to an edge
  virtual
  void edge_get_cells(const Entity_ID edgeid,
                      const Parallel_type ptype,
                      Entity_ID_List *cellids) const override {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }


  // Same level adjacencies
  //-----------------------

  // Face connected neighboring cells of given cell of a particular ptype
  // (e.g. a hex has 6 face neighbors)

  // The order in which the cellids are returned cannot be
  // guaranteed in general except when ptype = ALL, in which case
  // the cellids will correcpond to cells across the respective
  // faces given by cell_get_faces
  virtual
  void cell_get_face_adj_cells(const Entity_ID cellid,
          const Parallel_type ptype,
          Entity_ID_List *fadj_cellids) const override;


  //
  // Mesh entity geometry
  //--------------
  //

  // Node coordinates - 3 in 3D and 2 in 2D
  virtual
  void node_get_coordinates(const Entity_ID nodeid,
                            AmanziGeometry::Point *ncoord) const override;


  // Face coordinates - conventions same as face_to_nodes call
  // Number of nodes is the vector size divided by number of spatial dimensions
  virtual
  void face_get_coordinates(const Entity_ID faceid,
                            std::vector<AmanziGeometry::Point> *fcoords) const override;

  // Coordinates of cells in standard order (Exodus II convention)
  // STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES IN 3D
  // For a general polyhedron this will return the node coordinates in
  // arbitrary order
  // Number of nodes is vector size divided by number of spatial dimensions
  virtual
  void cell_get_coordinates(const Entity_ID cellid,
                            std::vector<AmanziGeometry::Point> *ccoords) const override;

  //
  // Mesh modification
  //-------------------
  // Set coordinates of node

  virtual
  void node_set_coordinates(const Entity_ID nodeid,
                            const AmanziGeometry::Point ncoord) override;


  virtual
  void node_set_coordinates(const Entity_ID nodeid,
                            const double *ncoord) override;


  // deformation not supported
  virtual
  int deform(const std::vector<double>& target_cell_volumes_in,
             const std::vector<double>& min_cell_volumes_in,
             const Entity_ID_List& fixed_nodes,
             const bool move_vertical) override;

  //
  // Epetra maps
  //------------
  virtual
  const Epetra_Map& cell_map(bool include_ghost) const override;

  virtual
  const Epetra_Map& face_map(bool include_ghost) const override;

  virtual
  const Epetra_Map& node_map(bool include_ghost) const override;

  virtual
  const Epetra_Map& exterior_face_map(bool include_ghost) const override;

  virtual
  const Epetra_Map& exterior_node_map(bool include_ghost) const override;

  // Epetra importer that will allow apps to import values from a
  // Epetra vector defined on all owned faces into an Epetra vector
  // defined only on exterior faces

  virtual
  const Epetra_Import& exterior_face_importer(void) const override;


  //
  // Mesh Sets for ICs, BCs, Material Properties and whatever else
  //--------------------------------------------------------------
  virtual
  bool valid_set_type(const AmanziGeometry::RegionType rtype, const Entity_kind kind) const override {
    return log_mesh_->valid_set_type(rtype, kind) ||
      bg_mesh_->valid_set_type(rtype, kind);
  }

  // Get list of entities of type 'category' in set
  virtual
  void get_set_entities(const std::string& setname,
                        const Entity_kind kind,
                        const Parallel_type ptype,
                        Entity_ID_List *entids) const override;

  virtual
  void get_set_entities_and_vofs(const std::string& setname,
                                 const Entity_kind kind,
                                 const Parallel_type ptype,
                                 Entity_ID_List *entids,
                                 std::vector<double> *vofs) const override;

  // Miscellaneous functions
  virtual
  void write_to_exodus_file(const std::string filename) const override;

 protected:
  // individual versions, if recompute is used
  virtual
  int compute_cell_geometry_(const Entity_ID cellid,
                             double *volume,
                             AmanziGeometry::Point *centroid) const override;
  virtual
  int compute_face_geometry_(const Entity_ID faceid,
                             double *area,
                             AmanziGeometry::Point *centroid,
                             std::vector<AmanziGeometry::Point> *normals) const override;

  // build maps
  void init_maps();


  // get faces of a cell and directions in which it is used - this function
  // is implemented in each mesh framework. The results are cached in
  // the base class
  virtual
  void cell_get_faces_and_dirs_internal_(const Entity_ID cellid,
          Entity_ID_List *faceids,
          std::vector<int> *face_dirs,
          const bool ordered=false) const override;

  // Cells connected to a face - this function is implemented in each
  // mesh framework. The results are cached in the base class
  virtual
  void face_get_cells_internal_(const Entity_ID faceid,
          const Parallel_type ptype,
          Entity_ID_List *cellids) const override;


  // edges of a face - this function is implemented in each mesh
  // framework. The results are cached in the base class
  virtual
  void face_get_edges_and_dirs_internal_(const Entity_ID faceid,
          Entity_ID_List *edgeids,
          std::vector<int> *edge_dirs,
          const bool ordered=true) const override;

  // edges of a cell - this function is implemented in each mesh
  // framework. The results are cached in the base class.
  virtual
  void cell_get_edges_internal_(const Entity_ID cellid,
          Entity_ID_List *edgeids) const override;

  virtual
  int build_columns_() const;

 protected:
  bool initialized_;

  std::map<Entity_kind,Entity_ID> num_entities_owned_;
  std::map<Entity_kind,Entity_ID> num_entities_used_;
  std::map<Entity_kind,Teuchos::RCP<Epetra_Map> > maps_owned_;
  std::map<Entity_kind,Teuchos::RCP<Epetra_Map> > maps_used_;
  std::map<Set_ID,std::vector<int> > regions_;
  Teuchos::RCP<Epetra_Import> exterior_face_importer_;
  std::vector<std::vector<AmanziGeometry::Point> > cell_face_bisectors_;

  Teuchos::RCP<Mesh> bg_mesh_;  // background mesh, typically a Mesh_MSTK
  Teuchos::RCP<Mesh> log_mesh_;  // embedded mesh, typically a MeshLogical
};

}  // namespace AmanziMesh
}  // namespace Amanzi

#endif

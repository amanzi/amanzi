/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//!

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

#ifndef AMANZI_LOGICAL_MESH_H_
#define AMANZI_LOGICAL_MESH_H_

#include <memory>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "AmanziMap.hh"
#include "AmanziComm.hh"

#include "VerboseObject.hh"
#include "dbc.hh"
#include "errors.hh"

#include "Region.hh"
#include "Mesh.hh"

namespace Amanzi {
namespace AmanziMesh {


class MeshLogical : public Mesh {
 public:
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
  //  - cell_centroids_          : (optional, for plotting) length ncell
  //                              array of centroids
  MeshLogical(
    const Comm_ptr_type& comm, const std::vector<Entity_ID_List>& face_cell_ids,
    const std::vector<AmanziGeometry::Point>& face_normals,
    const Teuchos::RCP<const Teuchos::ParameterList>& plist = Teuchos::null);

  MeshLogical(
    const Comm_ptr_type& comm, const std::vector<double>& cell_volumes,
    const std::vector<Entity_ID_List>& face_cell_ids,
    const std::vector<std::vector<double>>& face_cell_lengths,
    const std::vector<AmanziGeometry::Point>& face_area_normals,
    const std::vector<AmanziGeometry::Point>* cell_centroids = nullptr,
    const Teuchos::RCP<const Teuchos::ParameterList>& plist = Teuchos::null);


  void get_logical_geometry(
    std::vector<double>* const cell_volumes,
    std::vector<std::vector<double>>* const cell_face_lengths,
    std::vector<double>* const face_areas,
    std::vector<AmanziGeometry::Point>* const cell_centroids) const;

  void set_logical_geometry(
    std::vector<double> const* const cell_volumes,
    std::vector<std::vector<double>> const* const cell_face_lengths,
    std::vector<double> const* const face_areas,
    std::vector<AmanziGeometry::Point> const* const cell_centroids = NULL);


  // for testing
  bool operator==(const MeshLogical& other);

  // Get parallel type of entity - OWNED, GHOST, ALL (See MeshDefs.hh)
  virtual Parallel_type entity_get_ptype(const Entity_kind kind,
                                         const Entity_ID entid) const override;

  // Parent entity in the source mesh if mesh was derived from another mesh
  virtual Entity_ID entity_get_parent_type(const Entity_kind kind,
                                      const Entity_ID entid) const override;

  // Get cell type - UNKNOWN, TRI, QUAD, ... See MeshDefs.hh
  virtual Cell_type cell_get_type(const Entity_ID cellid) const override;

  //
  // General mesh information
  // -------------------------
  //
  // Number of entities of any kind (cell, face, node) and in a
  // particular category (OWNED, GHOST, ALL)
  virtual unsigned int num_entities(const Entity_kind kind,
                                    const Parallel_type ptype) const override;


  // Global ID of any entity
  virtual Entity_ID
  GID(const Entity_ID lid, const Entity_kind kind) const override;


  //
  // Mesh Entity Adjacencies
  //-------------------------

  // Downward Adjacencies
  //---------------------
  // Get nodes of a cell
  virtual void cell_get_nodes(const Entity_ID cellid,
                              Entity_ID_List& nodeids) const override;

  // Get the bisectors, i.e. vectors from cell centroid to face centroids.
  KOKKOS_INLINE_FUNCTION virtual void 
  cell_get_faces_and_bisectors(const Entity_ID cellid, 
                      Kokkos::View<Entity_ID*>& faceids,
                      Kokkos::View<AmanziGeometry::Point*>& bisectors) const override;

  // Get nodes of face
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal
  // In 2D, nfnodes is 2
  virtual void face_get_nodes(const Entity_ID faceid,
                              Entity_ID_List& nodeids) const override;


  // Get nodes of edge
  virtual void edge_get_nodes(const Entity_ID edgeid, Entity_ID* nodeid0,
                              Entity_ID* nodeid1) const override;

  // Upward adjacencies
  //-------------------

  // Cells of type 'ptype' connected to a node - The order of cells is
  // not guaranteed to be the same for corresponding nodes on
  // different processors
  virtual void node_get_cells(const Entity_ID nodeid, const Parallel_type ptype,
                              Entity_ID_List& cellids) const override;


  // Faces of type 'ptype' connected to a node - The order of faces is
  // not guarnateed to be the same for corresponding nodes on
  // different processors
  virtual void node_get_faces(const Entity_ID nodeid, const Parallel_type ptype,
                              Entity_ID_List& faceids) const override;

  // Get faces of ptype of a particular cell that are connected to the
  // given node - The order of faces is not guarnateed to be the same
  // for corresponding nodes on different processors
  virtual void
  node_get_cell_faces(const Entity_ID nodeid, const Entity_ID cellid,
                      const Parallel_type ptype,
                      Entity_ID_List& faceids) const override;

  // Cells of type 'ptype' connected to an edge - The order of cells is
  // not guaranteed to be the same for corresponding edges on
  // different processors
  virtual void edge_get_cells(const Entity_ID edgeid, const Parallel_type ptype,
                              Entity_ID_List& cellids) const override
  {
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
  virtual void
  cell_get_face_adj_cells(const Entity_ID cellid, const Parallel_type ptype,
                          Entity_ID_List& fadj_cellids) const override;

  // Node connected neighboring cells of given cell
  // (a hex in a structured mesh has 26 node connected neighbors)
  // The cells are returned in no particular order
  virtual void
  cell_get_node_adj_cells(const Entity_ID cellid, const Parallel_type ptype,
                          Entity_ID_List& nadj_cellids)  const override;

  //
  // Mesh entity geometry
  //--------------
  //

  // Node coordinates - 3 in 3D and 2 in 2D
  virtual void
  node_get_coordinates(const Entity_ID nodeid,
                       AmanziGeometry::Point* ncoord) const override;


  // Face coordinates - conventions same as face_to_nodes call
  // Number of nodes is the vector size divided by number of spatial dimensions
  virtual void
  face_get_coordinates(const Entity_ID faceid,
                       std::vector<AmanziGeometry::Point>& fcoords) const override;

  // Coordinates of cells in standard order (Exodus II convention)
  // STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES IN 3D
  // For a general polyhedron this will return the node coordinates in
  // arbitrary order
  // Number of nodes is vector size divided by number of spatial dimensions
  virtual void
  cell_get_coordinates(const Entity_ID cellid,
                       std::vector<AmanziGeometry::Point>& ccoords) const override;

  //
  // Mesh modification
  //-------------------
  // Set coordinates of node

  virtual void
  node_set_coordinates(const Entity_ID nodeid,
                       const AmanziGeometry::Point ncoord) override;

  virtual void
  node_set_coordinates(const Entity_ID nodeid, const double* ncoord) override;


  // deformation not supported
  virtual int
  deform(const std::vector<double>& target_cell_volumes_in,
         const std::vector<double>& min_cell_volumes_in,
         const Entity_ID_List& fixed_nodes, const bool move_vertical) override;

  //
  // Epetra maps
  //------------

  virtual Map_ptr_type cell_map(bool include_ghost) const override;

  virtual Map_ptr_type face_map(bool include_ghost) const override;

  virtual Map_ptr_type node_map(bool include_ghost) const override;

  virtual Map_ptr_type exterior_face_map(bool include_ghost) const override;

  // Epetra importer that will allow apps to import values from a
  // Epetra vector defined on all owned faces into an Epetra vector
  // defined only on exterior faces

  virtual Import_ptr_type exterior_face_importer(void) const override;


  //
  // Mesh Sets for ICs, BCs, Material Properties and whatever else
  //--------------------------------------------------------------

  // Get list of entities of type 'category' in set
  using Mesh::get_set_entities;

  virtual void get_set_entities(const std::string setname, const Entity_kind kind,
                   const Parallel_type ptype,
                   Entity_ID_List& entids) const override;

  virtual void
  get_set_entities_and_vofs(const std::string setname, const Entity_kind kind,
                            const Parallel_type ptype,
                            Entity_ID_List& entids,
                            std::vector<double>* vofs) const override;


  // Miscellaneous functions
  virtual void write_to_exodus_file(const std::string filename) const override;


 protected:
  // individual versions, if recompute is used
  virtual int
  compute_cell_geometry_(const Entity_ID cellid, double* volume,
                         AmanziGeometry::Point* centroid) const override;
  virtual int compute_face_geometry_(const Entity_ID faceid, double* area,
                         AmanziGeometry::Point* centroid,
                         std::vector<AmanziGeometry::Point>& normals) const override;
                
  // build the cache
  virtual int compute_cell_geometric_quantities_() const override;
  virtual int compute_face_geometric_quantities_() const override;

  // build maps
  void init_maps();


  // get faces of a cell and directions in which it is used - this function
  // is implemented in each mesh framework. The results are cached in
  // the base class

  virtual void
  cell_get_faces_and_dirs_internal_(const Entity_ID cellid,
                                    Entity_ID_List* faceids,
                                    std::vector<int>* face_dirs,
                                    const bool ordered = false) const override;

  // Cells connected to a face - this function is implemented in each
  // mesh framework. The results are cached in the base class

  virtual void
  face_get_cells_internal_(const Entity_ID faceid, const Parallel_type ptype,
                           Entity_ID_List* cellids) const override;


  // edges of a face - this function is implemented in each mesh
  // framework. The results are cached in the base class

  virtual void
  face_get_edges_and_dirs_internal_(const Entity_ID faceid,
                                    Entity_ID_List* edgeids,
                                    std::vector<int>* edge_dirs,
                                    const bool ordered = true) const override;

  // edges of a cell - this function is implemented in each mesh
  // framework. The results are cached in the base class.
  virtual void cell_get_edges_internal_(const Entity_ID cellid,
                                        Entity_ID_List* edgeids) const override;

  // edges and directions of a 2D cell - this function is implemented
  // in each mesh framework. The results are cached in the base class.
  virtual void cell_2D_get_edges_and_dirs_internal_(
    const Entity_ID cellid, Entity_ID_List* edgeids,
    std::vector<int>* edge_dirs) const override;

  // Cache connectivity info.
  virtual void cache_cell_face_info_() const override;


  virtual int build_columns_() const;

 protected:
  bool initialized_;

  std::map<Entity_kind, Entity_ID> num_entities_;
  std::map<Entity_kind, Map_ptr_type> maps_;
  std::map<Set_ID, std::vector<int>> regions_;
  Import_ptr_type exterior_face_importer_;

  std::vector<std::vector<AmanziGeometry::Point>> cell_face_bisectors_;
};

bool
viewMeshLogical(const Mesh& m, std::ostream& os = std::cout);


} // namespace AmanziMesh
} // namespace Amanzi


#endif

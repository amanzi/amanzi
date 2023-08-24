/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, Konstantin Lipnikov, others
*/

//! Implementation of the Mesh interface leveraging MOAB.
#ifndef AMANZI_MESH_MOAB_HH_
#define AMANZI_MESH_MOAB_HH_

// MOAB include files - for some reason they have to be first -
// if not the compiler cannot find MPI_COMM_WORLD

#include <memory>
#include <vector>

#include "moab_mpi.h"
#include "Types.hpp"
#include "Core.hpp"
#include "Range.hpp"
#include "MBTagConventions.hpp"
#include "ReadUtilIface.hpp"
#include "MBParallelConventions.h"
#include "ParallelComm.hpp"

#include "AmanziComm.hh"
#include "Epetra_Map.h"
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "GeometricModel.hh"
#include "Point.hh"

namespace Amanzi {
namespace AmanziMesh {

class Mesh_MOAB : public MeshFramework {
 public:
  // the request_faces and request_edges arguments have to be at the
  // end and not in the middle because if we omit them and specify a
  // pointer argument like gm or verbosity_obj, then there is implicit
  // conversion of the pointer to bool, thereby defeating the intent
  // of the call and making the pointer argument seem NULL. In C++11,
  // we could "delete" the illegal version of the call effectively
  // blocking the implicit conversion.
  Mesh_MOAB(const std::string& filename,
            const Comm_ptr_type& comm,
            const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm = Teuchos::null,
            const Teuchos::RCP<Teuchos::ParameterList>& plist = Teuchos::null);

  ~Mesh_MOAB();

  void update();

  // Get parallel type of entity
  Parallel_kind getEntityPtype(const Entity_kind kind, const Entity_ID entid) const override;


  // Get cell type
  Cell_kind getCellType(const Entity_ID cellid) const override;

  bool hasNodeFaces() const { return false; }

  //
  // General mesh information
  // ------------------------
  //

  // Number of entities of any kind (cell, face, node) and in a
  // particular category (OWNED, GHOST, ALL)
  std::size_t getNumEntities(const Entity_kind kind, const Parallel_kind ptype) const override;


  // Global ID of any entity
  Entity_ID GID(const Entity_ID lid, const Entity_kind kind) const;


  //
  // Mesh Entity Adjacencies
  //------------------------

  // Downward Adjacencies
  //---------------------

  // Get nodes of cell
  // On a distributed mesh, all nodes (OWNED or GHOST) of the cell
  // are returned
  // Nodes are returned in a standard order (Exodus II convention)
  // STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES in 3D
  // For a general polyhedron this will return the nodes in
  // arbitrary order
  // In 2D, the nodes of the polygon will be returned in ccw order
  // consistent with the face normal
  void getCellNodes(const Entity_ID cellid,
                    View_type<const Entity_ID, MemSpace_kind::HOST>& nodeids) const override;


  // Edges and edge directions of a face
  void getFaceEdgesAndDirs(
    const Entity_ID cellid,
    View_type<const Entity_ID, MemSpace_kind::HOST>& edgeids,
    View_type<const Direction_type, MemSpace_kind::HOST>* const dirs = nullptr) const override
  {
    Errors::Message mesg("Edges not implemented in this framework. Use MSTK");
    amanzi_throw(mesg);
  }

  // Get nodes of face
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal
  // In 2D, nfnodes is 2
  void getFaceNodes(const Entity_ID faceid,
                    View_type<const Entity_ID, MemSpace_kind::HOST>& nodeids) const override;


  // Get nodes of edge
  void getEdgeNodes(const Entity_ID edgeid,
                    View_type<const Entity_ID, MemSpace_kind::HOST>& nodes) const override
  {
    Errors::Message msg("Edges not implemented in this framework. Use MSTK");
    amanzi_throw(msg);
  }


  // Upward adjacencies
  //-------------------

  // Cells of type 'ptype' connected to a node
  void getNodeCells(const Entity_ID nodeid,
                    const Parallel_kind ptype,
                    View_type<const Entity_ID, MemSpace_kind::HOST>& cellids) const override;

  // Faces of type 'ptype' connected to a node
  void getNodeFaces(const Entity_ID nodeid,
                    const Parallel_kind ptype,
                    View_type<const Entity_ID, MemSpace_kind::HOST>& faceids) const override;

  // Cells of type 'ptype' connected to an edge
  void getEdgeCells(const Entity_ID edgeid,
                    const Parallel_kind ptype,
                    View_type<const Entity_ID, MemSpace_kind::HOST>& cellids) const override
  {
    Errors::Message msg("Edges not implemented in this framework. Use MSTK");
    amanzi_throw(msg);
  }


  // Same level adjacencies
  //-----------------------

  // Face connected neighboring cells of given cell of a particular ptype
  // (e.g. a hex has 6 face neighbors)

  // The order in which the cellids are returned cannot be
  // guaranteed in general except when ptype = ALL, in which case
  // the cellids will correcpond to cells across the respective
  // faces given by cell_get_faces
  void cell_get_face_adj_cells(const Entity_ID cellid,
                               const Parallel_kind ptype,
                               Entity_ID_View* fadj_cellids) const;

  //---------------------
  // Mesh entity geometry
  //---------------------

  // Node coordinates - 3 in 3D and 2 in 2D
  AmanziGeometry::Point getNodeCoordinate(const Entity_ID nodeid) const override;

  // Face coordinates - conventions same as face_to_nodes call
  // Number of nodes is the vector size divided by number of spatial dimensions
  cPoint_View getFaceCoordinates(const Entity_ID faceid) const override;

  // Coordinates of cells in standard order (Exodus II convention)
  // STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES IN 3D
  // For a general polyhedron this will return the node coordinates in
  // arbitrary order
  // Number of nodes is vector size divided by number of spatial dimensions
  cPoint_View getCellCoordinates(const Entity_ID cellid) const override;

  // Modify the coordinates of a node
  void setNodeCoordinates(const Entity_ID nodeid, const AmanziGeometry::Point& coords);

  //void node_set_coordinates(const Entity_ID nodeid, const double* coords);

  //
  // Boundary Conditions or Sets
  //----------------------------
  virtual bool isValidSetType(const AmanziGeometry::RegionType rtype, const Entity_kind kind) const
  {
    if (rtype == AmanziGeometry::RegionType::BOX || rtype == AmanziGeometry::RegionType::PLANE ||
        rtype == AmanziGeometry::RegionType::POINT ||
        rtype == AmanziGeometry::RegionType::POLYGON ||
        rtype == AmanziGeometry::RegionType::LABELEDSET ||
        rtype == AmanziGeometry::RegionType::COLORFUNCTION) {
      // while some LOGICAL operations are implemented, not all are...
      //rtype == AmanziGeometry::RegionType::LOGICAL) {
      return true;
    } else {
      return false;
    }
  }

  void getSetEntities(const AmanziGeometry::RegionLabeledSet& region,
                      const Entity_kind kind,
		      const Parallel_kind ptype,
		      View_type<const Entity_ID, MemSpace_kind::HOST>& setents) const override; 

  // Deform a mesh so that cell volumes conform as closely as possible
  // to target volumes without dropping below the minimum volumes.  If
  // move_vertical = true, nodes will be allowed to move only in the
  // vertical direction (right now arbitrary node movement is not allowed)
  int deform(const Double_List& target_cell_volumes__in,
             const Double_List& min_cell_volumes__in,
             const Entity_ID_View& fixed_nodes,
             const bool move_vertical);


  //
  // Miscellaneous
  // -------------

  void write_to_exodus_file(const std::string filename) const;

 private:
  void ErrorCheck_(int result, std::string msg) const
  {
    if (result != moab::MB_SUCCESS) {
      std::cerr << msg << std::endl;
      Errors::Message err(msg);
      Exceptions::amanzi_throw(err);
    }
  }

 private:
  std::unique_ptr<moab::Core> mbcore_;
  std::unique_ptr<moab::ParallelComm> mbcomm_;

  int serial_run;

  int celldim; // Topological dimension of highest level entities
  int facedim; // Topological dimension of 2nd highest level entities

  // Local handles to entity lists (Vertices, "Faces", "Cells")

  // For a surface mesh, "Faces" refers to mesh edges and "Cells"
  // refers to mesh faces
  //
  // For a solid mesh, "Faces" refers to mesh faces and "Cells"
  // refers to mesh regions

  // These are MOAB's definitions of types of parallel mesh entities
  // These definitions are slightly different from what Amanzi has defined
  //
  // There are 2 types of cells - Owned and Ghost

  // There are 4 types of lower dimensional entities
  //
  // 1. OWNED - owned by this processor
  // 2. NOTOWNED - not owned by this processor
  // 3. USED - connected to at least one cell owned by this
  // processor (may or may not be owned by this processor)
  // 4. GHOST - neither the entity nor a cell connected to the
  // entity is owned by the processor

  // UNFORTUNATELY, THE TERMINOLOGY USED BY THE API USES GHOST TO
  // MEAN NOTOWNED

  // ALL = OWNED + NOTOWNED or USED + GHOST

  moab::Range AllVerts, OwnedVerts, NotOwnedVerts;
  moab::Range AllFaces, OwnedFaces, NotOwnedFaces;
  moab::Range AllCells, OwnedCells, GhostCells;

  // tag handles
  moab::Tag lid_tag; // Local ID
  moab::Tag gid_tag; // Global ID
  moab::Tag cstag;   // Material tag
  moab::Tag sstag;   // Sideset tag
  moab::Tag nstag;   // Nodeset tag

  // Local ID to MOAB handle map
  std::vector<moab::EntityHandle> node_id_to_handle;
  std::vector<moab::EntityHandle> face_id_to_handle;
  std::vector<moab::EntityHandle> cell_id_to_handle;

  // Maps
  Epetra_Map *cell_map_wo_ghosts_, *face_map_wo_ghosts_, *node_map_wo_ghosts_;
  Epetra_Map *cell_map_w_ghosts_, *face_map_w_ghosts_, *node_map_w_ghosts_;
  Epetra_Map *extface_map_w_ghosts_,
    *extface_map_wo_ghosts_; // exterior faces (connected to only 1 cell)

  // Sets (material sets, sidesets, nodesets)
  // We store the number of sets in the whole problem regardless of whether
  // they are represented on this processor or not
  // We also store the IDs of the sets and the dimension of entities
  // in those sets
  int nsets;
  int *setids_, *setdims_;

  // Minimum and maximum global IDs of faces
  unsigned int minFGID, maxFGID;

  // flag whether to flip a face dir or not when returning nodes of a face
  bool* faceflip;

  // Private methods
  // ----------------------------
  bool valid_entity_kind_(int kind) const;

  void clear_internals_();

  void init_pvert_lists();
  void init_pface_lists();
  void init_pcell_lists();
  void init_pface_dirs();

  void init_id_handle_maps();
  void init_global_ids();

  void init_set_info();

  std::string internal_name_of_set(const Teuchos::RCP<const AmanziGeometry::Region>& r,
                                   const Entity_kind entity_kind) const;

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
  void
  getCellFacesAndDirs(const Entity_ID cellid,
                      View_type<const Entity_ID, MemSpace_kind::HOST>& faceids,
                      View_type<const Direction_type, MemSpace_kind::HOST>* const face_dirs) const;

  // Cells connected to a face
  void getFaceCells(const Entity_ID faceid,
                    const Parallel_kind ptype,
                    View_type<const Entity_ID, MemSpace_kind::HOST>& cellids) const;

  // Edges of a cell
  void getCellEdges(const Entity_ID cellid,
                    View_type<const Entity_ID, MemSpace_kind::HOST>& edgeids) const override
  {
    Errors::Message mesg("Edges not implemented in this framework. Use MSTK");
    Exceptions::amanzi_throw(mesg);
  }

  // Edges and edge directions of a face
  /*
  void face_get_edges_and_dirs_internal_(const Entity_ID cellid,
                                         Entity_ID_View *edgeids,
                                         std::vector<int> *edgedirs,
                                         bool ordered=true) const
  {
    Errors::Message mesg("Edges not implemented in this framework. Use MSTK");
    amanzi_throw(mesg);
  }
  */
};

} // namespace AmanziMesh
} // namespace Amanzi

#endif

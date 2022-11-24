/*
  Mesh

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Implementation of simple mesh.
*/

#ifndef AMANZI_MESH_SIMPLE_HH_
#define AMANZI_MESH_SIMPLE_HH_

#include <Epetra_Map.h>
#include <AmanziComm.hh>

#include <memory>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Mesh.hh"
#include "Region.hh"

#include "GeometricModel.hh"
#include "VerboseObject.hh"

namespace Amanzi {
namespace AmanziMesh {

class Mesh_simple : public Mesh {
 public:
  // the request_faces and request_edges arguments have to be at the
  // end and not in the middle because if we omit them and specify a
  // pointer argument like gm or verbosity_obj, then there is implicit
  // conversion of the pointer to bool, thereby defeating the intent
  // of the call and making the pointer argument seem NULL. In C++11,
  // we could "delete" the illegal version of the call effectively
  // blocking the implicit conversion.
  Mesh_simple(double x0,
              double y0,
              double z0,
              double x1,
              double y1,
              double z1,
              int nx,
              int ny,
              int nz,
              const Comm_ptr_type& comm,
              const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm = Teuchos::null,
              const Teuchos::RCP<const Teuchos::ParameterList>& plist = Teuchos::null,
              const bool request_faces = true,
              const bool request_edges = false);

  virtual ~Mesh_simple() = default;

  // Get parallel type of entity
  Parallel_type entity_get_ptype(const Entity_kind kind, const Entity_ID entid) const override
  {
    return Parallel_type::OWNED;
  }

  // Get cell type
  Cell_type cell_get_type(const Entity_ID cellid) const override { return HEX; }

  // Number of entities of any kind (cell, face, node) and in a
  // particular category (OWNED, GHOST, ALL)
  unsigned int num_entities(const Entity_kind kind, const Parallel_type ptype) const override;

  // Global ID of any entity (this is a serial code)
  Entity_ID GID(const Entity_ID lid, const Entity_kind kind) const override { return lid; }


  //---------------------
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
  void cell_get_nodes(const Entity_ID cellid, std::vector<Entity_ID>* nodeids) const override;

  // Get nodes of face
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal
  // In 2D, nfnodes is 2
  void face_get_nodes(const Entity_ID faceid, std::vector<Entity_ID>* nodeids) const override;

  // Get nodes of edge
  void
  edge_get_nodes(const Entity_ID edgeid, Entity_ID* nodeid0, Entity_ID* nodeid1) const override;

  //-------------------
  // Upward adjacencies
  //-------------------

  // Cells of type 'ptype' connected to a node
  void node_get_cells(const Entity_ID nodeid,
                      const Parallel_type ptype,
                      std::vector<Entity_ID>* cellids) const override;

  // Faces of type 'ptype' connected to a node
  void node_get_faces(const Entity_ID nodeid,
                      const Parallel_type ptype,
                      std::vector<Entity_ID>* faceids) const override;

  // Cells of type 'ptype' connected to an edge
  void edge_get_cells(const Entity_ID edgeid,
                      const Parallel_type ptype,
                      std::vector<Entity_ID>* cellids) const override
  {
    Errors::Message msg("Edge to cell connectivity is not implemented in this framework.");
    amanzi_throw(msg);
  }


  //-----------------------
  // Same level adjacencies
  //-----------------------

  // Face connected neighboring cells of given cell of a particular ptype
  // (e.g. a hex has 6 face neighbors)

  // The order in which the cellids are returned cannot be
  // guaranteed in general except when ptype = ALL, in which case
  // the cellids will correcpond to cells across the respective
  // faces given by cell_get_faces
  void cell_get_face_adj_cells(const Entity_ID cellid,
                               const Parallel_type ptype,
                               std::vector<Entity_ID>* fadj_cellids) const override;

  //---------------------
  // Mesh entity geometry
  //---------------------

  // Node coordinates - 3 in 3D and 2 in 2D
  void node_get_coordinates(const Entity_ID nodeid, AmanziGeometry::Point* ncoord) const override;

  // Modify the coordinates of a node
  void node_set_coordinates(const Entity_ID nodeid, const AmanziGeometry::Point coords) override;
  void node_set_coordinates(const Entity_ID nodeid, const double* coords) override;


  //------------
  // Epetra maps
  //------------
  const Epetra_Map& cell_map(bool include_ghost) const override { return *cell_map_; }
  const Epetra_Map& face_map(bool include_ghost) const override { return *face_map_; }
  const Epetra_Map& edge_map(bool include_ghost) const override { return *edge_map_; }
  const Epetra_Map& node_map(bool include_ghost) const override { return *node_map_; }

  const Epetra_Map& exterior_face_map(bool include_ghost) const override { return *extface_map_; }
  const Epetra_Map& exterior_node_map(bool include_ghost) const override;

  // Epetra importer that will allow apps to import values from a
  // Epetra vector defined on all owned faces into an Epetra vector
  // defined only on exterior faces
  const Epetra_Import& exterior_face_importer(void) const override;


  //----------------------------
  // Boundary Conditions or Sets
  //----------------------------
  virtual bool
  valid_set_type(const AmanziGeometry::RegionType rtype, const Entity_kind kind) const override
  {
    if (rtype == AmanziGeometry::RegionType::POINT && kind == Entity_kind::NODE) {
      return true;
    } else if (rtype == AmanziGeometry::RegionType::BOX ||
               rtype == AmanziGeometry::RegionType::PLANE ||
               rtype == AmanziGeometry::RegionType::POLYGON ||
               rtype == AmanziGeometry::RegionType::LABELEDSET ||
               rtype == AmanziGeometry::RegionType::COLORFUNCTION) {
      return true;
    } else {
      return false;
    }
  }

  // Get list of entities of type 'category' in set
  virtual void get_set_entities_and_vofs(const std::string& setname,
                                         const Entity_kind kind,
                                         const Parallel_type ptype,
                                         Entity_ID_List* entids,
                                         std::vector<double>* vofs) const override;


  // Miscellaneous
  void write_to_exodus_file(const std::string filename) const override;


  // this should be used with extreme caution:
  // modify coordinates
  void set_coordinate(Entity_ID local_node_id, double* source_begin, double* source_end);


  // Deform a mesh so that cell volumes conform as closely as possible
  // to target volumes without dropping below the minimum volumes.  If
  // move_vertical = true, nodes will be allowed to move only in the
  // vertical direction (right now arbitrary node movement is not allowed)
  int deform(const std::vector<double>& target_cell_volumes__in,
             const std::vector<double>& min_cell_volumes__in,
             const Entity_ID_List& fixed_nodes,
             const bool move_vertical) override;

 private:
  void CreateCache_();
  void BuildMaps_();

  AmanziGeometry::Point geometric_center_(std::vector<AmanziGeometry::Point>& vxyz) const;

  bool
  all_inside_(std::vector<AmanziGeometry::Point>& vxyz, const AmanziGeometry::Region& rgn) const;

  Teuchos::RCP<Epetra_Map> cell_map_, face_map_, edge_map_, node_map_;
  Teuchos::RCP<Epetra_Map> extface_map_;

  std::vector<double> coordinates_;

  unsigned int node_index_(int i, int j, int k) const;
  unsigned int cell_index_(int i, int j, int k) const;

  unsigned int xyface_index_(int i, int j, int k) const;
  unsigned int yzface_index_(int i, int j, int k) const;
  unsigned int xzface_index_(int i, int j, int k) const;

  unsigned int xedge_index_(int i, int j, int k) const;
  unsigned int yedge_index_(int i, int j, int k) const;
  unsigned int zedge_index_(int i, int j, int k) const;

  int nx_, ny_, nz_; // number of cells in the three coordinate directions
  double x0_, x1_, y0_, y1_, z0_,
    z1_; // coordinates of lower left front and upper right back of brick

  int num_cells_, num_faces_, num_edges_, num_nodes_;
  int num_faces_bnd_;

  // mesh connectivity arrays
  std::vector<Entity_ID> cell_to_face_;
  std::vector<Entity_ID> cell_to_edge_;
  std::vector<Entity_ID> cell_to_node_;
  std::vector<Entity_ID> face_to_edge_;
  std::vector<Entity_ID> face_to_node_;
  std::vector<Entity_ID> face_to_cell_;
  std::vector<Entity_ID> edge_to_node_;
  std::vector<Entity_ID> node_to_face_;
  std::vector<Entity_ID> node_to_cell_;

  // orientation arrays
  std::vector<int> cell_to_face_dirs_;
  std::vector<int> face_to_edge_dirs_;

  // The following are mutable because they have to be modified
  // after the class construction even though the class is instantiated
  // as a constant class
  mutable std::map<std::string, std::vector<Entity_ID>> sets_;

  // Get faces of a cell and directions in which the cell uses the face

  // In 3D, direction is 1 if face normal points out of cell
  // and -1 if face normal points into cell
  // In 2D, direction is 1 if face/edge is defined in the same
  // direction as the cell polygon, and -1 otherwise
  void cell_get_faces_and_dirs_internal_(const Entity_ID cellid,
                                         Entity_ID_List* faceids,
                                         std::vector<int>* face_dirs,
                                         const bool ordered = false) const override;

  // Cells connected to a face
  void face_get_cells_internal_(const Entity_ID faceid,
                                const Parallel_type ptype,
                                std::vector<Entity_ID>* cellids) const override;


  // Edges of a cell
  void cell_get_edges_internal_(const Entity_ID cellid, Entity_ID_List* edgeids) const override;

  // Get edges of a face and directions in which the face uses the edges.
  // In 3D, edge direction is 1 when it is oriented counter clockwise
  // with respect to the face natural normal.
  void face_get_edges_and_dirs_internal_(const Entity_ID cellid,
                                         Entity_ID_List* edgeids,
                                         std::vector<int>* edgedirs,
                                         bool ordered = true) const override;
};


// -------------------------
// Template & inline members
// ------------------------
inline unsigned int
Mesh_simple::node_index_(int i, int j, int k) const
{
  return i + j * (nx_ + 1) + k * (nx_ + 1) * (ny_ + 1);
}

inline unsigned int
Mesh_simple::cell_index_(int i, int j, int k) const
{
  return i + j * nx_ + k * nx_ * ny_;
}

inline unsigned int
Mesh_simple::xyface_index_(int i, int j, int k) const
{
  return i + j * nx_ + k * nx_ * ny_;
}

inline unsigned int
Mesh_simple::xzface_index_(int i, int j, int k) const
{
  return i + j * nx_ + k * nx_ * (ny_ + 1) + xyface_index_(0, 0, nz_ + 1);
}

inline unsigned int
Mesh_simple::yzface_index_(int i, int j, int k) const
{
  return i + j * (nx_ + 1) + k * (nx_ + 1) * ny_ + xzface_index_(0, 0, nz_);
}

inline unsigned int
Mesh_simple::xedge_index_(int i, int j, int k) const
{
  return i + j * nx_ + k * nx_ * (ny_ + 1);
}

inline unsigned int
Mesh_simple::yedge_index_(int i, int j, int k) const
{
  return i + j * (nx_ + 1) + k * (nx_ + 1) * ny_ + xedge_index_(0, 0, nz_ + 1);
}

inline unsigned int
Mesh_simple::zedge_index_(int i, int j, int k) const
{
  return i + j * (nx_ + 1) + k * (nx_ + 1) * (ny_ + 1) + yedge_index_(0, 0, nz_ + 1);
}

} // namespace AmanziMesh
} // namespace Amanzi

#endif

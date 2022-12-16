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
  1. nodes do not exist.
  2. cells are defined by a collection of cell-face bisectors.
  3. volumes can be defined by areas times the length of those bisectors

*/

#ifndef AMANZI_LOGICAL_MESH_H_
#define AMANZI_LOGICAL_MESH_H_

#include <memory>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Epetra_Map.h"
#include "AmanziComm.hh"
#include "Epetra_SerialComm.h"

#include "VerboseObject.hh"
#include "dbc.hh"
#include "errors.hh"

#include "Region.hh"
#include "Mesh.hh"

namespace Amanzi {
namespace AmanziMesh {


struct MeshLogicalAlgorithms : public MeshFrameworkAlgorithms {
  // lumped things for more efficient calculation
  virtual std::pair<double, AmanziGeometry::Point>
  computeCellGeometry(const MeshFramework& mesh, const Entity_ID c) const override;


  virtual std::tuple<double, AmanziGeometry::Point, Point_List>
  computeFaceGeometry(const MeshFramework& mesh, const Entity_ID f) const override;


  virtual std::pair<AmanziGeometry::Point, AmanziGeometry::Point>
  computeEdgeGeometry(const MeshFramework& mesh, const Entity_ID e) const override;

};


class MeshLogical : public MeshFramework {
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

  //
  // Topology only constructor
  // -----------------------------------------------------------------------------
  MeshLogical(const Comm_ptr_type& comm,
              const std::vector<Entity_ID_List>& face_cell_ids,
              const std::vector<std::vector<int> >& face_cell_dirs,
              const Teuchos::RCP<Teuchos::ParameterList>& plist=Teuchos::null);

  //
  // Topology and geometry constructor
  // -----------------------------------------------------------------------------
  MeshLogical(const Comm_ptr_type& comm,
              const std::vector<Entity_ID_List>& face_cell_ids,
              const std::vector<std::vector<int> >& face_cell_dirs,
              const std::vector<double>& cell_volumes,
              const std::vector<double>& face_areas,
              const std::vector<std::vector<AmanziGeometry::Point> >& face_cell_bisectors,
              const std::vector<AmanziGeometry::Point>* cell_centroids=nullptr,
              const Teuchos::RCP<Teuchos::ParameterList>& plist=Teuchos::null);


  void getLogicalGeometry(std::vector<double>* const cell_volumes,
                            std::vector<double>* const face_areas,
                            std::vector<std::vector<AmanziGeometry::Point> >* const face_cell_bisectors,
                            std::vector<AmanziGeometry::Point>* const cell_centroids) const;

  void setLogicalGeometry(std::vector<double> const* const cell_volumes,
                            std::vector<double> const* const face_areas,
                            std::vector<std::vector<AmanziGeometry::Point> > const* const face_cell_bisectors,
                            std::vector<AmanziGeometry::Point> const* const cell_centroids=nullptr);


  // for testing
  bool operator==(const MeshLogical& other);

  // disallow nodes
  virtual bool hasNodes() const override { return false; }

  // Some meshes are logical meshes and do not have coordinate info.
  virtual bool isLogical() const override { return true; }

  // Get cell type - UNKNOWN, TRI, QUAD, ... See MeshDefs.hh
  virtual
  Cell_type getCellType(const Entity_ID cellid) const override {
    return Cell_type::UNKNOWN;
  }

  // Number of entities of any kind (cell, face, node) and in a
  // particular category (OWNED, GHOST, ALL)
  virtual
  std::size_t getNumEntities(const Entity_kind kind,
          const Parallel_type ptype) const override;


  // All nodal methods throw -- there are no nodes in MeshLogical
  virtual AmanziGeometry::Point
  getNodeCoordinate(const Entity_ID node) const override;

  virtual void
  getFaceNodes(const Entity_ID f, Entity_ID_List& nodes) const override;

  virtual void
  getNodeFaces(const Entity_ID nodeid,
               const Parallel_type ptype,
               Entity_ID_List& faceids) const override;

  //
  // These are the important ones -- MeshLogical defines cell quantities
  //
  virtual void getCellFacesAndDirs(
    const Entity_ID c,
    Entity_ID_List& faces,
    Entity_Direction_List * const dirs) const override;

  // Get the bisectors, i.e. vectors from cell centroid to face centroids.
  virtual void getCellFacesAndBisectors(
          const Entity_ID cellid,
          Entity_ID_List& faceids,
          Point_List * const bisectors) const override;

  virtual double getCellVolume(const Entity_ID c) const override;
  virtual AmanziGeometry::Point getCellCentroid(const Entity_ID c) const override;

  //
  // MeshLogical defines face quantities
  //
  virtual void getFaceCells(const Entity_ID f,
                            const Parallel_type ptype,
                            Entity_ID_List& cells) const override;
  virtual double getFaceArea(const Entity_ID f) const override;
  virtual AmanziGeometry::Point getFaceCentroid(const Entity_ID f) const override;
  virtual AmanziGeometry::Point getFaceNormal(const Entity_ID f,
          const Entity_ID c, int * const orientation=nullptr) const override;

 protected:
  bool initialized_;

  Double_List cell_volumes_;
  Double_List face_areas_;
  Point_List cell_centroids_;
  RaggedArray_DualView<Entity_ID> face_cell_ids_;

  RaggedArray_DualView<Entity_ID> cell_face_ids_;
  RaggedArray_DualView<int> cell_face_dirs_;
  RaggedArray_DualView<AmanziGeometry::Point> cell_face_bisectors_;

};

bool viewMeshLogical(const Mesh& m, std::ostream& os=std::cout);

} // namespace AmanziMesh
} // namespace Amanzi

#endif

/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
*/

/*!

 Virtual mesh that can be modified and constructed on the fly.

 Virtual mesh is a logical mesh, i.e. it is a topologically defined mesh with no minimum 
 geometric data. It is designed for subgrid models and may work with a few spatial 
 discretizations, including MFD and FV.

 Unlike a logical mesh, a virtual mesh uses face normals which may not be aliged with 
 cell-face bisectors.

 Assumptions:
  1. nodes do not exist
  2. cells are defined by their centroids and volumes
  3. faces are defined by their centroids and total normals
  4. connectivy list (face -> cells) ordered from small to large cell index
  5. default orientation of face normal is also from small to large cell index
  6. face normal are exterior on the boundary

*/

#ifndef AMANZI_VIRTUAL_MESH_HH_
#define AMANZI_VIRTUAL_MESH_HH_

#include <memory>
#include <vector>

#include "Teuchos_ParameterList.hpp"

#include "AmanziComm.hh"
#include "Mesh.hh"
#include "MeshAlgorithms.hh"
#include "MeshFramework.hh"
#include "ViewUtils.hh"

#include "dbc.hh"
#include "errors.hh"

namespace Amanzi {
namespace AmanziMesh {

struct MeshVirtualAlgorithms : public MeshAlgorithms {
  virtual std::pair<double, AmanziGeometry::Point>
  computeCellGeometry(const Mesh& mesh, const Entity_ID c) const override;

  virtual std::tuple<double, AmanziGeometry::Point, cPoint_View>
  computeFaceGeometry(const Mesh& mesh, const Entity_ID f) const override;

  virtual std::pair<AmanziGeometry::Point, AmanziGeometry::Point>
  computeEdgeGeometry(const Mesh& mesh, const Entity_ID e) const override;

  virtual double computeCellVolume(const Mesh& mesh, const Entity_ID c) const override;
  virtual AmanziGeometry::Point
  computeCellCentroid(const Mesh& mesh, const Entity_ID c) const override;

  virtual double computeFaceArea(const Mesh& mesh, const Entity_ID f) const override;
  virtual AmanziGeometry::Point
  computeFaceCentroid(const Mesh& mesh, const Entity_ID f) const override;
  virtual AmanziGeometry::Point computeFaceNormal(const Mesh& mesh,
                                                  const Entity_ID f,
                                                  const Entity_ID c,
                                                  int* const orientation = nullptr) const override;
};


class MeshVirtual : public MeshFramework {
 public:
  MeshVirtual(const Comm_ptr_type& comm,
              const Teuchos::RCP<Teuchos::ParameterList>& plist,
              const std::vector<Entity_ID_List>& face_cells,
              const std::vector<AmanziGeometry::Point>& cell_centroids,
              const std::vector<double>& cell_volumes,
              const std::vector<AmanziGeometry::Point>& face_centroids,
              const std::vector<AmanziGeometry::Point>& face_normals);

  // General mesh information
  virtual bool hasNodes() const override { return false; }

  virtual Parallel_kind
  getEntityPtype(const Entity_kind kind, const Entity_ID entid) const override
  {
    return Parallel_kind::OWNED;
  }

  virtual Entity_ID
  getEntityParent(const Entity_kind kind, const Entity_ID entid) const override { return -1; }

  Cell_kind getCellType(const Entity_ID c) const override { return Cell_kind::UNKNOWN; }

  virtual std::size_t
  getNumEntities(const Entity_kind kind, const Parallel_kind ptype) const override
  {
    return num_entities_.at(kind);
  }


  // Downward Adjacencies
  virtual void
  getCellFacesAndDirs(const Entity_ID c,
                      cEntity_ID_View& faces,
                      cDirection_View* const dirs) const override;

  virtual void getCellNodes(const Entity_ID c, cEntity_ID_View& nodes) const override;
  virtual void getFaceNodes(const Entity_ID f, cEntity_ID_View& nodes) const override;


  // Upward adjacencies
  virtual void getNodeFaces(const Entity_ID n, cEntity_ID_View& faces) const override;
  virtual void getFaceCells(const Entity_ID f, cEntity_ID_View& cells) const override;


  // Mesh entity geometry
  virtual AmanziGeometry::Point getNodeCoordinate(const Entity_ID n) const override;

  AmanziGeometry::Point getFaceCentroid(const Entity_ID f) const { return face_centroids_[f]; }
  AmanziGeometry::Point getFaceNormal(const Entity_ID f, const Entity_ID c, int* orientation) const;
  double getFaceArea(const Entity_ID f) const { return norm(face_normals_[f]); }

  AmanziGeometry::Point getCellCentroid(const Entity_ID c) const { return cell_centroids_[c]; }
  double getCellVolume(const Entity_ID c) const { return cell_volumes_[c]; }

  virtual bool
  isValidSetType(const AmanziGeometry::RegionType rtype, const Entity_kind kind) const override
  {
    return (rtype == AmanziGeometry::RegionType::ALL) ||
           (rtype == AmanziGeometry::RegionType::ENUMERATED);
  }

 protected:
  const std::vector<AmanziGeometry::Point>& cell_centroids_;
  const std::vector<double> cell_volumes_;
  const std::vector<AmanziGeometry::Point> face_centroids_;
  const std::vector<AmanziGeometry::Point> face_normals_;

  std::map<Entity_kind, Entity_ID> num_entities_;

  RaggedArray_DualView<Entity_ID> face_cell_ids_;
  RaggedArray_DualView<Entity_ID> cell_face_ids_;
  RaggedArray_DualView<int> cell_face_dirs_;
};

} // namespace AmanziMesh
} // namespace Amanzi

#endif

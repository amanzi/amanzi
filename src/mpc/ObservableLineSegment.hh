/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef AMANZI_OBSERVABLE_LINE_SEGMENT_HH
#define AMANZI_OBSERVABLE_LINE_SEGMENT_HH

#include "ObservableAmanzi.hh"
#include "RegionPlane.hh"
#include "RegionPolygon.hh"
#include "RegionLineSegment.hh"
#include "Units.hh"

namespace Amanzi {

class ObservableLineSegment : public virtual Observable {
 public:
  ObservableLineSegment(std::string variable,
                        std::string region,
                        std::string functional,
                        Teuchos::ParameterList& plist,
                        Teuchos::ParameterList& units_plist,
                        Teuchos::RCP<const AmanziMesh::Mesh> mesh);

  virtual void
  ComputeObservation(State& S, double* value, double* volume, std::string& unit, double dt);
  virtual int ComputeRegionSize();
  void ComputeInterpolationPoints(Teuchos::RCP<const AmanziGeometry::Region> reg_ptr);

 protected:
  AmanziMesh::Double_View lofs_;
  std::string interpolation_;
  std::string weighting_;
  std::vector<AmanziGeometry::Point> line_points_;
  bool limiter_;
};


ObservableLineSegment::ObservableLineSegment(std::string variable,
                                             std::string region,
                                             std::string functional,
                                             Teuchos::ParameterList& plist,
                                             Teuchos::ParameterList& units_plist,
                                             Teuchos::RCP<const AmanziMesh::Mesh> mesh)
  : Observable(variable, region, functional, plist, units_plist, mesh)
{
  interpolation_ = plist.get<std::string>("interpolation", "linear");
  weighting_ = plist.get<std::string>("weighting", "none");
  limiter_ = plist.get<bool>("limiter", true);
};


int
ObservableLineSegment::ComputeRegionSize()
{
  //int mesh_block_size;
  Errors::Message msg;

  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm_ptr = mesh_->getGeometricModel();
  Teuchos::RCP<const AmanziGeometry::Region> reg_ptr = gm_ptr->FindRegion(region_);

  if (reg_ptr->get_type() != AmanziGeometry::RegionType::LINE_SEGMENT) {
    msg << "ObservableLineSegment works only with LineSegment region";
    Exceptions::amanzi_throw(msg);
  }

  // all others need cells
  region_size_ =
    mesh_->getSetSize(region_, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  Kokkos::resize(entity_ids_, region_size_);
  std::tie(entity_ids_, lofs_) = mesh_->getSetEntitiesAndVolumeFractions(
    region_, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  ComputeInterpolationPoints(reg_ptr);

  // find global meshblocksize
  int dummy = region_size_;
  int global_mesh_block_size(0);
  mesh_->getComm()->SumAll(&dummy, &global_mesh_block_size, 1);

  return global_mesh_block_size;
}


void
ObservableLineSegment::ComputeObservation(State& S,
                                          double* value,
                                          double* volume,
                                          std::string& unit,
                                          double dt)
{
  Errors::Message msg;
  msg << "Observation should be computed in classes inheritated from ObservableLineSegment\n";
  Exceptions::amanzi_throw(msg);
}


void
ObservableLineSegment::ComputeInterpolationPoints(
  Teuchos::RCP<const AmanziGeometry::Region> reg_ptr)
{
  std::vector<std::vector<int>> polytope_faces;

  line_points_.resize(entity_ids_.size());

  Teuchos::RCP<const AmanziGeometry::RegionLineSegment> line_segment_ptr =
    Teuchos::rcp_dynamic_cast<const AmanziGeometry::RegionLineSegment>(reg_ptr);

  // double sum = 0.0;

  for (int k = 0; k < entity_ids_.size(); k++) {
    int c = entity_ids_[k];
    auto polytope_nodes = mesh_->getCellCoordinates(c);

    if (mesh_->getSpaceDimension() == 3) {
      auto cnodes = mesh_->getCellNodes(c);
      auto [faces, dirs] = mesh_->getCellFacesAndDirections(c);
      int nfaces = faces.size();

      polytope_faces.clear();
      polytope_faces.resize(nfaces);
      for (int n = 0; n < nfaces; ++n) {
        auto fnodes = mesh_->getFaceNodes(faces[n]);
        int nnodes = fnodes.size();

        for (int i = 0; i < nnodes; ++i) {
          int j = (dirs[n] > 0) ? i : nnodes - i - 1;
          // Find node position of fnodes[j]
          int node_pos = -1;
          for (int m = 0; m < cnodes.size(); m++) {
            if (cnodes[m] == fnodes[j]) {
              node_pos = m;
              break;
            }
          }
          polytope_faces[n].push_back(node_pos);
        }
      }

      line_segment_ptr->ComputeInterLinePoints(
        asVector(polytope_nodes), polytope_faces, line_points_[k]);
      // sum += lofs_[k];
    }
  }
}

} // namespace Amanzi

#endif

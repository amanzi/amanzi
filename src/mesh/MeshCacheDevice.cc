#include "Teuchos_CommHelpers.hpp"

#include "Geometry.hh"
#include "Iterators.hh"
#include "MeshUtils.hh"
#include "MeshFramework.hh"
#include "MeshCacheDevice_impl.hh"
#include "MeshAlgorithms.hh"
#include "MeshSets.hh"

namespace Amanzi::AmanziMesh {

MeshCacheDevice::MeshCacheDevice(const MeshCacheHost& other) : MeshCacheBase(other)
{
  parent_ = Teuchos::RCP<const MeshCacheDevice>(
    other.getParentMesh() == Teuchos::null ? Teuchos::null : onMemDevice(other.getParentMesh()));
}

//
// TODO --etc
// This should be updated -- only cache the ALL set, then construct the
// OWNED or GHOST set on demand.  No need to save both OWNED and ALL.
//
// TODO --etc
// Don't cache "ALL" or "BOUNDARY" labeled sets.  Are there others that are
// just as fast to create as to cache?
//
typename MeshCacheDevice::cEntity_ID_View
MeshCacheDevice::getSetEntities(const std::string& region_name,
                                const Entity_kind kind,
                                const Parallel_kind ptype) const
{
  auto key = std::make_tuple(region_name, kind, ptype);
  auto key_all = std::make_tuple(region_name, kind, Parallel_kind::ALL);
  bool check_new = false;

  if (!sets_.count(key_all)) {
    auto region = getGeometricModel()->FindRegion(region_name);
    if (region == Teuchos::null) {
      Errors::Message msg;
      msg << "Cannot find region of name \"" << region_name << "\" in the geometric model.";
      Exceptions::amanzi_throw(msg);
    }

    auto this_on_host = onMemHost(*this);
    sets_[key_all] = asDualView(resolveMeshSet(*region, kind, Parallel_kind::ALL, this_on_host));
    check_new = true;
  }

  if (!sets_.count(key)) {
    if (ptype == Parallel_kind::OWNED) {
      auto v_all = view<MemSpace_kind::HOST>(sets_.at(key_all));
      size_t n_ents = getNumEntities(kind, Parallel_kind::OWNED);
      size_t i;
      for (i = 0; i != v_all.size(); ++i)
        if (v_all(i) >= n_ents) break;
      auto v_owned = Kokkos::subview(v_all, Kokkos::make_pair((size_t)0, i));
      sets_[key] = asDualView(v_owned);

    } else if (ptype == Parallel_kind::GHOST) {
      auto v_all = view<MemSpace_kind::HOST>(sets_.at(key_all));
      size_t n_ents = getNumEntities(kind, Parallel_kind::OWNED);
      size_t i;
      for (i = 0; i != v_all.size(); ++i)
        if (v_all(i) >= n_ents) break;
      auto v_ghosted = Kokkos::subview(v_all, Kokkos::make_pair(i, v_all.size()));
      sets_[key] = asDualView(v_ghosted);
    } else {
      AMANZI_ASSERT(false);
    }
    check_new = true;
  }

  if (check_new) {
    //   auto region = getGeometricModel()->FindRegion(region_name);
    //   // Error on zero -- some zeros already error internally (at the framework
    //   // level) but others don't.  This is the highest level we can catch these at.
    //   int lsize = view<MEM>(sets_.at(key)).size();
    //   int gsize = 0;
    //   getComm()->SumAll(&lsize, &gsize, 1);
    //   if (gsize == 0 && getComm()->MyPID() == 0) {
    //     Errors::Message msg;
    //     msg << "AmanziMesh::getSetEntities: Region \"" << region->get_name() << "\" of type \""
    //         << to_string(region->get_type()) << "\" and kind \"" << to_string(kind)
    //         << "\" is empty (globally).\n";
    //     std::cout << msg.what();
    //     // Exceptions::amanzi_throw(msg);
    //   }
  }

  return view<MEM>(sets_.at(key));
}

bool
MeshCacheDevice::isValidSetType(const AmanziGeometry::RegionType rtype,
                                const Entity_kind kind) const
{
  if (rtype == AmanziGeometry::RegionType::LABELEDSET && framework_mesh_.get()) {
    return framework_mesh_->isValidSetType(rtype, kind);
  }
  return true;
}


int
MeshCacheDevice::getSetSize(const std::string& region_name,
                            const Entity_kind kind,
                            const Parallel_kind ptype) const
{
  return getSetEntities(region_name, kind, ptype).size();
}


//===================
//    getFace*
//===================

// // Note that regions are cached on demand the first time they are requested,
// // but labeled sets must be pre-cached if the framework mesh is to be
// // destroyed.
// void MeshCacheDeviceprecacheLabeledSets()
// {
//   for (const auto& rgn : *getGeometricModel()) {
//     if (rgn->get_type() == AmanziGeometry::RegionType::LABELEDSET) {
//       auto rgn_lbl = Teuchos::rcp_dynamic_cast<const AmanziGeometry::RegionLabeledSet>(rgn);
//       AMANZI_ASSERT(rgn_lbl);

//       if (getParentMesh() != Teuchos::null) {
//         AMANZI_ASSERT(false); // not yet implemented lifted sets
//       } else {
//         auto entity_kind = createEntityKind(rgn_lbl->entity_str());
//         Entity_ID_View set_ents;
//         framework_mesh_->getSetEntities(*rgn_lbl, entity_kind, Parallel_kind::ALL, set_ents);
//         auto key = std::make_tuple(rgn->get_name(), entity_kind, Parallel_kind::ALL);
//         sets_[key] = asDualView(set_ents);
//       }
//     }
//   }
// }

// build the MeshMaps object
// void MeshCacheDevice::buildMaps_(bool natural)
// {
// }


// // common error messaging
// void MeshCacheDevicethrowAccessError_(const std::string& func_name) const
// {
//   Errors::Message msg;
//   msg << "MeshCacheDevice" << func_name << " cannot compute this quantity -- not cached and framework does not exist.";
//   Exceptions::amanzi_throw(msg);
// }

bool
MeshCacheDevice::isPointInCell(const AmanziGeometry::Point& p, const Entity_ID cellid) const
{
  if (manifold_dim_ == 3) {
    // 3D Elements with possibly curved faces
    // We have to build a description of the element topology
    // and send it into the polyhedron volume and centroid
    // calculation routine
    int nf;
    std::vector<std::size_t> nfnodes;
    std::vector<AmanziGeometry::Point> cfcoords;

    auto [faces, fdirs] = getCellFacesAndDirections(cellid);

    nf = faces.size();
    nfnodes.resize(nf);

    for (int j = 0; j < nf; j++) {
      auto fcoords = getFaceCoordinates(faces[j]);
      nfnodes[j] = fcoords.size();

      if (fdirs[j] == 1) {
        for (int k = 0; k < nfnodes[j]; k++) cfcoords.push_back(fcoords[k]);
      } else {
        for (int k = nfnodes[j] - 1; k >= 0; k--) cfcoords.push_back(fcoords[k]);
      }
    }

    auto ccoords = getCellCoordinates(cellid);
    return AmanziGeometry::point_in_polyhed(p, ccoords, nf, nfnodes, cfcoords);

  } else if (manifold_dim_ == 2) {
    auto ccoords = getCellCoordinates(cellid);
    return AmanziGeometry::point_in_polygon(p, ccoords);
  }

  return false;
}


void
MeshCacheDevice::PrintMeshStatistics() const
{
  // auto vo_ = Teuchos::rcp(new VerboseObject("Mesh Output", *plist_));
  // if (vo_.get() && vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
  //   int ncells = getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  //   int nfaces = getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::OWNED);
  //   int nnodes = getNumEntities(AmanziMesh::NODE, AmanziMesh::Parallel_kind::OWNED);
  //   int nedges(0);
  //   if (has_edges_) nedges = getNumEntities(AmanziMesh::EDGE, AmanziMesh::Parallel_kind::OWNED);

  //   int min_out[4], max_out[4], sum_out[4], tmp_in[4] = { ncells, nfaces, nedges, nnodes };
  //   getComm()->MinAll(tmp_in, min_out, 4);
  //   getComm()->MaxAll(tmp_in, max_out, 4);
  //   getComm()->SumAll(tmp_in, sum_out, 4);

  //   Teuchos::OSTab tab = vo_->getOSTab();
  //   *vo_->os() << "cells, tot/min/max: " << sum_out[0] << "/" << min_out[0] << "/" << max_out[0]
  //              << "\n";
  //   *vo_->os() << "faces, tot/min/max: " << sum_out[1] << "/" << min_out[1] << "/" << max_out[1]
  //              << "\n";
  //   if (has_edges_)
  //     *vo_->os() << "edges, tot/min/max: " << sum_out[2] << "/" << min_out[2] << "/" << max_out[2]
  //                << "\n";
  //   *vo_->os() << "nodes, tot/min/max: " << sum_out[3] << "/" << min_out[3] << "/" << max_out[3]
  //              << "\n\n";
  // }
}

} // namespace Amanzi::AmanziMesh

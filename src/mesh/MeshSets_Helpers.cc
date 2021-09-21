/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
// Helper functions for resolving regions on meshes.


#include "RegionLabeledSet.hh"
#include "MeshCache.hh"
#include "MeshSets_Helpers.hh"

namespace Amanzi {
namespace AmanziMesh {

Entity_ID_List
resolveMeshSet(const AmanziGeometry::Region& region,
               const Entity_kind kind,
               const Parallel_type ptype,
               const MeshCache& mesh)
{
  if (AmanziGeometry::RegionType::ENUMERATED == region.get_type()) {
    auto region_enumerated = dynamic_cast<const AmanziGeometry::RegionEnumerated*>(&region);
    AMANZI_ASSERT(region_enumerated);
    return resolveMeshSetEnumerated(*region_enumerated, kind, ptype, mesh);

  } else if (AmanziGeometry::RegionType::LOGICAL == region.get_type()) {
    auto region_logical = dynamic_cast<const AmanziGeometry::RegionLogical*>(&region);
    AMANZI_ASSERT(region_logical);
    return resolveMeshSetLogical(*region_logical, kind, ptype, mesh);

  } else if (AmanziGeometry::RegionType::LABELEDSET == region.get_type()) {
    auto region_ls = dynamic_cast<const AmanziGeometry::RegionLabeledSet*>(&region);
    AMANZI_ASSERT(region_ls);
    return resolveMeshSetLabeledSet(*region_ls, kind, ptype, mesh);

  } else if (AmanziGeometry::RegionType::ALL == region.get_type()) {
    return resolveMeshSetAll(region, kind, ptype, mesh);

  } else if (AmanziGeometry::RegionType::BOUNDARY == region.get_type()) {
    return resolveMeshSetBoundary(region, kind, ptype, mesh);

  } else {
    // geometric
    return resolveMeshSetGeometric(region, kind, ptype, mesh);
  }
}


Entity_ID_List
resolveMeshSetAll(const AmanziGeometry::Region& region,
        const Entity_kind kind,
        const Parallel_type ptype,
        const MeshCache& mesh)
{
  auto num_ents = mesh.getNumEntities(kind, ptype);
  Entity_ID_List ents(num_ents);
  for (Entity_ID i=0; i!=num_ents; ++i) ents[i] = i;
  return ents;
}


Entity_ID_List
resolveMeshSetBoundary(const AmanziGeometry::Region& region,
        const Entity_kind kind,
        const Parallel_type ptype,
        const MeshCache& mesh)
{
  AMANZI_ASSERT(AmanziGeometry::RegionType::BOUNDARY == region.get_type());
  Entity_ID_List ents;
  if (kind == Entity_kind::FACE) {
    ents = mesh.getBoundaryFaces();
  } else if (kind == Entity_kind::NODE) {
    ents = mesh.getBoundaryNodes();
  } else {
    Errors::Message msg;
    msg << "Developer Error: MeshCache::getSetEntities() on region \"" << region.get_name()
        << "\" of type BOUNDARY was requested with invalid type " << to_string(kind);
    Exceptions::amanzi_throw(msg);
  }
  return ents;
}


Entity_ID_List
resolveMeshSetEnumerated(const AmanziGeometry::RegionEnumerated& region,
                         const Entity_kind kind,
                         const Parallel_type ptype,
                         const MeshCache& mesh)
{
  AMANZI_ASSERT(region.get_space_dimension() == mesh.getSpaceDimension());
  AMANZI_ASSERT(false);
  return Entity_ID_List();
}


Entity_ID_List
resolveMeshSetGeometric(const AmanziGeometry::Region& region,
        const Entity_kind kind,
        const Parallel_type ptype,
        const MeshCache& mesh)
{
  AMANZI_ASSERT(region.is_geometric());
  AMANZI_ASSERT(region.get_space_dimension() == mesh.getSpaceDimension());

  // find the extent
  Entity_ID begin, end;
  switch(ptype) {
    case (Parallel_type::GHOST) :
      begin = mesh.getNumEntities(kind, Parallel_type::OWNED);
      end = mesh.getNumEntities(kind, Parallel_type::ALL);
      break;
    default :
      begin = 0;
      end = mesh.getNumEntities(kind, ptype);
  }

  // collect centroids
  std::function<AmanziGeometry::Point(Entity_ID)> getCentroid;
  switch(kind) {
    case (Entity_kind::CELL) :
      getCentroid = [&mesh](Entity_ID i) { return mesh.getCellCentroid(i); };
      break;
    case (Entity_kind::FACE) :
      getCentroid = [&mesh](Entity_ID i) { return mesh.getFaceCentroid(i); };
      break;
    case (Entity_kind::BOUNDARY_FACE) :
      getCentroid = [&mesh](Entity_ID i) {
        return mesh.getFaceCentroid(mesh.getBoundaryFaces()[i]); };
      break;
    case (Entity_kind::EDGE) :
      getCentroid = [&mesh](Entity_ID i) { return mesh.getEdgeCentroid(i); };
      break;
    case (Entity_kind::NODE) :
      getCentroid = [&mesh](Entity_ID i) { return mesh.getNodeCoordinate(i); };
      break;
    case (Entity_kind::BOUNDARY_NODE) :
      getCentroid = [&mesh](Entity_ID i) {
        return mesh.getNodeCoordinate(mesh.getBoundaryNodes()[i]); };
      break;
    default : {}
  }

  // check whether centroid is inside region
  Entity_ID_List entities(end - begin, -1);
  int lcv = 0;
  for (int i=begin; i!=end; ++i) {
    if (region.inside(getCentroid(i))) {
      entities[lcv++] = i;
    }
  }
  entities.resize(lcv);
  return entities;
}


Entity_ID_List
resolveMeshSetLabeledSet(const AmanziGeometry::RegionLabeledSet& region,
                         const Entity_kind kind,
                         const Parallel_type ptype,
                         const MeshCache& mesh)
{
  if (!mesh.getMeshFramework().get()) {
    Errors::Message msg;
    msg << "Developer Error: MeshCache::getSetEntities() on region \"" << region.get_name()
        << "\" of type LABLEDSET was requested for the first time after the framework mesh was deleted.";
    Exceptions::amanzi_throw(msg);
  }
  Entity_ID_List ents;
  mesh.getMeshFramework()->getSetEntities(region, kind, ptype, ents);
  return ents;
}


Entity_ID_List
resolveMeshSetLogical(const AmanziGeometry::RegionLogical& region,
                      const Entity_kind kind,
                      const Parallel_type ptype,
                      const MeshCache& mesh)
{
  Entity_ID_List result;
  switch(region.get_operation()) {
    case (AmanziGeometry::BoolOpType::COMPLEMENT) :
      // Get the set of ALL entities of the right kind and type.
      //
      // wow this is a fun hack.  Since an ALL region does not need any aspect
      // of region, we simply pass this region despite the fact that it is NOT
      // an ALL region!
      result = resolveMeshSetAll(region, kind, ptype, mesh);

      // then iterate and subtract
      for (const auto& rname : region.get_component_regions()) {
        auto comp_ents = mesh.getSetEntities(rname, kind, ptype);

        Entity_ID_List lresult;
        lresult.reserve(result.size());
        std::set_difference(result.begin(), result.end(),
                            comp_ents.begin(), comp_ents.end(),
                            std::back_inserter(lresult));
        result = std::move(lresult);
      }
      break;

    case(AmanziGeometry::BoolOpType::UNION) :
      for (const auto& rname : region.get_component_regions()) {
        auto comp_ents = mesh.getSetEntities(rname, kind, ptype);

        Entity_ID_List lresult;
        lresult.reserve(result.size() + comp_ents.size());

        std::set_union(result.begin(), result.end(),
                       comp_ents.begin(), comp_ents.end(),
                       std::back_inserter(lresult));
        result = std::move(lresult);
      }
      break;

    case(AmanziGeometry::BoolOpType::INTERSECT) :
      for (const auto& rname : region.get_component_regions()) {
        auto comp_ents = mesh.getSetEntities(rname, kind, ptype);

        Entity_ID_List lresult;
        lresult.reserve(std::max(result.size(), comp_ents.size()));

        std::set_intersection(result.begin(), result.end(),
                comp_ents.begin(), comp_ents.end(),
                std::back_inserter(lresult));
        result = std::move(lresult);
      }
      break;

    case(AmanziGeometry::BoolOpType::SUBTRACT) :
      for (const auto& rname : region.get_component_regions()) {
        auto comp_ents = mesh.getSetEntities(rname, kind, ptype);

        Entity_ID_List lresult;
        lresult.reserve(result.size());

        std::set_difference(result.begin(), result.end(),
                            comp_ents.begin(), comp_ents.end(),
                            std::back_inserter(lresult));
        result = std::move(lresult);
      }
      break;

    default : {
      // note this should have errored already!
      Errors::Message msg("RegionLogical: operation type not set");
      Exceptions::amanzi_throw(msg);
    }
  }
  return result;
}



} // namespace AmanziMesh
} // namespace Amanzi

/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
// Helper functions for resolving regions on meshes.


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
    return resolveMeshSet(*region_enumerated, kind, ptype, mesh);

  } else if (AmanziGeometry::RegionType::LOGICAL == region.get_type()) {
    auto region_logical = dynamic_cast<const AmanziGeometry::RegionLogical*>(&region);
    AMANZI_ASSERT(region_logical);
    return resolveMeshSet(*region_logical, kind, ptype, mesh);

  } else {
    // geometric
    return resolveMeshSetGeometric(region, kind, ptype, mesh);
  }
}


//
// The default implementation assumes that region is a geometric region.
//
Entity_ID_List
resolveMeshSetGeometric(const AmanziGeometry::Region& region,
        const Entity_kind kind,
        const Parallel_type ptype,
        const MeshCache& mesh)
{
  AMANZI_ASSERT(region.is_geometric());
  AMANZI_ASSERT(region.get_space_dimension() == mesh.getSpaceDimension());
  AMANZI_ASSERT(region.get_manifold_dimension() == mesh.getManifoldDimension());

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
  const Point_View* centroids = nullptr;
  switch(kind) {
    case (Entity_kind::CELL) :
      centroids = &mesh.cell_centroids;
      break;
    case (Entity_kind::FACE) :
      centroids = &mesh.face_centroids;
      break;
    case (Entity_kind::BOUNDARY_FACE) :
      centroids = &mesh.face_centroids;
      break;
    case (Entity_kind::EDGE) :
      centroids = &mesh.edge_centroids;
      break;
    case (Entity_kind::NODE) :
      centroids = &mesh.node_coordinates;
      break;
    case (Entity_kind::BOUNDARY_NODE) :
      centroids = &mesh.node_coordinates;
      break;
    default : {}
  }

  // check whether centroid is inside region
  Entity_ID_List entities(end - begin, -1);
  if (kind != Entity_kind::BOUNDARY_NODE &&
      kind != Entity_kind::BOUNDARY_FACE) {
    int lcv = 0;
    for (int i=begin; i!=end; ++i) {
      if (region.inside((*centroids)[i])) {
        entities[lcv] = i;
        lcv++;
      }
    }
    entities.resize(lcv);
  } else if (kind == Entity_kind::BOUNDARY_FACE) {
    auto& face_ids = mesh.getBoundaryFaces();
    int lcv = 0;
    for (int i=begin; i!=end; ++i) {
      if (region.inside((*centroids)[face_ids[i]])) {
        entities[lcv] = i;
        lcv++;
      }
    }
    entities.resize(lcv);
  } else if (kind == Entity_kind::BOUNDARY_NODE) {
    auto& node_ids = mesh.getBoundaryNodes();
    int lcv = 0;
    for (int i=begin; i!=end; ++i) {
      if (region.inside((*centroids)[node_ids[i]])) {
        entities[lcv] = i;
        lcv++;
      }
    }
    entities.resize(lcv);
  }
  return entities;
}


//
// Specialization for RegionLogical
//
Entity_ID_List
resolveMeshSet(const AmanziGeometry::RegionLogical& region,
               const Entity_kind kind,
               const Parallel_type ptype,
               const MeshCache& mesh)
{
  AMANZI_ASSERT(region.get_space_dimension() == mesh.getSpaceDimension());
  AMANZI_ASSERT(region.get_manifold_dimension() == mesh.getManifoldDimension());
  AMANZI_ASSERT(false);
  return Entity_ID_List();
}


//
// Specialization for RegionEnumerated
//
Entity_ID_List
resolveMeshSet(const AmanziGeometry::RegionEnumerated& region,
               const Entity_kind kind,
               const Parallel_type ptype,
               const MeshCache& mesh)
{
  AMANZI_ASSERT(region.get_space_dimension() == mesh.getSpaceDimension());
  AMANZI_ASSERT(region.get_manifold_dimension() == mesh.getManifoldDimension());
  AMANZI_ASSERT(false);
  return Entity_ID_List();
}


} // namespace AmanziMesh
} // namespace Amanzi

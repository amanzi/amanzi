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
#include "Mesh_Algorithms.hh"

namespace Amanzi {
namespace AmanziMesh {

//
// This function simply delegates to the other functions in Impl, depending
// upon region type, mesh type (extracted meshes) and entity type.
//
View_type<const Entity_ID, MemSpace_kind::HOST>
resolveMeshSet(const AmanziGeometry::Region& region,
               const Entity_kind kind,
               const Parallel_kind ptype,
               const MeshCache<MemSpace_kind::HOST>& mesh)
{
  // special function to deal with BOUNDARY_*, so no other function needs to
  // deal with them.
  if (kind == Entity_kind::BOUNDARY_NODE || kind == Entity_kind::BOUNDARY_FACE) {
    return Impl::resolveBoundaryEntityMeshSet(region, kind, ptype, mesh);
  }

  // If there is no parent, this is not an extracted mesh
  auto parent_mesh = mesh.getParentMesh();
  if (parent_mesh == Teuchos::null) { return Impl::resolveMeshSet_(region, kind, ptype, mesh); }

  // there is a parent..
  auto region_type = region.get_type();
  if (region_type == AmanziGeometry::RegionType::ENUMERATED ||
      region_type == AmanziGeometry::RegionType::LABELEDSET) {
    // enumerated and labeled sets can only refer to the parent mesh, so we
    // check the parent
    return Impl::resolveIDMeshSetFromParent(region, kind, ptype, mesh, *parent_mesh);

  } else if (region_type == AmanziGeometry::RegionType::LOGICAL ||
             region_type == AmanziGeometry::RegionType::ALL ||
             region_type == AmanziGeometry::RegionType::BOUNDARY) {
    // these are resolved on this mesh
    return Impl::resolveMeshSet_(region, kind, ptype, mesh);

  } else {
    // must be a geometric region -- check the dimensions
    if (mesh.getSpaceDimension() == region.get_space_dimension()) {
      // if the dimension is right, we should be able to resolve this
      return Impl::resolveMeshSet_(region, kind, ptype, mesh);
    } else {
      // try resolving it on the parent instead
      return Impl::resolveGeometricMeshSetFromParent(region, kind, ptype, mesh, *parent_mesh);
    }
  }
}

View_type<const Entity_ID, MemSpace_kind::HOST>
resolveMeshSetVolumeFractions(const AmanziGeometry::Region& region,
                              const Entity_kind kind,
                              const Parallel_kind ptype,
                              Double_View& vol_fracs,
                              const MeshCache<MemSpace_kind::HOST>& mesh)
{
  Kokkos::resize(vol_fracs, 0);
  Entity_ID_View ents;
  auto ents_cpt = 0;
  auto vol_fracs_cpt = 0;

  if ((AmanziGeometry::RegionType::BOX_VOF == region.get_type() ||
       AmanziGeometry::RegionType::LINE_SEGMENT == region.get_type()) &&
      (kind == Entity_kind::CELL || kind == Entity_kind::FACE)) {
    if (kind == Entity_kind::CELL) {
      auto ncells = mesh.getNumEntities(Entity_kind::CELL, ptype);
      Kokkos::resize(vol_fracs, ncells);
      Kokkos::resize(ents, ncells);

      for (int c = 0; c != ncells; ++c) {
        auto polytope_nodes = mesh.getCellCoordinates(c);
        std::vector<Entity_ID_List> polytope_faces;

        if (mesh.getSpaceDimension() == 3) {
          auto cnodes = mesh.getCellNodes(c);
          auto cfd = mesh.getCellFacesAndDirections(c);
          auto cfaces = cfd.first;
          auto cfdirs = cfd.second;
          polytope_faces.resize(cfaces.size());

          for (int n = 0; n < cfaces.size(); ++n) {
            auto fnodes = mesh.getFaceNodes(cfaces[n]);

            for (int i = 0; i != fnodes.size(); ++i) {
              int j = (cfdirs[n] > 0) ? i : fnodes.size() - i - 1;
              int pos =
                std::distance(cnodes.begin(), std::find(cnodes.begin(), cnodes.end(), fnodes[j]));
              polytope_faces[n].push_back(pos);
            }
          }
        }

        auto vpolytope_nodes = asVector(polytope_nodes);
        double volume = region.intersect(vpolytope_nodes, polytope_faces);
        if (volume > 0.0) {
          ents[ents_cpt++] = c;
          if (region.get_type() == AmanziGeometry::RegionType::LINE_SEGMENT)
            vol_fracs[vol_fracs_cpt++] = volume;
          else
            vol_fracs[vol_fracs_cpt++] = volume / mesh.getCellVolume(c);
        }
      }
      Kokkos::resize(ents, ents_cpt);
      Kokkos::resize(vol_fracs, vol_fracs_cpt);
    } else {
      // ind == FACE
      int nfaces = mesh.getNumEntities(Entity_kind::FACE, ptype);
      Kokkos::resize(vol_fracs, nfaces);
      Kokkos::resize(ents, nfaces);

      Point_List polygon;

      for (int f = 0; f != nfaces; ++f) {
        auto polygon = mesh.getFaceCoordinates(f);
        auto vpolygon = asVector(polygon);
        double area = region.intersect(vpolygon);
        if (area > 0.0) {
          ents[ents_cpt++] = f;
          vol_fracs[vol_fracs_cpt++] = area / mesh.getFaceArea(f);
        }
      }
      Kokkos::resize(ents, ents_cpt);
      Kokkos::resize(vol_fracs, vol_fracs_cpt);
    }

  } else {
    ents = mesh.getSetEntities(region.get_name(), kind, ptype);
  }
  return ents;
}


namespace Impl {
//
// This helper function resolves sets on this mesh, never a parent mesh.
//
View_type<const Entity_ID, MemSpace_kind::HOST>
resolveMeshSet_(const AmanziGeometry::Region& region,
                const Entity_kind kind,
                const Parallel_kind ptype,
                const MeshCache<MemSpace_kind::HOST>& mesh)
{
  View_type<const Entity_ID, MemSpace_kind::HOST> result;
  if (AmanziGeometry::RegionType::ENUMERATED == region.get_type()) {
    auto region_enumerated = dynamic_cast<const AmanziGeometry::RegionEnumerated*>(&region);
    AMANZI_ASSERT(region_enumerated);
    result = resolveMeshSetEnumerated(*region_enumerated, kind, ptype, mesh);

  } else if (AmanziGeometry::RegionType::LOGICAL == region.get_type()) {
    auto region_logical = dynamic_cast<const AmanziGeometry::RegionLogical*>(&region);
    AMANZI_ASSERT(region_logical);
    result = resolveMeshSetLogical(*region_logical, kind, ptype, mesh);

  } else if (AmanziGeometry::RegionType::LABELEDSET == region.get_type()) {
    auto region_ls = dynamic_cast<const AmanziGeometry::RegionLabeledSet*>(&region);
    AMANZI_ASSERT(region_ls);
    Entity_ID_View lresult;
    result = resolveMeshSetLabeledSet(*region_ls, kind, ptype, mesh);
    // labeled sets may not be sorted, though all other types are.  Sort labeled sets.
    lresult.fromConst(result);
    std::sort(lresult.begin(), lresult.end());
    result = lresult;

  } else if (AmanziGeometry::RegionType::ALL == region.get_type()) {
    result = resolveMeshSetAll(region, kind, ptype, mesh);
  } else if (AmanziGeometry::RegionType::BOUNDARY == region.get_type()) {
    result = resolveMeshSetBoundary(region, kind, ptype, mesh);
  } else if (AmanziGeometry::RegionType::LINE_SEGMENT == region.get_type()) {
    Double_View vofs;
    result = resolveMeshSetVolumeFractions(region, kind, ptype, vofs, mesh);
  } else if (AmanziGeometry::RegionType::POINT == region.get_type() && kind == Entity_kind::CELL) {
    auto region_point = dynamic_cast<const AmanziGeometry::RegionPoint*>(&region);
    result = resolveMeshSetPoint(*region_point, kind, ptype, mesh);
  } else {
    // geometric
    result = resolveMeshSetGeometric(region, kind, ptype, mesh);
  }
  return result;
}

cEntity_ID_View
resolveMeshSetPoint(const AmanziGeometry::RegionPoint& region,
                    const Entity_kind kind,
                    const Parallel_kind ptype,
                    const MeshCache<MemSpace_kind::HOST>& mesh)
{
  std::vector<Entity_ID> result;

  auto rgnpnt = region.point();

  int nnode = mesh.getNumEntities(NODE, Parallel_kind::ALL);
  double mindist2 = 1.e+16;
  int minnode = -1;
  int inode;
  for (inode = 0; inode < nnode; inode++) {
    auto vpnt = mesh.getNodeCoordinate(inode);
    double dist2 = (vpnt - rgnpnt) * (vpnt - rgnpnt);

    if (dist2 < mindist2) {
      mindist2 = dist2;
      minnode = inode;
      if (mindist2 <= 1.0e-14) break;
    }
  }

  auto cells = mesh.getNodeCells(minnode, ptype);

  int ncells = cells.size();
  for (int ic = 0; ic < ncells; ic++) {
    Entity_ID icell = cells[ic];

    // Check if point is contained in cell
    if (mesh.isPointInCell(rgnpnt, icell)) result.push_back(icell);
  }

  // finally check all cells, typical for anisotropic meshes
  if (result.empty()) {
    int ncells = mesh.getNumEntities(CELL, ptype);
    for (int c = 0; c < ncells; ++c)
      if (mesh.isPointInCell(rgnpnt, c)) result.push_back(c);
  }

  Entity_ID_View rv;
  vectorToView(rv, result);
  return rv;
}


// This helper function deals with BOUNDARY_* entity requests.
//
// This should probably be implemented, but it isn't currently supported by
// current master or by tests so left undone for now.
//
// Conceptually the implementation should just call with NODE instead of
// BOUNDARY_NODE, then filter to only BOUNDARY_NODE entities (respectively
// FACE).
View_type<const Entity_ID, MemSpace_kind::HOST>
resolveBoundaryEntityMeshSet(const AmanziGeometry::Region& region,
                             const Entity_kind kind,
                             const Parallel_kind ptype,
                             const MeshCache<MemSpace_kind::HOST>& parent_mesh)
{
  Errors::Message msg;
  msg << "Resolving " << to_string(kind) << " not yet implemented...";
  Exceptions::amanzi_throw(msg);
  return Entity_ID_View();
}


//
// Resolves sets that are discretely enumerated on the parent mesh.
//
// Remember, this could be volume extraction or surface extraction.
View_type<const Entity_ID, MemSpace_kind::HOST>
resolveIDMeshSetFromParent(const AmanziGeometry::Region& region,
                           const Entity_kind kind,
                           const Parallel_kind ptype,
                           const MeshCache<MemSpace_kind::HOST>& mesh,
                           const MeshCache<MemSpace_kind::HOST>& parent_mesh)
{
  // Enumerated or Labeled sets that have IDs, and so also have Entity_kind as
  // a region property.
  Entity_kind parent_kind;
  if (AmanziGeometry::RegionType::ENUMERATED == region.get_type()) {
    auto region_enumerated = dynamic_cast<const AmanziGeometry::RegionEnumerated*>(&region);
    AMANZI_ASSERT(region_enumerated);
    parent_kind = createEntityKind(region_enumerated->entity_str());
  } else if (AmanziGeometry::RegionType::LABELEDSET == region.get_type()) {
    auto region_ls = dynamic_cast<const AmanziGeometry::RegionLabeledSet*>(&region);
    AMANZI_ASSERT(region_ls);
    parent_kind = createEntityKind(region_ls->entity_str());
  } else {
    parent_kind = CELL;
    AMANZI_ASSERT(false); // this should not be reachable
  }

  // Get parent mesh entities
  View_type<const Entity_ID, MemSpace_kind::HOST> parent_entities =
    resolveMeshSet(region, parent_kind, ptype, parent_mesh);

  // filter, leaving only the extracted entities
  if (kind == Entity_kind::NODE) {
    if (parent_kind != Entity_kind::NODE) {
      Errors::Message msg;
      msg << "Mesh: cannot resolve entities of type " << to_string(kind)
          << " from a parent mesh region \"" << region.get_name()
          << "\" which supplies entities of type " << to_string(parent_kind);
      Exceptions::amanzi_throw(msg);
    }
    return Impl::filterParentEntities(mesh, Entity_kind::NODE, ptype, parent_entities);

  } else if (kind == Entity_kind::CELL) {
    if (parent_kind == Entity_kind::FACE) {
      // check -- this must be a surface extraction
      if (mesh.getManifoldDimension() != (parent_mesh.getManifoldDimension() - 1)) {
        Errors::Message msg;
        msg << "Mesh: cannot resolve entities of type CELL"
            << " from a parent mesh region FACE when the"
            << " extraction was not to a lower-dimension mesh.";
        Exceptions::amanzi_throw(msg);
      }
      return Impl::filterParentEntities(mesh, Entity_kind::CELL, ptype, parent_entities);

    } else if (parent_kind == Entity_kind::CELL) {
      if (mesh.getManifoldDimension() == parent_mesh.getManifoldDimension()) {
        // filter, direct parent from volumetric mesh
        return Impl::filterParentEntities(mesh, Entity_kind::CELL, ptype, parent_entities);
      } else {
        // filter, but look at the implied parent face's internal cell
        return Impl::filterParentEntities_SurfaceCellToCell(mesh, ptype, parent_entities);
      }
    } else {
      Errors::Message msg;
      msg << "Mesh: cannot resolve entities of type " << to_string(kind)
          << " from a parent mesh region \"" << region.get_name()
          << "\" which supplies entities of type " << to_string(parent_kind);
      Exceptions::amanzi_throw(msg);
    }

  } else if (kind == Entity_kind::FACE) {
    if (parent_kind == Entity_kind::FACE) {
      if (mesh.getManifoldDimension() == parent_mesh.getManifoldDimension()) {
        // filter, direct parent from volumetric mesh
        return Impl::filterParentEntities(mesh, Entity_kind::FACE, ptype, parent_entities);
      } else {
        // filter, how do we get this one?
        return Impl::filterParentEntities_SurfaceFaceToFace(mesh, ptype, parent_entities);
      }
    } else if (parent_kind == Entity_kind::EDGE) {
      // filter, direct parent?  Unless something on FACE-to-FACE is weird for
      // surface extraction?
      return Impl::filterParentEntities(mesh, Entity_kind::FACE, ptype, parent_entities);

    } else {
      Errors::Message msg;
      msg << "Mesh: cannot resolve entities of type " << to_string(kind)
          << " from a parent mesh region \"" << region.get_name()
          << "\" which supplies entities of type " << to_string(parent_kind);
      Exceptions::amanzi_throw(msg);
    }
  }

  Errors::Message msg;
  msg << "Mesh: cannot resolve entities of type " << to_string(kind)
      << " from a parent mesh region \"" << region.get_name()
      << "\" which supplies entities of type " << to_string(parent_kind);
  Exceptions::amanzi_throw(msg);
  return Entity_ID_View();
}


//
// Looking for entities on a region that is a different space_dimension than
// the mesh.  Try to resolve the region on the parent mesh, then collect the
// parent entity.  This is a lot like the ID version above, but we don't have
// the region's entity kind to help us guess...
//
View_type<const Entity_ID, MemSpace_kind::HOST>
resolveGeometricMeshSetFromParent(const AmanziGeometry::Region& region,
                                  const Entity_kind kind,
                                  const Parallel_kind ptype,
                                  const MeshCache<MemSpace_kind::HOST>& mesh,
                                  const MeshCache<MemSpace_kind::HOST>& parent_mesh)
{
  if (mesh.getManifoldDimension() == (parent_mesh.getManifoldDimension() - 1)) {
    // extracted surface mesh
    if (kind == Entity_kind::NODE) {
      // nodes are nodes -- resolve and filter
      View_type<const Entity_ID, MemSpace_kind::HOST> parent_entities =
        resolveMeshSet(region, Entity_kind::NODE, ptype, parent_mesh);
      return Impl::filterParentEntities(mesh, Entity_kind::NODE, ptype, parent_entities);
    } else if (kind == Entity_kind::CELL) {
      // cells are faces -- resolve and filter
      View_type<const Entity_ID, MemSpace_kind::HOST> parent_entities =
        resolveMeshSet(region, Entity_kind::FACE, ptype, parent_mesh);
      return Impl::filterParentEntities(mesh, Entity_kind::CELL, ptype, parent_entities);
    } else if (kind == Entity_kind::FACE) {
      // faces are edges, which are typically not formal entities -- how do we deal with this?
      Errors::Message msg;
      msg << "Mesh: cannot resolve set entities on parent mesh "
          << "for entities of kind \"" << to_string(kind) << "\"";
      Exceptions::amanzi_throw(msg);
    } else {
      Errors::Message msg;
      msg << "Mesh: cannot resolve set entities on parent mesh "
          << "for entities of kind \"" << to_string(kind) << "\"";
      Exceptions::amanzi_throw(msg);
    }

  } else if (mesh.getManifoldDimension() == parent_mesh.getManifoldDimension()) {
    // extracted volume mesh -- fail, this should already be dealt with or
    // something else is really weird.
    Errors::Message msg;
    msg << "Mesh: cannot resolve set entities on parent mesh.  "
        << "Likely this is developer error.";
    Exceptions::amanzi_throw(msg);

  } else {
    // who knows
    Errors::Message msg;
    msg << "Mesh: cannot resolve set entities on parent mesh, "
        << "meshes are questionably related for mesh of manifold dimension "
        << mesh.getManifoldDimension() << " and parent of manifold dimension "
        << parent_mesh.getManifoldDimension();
    Exceptions::amanzi_throw(msg);
  }
  return Entity_ID_View();
}

//
// Helper function to return ALL entities
//
View_type<const Entity_ID, MemSpace_kind::HOST>
resolveMeshSetAll(const AmanziGeometry::Region& region,
                  const Entity_kind kind,
                  const Parallel_kind ptype,
                  const MeshCache<MemSpace_kind::HOST>& mesh)
{
  auto num_ents = mesh.getNumEntities(kind, ptype);
  Entity_ID_View ents("ents", num_ents);
  for (Entity_ID i = 0; i != num_ents; ++i) ents[i] = i;
  return ents;
}


//
// Helper function to return BOUNDARY entities
//
View_type<const Entity_ID, MemSpace_kind::HOST>
resolveMeshSetBoundary(const AmanziGeometry::Region& region,
                       const Entity_kind kind,
                       const Parallel_kind ptype,
                       const MeshCache<MemSpace_kind::HOST>& mesh)
{
  AMANZI_ASSERT(AmanziGeometry::RegionType::BOUNDARY == region.get_type());
  Entity_ID_View ents;
  Entity_kind boundary_kind;
  if (kind == Entity_kind::FACE) {
    ents = mesh.getBoundaryFaces();
    boundary_kind = Entity_kind::BOUNDARY_FACE;
  } else if (kind == Entity_kind::NODE) {
    ents = mesh.getBoundaryNodes();
    boundary_kind = Entity_kind::BOUNDARY_NODE;
  } else {
    Errors::Message msg;
    msg << "Developer Error: MeshCache::getSetEntities() on region \"" << region.get_name()
        << "\" of type BOUNDARY was requested with invalid type " << to_string(kind);
    Exceptions::amanzi_throw(msg);
  }

  // the above calls always return ALL entities, adjust if needed
  if (ptype == Parallel_kind::OWNED) {
    // keep only the first OWNED ones
    Kokkos::resize(ents, mesh.getNumEntities(boundary_kind, ptype));
  } else if (ptype == Parallel_kind::GHOST) {
    // keep owned -- end
    int ct = 0;
    int numbk = mesh.getNumEntities(boundary_kind, Parallel_kind::OWNED);
    for (int i = numbk; i < ents.size(); ++i) ents[ct++] = ents[i];
    Kokkos::resize(ents, ct);
  }

  return ents;
}

View_type<const Entity_ID, MemSpace_kind::HOST>
resolveMeshSetEnumerated(const AmanziGeometry::RegionEnumerated& region,
                         const Entity_kind kind,
                         const Parallel_kind ptype,
                         const MeshCache<MemSpace_kind::HOST>& mesh)
{
  if (kind == createEntityKind(region.entity_str())) {
    Entity_ID_List region_entities = region.entities();
    bool ghosted = (ptype != Parallel_kind::OWNED);
    auto& mesh_map = mesh.getMap(kind, ghosted);

    std::size_t cpt = 0;
    Entity_ID_View result("SetEnumerated", mesh_map.NumMyElements());
    for (Entity_GID gid : region_entities) {
      Entity_ID lid = mesh_map.LID(gid);
      if (lid >= 0) result[cpt++] = lid;
    }
    Kokkos::resize(result, cpt);
    return result;
  } else {
    return Entity_ID_View();
  }
}

View_type<const Entity_ID, MemSpace_kind::HOST>
resolveMeshSetGeometric(const AmanziGeometry::Region& region,
                        const Entity_kind kind,
                        const Parallel_kind ptype,
                        const MeshCache<MemSpace_kind::HOST>& mesh)
{
  AMANZI_ASSERT(region.is_geometric());
  AMANZI_ASSERT(region.get_space_dimension() == mesh.getSpaceDimension());

  // find the extent
  Entity_ID begin, end;
  if (ptype == Parallel_kind::GHOST) {
    begin = mesh.getNumEntities(kind, Parallel_kind::OWNED);
    end = mesh.getNumEntities(kind, Parallel_kind::ALL);
  } else {
    begin = 0;
    end = mesh.getNumEntities(kind, ptype);
  }

  // check whether centroid is inside region
  Entity_ID_View entities("entities", end - begin);
  int lcv = 0;
  for (Entity_ID i = begin; i != end; ++i) {
    if (region.inside(mesh.getCentroid(kind, i))) { entities[lcv++] = i; }
  }
  Kokkos::resize(entities, lcv);
  return entities;
}


View_type<const Entity_ID, MemSpace_kind::HOST>
resolveMeshSetLabeledSet(const AmanziGeometry::RegionLabeledSet& region,
                         const Entity_kind kind,
                         const Parallel_kind ptype,
                         const MeshCache<MemSpace_kind::HOST>& mesh)
{
  if (!mesh.getMeshFramework().get()) {
    Errors::Message msg;
    msg << "Developer Error: MeshCache::getSetEntities() on region \"" << region.get_name()
        << "\" of type LABLEDSET was requested for the first time after the framework mesh was "
           "deleted.";
    Exceptions::amanzi_throw(msg);
  }
  View_type<const Entity_ID, MemSpace_kind::HOST> ents;
  mesh.getMeshFramework()->getSetEntities(region, kind, ptype, ents);
  return ents;
}


View_type<const Entity_ID, MemSpace_kind::HOST>
resolveMeshSetLogical(const AmanziGeometry::RegionLogical& region,
                      const Entity_kind kind,
                      const Parallel_kind ptype,
                      const MeshCache<MemSpace_kind::HOST>& mesh)
{
  View_type<const Entity_ID, MemSpace_kind::HOST> result;
  switch (region.get_operation()) {
  case (AmanziGeometry::BoolOpType::COMPLEMENT): {
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
      std::set_difference(result.begin(),
                          result.end(),
                          comp_ents.begin(),
                          comp_ents.end(),
                          std::back_inserter(lresult));
      vectorToConstView(result, lresult);
    }

  } break;
  case (AmanziGeometry::BoolOpType::UNION): {
    for (const auto& rname : region.get_component_regions()) {
      auto comp_ents = mesh.getSetEntities(rname, kind, ptype);

      Entity_ID_List lresult;
      lresult.reserve(result.size() + comp_ents.size());

      std::set_union(result.begin(),
                     result.end(),
                     comp_ents.begin(),
                     comp_ents.end(),
                     std::back_inserter(lresult));
      vectorToConstView(result, lresult);
    }

  } break;
  case (AmanziGeometry::BoolOpType::INTERSECT): {
    const auto& rnames = region.get_component_regions();
    AMANZI_ASSERT(rnames.size() > 1);
    result = mesh.getSetEntities(rnames[0], kind, ptype);
    for (int i = 1; i != rnames.size(); ++i) {
      auto comp_ents = mesh.getSetEntities(rnames[i], kind, ptype);

      Entity_ID_List lresult;
      lresult.reserve(std::max(result.size(), comp_ents.size()));

      std::set_intersection(result.begin(),
                            result.end(),
                            comp_ents.begin(),
                            comp_ents.end(),
                            std::back_inserter(lresult));
      vectorToConstView(result, lresult);
    }

  } break;
  case (AmanziGeometry::BoolOpType::SUBTRACT): {
    const auto& rnames = region.get_component_regions();
    AMANZI_ASSERT(rnames.size() > 1);
    result = mesh.getSetEntities(rnames[0], kind, ptype);
    for (int i = 1; i != rnames.size(); ++i) {
      auto comp_ents = mesh.getSetEntities(rnames[i], kind, ptype);

      Entity_ID_List lresult;
      lresult.reserve(result.size());

      std::set_difference(result.begin(),
                          result.end(),
                          comp_ents.begin(),
                          comp_ents.end(),
                          std::back_inserter(lresult));
      vectorToConstView(result, lresult);
    }
  } break;
  default: {
    // note this should have errored already!
    Errors::Message msg("RegionLogical: operation type not set");
    Exceptions::amanzi_throw(msg);
  }
  }
  return result;
}

//
// Filter for the direct parent entity.
View_type<const Entity_ID, MemSpace_kind::HOST>
filterParentEntities(const MeshCache<MemSpace_kind::HOST>& mesh,
                     Entity_kind kind,
                     Parallel_kind ptype,
                     const View_type<const Entity_ID, MemSpace_kind::HOST>& parent_entities)
{
  Entity_ID_View result;
  Entity_ID_List vresult;

  // This is possibly a poorly performing algorithm.  We do linear search for
  // each entry here, but parent_entities is gauranteed to be sorted, and I
  // believe that the list of all entities's parents are also sorted, so it
  // ought to be possible to do something fancy with iterators to make this a
  // much faster operation.  That said, we do this infrequently (likely only
  // once per set) and so it is not worth the complexity.
  int n_entities = mesh.getNumEntities(kind, ptype);
  for (int i = 0; i != n_entities; ++i) {
    Entity_ID parent_id = mesh.getEntityParent(kind, i);
    for (int j = 0; j != parent_entities.size(); ++j) {
      if (parent_entities[j] == parent_id)
        vresult.emplace_back(i);
      else if (parent_entities[j] > parent_id)
        break;
    }
  }
  vectorToView(result, vresult);
  return result;
}


View_type<const Entity_ID, MemSpace_kind::HOST>
filterParentEntities_SurfaceCellToCell(
  const MeshCache<MemSpace_kind::HOST>& mesh,
  Parallel_kind ptype,
  const View_type<const Entity_ID, MemSpace_kind::HOST>& parent_entities)
{
  Entity_ID_View result;
  Entity_ID_List vresult;

  // This is possibly a poorly performing algorithm.  We do linear search for
  // each entry here, but parent_entities is gauranteed to be sorted, and I
  // believe that the list of all entities's parents are also sorted, so it
  // ought to be possible to do something fancy with iterators to make this a
  // much faster operation.  That said, we do this infrequently (likely only
  // once per set) and so it is not worth the complexity.
  int n_entities = mesh.getNumEntities(Entity_kind::CELL, ptype);
  for (int i = 0; i != n_entities; ++i) {
    // cell's parent is a face
    Entity_ID face_parent_id = mesh.getEntityParent(Entity_kind::CELL, i);

    // better have one and only one internal cell
    Entity_ID parent_id = getFaceOnBoundaryInternalCell(*mesh.getParentMesh(), face_parent_id);

    for (int j = 0; j != parent_entities.size(); ++j) {
      if (parent_entities[j] == parent_id)
        vresult.emplace_back(i);
      else if (parent_entities[j] > parent_id)
        break;
    }
  }
  vectorToView(result, vresult);
  return result;
}


View_type<const Entity_ID, MemSpace_kind::HOST>
filterParentEntities_SurfaceFaceToFace(
  const MeshCache<MemSpace_kind::HOST>& mesh,
  Parallel_kind ptype,
  const View_type<const Entity_ID, MemSpace_kind::HOST>& parent_entities)
{
  Entity_ID_View result;
  Entity_ID_List vresult;

  // We wish to capture the relationship from surface faces (whose parent
  // entities are edges in the subsurface mesh) to subsurface faces whose
  // normal is horizontal (e.g. a column's lateral sides).  So we would like to
  // check that the parent entity (an edge) has a face that is in the set.
  // Note that this may not be possible if the parent mesh doesn't have edges,
  // so we may have to check nodes.
  if (mesh.getParentMesh()->hasEdges()) {
    // If we have edges, this is pretty straightforward
    auto num_faces = mesh.getNumEntities(Entity_kind::FACE, ptype);
    for (Entity_ID sf = 0; sf != num_faces; ++sf) {
      Entity_ID e = mesh.getEntityParent(Entity_kind::FACE, sf);

      for (Entity_ID f : parent_entities) {
        bool found = false;
        const auto fedges = mesh.getParentMesh()->getFaceEdges(f);
        for (int i = 0; i != fedges.size(); ++i) {
          if (e == fedges[i]) {
            found = true;
            break;
          }
        }
        if (found) {
          vresult.emplace_back(sf);
          break;
        }
      }
    }

  } else {
    // If we don't have edges, we have to compare nodes.
    auto num_faces = mesh.getNumEntities(Entity_kind::FACE, ptype);
    for (Entity_ID sf = 0; sf != num_faces; ++sf) {
      // for each surface face's nodes, get the parent node in the subsurface
      auto sfnodes = mesh.getFaceNodes(sf);
      AMANZI_ASSERT(sfnodes.size() == 2);
      Entity_ID_View parent_nodes("parent_nodes", sfnodes.size());
      for (int i = 0; i != sfnodes.size(); ++i)
        parent_nodes[i] = mesh.getEntityParent(Entity_kind::NODE, sfnodes[i]);

      // for each face in the set, are all of these nodes in the set of face
      // nodes?
      for (Entity_ID f : parent_entities) {
        const auto fnodes = mesh.getParentMesh()->getFaceNodes(f);
        bool found0 = false;
        for (int i = 0; i != fnodes.size(); ++i) {
          if (parent_nodes[0] == fnodes[i]) {
            found0 = true;
            break;
          }
        }
        if (!found0) continue;

        bool found1 = false;
        for (int i = 0; i != fnodes.size(); ++i) {
          if (parent_nodes[1] == fnodes[i]) {
            found1 = true;
            break;
          }
        }
        if (found0 && found1) {
          vresult.emplace_back(sf);
          break;
        }
      }
    }
  }
  vectorToView(result, vresult);
  return result;
}

} // namespace Impl
} // namespace AmanziMesh
} // namespace Amanzi

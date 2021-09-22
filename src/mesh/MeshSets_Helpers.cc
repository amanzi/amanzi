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
  std::cout << "Resolving set: " << region.get_name() << std::endl;
  Entity_ID_List result;
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
    result = resolveMeshSetLabeledSet(*region_ls, kind, ptype, mesh);
    // labeled sets may not be sorted, though all other types are.  Sort labeled sets.
    std::sort(result.begin(), result.end());

  } else if (AmanziGeometry::RegionType::ALL == region.get_type()) {
    result = resolveMeshSetAll(region, kind, ptype, mesh);

  } else if (AmanziGeometry::RegionType::BOUNDARY == region.get_type()) {
    result = resolveMeshSetBoundary(region, kind, ptype, mesh);

  } else {
    // geometric
    result = resolveMeshSetGeometric(region, kind, ptype, mesh);
  }

  int g_count = 0;
  int l_count = result.size();
  mesh.getComm()->SumAll(&l_count, &g_count, 1);
  if (g_count == 0) {
    // warn?  error?
    Errors::Message msg;
    msg << "AmanziMesh::resolveMeshSet: Region \"" << region.get_name() << "\" of type \"" << to_string(region.get_type()) << "\" is empty.";
    Exceptions::amanzi_throw(msg);
  }
  return result;
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
  if (ptype == Parallel_type::OWNED) {
    // keep only the first OWNED ones
    ents.resize(mesh.getNumEntities(boundary_kind, ptype));
  } else if (ptype == Parallel_type::GHOST) {
    // keep owned -- end
    Entity_ID_List ents2(ents.begin()+mesh.getNumEntities(boundary_kind, Parallel_type::OWNED),
                         ents.end());
    std::swap(ents, ents2);
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
  if (ptype == Parallel_type::GHOST) {
    begin = mesh.getNumEntities(kind, Parallel_type::OWNED);
    end = mesh.getNumEntities(kind, Parallel_type::ALL);
  } else {
    begin = 0;
    end = mesh.getNumEntities(kind, ptype);
  }

  // collect centroids
  std::function<AmanziGeometry::Point(Entity_ID)> getCentroid;
  switch(kind) {
    case (Entity_kind::CELL) : {
      getCentroid = [&mesh](Entity_ID i) { return mesh.getCellCentroid(i); };
    } break;
    case (Entity_kind::FACE) : {
      getCentroid = [&mesh](Entity_ID i) { return mesh.getFaceCentroid(i); };
    } break;
    case (Entity_kind::BOUNDARY_FACE) : {
      getCentroid = [&mesh](Entity_ID i) {
        return mesh.getFaceCentroid(mesh.getBoundaryFaces()[i]); };
    } break;
    case (Entity_kind::EDGE) : {
      getCentroid = [&mesh](Entity_ID i) { return mesh.getEdgeCentroid(i); };
    } break;
    case (Entity_kind::NODE) : {
      getCentroid = [&mesh](Entity_ID i) { return mesh.getNodeCoordinate(i); };
    } break;
    case (Entity_kind::BOUNDARY_NODE) : {
      getCentroid = [&mesh](Entity_ID i) {
        return mesh.getNodeCoordinate(mesh.getBoundaryNodes()[i]); };
    } break;
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
    case (AmanziGeometry::BoolOpType::COMPLEMENT) : {
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

    } break;
    case(AmanziGeometry::BoolOpType::UNION) : {
      for (const auto& rname : region.get_component_regions()) {
        auto comp_ents = mesh.getSetEntities(rname, kind, ptype);

        Entity_ID_List lresult;
        lresult.reserve(result.size() + comp_ents.size());

        std::set_union(result.begin(), result.end(),
                       comp_ents.begin(), comp_ents.end(),
                       std::back_inserter(lresult));
        result = std::move(lresult);
      }

    } break;
    case(AmanziGeometry::BoolOpType::INTERSECT) : {
      const auto& rnames = region.get_component_regions();
      AMANZI_ASSERT(rnames.size() > 1);
      result = mesh.getSetEntities(rnames[0], kind, ptype);
      for (int i=1; i!=rnames.size(); ++i) {
        auto comp_ents = mesh.getSetEntities(rnames[i], kind, ptype);

        Entity_ID_List lresult;
        lresult.reserve(std::max(result.size(), comp_ents.size()));

        std::set_intersection(result.begin(), result.end(),
                comp_ents.begin(), comp_ents.end(),
                std::back_inserter(lresult));
        result = std::move(lresult);
      }

    } break;
    case(AmanziGeometry::BoolOpType::SUBTRACT) : {
      const auto& rnames = region.get_component_regions();
      AMANZI_ASSERT(rnames.size() > 1);
      result = mesh.getSetEntities(rnames[0], kind, ptype);
      for (int i=1; i!=rnames.size(); ++i) {
        auto comp_ents = mesh.getSetEntities(rnames[i], kind, ptype);

        Entity_ID_List lresult;
        lresult.reserve(result.size());

        std::set_difference(result.begin(), result.end(),
                            comp_ents.begin(), comp_ents.end(),
                            std::back_inserter(lresult));
        result = std::move(lresult);
      }
    } break;
    default : {
      // note this should have errored already!
      Errors::Message msg("RegionLogical: operation type not set");
      Exceptions::amanzi_throw(msg);
    }
  }
  return result;
}


Entity_ID_List
resolveMeshSetVolumeFractions(const AmanziGeometry::Region& region,
        const Entity_kind kind,
        const Parallel_type ptype,
        Double_View& vol_fracs,
        const MeshCache& mesh)
{
  vol_fracs.resize(0);
  Entity_ID_List ents;

  if ((AmanziGeometry::RegionType::BOX_VOF == region.get_type() ||
       AmanziGeometry::RegionType::LINE_SEGMENT == region.get_type()) &&
      (kind == Entity_kind::CELL || kind == Entity_kind::FACE)) {

    if (kind == Entity_kind::CELL) {
      auto ncells = mesh.getNumEntities(Entity_kind::CELL, ptype);
      vol_fracs.reserve(ncells);
      ents.reserve(ncells);

      for (int c=0; c!=ncells; ++c) {
        auto polytope_nodes = mesh.getCellCoordinates(c);
        std::vector<Entity_ID_List> polytope_faces;

        if (mesh.getSpaceDimension() == 3) {
          auto cnodes = mesh.getCellNodes(c);
          Entity_ID_View cfaces;
          Entity_Direction_View cfdirs;
          std::tie(cfaces, cfdirs) = mesh.getCellFacesAndDirections(c);
          polytope_faces.resize(cfaces.size());

          for (int n = 0; n < cfaces.size(); ++n) {
            auto fnodes = mesh.getFaceNodes(cfaces[n]);

            for (int i=0; i!=fnodes.size(); ++i) {
              int j = (cfdirs[n] > 0) ? i : fnodes.size() - i - 1;
              int pos = std::distance(cnodes.begin(), std::find(cnodes.begin(), cnodes.end(), fnodes[j]));
              polytope_faces[n].push_back(pos);
            }
          }
        }

        double volume = region.intersect(polytope_nodes, polytope_faces);
        if (volume > 0.0) {
          ents.push_back(c);
          if (region.get_type()==AmanziGeometry::RegionType::LINE_SEGMENT) vol_fracs.push_back(volume);
          else vol_fracs.push_back(volume / mesh.getCellVolume(c));
        }
      }

    } else {
      // ind == FACE
      int nfaces = mesh.getNumEntities(Entity_kind::FACE, ptype);
      vol_fracs.reserve(nfaces);
      ents.reserve(nfaces);

      std::vector<AmanziGeometry::Point> polygon;

      for (int f=0; f!=nfaces; ++f) {
        auto polygon = mesh.getFaceCoordinates(f);
        double area = region.intersect(polygon);
        if (area > 0.0) {
          ents.push_back(f);
          vol_fracs.push_back(area / mesh.getFaceArea(f));
        }
      }
    }

  } else {
    ents = mesh.getSetEntities(region.get_name(), kind, ptype);
  }
  return ents;
}



} // namespace AmanziMesh
} // namespace Amanzi

/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Lipnikov Konstantin (lipnikov@lanl.gov)
      Rao Garimella (rao@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include <iterator>
#include <map>
#include <string>
#include <vector>

#include "Teuchos_CommHelpers.hpp"

#include "errors.hh"
#include "Point.hh"
#include "Region.hh"

#include "Mesh.hh"
#include "MeshDefs.hh"

namespace Amanzi {
namespace AmanziMesh {

void
Mesh::get_set_entities_box_vofs_(
  Teuchos::RCP<const AmanziGeometry::Region> region, const Entity_kind kind,
  const Parallel_type ptype, Kokkos::View<Entity_ID*>& setents,
  Kokkos::View<double*>* volume_fractions) const
{
  // AMANZI_ASSERT(volume_fractions != NULL);

  // setents->clear();
  // volume_fractions->clear();
  Kokkos::resize(*volume_fractions, 0);

  int sdim = space_dimension();
  int mdim = region->manifold_dimension();

  switch (kind) {
  case CELL: {
    std::string name = region->name() + "_cell";

    if (mesh_cache_.region_ids.find(name) !=  mesh_cache_.region_ids.end()) {
      Kokkos::resize(setents, mesh_cache_.region_ids[name].size());
      Kokkos::resize(*volume_fractions, mesh_cache_.region_vofs[name].size());

      assert(mesh_cache_.region_ids[name].size() == mesh_cache_.region_vofs.size());

      for (int i = 0; i < mesh_cache_.region_ids[name].size(); ++i) {
        setents(i) = mesh_cache_.region_ids[name][i];
        (*volume_fractions)(i) = mesh_cache_.region_vofs[name][i];
      }
    } else {
      int ncells = num_entities(CELL, ptype);
      double volume;

      Entity_ID_List cnodes, fnodes; 
      Kokkos::View<Entity_ID*> faces;
      Kokkos::View<int*> dirs;
      std::vector<AmanziGeometry::Point> polytope_nodes;
      std::vector<std::vector<int>> polytope_faces;

      size_t volume_fractions_size = 0;
      for (int c = 0; c < ncells; ++c) {
        cell_get_coordinates(c, polytope_nodes);

        if (space_dimension() == 3) {
          cell_get_nodes(c, cnodes);
          cell_get_faces_and_dirs(c, faces, dirs);
          int nfaces = faces.extent(0);

          polytope_faces.clear();
          polytope_faces.resize(nfaces);
          for (int n = 0; n < nfaces; ++n) {
            face_get_nodes(faces(n), fnodes);
            int nnodes = fnodes.size();

            for (int i = 0; i < nnodes; ++i) {
              int j = (dirs(n) > 0) ? i : nnodes - i - 1;
              int pos = 0;
              for (int p = 0; p < cnodes.size(); ++p) {
                if (cnodes[p] == fnodes[j]) {
                  pos = p;
                  break;
                }
              }
              polytope_faces[n].push_back(pos);
            }
          }
        }

        if ((volume = region->intersect(polytope_nodes, polytope_faces)) >
            0.0) {
          Kokkos::resize(setents, setents.extent(0) + 1);
          setents(setents.extent(0) - 1) = c;
          ;
          if (region->type() == AmanziGeometry::LINE_SEGMENT) {
            Kokkos::resize(*volume_fractions,
                           (*volume_fractions).extent(0) + 1);
            (*volume_fractions)(volume_fractions_size++) = volume;
          } else {
            Kokkos::resize(*volume_fractions,
                           (*volume_fractions).extent(0) + 1);
            (*volume_fractions)(volume_fractions_size++) =
              volume / cell_volume(c, false);
          }
        }
      }

      mesh_cache_.region_ids[name].resize(setents.extent(0));
      // assert(region_ids[name].size() == region_vofs.size());

      for (int i = 0; i < mesh_cache_.region_ids[name].size(); ++i) {
        mesh_cache_.region_ids[name][i] = setents(i);
      }
      for (int i = 0; i < mesh_cache_.region_vofs[name].size(); ++i) {
        mesh_cache_.region_vofs[name][i] = (*volume_fractions)(i);
      }
    }
    break;
  }

  case FACE: {
    std::string name = region->name() + "_face";

    if (mesh_cache_.region_ids.find(name) != mesh_cache_.region_ids.end()) {
      Kokkos::resize(setents, mesh_cache_.region_ids[name].size());
      assert(mesh_cache_.region_ids[name].size() == mesh_cache_.region_vofs.size());
      for (int i = 0; i < mesh_cache_.region_ids[name].size(); ++i) {
        setents(i) = mesh_cache_.region_ids[name][i];
        (*volume_fractions)(i) = mesh_cache_.region_vofs[name][i];
      }
    } else {
      int nfaces = num_entities(FACE, ptype);
      double area;

      std::vector<AmanziGeometry::Point> polygon;
      size_t volume_fractions_size = 0;
      for (int f = 0; f < nfaces; ++f) {
        face_get_coordinates(f, polygon);
        if ((area = region->intersect(polygon)) > 0.0) {
          Kokkos::resize(setents, setents.extent(0) + 1);
          setents(setents.extent(0) - 1) = f;
          (*volume_fractions)(volume_fractions_size++) = area / face_area(f);
        }
      }

      mesh_cache_.region_ids[name].resize(setents.extent(0));
      assert(mesh_cache_.region_ids[name].size() == mesh_cache_.region_vofs.size());
      for (int i = 0; i < mesh_cache_.region_ids[name].size(); ++i) {
        mesh_cache_.region_ids[name][i] = setents(i);
        mesh_cache_.region_vofs[name][i] = (*volume_fractions)(i);
      }
    }
    break;
  }

  case EDGE: {
    int nedges = num_entities(EDGE, ptype);

    for (int e = 0; e < nedges; ++e) {
      if (region->inside(mesh_cache_.edge_centroid(e))) {
        Kokkos::resize(setents, setents.extent(0) + 1);
        setents(setents.extent(0) - 1) = e;
      }
    }
    break;
  }

  case NODE: {
    int nnodes = num_entities(NODE, ptype);

    AmanziGeometry::Point xv(space_dimension());

    for (int v = 0; v < nnodes; ++v) {
      node_get_coordinates(v, &xv);
      if (region->inside(xv)) {
        Kokkos::resize(setents, setents.extent(0) + 1);
        setents(setents.extent(0) - 1) = v;
      }
    }
    break;
  }

  default:
    break;
  }

  // Check if no processor got any mesh entities
  int nents, nents_tmp = setents.extent(0);
  Teuchos::reduceAll(*get_comm(), Teuchos::REDUCE_SUM, 1, &nents_tmp, &nents);
  if (nents == 0) {
    Errors::Message msg;
    msg << "Could not retrieve any mesh entities for set \"" << region->name()
        << "\".\n";
    Exceptions::amanzi_throw(msg);
  }
}


//---------------------------------------------------------
// Generic implemnetation of set routines.
//---------------------------------------------------------
unsigned int
Mesh::get_set_size(const Set_ID setid, const Entity_kind kind,
                   const Parallel_type ptype) const
{
  Entity_ID_List ents;
  std::string setname = geometric_model()->FindRegion(setid)->name();

  get_set_entities(setname, kind, ptype, ents);

  return ents.size();
}


unsigned int
Mesh::get_set_size(const std::string setname, const Entity_kind kind,
                   const Parallel_type ptype) const
{
  Entity_ID_List setents;
  std::vector<double> vofs;

  get_set_entities_and_vofs(setname, kind, ptype, setents, &vofs);

  return setents.size();
}


void
Mesh::get_set_entities(const Set_ID setid, const Entity_kind kind,
                       const Parallel_type ptype,
                       std::vector<Entity_ID>& entids) const
{
  std::vector<double> vofs;
  std::string setname = geometric_model()->FindRegion(setid)->name();
  get_set_entities_and_vofs(setname, kind, ptype, entids, &vofs);
}

} // namespace AmanziMesh
} // namespace Amanzi

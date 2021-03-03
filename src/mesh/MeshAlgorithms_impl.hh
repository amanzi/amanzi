/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
      Ethan Coon (coonet@ornl.gov)
*/

#pragma once

#include "Point.hh"
#include "Geometry.hh"
#include "MeshDefs.hh"

namespace Amanzi {
namespace AmanziMesh {
namespace MeshAlgorithms {

template<class Mesh_type>
std::pair<double, AmanziGeometry::Point>
computeCellGeometry(const Mesh_type& mesh, const Entity_ID c)
{
  Entity_ID_List faces;
  std::vector<std::size_t> nfnodes;
  Entity_Direction_List fdirs;
  Point_List cfcoords;

  mesh.getCellFacesAndDirs(c, faces, &fdirs);
  nfnodes.resize(faces.size());

  for (int j = 0; j != faces.size(); ++j) {
    auto fcoords = mesh.face_coordinates(faces[j]);
    nfnodes[j] = fcoords.size();

    if (fdirs[j] == 1) {
      for (int k = 0; k != nfnodes[j]; ++k)
        cfcoords.emplace_back(fcoords[k]);
    } else {
      for (int k = nfnodes[j]-1; k >= 0; --k)
        cfcoords.emplace_back(fcoords[k]);
    }
  }

  auto ccoords = mesh.cell_coordinates(c);
  std::pair<double, AmanziGeometry::Point> vol_cent;
  AmanziGeometry::polyhed_get_vol_centroid(ccoords, faces.size(), nfnodes,
            cfcoords, &vol_cent.first, &vol_cent.second);
  return vol_cent;
}



template<class Mesh_type>
std::tuple<double,AmanziGeometry::Point,Point_List>
computeFaceGeometry(const Mesh_type& mesh, const Entity_ID f)
{
  if (mesh.get_manifold_dimension() == 3) {
    // 3D Elements with possibly curved faces
    Point_List fcoords = mesh.face_coordinates(f);

    double area;
    AmanziGeometry::Point centroid(3);
    AmanziGeometry::Point normal(3);
    AmanziGeometry::polygon_get_area_centroid_normal(fcoords, &area, &centroid, &normal);

    Entity_ID_List cellids;
    mesh.getFaceCells(f, Parallel_type::ALL, cellids);
    AMANZI_ASSERT(cellids.size() <= 2);

    Point_List normals(cellids.size());
    for (int i = 0; i < cellids.size(); i++) {
      Entity_ID_List cellfaceids;
      Entity_Direction_List cellfacedirs;
      int dir = 1;

      mesh.getCellFacesAndDirs(cellids[i], cellfaceids, &cellfacedirs);

      bool found = false;
      for (int j = 0; j < cellfaceids.size(); j++) {
        if (cellfaceids[j] == f) {
          found = true;
          dir = cellfacedirs[j];
          break;
        }
      }
      AMANZI_ASSERT(found);
      normals[i] = (dir == 1) ? normal : -normal;
    }
    return std::make_tuple(area, centroid, normals);

  } else if (mesh.get_manifold_dimension() == 2) {

    if (mesh.get_space_dimension() == 2) {
      // 2D mesh
      auto fcoords = mesh.face_coordinates(f);
      AMANZI_ASSERT(fcoords.size() == 2);

      AmanziGeometry::Point evec = fcoords[1] - fcoords[0];
      double area = std::sqrt(evec*evec);
      AmanziGeometry::Point centroid = 0.5 * (fcoords[0] + fcoords[1]);
      AmanziGeometry::Point normal(evec[1] - evec[0]);

      Entity_ID_List cellids;
      mesh.getFaceCells(f, Parallel_type::ALL, cellids);
      AMANZI_ASSERT(cellids.size() <= 2);

      Point_List normals(cellids.size());
      for (int i = 0; i < cellids.size(); i++) {
        Entity_ID_List cellfaceids;
        Entity_Direction_List cellfacedirs;
        int dir = 1;

        mesh.getCellFacesAndDirs(cellids[i], cellfaceids, &cellfacedirs);

        bool found = false;
        for (int j = 0; j < cellfaceids.size(); j++) {
          if (cellfaceids[j] == f) {
            found = true;
            dir = cellfacedirs[j];
            break;
          }
        }
        AMANZI_ASSERT(found);
        normals[i] = (dir == 1) ? normal : -normal;
      }
      return std::make_tuple(area, centroid, normals);

    } else {  // Surface mesh - cells are 2D, coordinates are 3D
      // Since the edge likely forms a discontinuity in the surface
      // (or may even be the intersection of several surfaces), we
      // have to compute an outward normal to the edge with respect to
      // each face
      auto fcoords = mesh.face_coordinates(f);
      AMANZI_ASSERT(fcoords.size() == 2);

      AmanziGeometry::Point evec = fcoords[1] - fcoords[0];
      double area = std::sqrt(evec*evec);
      AmanziGeometry::Point centroid = 0.5 * (fcoords[0] + fcoords[1]);

      Entity_ID_List cellids;
      mesh.getFaceCells(f, Parallel_type::ALL, cellids);
      Point_List normals(cellids.size());
      for (int i = 0; i < cellids.size(); i++) {
        Entity_ID_List cellfaceids;
        Entity_Direction_List cellfacedirs;
        mesh.getCellFacesAndDirs(cellids[i], cellfaceids, &cellfacedirs);

        bool found = false;
        for (int j = 0; j < cellfaceids.size(); j++) {
          if (cellfaceids[j] == f) {
            found = true;
            break;
          }
        }
        AMANZI_ASSERT(found);
        AmanziGeometry::Point cvec = fcoords[0] - mesh.getCellCentroid(cellids[i]);
        AmanziGeometry::Point trinormal = cvec^evec;
        AmanziGeometry::Point normal = evec^trinormal;

        double len = norm(normal);
        normal /= len;
        normal *= area;
        normals[i] = normal;  // Always an outward normal as calculated
      }
      return std::make_tuple(area, centroid, normals);
    }
  }

  Errors::Message msg;
  msg << "Invalid mesh argument to MeshAlgorithm: manifold_dim = " << mesh.get_manifold_dimension()
      << ", space_dim = " << mesh.get_space_dimension();
  Exceptions::amanzi_throw(msg);
}



template<class Mesh_type>
std::pair<AmanziGeometry::Point, AmanziGeometry::Point>
computeEdgeGeometry(const Mesh_type& mesh, const Entity_ID e)
{
  Entity_ID_List nodes;
  mesh.getEdgeNodes(e, nodes);
  AMANZI_ASSERT(nodes.size() == 2);
  auto x0 = mesh.getNodeCoordinate(nodes[0]);
  auto x1 = mesh.getNodeCoordinate(nodes[1]);
  return std::make_pair(x1-x0, (x1+x0)/2);
}

template<class Mesh_type>
void
computeBisectors(const Mesh_type& mesh, const Entity_ID c,
        const Entity_ID_List& faces, Point_List& bisectors)
{
  bisectors.resize(faces.size());
  for (int i = 0; i != faces.size(); ++i)
    bisectors[i] = mesh.getFaceCentroid(faces[i]) - mesh.getCellCentroid(c);
}


} // namespace MeshAlgorithms
} // namespace AmanziMesh
} // namspace Amanzi

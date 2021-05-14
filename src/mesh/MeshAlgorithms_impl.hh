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
Cell_type getCellType(const Mesh_type& mesh, const Entity_ID c)
{
  auto faces = mesh.getCellFaces(c);
  if (mesh.get_manifold_dimension() == 2) {
    switch (faces.size()) {
    case 3:
      return Cell_type::TRI;
      break;
    case 4:
      return Cell_type::QUAD;
      break;
    default:
      return Cell_type::POLYGON;
    }
  } else if (mesh.get_manifold_dimension() == 3) {
    int nquads = 0;
    for (const auto& f : faces) {
      Entity_ID_List fnodes;
      mesh.getFaceNodes(f, fnodes);
      if (fnodes.size() == 4) nquads++;
    }

    switch (faces.size()) {
      case 4:
        if (nquads == 0) return Cell_type::TET;
        else return Cell_type::POLYHED;
        break;
      case 5:
        if (nquads == 1) return Cell_type::PYRAMID;
        else if (nquads == 3) return Cell_type::PRISM;
        else return Cell_type::POLYHED;
        break;
      case 6:
        if (nquads == 6) return Cell_type::HEX;
        else return Cell_type::POLYHED;
        break;
      default:
        return Cell_type::POLYHED;
    }
  } else {
    Errors::Message msg;
    msg << "Mesh of manifold_dimension = " << mesh.get_manifold_dimension() << " not supported";
    Exceptions::amanzi_throw(msg);
  }
  return Cell_type::UNKNOWN;
}



template<class Mesh_type>
int getFaceDirectionInCell(const Mesh_type& mesh, const Entity_ID f, const Entity_ID c)
{
  Entity_ID_List cfaces;
  Entity_Direction_List dirs;
  mesh.getCellFacesAndDirs(c, cfaces, &dirs);
  int dir = 1;
  for (int j=0; j!=cfaces.size(); ++j) {
    if (cfaces[j] == f) {
      dir = dirs[j];
      break;
    }
  }
  return dir;
}

template<class Mesh_type>
std::pair<double, AmanziGeometry::Point>
computeCellGeometry(const Mesh_type& mesh, const Entity_ID c)
{
  if (mesh.get_manifold_dimension() == 2) {
    auto ccoords = mesh.getCellCoordinates(c);
    auto vol_cent = std::make_pair((double)0,
            AmanziGeometry::Point(mesh.get_space_dimension()));
    AmanziGeometry::Point normal(mesh.get_space_dimension());
    AmanziGeometry::polygon_get_area_centroid_normal(ccoords,
            &vol_cent.first, &vol_cent.second, &normal);
    return vol_cent;
  } else {

    Entity_ID_List faces;
    std::vector<std::size_t> nfnodes;
    Entity_Direction_List fdirs;
    Point_List cfcoords;

    mesh.getCellFacesAndDirs(c, faces, &fdirs);
    nfnodes.resize(faces.size());

    for (int j = 0; j != faces.size(); ++j) {
      auto fcoords = mesh.getFaceCoordinates(faces[j]);
      nfnodes[j] = fcoords.size();

      if (fdirs[j] == 1) {
        for (int k = 0; k != nfnodes[j]; ++k)
          cfcoords.emplace_back(fcoords[k]);
      } else {
        for (int k = nfnodes[j]-1; k >= 0; --k)
          cfcoords.emplace_back(fcoords[k]);
      }
    }

    auto ccoords = mesh.getCellCoordinates(c);
    auto vol_cent = std::make_pair((double)0,
            AmanziGeometry::Point(mesh.get_space_dimension()));
    AmanziGeometry::polyhed_get_vol_centroid(ccoords, faces.size(), nfnodes,
            cfcoords, &vol_cent.first, &vol_cent.second);
    return vol_cent;
  }
}


template<class Mesh_type>
std::tuple<double,AmanziGeometry::Point,Point_List>
computeFaceGeometry(const Mesh_type& mesh, const Entity_ID f)
{
  if (mesh.get_manifold_dimension() == 3) {
    // 3D Elements with possibly curved faces
    Point_List fcoords = mesh.getFaceCoordinates(f);

    double area;
    AmanziGeometry::Point centroid(3);
    AmanziGeometry::Point normal(3);
    AmanziGeometry::polygon_get_area_centroid_normal(fcoords, &area, &centroid, &normal);

    Entity_ID_List fcells;
    mesh.getFaceCells(f, Parallel_type::ALL, fcells);
    Point_List normals(fcells.size(), normal);

    for (int i=0; i!=fcells.size(); ++i) {
      int dir = MeshAlgorithms::getFaceDirectionInCell(mesh, f, fcells[i]);
      normals[i] = dir * normals[i];
    }
    return std::make_tuple(area, centroid, normals);

  } else if (mesh.get_manifold_dimension() == 2) {

    if (mesh.get_space_dimension() == 2) {
      // 2D mesh
      auto fcoords = mesh.getFaceCoordinates(f);
      AMANZI_ASSERT(fcoords.size() == 2);

      AmanziGeometry::Point evec = fcoords[1] - fcoords[0];
      double area = std::sqrt(evec*evec);
      AmanziGeometry::Point centroid = 0.5 * (fcoords[0] + fcoords[1]);
      AmanziGeometry::Point normal(evec[1], -evec[0]);

      // only one normal needed
      Entity_ID_List fcells;
      mesh.getFaceCells(f, Parallel_type::ALL, fcells);
      Point_List normals(fcells.size(), normal);
      for (int i=0; i!=fcells.size(); ++i) {
        int dir = MeshAlgorithms::getFaceDirectionInCell(mesh, f, fcells[i]);
        normals[i] = dir * normals[i];
      }
      return std::make_tuple(area, centroid, normals);

    } else {
      // Surface mesh - cells are 2D, coordinates are 3D
      // Since the face likely forms a discontinuity in the surface
      // (or may even be the intersection of several surfaces), we
      // have to compute an outward normal to the face with respect to
      // each cell.
      auto fcoords = mesh.getFaceCoordinates(f);
      AMANZI_ASSERT(fcoords.size() == 2);

      AmanziGeometry::Point evec = fcoords[1] - fcoords[0];
      double area = AmanziGeometry::norm(evec);
      AmanziGeometry::Point centroid = 0.5 * (fcoords[0] + fcoords[1]);

      Entity_ID_List cellids;
      mesh.getFaceCells(f, Parallel_type::ALL, cellids);
      Point_List normals(cellids.size());
      for (int i = 0; i < cellids.size(); i++) {
        AmanziGeometry::Point cvec = fcoords[0] - mesh.getCellCentroid(cellids[i]);
        AmanziGeometry::Point trinormal = cvec^evec;
        AmanziGeometry::Point normal = evec^trinormal;

        double len = norm(normal);
        normal *= (area/len);
      }
      return std::make_tuple(area, centroid, normals);
    }
  }

  Errors::Message msg;
  msg << "Invalid mesh argument to MeshAlgorithm: manifold_dim = " << mesh.get_manifold_dimension()
      << ", space_dim = " << mesh.get_space_dimension();
  Exceptions::amanzi_throw(msg);
  return std::make_tuple(0, AmanziGeometry::Point(), std::vector<AmanziGeometry::Point>());
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


template<class Mesh_type>
void debugCell(const Mesh_type& mesh, const Entity_ID c)
{
  auto cc = mesh.getCellCentroid(c);
  std::stringstream stream;
  stream << "Debug cell: LID " << c << " Rank " << mesh.get_comm()->MyPID()
         << " GID " << mesh.getEntityGID(Entity_kind::CELL, c)
         << " of type " << to_string(mesh.getCellType(c)) << std::endl
          << "  centroid = " << mesh.getCellCentroid(c) << std::endl;
  Entity_ID_List faces;
  Entity_Direction_List f_dirs;
  mesh.getCellFacesAndDirs(c, faces, &f_dirs);
  for (int fi=0; fi!=faces.size(); ++fi) {
    auto fc = mesh.getFaceCentroid(faces[fi]);
    stream << "  face " << faces[fi]
           << " GID " << mesh.getEntityGID(Entity_kind::FACE, faces[fi])
           << " dir " << f_dirs[fi] << " centroid " << fc << " normal " << mesh.getFaceNormal(faces[fi]) << " bisector " <<  fc - cc << std::endl;
  }
  std::cout << stream.str();
}



} // namespace MeshAlgorithms
} // namespace AmanziMesh
} // namspace Amanzi

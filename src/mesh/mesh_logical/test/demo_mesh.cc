/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

//
// Set of functions that create demo meshes for testing.
//
// ------------------------------------------------------------------

#include <cmath>

#include "AmanziComm.hh"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "Mesh_MSTK.hh"
#include "MeshLogicalFactory.hh"
#include "Geometry.hh"

#include "demo_mesh.hh"

namespace Amanzi {
namespace Testing {

// Single segment of plant, manually generated
Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical>
demoMeshLogicalSegmentRegularManual()
{
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  auto comm = Amanzi::getDefaultComm();
  const int nproc(comm->getSize());
  const int me(comm->getRank());

  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3));

  // Create the mesh:
  std::vector<double> cell_volumes_;
  cell_volumes_.push_back(.25);
  cell_volumes_.push_back(.25);
  cell_volumes_.push_back(.25);
  cell_volumes_.push_back(.25);

  std::vector<std::vector<int>> face_cell_list(5);
  face_cell_list[0].push_back(0);
  face_cell_list[1].push_back(0);
  face_cell_list[1].push_back(1);
  face_cell_list[2].push_back(1);
  face_cell_list[2].push_back(2);
  face_cell_list[3].push_back(2);
  face_cell_list[3].push_back(3);
  face_cell_list[4].push_back(3);

  std::vector<std::vector<double>> face_cell_lengths;
  face_cell_lengths.resize(5);
  face_cell_lengths[0].push_back(0.125);
  face_cell_lengths[1].resize(2, 0.125);
  face_cell_lengths[2].resize(2, 0.125);
  face_cell_lengths[3].resize(2, 0.125);
  face_cell_lengths[4].push_back(0.125);

  std::vector<Amanzi::AmanziGeometry::Point> face_area_normals;
  Amanzi::AmanziGeometry::Point normal(3);
  normal[0] = 1.0;
  normal[1] = 0.0;
  normal[2] = 0.0;
  face_area_normals.resize(5, normal);

  Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical> mesh =
    Teuchos::rcp(new MeshLogical(comm,
                                 cell_volumes_,
                                 face_cell_list,
                                 face_cell_lengths,
                                 face_area_normals));
  mesh->set_geometric_model(gm);
  return mesh;
};


// Single segment of plant with irregular mesh spacing
Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical>
demoMeshLogicalSegmentIrregularManual()
{
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = Amanzi::getDefaultComm();
  const int nproc(comm->getSize());
  const int me(comm->getRank());

  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3));

  // Create the mesh:
  std::vector<double> cell_volumes_;
  cell_volumes_.push_back(2);
  cell_volumes_.push_back(3);
  cell_volumes_.push_back(4);

  std::vector<std::vector<int>> face_cell_list;
  face_cell_list.resize(4);
  face_cell_list[0].push_back(0);
  face_cell_list[1].push_back(0);
  face_cell_list[1].push_back(1);
  face_cell_list[2].push_back(1);
  face_cell_list[2].push_back(2);
  face_cell_list[3].push_back(2);

  std::vector<std::vector<double>> face_cell_lengths;
  face_cell_lengths.resize(4);
  face_cell_lengths[0].push_back(0.5);
  face_cell_lengths[1].push_back(0.5);
  face_cell_lengths[1].push_back(.75);
  face_cell_lengths[2].push_back(.75);
  face_cell_lengths[2].push_back(1.);
  face_cell_lengths[3].push_back(1.);

  std::vector<Amanzi::AmanziGeometry::Point> face_area_normals;
  Amanzi::AmanziGeometry::Point normal(3);
  normal[0] = 2.0;
  normal[1] = 0.0;
  normal[2] = 0.0;
  face_area_normals.resize(4, normal);

  Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical> mesh =
    Teuchos::rcp(new MeshLogical(comm,
                                 cell_volumes_,
                                 face_cell_list,
                                 face_cell_lengths,
                                 face_area_normals));
  mesh->set_geometric_model(gm);
  return mesh;
}

// Single coarse root to 1 meter, then branches to 4 fine roots
Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical>
demoMeshLogicalYManual()
{
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = Amanzi::getDefaultComm();
  const int nproc(comm->getSize());
  const int me(comm->getRank());

  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3));

  std::vector<double> cell_vols;
  std::vector<std::vector<Entity_ID>> face_cells;
  std::vector<std::vector<double>> face_cell_lengths;
  std::vector<Point> face_normals;
  std::vector<Point> cell_centroids_;

  double half_cell_coarse = 0.25;
  double half_cell_fine = 0.25 * std::sqrt(2);
  double cross_section_coarse = 1.e-4;
  double cross_section_fine = cross_section_coarse / 4.;

  // Structure: set up centroids
  // ---------------------------------
  // coarse root
  Point stem_root_intersection(0., 0., 0.);
  Point down(0., 0., -half_cell_coarse);
  cell_centroids_.push_back(stem_root_intersection + down); // c0
  cell_centroids_.push_back(cell_centroids_[0] + 2 * down); // c1
  cell_centroids_.push_back(cell_centroids_[1] + 2 * down); // c2

  // branch cell
  Point branch = cell_centroids_[cell_centroids_.size() - 1];

  // fine root 1
  Point root_normal1(1., 0., -1.);
  root_normal1 *= half_cell_fine / norm(root_normal1);
  cell_centroids_.push_back(branch + 2 * root_normal1);             // c3
  cell_centroids_.push_back(cell_centroids_[3] + 2 * root_normal1); // c4

  // fine root 2
  Point root_normal2(-1., 0., -1.);
  root_normal2 *= half_cell_fine / norm(root_normal2);
  cell_centroids_.push_back(branch + 2 * root_normal2);             // c5
  cell_centroids_.push_back(cell_centroids_[5] + 2 * root_normal2); // c6

  // fine root 3
  Point root_normal3(0., 1., -1.);
  root_normal3 *= half_cell_fine / norm(root_normal3);
  cell_centroids_.push_back(branch + 2 * root_normal3);             // c7
  cell_centroids_.push_back(cell_centroids_[7] + 2 * root_normal3); // c8

  // fine root 4
  Point root_normal4(0., -1., -1.);
  root_normal4 *= half_cell_fine / norm(root_normal4);
  cell_centroids_.push_back(branch + 2 * root_normal4);             // c9
  cell_centroids_.push_back(cell_centroids_[9] + 2 * root_normal4); // c10


  // cell volumes
  // ---------------------------------
  for (int lcv = 0; lcv != 2; ++lcv) {
    cell_vols.push_back(cross_section_coarse * 2 * half_cell_coarse); // coarse
  }
  cell_vols.push_back(cross_section_coarse * half_cell_coarse +
                      4 * cross_section_fine * half_cell_fine); // junction cell
  for (int lcv = 0; lcv != 8; ++lcv) {
    cell_vols.push_back(cross_section_fine * 2 * half_cell_fine); // coarse
  }

  // faces
  face_cells.resize(15);
  face_cell_lengths.resize(15);
  face_normals.resize(15);

  // renormalize by area
  down *= cross_section_coarse / half_cell_coarse;
  root_normal1 *= cross_section_fine / half_cell_fine;
  root_normal2 *= cross_section_fine / half_cell_fine;
  root_normal3 *= cross_section_fine / half_cell_fine;
  root_normal4 *= cross_section_fine / half_cell_fine;

  // coarse root -- 3 faces
  // -- f0
  face_cells[0].push_back(0);
  face_cell_lengths[0].push_back(half_cell_coarse);
  face_normals[0] = -down;
  // -- f1
  face_cells[1].push_back(0);
  face_cells[1].push_back(1);
  face_cell_lengths[1].push_back(half_cell_coarse);
  face_cell_lengths[1].push_back(half_cell_coarse);
  face_normals[1] = down;
  // -- f2
  face_cells[2].push_back(1);
  face_cells[2].push_back(2);
  face_cell_lengths[2].push_back(half_cell_coarse);
  face_cell_lengths[2].push_back(half_cell_coarse);
  face_normals[2] = down;

  // fine root 1, 3 faces
  // -- f3
  face_cells[5].push_back(2);
  face_cells[5].push_back(3);
  face_cell_lengths[5].push_back(half_cell_fine);
  face_cell_lengths[5].push_back(half_cell_fine);
  face_normals[5] = root_normal1;
  // -- f4
  face_cells[3].push_back(3);
  face_cells[3].push_back(4);
  face_cell_lengths[3].push_back(half_cell_fine);
  face_cell_lengths[3].push_back(half_cell_fine);
  face_normals[3] = root_normal1;
  // -- f5
  face_cells[4].push_back(4);
  face_cell_lengths[4].push_back(half_cell_fine);
  face_normals[4] = root_normal1;

  // fine root 2, 3 faces
  // -- f6
  face_cells[8].push_back(2);
  face_cells[8].push_back(5);
  face_cell_lengths[8].push_back(half_cell_fine);
  face_cell_lengths[8].push_back(half_cell_fine);
  face_normals[8] = root_normal2;
  // -- f7
  face_cells[6].push_back(5);
  face_cells[6].push_back(6);
  face_cell_lengths[6].push_back(half_cell_fine);
  face_cell_lengths[6].push_back(half_cell_fine);
  face_normals[6] = root_normal2;
  // -- f8
  face_cells[7].push_back(6);
  face_cell_lengths[7].push_back(half_cell_fine);
  face_normals[7] = root_normal2;

  // fine root 1, 3 faces
  // -- f9
  face_cells[11].push_back(2);
  face_cells[11].push_back(7);
  face_cell_lengths[11].push_back(half_cell_fine);
  face_cell_lengths[11].push_back(half_cell_fine);
  face_normals[11] = root_normal3;
  // -- f10
  face_cells[9].push_back(7);
  face_cells[9].push_back(8);
  face_cell_lengths[9].push_back(half_cell_fine);
  face_cell_lengths[9].push_back(half_cell_fine);
  face_normals[9] = root_normal3;
  // -- f11
  face_cells[10].push_back(8);
  face_cell_lengths[10].push_back(half_cell_fine);
  face_normals[10] = root_normal3;

  // fine root 1, 3 faces
  // -- f12
  face_cells[14].push_back(2);
  face_cells[14].push_back(9);
  face_cell_lengths[14].push_back(half_cell_fine);
  face_cell_lengths[14].push_back(half_cell_fine);
  face_normals[14] = root_normal4;
  // -- f13
  face_cells[12].push_back(9);
  face_cells[12].push_back(10);
  face_cell_lengths[12].push_back(half_cell_fine);
  face_cell_lengths[12].push_back(half_cell_fine);
  face_normals[12] = root_normal4;
  // -- f14
  face_cells[13].push_back(10);
  face_cell_lengths[13].push_back(half_cell_fine);
  face_normals[13] = root_normal4;

  // make the mesh
  Teuchos::RCP<MeshLogical> m = Teuchos::rcp(new MeshLogical(comm,
                                                             cell_vols,
                                                             face_cells,
                                                             face_cell_lengths,
                                                             face_normals,
                                                             &cell_centroids_));

  // make sets
  // -- coarse roots
  Entity_ID_List coarse_cells;
  for (int lcv = 0; lcv != 3; ++lcv) coarse_cells.push_back(lcv);
  Teuchos::RCP<AmanziGeometry::RegionEnumerated> coarse_rgn =
    Teuchos::rcp(new AmanziGeometry::RegionEnumerated(
      "coarse_root", gm->RegionSize(), "CELL", coarse_cells));
  gm->AddRegion(coarse_rgn);

  Entity_ID_List fine_cells;
  for (int lcv = 3; lcv != 11; ++lcv) fine_cells.push_back(lcv);
  Teuchos::RCP<AmanziGeometry::RegionEnumerated> fine_rgn =
    Teuchos::rcp(new AmanziGeometry::RegionEnumerated(
      "fine_root", gm->RegionSize(), "CELL", fine_cells));
  gm->AddRegion(fine_rgn);

  m->set_geometric_model(gm);
  return m;
}


// Single coarse root to 1 meter, then branches to 4 fine roots
Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical>
demoMeshLogicalYFromXML(const std::string& meshname)
{
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = Amanzi::getDefaultComm();
  const int nproc(comm->getSize());
  const int me(comm->getRank());

  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3));
  MeshLogicalFactory fac(comm, gm);

  // load the xml
  std::string xmlFileName = "test/demo_mesh_Y.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  return fac.Create(plist.sublist(meshname));
}

// helper class to get indices from a regular mesh
RegularMeshCellFromCoordFunctor::RegularMeshCellFromCoordFunctor(
  const Amanzi::AmanziGeometry::Point& X0,
  const Amanzi::AmanziGeometry::Point& X1, int nx, int ny, int nz)
  : X0_(X0), X1_(X1), nx_(nx), ny_(ny), nz_(nz), dX_(3)
{
  dX_[0] = (X1_[0] - X0_[0]) / nx_;
  dX_[1] = (X1_[1] - X0_[1]) / ny_;
  dX_[2] = (X1_[2] - X0_[2]) / nz_;
}

Amanzi::AmanziMesh::Entity_ID
RegularMeshCellFromCoordFunctor::
operator()(const Amanzi::AmanziGeometry::Point& p)
{
  Amanzi::AmanziGeometry::Point dp = p - X0_;

  int i = round(dp[0] / dX_[0] - .5);
  int j = round(dp[1] / dX_[1] - .5);
  int k = round(dp[2] / dX_[2] - .5);
  return k + j * nz_ + i * nz_ * ny_;
};


Teuchos::RCP<Amanzi::AmanziMesh::MeshEmbeddedLogical>
demoMeshLogicalYEmbedded()
{
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = Amanzi::getDefaultComm();
  const int nproc(comm->getSize());
  const int me(comm->getRank());

  Teuchos::RCP<MeshLogical> m_log = demoMeshLogicalYManual();

  Point X0(-1.75, -1.75, -3.);
  Point X1(1.75, 1.75, 0.);
  int nx = 7, ny = 7, nz = 6;

  Teuchos::RCP<Mesh_MSTK> m_bg =
    Teuchos::rcp(new Mesh_MSTK(X0[0],
                               X0[1],
                               X0[2],
                               X1[0],
                               X1[1],
                               X1[2],
                               nx,
                               ny,
                               nz,
                               comm,
                               m_log->geometric_model(),
                               Teuchos::null,
                               true,
                               false));

  // make the new connections, 1 per logical cell
  int ncells_log = m_log->num_entities(CELL, Parallel_type::ALL);
  std::vector<Entity_ID_List> face_cell_ids_(ncells_log);

  RegularMeshCellFromCoordFunctor cellID(X0, X1, nx, ny, nz);
  for (int c = 0; c != ncells_log; ++c) {
    int bg_c = cellID(m_log->cell_centroid(c));
    face_cell_ids_[c].push_back(c);
    face_cell_ids_[c].push_back(bg_c);
  }

  // cross-sectional areas and lengths given by root class
  // relatively made up numbers, this could be done formally if we wanted to.
  Entity_ID_List coarse_roots, fine_roots;
  m_log->get_set_entities(
    "coarse_root", CELL, Parallel_type::ALL, &coarse_roots);
  m_log->get_set_entities("fine_root", CELL, Parallel_type::ALL, &fine_roots);

  std::vector<std::vector<double>> face_cell_lengths(ncells_log);
  std::vector<AmanziGeometry::Point> face_area_normals(ncells_log);
  Point ihat(1., 0., 0.); // just must have 0 z to avoid gravity
  for (Entity_ID_List::iterator c = coarse_roots.begin();
       c != coarse_roots.end();
       ++c) {
    face_area_normals[*c] = 2.e-2 * ihat;
    face_cell_lengths[*c].push_back(1.e-2);
    face_cell_lengths[*c].push_back(0.5 / 4.);
  }
  for (Entity_ID_List::iterator c = fine_roots.begin(); c != fine_roots.end();
       ++c) {
    face_area_normals[*c] = 1.e-2 * ihat;
    face_cell_lengths[*c].push_back(5.e-3);
    face_cell_lengths[*c].push_back(0.5 / 4.);
  }

  Teuchos::RCP<MeshEmbeddedLogical> m = Teuchos::rcp(new MeshEmbeddedLogical(
    comm, m_bg, m_log, face_cell_ids_, face_cell_lengths, face_area_normals));
  return m;
}


} // namespace Testing
} // namespace Amanzi

/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//
// Set of functions that create demo meshes for testing.
//
// ------------------------------------------------------------------

#include <cmath>

#include "AmanziComm.hh"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "Mesh_MSTK.hh"
#include "Mesh.hh"
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

  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3));

  // Create the mesh:
  std::vector<double> cell_volumes(4, 0.25);
  std::vector<double> face_areas(5, 1.);

  std::vector<std::vector<int> > face_cell_list(5);
  face_cell_list[0] = { 0 };
  face_cell_list[1] = { 0, 1 };
  face_cell_list[2] = { 1, 2 };
  face_cell_list[3] = { 2, 3 };
  face_cell_list[4] = { 3 };

  std::vector<std::vector<int> > face_cell_dirs(5);
  face_cell_dirs[0] = { 1 };
  face_cell_dirs[1] = { -1, 1 };
  face_cell_dirs[2] = { -1, 1 };
  face_cell_dirs[3] = { -1, 1 };
  face_cell_dirs[4] = { -1 };

  Point bisector(3);
  bisector[0] = 0.125;
  bisector[1] = 0.0;
  bisector[2] = 0.0;

  std::vector<std::vector<Point> > face_cell_bisectors(5);
  face_cell_bisectors[0] = std::vector<Point>{ bisector };
  face_cell_bisectors[1] = std::vector<Point>{ -bisector, bisector };
  face_cell_bisectors[2] = std::vector<Point>{ -bisector, bisector };
  face_cell_bisectors[3] = std::vector<Point>{ -bisector, bisector };
  face_cell_bisectors[4] = std::vector<Point>{ -bisector };

  Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical> mesh =
    Teuchos::rcp(new MeshLogical(comm,
            face_cell_list,
            face_cell_dirs,
            cell_volumes,
            face_areas,
            face_cell_bisectors));
  mesh->setGeometricModel(gm);
  return mesh;
};



// Single segment of plant with irregular mesh spacing
Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical>
demoMeshLogicalSegmentIrregularManual() {
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = Amanzi::getDefaultComm();

  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3));

  // Create the mesh:
  std::vector<double> cell_volumes = { 2, 3, 4 };
  std::vector<double> face_areas = { 2, 2, 2, 2 };

  std::vector<std::vector<int> > face_cell_list(4);
  face_cell_list[0] = { 0 };
  face_cell_list[1] = { 0, 1 };
  face_cell_list[2] = { 1, 2 };
  face_cell_list[3] = { 2 };

  std::vector<std::vector<int> > face_cell_dirs(4);
  face_cell_dirs[0] = { 1 };
  face_cell_dirs[1] = { -1, 1 };
  face_cell_dirs[2] = { -1, 1 };
  face_cell_dirs[3] = { -1 };

  Amanzi::AmanziGeometry::Point normal(3);
  normal[0] = 1.0;
  normal[1] = 0.0;
  normal[2] = 0.0;

  std::vector<std::vector<Point> > face_cell_bisectors;
  face_cell_bisectors.resize(4);
  face_cell_bisectors[0] = { .5 * normal };
  face_cell_bisectors[1] = { -.5 * normal, .75 * normal };
  face_cell_bisectors[2] = { -.75 * normal, normal };
  face_cell_bisectors[3] = { -normal };

  Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical> mesh =
    Teuchos::rcp(new MeshLogical(comm,
            face_cell_list,
            face_cell_dirs,
            cell_volumes,
            face_areas,
            face_cell_bisectors));
  mesh->setGeometricModel(gm);
  return mesh;
}

// Single coarse root to 1 meter, then branches to 4 fine roots
Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical>
demoMeshLogicalYManual() {
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = Amanzi::getDefaultComm();

  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3));

  std::vector<double> cell_vols;
  std::vector<double> face_areas;
  std::vector<std::vector<Entity_ID> > face_cells;
  std::vector<std::vector<int> > face_cell_dirs;
  std::vector<std::vector<Point> > face_cell_bisectors;
  std::vector<Point> cell_centroids;

  double half_cell_coarse = 0.5;
  double half_cell_fine = 0.5 * 0.75;
  double cross_section_coarse = 1.e-4;
  double cross_section_fine = cross_section_coarse / 4.;

  // Structure: set up centroids
  // ---------------------------------
  // coarse root
  Point stem_root_intersection(0., 0., 0.);
  Point down(0.,0.,-half_cell_coarse);
  cell_centroids.push_back(stem_root_intersection + down); // c0
  cell_centroids.push_back(cell_centroids[0] + 2*down); // c1
  cell_centroids.push_back(cell_centroids[1] + 2*down); // c2

  // branch cell
  Point branch = cell_centroids[cell_centroids.size()-1] + down;

  // fine root 1
  Point root_normal1(1., 0., -1.);
  root_normal1 *= half_cell_fine / norm(root_normal1);
  cell_centroids.push_back(branch + 2*root_normal1); // c3
  cell_centroids.push_back(cell_centroids[3] + 2*root_normal1); // c4

  // fine root 2
  Point root_normal2(-1., 0., -1.);
  root_normal2 *= half_cell_fine / norm(root_normal2);
  cell_centroids.push_back(branch + 2*root_normal2); // c5
  cell_centroids.push_back(cell_centroids[5] + 2*root_normal2); // c6

  // fine root 3
  Point root_normal3(0., 1., -1.);
  root_normal3 *= half_cell_fine / norm(root_normal3);
  cell_centroids.push_back(branch + 2*root_normal3); // c7
  cell_centroids.push_back(cell_centroids[7] + 2*root_normal3); // c8

  // fine root 4
  Point root_normal4(0., -1., -1.);
  root_normal4 *= half_cell_fine / norm(root_normal4);
  cell_centroids.push_back(branch + 2*root_normal4); // c9
  cell_centroids.push_back(cell_centroids[9] + 2*root_normal4); // c10


  // cell volumes
  // ---------------------------------
  for (int lcv=0; lcv!=2; ++lcv) {
    cell_vols.push_back(cross_section_coarse * 2*half_cell_coarse); // coarse
  }
  cell_vols.push_back(cross_section_coarse*half_cell_coarse + 4*cross_section_fine*half_cell_coarse); // junction cell
  for (int lcv=0; lcv!=8; ++lcv) {
    cell_vols.push_back(cross_section_fine * 2*half_cell_fine); // coarse
  }

  // faces
  face_cells.resize(15);
  face_cell_dirs.resize(15);
  face_cell_bisectors.resize(15);
  face_areas.resize(15);

  // coarse root -- 3 faces
  // -- f0
  face_cells[0] = { 0 };
  face_cell_bisectors[0] = { -down };
  face_areas[0] = cross_section_coarse;
  face_cell_dirs[0] = { 1 };
  // -- f1
  face_cells[1] = { 0, 1 };
  face_cell_bisectors[1] = { down, -down };
  face_areas[1] = cross_section_coarse;
  face_cell_dirs[1] = { -1, 1 };
  // -- f2
  face_cells[2] = { 1, 2 };
  face_cell_bisectors[2] = { down, -down };
  face_areas[2] = cross_section_coarse;
  face_cell_dirs[2] = { -1, 1 };

  // fine root 1, 3 faces
  // -- f3
  face_cells[3] = { 2, 3 };
  face_cell_bisectors[3] = { down, -root_normal1 };
  face_areas[3] = cross_section_fine;
  face_cell_dirs[3] = { -1, 1 };
  // -- f4
  face_cells[4] = { 3, 4 };
  face_cell_bisectors[4] = { root_normal1, -root_normal1 };
  face_areas[4] = cross_section_fine;
  face_cell_dirs[4] = { -1, 1 };
  // -- f5
  face_cells[5] = { 4 };
  face_cell_bisectors[5] = { root_normal1 };
  face_areas[5] = cross_section_fine;
  face_cell_dirs[5] = { -1 };

  // fine root 2, 3 faces
  // -- f6
  face_cells[6] = { 2, 5 };
  face_cell_bisectors[6] = { down, -root_normal2 };
  face_areas[6] = cross_section_fine;
  face_cell_dirs[6] = { -1, 1 };
  // -- f7
  face_cells[7] = { 5, 6 };
  face_cell_bisectors[7] = { root_normal2, -root_normal2 };
  face_areas[7] = cross_section_fine;
  face_cell_dirs[7] = { -1, 1 };
  // -- f8
  face_cells[8] = { 6 };
  face_cell_bisectors[8] = { root_normal2 };
  face_areas[8] = cross_section_fine;
  face_cell_dirs[8] = { -1 };

  // fine root 3, 3 faces
  // -- f9
  face_cells[9] = { 2, 7 };
  face_cell_bisectors[9] = { down, -root_normal3 };
  face_areas[9] = cross_section_fine;
  face_cell_dirs[9] = { -1, 1 };
  // -- f10
  face_cells[10] = { 7, 8 };
  face_cell_bisectors[10] = { root_normal3, -root_normal3 };
  face_areas[10] = cross_section_fine;
  face_cell_dirs[10] = { -1, 1 };
  // -- f11
  face_cells[11] = { 8 };
  face_cell_bisectors[11] = { root_normal3 };
  face_areas[11] = cross_section_fine;
  face_cell_dirs[11] = { -1 };

  // fine root 4, 3 faces
  // -- f12
  face_cells[12] = { 2, 9 };
  face_cell_bisectors[12] = { down, -root_normal4 };
  face_areas[12] = cross_section_fine;
  face_cell_dirs[12] = { -1, 1 };
  // -- f13
  face_cells[13] = { 9, 10 };
  face_cell_bisectors[13] = { root_normal4, -root_normal4 };
  face_areas[13] = cross_section_fine;
  face_cell_dirs[13] = { -1, 1 };
  // -- f14
  face_cells[14] = { 10 };
  face_cell_bisectors[14] = { root_normal4 };
  face_areas[14] = cross_section_fine;
  face_cell_dirs[14] = { -1 };

  // make the mesh
  Teuchos::RCP<MeshLogical> m = Teuchos::rcp(new MeshLogical(comm, face_cells, face_cell_dirs,
          cell_vols, face_areas, face_cell_bisectors, &cell_centroids));

  // make sets
  // -- coarse roots
  Entity_ID_List coarse_cells;
  for (int lcv=0; lcv!=3; ++lcv) coarse_cells.push_back(lcv);
  Teuchos::RCP<AmanziGeometry::RegionEnumerated> coarse_rgn =
      Teuchos::rcp(new AmanziGeometry::RegionEnumerated("coarse_root",
              gm->size(), "CELL", coarse_cells));
  gm->AddRegion(coarse_rgn);

  Entity_ID_List fine_cells;
  for (int lcv=3; lcv!=11; ++lcv) fine_cells.push_back(lcv);
  Teuchos::RCP<AmanziGeometry::RegionEnumerated> fine_rgn =
      Teuchos::rcp(new AmanziGeometry::RegionEnumerated("fine_root",
              gm->size(), "CELL", fine_cells));
  gm->AddRegion(fine_rgn);

  m->setGeometricModel(gm);
  return m;
}


// Single coarse root to 1 meter, then branches to 4 fine roots
Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical>
demoMeshLogicalFromXML(const std::string& meshname) {
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = Amanzi::getDefaultComm();

  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3));
  MeshLogicalFactory fac(comm, gm);

  // load the xml
  std::string xmlFileName = "test/demo_mesh.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  return fac.Create(plist.sublist(meshname));
}


// helper class to get indices from a regular mesh
RegularMeshCellFromCoordFunctor::RegularMeshCellFromCoordFunctor(
				 const Amanzi::AmanziGeometry::Point& X0,
				 const Amanzi::AmanziGeometry::Point& X1,
				 int nx, int ny, int nz)
  : X0_(X0), X1_(X1), nx_(nx), ny_(ny), nz_(nz), dX_(3) {
  dX_[0] = (X1_[0] - X0_[0]) / nx_;
  dX_[1] = (X1_[1] - X0_[1]) / ny_;
  dX_[2] = (X1_[2] - X0_[2]) / nz_;
}

Amanzi::AmanziMesh::Entity_ID
RegularMeshCellFromCoordFunctor::operator()(const Amanzi::AmanziGeometry::Point& p) {
  Amanzi::AmanziGeometry::Point dp = p - X0_;

  int i = round(dp[0]/dX_[0] - .5);
  int j = round(dp[1]/dX_[1] - .5);
  int k = round(dp[2]/dX_[2] - .5);
  return k + j*nz_ + i*nz_*ny_;
};


Teuchos::RCP<Amanzi::AmanziMesh::MeshEmbeddedLogical>
demoMeshLogicalYEmbedded() {
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = Amanzi::getDefaultComm();

  Teuchos::RCP<MeshFramework> m_log_fw = demoMeshLogicalYManual();
  Teuchos::RCP<Mesh> m_log = Teuchos::rcp(new Mesh(m_log_fw));

  Point X0(-1.75, -1.75, -3.);
  Point X1(1.75, 1.75, 0.);
  int nx = 7, ny = 7, nz = 6;

  Teuchos::RCP<Mesh_MSTK> m_bg =
    Teuchos::rcp(new Mesh_MSTK(X0[0], X0[1],X0[2], X1[0], X1[1], X1[2],
			       nx, ny, nz, comm, m_log->getGeometricModel(), Teuchos::null));

  // make the new connections, 1 per logical cell
  int ncells_log = m_log->getNumEntities(Entity_kind::CELL, Parallel_type::ALL);
  std::vector<Entity_ID_List> face_cell_ids_(ncells_log);

  RegularMeshCellFromCoordFunctor cellID(X0, X1, nx, ny, nz);
  for (int c=0; c!=ncells_log; ++c) {
    int bg_c = cellID(m_log->getCellCentroid(c));
    face_cell_ids_[c].push_back(c);
    face_cell_ids_[c].push_back(bg_c);
  }

  // cross-sectional areas and lengths given by root class
  // relatively made up numbers, this could be done formally if we wanted to.
  Entity_ID_List coarse_roots, fine_roots;
  coarse_roots = m_log->getSetEntities("coarse_root", Entity_kind::CELL, Parallel_type::ALL);
  fine_roots = m_log->getSetEntities("fine_root", Entity_kind::CELL, Parallel_type::ALL);

  std::vector<std::vector<double> > face_cell_lengths(ncells_log);
  std::vector<AmanziGeometry::Point> face_area_normals(ncells_log);
  Point ihat(1.,0.,0.);  // just must have 0 z to avoid gravity
  for (Entity_ID_List::iterator c=coarse_roots.begin();
       c!=coarse_roots.end(); ++c) {
    face_area_normals[*c] = 2.e-2*ihat;
    face_cell_lengths[*c].push_back(1.e-2);
    face_cell_lengths[*c].push_back(0.5/4.);
  }
  for (Entity_ID_List::iterator c=fine_roots.begin();
       c!=fine_roots.end(); ++c) {
    face_area_normals[*c] = 1.e-2*ihat;
    face_cell_lengths[*c].push_back(5.e-3);
    face_cell_lengths[*c].push_back(0.5/4.);
  }

  Teuchos::RCP<MeshEmbeddedLogical> m =
    Teuchos::rcp(new MeshEmbeddedLogical(comm, m_bg, m_log->getMeshFramework(),
  					 face_cell_ids_, face_cell_lengths,
  					 face_area_normals));
  return m;
}



} // namespace Testing
} // namespace Amanzi

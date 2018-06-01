#include <fstream>
#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"

#include "boost/math/constants/constants.hpp"
#include "boost/math/special_functions/round.hpp"

#include "plant_1D_mesh.hh"
#include "MeshLogicalFactory.hh"
#include "MeshLogicalAudit.hh"
#include "Mesh_MSTK.hh"
#include "RegionPlane.hh"
#include "depth_model.hh"


const double PI = boost::math::constants::pi<double>();

TEST(PLANT_CONSTRUCTION) {
  using namespace ATS::Testing;
  using namespace Amanzi;

  Epetra_MpiComm comm(MPI_COMM_SELF);
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3));

  AmanziMesh::MeshLogicalFactory fac(&comm, gm, false);
  double surface_area = 9;  // 3m by 3m single tree?
  
  // topology
  int nstem = 3;
  int nshell = 3;
  int nsoil = 6;

  // geometry
  HydraulicMeshParams p("plant");
  p.height = 3.;
  p.leaf_height = 0.1;
  p.stem_length = p.height - p.leaf_height;
  
  double leaf_density = -2.3231*30 / 2 + 781.899;  // sla = 30
  p.leaf_volume = 5.0 * surface_area / leaf_density; // lai = 5, surface area = 9 m^2
  p.sapwood_area = 0.01/2; // one stem, 10cm square, /2 for area -> sapwood
  p.troot_length = p.height / 4.;  // totally made up
  p.troot_volume = p.troot_length * p.sapwood_area;
  p.aroot_radius = 0.003; // totally made up
  p.aroot_volume = p.troot_volume;
  p.aroot_length = p.troot_volume / (PI * std::pow(p.aroot_radius,2));
  
  std::vector<double> rooting_fraction = { 0.5, 0.3, 0.15, 0.05, 0., 0. };
  std::vector<double> dz = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };
  std::vector<double> rheiz_outer_radius(dz.size(), 0.2);
  
  AMANZI_ASSERT(rooting_fraction.size() == nsoil);
  AMANZI_ASSERT(dz.size() == nsoil);
  AmanziGeometry::Point centroid(0.,0.,0.);
  AmanziGeometry::Point displacement(std::sqrt(surface_area/PI)/2.,0.,0.);
  
  std::vector<std::vector<int>> soil_rheiz_conn_cells;
  std::vector<std::vector<double>> soil_rheiz_conn_lengths;
  std::vector<AmanziGeometry::Point> soil_rheiz_conn_area_normals;
  std::vector<double> belowground_volume;

  plantMesh(fac, p, nstem, nshell,
            rooting_fraction, dz, centroid, displacement,
            rheiz_outer_radius,
            soil_rheiz_conn_cells, soil_rheiz_conn_lengths,
            soil_rheiz_conn_area_normals, belowground_volume);

  auto mesh = fac.Create();
  Teuchos::RCP<AmanziMesh::Mesh> m2(mesh);
  
  AmanziMesh::MeshLogicalAudit audit(m2, std::cout);
  CHECK(!audit.Verify());


  int nroot = 0;
  for (auto r : rooting_fraction) {
    if (r == 0.) break;
    nroot++;
  }
  
  CHECK_EQUAL(1 + nstem + (nshell+2)*nroot,
              m2->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL));
  CHECK_EQUAL(nstem + 2 + (nroot-1) + (nshell+1)*nroot, 
              m2->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL));

  std::ofstream f("plant_mesh.txt");
  CHECK(!viewMeshLogical(*m2,f));
  f.close();

}


TEST(PLANT_CONSTRUCTION_NO_SHELL) {
  using namespace ATS::Testing;
  using namespace Amanzi;

  Epetra_MpiComm comm(MPI_COMM_SELF);
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3));

  AmanziMesh::MeshLogicalFactory fac(&comm, gm, false);
  double surface_area = 9;  // 3m by 3m single tree?
  
  // topology
  int nstem = 3;
  int nshell = 0;
  int nsoil = 6;

  // geometry
  HydraulicMeshParams p("plant");
  p.height = 3.;
  p.leaf_height = 0.1;
  p.stem_length = p.height - p.leaf_height;
  
  double leaf_density = -2.3231*30 / 2 + 781.899;  // sla = 30
  p.leaf_volume = 5.0 * surface_area / leaf_density; // lai = 5, surface area = 9 m^2
  p.sapwood_area = 0.01/2; // one stem, 10cm square, /2 for area -> sapwood
  p.troot_length = p.height / 4.;  // totally made up
  p.troot_volume = p.troot_length * p.sapwood_area;
  p.aroot_radius = 0.003; // totally made up
  p.aroot_volume = p.troot_volume;
  p.aroot_length = p.troot_volume / (PI * std::pow(p.aroot_radius,2));
  
  std::vector<double> rooting_fraction = { 0.5, 0.3, 0.15, 0.05, 0., 0. };
  std::vector<double> dz = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };
  std::vector<double> rheiz_outer_radius;

  AMANZI_ASSERT(rooting_fraction.size() == nsoil);
  AMANZI_ASSERT(dz.size() == nsoil);
  AmanziGeometry::Point centroid(0.,0.,0.);
  AmanziGeometry::Point displacement(std::sqrt(surface_area/PI)/2.,0.,0.);
  
  std::vector<std::vector<int>> soil_rheiz_conn_cells;
  std::vector<std::vector<double>> soil_rheiz_conn_lengths;
  std::vector<AmanziGeometry::Point> soil_rheiz_conn_area_normals;
  std::vector<double> belowground_volume;

  plantMesh(fac, p, nstem, nshell,
            rooting_fraction, dz, centroid, displacement,
            rheiz_outer_radius,
            soil_rheiz_conn_cells, soil_rheiz_conn_lengths,
            soil_rheiz_conn_area_normals, belowground_volume);

  auto mesh = fac.Create();
  Teuchos::RCP<AmanziMesh::Mesh> m2(mesh);
  
  AmanziMesh::MeshLogicalAudit audit(m2, std::cout);
  CHECK(!audit.Verify());


  int nroot = 0;
  for (auto r : rooting_fraction) {
    if (r == 0.) break;
    nroot++;
  }
  
  CHECK_EQUAL(1 + nstem + (nshell+2)*nroot,
              m2->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL));
  CHECK_EQUAL(nstem + 2 + (nroot-1) + (nshell+1)*nroot, 
              m2->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL));

  std::ofstream f("plant_mesh_no_shell.txt");
  CHECK(!viewMeshLogical(*m2,f));
  f.close();

}


TEST(COLUMN_CONSTRUCTION_ONE_PLANT) {
  using namespace ATS::Testing;
  using namespace Amanzi;

  auto comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  // make a geometric model including the surface region
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3));
  auto surface = Teuchos::rcp(new AmanziGeometry::RegionPlane("surface", 1, AmanziGeometry::Point(0.,0.,0.), AmanziGeometry::Point(0.,0.,1.)));
  gm->AddRegion(surface);

  // make a subsurface mesh
  auto subsurface_mesh = Teuchos::rcp(new AmanziMesh::Mesh_MSTK(-1.5, -1.5, -10, 1.5, 1.5, 0, 1, 1, 100, comm, gm));
  subsurface_mesh->build_columns();

  // extract the surface mesh
  auto surface_mesh = Teuchos::rcp(new AmanziMesh::Mesh_MSTK(&*subsurface_mesh, std::vector<std::string>(1,"surface"), AmanziMesh::FACE, true));

  // plant mesh
  // -- rooting fraction
  Epetra_MultiVector root_fraction(subsurface_mesh->cell_map(false),1);
  double sa = surface_mesh->cell_volume(0);
  for (int c=0; c!=subsurface_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED); ++c) {
    auto cc = subsurface_mesh->cell_centroid(c);
    double depth = 0 - cc[2];
    double thickness = subsurface_mesh->cell_volume(c) / sa;
    root_fraction[0][c] = depth > 1 ? 0. : 2.0 * (1 - depth) * thickness;
  }

  // -- pfts
  std::vector<std::vector<BGC::PFT>> pfts(1);
  pfts[0].emplace_back(BGC::PFT("plant", 100));
  Teuchos::ParameterList pft_plist;
  pfts[0][0].Init(pft_plist, sa);

  pfts[0][0].lai = 5;
  pfts[0][0].Bleaf = pfts[0][0].lai * sa / pfts[0][0].SLA;
  pfts[0][0].Bstem = pfts[0][0].Bleaf / pfts[0][0].leaf2stemratio;
  pfts[0][0].Broot = pfts[0][0].Bleaf / pfts[0][0].leaf2rootratio;

  // -- constructor
  auto embedded = plantMeshes(subsurface_mesh, gm, *surface_mesh, root_fraction, pfts, 5, 5, 0.1);

  // verify
  AmanziMesh::MeshLogicalAudit audit(embedded, std::cout);
  CHECK(!audit.Verify());
  
  std::ofstream f("soil_embedded_plant.txt");
  CHECK(!viewMeshLogical(*embedded,f));
  f.close();
  

}


TEST(COLUMN_CONSTRUCTION_TREE_SHRUB) {
  using namespace ATS::Testing;
  using namespace Amanzi;

  auto comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  // make a geometric model including the surface region
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3));
  auto surface = Teuchos::rcp(new AmanziGeometry::RegionPlane("surface", 1, AmanziGeometry::Point(0.,0.,0.), AmanziGeometry::Point(0.,0.,1.)));
  gm->AddRegion(surface);

  // make a subsurface mesh
  auto subsurface_mesh = Teuchos::rcp(new AmanziMesh::Mesh_MSTK(-1.5, -1.5, -10, 1.5, 1.5, 0, 1, 1, 100, comm, gm));
  subsurface_mesh->build_columns();

  // extract the surface mesh
  auto surface_mesh = Teuchos::rcp(new AmanziMesh::Mesh_MSTK(&*subsurface_mesh, std::vector<std::string>(1,"surface"), AmanziMesh::FACE, true));

  // plant mesh
  // -- rooting fraction
  Epetra_MultiVector root_fraction(subsurface_mesh->cell_map(false),2);
  double sa = surface_mesh->cell_volume(0);
  for (int c=0; c!=subsurface_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED); ++c) {
    auto cc = subsurface_mesh->cell_centroid(c);
    double depth = 0 - cc[2];
    double thickness = subsurface_mesh->cell_volume(c) / sa;
    root_fraction[0][c] = depth > 1 ? 0. : 2.0 * (1 - depth) * thickness;
    root_fraction[1][c] = depth > 0.5 ? 0. : 8 * (0.5-depth) * thickness;
  }

  // -- pfts
  std::vector<std::vector<BGC::PFT>> pfts(2);
  pfts[0].emplace_back(BGC::PFT("tree", 100));
  pfts[0].emplace_back(BGC::PFT("shrub", 100));
  Teuchos::ParameterList pft_plist;
  pfts[0][0].Init(pft_plist, sa);
  pfts[0][1].Init(pft_plist, sa);

  pfts[0][0].lai = 5;
  pfts[0][0].Bleaf = pfts[0][0].lai * sa / pfts[0][0].SLA;
  pfts[0][0].Bstem = pfts[0][0].Bleaf / pfts[0][0].leaf2stemratio;
  pfts[0][0].Broot = pfts[0][0].Bleaf / pfts[0][0].leaf2rootratio;

  pfts[0][1].lai = 2;
  pfts[0][1].height = 1;
  pfts[0][1].Bleaf = pfts[0][0].lai * sa / pfts[0][0].SLA;
  pfts[0][1].leaf2stemratio = 10;
  pfts[0][1].leaf2rootratio = 1;
  pfts[0][1].Bstem = pfts[0][0].Bleaf / pfts[0][0].leaf2stemratio;
  pfts[0][1].Broot = pfts[0][0].Bleaf / pfts[0][0].leaf2rootratio;
  
  // -- constructor
  auto embedded = plantMeshes(subsurface_mesh, gm, *surface_mesh, root_fraction, pfts, 5, 5, 0.1);

  // verify
  AmanziMesh::MeshLogicalAudit audit(embedded, std::cout);
  CHECK(!audit.Verify());
  
  std::ofstream f("soil_embedded_tree_shrub.txt");
  CHECK(!viewMeshLogical(*embedded,f));
  f.close();
  

}




TEST(COLUMN_CONSTRUCTION_TREE_SHRUB_HILLSLOPE) {
  using namespace ATS::Testing;
  using namespace Amanzi;

  auto comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  // read input to get regions, mesh
  std::string xmlFileName = "test/open_book.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  
  // make a geometric model including the surface region
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, region_list, comm));

  // make a subsurface mesh
  auto subsurface_mesh = Teuchos::rcp(new AmanziMesh::Mesh_MSTK("test/open_book.exo", comm, gm));
  subsurface_mesh->build_columns("surface");

  // extract the surface mesh
  std::vector<std::string> surface_ss = { "surface" };
  auto surface_mesh = Teuchos::rcp(new AmanziMesh::Mesh_MSTK(&*subsurface_mesh, surface_ss, AmanziMesh::FACE, true));

  // plant mesh
  // -- depth of each cell
  Epetra_Vector depth(subsurface_mesh->cell_map(false));
  Flow::DepthModel(*subsurface_mesh, depth);

  // -- cumulative rooting fraction
  Epetra_MultiVector root_fraction(subsurface_mesh->cell_map(false),2);
  root_fraction.PutScalar(-1.0);
  for (int sc=0; sc!=surface_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED); ++sc) {
    auto scc = subsurface_mesh->face_centroid(surface_mesh->entity_get_parent(AmanziMesh::CELL, sc));
    double sa = surface_mesh->cell_volume(sc);
    
    const auto& faces_of_col = subsurface_mesh->faces_of_column(sc);
    const auto& cells_of_col = subsurface_mesh->cells_of_column(sc);

    for (int i_pft=0; i_pft!=2; ++i_pft) {
      std::vector<double> cum_rooting_frac(faces_of_col.size(), 0.);
      for (int i=0; i!=faces_of_col.size(); ++i) {
        auto fc = subsurface_mesh->face_centroid(faces_of_col[i]);
        double depth = scc[2] - fc[2];
        AMANZI_ASSERT(depth >= 0.);

        cum_rooting_frac[i] = i_pft == 0 ? (depth > 3 ? 1 : depth/3) : (depth > 1 ? 1 : depth);
      }

      for (int cc=0; cc!=faces_of_col.size()-1; ++cc) {
        root_fraction[i_pft][cells_of_col[cc]] = cum_rooting_frac[cc+1] - cum_rooting_frac[cc];
      }
    }
  }

  // -- pfts
  std::vector<std::vector<BGC::PFT>> pfts;
  for (int sc=0; sc!=surface_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED); ++sc) {
    double sa = surface_mesh->cell_volume(sc);

    std::vector<BGC::PFT> cell_pfts = { BGC::PFT("tree", 100), BGC::PFT("shrub", 100) };
    Teuchos::ParameterList pft_plist;
    cell_pfts[0].Init(pft_plist, sa);
    cell_pfts[1].Init(pft_plist, sa);

    cell_pfts[0].lai = 5;
    cell_pfts[0].Bleaf = cell_pfts[0].lai * sa / cell_pfts[0].SLA;
    cell_pfts[0].Bstem = cell_pfts[0].Bleaf / cell_pfts[0].leaf2stemratio;
    cell_pfts[0].Broot = cell_pfts[0].Bleaf / cell_pfts[0].leaf2rootratio;
    
    cell_pfts[1].lai = 2;
    cell_pfts[1].height = 1;
    cell_pfts[1].Bleaf = cell_pfts[0].lai * sa / cell_pfts[0].SLA;
    cell_pfts[1].leaf2stemratio = 10;
    cell_pfts[1].leaf2rootratio = 1;
    cell_pfts[1].Bstem = cell_pfts[0].Bleaf / cell_pfts[0].leaf2stemratio;
    cell_pfts[1].Broot = cell_pfts[0].Bleaf / cell_pfts[0].leaf2rootratio;

    pfts.emplace_back(cell_pfts);
  }
  
  // -- constructor
  auto embedded = plantMeshes(subsurface_mesh, gm, *surface_mesh, root_fraction, pfts, 5, 5, 0.1);

  // verify
  AmanziMesh::MeshLogicalAudit audit(embedded, std::cout);
  CHECK(!audit.Verify());
  
  std::ofstream f("soil_embedded_tree_shrub.txt");
  CHECK(!viewMeshLogical(*embedded,f));
  f.close();
  

}

  

  

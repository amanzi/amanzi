/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
// -------------------------------------------------------------
/**
 * @file   test_column_mesh.cc
 * @author Rao V. Garimella
 * @date   Thu Mar 18, 2015
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------

#include <UnitTest++.h>

#include <mpi.h>
#include <fstream>

#include "AmanziComm.hh"

#include "Geometry.hh"
#include "../Mesh_MSTK.hh"
#include "../MeshColumn.hh"
#include "../MeshSurfaceCell.hh"
#include "RegionBox.hh"
#include "RegionLabeledSet.hh"
#include "GeometricModel.hh"


TEST(SURFACE_COLUMN_MESH_3D)
{
  auto comm = Amanzi::getDefaultComm();

  int nx = 4, ny = 4, nz = 4;
  double lx = 4, ly = 4, lz = 4;

  // create a geometric model with surface region
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3));

  Amanzi::AmanziGeometry::Point p0(3), p1(3);
  p0[2] = 3.999;
  p1[0] = 4.; p1[1] = 4.; p1[2] = 5.;
  auto r0 = Teuchos::rcp(new Amanzi::AmanziGeometry::RegionBox("surface", -1, p0, p1));
  gm->AddRegion(r0);

  Amanzi::AmanziGeometry::Point p2(2), p3(2);
  p3[0] = 4.; p3[1] = 4.;
  auto r1 = Teuchos::rcp(new Amanzi::AmanziGeometry::RegionBox("surface_domain", -1, p2, p3));
  gm->AddRegion(r1);
  
  // Create the mesh
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh =
      Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(0.0,0.0,0.0,
              lx,ly,lz,nx,ny,nz,
              comm, gm));

  CHECK_EQUAL(1, mesh->build_columns());
  
  // Perturb the nodes above the base layer just a bit
  int nnodes = mesh->num_entities(Amanzi::AmanziMesh::NODE,
          Amanzi::AmanziMesh::Parallel_type::OWNED);
    
  for (int n = 0; n < nnodes; n++) {
    Amanzi::AmanziGeometry::Point xyz(3);
    mesh->node_get_coordinates(n,&xyz);
    xyz[2] += 0.005*xyz[0]*xyz[1]*xyz[2];
    mesh->node_set_coordinates(n,xyz);
  }

  // Create a column mesh from one of the columns
  auto colmesh = Teuchos::rcp(new Amanzi::AmanziMesh::MeshColumn(mesh,10));

  // Extract the surface from this column
  Amanzi::AmanziMesh::MeshSurfaceCell col_surf(colmesh, "surface");

  // -- check basic mesh structure
  CHECK_EQUAL(1, col_surf.num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED));
  CHECK_EQUAL(4, col_surf.num_entities(Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::OWNED));
  CHECK_EQUAL(4, col_surf.num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::OWNED));

  // -- check flattened
  Amanzi::AmanziGeometry::Point node;
  col_surf.node_get_coordinates(0, &node);
  CHECK_EQUAL(2, node.dim());

  // -- check sets
  Amanzi::AmanziMesh::Entity_ID_List cells_in_surf;
  col_surf.get_set_entities_and_vofs("surface", Amanzi::AmanziMesh::CELL,
                                     Amanzi::AmanziMesh::Parallel_type::OWNED,
                                     &cells_in_surf, NULL);
  CHECK_EQUAL(1, cells_in_surf.size());
  CHECK_EQUAL(0, cells_in_surf[0]);

  Amanzi::AmanziMesh::Entity_ID_List cells_in_surf_2D;
  col_surf.get_set_entities_and_vofs("surface_domain", Amanzi::AmanziMesh::CELL,
                                     Amanzi::AmanziMesh::Parallel_type::OWNED,
                                     &cells_in_surf_2D, NULL);
  CHECK_EQUAL(1, cells_in_surf.size());
  CHECK_EQUAL(0, cells_in_surf[0]);

  // -- check volumes  
  CHECK_CLOSE(1.0, col_surf.cell_volume(0), 1.e-9);
  CHECK_CLOSE(1.0, col_surf.face_area(3), 1.e-9);
  
}

TEST(SURFACE_COLUMN_MESH_3D_UNSTRUCTURED)
{
  auto comm = Amanzi::getDefaultComm();

  // create a geometric model with surface region
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3));

  Teuchos::RCP<Amanzi::AmanziGeometry::RegionLabeledSet> r0 =
      Teuchos::rcp(new Amanzi::AmanziGeometry::RegionLabeledSet("surface", -1, "FACE",
              "test/slab-0.05-5x4x25.exo", "Exodus II", "1"));
  gm->AddRegion(r0);

  Amanzi::AmanziGeometry::Point p2(2), p3(2);
  p3[0] = 4.; p3[1] = 4.;
  auto r1 = Teuchos::rcp(new Amanzi::AmanziGeometry::RegionBox("surface_domain", -1, p2, p3));
  gm->AddRegion(r1);

  // Create the mesh
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh =
      Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK("test/slab-0.05-5x4x25.exo", comm, gm));

  CHECK_EQUAL(20, mesh->get_set_size("surface",Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::ALL));

  // Build columns in the mesh
  CHECK_EQUAL(1, mesh->build_columns());
  
  // Create a column mesh from one of the columns
  auto colmesh = Teuchos::rcp(new Amanzi::AmanziMesh::MeshColumn(mesh,10));
  CHECK_EQUAL(1, colmesh->get_set_size("surface",Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::ALL));

  // Extract the surface from this column
  Amanzi::AmanziMesh::MeshSurfaceCell col_surf(colmesh, "surface");

  CHECK_EQUAL(1, col_surf.num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED));
  CHECK_EQUAL(4, col_surf.num_entities(Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::OWNED));
  CHECK_EQUAL(4, col_surf.num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::OWNED));

  Amanzi::AmanziMesh::Entity_ID_List cells_in_surf;
  col_surf.get_set_entities_and_vofs("surface", Amanzi::AmanziMesh::CELL,
                                     Amanzi::AmanziMesh::Parallel_type::OWNED,
                                     &cells_in_surf, NULL);
  CHECK_EQUAL(1, cells_in_surf.size());
  CHECK_EQUAL(0, cells_in_surf[0]);

  Amanzi::AmanziMesh::Entity_ID_List cells_in_surf_2D;
  col_surf.get_set_entities_and_vofs("surface_domain", Amanzi::AmanziMesh::CELL,
                                     Amanzi::AmanziMesh::Parallel_type::OWNED,
                                     &cells_in_surf_2D, NULL);
  CHECK_EQUAL(1, cells_in_surf.size());
  CHECK_EQUAL(0, cells_in_surf[0]);
  
  CHECK_CLOSE(6400.0, col_surf.cell_volume(0), 1.e-9);
  CHECK_CLOSE(80.0, col_surf.face_area(3), 1.e-9);
  
}


TEST(SURFACE_COLUMN_MESH_3D_UNSTRUCTURED_SETS)
{
  std::string filename("test/hex_3x3x3_sets.exo");
  auto comm = Amanzi::getDefaultComm();
  Teuchos::ParameterList parameterlist;
 
  // create a sublist name Regions and put a reference to it in
  // reg_spec and other sublists as references. Turns out it is
  // important to define reg_spec and other lists below as references
  // - otherwise, a new copy is made of the sublist that is retrieved
  Teuchos::ParameterList& reg_spec = parameterlist.sublist("regions"); 
  
  Teuchos::ParameterList& top_surface = reg_spec.sublist("Top Surface");
  Teuchos::ParameterList& top_surface_def = top_surface.sublist("region: labeled set");
  top_surface_def.set<std::string>("label","106");
  top_surface_def.set<std::string>("file",filename.c_str());
  top_surface_def.set<std::string>("format","Exodus II");
  top_surface_def.set<std::string>("entity","face");
  
  Teuchos::ParameterList& side_surface = reg_spec.sublist("Side Surface");
  Teuchos::ParameterList& side_surface_def = side_surface.sublist("region: labeled set");
  side_surface_def.set<std::string>("label","102");
  side_surface_def.set<std::string>("file",filename.c_str());
  side_surface_def.set<std::string>("format","Exodus II");
  side_surface_def.set<std::string>("entity","face");
  
  Teuchos::ParameterList& r1_surface = reg_spec.sublist("Region 1");
  Teuchos::ParameterList& r1_surface_def = r1_surface.sublist("region: labeled set");
  r1_surface_def.set<std::string>("label","30000");
  r1_surface_def.set<std::string>("file",filename.c_str());
  r1_surface_def.set<std::string>("format","Exodus II");
  r1_surface_def.set<std::string>("entity","cell");
  
  Teuchos::ParameterList& r2_surface = reg_spec.sublist("Region 2");
  Teuchos::ParameterList& r2_surface_def = r2_surface.sublist("region: labeled set");
  r2_surface_def.set<std::string>("label","20000");
  r2_surface_def.set<std::string>("file",filename.c_str());
  r2_surface_def.set<std::string>("format","Exodus II");
  r2_surface_def.set<std::string>("entity","cell");

  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));
  auto mesh = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(filename.c_str(), comm, gm));

  CHECK_EQUAL(9, mesh->get_set_size("Top Surface",Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::ALL));
  CHECK_EQUAL(1, mesh->build_columns());
  
  // Create a column mesh from one of the columns
  auto colmesh = Teuchos::rcp(new Amanzi::AmanziMesh::MeshColumn(mesh,0));
  CHECK_EQUAL(1, colmesh->get_set_size("Top Surface",Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::ALL));
  CHECK_EQUAL(1, colmesh->get_set_size("Region 1",Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::ALL));

  // Extract the surface from this column
  Amanzi::AmanziMesh::MeshSurfaceCell col_surf(colmesh, "Top Surface");

  CHECK_EQUAL(1, col_surf.num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED));
  CHECK_EQUAL(4, col_surf.num_entities(Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::OWNED));
  CHECK_EQUAL(4, col_surf.num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::OWNED));

  CHECK(col_surf.valid_set_name("Top Surface", Amanzi::AmanziMesh::CELL));
  CHECK(col_surf.valid_set_name("Region 1", Amanzi::AmanziMesh::CELL));
  CHECK(col_surf.valid_set_name("Top Surface", Amanzi::AmanziMesh::FACE));
  CHECK(!col_surf.valid_set_name("Region 1", Amanzi::AmanziMesh::FACE));
  
  Amanzi::AmanziMesh::Entity_ID_List cells_in_surf;
  col_surf.get_set_entities_and_vofs("Top Surface", Amanzi::AmanziMesh::CELL,
                                     Amanzi::AmanziMesh::Parallel_type::OWNED,
                                     &cells_in_surf, NULL);
  CHECK_EQUAL(1, cells_in_surf.size());
  CHECK_EQUAL(0, cells_in_surf[0]);

  Amanzi::AmanziMesh::Entity_ID_List cells_in_surf2;
  col_surf.get_set_entities_and_vofs("Region 1", Amanzi::AmanziMesh::CELL,
                                     Amanzi::AmanziMesh::Parallel_type::OWNED,
                                     &cells_in_surf2, NULL);
  CHECK_EQUAL(1, cells_in_surf.size());
  CHECK_EQUAL(0, cells_in_surf[0]);
}

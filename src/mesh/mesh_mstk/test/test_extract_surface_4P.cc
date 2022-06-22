/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <UnitTest++.h>
#include <fstream>

#include "../Mesh_MSTK.hh"


#include "AmanziMap.hh"
#include "AmanziComm.hh"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Array.hpp"

// Extract some surfaces as-is from 3D mesh

TEST(Extract_Surface_MSTK1_4P)
{
  auto comm = Amanzi::getDefaultComm();


  Teuchos::ParameterList parameterlist;


  // create a sublist name Regions and put a reference to it in
  // reg_spec and other sublists as references. Turns out it is
  // important to define reg_spec and other lists below as references
  // - otherwise, a new copy is made of the sublist that is retrieved

  Teuchos::ParameterList& reg_spec = parameterlist.sublist("regions");

  Teuchos::ParameterList& top_surface = reg_spec.sublist("Top Surface");
  Teuchos::ParameterList& top_surface_def =
    top_surface.sublist("region: plane");
  Teuchos::Array<double> loc1 = Teuchos::tuple(0.0, 0.0, 1.0);
  Teuchos::Array<double> dir1 = Teuchos::tuple(0.0, 0.0, -1.0);
  top_surface_def.set<Teuchos::Array<double>>("point", loc1);
  top_surface_def.set<Teuchos::Array<double>>("normal", dir1);

  Teuchos::ParameterList& right_surface = reg_spec.sublist("Right Surface");
  Teuchos::ParameterList& right_surface_def =
    right_surface.sublist("region: plane");
  Teuchos::Array<double> loc2 = Teuchos::tuple(1.0, 0.0, 0.0);
  Teuchos::Array<double> dir2 = Teuchos::tuple(1.0, 0.0, 0.0);
  right_surface_def.set<Teuchos::Array<double>>("point", loc2);
  right_surface_def.set<Teuchos::Array<double>>("normal", dir2);

  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm = Teuchos::rcp(
    new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // Generate a mesh consisting of 4x4x4 elements
  auto mesh = Teuchos::rcp(
    new Amanzi::AmanziMesh::Mesh_MSTK(0, 0, 0, 1, 1, 1, 4, 4, 4, comm, gm));

  std::vector<std::string> setnames;
  setnames.push_back(std::string("Top Surface"));
  setnames.push_back(std::string("Right Surface"));

  Amanzi::AmanziMesh::Entity_ID_List ids1, ids2;
  mesh->get_set_entities(setnames[0],
                         Amanzi::AmanziMesh::FACE,
                         Amanzi::AmanziMesh::Parallel_type::OWNED,
                         ids1);
  mesh->get_set_entities(setnames[1],
                         Amanzi::AmanziMesh::FACE,
                         Amanzi::AmanziMesh::Parallel_type::OWNED,
                         ids2);

  ids1 = ids2;

  Amanzi::AmanziMesh::Mesh_MSTK surfmesh(
    mesh, ids1, Amanzi::AmanziMesh::FACE, false, mesh->get_comm());

  // Number of cells (quadrilaterals) in surface mesh
  int ncells_surf = surfmesh.num_entities(
    Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::ALL);

  // Check if centroid of the surface mesh cell is the same as its
  // parent (face) in the volume mesh
  for (int k = 0; k < ncells_surf; k++) {
    Amanzi::AmanziMesh::Entity_ID parent =
      surfmesh.entity_get_parent_host(Amanzi::AmanziMesh::CELL, k);

    Amanzi::AmanziGeometry::Point centroid1 = mesh->face_centroid_host(parent);
    Amanzi::AmanziGeometry::Point centroid2 = surfmesh.cell_centroid_host(k);
    CHECK_ARRAY_EQUAL(centroid1, centroid2, 3);
  }

  // Number of nodes in surface mesh
  int nnodes_surf = surfmesh.num_entities(
    Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::ALL);
  int nnodes_surf_owned = surfmesh.num_entities(
    Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::OWNED);

  // Check if coordinates of surface mesh node is the same as its
  // parent node in the volume mesh
  for (int k = 0; k < nnodes_surf; k++) {
    Amanzi::AmanziMesh::Entity_ID parent =
      surfmesh.entity_get_parent_host(Amanzi::AmanziMesh::NODE, k);
    Amanzi::AmanziGeometry::Point coord1;
    mesh->node_get_coordinates(parent, &coord1);
    Amanzi::AmanziGeometry::Point coord2;
    surfmesh.node_get_coordinates(k, &coord2);
    CHECK_ARRAY_EQUAL(coord1, coord2, 3);
  }
}


// Extract a surface of a 3D mesh and flatten it to 2D to make new mesh

TEST(Extract_Surface_MSTK2_4P)
{
  auto comm = Amanzi::getDefaultComm();


  Teuchos::ParameterList parameterlist;


  // create a sublist name Regions and put a reference to it in
  // reg_spec and other sublists as references. Turns out it is
  // important to define reg_spec and other lists below as references
  // - otherwise, a new copy is made of the sublist that is retrieved

  Teuchos::ParameterList& reg_spec = parameterlist.sublist("regions");

  Teuchos::ParameterList& top_surface = reg_spec.sublist("Top Surface");
  Teuchos::ParameterList& top_surface_def =
    top_surface.sublist("region: plane");
  Teuchos::Array<double> loc1 = Teuchos::tuple(0.0, 0.0, 1.0);
  Teuchos::Array<double> dir1 =
    Teuchos::tuple(-1 / sqrt(2.0), 0.0, 1 / sqrt(2.0));
  top_surface_def.set<Teuchos::Array<double>>("point", loc1);
  top_surface_def.set<Teuchos::Array<double>>("normal", dir1);

  //  Teuchos::writeParameterListToXmlOStream(parameterlist,std::cout);
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm = Teuchos::rcp(
    new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // Generate a mesh consisting of 4x4x4 elements
  auto mesh = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(
    0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4, comm, gm));


  // Perturb some nodes

  int nv = mesh->num_entities(Amanzi::AmanziMesh::NODE,
                              Amanzi::AmanziMesh::Parallel_type::OWNED);

  for (int i = 0; i < nv; i++) {
    Amanzi::AmanziGeometry::Point pt;
    mesh->node_get_coordinates(i, &pt);
    if (pt[2] == 1.0) {
      pt[2] = pt[2] + pt[0];
      mesh->node_set_coordinates(i, pt);
    }
  }


  std::vector<std::string> setnames;
  setnames.push_back(std::string("Top Surface"));

  Amanzi::AmanziMesh::Entity_ID_List ids1;
  mesh->get_set_entities(setnames[0],
                         Amanzi::AmanziMesh::FACE,
                         Amanzi::AmanziMesh::Parallel_type::OWNED,
                         ids1);

  // Extract surface mesh while projecting to 2D
  Amanzi::AmanziMesh::Mesh_MSTK surfmesh(
    mesh, ids1, Amanzi::AmanziMesh::FACE, true, mesh->get_comm());

  CHECK_EQUAL(surfmesh.space_dimension(), 2);

  // Number of cells (quadrilaterals) in surface mesh

  int ncells_surf = surfmesh.num_entities(
    Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::ALL);

  // Check if centroid of a surface mesh cell is the same as its
  // parent (face) in the volume mesh

  for (int k = 0; k < ncells_surf; k++) {
    Amanzi::AmanziMesh::Entity_ID parent =
      surfmesh.entity_get_parent_host(Amanzi::AmanziMesh::CELL, k);

    Amanzi::AmanziGeometry::Point centroid1 = mesh->face_centroid_host(parent);
    Amanzi::AmanziGeometry::Point centroid2 = surfmesh.cell_centroid_host(k);
    CHECK_EQUAL(2, centroid2.dim());
    CHECK_CLOSE(centroid1[0], centroid1[0], 1.0e-10);
    CHECK_CLOSE(centroid1[1], centroid1[1], 1.0e-10);
  }


  // Number of nodes in surface mesh

  int nnodes_surf = surfmesh.num_entities(
    Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::ALL);

  // Check if coordinate of a node in a surface mesh cell is the same
  // as its parent in the volume mesh

  for (int k = 0; k < nnodes_surf; k++) {
    //    Amanzi::AmanziMesh::Amanzi::AmanziMesh::Entity_ID_List nodecells;
    //    surfmesh.node_get_cells(k,Amanzi::AmanziMesh::Parallel_type::OWNED,&nodecells);
    //    if (nodecells.size() == 0) continue;

    Amanzi::AmanziMesh::Entity_ID parent =
      surfmesh.entity_get_parent_host(Amanzi::AmanziMesh::NODE, k);
    Amanzi::AmanziGeometry::Point coord1;
    mesh->node_get_coordinates(parent, &coord1);
    Amanzi::AmanziGeometry::Point coord2;
    surfmesh.node_get_coordinates(k, &coord2);
    CHECK_CLOSE(coord1[0], coord2[0], 1.0e-10);
    CHECK_CLOSE(coord1[1], coord2[1], 1.0e-10);
  }
}


// Extract a surface defined by a labeled side set

TEST(Extract_Surface_MSTK3_4P)
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
  Teuchos::ParameterList& top_surface_def =
    top_surface.sublist("region: labeled set");
  top_surface_def.set<std::string>("label", "106");
  top_surface_def.set<std::string>("file", filename.c_str());
  top_surface_def.set<std::string>("format", "Exodus II");
  top_surface_def.set<std::string>("entity", "face");


  //  Teuchos::writeParameterListToXmlOStream(parameterlist,std::cout);

  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm = Teuchos::rcp(
    new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // Read a mesh from the file
  auto mesh =
    Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(filename.c_str(), comm, gm));


  std::vector<std::string> setnames;
  setnames.push_back(std::string("Top Surface"));
  Amanzi::AmanziMesh::Entity_ID_List ids1;
  mesh->get_set_entities(setnames[0],
                         Amanzi::AmanziMesh::FACE,
                         Amanzi::AmanziMesh::Parallel_type::OWNED,
                         ids1);

  Amanzi::AmanziMesh::Mesh_MSTK surfmesh(
    mesh, ids1, Amanzi::AmanziMesh::FACE, false, mesh->get_comm());


  // Number of cells (quadrilaterals) in surface mesh

  int ncells_surf = surfmesh.num_entities(
    Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::ALL);

  // Check if centroid of the surface mesh cell is the same as its
  // parent (face) in the volume mesh

  for (int k = 0; k < ncells_surf; k++) {
    Amanzi::AmanziMesh::Entity_ID parent =
      surfmesh.entity_get_parent_host(Amanzi::AmanziMesh::CELL, k);

    Amanzi::AmanziGeometry::Point centroid1 = mesh->face_centroid_host(parent);
    Amanzi::AmanziGeometry::Point centroid2 = surfmesh.cell_centroid_host(k);
    CHECK_ARRAY_EQUAL(centroid1, centroid2, 3);
  }

  // Number of nodes in surface mesh

  int nnodes_surf = surfmesh.num_entities(
    Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::ALL);

  // Check if coordinate of a node in a surface mesh cell is the same
  // as its parent in the volume mesh

  for (int k = 0; k < nnodes_surf; k++) {
    //    Amanzi::AmanziMesh::Entity_ID_List nodecells;
    //    surfmesh.node_get_cells(k,Amanzi::AmanziMesh::Parallel_type::OWNED,&nodecells);
    //    if (nodecells.size() == 0) continue;

    Amanzi::AmanziMesh::Entity_ID parent =
      surfmesh.entity_get_parent_host(Amanzi::AmanziMesh::NODE, k);
    Amanzi::AmanziGeometry::Point coord1;
    mesh->node_get_coordinates(parent, &coord1);
    Amanzi::AmanziGeometry::Point coord2;
    surfmesh.node_get_coordinates(k, &coord2);
    CHECK_CLOSE(coord1[0], coord2[0], 1.0e-10);
    CHECK_CLOSE(coord1[1], coord2[1], 1.0e-10);
  }
}

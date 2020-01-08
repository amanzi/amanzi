#include <UnitTest++.h>
#include <fstream>

#include "../Mesh_MSTK.hh"


#include "Epetra_Map.h"
#include "AmanziComm.hh"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Array.hpp"

#include "mpi.h"


// Extract some surfaces as-is from 3D mesh

TEST(Extract_Surface_MSTK1)
{

  auto comm = Amanzi::getDefaultComm();


  Teuchos::ParameterList parameterlist;

 
  // create a sublist name Regions and put a reference to it in
  // reg_spec and other sublists as references. Turns out it is
  // important to define reg_spec and other lists below as references
  // - otherwise, a new copy is made of the sublist that is retrieved

  Teuchos::ParameterList& reg_spec = parameterlist.sublist("regions"); 
  
  Teuchos::ParameterList& top_surface = reg_spec.sublist("Top Surface");
  Teuchos::ParameterList& top_surface_def = top_surface.sublist("region: plane");
  Teuchos::Array<double> loc1 = Teuchos::tuple(0.0,0.0,1.0);
  Teuchos::Array<double> dir1 = Teuchos::tuple(0.0,0.0,-1.0);
  top_surface_def.set< Teuchos::Array<double> >("point",loc1);
  top_surface_def.set< Teuchos::Array<double> >("normal",dir1);

  Teuchos::ParameterList& right_surface = reg_spec.sublist("Right Surface");
  Teuchos::ParameterList& right_surface_def = right_surface.sublist("region: plane");
  Teuchos::Array<double> loc2 = Teuchos::tuple(1.0,0.0,0.0);
  Teuchos::Array<double> dir2 = Teuchos::tuple(1.0,0.0,0.0);
  right_surface_def.set< Teuchos::Array<double> >("point",loc2);
  right_surface_def.set< Teuchos::Array<double> >("normal",dir2);


  //  Teuchos::writeParameterListToXmlOStream(parameterlist,std::cout);
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // Generate a mesh consisting of 3x3x3 elements 
  auto mesh = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(0,0,0,1,1,1,3,3,3,comm,gm));

  std::vector<std::string> setnames;
  setnames.push_back(std::string("Top Surface"));
  setnames.push_back(std::string("Right Surface"));

  Amanzi::AmanziMesh::Entity_ID_List ids1, ids2;
  mesh->get_set_entities(setnames[0], Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::OWNED, &ids1);
  mesh->get_set_entities(setnames[1], Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::OWNED, &ids2);
  ids1.insert(ids1.end(), ids2.begin(), ids2.end());
  Amanzi::AmanziMesh::Mesh_MSTK surfmesh(mesh,ids1,Amanzi::AmanziMesh::FACE, false, mesh->get_comm());


  // Number of cells (quadrilaterals) in surface mesh
  int ncells_surf = surfmesh.get_set_size(setnames[0],Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::Parallel_type::OWNED);
  ncells_surf += surfmesh.get_set_size(setnames[1],Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::Parallel_type::OWNED);
  //  int ncells_surf = surfmesh.num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(18,ncells_surf);
      
  // Number of "faces" (edges) in surface mesh

  int nfaces_surf = surfmesh.num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(45,nfaces_surf);

  // Number of nodes in surface mesh

  int nnodes_surf = surfmesh.num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(28,nnodes_surf);


  int exp_parent_faces[18] = {27,28,29,30,31,32,33,34,35,99,100,101,102,103,104,105,106,107};

  int *found = new int[ncells_surf];
  for (int k = 0; k < ncells_surf; k++) {
    Amanzi::AmanziMesh::Entity_ID parent = surfmesh.entity_get_parent(Amanzi::AmanziMesh::CELL,k);
        
    found[k] = 0;
    for (int kk = 0; kk < ncells_surf; kk++) {
      if (exp_parent_faces[kk] == parent) {        
        Amanzi::AmanziGeometry::Point centroid1 = mesh->face_centroid(parent);
        Amanzi::AmanziGeometry::Point centroid2 = surfmesh.cell_centroid(k);
        CHECK_ARRAY_EQUAL(centroid1,centroid2,3);
        found[k]++;
      }
    }
  }
  
  for (int k = 0; k < ncells_surf; k++)
    CHECK_EQUAL(1,found[k]);
  delete [] found;
}


// Extract a surface of a 3D mesh and flatten it to 2D to make new mesh

TEST(Extract_Surface_MSTK2)
{

  auto comm = Amanzi::getDefaultComm();


  Teuchos::ParameterList parameterlist;

 
  // create a sublist name Regions and put a reference to it in
  // reg_spec and other sublists as references. Turns out it is
  // important to define reg_spec and other lists below as references
  // - otherwise, a new copy is made of the sublist that is retrieved

  Teuchos::ParameterList& reg_spec = parameterlist.sublist("regions"); 

  Teuchos::ParameterList& top_surface = reg_spec.sublist("Top Surface");
  Teuchos::ParameterList& top_surface_def = top_surface.sublist("region: plane");
  Teuchos::Array<double> loc1 = Teuchos::tuple(0.0,0.0,1.0);
  Teuchos::Array<double> dir1 = Teuchos::tuple(-1/sqrt(2.0),0.0,1/sqrt(2.0));
  top_surface_def.set< Teuchos::Array<double> >("point",loc1);
  top_surface_def.set< Teuchos::Array<double> >("normal",dir1);

  

  //  Teuchos::writeParameterListToXmlOStream(parameterlist,std::cout);
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // Generate a mesh consisting of 3x3x3 elements 
  auto mesh = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(0.0,0.0,0.0,1.0,1.0,1.0,4,4,4,comm,gm));


  // Perturb some nodes

  int nv = mesh->num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::Parallel_type::OWNED);

  for (int i = 0; i < nv; i++) {
    Amanzi::AmanziGeometry::Point pt;
    mesh->node_get_coordinates(i,&pt);
    if (pt[2] == 1.0) {
      pt[2] = pt[2]+pt[0];
      mesh->node_set_coordinates(i,pt);
    }
  }
  

  std::vector<std::string> setnames;
  setnames.push_back(std::string("Top Surface"));

  Amanzi::AmanziMesh::Entity_ID_List ids1;
  mesh->get_set_entities(setnames[0], Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::OWNED, &ids1);

  // Extract surface mesh while projecting to 2D
  Amanzi::AmanziMesh::Mesh_MSTK surfmesh(mesh,ids1,Amanzi::AmanziMesh::FACE,true,mesh->get_comm());

  CHECK_EQUAL(surfmesh.space_dimension(),2);

  // Number of cells (quadrilaterals) in surface mesh

  int ncells_surf = surfmesh.num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(16,ncells_surf);
      
  // Number of "faces" (edges) in surface mesh

  int nfaces_surf = surfmesh.num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(40,nfaces_surf);

  // Number of nodes in surface mesh

  int nnodes_surf = surfmesh.num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(25,nnodes_surf);


  int exp_parent_faces[16] = {224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239};

  int *found = new int[ncells_surf];
  for (int k = 0; k < ncells_surf; k++) {
    Amanzi::AmanziMesh::Entity_ID parent = surfmesh.entity_get_parent(Amanzi::AmanziMesh::CELL,k);
        
    found[k] = 0;
    for (int kk = 0; kk < ncells_surf; kk++) {
      if (exp_parent_faces[kk] == parent) {        
        Amanzi::AmanziGeometry::Point centroid1 = mesh->face_centroid(parent);
        Amanzi::AmanziGeometry::Point centroid2 = surfmesh.cell_centroid(k);
        CHECK_EQUAL(2,centroid2.dim());
        CHECK_CLOSE(centroid1[0],centroid1[0],1.0e-10);
        CHECK_CLOSE(centroid1[1],centroid1[1],1.0e-10);
        found[k]++;
      }
    }
  }
  
  for (int k = 0; k < ncells_surf; k++)
    CHECK_EQUAL(1,found[k]);
  delete [] found;
}


// Extract a surface defined by a labeled side set

TEST(Extract_Surface_MSTK3)
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
  r1_surface_def.set<std::string>("label","10000");
  r1_surface_def.set<std::string>("file",filename.c_str());
  r1_surface_def.set<std::string>("format","Exodus II");
  r1_surface_def.set<std::string>("entity","cell");

  Teuchos::ParameterList& r2_surface = reg_spec.sublist("Region 2");
  Teuchos::ParameterList& r2_surface_def = r2_surface.sublist("region: labeled set");
  r2_surface_def.set<std::string>("label","20000");
  r2_surface_def.set<std::string>("file",filename.c_str());
  r2_surface_def.set<std::string>("format","Exodus II");
  r2_surface_def.set<std::string>("entity","cell");
  
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // Read a mesh from the file
  auto mesh = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(filename.c_str(),comm,gm));

  // get the surface set
  std::vector<std::string> setnames;
  setnames.push_back(std::string("Top Surface"));
  Amanzi::AmanziMesh::Entity_ID_List ids1;
  mesh->get_set_entities(setnames[0], Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::OWNED, &ids1);

  // extract a surface mesh from the 3D mesh using the set
  Amanzi::AmanziMesh::Mesh_MSTK surfmesh(mesh,ids1,Amanzi::AmanziMesh::FACE, false, mesh->get_comm());

  // Number of cells (quadrilaterals) in surface mesh
  int ncells_surf = surfmesh.num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(9,ncells_surf);
      
  // Number of "faces" (edges) in surface mesh
  int nfaces_surf = surfmesh.num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(24,nfaces_surf);

  // Number of nodes in surface mesh
  int nnodes_surf = surfmesh.num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(16,nnodes_surf);

  // parent ids are the GIDs of the parent mesh's corresponding entities
  int exp_parent_faces[9] = {79,83,91,87,94,101,97,104,107};

  int *found = new int[ncells_surf];
  for (int k = 0; k < ncells_surf; k++) {
    Amanzi::AmanziMesh::Entity_ID parent = surfmesh.entity_get_parent(Amanzi::AmanziMesh::CELL,k);
        
    found[k] = 0;
    for (int kk = 0; kk < ncells_surf; kk++) {
      if (exp_parent_faces[kk] == parent) {        
        Amanzi::AmanziGeometry::Point centroid1 = mesh->face_centroid(parent);
        Amanzi::AmanziGeometry::Point centroid2 = surfmesh.cell_centroid(k);
        CHECK_ARRAY_EQUAL(centroid1,centroid2,3);
        found[k]++;
      }
    }
  }
  
  for (int k = 0; k < ncells_surf; k++)
    CHECK_EQUAL(1,found[k]);
  delete [] found;

  // Test if the labeled set was inherited correctly and if we get the
  // right entities for this set
  //
  // In this case, a face set in the parent becomes a face set in the surface
  Amanzi::AmanziMesh::Entity_ID_List setents;
  surfmesh.get_set_entities("Side Surface",Amanzi::AmanziMesh::FACE,
                            Amanzi::AmanziMesh::Parallel_type::OWNED,&setents);
  CHECK(setents.size() == 3);

  // In this case, a face set in the parent becomes a cell set in the surface
  setents.clear();
  surfmesh.get_set_entities("Top Surface",Amanzi::AmanziMesh::CELL,
                            Amanzi::AmanziMesh::Parallel_type::OWNED,&setents);
  CHECK(setents.size() == 9);

  // In this case, a cell set in the parent becomes a cell set in the surface
  setents.clear();
  surfmesh.get_set_entities("Region 1",Amanzi::AmanziMesh::CELL,
                            Amanzi::AmanziMesh::Parallel_type::OWNED,&setents);
  CHECK(setents.size() == 9);

  // In this case, a cell set in the parent does not exist in the parent
  setents.clear();
  surfmesh.get_set_entities("Region 2",Amanzi::AmanziMesh::CELL,
                            Amanzi::AmanziMesh::Parallel_type::OWNED,&setents);
  CHECK(setents.size() == 0);

  // make sure bad combos dont seg fault?
  setents.clear();
  surfmesh.get_set_entities("Top Surface", Amanzi::AmanziMesh::FACE,
                            Amanzi::AmanziMesh::Parallel_type::OWNED,&setents);
  CHECK(setents.size() == 0);

  setents.clear();
  surfmesh.get_set_entities("Side Surface", Amanzi::AmanziMesh::CELL,
                            Amanzi::AmanziMesh::Parallel_type::OWNED,&setents);
  CHECK(setents.size() == 0);
}




// Extract cells of a surface mesh and flatten to 2D to make new mesh
// Also tests if we can use the same region definition in multiple meshes

TEST(Extract_Surface_MSTK4)
{

  auto comm = Amanzi::getDefaultComm();


  Teuchos::ParameterList parameterlist;

 
  // create a sublist name Regions and put a reference to it in
  // reg_spec and other sublists as references. Turns out it is
  // important to define reg_spec and other lists below as references
  // - otherwise, a new copy is made of the sublist that is retrieved

  Teuchos::ParameterList& reg_spec = parameterlist.sublist("regions"); 

  Teuchos::ParameterList& top_surface = reg_spec.sublist("Top Surface");
  Teuchos::ParameterList& top_surface_def = top_surface.sublist("region: plane");
  Teuchos::Array<double> loc1 = Teuchos::tuple(0.0,0.0,1.0);
  Teuchos::Array<double> dir1 = Teuchos::tuple(-1/sqrt(2.0),0.0,1/sqrt(2.0));
  top_surface_def.set< Teuchos::Array<double> >("point",loc1);
  top_surface_def.set< Teuchos::Array<double> >("normal",dir1);


  Teuchos::ParameterList& whole_surface = reg_spec.sublist("Whole Surface");
  Teuchos::ParameterList& whole_surface_def = whole_surface.sublist("region: box");
  Teuchos::Array<double> corner1 = Teuchos::tuple(0.0,0.0,1.0);
  Teuchos::Array<double> corner2 = Teuchos::tuple(1.0,1.0,3.0);
  whole_surface_def.set< Teuchos::Array<double> >("low coordinate",corner1);
  whole_surface_def.set< Teuchos::Array<double> >("high coordinate",corner2);

  Teuchos::ParameterList& top_surface2D = reg_spec.sublist("Top Surface 2D");
  Teuchos::ParameterList& top_surface2D_def = top_surface2D.sublist("region: plane");
  Teuchos::Array<double> loc2 = Teuchos::tuple(0.0,1.0);
  Teuchos::Array<double> dir2 = Teuchos::tuple(0.0,1.0);
  top_surface2D_def.set< Teuchos::Array<double> >("point",loc2);
  top_surface2D_def.set< Teuchos::Array<double> >("normal",dir2);

  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // Generate a mesh consisting of 3x3x3 elements 
  auto mesh = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(0.0,0.0,0.0,1.0,1.0,1.0,4,4,4,comm,gm));


  // Perturb some nodes
  int nv = mesh->num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  for (int i = 0; i < nv; i++) {
    Amanzi::AmanziGeometry::Point pt;
    mesh->node_get_coordinates(i,&pt);
    if (pt[2] == 1.0) {
      pt[2] = pt[2]+pt[0];
      mesh->node_set_coordinates(i,pt);
    }
  }

  // Extract top surface mesh 
  std::vector<std::string> setnames;
  setnames.push_back(std::string("Top Surface"));
  Amanzi::AmanziMesh::Entity_ID_List ids1;
  mesh->get_set_entities(setnames[0], Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::OWNED, &ids1);

  // construct the surface mesh (as a 3D submanifold)
  auto surfmesh = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(mesh,ids1,Amanzi::AmanziMesh::FACE,false,mesh->get_comm()));

  // "Extract" all cells of this surface mesh into another 2D mesh
  setnames.clear();
  setnames.push_back("Whole Surface");
  Amanzi::AmanziMesh::Entity_ID_List ids2;
  surfmesh->get_set_entities(setnames[0], Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED, &ids2);

  Amanzi::AmanziMesh::Mesh_MSTK surfmesh2D(surfmesh,ids2,Amanzi::AmanziMesh::CELL,true,surfmesh->get_comm());
  CHECK_EQUAL(surfmesh2D.space_dimension(),2);

  // Number of cells (quadrilaterals) in surface mesh
  int ncells_surf2D = surfmesh2D.num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(16,ncells_surf2D);
      
  // Number of "faces" (edges) in surface mesh
  int nfaces_surf2D = surfmesh2D.num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(40,nfaces_surf2D);

  // Number of nodes in surface mesh
  int nnodes_surf2D = surfmesh2D.num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(25,nnodes_surf2D);

  int exp_parent_faces[16] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};

  int *found = new int[ncells_surf2D];
  for (int k = 0; k < ncells_surf2D; k++) {
    Amanzi::AmanziMesh::Entity_ID parent = surfmesh2D.entity_get_parent(Amanzi::AmanziMesh::CELL,k);
        
    found[k] = 0;
    for (int kk = 0; kk < ncells_surf2D; kk++) {
      if (exp_parent_faces[kk] == parent) {        
        Amanzi::AmanziGeometry::Point centroid1 = surfmesh->cell_centroid(parent);
        Amanzi::AmanziGeometry::Point centroid2 = surfmesh2D.cell_centroid(k);
        CHECK_EQUAL(2,centroid2.dim());
        CHECK_CLOSE(centroid1[0],centroid2[0],1.0e-10);
        CHECK_CLOSE(centroid1[1],centroid2[1],1.0e-10);
        found[k]++;
      }
    }
  }
  
  for (int k = 0; k < ncells_surf2D; k++)
    CHECK_EQUAL(1,found[k]);
  delete [] found;

  // check the sets
  Amanzi::AmanziMesh::Entity_ID_List setents;  
  surfmesh2D.get_set_entities("Top Surface 2D",Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::Parallel_type::OWNED,&setents);
  
  Amanzi::AmanziMesh::Entity_ID expents[4] = {11,21,30,39};
  double exp_centroid[4][2] = {{0.125,1.0},{0.375,1.0},{0.625,1.0},{0.875,1.0}};

  for (int k = 0; k < setents.size(); k++) {
    Amanzi::AmanziGeometry::Point centroid1 = surfmesh2D.face_centroid(setents[k]);
    CHECK_EQUAL(expents[k],setents[k]);
    CHECK_CLOSE(exp_centroid[k][0],centroid1[0],1.0e-10);
    CHECK_CLOSE(exp_centroid[k][1],centroid1[1],1.0e-10);
  } 
}


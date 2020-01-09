#include <UnitTest++.h>
#include <fstream>

#include "../Mesh_MSTK.hh"


#include "Epetra_Map.h"
#include "AmanziComm.hh"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Array.hpp"

#include "mpi.h"

// startup function


// Extract some surfaces as-is from 3D mesh, generated mesh.
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

  // get a list of the top surface
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

  // parent entities
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


// Doing these as a test suite to really see what breaks.   Failures all throw unfortunately.
SUITE(Extract_Surface_MSTK_Sets) {

struct test_fixture {
  test_fixture() {
    std::string filename("test/hex_3x3x3_sets.exo");
    comm = Amanzi::getDefaultComm();
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

    gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));
    mesh = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(filename.c_str(), comm, gm));
  }

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> extract(bool flatten) {
    // get the surface set
    std::vector<std::string> setnames;
    setnames.push_back(std::string("Top Surface"));
    Amanzi::AmanziMesh::Entity_ID_List ids1;
    mesh->get_set_entities(setnames[0], Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::OWNED, &ids1);

    // extract a surface mesh from the 3D mesh using the set
    surfmesh = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(mesh,ids1,Amanzi::AmanziMesh::FACE, flatten, mesh->get_comm()));
    return surfmesh;
  }
    
  
  Amanzi::Comm_ptr_type comm;
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm;
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> surfmesh;
};

// Check the basic extraction geometry, topology.
TEST_FIXTURE(test_fixture, Extract_Surface_MSTK3_BASIC_EXTRACTION) {
  extract(false);
  
  // Number of cells (quadrilaterals) in surface mesh
  int ncells_surf = surfmesh->num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(9,ncells_surf);
      
  // Number of "faces" (edges) in surface mesh
  int nfaces_surf = surfmesh->num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(24,nfaces_surf);

  // Number of nodes in surface mesh
  int nnodes_surf = surfmesh->num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(16,nnodes_surf);

  // parent ids are the GIDs of the parent mesh's corresponding entities
  int exp_parent_faces[9] = {79,83,91,87,94,101,97,104,107};

  int *found = new int[ncells_surf];
  for (int k = 0; k < ncells_surf; k++) {
    Amanzi::AmanziMesh::Entity_ID parent = surfmesh->entity_get_parent(Amanzi::AmanziMesh::CELL,k);
        
    found[k] = 0;
    for (int kk = 0; kk < ncells_surf; kk++) {
      if (exp_parent_faces[kk] == parent) {        
        Amanzi::AmanziGeometry::Point centroid1 = mesh->face_centroid(parent);
        Amanzi::AmanziGeometry::Point centroid2 = surfmesh->cell_centroid(k);
        CHECK_ARRAY_EQUAL(centroid1,centroid2,3);
        found[k]++;
      }
    }
  }
  
  for (int k = 0; k < ncells_surf; k++)
    CHECK_EQUAL(1,found[k]);
  delete [] found;
}

TEST_FIXTURE(test_fixture, Extract_Surface_MSTK3_SIDE_FACES) {
  extract(false);

  // Test if the labeled set was inherited correctly and if we get the
  // right entities for this set
  //
  bool is_valid = surfmesh->valid_set_name("Side Surface", Amanzi::AmanziMesh::FACE);
  CHECK(is_valid);
  Amanzi::AmanziMesh::Entity_ID_List setents;
  // In this case, a face set in the parent becomes a face set in the surface
  surfmesh->get_set_entities("Side Surface",Amanzi::AmanziMesh::FACE,
                            Amanzi::AmanziMesh::Parallel_type::OWNED,&setents);
  CHECK_EQUAL(3, setents.size());
}

TEST_FIXTURE(test_fixture, Extract_Surface_MSTK3_SIDE_FACES_NO_CELLS) {
  extract(false);

  // The side surface does not intersect with the top surface, so there are no
  // cells.  This is either valid, and get_set_entities() does not throw but
  // returns an empty list, or is not valid.
  bool is_valid = surfmesh->valid_set_name("Side Surface", Amanzi::AmanziMesh::CELL);
  if (is_valid) {  
    Amanzi::AmanziMesh::Entity_ID_List setents;
    surfmesh->get_set_entities("Side Surface",Amanzi::AmanziMesh::CELL,
            Amanzi::AmanziMesh::Parallel_type::OWNED,&setents);
    CHECK_EQUAL(0, setents.size());
  }
}

TEST_FIXTURE(test_fixture, Extract_Surface_MSTK3_TOP_FACES) {
  extract(false);

  // In this case, a face set in the parent becomes a cell set in the surface
  bool is_valid = surfmesh->valid_set_name("Top Surface", Amanzi::AmanziMesh::CELL);
  CHECK(is_valid);
  Amanzi::AmanziMesh::Entity_ID_List setents;
  surfmesh->get_set_entities("Top Surface",Amanzi::AmanziMesh::CELL,
                            Amanzi::AmanziMesh::Parallel_type::OWNED,&setents);
  CHECK_EQUAL(9, setents.size());
}

TEST_FIXTURE(test_fixture, Extract_Surface_MSTK3_TOP_FACES_NO_FACES) {
  extract(false);

  // In this case, a face set in the parent becomes a face set in the surface
  //
  // The top surface is only cells, so there are no faces.  This is either
  // valid, and get_set_entities() does not throw but returns an empty list, or
  // is not valid.
  bool is_valid = surfmesh->valid_set_name("Side Surface", Amanzi::AmanziMesh::FACE);
  if (is_valid) {  
    Amanzi::AmanziMesh::Entity_ID_List setents;
    surfmesh->get_set_entities("Top Surface",Amanzi::AmanziMesh::FACE,
            Amanzi::AmanziMesh::Parallel_type::OWNED,&setents);
    CHECK_EQUAL(0, setents.size());
  }
}

TEST_FIXTURE(test_fixture, Extract_Surface_MSTK3_TOP_CELLS) {
  extract(false);

  // In this case, a cell set in the parent becomes a cell set in the surface
  bool is_valid = surfmesh->valid_set_name("Region 1", Amanzi::AmanziMesh::CELL);
  CHECK(is_valid);
  Amanzi::AmanziMesh::Entity_ID_List setents;
  surfmesh->get_set_entities("Region 1",Amanzi::AmanziMesh::CELL,
                            Amanzi::AmanziMesh::Parallel_type::OWNED,&setents);
  CHECK_EQUAL(9, setents.size());
}

TEST_FIXTURE(test_fixture, Extract_Surface_MSTK3_BOTTOM_CELLS) {
  extract(false);

  // In this case, a cell set in the parent becomes a cell set in the surface,
  // but there is no intersection.  This is either valid, but doesn't throw and
  // returns an empty set, or is not valid.  Either is OK?
  bool is_valid = surfmesh->valid_set_name("Region 2", Amanzi::AmanziMesh::CELL);
  if (is_valid) {  
    Amanzi::AmanziMesh::Entity_ID_List setents;
    surfmesh->get_set_entities("Region 2",Amanzi::AmanziMesh::CELL,
            Amanzi::AmanziMesh::Parallel_type::OWNED,&setents);
    CHECK_EQUAL(0, setents.size());
  }
}


// Check the basic extraction geometry, topology.
TEST_FIXTURE(test_fixture, Extract_Flatten_Surface_MSTK3_BASIC_EXTRACTION) {
  extract(true);
  
  // Number of cells (quadrilaterals) in surface mesh
  int ncells_surf = surfmesh->num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(9,ncells_surf);
      
  // Number of "faces" (edges) in surface mesh
  int nfaces_surf = surfmesh->num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(24,nfaces_surf);

  // Number of nodes in surface mesh
  int nnodes_surf = surfmesh->num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(16,nnodes_surf);

  // parent ids are the GIDs of the parent mesh's corresponding entities
  int exp_parent_faces[9] = {79,83,91,87,94,101,97,104,107};

  int *found = new int[ncells_surf];
  for (int k = 0; k < ncells_surf; k++) {
    Amanzi::AmanziMesh::Entity_ID parent = surfmesh->entity_get_parent(Amanzi::AmanziMesh::CELL,k);
        
    found[k] = 0;
    for (int kk = 0; kk < ncells_surf; kk++) {
      if (exp_parent_faces[kk] == parent) {        
        Amanzi::AmanziGeometry::Point centroid1 = mesh->face_centroid(parent);
        Amanzi::AmanziGeometry::Point centroid2 = surfmesh->cell_centroid(k);
        CHECK_ARRAY_EQUAL(centroid1,centroid2,2);
        found[k]++;
      }
    }
  }
  
  for (int k = 0; k < ncells_surf; k++)
    CHECK_EQUAL(1,found[k]);
  delete [] found;
}

TEST_FIXTURE(test_fixture, Extract_Flatten_Surface_MSTK3_SIDE_FACES) {
  extract(true);

  // Test if the labeled set was inherited correctly and if we get the
  // right entities for this set
  //
  bool is_valid = surfmesh->valid_set_name("Side Surface", Amanzi::AmanziMesh::FACE);
  CHECK(is_valid);
  Amanzi::AmanziMesh::Entity_ID_List setents;
  // In this case, a face set in the parent becomes a face set in the surface
  surfmesh->get_set_entities("Side Surface",Amanzi::AmanziMesh::FACE,
                            Amanzi::AmanziMesh::Parallel_type::OWNED,&setents);
  CHECK_EQUAL(3, setents.size());
}

TEST_FIXTURE(test_fixture, Extract_Flatten_Surface_MSTK3_SIDE_FACES_NO_CELLS) {
  extract(true);

  // The side surface does not intersect with the top surface, so there are no
  // cells.  This is either valid, and get_set_entities() does not throw but
  // returns an empty list, or is not valid.
  bool is_valid = surfmesh->valid_set_name("Side Surface", Amanzi::AmanziMesh::CELL);
  if (is_valid) {  
    Amanzi::AmanziMesh::Entity_ID_List setents;
    surfmesh->get_set_entities("Side Surface",Amanzi::AmanziMesh::CELL,
            Amanzi::AmanziMesh::Parallel_type::OWNED,&setents);
    CHECK_EQUAL(0, setents.size());
  }
}

TEST_FIXTURE(test_fixture, Extract_Flatten_Surface_MSTK3_TOP_FACES) {
  extract(true);

  // In this case, a face set in the parent becomes a cell set in the surface
  bool is_valid = surfmesh->valid_set_name("Top Surface", Amanzi::AmanziMesh::CELL);
  CHECK(is_valid);
  Amanzi::AmanziMesh::Entity_ID_List setents;
  surfmesh->get_set_entities("Top Surface",Amanzi::AmanziMesh::CELL,
                            Amanzi::AmanziMesh::Parallel_type::OWNED,&setents);
  CHECK_EQUAL(9, setents.size());
}

TEST_FIXTURE(test_fixture, Extract_Flatten_Surface_MSTK3_TOP_FACES_NO_FACES) {
  extract(true);

  // In this case, a face set in the parent becomes a face set in the surface
  //
  // The top surface is only cells, so there are no faces.  This is either
  // valid, and get_set_entities() does not throw but returns an empty list, or
  // is not valid.
  bool is_valid = surfmesh->valid_set_name("Side Surface", Amanzi::AmanziMesh::FACE);
  if (is_valid) {  
    Amanzi::AmanziMesh::Entity_ID_List setents;
    surfmesh->get_set_entities("Top Surface",Amanzi::AmanziMesh::FACE,
            Amanzi::AmanziMesh::Parallel_type::OWNED,&setents);
    CHECK_EQUAL(0, setents.size());
  }
}

TEST_FIXTURE(test_fixture, Extract_Flatten_Surface_MSTK3_TOP_CELLS) {
  extract(true);

  // In this case, a cell set in the parent becomes a cell set in the surface
  bool is_valid = surfmesh->valid_set_name("Region 1", Amanzi::AmanziMesh::CELL);
  CHECK(is_valid);
  Amanzi::AmanziMesh::Entity_ID_List setents;
  surfmesh->get_set_entities("Region 1",Amanzi::AmanziMesh::CELL,
                            Amanzi::AmanziMesh::Parallel_type::OWNED,&setents);
  CHECK_EQUAL(9, setents.size());
}

TEST_FIXTURE(test_fixture, Extract_Flatten_Surface_MSTK3_BOTTOM_CELLS) {
  extract(true);

  // In this case, a cell set in the parent becomes a cell set in the surface,
  // but there is no intersection.  This is either valid, but doesn't throw and
  // returns an empty set, or is not valid.  Either is OK?
  bool is_valid = surfmesh->valid_set_name("Region 2", Amanzi::AmanziMesh::CELL);
  if (is_valid) {  
    Amanzi::AmanziMesh::Entity_ID_List setents;
    surfmesh->get_set_entities("Region 2",Amanzi::AmanziMesh::CELL,
            Amanzi::AmanziMesh::Parallel_type::OWNED,&setents);
    CHECK_EQUAL(0, setents.size());
  }
}


TEST_FIXTURE(test_fixture, Extract_Flatten_Surface_MSTK3_BAD_NAME) {
  extract(false);

  // In this case, a cell set in the parent becomes a cell set in the surface,
  // but there is no intersection.  This is either valid, but doesn't throw and
  // returns an empty set, or is not valid.  Either is OK?
  bool is_valid = surfmesh->valid_set_name("Region Not Here", Amanzi::AmanziMesh::CELL);
  CHECK(!is_valid);
  is_valid = surfmesh->valid_set_name("Region Not Here", Amanzi::AmanziMesh::FACE);
  CHECK(!is_valid);
}

}


#include <UnitTest++.h>
#include <iostream>

#include "../Mesh_MSTK.hh"


#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Array.hpp"

#include "mpi.h"


TEST(Extract_Surface_MSTK1)
{

  Teuchos::RCP<Epetra_MpiComm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));


  Teuchos::ParameterList parameterlist;

 
  // create a sublist name Regions and put a reference to it in
  // reg_spec and other sublists as references. Turns out it is
  // important to define reg_spec and other lists below as references
  // - otherwise, a new copy is made of the sublist that is retrieved

  Teuchos::ParameterList& reg_spec = parameterlist.sublist("Regions"); 
  
  Teuchos::ParameterList& top_surface = reg_spec.sublist("Top Surface");
  Teuchos::ParameterList& top_surface_def = top_surface.sublist("Region: Plane");
  Teuchos::Array<double> loc1 = Teuchos::tuple(0.0,0.0,1.0);
  Teuchos::Array<double> dir1 = Teuchos::tuple(0.0,0.0,-1.0);
  top_surface_def.set< Teuchos::Array<double> >("Location",loc1);
  top_surface_def.set< Teuchos::Array<double> >("Direction",dir1);

  Teuchos::ParameterList& right_surface = reg_spec.sublist("Right Surface");
  Teuchos::ParameterList& right_surface_def = right_surface.sublist("Region: Plane");
  Teuchos::Array<double> loc2 = Teuchos::tuple(1.0,0.0,0.0);
  Teuchos::Array<double> dir2 = Teuchos::tuple(1.0,0.0,0.0);
  right_surface_def.set< Teuchos::Array<double> >("Location",loc2);
  right_surface_def.set< Teuchos::Array<double> >("Direction",dir2);


  Teuchos::writeParameterListToXmlOStream(parameterlist,std::cout);


  Amanzi::AmanziGeometry::GeometricModelPtr gm = new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, comm.get());

  // Generate a mesh consisting of 3x3x3 elements 

  Amanzi::AmanziMesh::Mesh_MSTK mesh(0,0,0,1,1,1,3,3,3,comm.get(),gm);


  std::vector<std::string> setnames;
  setnames.push_back(std::string("Top Surface"));
  setnames.push_back(std::string("Right Surface"));

  Amanzi::AmanziMesh::Mesh_MSTK surfmesh(mesh,setnames,Amanzi::AmanziMesh::FACE);


  // Number of cells (quadrilaterals) in surface mesh

  int ncells_surf = surfmesh.num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);
  CHECK_EQUAL(18,ncells_surf);
      
  // Number of "faces" (edges) in surface mesh

  int nfaces_surf = surfmesh.num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED);
  CHECK_EQUAL(45,nfaces_surf);

  // Number of nodes in surface mesh

  int nnodes_surf = surfmesh.num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::OWNED);
  CHECK_EQUAL(28,nnodes_surf);


  int exp_parent_faces[18] = {27,28,29,30,31,32,33,34,35,99,100,101,102,103,104,105,106,107};

  int *found = new int[ncells_surf];
  for (int k = 0; k < ncells_surf; k++) {
    Amanzi::AmanziMesh::Entity_ID parent = surfmesh.entity_get_parent(Amanzi::AmanziMesh::CELL,k);
        
    found[k] = 0;
    for (int kk = 0; kk < ncells_surf; kk++) {
      if (exp_parent_faces[kk] == parent) {        
        Amanzi::AmanziGeometry::Point centroid1 = mesh.face_centroid(parent);
        Amanzi::AmanziGeometry::Point centroid2 = surfmesh.cell_centroid(k);
        CHECK_ARRAY_EQUAL(centroid1,centroid2,3);
        found[k]++;
      }
    }
  }
  
  for (int k = 0; k < ncells_surf; k++)
    CHECK_EQUAL(1,found[k]);


  // Once we can make RegionFactory work with reference counted pointers 
  // we can get rid of this code

  for (int i = 0; i < gm->Num_Regions(); i++)
    delete (gm->Region_i(i));
  delete gm;

  
}

TEST(Extract_Surface_MSTK2)
{

  Teuchos::RCP<Epetra_MpiComm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));


  Teuchos::ParameterList parameterlist;

 
  // create a sublist name Regions and put a reference to it in
  // reg_spec and other sublists as references. Turns out it is
  // important to define reg_spec and other lists below as references
  // - otherwise, a new copy is made of the sublist that is retrieved

  Teuchos::ParameterList& reg_spec = parameterlist.sublist("Regions"); 

  Teuchos::ParameterList& top_surface = reg_spec.sublist("Top Surface");
  Teuchos::ParameterList& top_surface_def = top_surface.sublist("Region: Plane");
  Teuchos::Array<double> loc1 = Teuchos::tuple(0.0,0.0,1.0);
  Teuchos::Array<double> dir1 = Teuchos::tuple(-1/sqrt(2.0),0.0,1/sqrt(2.0));
  top_surface_def.set< Teuchos::Array<double> >("Location",loc1);
  top_surface_def.set< Teuchos::Array<double> >("Direction",dir1);

  

  Teuchos::writeParameterListToXmlOStream(parameterlist,std::cout);


  Amanzi::AmanziGeometry::GeometricModelPtr gm = new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, comm.get());

  // Generate a mesh consisting of 3x3x3 elements 

  Amanzi::AmanziMesh::Mesh_MSTK mesh(0.0,0.0,0.0,1.0,1.0,1.0,4,4,4,comm.get(),gm);


  // Perturb some nodes

  int nv = mesh.num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::OWNED);

  for (int i = 0; i < nv; i++) {
    Amanzi::AmanziGeometry::Point pt;
    mesh.node_get_coordinates(i,&pt);
    if (pt[2] == 1.0) {
      pt[2] = pt[2]+pt[0];
      mesh.node_set_coordinates(i,pt);
    }
  }
  

  std::vector<std::string> setnames;
  setnames.push_back(std::string("Top Surface"));

  // Extract surface mesh while projecting to 2D
  Amanzi::AmanziMesh::Mesh_MSTK surfmesh(mesh,setnames,Amanzi::AmanziMesh::FACE,true,false);

  CHECK_EQUAL(surfmesh.space_dimension(),2);

  // Number of cells (quadrilaterals) in surface mesh

  int ncells_surf = surfmesh.num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);
  CHECK_EQUAL(16,ncells_surf);
      
  // Number of "faces" (edges) in surface mesh

  int nfaces_surf = surfmesh.num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED);
  CHECK_EQUAL(40,nfaces_surf);

  // Number of nodes in surface mesh

  int nnodes_surf = surfmesh.num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::OWNED);
  CHECK_EQUAL(25,nnodes_surf);


  int exp_parent_faces[16] = {224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239};

  int *found = new int[ncells_surf];
  for (int k = 0; k < ncells_surf; k++) {
    Amanzi::AmanziMesh::Entity_ID parent = surfmesh.entity_get_parent(Amanzi::AmanziMesh::CELL,k);
        
    found[k] = 0;
    for (int kk = 0; kk < ncells_surf; kk++) {
      if (exp_parent_faces[kk] == parent) {        
        Amanzi::AmanziGeometry::Point centroid1 = mesh.face_centroid(parent);
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


  // Once we can make RegionFactory work with reference counted pointers 
  // we can get rid of this code

  for (int i = 0; i < gm->Num_Regions(); i++)
    delete (gm->Region_i(i));
  delete gm;

  
}


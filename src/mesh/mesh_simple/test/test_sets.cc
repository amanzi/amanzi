#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include "../Mesh_simple.hh"
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"
#include "GenerationSpec.hh"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

SUITE (MeshSimple) {
TEST(SETS) {
  
  using namespace std;

#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  std::string expcsetnames[3] = {"Bottom Box", "Bottom+Middle Box",
                                 "Vertical Box"};
  unsigned int csetsize, expcsetsizes[5] = {9,18,9};
  
  unsigned int expcsetcells[3][18] = {{0,3,6,1,4,7,2,5,8},
				      {0,9,3,12,6,15,1,10,4,13,7,16,2,11,5,14,8,17},
				      {1,10,19,4,13,22,7,16,25}};

  std::string expfsetnames[2] = {"ZLO FACE Plane", "YLO FACE Box"};

  unsigned int fsetsize, expfsetsizes[2] = {9,3};


  unsigned int expfsetfaces[6][9] = {{0,3,6,1,4,7,2,5,8},
				     {60,61,62,-1,-1,-1,-1,-1,-1}};


  std::string infilename = "test/hex_3x3x3.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());

  Amanzi::AmanziGeometry::GeometricModelPtr gm = new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec);

  // Create a mesh consisting of 3x3x3 elements (4x4x4 nodes)

  Amanzi::AmanziMesh::Mesh_simple mesh(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3, comm, gm); 

  Teuchos::ParameterList::ConstIterator i;
  for (i = reg_spec.begin(); i != reg_spec.end(); i++) {
        const std::string reg_name = reg_spec.name(i);     

    Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

    // See if the geometric model has a region by this name
  
    Amanzi::AmanziGeometry::RegionPtr reg = gm->FindRegion(reg_name);

    CHECK(reg != NULL);

    // Do their names match ?

    CHECK_EQUAL(reg->name(),reg_name);


    // Get the region info directly from the XML and compare
  
    Teuchos::ParameterList::ConstIterator j = reg_params.begin(); 

    std::string shape = reg_params.name(j);

    if (shape == "Region: Plane") {

      // Do we have a valid sideset by this name

      CHECK(mesh.valid_set_name(reg_name,Amanzi::AmanziMesh::FACE));

      int j;
      for (j = 0; j < 2; j++) {
        if (expfsetnames[j] == reg_name) break;
      }

      CHECK(j < 2);


      // Verify that we can get the right number of entities in the set

      int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::USED);

      CHECK_EQUAL(expfsetsizes[j],set_size);


      // Verify that we can get the correct set entities
     
      Amanzi::AmanziMesh::Entity_ID_List setents;
      mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::USED,&setents);

      CHECK_ARRAY_EQUAL(expfsetfaces[j],setents,set_size);

    }
    else if (shape == "Region: Box") {

      Teuchos::ParameterList box_params = reg_params.sublist(shape);
      Teuchos::Array<double> pmin = box_params.get< Teuchos::Array<double> >("Low Coordinate");
      Teuchos::Array<double> pmax = box_params.get< Teuchos::Array<double> >("High Coordinate");

      if (pmin[0] == pmax[0] || pmin[1] == pmax[1] || pmin[2] == pmax[2])
	{
	  // This is a reduced dimensionality box - request a faceset

	  // Do we have a valid sideset by this name

	  CHECK(mesh.valid_set_name(reg_name,Amanzi::AmanziMesh::FACE));
	  
	  int j;
	  for (j = 0; j < 2; j++) {
	    if (expfsetnames[j] == reg_name) break;
	  }
	  
	  CHECK(j < 2);
	  
	  
	  // Verify that we can get the right number of entities in the set
	  
	  int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::USED);
	  
	  CHECK_EQUAL(expfsetsizes[j],set_size);
	  
	  
	  // Verify that we can get the correct set entities
	  
	  Amanzi::AmanziMesh::Entity_ID_List setents;
	  mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::USED,&setents);
	  
	  CHECK_ARRAY_EQUAL(expfsetfaces[j],setents,set_size);	  
	}
      else 
	{
	  // Do we have a valid cellset by this name
	  
	  CHECK(mesh.valid_set_name(reg_name,Amanzi::AmanziMesh::CELL));
	  
	  // Find the expected cell set info corresponding to this name 
	  
	  int j;
	  for (j = 0; j < 3; j++)
	    if (reg_name == expcsetnames[j]) break;
	  
	  CHECK(j < 3);
	  
	  // Verify that we can get the right number of entities in the set
	  
	  int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::USED);
	  
	  CHECK_EQUAL(expcsetsizes[j],set_size);
	  
	  // Verify that we can get the correct set entities
	  
	  Amanzi::AmanziMesh::Entity_ID_List setents;
	  mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::USED,&setents);
	  
	  CHECK_ARRAY_EQUAL(expcsetcells[j],setents,set_size);
	}
    }
    else if (shape == "Region: Labeled Set") {

      std::cerr << "Mesh framework cannot do labeled sets" << std::endl;

      CHECK(false); 

    }
  }

}
}

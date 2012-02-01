#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"


#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "mpi.h"

// Unless this example is enhanced, it does lesser testing than test_hex_3x3x2.cc

TEST(MSTK_HEX_4x4x4_SETS_4P)
{
  int rank, size;

  std::string expcsetnames[8] = {"Bottom LS", "Middle LS", "Top LS", 
                                 "Bottom+Middle Box", "Top Box",
                                 "Bottom ColFunc", "Middle ColFunc", "Top ColFunc"};

  int csetsize, expcsetsizes[4][8] = {{0,0,6,0,6,0,0,6},
				      {1,3,3,4,3,1,3,3},
				      {3,4,0,7,0,3,4,0},
				      {5,2,0,7,0,5,2,0}};
  
  int expcsetcells[4][8][9] = 
    {
      {{-1,-1,-1,-1,-1,-1,-1,-1,-1},
       {-1,-1,-1,-1,-1,-1,-1,-1,-1},
       {0,1,2,3,4,5,-1,-1,-1},
       {-1,-1,-1,-1,-1,-1,-1,-1,-1},
       {0,1,2,3,4,5,-1,-1,-1},
       {-1,-1,-1,-1,-1,-1,-1,-1,-1},
       {-1,-1,-1,-1,-1,-1,-1,-1,-1},
       {0,1,2,3,4,5,-1,-1,-1}},

      {{0,-1,-1,-1,-1,-1,-1,-1,-1},
       {1,2,3,-1,-1,-1,-1,-1,-1},
       {4,5,6,-1,-1,-1,-1,-1,-1},
       {0,1,2,3,-1,-1,-1,-1,-1},
       {4,5,6,-1,-1,-1,-1,-1,-1},
       {0,-1,-1,-1,-1,-1,-1,-1,-1},
       {1,2,3,-1,-1,-1,-1,-1,-1},
       {4,5,6,-1,-1,-1,-1,-1,-1}},
      
      {{0,1,2,-1,-1,-1,-1,-1,-1},
       {3,4,5,6,-1,-1,-1,-1,-1},
       {-1,-1,-1,-1,-1,-1,-1,-1,-1},
       {0,1,2,3,4,5,6,-1,-1},
       {-1,-1,-1,-1,-1,-1,-1,-1},
       {0,1,2,-1,-1,-1,-1,-1,-1},
       {3,4,5,6,-1,-1,-1,-1,-1},
       {-1,-1,-1,-1,-1,-1,-1,-1,-1}},

      {{0,1,2,3,4,-1,-1,-1,-1},
       {5,6,-1,-1,-1,-1,-1,-1,-1},
       {-1,-1,-1,-1,-1,-1,-1,-1,-1},
       {0,1,2,3,4,5,6,-1,-1},
       {-1,-1,-1,-1,-1,-1,-1,-1,-1},
       {0,1,2,3,4,-1,-1,-1,-1},
       {5,6,-1,-1,-1,-1,-1,-1,-1},
       {-1,-1,-1,-1,-1,-1,-1,-1,-1}}
    };


  std::string expfsetnames[4] = {"Face 101",  
				  "Face 30004",
                                  "ZLO FACE Plane", 
				 "YLO FACE Box"};

  unsigned int expfsetids[4]={101,30004,0,0};

  int fsetsize, expfsetsizes[4][4] = {{0,2,0,2},
				      {1,1,1,2},
				      {3,0,3,3},
				      {5,0,5,2}};
  
  int expfsetfaces[4][4][5] = {{{-1,-1,-1,-1,-1},
				{21,26,-1,-1,-1},
				{-1,-1,-1,-1,-1},
				{0,6,-1,-1,-1}},
			       
			       {{4,-1,-1,-1,-1},
				{28,-1,-1,-1,-1},
				{4,-1,-1,-1,-1},
				{6,21,-1,-1,-1}},
			       
			       {{4,9,14,-1,-1},
				{-1,-1,-1,-1,-1},
				{4,9,14,-1,-1},
				{0,16,20,-1,-1}},
			       
			       {{4,10,7,12,15},
				{-1,-1,-1,-1,-1},
				{4,7,10,12,15},
				{0,5,-1,-1,-1}}};

  
  int initialized;
  MPI_Initialized(&initialized);
  
  if (!initialized)
    MPI_Init(NULL,NULL);

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  if (size != 4) {
    cerr << "Test must be run with 4 processors" << std::endl;
  }
  CHECK_EQUAL(4,size);


  std::string infilename = "test/hex_4x4x4_4P.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());

  Epetra_MpiComm ecomm(MPI_COMM_WORLD);

  Amanzi::AmanziGeometry::GeometricModelPtr gm = new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, &ecomm);

  // Load a mesh consisting of 3x3x3 elements (4x4x4 nodes)

  Amanzi::AmanziMesh::Mesh_MSTK mesh("test/hex_4x4x4_ss.exo",MPI_COMM_WORLD,3,gm);

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
      for (j = 0; j < 4; j++) {
        if (expfsetnames[j] == reg_name) break;
      }

      CHECK(j < 4);


      // Verify that we can get the right number of entities in the set

      int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED);

      CHECK_EQUAL(expfsetsizes[rank][j],set_size);


      // Verify that we can get the correct set entities
     
      Amanzi::AmanziMesh::Entity_ID_List setents;
      mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED,&setents);

      CHECK_ARRAY_EQUAL(expfsetfaces[rank][j],setents,set_size);

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
	  for (j = 0; j < 4; j++) {
	    if (expfsetnames[j] == reg_name) break;
	  }
	  
	  CHECK(j < 4);
	  
	  
	  // Verify that we can get the right number of entities in the set
	  
	  int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED);

	  if (expfsetsizes[rank][j] != set_size)
	    std::cout << "Wrong set size for " << reg_name << " on proc " << rank << std::endl;
	  CHECK_EQUAL(expfsetsizes[rank][j],set_size);
	  
	  
	  // Verify that we can get the correct set entities
	  
	  Amanzi::AmanziMesh::Entity_ID_List setents;
	  mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED,&setents);
	  
	  CHECK_ARRAY_EQUAL(expfsetfaces[rank][j],setents,set_size);	  
	}
      else 
	{
	  // Do we have a valid cellset by this name
	  
	  CHECK(mesh.valid_set_name(reg_name,Amanzi::AmanziMesh::CELL));
	  
	  // Find the expected cell set info corresponding to this name 
	  
	  int j;
	  for (j = 0; j < 8; j++)
	    if (reg_name == expcsetnames[j]) break;
	  
	  CHECK(j < 8);
	  
	  // Verify that we can get the right number of entities in the set
	  
	  int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);
	  
	  CHECK_EQUAL(expcsetsizes[rank][j],set_size);
	  
	  // Verify that we can get the correct set entities
	  
	  Amanzi::AmanziMesh::Entity_ID_List setents;
	  mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED,&setents);
	  
	  CHECK_ARRAY_EQUAL(expcsetcells[rank][j],setents,set_size);
	}
    }
    else if (shape == "Region: Labeled Set") {

      Teuchos::ParameterList lsparams = reg_params.sublist(shape);

      // Find the entity type in this parameter list

      std::string entity_type = lsparams.get<std::string>("Entity");

      if (entity_type == "Face") {

	// Do we have a valid sideset by this name

	CHECK(mesh.valid_set_name(reg_name,Amanzi::AmanziMesh::FACE));

        // Find the expected face set info corresponding to this name

        int j;
        for (j = 0; j < 4; j++)
          if (reg_name == expfsetnames[j]) break;

	if (j >= 4) 
	  std::cerr << "Cannot find regname " << reg_name << "on processor " << rank << std::endl;
        CHECK(j < 4);
	
	// Verify that we can get the right number of entities in the set
	
	int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED);
		
        CHECK_EQUAL(expfsetsizes[rank][j],set_size);

	// Verify that we can get the correct set entities
	
        Amanzi::AmanziMesh::Entity_ID_List setents;
	mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED,&setents);

        CHECK_ARRAY_EQUAL(expfsetfaces[rank][j],setents,set_size);

      }
      else if (entity_type == "Cell") {

	// Do we have a valid sideset by this name

	CHECK(mesh.valid_set_name(reg_name,Amanzi::AmanziMesh::CELL));
	
        // Find the expected face set info corresponding to this name

        int j;
        for (j = 0; j < 8; j++)
          if (reg_name == expcsetnames[j]) break;

        CHECK(j < 8);
	
	// Verify that we can get the right number of entities in the set
	
	int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);

        CHECK_EQUAL(expcsetsizes[rank][j],set_size);
	
	// Verify that we can get the correct set entities
	
        Amanzi::AmanziMesh::Entity_ID_List setents;
	mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED,&setents);

        CHECK_ARRAY_EQUAL(expcsetcells[rank][j],setents,set_size);
      }

    }
    else if (shape == "Region: Color Function") {

      // Do we have a valid cellset by this name

      CHECK(mesh.valid_set_name(reg_name,Amanzi::AmanziMesh::CELL));
	
      // Find the expected cell set info corresponding to this name

      int j;
      for (j = 0; j < 8; j++)
        if (reg_name == expcsetnames[j]) break;

      CHECK(j < 8);
	
      // Verify that we can get the right number of entities in the set
	
      int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);

      CHECK_EQUAL(expcsetsizes[rank][j],set_size);
	
      // Verify that we can get the correct set entities
	
      Amanzi::AmanziMesh::Entity_ID_List setents;
      mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED,&setents);
      
      CHECK_ARRAY_EQUAL(expcsetcells[rank][j],setents,set_size);

    }
  }


  // Once we can make RegionFactory work with reference counted pointers 
  // we can get rid of this code

  for (int i = 0; i < gm->Num_Regions(); i++)
    delete (gm->Region_i(i));
  delete gm;

}


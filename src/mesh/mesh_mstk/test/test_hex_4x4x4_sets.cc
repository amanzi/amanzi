#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"


#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "mpi.h"

// Unless this example is enhanced, it does lesser testing than test_hex_3x3x2.cc

TEST(MSTK_HEX_4x4x4_SETS)
{

  std::string expcsetnames[9] = {"Bottom LS", "Middle LS", "Top LS", 
                                 "Bottom+Middle Box", "Top Box",
                                 "Sample Point InCell", "Sample Point OnFace",
                                 "Sample Point OnEdge", "Sample Point OnVertex"};
  unsigned int csetsize, expcsetsizes[9] = {9,9,9,18,9,1,2,4,8};
  
  unsigned int expcsetcells[9][18] = {{0,1,2,3,4,5,6,7,8,0,0,0,0,0,0,0,0},
				      {9,10,11,12,13,14,15,16,17,0,0,0,0,0,0,0,0,0},
				      {18,19,20,21,22,23,24,25,26,0,0,0,0,0,0,0,0,0},
				      {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17},
				      {18,19,20,21,22,23,24,25,26,0,0,0,0,0,0,0,0,0},
                                      {13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                      {12,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                      {9,10,12,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                      {0,1,3,4,9,10,12,13,0,0,0,0,0,0,0,0,0,0}};

  std::string expfsetnames[7] = {"Face 101", "Face 102", 
				  "Face 10005", "Face 20004", "Face 30004",
                                  "ZLO FACE Plane", "YLO FACE Box"};

  unsigned int expfsetids[7]={101,102,10005,20004,30004,0,0};
  
  unsigned int fsetsize, expfsetsizes[7] = {9,9,3,3,3,9,9};


  unsigned int expfsetfaces[7][9] = {{4,19,9,32,23,14,36,27,40},
   				      {0,6,42,11,47,75,51,80,84},
   				      {30,35,39,0,0,0,0,0,0},
   				      {66,70,73,0,0,0,0,0,0},
   				      {99,103,106,0,0,0,0,0,0},
				      {4,9,14,19,23,27,32,36,40},
				      {0,6,11,42,47,51,75,80,84}};

  std::string expnsetnames[2] = {"INTERIOR XY PLANE", "TOP BOX"};

  unsigned int nsetsize, expnsetsizes[2] = {16, 4};
  unsigned int expnsetnodes[2][16] = {{16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},
                                      {53,54,57,58,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}};

			   

  std::string infilename = "test/hex_4x4x4.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());

  Amanzi::AmanziGeometry::GeometricModelPtr gm = new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec);

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

      if (reg_name == "ZLO FACE Plane") {

        // Do we have a valid sideset by this name
        
        CHECK(mesh.valid_set_name(reg_name,Amanzi::AmanziMesh::FACE));
        
        int j;
        for (j = 0; j < 7; j++) {
          if (expfsetnames[j] == reg_name) break;
        }
        
        CHECK(j < 7);
        
        
        // Verify that we can get the right number of entities in the set
        
        int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED);
        
        CHECK_EQUAL(expfsetsizes[j],set_size);
        
        
        // Verify that we can get the correct set entities
        
        Amanzi::AmanziMesh::Entity_ID_List setents;
        mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED,&setents);
        
        CHECK_ARRAY_EQUAL(expfsetfaces[j],setents,set_size);
      }
      else if (reg_name == "INTERIOR XY PLANE") {

        // Do we have a valid nodeset by this name
        
        CHECK(mesh.valid_set_name(reg_name,Amanzi::AmanziMesh::NODE));
        
        int j;
        for (j = 0; j < 2; j++) {
          if (expnsetnames[j] == reg_name) break;
        }

        CHECK(j < 2);
        
        
        // Verify that we can get the right number of entities in the set
        
        int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::USED);
        
        CHECK_EQUAL(expnsetsizes[j],set_size);
        
        
        // Verify that we can get the correct set entities
        
        Amanzi::AmanziMesh::Entity_ID_List setents;
        mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::USED,&setents);
        
        CHECK_ARRAY_EQUAL(expnsetnodes[j],setents,set_size);
        
      }

    }
    else if (shape == "Region: Box") {

      Teuchos::ParameterList box_params = reg_params.sublist(shape);
      Teuchos::Array<double> pmin = box_params.get< Teuchos::Array<double> >("Low Coordinate");
      Teuchos::Array<double> pmax = box_params.get< Teuchos::Array<double> >("High Coordinate");

      if (pmin[0] == pmax[0] || pmin[1] == pmax[1] || pmin[2] == pmax[2])
	{          

	  // This is a reduced dimensionality box - request a faceset or nodeset

          if (reg_name == "YLO FACE BOX") {

            // Do we have a valid sideset by this name
            
            CHECK(mesh.valid_set_name(reg_name,Amanzi::AmanziMesh::FACE));
            
            int j;
            for (j = 0; j < 7; j++) {
              if (expfsetnames[j] == reg_name) break;
            }
            
            CHECK(j < 7);
            
            
            // Verify that we can get the right number of entities in the set
            
            int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED);
            
            CHECK_EQUAL(expfsetsizes[j],set_size);
            
            
            // Verify that we can get the correct set entities
            
            Amanzi::AmanziMesh::Entity_ID_List setents;
            mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED,&setents);
            
            CHECK_ARRAY_EQUAL(expfsetfaces[j],setents,set_size);	  
            
          }
          else if (reg_name == "TOP BOX") {

            // Do we have a valid set by this name
            
            CHECK(mesh.valid_set_name(reg_name,Amanzi::AmanziMesh::NODE));
            
            int j;
            for (j = 0; j < 2; j++) {
              if (expnsetnames[j] == reg_name) break;
            }
            
            CHECK(j < 2);
            
            
            // Verify that we can get the right number of entities in the set
            
            int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::USED);
            
            CHECK_EQUAL(expnsetsizes[j],set_size);
            
            
            // Verify that we can get the correct set entities
            
            Amanzi::AmanziMesh::Entity_ID_List setents;
            mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::USED,&setents);
            
            CHECK_ARRAY_EQUAL(expnsetnodes[j],setents,set_size);	  
                        
          }

	}
      else 
	{
	  // Do we have a valid cellset by this name
	  
	  CHECK(mesh.valid_set_name(reg_name,Amanzi::AmanziMesh::CELL));
	  
	  // Find the expected cell set info corresponding to this name 
	  
	  int j;
	  for (j = 0; j < 9; j++)
	    if (reg_name == expcsetnames[j]) break;
	  
	  CHECK(j < 9);
	  
	  // Verify that we can get the right number of entities in the set
	  
	  int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);
	  
	  CHECK_EQUAL(expcsetsizes[j],set_size);
	  
	  // Verify that we can get the correct set entities
	  
	  Amanzi::AmanziMesh::Entity_ID_List setents;
	  mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED,&setents);
	  
	  CHECK_ARRAY_EQUAL(expcsetcells[j],setents,set_size);
	}
    }
    else if (shape == "Region: Point") {

      // Do we have a valid cell set by this name
      
      CHECK(mesh.valid_set_name(reg_name,Amanzi::AmanziMesh::CELL));
      
      int j;
      for (j = 0; j < 9; j++) {
        if (expcsetnames[j] == reg_name) break;
      }
      
      CHECK(j < 9);
            
            
      // Verify that we can get the right number of entities in the set
      
      int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::USED);
      
      CHECK_EQUAL(expcsetsizes[j],set_size);
      
      
      // Verify that we can get the correct set entities
      
      Amanzi::AmanziMesh::Entity_ID_List setents;
      mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::USED,&setents);
      
      CHECK_ARRAY_EQUAL(expcsetcells[j],setents,set_size);	  
                        
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
        for (j = 0; j < 7; j++)
          if (reg_name == expfsetnames[j]) break;

        CHECK(j < 7);
	
	// Verify that we can get the right number of entities in the set
	
	int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED);
		
        CHECK_EQUAL(expfsetsizes[j],set_size);

	// Verify that we can get the correct set entities
	
        Amanzi::AmanziMesh::Entity_ID_List setents;
	mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED,&setents);

        CHECK_ARRAY_EQUAL(expfsetfaces[j],setents,set_size);

      }
      else if (entity_type == "Cell") {

	// Do we have a valid sideset by this name

	CHECK(mesh.valid_set_name(reg_name,Amanzi::AmanziMesh::CELL));
	
        // Find the expected face set info corresponding to this name

        int j;
        for (j = 0; j < 5; j++)
          if (reg_name == expcsetnames[j]) break;

        CHECK(j < 5);
	
	// Verify that we can get the right number of entities in the set
	
	int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);

        CHECK_EQUAL(expcsetsizes[j],set_size);
	
	// Verify that we can get the correct set entities
	
        Amanzi::AmanziMesh::Entity_ID_List setents;
	mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED,&setents);

        CHECK_ARRAY_EQUAL(expcsetcells[j],setents,set_size);
      }

    }
  }

}


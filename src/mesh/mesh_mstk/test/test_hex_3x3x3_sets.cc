#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"


#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"


TEST(MSTK_HEX_3x3x3_SETS)
{

  std::string expcsetnames[12] = {"Bottom LS", "Middle LS", "Top LS", 
                                 "Bottom+Middle Box", "Top Box",
                                 "Sample Point InCell", "Sample Point OnFace",
                                 "Sample Point OnEdge", "Sample Point OnVertex",
                                 "Bottom ColFunc", "Middle ColFunc", "Top ColFunc"};
  Amanzi::AmanziMesh::Set_ID csetsize;
  
  std::string expfsetnames[7] = {"Face 101", "Face 102", 
				  "Face 10005", "Face 20004", "Face 30004",
                                  "ZLO FACE Plane", "YLO FACE Box"};


  std::string expnsetnames[2] = {"INTERIOR XY PLANE", "TOP BOX"};

			   
  Teuchos::RCP<Epetra_MpiComm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));


  std::string infilename = "test/hex_3x3x3.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());

  Amanzi::AmanziGeometry::GeometricModelPtr gm = new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, comm.get());

  // Load a mesh consisting of 3x3x3 elements

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh(new Amanzi::AmanziMesh::Mesh_MSTK("test/hex_3x3x3_ss.exo",comm.get(),3,gm));


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
        
        CHECK(mesh->valid_set_name(reg_name,Amanzi::AmanziMesh::FACE));
        
        int j;
        for (j = 0; j < 7; j++) {
          if (expfsetnames[j] == reg_name) break;
        }
        
        CHECK(j < 7);
        
        
        // Verify that we can get the number of entities in the set
        
        int set_size = mesh->get_set_size(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED);
        
        
        // Verify that we can retrieve the set entities
        
        Amanzi::AmanziMesh::Entity_ID_List setents;
        mesh->get_set_entities(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED,&setents);
        
      }
      else if (reg_name == "INTERIOR XY PLANE") {

        // Do we have a valid nodeset by this name
        
        CHECK(mesh->valid_set_name(reg_name,Amanzi::AmanziMesh::NODE));
        
        int j;
        for (j = 0; j < 2; j++) {
          if (expnsetnames[j] == reg_name) break;
        }

        CHECK(j < 2);
        
        
        // Verify that we can get the number of entities in the set
        
        int set_size = mesh->get_set_size(reg_name,Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::USED);
        
        
        // Verify that we can retrieve the set entities
        
        Amanzi::AmanziMesh::Entity_ID_List setents;
        mesh->get_set_entities(reg_name,Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::USED,&setents);
        
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
            
            CHECK(mesh->valid_set_name(reg_name,Amanzi::AmanziMesh::FACE));
            
            int j;
            for (j = 0; j < 7; j++) {
              if (expfsetnames[j] == reg_name) break;
            }
            
            CHECK(j < 7);
            
            
            // Verify that we can get the number of entities in the set
            
            int set_size = mesh->get_set_size(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED);
            
            // Verify that we can retrieve the set entities
            
            Amanzi::AmanziMesh::Entity_ID_List setents;
            mesh->get_set_entities(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED,&setents);
            
          }
          else if (reg_name == "TOP BOX") {

            // Do we have a valid set by this name
            
            CHECK(mesh->valid_set_name(reg_name,Amanzi::AmanziMesh::NODE));
            
            int j;
            for (j = 0; j < 2; j++) {
              if (expnsetnames[j] == reg_name) break;
            }
            
            CHECK(j < 2);
            
            
            // Verify that we can get the number of entities in the set
            
            int set_size = mesh->get_set_size(reg_name,Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::USED);
            
            
            // Verify that we can retrieve the set entities
            
            Amanzi::AmanziMesh::Entity_ID_List setents;
            mesh->get_set_entities(reg_name,Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::USED,&setents);
            
          }

	}
      else 
	{
	  // Do we have a valid cellset by this name
	  
	  CHECK(mesh->valid_set_name(reg_name,Amanzi::AmanziMesh::CELL));
	  
	  // Find the expected cell set info corresponding to this name 
	  
	  int j;
	  for (j = 0; j < 9; j++)
	    if (reg_name == expcsetnames[j]) break;
	  
	  CHECK(j < 9);
	  
	  // Verify that we can get the number of entities in the set
	  
	  int set_size = mesh->get_set_size(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);
	  
	  // Verify that we can retrieve the set entities
	  
	  Amanzi::AmanziMesh::Entity_ID_List setents;
	  mesh->get_set_entities(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED,&setents);
	}
    }
    else if (shape == "Region: Point") {

      // Do we have a valid cell set by this name
      
      CHECK(mesh->valid_set_name(reg_name,Amanzi::AmanziMesh::CELL));
      
      int j;
      for (j = 0; j < 9; j++) {
        if (expcsetnames[j] == reg_name) break;
      }
      
      CHECK(j < 9);
            
            
      // Verify that we can get the number of entities in the set
      
      int set_size = mesh->get_set_size(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::USED);
      
      
      // Verify that we can retrieve the set entities
      
      Amanzi::AmanziMesh::Entity_ID_List setents;
      mesh->get_set_entities(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::USED,&setents);
      
    }
    else if (shape == "Region: Labeled Set") {

      Teuchos::ParameterList lsparams = reg_params.sublist(shape);

      // Find the entity type in this parameter list

      std::string entity_type = lsparams.get<std::string>("Entity");

      if (entity_type == "Face") {

	// Do we have a valid sideset by this name

	CHECK(mesh->valid_set_name(reg_name,Amanzi::AmanziMesh::FACE));

        // Find the expected face set info corresponding to this name

        int j;
        for (j = 0; j < 7; j++)
          if (reg_name == expfsetnames[j]) break;

        CHECK(j < 7);
	
	// Verify that we can get the number of entities in the set
	
	int set_size = mesh->get_set_size(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED);
		

	// Verify that we can retrieve the set entities
	
        Amanzi::AmanziMesh::Entity_ID_List setents;
	mesh->get_set_entities(reg_name,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED,&setents);

      }
      else if (entity_type == "Cell") {

	// Do we have a valid sideset by this name

	CHECK(mesh->valid_set_name(reg_name,Amanzi::AmanziMesh::CELL));
	
        // Find the expected face set info corresponding to this name

        int j;
        for (j = 0; j < 5; j++)
          if (reg_name == expcsetnames[j]) break;

        CHECK(j < 5);
	
	// Verify that we can get the number of entities in the set
	
	int set_size = mesh->get_set_size(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);

	// Verify that we can retrieve the set entities
	
        Amanzi::AmanziMesh::Entity_ID_List setents;
	mesh->get_set_entities(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED,&setents);
      }

    }
    else if (shape == "Region: Color Function") {

      // Do we have a valid sideset by this name
      
      CHECK(mesh->valid_set_name(reg_name,Amanzi::AmanziMesh::CELL));
      
      // Find the expected face set info corresponding to this name
      
      int j;
      for (j = 0; j < 12; j++)
        if (reg_name == expcsetnames[j]) break;
      
      CHECK(j < 12);
      
      // Verify that we can get the number of entities in the set
      
      int set_size = mesh->get_set_size(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);
      
      // Verify that we can retrieve the set entities
      
      Amanzi::AmanziMesh::Entity_ID_List setents;
      mesh->get_set_entities(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED,&setents);
      
    }
  }


  // Once we can make RegionFactory work with reference counted pointers 
  // we can get rid of this code

  for (int i = 0; i < gm->Num_Regions(); i++)
    delete (gm->Region_i(i));
  delete gm;

}


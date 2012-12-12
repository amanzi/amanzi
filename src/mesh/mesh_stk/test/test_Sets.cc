// -------------------------------------------------------------
/**
 * @file   test_Sets.cc
 * @author Rao Garimella 
 * @date 
 * 
 * @brief Some unit tests for generating a Mesh_STK instance,
 * creating the right sets from an input specification, and
 * responding correctly to set queries
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 22, 2010 by William A. Perkins
// Last Change: Wed May 18 14:15:54 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <sstream>
#include <UnitTest++.h>
#include <Epetra_MpiComm.h>

#include "Exodus_readers.hh"
#include "Parallel_Exodus_file.hh"
#include "../Mesh_STK_Impl.hh"
#include "../Mesh_STK.hh"
#include "Auditor.hh"

#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "mpi.h"


SUITE (STK_SETS)
{
  TEST (SETS_READ)
  {

    std::string expcsetnames[12] = {"Bottom LS", "Middle LS", "Top LS", 
                                    "Bottom+Middle Box", "Top Box",
                                    "Sample Point InCell", "Sample Point OnFace",
                                    "Sample Point OnEdge", "Sample Point OnVertex",
                                    "Bottom ColFunc", "Middle ColFunc", "Top ColFunc"};  
  
    unsigned int csetsize, expcsetsizes[12] = {9,9,9,18,9,1,2,4,8,9,9,9};
  
    int expcsetcells[12][18] = {{0,1,2,3,4,5,6,7,8,-1,-1,-1,-1,-1,-1,-1,-1},
				{9,10,11,12,13,14,15,16,17,-1,-1,-1,-1,-1,-1,-1,-1,-1},
				{18,19,20,21,22,23,24,25,26,-1,-1,-1,-1,-1,-1,-1,-1,-1},
				{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17},
				{18,19,20,21,22,23,24,25,26,-1,-1,-1,-1,-1,-1,-1,-1,-1},
				{13,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
				{12,13,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
				{9,10,12,13,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
				{0,1,3,4,9,10,12,13,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
				{0,1,2,3,4,5,6,7,8,-1,-1,-1,-1,-1,-1,-1,-1},
				{9,10,11,12,13,14,15,16,17,-1,-1,-1,-1,-1,-1,-1,-1,-1},
				{18,19,20,21,22,23,24,25,26,-1,-1,-1,-1,-1,-1,-1,-1,-1}};
    
    std::string expfsetnames[7] = {"Face 101", "Face 102", 
                                   "Face 10005", "Face 20004", "Face 30004",
                                   "ZLO FACE Plane", "YLO FACE Box"};

    unsigned int expfsetids[7]={101,102,10005,20004,30004,0,0};
  
    unsigned int fsetsize, expfsetsizes[7] = {9,9,3,3,3,9,9};


    int expfsetfaces[7][9] = {{4,9,14,19,23,27,32,36,40},
			      {0,6,11,42,47,51,75,80,84},
			      {30,35,39,-1,-1,-1,-1,-1,-1},
			      {66,70,73,-1,-1,-1,-1,-1,-1},
			      {99,103,106,-1,-1,-1,-1,-1,-1},
			      {4,9,14,19,23,27,32,36,40},
			      {0,6,11,42,47,51,75,80,84}};
    

    std::string infilename = "test/hex_3x3x3_read.xml";
    Teuchos::ParameterXMLFileReader xmlreader(infilename);

    Teuchos::ParameterList reg_spec(xmlreader.getParameters());

    Epetra_MpiComm *comm(new Epetra_MpiComm(MPI_COMM_WORLD));

    Amanzi::AmanziGeometry::GeometricModelPtr gm = new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, comm);

    // Load a mesh consisting of 3x3x3 elements

    Amanzi::AmanziMesh::Mesh_STK mesh("test/hex_3x3x3_ss.exo",comm,gm);


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
        else 
          {
            // Do we have a valid cellset by this name
	  
            CHECK(mesh.valid_set_name(reg_name,Amanzi::AmanziMesh::CELL));
	  
            // Find the expected cell set info corresponding to this name 
	  
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
      else if (shape == "Region: Point") {

        // Get cells around this point

        Teuchos::ParameterList point_params = reg_params.sublist(shape);
        Teuchos::Array<double> p_vec = point_params.get< Teuchos::Array<double> >("Coordinate");

        // Do we have a valid set by this name
      
        CHECK(mesh.valid_set_name(reg_name,Amanzi::AmanziMesh::CELL));
	  
        int j;
        for (j = 0; j < 12; j++) {
          if (expcsetnames[j] == reg_name) break;
        }
	  
        CHECK(j < 12);
	  
	  
        // Verify that we can get the right number of entities in the set
	  
        int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);
	  
        CHECK_EQUAL(expcsetsizes[j],set_size);
	  
	  
        // Verify that we can get the correct set entities
      
        Amanzi::AmanziMesh::Entity_ID_List setents;
        mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED,&setents);
	  
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
          for (j = 0; j < 12; j++)
            if (reg_name == expcsetnames[j]) break;

          CHECK(j < 12);
	
          // Verify that we can get the right number of entities in the set
	
          int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);

          CHECK_EQUAL(expcsetsizes[j],set_size);
	
          // Verify that we can get the correct set entities
	
          Amanzi::AmanziMesh::Entity_ID_List setents;
          mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED,&setents);

          CHECK_ARRAY_EQUAL(expcsetcells[j],setents,set_size);
        }

      }
      else if (shape == "Region: Color Function") {

        // Do we have a valid sideset by this name
      
        CHECK(mesh.valid_set_name(reg_name,Amanzi::AmanziMesh::CELL));
      
        // Find the expected face set info corresponding to this name
      
        int j;
        for (j = 0; j < 12; j++)
          if (reg_name == expcsetnames[j]) break;
      
        CHECK(j < 12);
	
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



  TEST (SETS_GEN)
  {
    using namespace std;

#ifdef HAVE_MPI
    Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

    std::string expcsetnames[4] = {"Bottom Box", "Bottom+Middle Box",
				   "Vertical Box", "Sample Point 1"};
    unsigned int csetsize, expcsetsizes[4] = {9,18,9,8};
  
    int expcsetcells[4][18] = {{0,1,2,3,4,5,6,7,8,-1,-1,-1,-1,-1,-1,-1,-1,-1},
			       {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17},
			       {1,4,7,10,13,16,19,22,25,-1,-1,-1,-1,-1,-1,-1,-1,-1},
			       {0,1,3,4,9,10,12,13,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}
    };

    std::string expfsetnames[2] = {"ZLO FACE Plane", "YLO FACE Box"};

    unsigned int fsetsize, expfsetsizes[2] = {9,3};


    int expfsetfaces[6][9] = {{3,8,13,19,23,27,32,36,40},
			      {78,82,86,-1,-1,-1,-1,-1,-1}};


    std::string infilename = "test/hex_3x3x3_gen.xml";
    Teuchos::ParameterXMLFileReader xmlreader(infilename);

    Teuchos::ParameterList reg_spec(xmlreader.getParameters());

    Epetra_MpiComm ecomm(MPI_COMM_WORLD);

    Amanzi::AmanziGeometry::GeometricModelPtr gm = new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, &ecomm);

    // Create a mesh consisting of 3x3x3 elements

    Amanzi::AmanziMesh::Mesh_STK mesh(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3, comm, gm); 

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
      else if (shape == "Region: Point") {
        
        Teuchos::ParameterList point_params = reg_params.sublist(shape);
        Teuchos::Array<double> p_vec = point_params.get< Teuchos::Array<double> >("Coordinate");
        
        // Do we have a valid set by this name
        
        CHECK(mesh.valid_set_name(reg_name,Amanzi::AmanziMesh::CELL));
        
        int j;
        for (j = 0; j < 4; j++) {
          if (expcsetnames[j] == reg_name) break;
        }
        
        CHECK(j < 4);
        
	
        // Verify that we can get the right number of entities in the set
        
        int set_size = mesh.get_set_size(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);
        
        CHECK_EQUAL(expcsetsizes[j],set_size);
        
	
        // Verify that we can get the correct set entities
        
        Amanzi::AmanziMesh::Entity_ID_List setents;
        mesh.get_set_entities(reg_name,Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED,&setents);
        
        CHECK_ARRAY_EQUAL(expcsetcells[j],setents,set_size);	  
        
      }
      else if (shape == "Region: Labeled Set") {

	std::cerr << "Mesh framework cannot do labeled sets" << std::endl;

	CHECK(false); 

      }
    }

  }
}





#include "UnitTest++.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "State.hpp"
#include "Mesh_STK.hh"
#include "Restart.hpp"


SUITE(RESTART) {

  TEST(RESTART_DUMP_REQUIRED_INTERVAL) {
    
    Teuchos::ParameterList plist;

    plist.set<string>("File Name Base","restartdump");

    Teuchos::ParameterList& i3_ = plist.sublist("Cycle Data");
    
    i3_.set<int>("Start",0);
    i3_.set<int>("End",10);
    i3_.set<int>("Interval",3);
    

    Epetra_MpiComm comm(MPI_COMM_WORLD);
    Amanzi::Restart R(plist, &comm);

    
    // test the cycle stuff, the expected result is in cycles_ and 
    // we store the computed result in cycles
    
    int cycles_[31] = { 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
                        0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  };
    int cycles [31];
    for (int ic = 0; ic<=30; ic++)
      {
  	cycles[ic] = R.dump_requested(ic);
      }
    CHECK_ARRAY_EQUAL(cycles_, cycles, 31);

  }

  TEST(RESTART_DUMP_REQUIRED_INTERVAL_OPENENDED1) {
    
    Teuchos::ParameterList plist;

    plist.set<string>("File Name Base","restartdump");

    Teuchos::ParameterList& i3_ = plist.sublist("Cycle Data");
    
    i3_.set<int>("Start",0);
    i3_.set<int>("End",-1);
    i3_.set<int>("Interval",3);
    

    Epetra_MpiComm comm(MPI_COMM_WORLD);
    Amanzi::Restart R(plist, &comm);

    
    // test the cycle stuff, the expected result is in cycles_ and 
    // we store the computed result in cycles
    
    int cycles_[31] = { 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 
                        1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1  };
    int cycles [31];
    for (int ic = 0; ic<=30; ic++)
      {
  	cycles[ic] = R.dump_requested(ic);
      }
    CHECK_ARRAY_EQUAL(cycles_, cycles, 31);

  }


  TEST(RESTART_DUMP_REQUIRED_INTERVAL_OPENENDED2) {
    
    Teuchos::ParameterList plist;

    plist.set<string>("File Name Base","restartdump");

    Teuchos::ParameterList& i3_ = plist.sublist("Cycle Data");
    
    i3_.set<int>("Start",5);
    i3_.set<int>("End",-1);
    i3_.set<int>("Interval",3);
    

    Epetra_MpiComm comm(MPI_COMM_WORLD);
    Amanzi::Restart R(plist, &comm);

    
    // test the cycle stuff, the expected result is in cycles_ and 
    // we store the computed result in cycles
    
    int cycles_[31] = { 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 
                        0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0  };
    int cycles [31];
    for (int ic = 0; ic<=30; ic++)
      {
  	cycles[ic] = R.dump_requested(ic);
      }
    CHECK_ARRAY_EQUAL(cycles_, cycles, 31);

  }



  TEST(RESTART_DUMP_REQUIRED_STEPS) {
    
    Teuchos::ParameterList plist;

    plist.set<string>("File Name Base","restartdump");

    Teuchos::ParameterList& i3_ = plist.sublist("Cycle Data");
    
    i3_.set<int>("Start",0);
    i3_.set<int>("End",10);
    i3_.set<int>("Interval",3);
    
    Teuchos::Array<int> steps(3);
    steps[0] = 2;
    steps[1] = 3;
    steps[2] = 4;
    i3_.set<Teuchos::Array<int> >("Steps",steps);
    

    Epetra_MpiComm comm(MPI_COMM_WORLD);
    Amanzi::Restart R(plist, &comm);

    
    // test the cycle stuff, the expected result is in cycles_ and 
    // we store the computed result in cycles
    
    int cycles_[31] = { 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  };
    int cycles [31];
    for (int ic = 0; ic<=30; ic++)
      {
  	cycles[ic] = R.dump_requested(ic);
      }
    CHECK_ARRAY_EQUAL(cycles_, cycles, 31);

  }




  TEST(DUMP_DATA) 
  {
    Teuchos::ParameterList plist;
    plist.set<string>("File Name Base","restart_dump");
    plist.set<int>("File Name Digits",4);
    Teuchos::ParameterList& i1_ = plist.sublist("Cycle Data");
    
    i1_.set<int>("Start",0);
    i1_.set<int>("End",10);
    i1_.set<int>("Interval",1);

    Epetra_MpiComm comm(MPI_COMM_WORLD);
    Amanzi::Restart R(plist, &comm);   

    // make a simple mesh
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh_STK> Mesh
      = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_STK(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 1, 1, &comm));

    // create a state object with some data in it
    int number_of_components = 2;
    int number_of_minerals = 2;
    State S0(number_of_components, number_of_minerals, Mesh);

    S0.set_time(1.02);

    // make sure that the parameter list is such that we will actually dump something (see above)
    S0.set_cycle(4);

    double g[3] = { 1.0, 2.0, 3.0 }; // totally bogus gravity
    S0.set_gravity(g);

    S0.set_water_density(99.9);
    S0.set_viscosity(0.22);

    Epetra_Vector* face_vector = new Epetra_Vector(S0.get_mesh().face_epetra_map(false));
    face_vector->Random();
    S0.set_darcy_flux(*face_vector);
    delete face_vector;

    Epetra_Vector* cell_vector = new Epetra_Vector(S0.get_mesh().cell_epetra_map(false));
    cell_vector->Random();
    S0.set_water_saturation(*cell_vector);
    delete cell_vector;

    cell_vector = new Epetra_Vector(S0.get_mesh().cell_epetra_map(false));
    cell_vector->Random();
    S0.set_water_density(*cell_vector);
    delete cell_vector;
    
    cell_vector = new Epetra_Vector(S0.get_mesh().cell_epetra_map(false));
    cell_vector->Random();
    S0.set_pressure(*cell_vector);
    delete cell_vector;

    face_vector = new Epetra_Vector(S0.get_mesh().face_epetra_map(false));
    face_vector->Random();
    S0.set_lambda(*face_vector);
    delete face_vector;

    cell_vector = new Epetra_Vector(S0.get_mesh().cell_epetra_map(false));
    cell_vector->Random();
    S0.set_porosity(*cell_vector);
    delete cell_vector;    
    
    cell_vector = new Epetra_Vector(S0.get_mesh().cell_epetra_map(false));
    cell_vector->Random();
    S0.set_permeability(*cell_vector);
    delete cell_vector;    
    
    Epetra_MultiVector* cell_multivector = new Epetra_MultiVector(S0.get_mesh().cell_epetra_map(false), number_of_components);
    cell_multivector->Random();
    S0.set_total_component_concentration(*cell_multivector);
    delete cell_multivector;       

    std::vector<std::string> mineral_names(number_of_minerals);
    mineral_names.at(0) = "Aoeui";
    mineral_names.at(1) = "Snthd";
    S0.set_mineral_names(mineral_names);
 
    cell_multivector = new Epetra_MultiVector(S0.get_mesh().cell_epetra_map(false),
                                              S0.number_of_minerals());
    cell_multivector->Random();
    S0.set_mineral_volume_fractions(*cell_multivector);
    cell_multivector->Random();
    S0.set_mineral_specific_surface_area(*cell_multivector);
    delete cell_multivector;

    R.dump_state(S0);

    // now read the file into a new state object

    State S1(number_of_components, number_of_minerals, Mesh);

    std::string filename = "restart_dump0004.h5";

    R.read_state(S1, filename);

    // and compare with the original

    CHECK_EQUAL(S0.get_time(), S1.get_time());
    CHECK_EQUAL(S0.get_cycle(), S1.get_cycle());
    CHECK_EQUAL((*S0.get_gravity())[0],(*S1.get_gravity())[0]);
    CHECK_EQUAL((*S0.get_gravity())[1],(*S1.get_gravity())[1]);
    CHECK_EQUAL((*S0.get_gravity())[2],(*S1.get_gravity())[2]);
    CHECK_EQUAL(*S0.get_density(),*S1.get_density());
    CHECK_EQUAL(*S0.get_viscosity(),*S1.get_viscosity());
    CHECK_EQUAL(S0.get_number_of_components(),S1.get_number_of_components());
    
    int s0_size = S0.get_darcy_flux()->MyLength();
    int s1_size = S1.get_darcy_flux()->MyLength();
    CHECK_EQUAL(s0_size, s1_size);

    for (int i=0; i<s0_size; i++) 
      {
  	CHECK_EQUAL((*S0.get_darcy_flux())[i],(*S1.get_darcy_flux())[i]);
      }


    s0_size = S0.get_water_saturation()->MyLength();
    s1_size = S1.get_water_saturation()->MyLength();
    CHECK_EQUAL(s0_size, s1_size);

    for (int i=0; i<s0_size; i++) 
      {
  	CHECK_EQUAL((*S0.get_water_saturation())[i],(*S1.get_water_saturation())[i]);
      }    

    
    s0_size = S0.get_water_density()->MyLength();
    s1_size = S1.get_water_density()->MyLength();
    CHECK_EQUAL(s0_size, s1_size);

    for (int i=0; i<s0_size; i++) 
      {
  	CHECK_EQUAL((*S0.get_water_density())[i],(*S1.get_water_density())[i]);
      }    
    
    s0_size = S0.get_pressure()->MyLength();
    s1_size = S1.get_pressure()->MyLength();
    CHECK_EQUAL(s0_size, s1_size);

    for (int i=0; i<s0_size; i++) 
      {
  	CHECK_EQUAL((*S0.get_pressure())[i],(*S1.get_pressure())[i]);
      }    

    s0_size = S0.get_lambda()->MyLength();
    s1_size = S1.get_lambda()->MyLength();
    CHECK_EQUAL(s0_size, s1_size);

    for (int i=0; i<s0_size; i++) 
      {
  	CHECK_EQUAL((*S0.get_lambda())[i],(*S1.get_lambda())[i]);
      }    

    s0_size = S0.get_porosity()->MyLength();
    s1_size = S1.get_porosity()->MyLength();
    CHECK_EQUAL(s0_size, s1_size);

    for (int i=0; i<s0_size; i++) 
      {
  	CHECK_EQUAL((*S0.get_porosity())[i],(*S1.get_porosity())[i]);
      }    

    s0_size = S0.get_vertical_permeability()->MyLength();
    s1_size = S1.get_vertical_permeability()->MyLength();
    CHECK_EQUAL(s0_size, s1_size);

    s0_size = S0.get_horizontal_permeability()->MyLength();
    s1_size = S1.get_horizontal_permeability()->MyLength();
    CHECK_EQUAL(s0_size, s1_size);

    for (int i=0; i<s0_size; i++) 
      {
  	CHECK_EQUAL((*S0.get_vertical_permeability())[i],(*S1.get_vertical_permeability())[i]);
      }    

    for (int i=0; i<s0_size; i++) 
      {
  	CHECK_EQUAL((*S0.get_horizontal_permeability())[i],(*S1.get_horizontal_permeability())[i]);
      }    


    s0_size = S0.get_total_component_concentration()->MyLength();
    s1_size = S1.get_total_component_concentration()->MyLength();
    CHECK_EQUAL(s0_size, s1_size);

    for (int i=0; i<s0_size; i++) 
      {
  	for (int n=0; n<S0.get_number_of_components(); n++)
  	  {
  	    CHECK_EQUAL( (*(*S0.get_total_component_concentration())(n))[i], 
  			 (*(*S1.get_total_component_concentration())(n))[i]);
  	  }
      }    


    CHECK_EQUAL(S0.number_of_minerals(), S1.number_of_minerals());
    CHECK_EQUAL(S0.mineral_names().size(), S1.mineral_names().size());
    CHECK_EQUAL(S1.number_of_minerals(), S1.mineral_names().size());
    for (int m = 0; m < S0.mineral_names().size(); ++m) {
      CHECK_EQUAL(S0.mineral_names().at(m), S1.mineral_names().at(m));
    }
    
    int num_cells = S0.mineral_volume_fractions()->MyLength();
    int num_cells_2 = S1.mineral_volume_fractions()->MyLength();
    CHECK_EQUAL(num_cells, num_cells_2);
    for (int m = 0; m < S0.number_of_minerals(); ++m) {
      for (int cell = 0; cell < num_cells; ++cell) {
        double v1 = (*(*S0.mineral_volume_fractions())(m))[cell];
        double v2 = (*(*S1.mineral_volume_fractions())(m))[cell];
        CHECK_EQUAL(v1, v2);
      }
    }
    num_cells = S0.mineral_specific_surface_area()->MyLength();
    num_cells_2 = S1.mineral_specific_surface_area()->MyLength();
    CHECK_EQUAL(num_cells, num_cells_2);
    for (int m = 0; m < S0.number_of_minerals(); ++m) {
      for (int cell = 0; cell < num_cells; ++cell) {
        double v1 = (*(*S0.mineral_specific_surface_area())(m))[cell];
        double v2 = (*(*S1.mineral_specific_surface_area())(m))[cell];
        CHECK_EQUAL(v1, v2);
      }
    }
  }
  
}

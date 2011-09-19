#include "UnitTest++.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "State.hpp"
#include "Mesh_STK.hh"
#include "Restart.hpp"


SUITE(RESTART) {

  TEST(RESTART_DUMP_REQUIRED) {
    
    Teuchos::ParameterList plist;

    plist.set<string>("file base name","restartdump");

    Teuchos::ParameterList& i3 = plist.sublist("Interval 1");
    Teuchos::ParameterList& i3_ = i3.sublist("cycle range");
    
    i3_.set<int>("start cycle",0);
    i3_.set<int>("end cycle",10);
    i3_.set<int>("cycle frequency",3);

    Teuchos::ParameterList& i4 = plist.sublist("Interval 2");
    Teuchos::ParameterList& i4_ = i4.sublist("cycle range");
    
    i4_.set<int>("start cycle",15);
    i4_.set<int>("end cycle",30);
    i4_.set<int>("cycle frequency",10);
   
    
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    Amanzi::Restart R(plist, &comm);

    
    // test the cycle stuff, the expected result is in cycles_ and 
    // we store the computed result in cycles
    
    int cycles_[31] = { 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
                        0, 0, 0, 0,
                        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0  };
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
    plist.set<string>("file base name","restart_dump");
    Teuchos::ParameterList& i1 = plist.sublist("Interval 1");
    Teuchos::ParameterList& i1_ = i1.sublist("cycle range");
    
    i1_.set<int>("start cycle",0);
    i1_.set<int>("end cycle",10);
    i1_.set<int>("cycle frequency",1);

    Epetra_MpiComm comm(MPI_COMM_WORLD);
    Amanzi::Restart R(plist, &comm);   

    // make a simple mesh
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh_STK> Mesh
      = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_STK(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 1, 1, &comm));

    R.create_files(); 

    // create a state object with some data in it
    int number_of_components = 2;
    State S0(number_of_components, Mesh);

    S0.set_time(1.02);

    // make sure that the parameter list is such that we will actually dump something (see above)
    S0.set_cycle(4);

    double g[3] = { 1.0, 2.0, 3.0 }; // totally bogus gravity
    S0.set_gravity(g);

    S0.set_water_density(99.9);
    S0.set_viscosity(0.22);

    Epetra_Vector face_vector(S0.get_mesh().face_epetra_map(false));
    face_vector.Random();
    S0.set_darcy_flux(face_vector);

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

    cell_vector = new Epetra_Vector(S0.get_mesh().cell_epetra_map(false));
    cell_vector->Random();
    S0.set_porosity(*cell_vector);
    delete cell_vector;    
    
    cell_vector = new Epetra_Vector(S0.get_mesh().cell_epetra_map(false));
    cell_vector->Random();
    S0.set_permeability(*cell_vector);
    delete cell_vector;    
    
    Epetra_MultiVector* cell_multivector = new Epetra_MultiVector(S0.get_mesh().cell_epetra_map(false), 3);
    cell_multivector->Random();
    S0.set_darcy_velocity(*cell_multivector);
    delete cell_multivector;        

    cell_multivector = new Epetra_MultiVector(S0.get_mesh().cell_epetra_map(false), number_of_components);
    cell_multivector->Random();
    S0.set_total_component_concentration(*cell_multivector);
    delete cell_multivector;       



    R.dump_state(S0);


    // now read the file into a new state object

    State S1(Mesh);
    R.read_state(S1);

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

    s0_size = S0.get_porosity()->MyLength();
    s1_size = S1.get_porosity()->MyLength();
    CHECK_EQUAL(s0_size, s1_size);

    for (int i=0; i<s0_size; i++) 
      {
	CHECK_EQUAL((*S0.get_porosity())[i],(*S1.get_porosity())[i]);
      }    

    s0_size = S0.get_permeability()->MyLength();
    s1_size = S1.get_permeability()->MyLength();
    CHECK_EQUAL(s0_size, s1_size);

    for (int i=0; i<s0_size; i++) 
      {
	CHECK_EQUAL((*S0.get_permeability())[i],(*S1.get_permeability())[i]);
      }    

    s0_size = S0.get_darcy_velocity()->MyLength();
    s1_size = S1.get_darcy_velocity()->MyLength();
    CHECK_EQUAL(s0_size, s1_size);

    for (int i=0; i<s0_size; i++) 
      {
	CHECK_EQUAL( (*(*S0.get_darcy_velocity())(0))[i], (*(*S1.get_darcy_velocity())(0))[i]);
	CHECK_EQUAL( (*(*S0.get_darcy_velocity())(1))[i], (*(*S1.get_darcy_velocity())(1))[i]);
	CHECK_EQUAL( (*(*S0.get_darcy_velocity())(2))[i], (*(*S1.get_darcy_velocity())(2))[i]);

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

  }



}

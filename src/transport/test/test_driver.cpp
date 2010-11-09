#include <iostream>
#include <cstdlib>
#include <cmath>
#include "UnitTest++.h"
#include <vector>
#include <cstring>
#include <string>

#include "moab_mesh/Mesh_maps_moab.hh"
#include "simple_mesh/Mesh_maps_simple.hh"
#include "mpc/State.hpp"
#include "mpc/MPC.hpp"
#include "transport/Transport_PK.hpp"
#include "transport/Transport_State.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"


TEST(DRIVER) {

  using namespace Teuchos;
  using namespace std;

#ifdef HAVE_MPI
  Epetra_MpiComm     *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm  *comm = new Epetra_SerialComm();
#endif
  
  string  xmlInFileName = "test/test_driver.xml";


  /* read the main parameter list */
  ParameterList  driver_list;
  updateParametersFromXmlFile( xmlInFileName, &driver_list );
  
  /* get the Mesh sublist */
  ParameterList  mesh_list = driver_list.sublist( "Mesh" );

  string  mesh_class = mesh_list.get<string>( "Mesh Class" );

  if ( mesh_class == "Simple" ) {
cout << "Simple mesh" << endl;
     RCP<Mesh_maps_simple>  MMS;
     ParameterList  simple_mesh_list = mesh_list.sublist( "Simple Mesh Parameters" );

     MMS = rcp( new Mesh_maps_simple( simple_mesh_list, comm ) );

     MPC  mpc( driver_list, MMS );
     mpc.cycle_driver();
  } 
  else if ( mesh_class == "MOAB" ) {
cout << "MOAB mesh" << endl;
     RCP<Mesh_maps_moab>  MMB;
     ParameterList  moab_mesh_list = mesh_list.sublist( "MOAB Mesh Parameters" );

     string  file_name = moab_mesh_list.get<string>( "Mesh file name" );

     MMB = rcp( new Mesh_maps_moab( file_name.c_str(), MPI_COMM_WORLD ) );

     MPC  mpc( driver_list, MMB );
     mpc.cycle_driver();
  }

  delete comm;
}



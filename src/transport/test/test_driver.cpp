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
  Epetra_MpiComm     *comm = new Epetra_MpiComm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm  *comm = new Epetra_SerialComm();
#endif
  
  /* read the main parameter list */
  ParameterList  driver_list;
  string  xmlFileName = "test/test_driver.xml";

  updateParametersFromXmlFile( xmlFileName, &driver_list );
  
  /* get the Mesh sublist */
  ParameterList  mesh_list = driver_list.sublist( "Mesh" );

  string  mesh_class = mesh_list.get<string>( "Mesh Class" );

  RCP<Mesh_maps_base>  mesh;

  if ( mesh_class == "Simple" ) {
     ParameterList  simple_mesh_list = mesh_list.sublist( "Simple Mesh Parameters" );

     RCP<Mesh_maps_simple> MMS = rcp( new Mesh_maps_simple( simple_mesh_list, comm ) );

     mesh = MMS;
  } 
  else if ( mesh_class == "MOAB" ) {
     ParameterList  moab_mesh_list = mesh_list.sublist( "MOAB Mesh Parameters" );

     string  file_name = moab_mesh_list.get<string>( "Mesh file name" );

     RCP<Mesh_maps_moab>  MMM = rcp( new Mesh_maps_moab( file_name.c_str(), MPI_COMM_WORLD ) );

     mesh = MMM;
  }

  MPC  mpc( driver_list, mesh );
  mpc.cycle_driver();

 
  /* print influx data */
  /*
  int  i, n;
  for( n=0; n<TPK.bcs.size(); n++ ) {
     cout << "Side set " << n << " of type " << TPK.bcs[n].type << endl;
     cout << "   influx: ";
     for( i=0; i<number_components; i++ ) cout << TPK.bcs[n].influx[i] << " ";
     cout << endl;
  }
  */

  delete comm;
}



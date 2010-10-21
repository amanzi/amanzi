#include "Epetra_MultiVector.h"

#include "../simple_mesh/Mesh_maps_simple.hh"
#include "Transport_PK.hpp"


using namespace Teuchos;


/* 
   Constructor for initializing the transport PK.
   Its call is made usually at time T=0 by the MPC.
 */
Transport_PK::Transport_PK ( RCP<Transport_State> TS_MPC )
{ 
  /* make a copy of the transport state object */
  TS = TS_MPC;

  /* copy pointers for state variables that will remain unchanged */
  TS_next = rcp(new Transport_State());
  TS_next->get_porosity () = TS->get_porosity (); 
  TS_next->get_water_saturation () = TS->get_water_saturation (); 
  TS_next->get_darcy_flux () = TS->get_darcy_flux (); 

  /* allocate memory for state variables that will be changed */
  RCP<Epetra_MultiVector>  tcc      = TS->get_total_component_concentration ();
  RCP<Epetra_MultiVector>  tcc_next = TS_next->get_total_component_concentration ();
  tcc_next = rcp( new Epetra_MultiVector( *tcc ) );

  /* set null/zero values to all internal parameters */
  dT = 0.0;
  status = 0;
};



/* null destructor */
Transport_PK::~Transport_PK () {};



void Transport_PK::advance_transport_state ()
{
  cout << "advancing the state of the transport process model here" << endl;

  // the MPC will call this function to advance the state 
  // with this particular process kernel

  // this is how to get the total component concentration
  // TS->get_total_component_concentration()

  // see the Transport_State for the other available 
  // data in the transport specific state

};




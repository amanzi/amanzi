#include "Epetra_MultiVector.h"

#include "../simple_mesh/Mesh_maps_simple.hh"
#include "Transport_PK.hpp"


/* constructor initializes the transport PK */
Transport_PK::Transport_PK ( Teuchos::RCP<Transport_State> TS_MPC )
{ 
  // use the constructor to initialize the transport process kernel

  // make a deep copy of the transport state object

  TS->copy(TS_MPC);

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




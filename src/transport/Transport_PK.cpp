#include "Transport_PK.hpp"
#include "Epetra_MultiVector.h"


/* constructor initializes the transport PK */
Transport_PK::Transport_PK ( Teuchos::RCP<Transport_State> TS_MPC )
{ 
  TS = TS_MPC;
  
  status = NULL;
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




#include "Transport_PK.hpp"
#include "Epetra_MultiVector.h"


Transport_PK::Transport_PK( Teuchos::RCP<Transport_State> TS_):
  TS(TS_)
{ 
  // use the constructor to initialize the transport process kernel

};

Transport_PK::~Transport_PK()
{ 

};


void Transport_PK::advance( )
{
  cout << "advancing the state of the transport process model here" << endl;

  // the MPC will call this function to advance the state 
  // with this particular process kernel

  // this is how to get the total component concentration
  // TS->get_total_component_concentration()

  // see the Transport_State for the other available 
  // data in the transport specific state

};


void Transport_PK::commit_state( Teuchos::RCP<Transport_State> )
{
  cout << "committing the internal state of the transport process model" << endl;

  // the MPC will call this function to signal to the 
  // process kernel that it has accepted the 
  // state update, thus, the PK should update
  // possible auxilary state variables here 

  // use this function to commit the internal state
  // of the transport process kernel

};



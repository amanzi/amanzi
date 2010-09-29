#include "Flow_PK.hpp"
#include "Epetra_MultiVector.h"


Flow_PK::Flow_PK( Teuchos::RCP<Flow_State> CS_):
  CS(CS_)
{ 
  // use the constructor to initialize the flow process kernel

};

Flow_PK::~Flow_PK()
{ 

};


void Flow_PK::advance( Teuchos::RCP<Epetra_MultiVector> tcc_star )
{
  cout << "advancing the state of the flow process model here" << endl;

  // this is how to get the element volumes...
  // mesh_wrapper->get_element_volumes() 

  // this is how to get the total component concentration
  // CS->get_total_component_concentration()

  // please update the argument to this function called tcc_star
  // with the result of your flow computation which is
  // the total component concentration ^star

  // see the Flow_State for the other available 
  // data in the flow specific state

};


void Flow_PK::commit_state( Teuchos::RCP<Flow_State> )
{
  cout << "committing the internal state of the flow process model" << endl;

  // use this function to commit the internal state
  // of the flow process kernel

};



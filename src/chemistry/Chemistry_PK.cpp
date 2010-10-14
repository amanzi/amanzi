#include "Chemistry_PK.hpp"
#include "Epetra_MultiVector.h"


Chemistry_PK::Chemistry_PK( Teuchos::RCP<Chemistry_State> CS_):
  CS(CS_)
{ 
  // use the constructor to initialize the chemistry process kernel

};

Chemistry_PK::~Chemistry_PK()
{ 

};


void Chemistry_PK::advance( Teuchos::RCP<Epetra_MultiVector> tcc_star )
{
  cout << "advancing the state of the chemistry process model here" << endl;

  // this is how to get the element volumes...
  // mesh_wrapper->get_element_volumes() 

  // this is how to get the total component concentration
  // CS->get_total_component_concentration()

  // please update the argument to this function called tcc_star
  // with the result of your chemistry computation which is
  // the total component concentration ^star

  // see the Chemistry_State for the other available 
  // data in the chemistry specific state

};


void Chemistry_PK::commit_state( Teuchos::RCP<Chemistry_State> )
{
  cout << "committing the internal state of the chemistry process model" << endl;

  // use this function to commit the internal state
  // of the chemistry process kernel

};



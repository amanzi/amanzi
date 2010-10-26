#include "Chemistry_PK.hpp"
#include "Epetra_MultiVector.h"


Chemistry_PK::Chemistry_PK( Teuchos::ParameterList &parameter_list_,
			    Teuchos::RCP<Chemistry_State> CS_):
  CS(CS_), 
  parameter_list(parameter_list_)
{ 
  // use the constructor to initialize the chemistry process kernel

};

Chemistry_PK::~Chemistry_PK()
{ 

};


void Chemistry_PK::advance( )
{
  cout << "advancing the state of the chemistry process model here" << endl;

  // the MPC will call this function to advance the state 
  // with this particular process kernel

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

  // the MPC will call this function to signal to the 
  // process kernel that it has accepted the 
  // state update, thus, the PK should update
  // possible auxilary state variables here 

};



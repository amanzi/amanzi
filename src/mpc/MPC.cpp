#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "MPC.hpp"
#include "State.hpp"
#include "../chemistry/Chemistry_State.hpp"
#include "../flow/Flow_State.hpp"
#include "../transport/Transport_State.hpp"



MPC::MPC(Teuchos::RCP<Teuchos::ParameterList> Parameters_,
		   Teuchos::RCP<Mesh_maps_simple> mesh_maps_):
  Parameters(Parameters_),
  mesh_maps(mesh_maps_)
  
 {
   int number_of_components = 10; // just a wild guess, should probably come in from input
   
   // create the state object
   S = Teuchos::rcp( new State( number_of_components, mesh_maps) );

  // create auxilary state objects for the process models
  // chemistry...
  CS = Teuchos::rcp( new Chemistry_State( S ) );

  // ... 
  // done creating auxilary state objects for the process models

  
  // create the individual process models
  // chemistry...

  CPK = Teuchos::rcp( new Chemistry_PK(CS) );

  // ...
  // done creating the individual process models

  // chemistry computes new total_component_concentration, so
  // we create storage for that return multi vector

  total_component_concentration_star = Teuchos::rcp(new Epetra_MultiVector( *CS->get_total_component_concentration() ));


}


void MPC::cycle_driver () {
  
  cout << "this is the Chem_MPC::cycle_driver" << endl; 

  CPK->advance( total_component_concentration_star );

  // let's ponder if we can accept the return value...
  // ...
  // o.k. 
  
  // accept total_component_concentration_star
  Teuchos::RCP<Epetra_MultiVector> tcc = S->get_total_component_concentration();
  *tcc = *total_component_concentration_star;

  // make the process model commit its state
  CPK->commit_state( CS );

}

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "MPC.hpp"
#include "State.hpp"
#include "../chemistry/Chemistry_State.hpp"
#include "../chemistry/Chemistry_PK.hpp"
#include "../flow/Flow_State.hpp"
#include "../flow/Flow_PK.hpp"
#include "../transport/Transport_State.hpp"
#include "../transport/Transport_PK.hpp"


MPC::MPC(Teuchos::ParameterList parameter_list_,
	 Teuchos::RCP<Mesh_maps_base> mesh_maps_):
  parameter_list(parameter_list_),
  mesh_maps(mesh_maps_)
  
 {
   
   Teuchos::ParameterList state_parameter_list = 
     parameter_list.sublist("State");

   // create the state object
   S = Teuchos::rcp( new State( state_parameter_list, mesh_maps) );
   
   // create auxilary state objects for the process models
   // chemistry...
   
   CS = Teuchos::rcp( new Chemistry_State( S ) );
   
   TS = Teuchos::rcp( new Transport_State( *S ) );

   FS = Teuchos::rcp( new Flow_State( S ) ); 
   // done creating auxilary state objects for the process models

  
   // create the individual process models
   // chemistry...
   Teuchos::ParameterList chemistry_parameter_list = 
     parameter_list.sublist("Chemistry");
   
   CPK = Teuchos::rcp( new Chemistry_PK(chemistry_parameter_list, CS) );
   
   // transport...
   Teuchos::ParameterList transport_parameter_list = 
     parameter_list.sublist("Transport");
   
   TPK = Teuchos::rcp( new Transport_PK(transport_parameter_list, TS) );
   
   // flow...
   Teuchos::ParameterList flow_parameter_list = 
     parameter_list.sublist("Flow");
   
   FPK = Teuchos::rcp( new Flow_PK(flow_parameter_list, FS) );
   // done creating the individual process models



   
  // chemistry computes new total_component_concentration, so
  // we create storage for that return multi vector

  total_component_concentration_star = Teuchos::rcp(new Epetra_MultiVector( *CS->get_total_component_concentration() ));


}


void MPC::cycle_driver () {
  
  FPK->advance();
  FPK->commit_state(FS);

  TPK->advance();
  TPK->commit_state(TS);

  CPK->advance();
  CPK->commit_state(CS);

}



void MPC::write_mesh()
{

  Teuchos::ParameterList gmv_parameter_list = parameter_list.sublist("GMV");
  std::string gmv_filename = gmv_parameter_list.get<string>("File Name");

  S->write_gmv(gmv_filename);

}

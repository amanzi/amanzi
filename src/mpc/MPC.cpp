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
   
   mpc_parameter_list =  parameter_list.sublist("MPC");
   read_parameter_list();

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

}


void MPC::read_parameter_list()
{
  T0 = mpc_parameter_list.get<double>("Start Time");
  T1 = mpc_parameter_list.get<double>("End Time");
}


void MPC::cycle_driver () {
  
  // so far we only have transport working


  TS->analytic_total_component_concentration();
  TS->analytic_porosity();
  TS->analytic_darcy_flux();
  TS->analytic_water_saturation();
  TS->analytic_water_density();

  // start at time T=T0;
  S->set_time(T0);


  // dump the initial state to gmv
  std::cout << "Time = " <<  S->get_time() << std::endl;
  for( int k=0; k<20; k++ ) printf("%7.4f", (*TPK->get_transport_state_next()->get_total_component_concentration())[0][k]); cout << endl;



  while (S->get_time() <= T1) {
    TPK->advance();
    double transport_dT = TPK->get_transport_dT();

    std::cout << "Transport dT = " << transport_dT << std::endl;

    if (TPK->get_transport_status() == Amanzi_Transport::TRANSPORT_STATE_COMPLETE) 
      {
	// get the transport state and commit it to the state

	RCP<Transport_State> TS_next = TPK->get_transport_state_next();
	
	S->update_total_component_concentration(TS_next->get_total_component_concentration());

	TPK->commit_state(TS);
      }
    else
      {
	// something went wrong
	throw std::exception();
      }

    // update the state
    S->advance_time(transport_dT);
    
    for( int k=0; k<20; k++ ) printf("%7.4f", (*TPK->get_transport_state_next()->get_total_component_concentration())[0][k]); cout << endl;
   
    std::cout << "Time = " <<  S->get_time() << std::endl;

    
  }



}



void MPC::write_mesh()
{

  Teuchos::ParameterList gmv_parameter_list = parameter_list.sublist("GMV");
  std::string gmv_filename = gmv_parameter_list.get<string>("File Name");

  S->write_gmv(gmv_filename);

}

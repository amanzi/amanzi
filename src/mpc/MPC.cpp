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
#include "../utils/gmv_mesh.hh"


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

  // get the GMV data from the parameter list
  Teuchos::ParameterList gmv_parameter_list =  mpc_parameter_list.sublist("GMV");
  
  string gmv_meshfile = gmv_parameter_list.get<string>("Mesh file name");
  string gmv_datafile = gmv_parameter_list.get<string>("Data file name");
  const int gmv_cycle_freq = gmv_parameter_list.get<int>("Dump cycle frequency");
  
  

  GMV::create_mesh_file(*mesh_maps,gmv_meshfile);
  int iter = 0;

  GMV::open_data_file(*mesh_maps, gmv_datafile, iter, 6); 
  GMV::write_time(T0);
  GMV::write_cycle(iter);
  GMV::start_data();
  GMV::write_cell_data( *(*S->get_total_component_concentration())(0), "concentration0");
  GMV::close_data_file();

  while (S->get_time() <= T1) {

    TPK->advance();
    double transport_dT = TPK->get_transport_dT();
    
    std::cout << "MPC: ";
    std::cout << "Cycle = " << iter;
    std::cout << ",  Transport dT = " << transport_dT << std::endl;

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
    
    // advance the 
    iter++;
   
    if (  iter % gmv_cycle_freq   == 0 ) {
      GMV::open_data_file(gmv_meshfile, gmv_datafile,
			  mesh_maps->count_entities(Mesh_data::NODE, OWNED),
			  mesh_maps->count_entities(Mesh_data::NODE, OWNED),
			  iter, 6);
      GMV::write_time(S->get_time());
      GMV::write_cycle(iter);
      GMV::start_data();
      GMV::write_cell_data( *(*S->get_total_component_concentration())(0), "concentration0");
      GMV::close_data_file();    
    }
  }



}



void MPC::write_mesh()
{

  Teuchos::ParameterList gmv_parameter_list = parameter_list.sublist("GMV");
  std::string gmv_filename = gmv_parameter_list.get<string>("File Name");

  S->write_gmv(gmv_filename);

}

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "MPC.hpp"
#include "State.hpp"
#include "Chemistry_State.hpp"
#include "Chemistry_PK.hpp"
#include "Flow_State.hpp"
#include "Flow_PK.hpp"
#include "Transport_State.hpp"
#include "Transport_PK.hpp"
#include "gmv_mesh.hh"
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"

MPC::MPC(Teuchos::ParameterList parameter_list_,
	 Teuchos::RCP<Mesh_maps_base> mesh_maps_):
  parameter_list(parameter_list_),
  mesh_maps(mesh_maps_)
  
{
   
   mpc_parameter_list =  parameter_list.sublist("MPC");
   read_parameter_list();


   // let users selectively disable individual process kernels
   // to allow for testing of the process kernels separately
   transport_enabled = 
     (mpc_parameter_list.get<string>("disable Transport_PK","no") == "no");
   chemistry_enabled =
     (mpc_parameter_list.get<string>("disable Chemistry_PK","no") == "no");
   flow_enabled =
     (mpc_parameter_list.get<string>("disable Flow_PK","no") == "no");
   
   cout << "MPC: The following process kernels are enabled: ";
   if (flow_enabled) {
     cout << "Flow ";
   }
   if (transport_enabled) {
     cout << "Transport "; 
   }
   if (chemistry_enabled) {
     cout << "Chemistry ";
   }
   cout << endl;


   Teuchos::ParameterList state_parameter_list = 
     parameter_list.sublist("State");

   // create the state object
   S = Teuchos::rcp( new State( state_parameter_list, mesh_maps) );
   
   // create auxilary state objects for the process models
   // chemistry...
   
   if (chemistry_enabled) {
     CS = Teuchos::rcp( new Chemistry_State( S ) );

     Teuchos::ParameterList chemistry_parameter_list = 
       parameter_list.sublist("Chemistry");
     
     CPK = Teuchos::rcp( new Chemistry_PK(chemistry_parameter_list, CS) );
   }
   
   if (transport_enabled) {
     TS = Teuchos::rcp( new Transport_State( *S ) );

     Teuchos::ParameterList transport_parameter_list = 
       parameter_list.sublist("Transport");
     
     TPK = Teuchos::rcp( new Transport_PK(transport_parameter_list, TS) );
   }

   if (flow_enabled) {
     FS = Teuchos::rcp( new Flow_State( S ) );

     Teuchos::ParameterList flow_parameter_list = 
       parameter_list.sublist("Flow");
     
     FPK = Teuchos::rcp( new Flow_PK(flow_parameter_list, FS) );      
   }
   // done creating auxilary state objects and  process models

}


void MPC::read_parameter_list()
{
  T0 = mpc_parameter_list.get<double>("Start Time");
  T1 = mpc_parameter_list.get<double>("End Time");
}


void MPC::cycle_driver () {
  
  // so far we only have transport working

  // start at time T=T0;
  S->set_time(T0);


  // get the GMV data from the parameter list
  Teuchos::ParameterList gmv_parameter_list =  mpc_parameter_list.sublist("GMV");
  
  string gmv_meshfile_ = gmv_parameter_list.get<string>("Mesh file name");
  string gmv_datafile_ = gmv_parameter_list.get<string>("Data file name");
  string gmv_prefix = gmv_parameter_list.get<string>("GMV prefix","./");
  
  // create the GMV subdirectory
  boost::filesystem::path prefix(gmv_prefix);
  if (!boost::filesystem::is_directory(prefix.directory_string())) {
    boost::filesystem::create_directory(prefix.directory_string());
  }
  boost::filesystem::path fmesh(gmv_meshfile_);
  boost::filesystem::path fdata(gmv_datafile_);
  boost::filesystem::path slash("/");

  std::stringstream  mesh_filename;
  mesh_filename << prefix.directory_string(); 
  mesh_filename << slash.directory_string();
  mesh_filename << fmesh.directory_string();

  string gmv_meshfile = mesh_filename.str();

  std::stringstream  data_filename;
  data_filename << prefix.directory_string(); 
  data_filename << slash.directory_string();
  data_filename << fdata.directory_string();

  string gmv_datafile = data_filename.str();

  const int gmv_cycle_freq = gmv_parameter_list.get<int>("Dump cycle frequency",100000);
  const double gmv_time_freq = gmv_parameter_list.get<double>("Dump time frequency",1.0e99);
    
  // write the GMV mesh file
  GMV::create_mesh_file(*mesh_maps,gmv_meshfile);
  
  int iter = 0;
  
  int gmv_freq_dump = 0;
  double gmv_time_dump = 0.0;
  int gmv_time_dump_int = 0;

  // first solve the flow equation
  if (flow_enabled) {
    FPK->advance();
    S->update_darcy_flux(FPK->DarcyFlux());
    S->update_pressure(FPK->Pressure());
    FPK->commit_state(FS);
  }
  
  
  // write the GMV data file
  write_mesh_data(gmv_meshfile, gmv_datafile, iter, 6);
  gmv_time_dump_int++;
  
  if (chemistry_enabled || transport_enabled) {
    
    // then iterate transport and chemistry
    while (S->get_time() <= T1) {
      double mpc_dT, chemistry_dT=1e+99, transport_dT=1e+99;

      if (transport_enabled) transport_dT = TPK->calculate_transport_dT();

      mpc_dT = min( transport_dT, chemistry_dT );
      
      std::cout << "MPC: ";
      std::cout << "Cycle = " << iter; 
      std::cout << ",  Time = "<< S->get_time();
      std::cout << ",  Transport dT = " << transport_dT << std::endl;
      
      if (transport_enabled) {
	// now advance transport
	TPK->advance( mpc_dT );	
	if (TPK->get_transport_status() == Amanzi_Transport::TRANSPORT_STATE_COMPLETE) 
	  {
	    // get the transport state and commit it to the state
	    RCP<Transport_State> TS_next = TPK->get_transport_state_next();
	    S->update_total_component_concentration(TS_next->get_total_component_concentration());
	  }
	else
	  {
	    // something went wrong
	    throw std::exception();
	  }
      }
      
      if (chemistry_enabled) {
	// now advance chemistry
	chemistry_dT = transport_dT; // units?
	CPK->advance(chemistry_dT, total_component_concentration_star);
	Chemistry_PK::ChemistryStatus cpk_status = CPK->status();
      }

      // update the time in the state object
      S->advance_time(mpc_dT);
 

      // we're done with this time step, commit the state 
      // in the process kernels
      if (transport_enabled) TPK->commit_state(TS);
      if (chemistry_enabled) CPK->commit_state(CS, chemistry_dT);
	
     
      // advance the iteration count
      iter++;

      gmv_freq_dump = iter;
      gmv_time_dump = S->get_time();
     
      if (  gmv_freq_dump % gmv_cycle_freq   ==  0 ) {

	cout << "Writing GMV file at cycle " << gmv_freq_dump << endl;
	write_mesh_data(gmv_meshfile, gmv_datafile, iter, 6);      

      } else if ( (gmv_time_dump+ mpc_dT/1000.0) / gmv_time_freq  >= gmv_time_dump_int ) {

	gmv_time_dump_int ++;
	

	cout << "Writing GMV file at time T=" << gmv_time_dump << endl;
	write_mesh_data(gmv_meshfile, gmv_datafile, iter, 6);      
	 
      }
    }
    
  }

}



void MPC::write_mesh_data(std::string gmv_meshfile, std::string gmv_datafile, 
			  const int iter, const int digits)
{
  
  GMV::open_data_file(gmv_meshfile, gmv_datafile,
		      mesh_maps->count_entities(Mesh_data::NODE, OWNED),
		      mesh_maps->count_entities(Mesh_data::CELL, OWNED),
		      iter, digits);
  GMV::write_time(S->get_time());
  GMV::write_cycle(iter);
  GMV::start_data();
  
  string basestring = "concentration";
  string suffix = ".00";
  
  for (int nc=0; nc<S->get_number_of_components(); nc++) {
    string concstring(basestring);
    GMV::suffix_no(suffix,nc);
    concstring.append(suffix);
    
    GMV::write_cell_data( *(*S->get_total_component_concentration())(nc), concstring);
  }
  
  GMV::write_cell_data(*S->get_pressure(), "pressure");
  GMV::write_cell_data(*S->get_permeability(), "permeability");
  GMV::write_cell_data(*S->get_porosity(),"porosity");
  
  GMV::close_data_file();     
 
}

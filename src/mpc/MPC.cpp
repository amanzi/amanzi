#include "errors.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Epetra_Comm.h"
#include "Epetra_MpiComm.h"
#include "MPC.hpp"
#include "State.hpp"
#include "chemistry_state.hh"
#include "chemistry_pk.hh"
#include "Flow_State.hpp"
#include "Darcy_PK.hpp"
//#include "SteadyState_Richards_PK.hpp"
#include "Transient_Richards_PK.hpp"
#include "Transport_State.hpp"
#include "Transport_PK.hpp"
// TODO: We are using depreciated parts of boost::filesystem
#define BOOST_FILESYSTEM_VERSION 2
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"


namespace Amanzi
{

  using amanzi::chemistry::Chemistry_State;
  using amanzi::chemistry::Chemistry_PK;
  using amanzi::chemistry::ChemistryException;
  

  MPC::MPC(Teuchos::ParameterList parameter_list_,
	   Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_maps_,
	   Epetra_MpiComm* comm_,
	   Amanzi::ObservationData& output_observations_):
    parameter_list(parameter_list_),
    mesh_maps(mesh_maps_),
    comm(comm_),
    output_observations(output_observations_)
  {
    mpc_init();
  }
  
  
void MPC::mpc_init()
{
  // set the line prefix for output
  this->setLinePrefix("Amanzi::MPC         ");
  // make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);
  
  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&parameter_list,this);    

  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab
  

  
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
   
  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true))	  
    {
      *out << "The following process kernels are enabled: ";
      
      if (flow_enabled) {
	*out << "Flow ";
      }
      if (transport_enabled) {
	*out << "Transport "; 
      }
      if (chemistry_enabled) {
	*out << "Chemistry ";
      }
      *out << std::endl;
    }
     
   if (transport_enabled || flow_enabled || chemistry_enabled) {
     Teuchos::ParameterList state_parameter_list = 
       parameter_list.sublist("State");
     
     // create the state object
     S = Teuchos::rcp( new State( state_parameter_list, mesh_maps) );
   }

   // create auxilary state objects for the process models
   // chemistry...
   
   if (chemistry_enabled) {
     try {
       CS = Teuchos::rcp( new Chemistry_State( S ) );

       Teuchos::ParameterList chemistry_parameter_list = 
           parameter_list.sublist("Chemistry");
     
       CPK = Teuchos::rcp( new Chemistry_PK(chemistry_parameter_list, CS) );
     } catch (ChemistryException& chem_error) {
       std::cout << "MPC: Chemistry_PK constructor returned an error: " 
                 << std::endl << chem_error.what() << std::endl;
       amanzi_throw(chem_error);
     }
   }
   
   // transport...
   if (transport_enabled) {
     TS = Teuchos::rcp( new AmanziTransport::Transport_State( *S ) );

     Teuchos::ParameterList transport_parameter_list = 
       parameter_list.sublist("Transport");
     
     TPK = Teuchos::rcp( new AmanziTransport::Transport_PK(transport_parameter_list, TS) );
   }

   // flow...
   if (flow_enabled) {
     FS = Teuchos::rcp( new Flow_State( S ) );

     Teuchos::ParameterList flow_parameter_list = 
       parameter_list.sublist("Flow");
     
     flow_model = mpc_parameter_list.get<string>("Flow model","Darcy");
     if (flow_model == "Darcy")
       FPK = Teuchos::rcp( new Darcy_PK(flow_parameter_list, FS) );  
     else if (flow_model == "Richards")
       FPK = Teuchos::rcp( new Transient_Richards_PK(flow_parameter_list, FS) );
     else {
       cout << "MPC: unknown flow model: " << flow_model << endl;
       throw std::exception();
     } 
   }
   // done creating auxilary state objects and  process models

   // create the observations
   Teuchos::ParameterList observation_plist = parameter_list.sublist("Observation Data");
   observations = new Amanzi::Unstructured_observations(observation_plist, output_observations);


   // create the visualization object
   if (parameter_list.isSublist("Visualization Data"))
     {
       
       Teuchos::ParameterList vis_parameter_list = 
	 parameter_list.sublist("Visualization Data");
       visualization = new Amanzi::Vis(vis_parameter_list, comm);
       visualization->create_files(*mesh_maps);
     }
   else
     {
       visualization = new Amanzi::Vis();
     }


   // create the restart object
   if (parameter_list.isSublist("Checkpoint Data"))
     {
       
       Teuchos::ParameterList checkpoint_parameter_list = 
	 parameter_list.sublist("Checkpoint Data");
       restart = new Amanzi::Restart(checkpoint_parameter_list, comm);
     }
   else
     {
       restart = new Amanzi::Restart();
     }   


   // are we restarting from a file?
   // assume we're not
   restart_requested = false;
   
   // then check if indeed we are
   if (parameter_list.isSublist("Execution Control"))
     {
       if (parameter_list.sublist("Execution Control").isSublist("Restart from Checkpoint File"))
	 {
	   restart_requested = true;
	   
	   Teuchos::ParameterList restart_parameter_list = 
	     parameter_list.sublist("Execution Control").sublist("Restart from Checkpoint File");
	   
	   restart_from_filename = restart_parameter_list.get<string>("Checkpoint File Name");
	 }
     }
}

void MPC::read_parameter_list()
{
  T0 = mpc_parameter_list.get<double>("Start Time");
  T1 = mpc_parameter_list.get<double>("End Time");
  dT0 = mpc_parameter_list.get<double>("Initial time step",0.0);
  end_cycle = mpc_parameter_list.get<int>("End Cycle",-1);
}


void MPC::cycle_driver ()
{
  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab  

  if (transport_enabled || flow_enabled || chemistry_enabled) {
    // start at time T=T0;
    S->set_time(T0);
  }



  if (chemistry_enabled) {
    try {
      // these are the vectors that chemistry will populate with
      // the names for the auxillary output vectors and the
      // names of components
      std::vector<string> compnames; 
      
      // total view needs this to be outside the constructor 
      CPK->InitializeChemistry();
      CPK->set_chemistry_output_names(&auxnames);
      CPK->set_component_names(&compnames);

      // set the names in the visualization object
      S->set_compnames(compnames);

    } catch (ChemistryException& chem_error) {
      std::cout << "MPC: Chemistry_PK.InitializeChemistry returned an error " 
                << std::endl << chem_error.what() << std::endl;
      Exceptions::amanzi_throw(chem_error);
    }
  }

  int iter = 0;
  S->set_cycle(iter);

  // we cannot at the moment restart in the middle of the 
  // steady state flow calculation, so we check whether a
  // restart was requested, and if not we do the steady
  // state flow calculation

  // if (restart_requested == false)
  //   {
      
  //     // first solve the flow equation to steady state
  //     if (flow_enabled) {
  // 	FPK->advance_to_steady_state();
	
  // 	// reset the time after the steady state solve
  // 	S->set_time(T0);
	
  // 	S->update_darcy_flux(FPK->Flux());
  // 	S->update_pressure(FPK->Pressure());
  // 	FPK->commit_state(FS);
  // 	FPK->GetVelocity(*S->get_darcy_velocity());
	
  // 	if ( flow_model == "Richards") 
  // 	  {
  // 	    Transient_Richards_PK *RPK = dynamic_cast<Transient_Richards_PK*> (&*FPK); 
	    
  // 	    RPK->GetSaturation(*S->get_water_saturation()); 
  // 	  }

  //     }
  //   }
  // else
  if (restart_requested == true)
    {
      // re-initialize the state object
      restart->read_state( *S, restart_from_filename );

      iter = S->get_cycle();
    }
  
  // write visualization output
  if (chemistry_enabled) 
    {
      // get the auxillary data
      Teuchos::RCP<Epetra_MultiVector> aux = CPK->get_extra_chemistry_output_data();
      
      // write visualization data for timestep if requested
      S->write_vis(*visualization, &(*aux), auxnames);
    }
  else
    {
      S->write_vis(*visualization);
    }

  // write a restart dump if requested
  restart->dump_state(*S);

  if (flow_enabled)
    {
      FPK->init_transient(T0, dT0);
    }


  if (flow_enabled || transport_enabled || chemistry_enabled) 
    {
      // make observations
      observations->make_observations(*S);
	
      // we need to create an EpetraMulitVector that will store the 
      // intermediate value for the total component concentration
      total_component_concentration_star =
	Teuchos::rcp(new Epetra_MultiVector(*S->get_total_component_concentration()));
      
      // then iterate transport and chemistry
      while (  (S->get_time() <= T1)  &&   ((end_cycle == -1) || (iter <= end_cycle)) ) {
	double mpc_dT, chemistry_dT=1e+99, transport_dT=1e+99, flow_dT=1e+99;
	
	if (flow_enabled && flow_model == "Richards") 
	  {
	    flow_dT = FPK->get_flow_dT();
	  }
	
	if (transport_enabled) transport_dT = TPK->calculate_transport_dT();

	if (chemistry_enabled) {
	  chemistry_dT = CPK->max_time_step();
	}
	
	mpc_dT = std::min( std::min(flow_dT, transport_dT), chemistry_dT );
	
	if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true))	  
	  {
	    *out << "Cycle = " << iter; 
	    *out << ",  Time = "<< S->get_time() / (60*60*24);
	    *out << ",  dT = " << mpc_dT / (60*60*24)  << std::endl;
	  }

	
	if (flow_enabled)
	  {
	    FPK->advance_transient(mpc_dT);
	  }


	if (transport_enabled) {
	  
	  // now advance transport
	  TPK->advance( mpc_dT );	
	  if (TPK->get_transport_status() == AmanziTransport::TRANSPORT_STATE_COMPLETE) 
	    {
	      // get the transport state and commit it to the state
	      Teuchos::RCP<AmanziTransport::Transport_State> TS_next = TPK->get_transport_state_next();
	      *total_component_concentration_star = *TS_next->get_total_component_concentration();
	    }
	  else
	    {
	      Errors::Message message("MPC: error... Transport_PK.advance returned an error status"); 
	      Exceptions::amanzi_throw(message);
	    }
	} else {
	  
	  *total_component_concentration_star = *S->get_total_component_concentration();
	}
	
	
	if (chemistry_enabled) {
	  try {
	    // now advance chemistry
	    CPK->advance(mpc_dT, total_component_concentration_star);
	    S->update_total_component_concentration(CPK->get_total_component_concentration());
	  } catch (ChemistryException& chem_error) {
	    std::ostringstream error_message;
	    error_message << "MPC: error... Chemistry_PK.advance returned an error status";
	    error_message << chem_error.what();
	    Errors::Message message(error_message.str());
	    Exceptions::amanzi_throw(message);
	    // dump data and give up...
	    // 	  S->update_total_component_concentration(CPK->get_total_component_concentration());
	    
	    // 	  S->advance_time(mpc_dT);
	    // 	  iter++;
	  }
	} else {
	  // commit total_component_concentration_star to the state
	  
	  S->update_total_component_concentration(*total_component_concentration_star);
	}
	
	
	
	// update the time in the state object
	S->advance_time(mpc_dT);
	
	// we're done with this time step, commit the state 
	// in the process kernels
	if (transport_enabled) TPK->commit_state(TS);
	if (chemistry_enabled) CPK->commit_state(CS, mpc_dT);
	
	// advance the iteration count
	iter++;
	S->set_cycle(iter);


	// make observations
	observations->make_observations(*S);

	
	// write visualization if requested
	if (chemistry_enabled) 
	  {
	    // get the auxillary data
	    Teuchos::RCP<Epetra_MultiVector> aux = CPK->get_extra_chemistry_output_data();
	    
	    // write visualization data for timestep if requested
	    S->write_vis(*visualization, &(*aux), auxnames);
	  }
	else
	  {
	    S->write_vis(*visualization);
	  }

	// write restart dump if requested
	restart->dump_state(*S);

      }
      
    }
 


  // dump observations
  output_observations.print(std::cout);
  
  
}


} // close namespace Amanzi

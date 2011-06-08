#include "errors.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "MPC.hpp"
#include "State.hpp"
#include "chemistry-state.hh"
#include "chemistry-pK.hh"
#include "Flow_State.hpp"
#include "Darcy_PK.hpp"
//#include "SteadyState_Richards_PK.hpp"
#include "Transient_Richards_PK.hpp"
#include "Transport_State.hpp"
#include "Transport_PK.hpp"
#include "gmv_mesh.hh"
#ifdef ENABLE_CGNS
#include "cgns_mesh_par.hh"
#include "cgns_mesh.hh"
#endif
// TODO: We are using depreciated parts of boost::filesystem
#define BOOST_FILESYSTEM_VERSION 2
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"


#ifdef ENABLE_CGNS
using namespace CGNS_PAR;
#endif

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
     TS = Teuchos::rcp( new Transport_State( *S ) );

     Teuchos::ParameterList transport_parameter_list = 
       parameter_list.sublist("Transport");
     
     TPK = Teuchos::rcp( new Transport_PK(transport_parameter_list, TS) );
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

}


void MPC::read_parameter_list()
{
  T0 = mpc_parameter_list.get<double>("Start Time");
  T1 = mpc_parameter_list.get<double>("End Time");
  end_cycle = mpc_parameter_list.get<int>("End Cycle",-1);
}


void MPC::cycle_driver () {
  
  if (transport_enabled || flow_enabled || chemistry_enabled) {
    // start at time T=T0;
    S->set_time(T0);
  }



  if (chemistry_enabled) {
    try {
      // total view needs this to be outside the constructor 
      CPK->InitializeChemistry();
      CPK->set_chemistry_output_names(&auxnames);
      CPK->set_component_names(&compnames);
    } catch (ChemistryException& chem_error) {
      std::cout << "MPC: Chemistry_PK.InitializeChemistry returned an error " 
                << std::endl << chem_error.what() << std::endl;
      Exceptions::amanzi_throw(chem_error);
    }
  }
  
  bool gmv_output = mpc_parameter_list.isSublist("GMV");
#ifdef ENABLE_CGNS
  bool cgns_output = mpc_parameter_list.isSublist("CGNS");
  if (cgns_output) {
    cout << "MPC: will write cgns output" << endl;
  }

  // bandre: moved up to the previous chemistry try block
  // if (chemistry_enabled) {
  //   try {
  //     CPK->set_chemistry_output_names(&auxnames);
  //     CPK->set_component_names(compnames);
  //   } catch (ChemistryException& chem_error) {
  //     std::cout << chem_error.what() << std::endl;
  //     Exceptions::amanzi_throw(chem_error);
  //   }
  // }
#endif
  bool gnuplot_output = mpc_parameter_list.get<bool>("Gnuplot output",false);
  if ( mesh_maps->get_comm()->NumProc() != 1 ) {
    gnuplot_output = false;
  }
  if (gnuplot_output) {
    cout << "MPC: will write gnuplot output" << endl;
  }

  const int vizdump_cycle_freq = mpc_parameter_list.get<int>("Viz dump cycle frequency",-1);
  const double vizdump_time_freq = mpc_parameter_list.get<double>("Viz dump time frequency",-1);

  string gmv_mesh_filename_str;
  string gmv_data_filename_path_str;
  string gmv_mesh_filename_path_str;  

  if (gmv_output) {
    // get the GMV data from the parameter list  
    Teuchos::ParameterList gmv_parameter_list =  mpc_parameter_list.sublist("GMV");
    
    create_gmv_paths(gmv_mesh_filename_path_str, gmv_data_filename_path_str,
		     gmv_mesh_filename_str, gmv_parameter_list);
    
    // write the GMV mesh file
    GMV::create_mesh_file(*mesh_maps, gmv_mesh_filename_path_str);
  }
  
#ifdef ENABLE_CGNS
  std::string cgns_filename;
  if (cgns_output) {
    Teuchos::ParameterList cgns_parameter_list =  mpc_parameter_list.sublist("CGNS");
    
    cgns_filename = cgns_parameter_list.get<string>("File name");
    create_mesh_file(*mesh_maps, cgns_filename);


    // print out the parallel distribution
    int rank = mesh_maps->get_comm()->MyPID(); 
    
    Epetra_Vector RNK(mesh_maps->cell_map(false));
    RNK.PutScalar((double)rank);

    open_data_file(cgns_filename);
    create_timestep(0.0, 0, Mesh_data::CELL);
    write_field_data(RNK,"PE");


    
    if (!flow_enabled && !transport_enabled && !chemistry_enabled) close_data_file();
    
  }
#endif  


  int iter = 0;
  
  int vizdump_cycle = 0;
  double vizdump_time = 0.0;
  int vizdump_time_count = 0;

  // first solve the flow equation
  if (flow_enabled) {
    FPK->advance();
    S->update_darcy_flux(FPK->Flux());
    S->update_pressure(FPK->Pressure());
    FPK->commit_state(FS);
    FPK->GetVelocity(*S->get_darcy_velocity());
    
    if ( flow_model == "Richards") 
      {
	Transient_Richards_PK *RPK = dynamic_cast<Transient_Richards_PK*> (&*FPK); 
	
	RPK->GetSaturation(*S->get_water_saturation()); 
      }
  }
  
  if (flow_enabled || transport_enabled || chemistry_enabled) {
    
    if (gmv_output) {
      // write the GMV data file
      write_gmv_data(gmv_data_filename_path_str, gmv_mesh_filename_str, iter, 6);
    }
#ifdef ENABLE_CGNS
    if (cgns_output) {
      // data file is already open, see above

      write_field_data(*S->get_pressure(), "pressure");
      write_field_data(*S->get_permeability(), "permeability");
      write_field_data(*S->get_porosity(),"porosity");
   
      if (flow_model == "Richards")
	{
	  write_field_data(*S->get_water_saturation(),"water saturation");
	}

      for (int nc=0; nc<S->get_number_of_components(); nc++) {
	std::stringstream cname;
	cname << "concentration " << nc;
	
	if (chemistry_enabled) {
	  cname << " " << compnames[nc];
	}

	write_field_data( *(*S->get_total_component_concentration())(nc), cname.str());
      }

      if (flow_enabled) {
	const Epetra_MultiVector &DV = *S->get_darcy_velocity();
	
	write_field_data( *DV(0), "darcy velocity x");
	write_field_data( *DV(1), "darcy velocity y");
	write_field_data( *DV(2), "darcy velocity z");
      }    
      
      if (chemistry_enabled) {
        try {
          // get the auxillary data
          Teuchos::RCP<Epetra_MultiVector> aux = CPK->get_extra_chemistry_output_data();

          // how much of it is there?
          int naux = aux->NumVectors();
          for (int i=0; i<naux; i++) {
            std::stringstream name;
            name << auxnames[i];

            write_field_data( *(*aux)(i), name.str());
          }
        } catch (ChemistryException& chem_error) {
          std::cout << chem_error.what() << std::endl;
          amanzi_throw(chem_error);
	}

      }
      
      close_data_file();
    }
#endif
    if (gnuplot_output) write_gnuplot_data(0, 0.0);


    vizdump_time_count ++;
  }

  if (chemistry_enabled || transport_enabled) {
    // we need to create an EpetraMulitVector that will store the 
    // intermediate value for the total component concentration
    total_component_concentration_star =
      Teuchos::rcp(new Epetra_MultiVector(*S->get_total_component_concentration()));
    
    // then iterate transport and chemistry
    while ( (S->get_time() <= T1)  &&  ((end_cycle == -1) || (iter <= end_cycle)) ) {
      double mpc_dT, chemistry_dT=1e+99, transport_dT=1e+99;
      
      if (transport_enabled) transport_dT = TPK->calculate_transport_dT();
      
      if (chemistry_enabled) {
        chemistry_dT = CPK->max_time_step();
      }
      
      mpc_dT = std::min( transport_dT, chemistry_dT );
      
      std::cout << "MPC: ";
      std::cout << "Cycle = " << iter; 
      std::cout << ",  Time = "<< S->get_time() / (60*60*24);
      std::cout << ",  dT = " << mpc_dT / (60*60*24)  << std::endl;
      
      if (transport_enabled) {
	// now advance transport
	TPK->advance( mpc_dT );	
	if (TPK->get_transport_status() == Amanzi_Transport::TRANSPORT_STATE_COMPLETE) 
	  {
	    // get the transport state and commit it to the state
	    Teuchos::RCP<Transport_State> TS_next = TPK->get_transport_state_next();
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
	  
// #ifdef ENABLE_CGNS
// 	  if (cgns_output) {
// 	    cout << "MPC: Writing to CGNS file at cycle "<< vizdump_cycle << endl;
// 	    write_cgns_data(cgns_filename, iter);
// 	  }
// #endif	
	  
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
      
      // TODO: ask the CPK to dump its data someplace....

      vizdump_cycle = iter;
      vizdump_time = S->get_time();
      
      if (  (vizdump_cycle_freq > 0) && (vizdump_cycle % vizdump_cycle_freq == 0) ) {
	if (gmv_output) {
	  cout << "MPC: Writing GMV file at cycle " << vizdump_cycle << endl;
	  write_gmv_data(gmv_data_filename_path_str, 
			  gmv_mesh_filename_str, iter, 6);
	}
#ifdef ENABLE_CGNS
	if (cgns_output) {
	  cout << "MPC: Writing to CGNS file at cycle "<< vizdump_cycle << endl;
	  write_cgns_data(cgns_filename, iter);
	}
#endif	
	if (gnuplot_output) write_gnuplot_data(iter, vizdump_time + mpc_dT);

      } else if ( (vizdump_time_freq > 0) && ((vizdump_time + mpc_dT/1000.0) / vizdump_time_freq  > vizdump_time_count) ) {
	
	vizdump_time_count ++;
	
	if (gmv_output) {
	  cout << "MPC: Writing GMV file at time T=" << vizdump_time << endl;
	  write_gmv_data(gmv_data_filename_path_str, 
			  gmv_mesh_filename_str, iter, 6);
	}
#ifdef ENABLE_CGNS
	if (cgns_output) {
	  cout << "MPC: Writing to CGNS file at time T=" << vizdump_time << endl;
	  write_cgns_data(cgns_filename, iter);
	}
#endif
	if (gnuplot_output) write_gnuplot_data(iter, vizdump_time + mpc_dT);

      }
    }
    
  }
}




void MPC::write_gmv_data(std::string gmv_datafile_path,
			 std::string gmv_meshfile, const int iter, const int digits)
{
  
  GMV::open_data_file(gmv_meshfile, gmv_datafile_path,
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

#ifdef ENABLE_CGNS
void MPC::write_cgns_data(std::string filename, int iter)
{
  open_data_file(filename);
  
  create_timestep(S->get_time(), iter, Mesh_data::CELL);

  for (int nc=0; nc<S->get_number_of_components(); nc++) {
    
    std::stringstream cname;
    cname << "concentration " << nc;
    
    if (chemistry_enabled) {
      cname << " " << compnames[nc];
    }

    write_field_data( *(*S->get_total_component_concentration())(nc), cname.str());
  }
  
  if (chemistry_enabled) {
    // get the auxillary data
    Teuchos::RCP<Epetra_MultiVector> aux = CPK->get_extra_chemistry_output_data();
    
    // how much of it is there?
    int naux = aux->NumVectors();
    for (int i=0; i<naux; i++) {
      std::stringstream name;
      name << auxnames[i];
      
      write_field_data( *(*aux)(i), name.str()); 
      
    }
    
  }
  close_data_file();     
}
#endif


void MPC::write_gnuplot_data(int iter, double time)
{
  // write each variable to a separate file
  // only works on one PE, not on a truly parallel run
  
  cout << " ...writing gnuplot output" << endl;
  cout << " ...number of components : " << S->get_number_of_components() << endl;
  
  for (int nc=0; nc<S->get_number_of_components(); nc++) {
    std::stringstream fname;
    fname << "conc_" << nc;
    
    if (chemistry_enabled) {
      fname << "_" << compnames[nc];
    }    
    
    fname << "_" << std::setfill('0') << std::setw(5) << iter << ".dat";

    std::filebuf fb;
    fb.open (fname.str().c_str(), std::ios::out);
    ostream os(&fb);
    os << "# time = " << time / (60.0 * 60.0 * 24.0 * 365.25) << " years" << std::endl;

    // now dump the Epetra Vector
    for (int i=0; i< (S->get_total_component_concentration())->MyLength(); i++) 
      os << (*(*S->get_total_component_concentration())(nc))[i] << endl;

    fb.close();
  }
  
  if (chemistry_enabled) {
    // get the auxillary data
    Teuchos::RCP<Epetra_MultiVector> aux = CPK->get_extra_chemistry_output_data();
    
    // how much of it is there?
    int naux = aux->NumVectors();
    for (int n=0; n<naux; n++) {
      std::stringstream fname;
      fname << auxnames[n];
      
      fname << "_" << std::setfill('0') << std::setw(5) << iter << ".dat";
      
      std::filebuf fb;
      fb.open (fname.str().c_str(), std::ios::out);
      ostream os(&fb);      
      os << "# time = " << time / (60.0 * 60.0 * 24.0 * 365.25) << " years" << std::endl;


      // now dump the Epetra Vector
      for (int i=0; i< aux->MyLength(); i++) 
	os << (*(*aux)(n))[i] << endl;
      fb.close();
    }
  }
}



void MPC::create_gmv_paths(std::string  &gmv_mesh_filename_path_str,
			   std::string  &gmv_data_filename_path_str,
			   std::string  &gmv_mesh_filename_str,
			   Teuchos::ParameterList &gmv_parameter_list)
{
  string gmv_meshfile_str = gmv_parameter_list.get<string>("Mesh file name");
  string gmv_datafile_str = gmv_parameter_list.get<string>("Data file name");
  string gmv_prefix_str = gmv_parameter_list.get<string>("GMV prefix","./");
  
  // create the GMV subdirectory if does not exist already
  boost::filesystem::path gmv_prefix_path(gmv_prefix_str);
  if (!boost::filesystem::is_directory(gmv_prefix_path.directory_string())) {
    boost::filesystem::create_directory(gmv_prefix_path.directory_string());
  }
  
  // convert strings to paths
  boost::filesystem::path gmv_meshfile_path(gmv_meshfile_str);
  boost::filesystem::path gmv_datafile_path(gmv_datafile_str);
  boost::filesystem::path slash("/");
  
  // create a portable mesh file path string 
  std::stringstream  gmv_mesh_filename_path_sstr;
  gmv_mesh_filename_path_sstr << gmv_prefix_path.directory_string(); 
  gmv_mesh_filename_path_sstr << slash.directory_string();
  gmv_mesh_filename_path_sstr << gmv_meshfile_path.directory_string();
  gmv_mesh_filename_path_str = gmv_mesh_filename_path_sstr.str();
  
  // create a portable data file path string
  std::stringstream  gmv_data_filename_path_sstr;
  gmv_data_filename_path_sstr << gmv_prefix_path.directory_string(); 
  gmv_data_filename_path_sstr << slash.directory_string();
  gmv_data_filename_path_sstr << gmv_datafile_path.directory_string();
  gmv_data_filename_path_str = gmv_data_filename_path_sstr.str();
  
  // create a portable mesh file string
  std::stringstream  gmv_mesh_filename_sstr;
  gmv_mesh_filename_sstr << gmv_meshfile_path.directory_string();
  gmv_mesh_filename_str = gmv_mesh_filename_sstr.str();  
}

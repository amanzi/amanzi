#ifdef ENABLE_Unstructured
  #include "amanzi_unstructured_grid_simulation_driver.hpp"
#endif
#ifdef ENABLE_Structured
  #include "amanzi_structured_grid_simulation_driver.H"
#endif

#include <iostream>

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "TimerManager.hh"

// include fenv if it exists
#include "boost/version.hpp"
#if (BOOST_VERSION / 100 % 1000 >= 46)
 #include "boost/config.hpp"
 #ifndef BOOST_NO_FENV_H
  #ifdef _GNU_SOURCE
   #define AMANZI_USE_FENV
   #include "boost/detail/fenv.hpp"
  #endif
 #endif
#endif


int main(int argc, char *argv[]) {

#ifdef AMANZI_USE_FENV
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  try {
    Teuchos::GlobalMPISession mpiSession(&argc,&argv,0);
    
    // make sure only PE0 can write to std::cout
    int rank;
    rank = mpiSession.getRank();
    if (rank != 0) cout.rdbuf(0);

    Teuchos::CommandLineProcessor CLP;
    
    CLP.setDocString("\nThe Amanzi driver reads an XML input file and\n"
		     "runs a reactive flow and transport simulation.\n");
    
    std::string xmlInFileName = "options.xml";
    CLP.setOption("xml_file", &xmlInFileName, "XML options file");
    
    CLP.throwExceptions(true);
    
    Teuchos::CommandLineProcessor::EParseCommandLineReturn
      parseReturn = CLP.parse(argc, argv);
       

    // read the main parameter list
    Teuchos::ParameterList driver_parameter_list;
    Teuchos::updateParametersFromXmlFile(xmlInFileName,&driver_parameter_list);
    const Teuchos::ParameterList& mesh_parameter_list = driver_parameter_list.sublist("Mesh");
    driver_parameter_list.set<string>("input file name", xmlInFileName);
    
    // The Mesh list contains a "Structured" sublist or a "Unstructured" sublist, and will 
    // determine which simulation driver to call
    std::string framework;
    if (mesh_parameter_list.isSublist("Structured")) {
      framework = "Structured";
    } else if (mesh_parameter_list.isSublist("Unstructured")) {
      framework = "Unstructured";
    } else {
      amanzi_throw(Errors::Message("The Mesh parameter list must contain one sublist: \"Structured\" or \"Unstructured\""));
    }
    
    Amanzi::Simulator* simulator = NULL;
    
    Amanzi::timer_manager.add("Full Simulation", Amanzi::Timer::ONCE);
    Amanzi::timer_manager.start("Full Simulation");
    
    if (framework=="Structured") {
#ifdef ENABLE_Structured
      simulator = new AmanziStructuredGridSimulationDriver();
#else
      amanzi_throw(Errors::Message("Structured not supported in current build"));
#endif
    } else {
#ifdef ENABLE_Unstructured
      simulator = new AmanziUnstructuredGridSimulationDriver();
#else
      amanzi_throw(Errors::Message("Unstructured not supported in current build"));
#endif
    }
    
    MPI_Comm mpi_comm(MPI_COMM_WORLD);
    Amanzi::ObservationData output_observations;  
    Amanzi::Simulator::ReturnType ret = simulator->Run(mpi_comm, driver_parameter_list, output_observations);

    if (ret == Amanzi::Simulator::FAIL) {
      amanzi_throw(Errors::Message("The amanzi simulator returned an error code, this is most likely due to an error in the mesh creation."));
    }

    // print out observation file in ASCII format 
    Teuchos::ParameterList obs_list;
    if (driver_parameter_list.get<bool>("Native Unstructured Input",true)) {
      obs_list = driver_parameter_list.sublist("Observation Data");
    } else {
      obs_list = driver_parameter_list.sublist("Output").sublist("Observation Data");
    }

    if (obs_list.isParameter("Observation Output Filename")) {
      std::string obs_file = obs_list.get<std::string>("Observation Output Filename");

      if (rank == 0) {
        std::ofstream out; out.open(obs_file.c_str(),std::ios::out);
        if (!out.good()) {
            std::cout << "OPEN PROBLEM" << endl;
        }
      
        out.precision(16);
        out.setf(std::ios::scientific);
	
        out << "Observation Name, Region, Functional, Variable, Time, Value\n";
        out << "===========================================================\n";

	for (Teuchos::ParameterList::ConstIterator i=obs_list.begin(); i!=obs_list.end(); ++i) {
	  std::string label  = obs_list.name(i);
	  std::string _label = label;
#ifdef ENABLE_Structured
	  if (framework=="Structured") _label = Amanzi::AmanziInput::underscore(label);
#endif
	  const Teuchos::ParameterEntry& entry = obs_list.getEntry(label);
	  if (entry.isList()) {
            const Teuchos::ParameterList& ind_obs_list = obs_list.sublist(label);
            for (int j = 0; j < output_observations[_label].size(); j++) {

              if (output_observations[_label][j].is_valid) {
                  if (!out.good()) {
                      std::cout << "PROBLEM BEFORE" << endl;
                  }
                  out << label << ", " 
                    << ind_obs_list.get<std::string>("Region") << ", "
                    << ind_obs_list.get<std::string>("Functional") << ", "
                    << ind_obs_list.get<std::string>("Variable") << ", "
                    << output_observations[_label][j].time << ", "
                      << output_observations[_label][j].value << '\n';
                  if (!out.good()) {
                      std::cout << "PROBLEM AFTER" << endl;
                  }
              }
            }
          }
        }
        out.close();
      }
    }
 
    Amanzi::timer_manager.stop( "Full Simulation" );
    std::cout << "Amanzi::SIMULATION_SUCCESSFUL\n\n";
    
    std::cout << Amanzi::timer_manager << std::endl;
    
    delete simulator;
  }

  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    std::cout << "Amanzi::SIMULATION_FAILED\n";
  }

  // catch all
  catch (...) {
    std::cout << "Amanzi::SIMULATION_FAILED\n";    
  }

}


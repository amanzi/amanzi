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

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,0);

  // make sure only PE0 can write to std::cout
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if (rank!=0) {
    cout.rdbuf(0);
  } 

  Teuchos::CommandLineProcessor CLP;
  
  CLP.setDocString("\nThe Amanzi driver reads an XML input file and\n"
		   "runs a reactive flow and transport simulation.\n");
  
  std::string xmlInFileName = "options.xml";
  CLP.setOption("xml_file", &xmlInFileName, "XML options file");
  
  CLP.throwExceptions(false);
  
  Teuchos::CommandLineProcessor::EParseCommandLineReturn
    parseReturn = CLP.parse(argc, argv);
  
  if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
    return 0;
  }
  if (parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return 1;
  }
  
  MPI_Comm mpi_comm(MPI_COMM_WORLD);
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(mpi_comm);
#else  
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif
  
  // read the main parameter list
  Teuchos::ParameterList driver_parameter_list;
  Teuchos::updateParametersFromXmlFile(xmlInFileName,&driver_parameter_list);
  const Teuchos::ParameterList& mesh_parameter_list = driver_parameter_list.sublist("Mesh");
  
  // Framework is now determined by noting whether the Mesh list contains a "Structured" sublist or
  // a "Unstructured" sublist.
  // Read the "Framework" from the "Mesh" parameter list so that we know which driver to call
  // The available options are "Structured" and the unstructured options "SimpleMesh", "stk::mesh"
  std::string framework;
  if (mesh_parameter_list.isParameter("Framework")) {
    typedef Teuchos::StringToIntegralParameterEntryValidator<int> StrValidator;
    Teuchos::RCP<StrValidator> frameworkValidator
      = Teuchos::rcp(new StrValidator(Teuchos::tuple<std::string>( "Structured", "SimpleMesh", "stk::mesh" )
				      ,"Framework") );
    
    framework
      = frameworkValidator->validateString(Teuchos::getParameter<std::string>(mesh_parameter_list,"Framework"));
  }
  else {
    if (mesh_parameter_list.isSublist("Structured"))
      framework = "Structured";
    else
      framework = "Unstructured";
  }
    
  Amanzi::Simulator* simulator = 0;
  
  if (framework=="Structured")
    {
#ifdef ENABLE_Structured
      simulator = new AmanziStructuredGridSimulationDriver();
#else
      amanzi_throw(Errors::Message("Structured not supported in current build"));
#endif
    }
  else
    {
#ifdef ENABLE_Unstructured
      simulator = new AmanziUnstructuredGridSimulationDriver();
#else
      amanzi_throw(Errors::Message("Unstructured not supported in current build"));
#endif
    }

  Amanzi::ObservationData output_observations;  
  Amanzi::Simulator::ReturnType ret = simulator->Run(mpi_comm,driver_parameter_list,output_observations);

  // print out observation file in ASCII format 
  const Teuchos::ParameterList& obs_list = driver_parameter_list.sublist("Output").sublist("Observation Data");
  if (obs_list.isParameter("Observation Output Filename"))
    {
      std::string obs_file = obs_list.get<std::string>("Observation Output Filename");

      if (rank == 0)
	{
	  std::ofstream out(obs_file.c_str());
	  out.precision(16);
	  out.setf(std::ios::scientific);
	  std::vector<std::string> labellist = output_observations.observationLabels();

	  out << "Observation Name, Region, Functional, Variable, Time, Value\n";
	  out << "===========================================================\n";
	  for (int i = 0; i < labellist.size(); i++)
	    {	  
	      const Teuchos::ParameterList& ind_obs_list = 
		obs_list.sublist(labellist[i]);

	      for (int j = 0; j < output_observations[labellist[i]].size(); j++)
		{
		  if (output_observations[labellist[i]][j].is_valid)
		    out << labellist[i] << ", " 
		        << ind_obs_list.get<std::string>("Region") << ", "
			<< ind_obs_list.get<std::string>("Functional") << ", "
			<< ind_obs_list.get<Teuchos::Array<std::string> >("Variable Macro") << ", "
			<< output_observations[labellist[i]][j].time << ", "
			<< output_observations[labellist[i]][j].value << std::endl;
		}
	    }
	  out.close();
	}
    }

  std::cout << "SIMULATION COMPLETED\n";

  delete simulator;
  delete comm;
}



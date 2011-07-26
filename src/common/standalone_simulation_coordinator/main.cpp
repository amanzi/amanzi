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
  
  // Read the "grid_option" from the parameter list so that we know which driver to call
  // The available options are "Structured" and "Unstructured"
  Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
    gridOptionValidator = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<int>(
        Teuchos::tuple<std::string>( "Structured", "Unstructured" )
        ,"grid_option"
        )
      );

  std::string
    grid_option = gridOptionValidator->validateString(
      Teuchos::getParameter<std::string>(driver_parameter_list,"grid_option")
      );

  Amanzi::Simulator* simulator = 0;
  
  if (grid_option=="Structured")
    {
#ifdef ENABLE_Structured
      simulator = new AmanziStructuredGridSimulationDriver();
#else
      amanzi_throw(Errors::Message("Structured not supported in current build"));
#endif
    }
  else if (grid_option=="Unstructured")
    {
#ifdef ENABLE_Unstructured
      simulator = new AmanziUnstructuredGridSimulationDriver();
#else
      amanzi_throw(Errors::Message("Unstructured not supported in current build"));
#endif
    }

  Amanzi::ObservationData output_observations;  
  Amanzi::Simulator::ReturnType ret = simulator->Run(mpi_comm,driver_parameter_list,output_observations);

  delete simulator;
  delete comm;
}



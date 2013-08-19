#include <iostream>

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "global_verbosity.hh"

#include "dbc.hh"
#include "errors.hh"
#include "amanzi_unstructured_grid_simulation_driver.hh"

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


Teuchos::EVerbosityLevel Amanzi::VerbosityLevel::level_ = Teuchos::VERB_MEDIUM;

int main(int argc, char *argv[])
{

#ifdef AMANZI_USE_FENV
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif


  Teuchos::GlobalMPISession mpiSession(&argc,&argv,0);

  // make sure only PE0 can write to std::cout
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  //if (rank!=0) {
  //  cout.rdbuf(0);
  //} 

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
  Teuchos::ParameterXMLFileReader xmlreader(xmlInFileName);
  driver_parameter_list = xmlreader.getParameters();
  Teuchos::RCP<Teuchos::FancyOStream> fos;
  Teuchos::readVerboseObjectSublist(&driver_parameter_list,&fos,&Amanzi::VerbosityLevel::level_);

  const Teuchos::ParameterList& mesh_parameter_list = driver_parameter_list.sublist("Mesh");

  // Read the "Framework" from the "Mesh" parameter list so that we know which
  // driver to call.  The available options are "Structured" and the
  // unstructured options "SimpleMesh", "stk::mesh"
  typedef Teuchos::StringToIntegralParameterEntryValidator<int> StrValidator;
  Teuchos::RCP<StrValidator> frameworkValidator = Teuchos::rcp(
          new StrValidator(Teuchos::tuple<std::string>( "Structured", "SimpleMesh",
                  "stk::mesh", "MSTK" ) ,"Framework") );

  std::string framework = frameworkValidator->validateString(
          Teuchos::getParameter<std::string>(mesh_parameter_list,"Framework"));

  AmanziUnstructuredGridSimulationDriver* simulator = new AmanziUnstructuredGridSimulationDriver();
  Amanzi::Simulator::ReturnType ret = simulator->Run(mpi_comm, driver_parameter_list);

  delete simulator;
  delete comm;
}



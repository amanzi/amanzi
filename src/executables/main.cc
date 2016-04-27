#include <iostream>

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "VerboseObject_objs.hh"

#include "dbc.hh"
#include "errors.hh"
#include "amanzi_unstructured_grid_simulation_driver.hh"

#include "state_evaluators_registration.hh"

#include "constitutive_relations_eos_registration.hh"
#include "constitutive_relations_surface_subsurface_fluxes_registration.hh"
#include "constitutive_relations_generic_evaluators_registration.hh"

#include "flow_relations_registration.hh"
#include "flow_icy_overland_registration.hh"
#include "flow_overland_pressure_registration.hh"
#include "flow_permafrost_registration.hh"
#include "flow_richards_registration.hh"
#include "transport_amanzi_registration.hh"
//#include "pks_transport_registration.hh"
#include "multiscale_transport_registration.hh"
#include "mdm_transport_registration.hh"

#include "chemistry_amanzi_registration.hh"

//#include "deform_constitutive_relations_porosity_registration.hh"
//#include "deform_prescribed_deformation_registration.hh"
//#include "deform_volumetric_deformation_registration.hh"

#include "energy_advection_diffusion_registration.hh"
// #include "energy_constant_temperature_registration.hh"
#include "energy_constitutive_relations_internal_energy_registration.hh"
#include "energy_constitutive_relations_source_terms_registration.hh"
#include "energy_constitutive_relations_thermal_conductivity_registration.hh"
#include "energy_surface_ice_registration.hh"
#include "energy_two_phase_registration.hh"
#include "energy_three_phase_registration.hh"

//#include "surface_balance_SEB_registration.hh"
//#include "BGC_registration.hh"

// #include "test_pks_registration.hh"

#include "mpc_registration.hh"



// include fenv if it exists
#include "boost/version.hpp"
#if (BOOST_VERSION / 100 % 1000 >= 46)
 #include "boost/config.hpp"
 #ifndef BOOST_NO_FENV_H
//  #ifdef _GNU_SOURCE
   #define AMANZI_USE_FENV
   #include "boost/detail/fenv.hpp"
//  #endif
 #endif
#endif


Teuchos::EVerbosityLevel Amanzi::VerbosityLevel::level_ = Teuchos::VERB_MEDIUM;

int main(int argc, char *argv[])
{

#ifdef AMANZI_USE_FENV
  //  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  feraiseexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
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
  Epetra_MpiComm *comm = new Epetra_MpiComm(mpi_comm);

  // read the main parameter list
  Teuchos::ParameterList driver_parameter_list;
  // Teuchos::ParameterXMLFileReader xmlreader(xmlInFileName);
  // driver_parameter_list = xmlreader.getParameters();
  {
    Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName); 
    driver_parameter_list = *plist;
  }
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
  int ret = simulator->Run(mpi_comm, driver_parameter_list);

  delete simulator;
  delete comm;
}



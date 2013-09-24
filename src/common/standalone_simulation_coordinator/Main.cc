#ifdef ENABLE_Unstructured
#include "AmanziUnstructuredGridSimulationDriver.hh"
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
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/framework/StdOutFormatTarget.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>
#include "DOMTreeErrorReporter.hpp"
#include "ErrorHandler.hpp"
#include "InputTranslator.hh"
//#include "DOMPrintErrorHandler.hpp"
#include "XMLParameterListWriter.hh"

#include "amanzi_version.hh"
#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "TimerManager.hh"
#include "VerboseObject_objs.hh"

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


#ifdef ENABLE_Unstructured
#include "evaluator_reg.hh"
#endif

#include "tpl_versions.h"

#include <iostream>
#include <boost/filesystem.hpp>
using namespace boost::filesystem;



struct RunLog
    : public std::ostream
{
  RunLog(std::ostream& _os);
 protected:
  int rank;
};

int main(int argc, char *argv[]) {

#ifdef AMANZI_USE_FENV
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  Teuchos::GlobalMPISession mpiSession(&argc,&argv,0);
  int rank = mpiSession.getRank();

  try {
    Teuchos::CommandLineProcessor CLP;

    CLP.setDocString("The Amanzi driver reads an XML input file and\n"
                     "runs a reactive flow and transport simulation.\n");

    std::string xmlInFileName = "options.xml";
    CLP.setOption("xml_file", &xmlInFileName, "XML options file");

    std::string xmlSchema = "doc/input_spec/amanzi.xsd";
    CLP.setOption("xml_schema", &xmlSchema, "XML Schema File"); 

    bool print_version(false);
    CLP.setOption("version", "no_version", &print_version, "Print version number and exit.");

    bool print_tpl_versions(false);
    CLP.setOption("tplversions", "no_tplversions", &print_tpl_versions, "Print version numbers of third party libraries and exit.");

    CLP.throwExceptions(false);
    CLP.recogniseAllOptions(true);

    Teuchos::CommandLineProcessor::EParseCommandLineReturn
        parseReturn = CLP.parse(argc, argv);

    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
      // do nothing
    }    
    
    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION) {
      throw 0;
    }

    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_ERROR) {
      throw 0;
    }
    

    // strinigy magic
#define XSTR(s) STR(s)
#define STR(s) #s

    // check for verbose option
    if (print_version) {
      if (rank == 0) {
        std::cout << "Amanzi Version " << XSTR(AMANZI_VERSION) << std::endl;
	std::cout << "HG branch      " << XSTR(AMANZI_HG_BRANCH) << std::endl;
	std::cout << "HG global hash " << XSTR(AMANZI_HG_GLOBAL_HASH) << std::endl;
	std::cout << "HG local id    " << XSTR(AMANZI_HG_LOCAL_ID) << std::endl;
      }
    }

    if (print_tpl_versions) {
      if (rank == 0) {
	std::cout << "Amanzi TPL collection version "<<  XSTR(AMANZI_MAJOR) << "." << XSTR(AMANZI_MINOR) << "." << XSTR(AMANZI_PATCH) << std::endl;
	std::cout << "Third party libraries that this amanzi binary is linked against:" << std::endl;
	std::cout << "ALQUIMIA       " << XSTR(ALQUIMIA_MAJOR) << "." << XSTR(ALQUIMIA_MINOR) << "." << XSTR(ALQUIMIA_PATCH) << std::endl;
	std::cout << "ASCEMIO        " << XSTR(ASCEMIO_MAJOR) << "." << XSTR(ASCEMIO_MINOR) << "." << XSTR(ASCEMIO_PATCH) << std::endl;
	std::cout << "Boost          " << XSTR(Boost_MAJOR) << "." << XSTR(Boost_MINOR) << "." << XSTR(Boost_PATCH) << std::endl;
	std::cout << "BoostCmake     " << XSTR(BoostCmake_MAJOR) << "." << XSTR(BoostCmake_MINOR) << "." << XSTR(BoostCmake_PATCH) << std::endl;
	std::cout << "CCSE           " << XSTR(CCSE_MAJOR) << "." << XSTR(CCSE_MINOR) << "." << XSTR(CCSE_PATCH) << std::endl;
	std::cout << "CURL           " << XSTR(CURL_MAJOR) << "." << XSTR(CURL_MINOR) << "." << XSTR(CURL_PATCH) << std::endl;
	std::cout << "ExodusII       " << XSTR(ExodusII_MAJOR) << "." << XSTR(ExodusII_MINOR) << "." << XSTR(ExodusII_PATCH) << std::endl;
	std::cout << "HDF5           " << XSTR(HDF5_MAJOR) << "." << XSTR(HDF5_MINOR) << "." << XSTR(HDF5_PATCH) << std::endl;
	std::cout << "HYPRE          " << XSTR(HYPRE_MAJOR) << "." << XSTR(HYPRE_MINOR) << "." << XSTR(HYPRE_PATCH) << std::endl;
	std::cout << "METIS          " << XSTR(METIS_MAJOR) << "." << XSTR(METIS_MINOR) << "." << XSTR(METIS_PATCH) << std::endl;	
	std::cout << "MOAB           " << XSTR(MOAB_MAJOR) << "." << XSTR(MOAB_MINOR) << "." << XSTR(MOAB_PATCH) << std::endl;
	std::cout << "MSTK           " << XSTR(MSTK_MAJOR) << "." << XSTR(MSTK_MINOR) << "." << XSTR(MSTK_PATCH) << std::endl;
	std::cout << "NetCDF         " << XSTR(NetCDF_MAJOR) << "." << XSTR(NetCDF_MINOR) << "." << XSTR(NetCDF_PATCH) << std::endl;
	std::cout << "NetCDF_Fortran " << XSTR(NetCDF_Fortran_MAJOR) << "." << XSTR(NetCDF_Fortran_MINOR) << "." << XSTR(NetCDF_Fortran_PATCH) << std::endl;
	std::cout << "OpenMPI        " << XSTR(OpenMPI_MAJOR) << "." << XSTR(OpenMPI_MINOR) << "." << XSTR(OpenMPI_PATCH) << std::endl;	
	std::cout << "ParMetis       " << XSTR(ParMetis_MAJOR) << "." << XSTR(ParMetis_MINOR) << "." << XSTR(ParMetis_PATCH) << std::endl;
	std::cout << "PFLOTRAN       " << XSTR(PFLOTRAN_MAJOR) << "." << XSTR(PFLOTRAN_MINOR) << "." << XSTR(PFLOTRAN_PATCH) << std::endl;
	std::cout << "SEACAS         " << XSTR(SEACAS_MAJOR) << "." << XSTR(SEACAS_MINOR) << "." << XSTR(SEACAS_PATCH) << std::endl;
	std::cout << "SuperLU        " << XSTR(SuperLU_MAJOR) << "." << XSTR(SuperLU_MINOR) << "." << XSTR(SuperLU_PATCH) << std::endl;
	std::cout << "SuperLUDist    " << XSTR(SuperLUDist_MAJOR) << "." << XSTR(SuperLUDist_MINOR) << "." << XSTR(SuperLUDist_PATCH) << std::endl;
	std::cout << "Trilinos       " << XSTR(Trilinos_MAJOR) << "." << XSTR(Trilinos_MINOR) << "." << XSTR(Trilinos_PATCH) << std::endl;
	std::cout << "UnitTest       " << XSTR(UnitTest_MAJOR) << "." << XSTR(UnitTest_MINOR) << "." << XSTR(UnitTest_PATCH) << std::endl;
	std::cout << "XERCES         " << XSTR(XERCES_MAJOR) << "." << XSTR(XERCES_MINOR) << "." << XSTR(XERCES_PATCH) << std::endl;
	std::cout << "ZLIB           " << XSTR(ZLIB_MAJOR) << "." << XSTR(ZLIB_MINOR) << "." << XSTR(ZLIB_PATCH) << std::endl;
      }
    }

    // check if the input file actually exists
    if (!exists(xmlInFileName)) {
      Exceptions::amanzi_throw(Errors::Message("The input file " + xmlInFileName + " does not exist."));
    }


    // EIB - this is the new piece which reads either the new or old input
    /***************************************/
    MPI_Comm mpi_comm(MPI_COMM_WORLD);
    Teuchos::ParameterList driver_parameter_list;
    try {
      xercesc::XMLPlatformUtils::Initialize();
      xercesc::XercesDOMParser *parser = new xercesc::XercesDOMParser;
      parser->setValidationScheme(xercesc::XercesDOMParser::Val_Never);
      bool errorsOccured = false;
      try{
        parser->parse(xmlInFileName.c_str());
      }
      catch (const xercesc::OutOfMemoryException&)
      {
	std::cerr << "OutOfMemoryException" << std::endl;
        errorsOccured = true;
      }
      xercesc::DOMDocument *doc = parser->getDocument();
      xercesc::DOMElement *root = doc->getDocumentElement();
      char* temp2 = xercesc::XMLString::transcode(root->getTagName());
      //DOMImplementation* impl =  DOMImplementationRegistry::getDOMImplementation(X("Core"));
      if (strcmp(temp2,"amanzi_input")==0) {

	if (!exists(xmlSchema)) {
	  Exceptions::amanzi_throw(Errors::Message("The schema file " + xmlSchema + " does not exist."));
	}
	
	//amanzi_throw(Errors::Message("Translation for new input spec is not yet complete, please use old input spec"));
	driver_parameter_list = Amanzi::AmanziNewInput::translate(xmlInFileName, xmlSchema);
	
	//driver_parameter_list.print(std::cout,true,false);
	std::string new_filename(xmlInFileName);
        std::string new_extension("_oldspec.xml");
        size_t pos = new_filename.find(".xml");
        new_filename.replace(pos, (size_t)4, new_extension, (size_t)0, (size_t)12);
        if (rank == 0) {
          printf("Amanzi: writing the translated parameter list to file %s...\n", new_filename.c_str());
	  Teuchos::Amanzi_XMLParameterListWriter XMLWriter;
          Teuchos::XMLObject XMLobj = XMLWriter.toXML(driver_parameter_list);
          std::ofstream xmlfile;
          xmlfile.open(new_filename.c_str());
          xmlfile << XMLobj;
	}

      } else if(strcmp(temp2,"ParameterList")==0) {
	Teuchos::ParameterXMLFileReader xmlreader(xmlInFileName);
        driver_parameter_list = xmlreader.getParameters();
      }
      else {
	amanzi_throw(Errors::Message("Unexpected Error reading input file"));
      }

      // check root tag 
      // if ParameterList - do old and pass thru
      // if amanzi_input  - validate, convert to old

      xercesc::XMLString::release( &temp2) ;
      doc->release();
      xercesc::XMLPlatformUtils::Terminate();
    }
    catch (std::exception& e)
    {
      if (rank == 0) {
        std::cout << e.what() << std::endl;
        std::cout << "Amanzi::XERCES-INPUT_FAILED\n";
      }
    }
    /***************************************/
    // EIB - this is the old stuff I'm replacing
    // *************************************//
    //// read the main parameter list
    //Teuchos::ParameterList driver_parameter_list;
    //// DEPRECATED    Teuchos::updateParametersFromXmlFile(xmlInFileName,&driver_parameter_list);
    //Teuchos::ParameterXMLFileReader xmlreader(xmlInFileName);
    //driver_parameter_list = xmlreader.getParameters();
    // *************************************//

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

    //MPI_Comm mpi_comm(MPI_COMM_WORLD);
    Amanzi::ObservationData output_observations;
    Amanzi::Simulator::ReturnType ret = simulator->Run(mpi_comm, driver_parameter_list, output_observations);

    if (ret == Amanzi::Simulator::FAIL) {
      amanzi_throw(Errors::Message("The amanzi simulator returned an error code, this is most likely due to an error in the mesh creation."));
    }

    // print out observation file in ASCII format
    Teuchos::ParameterList obs_list;

    bool is_native_unstructured = false;
    std::string nui_str = "Native Unstructured Input";
    if (driver_parameter_list.isParameter(nui_str)) {
      is_native_unstructured = driver_parameter_list.get<bool>(nui_str);
    }

    if (is_native_unstructured) {
      obs_list = driver_parameter_list.sublist("Observation Data");
    } else {
      obs_list = driver_parameter_list.sublist("Output").sublist("Observation Data");
    }


    if (obs_list.isParameter("Observation Output Filename")) {
      std::string obs_file = obs_list.get<std::string>("Observation Output Filename");
      if (rank == 0) {
        std::ofstream out;
        out.open(obs_file.c_str(),std::ios::out);
        if (!out.good()) {
          std::cout << "OPEN PROBLEM" << endl;
        }

        out.precision(16);
        out.setf(std::ios::scientific);

        out << "Observation Name, Region, Functional, Variable, Time, Value\n";
        out << "===========================================================\n";

        for (Teuchos::ParameterList::ConstIterator i=obs_list.begin(); i!=obs_list.end(); ++i) {
          std::string label  = obs_list.name(i);
          const Teuchos::ParameterEntry& entry = obs_list.getEntry(label);
          if (entry.isList()) {
            const Teuchos::ParameterList& ind_obs_list = obs_list.sublist(label);
            for (int j = 0; j < output_observations[label].size(); j++) {

              if (output_observations[label][j].is_valid) {
                if (!out.good()) {
                  std::cout << "PROBLEM BEFORE" << endl;
                }
                out << label << ", "
                    << ind_obs_list.get<std::string>("Region") << ", "
                    << ind_obs_list.get<std::string>("Functional") << ", "
                    << ind_obs_list.get<std::string>("Variable") << ", "
                    << output_observations[label][j].time << ", "
                    << output_observations[label][j].value << '\n';
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

    Amanzi::timer_manager.parSync(mpi_comm);

    if (rank == 0) {
      std::cout << "Amanzi::SIMULATION_SUCCESSFUL\n\n";
      std::cout << Amanzi::timer_manager << std::endl;
    }
    delete simulator;
  }

  catch (std::exception& e) {
    if (rank == 0) {
      std::cout << e.what() << std::endl;
      std::cout << "Amanzi::SIMULATION_FAILED\n";
    }
  }

  // catch all
  catch (...) {
    if (rank == 0) {
      std::cout << "Amanzi::SIMULATION_FAILED\n";
    }
  }

}

/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "ErrorHandler.hpp"
#include "SimulatorFactory.hh"
#include "XMLParameterListWriter.hh"

#include "amanzi_version.hh"
#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "VerboseObject_objs.hh"
#include "AmanziComm.hh"

// include fenv if it exists
#include "boost/version.hpp"
#if (BOOST_VERSION / 100 % 1000 >= 46)
#  include "boost/config.hpp"
#  ifndef BOOST_NO_FENV_H
#    ifdef _GNU_SOURCE
#      define AMANZI_USE_FENV
#      include "boost/detail/fenv.hpp"
#    endif
#  endif
#endif

#ifdef ENABLE_Unstructured
#  include "state_evaluators_registration.hh"
#endif

#include "tpl_versions.h"

#include <iostream>
#include <filesystem>


int
main(int argc, char* argv[])
{
#ifdef AMANZI_USE_FENV
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, 0);
  Kokkos::initialize();
  int rank = mpiSession.getRank();

  try {
    Teuchos::CommandLineProcessor CLP;

    CLP.setDocString("The Amanzi driver reads an XML input file and\n"
                     "runs a reactive flow and transport simulation.\n");

    std::string xmlInFileName = "";
    CLP.setOption("xml_file", &xmlInFileName, "XML options file");

    std::string xmlSchema = "";
    CLP.setOption("xml_schema", &xmlSchema, "XML schema file");

    std::string outputPrefix = "";
    CLP.setOption("output_prefix", &outputPrefix, "Path to output data");

    bool print_version(false);
    CLP.setOption(
      "print_version", "no_print_version", &print_version, "Print version number and exit.");

    bool print_tpl_versions(false);
    CLP.setOption("print_tplversions",
                  "no_print_tplversions",
                  &print_tpl_versions,
                  "Print version numbers of third party libraries and exit.");

    bool print_all(false);
    CLP.setOption("print_all", "no_print_all", &print_all, "Print all pre-run information.");

    bool print_paths(false);
    CLP.setOption("print_paths",
                  "no_print_paths",
                  &print_paths,
                  "Print paths of the xml input file and the xml schema file.");

    CLP.throwExceptions(false);
    CLP.recogniseAllOptions(true);

    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = CLP.parse(argc, argv);

    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
      throw std::string("amanzi not run");
    }
    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION) {
      throw std::string("amanzi not run");
    }
    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_ERROR) {
      throw std::string("amanzi not run");
    }

    if (print_all) {
      print_paths = true;
      print_tpl_versions = true;
      print_version = true;
    }

    // strinigy magic
#define XSTR(s) STR(s)
#define STR(s) #s

    // check for verbose option
    if (print_version) {
      if (rank == 0) {
        std::cout << "Amanzi Version  " << XSTR(AMANZI_VERSION) << std::endl;
        std::cout << "GIT branch      " << XSTR(AMANZI_GIT_BRANCH) << std::endl;
        std::cout << "GIT global hash " << XSTR(AMANZI_GIT_GLOBAL_HASH) << std::endl;
      }
    }

    if (print_tpl_versions) {
      if (rank == 0) {
#ifdef AMANZI_MAJOR
        std::cout << "Amanzi TPL collection version " << XSTR(AMANZI_MAJOR) << "."
                  << XSTR(AMANZI_MINOR) << "." << XSTR(AMANZI_PATCH) << std::endl;
#endif
        std::cout << "Third party libraries that this amanzi binary is linked against:"
                  << std::endl;
#ifdef ALQUIMIA_MAJOR
        std::cout << "  ALQUIMIA       " << XSTR(ALQUIMIA_MAJOR) << "." << XSTR(ALQUIMIA_MINOR)
                  << "." << XSTR(ALQUIMIA_PATCH) << std::endl;
#endif
#ifdef ASCEMIO_MAJOR
        std::cout << "  ASCEMIO        " << XSTR(ASCEMIO_MAJOR) << "." << XSTR(ASCEMIO_MINOR) << "."
                  << XSTR(ASCEMIO_PATCH) << std::endl;
#endif
#ifdef Boost_MAJOR
        std::cout << "  Boost          " << XSTR(Boost_MAJOR) << "." << XSTR(Boost_MINOR) << "."
                  << XSTR(Boost_PATCH) << std::endl;
#endif
#ifdef CCSE_MAJOR
        std::cout << "  CCSE           " << XSTR(CCSE_MAJOR) << "." << XSTR(CCSE_MINOR) << "."
                  << XSTR(CCSE_PATCH) << std::endl;
#endif
#ifdef CRUNCHTOPE_MAJOR
        std::cout << "  CRUNCHTOPE     " << XSTR(CRUNCHTOPE_MAJOR) << "." << XSTR(CRUNCHTOPE_MINOR)
                  << "." << XSTR(CRUNCHTOPE_PATCH) << std::endl;
#endif
#ifdef ExodusII_MAJOR
        std::cout << "  ExodusII       " << XSTR(ExodusII_MAJOR) << "." << XSTR(ExodusII_MINOR)
                  << "." << XSTR(ExodusII_PATCH) << std::endl;
#endif
#ifdef HDF5_MAJOR
        std::cout << "  HDF5           " << XSTR(HDF5_MAJOR) << "." << XSTR(HDF5_MINOR) << "."
                  << XSTR(HDF5_PATCH) << std::endl;
#endif
#ifdef HYPRE_MAJOR
        std::cout << "  HYPRE          " << XSTR(HYPRE_MAJOR) << "." << XSTR(HYPRE_MINOR) << "."
                  << XSTR(HYPRE_PATCH) << std::endl;
#endif
#ifdef METIS_MAJOR
        std::cout << "  METIS          " << XSTR(METIS_MAJOR) << "." << XSTR(METIS_MINOR) << "."
                  << XSTR(METIS_PATCH) << std::endl;
#endif
#ifdef MOAB_MAJOR
        std::cout << "  MOAB           " << XSTR(MOAB_MAJOR) << "." << XSTR(MOAB_MINOR) << "."
                  << XSTR(MOAB_PATCH) << std::endl;
#endif
#ifdef MSTK_MAJOR
        std::cout << "  MSTK           " << XSTR(MSTK_MAJOR) << "." << XSTR(MSTK_MINOR) << "."
                  << XSTR(MSTK_PATCH) << std::endl;
#endif
#ifdef Nanoflann_MAJOR
        std::cout << "  Nanoflann      " << XSTR(Nanoflann_MAJOR) << "." << XSTR(Nanoflann_MINOR)
                  << "." << XSTR(Nanoflann_PATCH) << std::endl;
#endif
#ifdef NetCDF_MAJOR
        std::cout << "  NetCDF         " << XSTR(NetCDF_MAJOR) << "." << XSTR(NetCDF_MINOR) << "."
                  << XSTR(NetCDF_PATCH) << std::endl;
#endif
#ifdef NetCDF_Fortran_MAJOR
        std::cout << "  NetCDF_Fortran " << XSTR(NetCDF_Fortran_MAJOR) << "."
                  << XSTR(NetCDF_Fortran_MINOR) << "." << XSTR(NetCDF_Fortran_PATCH) << std::endl;
#endif
#ifdef ParMetis_MAJOR
        std::cout << "  ParMetis       " << XSTR(ParMetis_MAJOR) << "." << XSTR(ParMetis_MINOR)
                  << "." << XSTR(ParMetis_PATCH) << std::endl;
#endif
#ifdef PETSc_MAJOR
        std::cout << "  PETSc          " << XSTR(PETSc_MAJOR) << "." << XSTR(PETSc_MINOR) << "."
                  << XSTR(PETSc_PATCH) << std::endl;
#endif
#ifdef PFLOTRAN_MAJOR
        std::cout << "  PFLOTRAN       " << XSTR(PFLOTRAN_MAJOR) << "." << XSTR(PFLOTRAN_MINOR)
                  << "." << XSTR(PFLOTRAN_PATCH) << std::endl;
#endif
#ifdef SEACAS_MAJOR
        std::cout << "  SEACAS         " << XSTR(SEACAS_MAJOR) << "." << XSTR(SEACAS_MINOR) << "."
                  << XSTR(SEACAS_PATCH) << std::endl;
#endif
#ifdef Silo_MAJOR
        std::cout << "  Silo           " << XSTR(Silo_MAJOR) << "." << XSTR(Silo_MINOR) << "."
                  << XSTR(Silo_PATCH) << std::endl;
#endif
#ifdef SuperLU_MAJOR
        std::cout << "  SuperLU        " << XSTR(SuperLU_MAJOR) << "." << XSTR(SuperLU_MINOR) << "."
                  << XSTR(SuperLU_PATCH) << std::endl;
#endif
#ifdef SuperLUDist_MAJOR
        std::cout << "  SuperLUDist    " << XSTR(SuperLUDist_MAJOR) << "."
                  << XSTR(SuperLUDist_MINOR) << "." << XSTR(SuperLUDist_PATCH) << std::endl;
#endif
#ifdef Trilinos_MAJOR
        std::cout << "  Trilinos       " << XSTR(Trilinos_MAJOR) << "." << XSTR(Trilinos_MINOR)
                  << "." << XSTR(Trilinos_PATCH) << std::endl;
#endif
#ifdef UnitTest_MAJOR
        std::cout << "  UnitTest       " << XSTR(UnitTest_MAJOR) << "." << XSTR(UnitTest_MINOR)
                  << "." << XSTR(UnitTest_PATCH) << std::endl;
#endif
#ifdef XERCES_MAJOR
        std::cout << "  XERCES         " << XSTR(XERCES_MAJOR) << "." << XSTR(XERCES_MINOR) << "."
                  << XSTR(XERCES_PATCH) << std::endl;
#endif
#ifdef ZLIB_MAJOR
        std::cout << "  ZLIB           " << XSTR(ZLIB_MAJOR) << "." << XSTR(ZLIB_MINOR) << "."
                  << XSTR(ZLIB_PATCH) << std::endl;
#endif
      }
    }

    if (print_paths) {
      if (rank == 0) { std::cout << "xml input file:  " << xmlInFileName << std::endl; }
    }

    if (xmlInFileName.size() == 0) {
      if (rank == 0) {
        std::cout << "ERROR: No xml input file was specified. Use the command line option "
                     "--xml_file to specify one."
                  << std::endl;
      }
      throw std::string("Amanzi not run");
    }

    // check if the input file actually exists
    if (!std::filesystem::exists(xmlInFileName)) {
      if (rank == 0) {
        std::cout << "ERROR: The xml input file \"" << xmlInFileName
                  << "\" specified with the command line option --xml_file does not exist."
                  << std::endl;
      }
      throw std::string("Amanzi not run");
    }

    // Create a simulator that corresponds to our input file.
    auto simulator = Amanzi::SimulatorFactory::Create(xmlInFileName, outputPrefix);

    auto comm = Amanzi::getDefaultComm();
    Amanzi::ObservationData observations_data;
    Amanzi::Simulator::ReturnType ret = simulator->Run(comm, observations_data);
    Teuchos::TimeMonitor::summarize();

    if (ret == Amanzi::Simulator::FAIL) {
      amanzi_throw(Errors::Message("The amanzi simulator returned an error code, this is most "
                                   "likely due to an error in the mesh creation."));
    }

    if (rank == 0) { std::cout << "Amanzi::SIMULATION_SUCCESSFUL\n\n"; }
  } catch (std::string& s) {
    if (rank == 0) {
      if (s == "Amanzi not run") { std::cout << "Amanzi::SIMULATION_DID_NOT_RUN\n"; }
    }
    Kokkos::finalize();
    return 1;
  } catch (std::exception& e) {
    if (rank == 0) {
      if (!strcmp(e.what(), "Amanzi not run")) {
        std::cout << "Amanzi::SIMULATION_DID_NOT_RUN\n";
      } else {
        std::cout << e.what() << std::endl;
        std::cout << "Amanzi::SIMULATION_FAILED\n";
      }
    }
    Kokkos::finalize();
    return 1;
  } catch (int& ierr) {
    if (rank == 0) {
      std::cout << "Caught unknown exception with integer code " << ierr
                << ". Known sources: Epetra_MultiVector::AllocateForCopy" << std::endl;
      std::cout << "Amanzi::SIMULATION_FAILED\n";
    }
    Kokkos::finalize();
    return 1;
  }

  // catch all
  catch (...) {
    if (rank == 0) {
      std::cout << "Unknown exception" << std::endl;
      std::cout << "Amanzi::SIMULATION_FAILED\n";
    }
    Kokkos::finalize();
    return 1;
  }
  Kokkos::finalize();
  return 0;
}

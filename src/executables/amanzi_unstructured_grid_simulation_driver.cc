/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $AMANZI_DIR/COPYRIGHT
Author: ??

Effectively stolen from Amanzi, with few modifications.
------------------------------------------------------------------------- */

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "Exodus_readers.hh"
#include "Parallel_Exodus_file.hh"

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "MeshFactory.hh"
#include "Domain.hh"
#include "GeometricModel.hh"
#include "coordinator.hh"

#include "errors.hh"
#include "exceptions.hh"

#include "amanzi_unstructured_grid_simulation_driver.hh"
#include "InputParser.H"
#include "global_verbosity.hh"

Amanzi::Simulator::ReturnType AmanziUnstructuredGridSimulationDriver::Run(
        const MPI_Comm& mpi_comm, Teuchos::ParameterList& input_parameter_list) {

  using Teuchos::OSTab;
  setDefaultVerbLevel(ATS::VerbosityLevel::level_);
  Teuchos::EVerbosityLevel verbLevel = getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = getOStream();
  OSTab tab = getOSTab(); // This sets the line prefix and adds one tab

#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(mpi_comm);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  int rank;
  MPI_Comm_rank(mpi_comm,&rank);

  ParameterList params_copy;
  bool native = input_parameter_list.get<bool>("Native Unstructured Input",false);
  if (! native) {
    params_copy = Amanzi::AmanziInput::translate_state_sublist(input_parameter_list);
  } else {
    params_copy = input_parameter_list;
  }

  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true)) {
    // print parameter list
    *out << "======================> dumping parameter list <======================" <<
      std::endl;
    Teuchos::writeParameterListToXmlOStream(params_copy, *out);
    *out << "======================> done dumping parameter list. <================" <<
      std::endl;
  }

  // ------  domain and geometric model  ------
  Teuchos::ParameterList domain_params = params_copy.sublist("Domain");
  unsigned int spdim = domain_params.get<int>("Spatial Dimension");
  Amanzi::AmanziGeometry::Domain *simdomain_ptr = new Amanzi::AmanziGeometry::Domain(spdim);

  // Under the simulation domain we have different geometric
  // models. We can also have free geometric regions not associated
  // with a geometric model.

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList reg_params = params_copy.sublist("Regions");

  Amanzi::AmanziGeometry::GeometricModelPtr
    geom_model_ptr( new Amanzi::AmanziGeometry::GeometricModel(spdim, reg_params, comm) );

  // Add the geometric model to the domain
  simdomain_ptr->Add_Geometric_Model(geom_model_ptr);

  // If we had geometric models and free regions coexisting then we would 
  // create the free regions here and add them to the simulation domain
  // Nothing to do for now

  // ------ mesh ------
  Amanzi::AmanziMesh::MeshFactory factory(comm);
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

  // select the mesh framework
  Teuchos::ParameterList mesh_plist = params_copy.sublist("Mesh");

  int ierr = 0;
  try {
    std::string framework = mesh_plist.get<string>("Framework");
    Amanzi::AmanziMesh::FrameworkPreference prefs(factory.preference());
    if (framework == Amanzi::AmanziMesh::framework_name(Amanzi::AmanziMesh::Simple)) {
      prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::Simple);
    } else if (framework == Amanzi::AmanziMesh::framework_name(Amanzi::AmanziMesh::MOAB)) {
      prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::MOAB);
    } else if (framework == Amanzi::AmanziMesh::framework_name(Amanzi::AmanziMesh::STKMESH)) {
      prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::STKMESH);
    } else if (framework == Amanzi::AmanziMesh::framework_name(Amanzi::AmanziMesh::MSTK)) {
      prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::MSTK);
    } else if (framework == "") {
      // do nothing
    } else {
      std::string s(framework);
      s += ": specified mesh framework preference not understood";
      amanzi_throw(Errors::Message(s));
    }
    factory.preference(prefs);
  } catch (const Teuchos::Exceptions::InvalidParameterName& e) {
    // do nothing, this means that the "Framework" parameter was not in the input
  } catch (const std::exception& e) {
    std::cerr << rank << ": error: " << e.what() << std::endl;
    ierr++;
  }

  int aerr = 0;
  comm->SumAll(&ierr, &aerr, 1);
  if (aerr > 0) {
    return Amanzi::Simulator::FAIL;
  }

  // Create the mesh
  std::string file(""), format("");

  if (mesh_plist.isSublist("Read Mesh File")) {
    // try to read mesh from file
    Teuchos::ParameterList read_params = mesh_plist.sublist("Read Mesh File");

    if (read_params.isParameter("File")) {
      file = read_params.get<string>("File");
    } else {
      std::cerr << "Must specify File parameter for Read option under Mesh" << std::endl;
      throw std::exception();
    }

    if (read_params.isParameter("Format")) {
      // Is the format one that we can read?
      format = read_params.get<string>("Format");
      if (format != "Exodus II") {
        std::cerr << "Can only read files in Exodus II format" << std::endl;
        throw std::exception();
      }
    } else {
      std::cerr << "Must specify Format parameter for Read option under Mesh" << std::endl;
      throw std::exception();
    }

    if (!file.empty()) {
      ierr = 0;
      try {
        // create the mesh from the file
        Teuchos::RCP<Teuchos::Time> volmeshtime = Teuchos::TimeMonitor::getNewCounter("volume mesh creation");
        Teuchos::TimeMonitor timer(*volmeshtime);
        mesh = factory.create(file, geom_model_ptr);
      } catch (const std::exception& e) {
        std::cerr << rank << ": error: " << e.what() << std::endl;
        ierr++;
      }

      comm->SumAll(&ierr, &aerr, 1);
      if (aerr > 0) return Amanzi::Simulator::FAIL;
    }
  } else if (mesh_plist.isSublist("Generate Mesh")) {
    // try to generate the mesh from data in plist
    Teuchos::ParameterList gen_params = mesh_plist.sublist("Generate Mesh");
    ierr = 0;

    try {
      mesh = factory.create(gen_params, geom_model_ptr);
    } catch (const std::exception& e) {
      std::cerr << rank << ": error: " << e.what() << std::endl;
      ierr++;
    }

    comm->SumAll(&ierr, &aerr, 1);
    if (aerr > 0) return Amanzi::Simulator::FAIL;

  } else {
    std::cerr << rank << ": error: "
              << "Neither Read nor Generate options specified for mesh" << std::endl;
    throw std::exception();
  }

  ASSERT(!mesh.is_null());
  Teuchos::TimeMonitor::summarize();
  Teuchos::TimeMonitor::zeroOutTimers();

  // create the top level Coordinator
  Amanzi::Coordinator coordinator(params_copy, mesh, comm);

  // run the simulation
  coordinator.cycle_driver();

  mesh.reset();
  delete comm;
  delete simdomain_ptr;
  delete geom_model_ptr;
  return Amanzi::Simulator::SUCCESS;
}



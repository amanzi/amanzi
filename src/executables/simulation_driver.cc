/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see ATS_DIR/COPYRIGHT
Author: Ethan Coon (ecoon@lanl.gov)

------------------------------------------------------------------------- */

#include <iostream>

#include <Epetra_MpiComm.h>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "GlobalVerbosity.hh"
#include "VerboseObject.hh"

#include "GeometricModel.hh"
#include "coordinator.hh"
#include "State.hh"

#include "errors.hh"
#include "exceptions.hh"

#include "ats_mesh_factory.hh"
#include "simulation_driver.hh"


int SimulationDriver::Run(
    const MPI_Comm& mpi_comm, Teuchos::ParameterList& plist) {

  Teuchos::RCP<Epetra_MpiComm> comm = Teuchos::rcp(new Epetra_MpiComm(mpi_comm));
  
  // verbosity settings
  setDefaultVerbLevel(Amanzi::VerbosityLevel::level_);
  Teuchos::EVerbosityLevel verbLevel = getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = getOStream();
  Teuchos::OSTab tab = getOSTab(); // This sets the line prefix and adds one tab

  // size, rank
  int rank = comm->MyPID();
  int size = comm->NumProc();

  // print header material
  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true)) {
    // print parameter list
    *out << "======================> dumping parameter list <======================" <<
      std::endl;
    Teuchos::writeParameterListToXmlOStream(plist, *out);
    *out << "======================> done dumping parameter list. <================" <<
      std::endl;
  }

  // create the geometric model and regions
  Teuchos::ParameterList reg_params = plist.sublist("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_params, comm.get()));

  // Create the state.
  Teuchos::ParameterList state_plist = plist.sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));

  // create and register meshes
  //ATS::createMeshes(plist.sublist("mesh"), comm, gm, *S);
  ATS::createMeshes(plist, comm, gm, *S);
  
  // create the top level Coordinator
  ATS::Coordinator coordinator(plist, S, comm.get());

  // run the simulation
  coordinator.cycle_driver();
  return 0;
}



/*
  Simulator

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Markus Brendt (brendt@lanl.gov)
*/

#include <iostream>
#include <fstream>

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "XMLParameterListWriter.hh"

#include "AmanziUnstructuredGridSimulationDriver.hh"
#include "CycleDriver.hh"
#include "MeshInfo.hh"
#include "Domain.hh"
#include "GeometricModel.hh"
#include "InputConverterU.hh"
#include "InputAnalysis.hh"
#include "MeshAudit.hh"
#include "MeshFactory.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "State.hh"
#include "TimerManager.hh"

#include "dbc.hh"
#include "energy_tcm_registration.hh"
#include "energy_iem_registration.hh"
#include "eos_registration.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "mdm_transport_registration.hh"
#include "multiscale_flow_registration.hh"
#include "multiscale_transport_registration.hh"
#include "mpc_pks_registration.hh"
#include "pks_chemistry_registration.hh"
#include "pks_flow_registration.hh"
#include "pks_transport_registration.hh"
#include "pks_energy_registration.hh"
#include "wrm_flow_registration.hh"

// v1 spec constructor -- delete when we get rid of v1.2 spec.
AmanziUnstructuredGridSimulationDriver::AmanziUnstructuredGridSimulationDriver(const std::string& xmlInFileName)
{
  Teuchos::RCP<Teuchos::ParameterList> native = Teuchos::getParametersFromXmlFile(xmlInFileName);
  plist_ = new Teuchos::ParameterList(*native);
}


AmanziUnstructuredGridSimulationDriver::AmanziUnstructuredGridSimulationDriver(const std::string& xmlInFileName,
                                                                               xercesc::DOMDocument* input,
                                                                               const std::string& output_prefix)
{
  int rank = Teuchos::GlobalMPISession::getRank();
  int num_proc = Teuchos::GlobalMPISession::getNProc();

  Amanzi::AmanziInput::InputConverterU converter(xmlInFileName, input, output_prefix);
  plist_ = new Teuchos::ParameterList(converter.Translate(rank, num_proc));
}


AmanziUnstructuredGridSimulationDriver::~AmanziUnstructuredGridSimulationDriver()
{
  delete plist_;
}


Amanzi::Simulator::ReturnType
AmanziUnstructuredGridSimulationDriver::Run(const MPI_Comm& mpi_comm,
                                            Amanzi::ObservationData& observations_data)
{
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  // Teuchos::OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(mpi_comm);
#else  
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  int rank, ierr, aerr, size;
  MPI_Comm_rank(mpi_comm,&rank);
  MPI_Comm_size(mpi_comm,&size);

  //------------ DOMAIN, GEOMETRIC MODEL, ETC ----------------------------
  // Create the simulation domain
  Amanzi::timer_manager.add("Geometric Model creation", Amanzi::Timer::ONCE);
  Amanzi::timer_manager.start("Geometric Model creation");

  Teuchos::ParameterList domain_params = plist_->sublist("domain");
  unsigned int spdim = domain_params.get<int>("spatial dimension");
  
  Amanzi::AmanziGeometry::Domain *simdomain_ptr = new Amanzi::AmanziGeometry::Domain(spdim);

  // Parse the domain description and create.

  // Under the simulation domain we have different geometric
  // models. We can also have free geometric regions not associated
  // with a geometric model.

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList reg_params = plist_->sublist("regions");

  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> geom_model =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(spdim, reg_params, comm));

  // Add the geometric model to the domain
  simdomain_ptr->Add_Geometric_Model(geom_model);

  Amanzi::timer_manager.stop("Geometric Model creation");

  // If we had geometric models and free regions coexisting then we would 
  // create the free regions here and add them to the simulation domain
  // Nothing to do for now


  // ---------------- MESH -----------------------------------------------
  Amanzi::timer_manager.add("Mesh creation",Amanzi::Timer::ONCE);
  Amanzi::timer_manager.start("Mesh creation");

  // Create a verbose object to pass to the mesh_factory and mesh
  Teuchos::ParameterList mesh_params = plist_->sublist("mesh");

  Teuchos::RCP<Amanzi::VerboseObject> mesh_vo =
      Teuchos::rcp(new Amanzi::VerboseObject("Mesh", mesh_params));

  // Create a mesh factory for this geometric model
  Amanzi::AmanziMesh::MeshFactory factory(comm, mesh_vo);

  // Prepare to read/create the mesh specification
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

  // Make sure the unstructured mesh option was chosen
  ierr = 0;
  bool unstructured_option = mesh_params.isSublist("unstructured");

  if (!unstructured_option) {
    std::cerr << "Unstructured simulator invoked for structured mesh request" << std::endl;
    throw std::exception();
  }

  // Read and initialize the unstructured mesh parameters
  Teuchos::ParameterList unstr_mesh_params = mesh_params.sublist("unstructured");

  // Decide on which mesh framework to use
  bool expert_params_specified = unstr_mesh_params.isSublist("expert");

  try {
    Amanzi::AmanziMesh::FrameworkPreference prefs(Amanzi::AmanziMesh::default_preference());

    if (expert_params_specified) {
      Teuchos::ParameterList expert_mesh_params = unstr_mesh_params.sublist("expert");  

      bool framework_specified = expert_mesh_params.isParameter("framework");

      // If caller has specified a particular framework to use, make
      // that the primary framework. Otherwise, use default framework
      // preferences

      if (framework_specified) {
	std::string framework = expert_mesh_params.get<std::string>("framework");

	if (framework == Amanzi::AmanziMesh::framework_name(Amanzi::AmanziMesh::Simple)) {
	  prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::Simple);
	} else if (framework == Amanzi::AmanziMesh::framework_name(Amanzi::AmanziMesh::MSTK)) {
	  prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::MSTK);
	} else if (framework == Amanzi::AmanziMesh::framework_name(Amanzi::AmanziMesh::STKMESH)) {
	  prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::STKMESH);
	} else if (framework == Amanzi::AmanziMesh::framework_name(Amanzi::AmanziMesh::MOAB)) {
	  prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::MOAB);
	// } else if (framework == "") {
	} else {
          std::string s(framework);
          s += ": specified mesh framework preference not understood";
          amanzi_throw(Errors::Message(s));
	}
      }

      bool partitioner_specified = expert_mesh_params.isParameter("partitioner");
      if (partitioner_specified) {
        std::string partitioner = expert_mesh_params.get<std::string>("partitioner");
        if (partitioner == "METIS" || partitioner == "metis")
          factory.set_partitioner(Amanzi::AmanziMesh::Partitioner_type::METIS);
        else if (partitioner == "ZOLTAN_GRAPH" || partitioner == "zoltan_graph")
          factory.set_partitioner(Amanzi::AmanziMesh::Partitioner_type::ZOLTAN_GRAPH);
        else if (partitioner == "ZOLTAN_RCB" || partitioner == "zoltan_rcb")
          factory.set_partitioner(Amanzi::AmanziMesh::Partitioner_type::ZOLTAN_RCB);
      }
    }

    // Create a mesh factory with default or user preferences for a
    // mesh framework
    // prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::MSTK);
    factory.preference(prefs);

  } catch (const Teuchos::Exceptions::InvalidParameterName& e) {
    // do nothing, this means that the "framework" parameter was 
    // not in the input

  } catch (const std::exception& e) {
    std::cerr << rank << ": error: " << e.what() << std::endl;
    ierr++;
  }

  comm->SumAll(&ierr, &aerr, 1);
  if (aerr > 0) {
    return Amanzi::Simulator::FAIL;
  }

  // Read or generate the mesh
  std::string file(""), format("");

  if (unstr_mesh_params.isSublist("read mesh file")) {
    Teuchos::ParameterList read_params = unstr_mesh_params.sublist("read mesh file");
    
    if (read_params.isParameter("file")) {
      file = read_params.get<std::string>("file");
    } else {
      std::cerr << "Must specify File parameter for Read option under mesh" << std::endl;
      throw std::exception();
    }

    if (read_params.isParameter("format")) {
      // Is the format one that we can read?
      format = read_params.get<std::string>("format");

      if (format != "Exodus II" && format != "H5M") {	    
	std::cerr << "Can only read files in Exodus II or H5M format" << std::endl;
	throw std::exception();
      }
    } else {
      std::cerr << "Must specify 'format' parameter for Read option under mesh" << std::endl;
      throw std::exception();
    }

    if (!file.empty()) {
      ierr = 0;
      try {
        // create the mesh from the file
	mesh = factory.create(file, geom_model);
	    
      } catch (const std::exception& e) {
	std::cerr << rank << ": error: " << e.what() << std::endl;
	ierr++;
      }
	  
      comm->SumAll(&ierr, &aerr, 1);
      if (aerr > 0) {
	return Amanzi::Simulator::FAIL;
      }
    }

  } else if (unstr_mesh_params.isSublist("generate mesh")) {  // If Read parameters are specified
    Teuchos::ParameterList gen_params = unstr_mesh_params.sublist("generate mesh");
    ierr = 0;
    
    try {
      // create the mesh by internal generation
      mesh = factory.create(gen_params, geom_model);

    } catch (const std::exception& e) {
      std::cerr << rank << ": error: " << e.what() << std::endl;
      ierr++;
    }
  
    comm->SumAll(&ierr, &aerr, 1);
    if (aerr > 0) {
      return Amanzi::Simulator::FAIL;
    }

  } else {  // If Generate parameters are specified
    std::cerr << rank << ": error: " << "Neither Read nor Generate options specified for mesh" << std::endl;
    throw std::exception();
  }
  Amanzi::timer_manager.stop("Mesh creation");

  AMANZI_ASSERT(!mesh.is_null());

  // Verify mesh and geometric model compatibility
  if (geom_model->dimension() != mesh->space_dimension()) {
    amanzi_throw(Errors::Message("Geometric model and mesh have different dimensions."));
  }

  if (expert_params_specified) {
    Teuchos::ParameterList expert_mesh_params = unstr_mesh_params.sublist("expert");  
    bool verify_mesh_param = expert_mesh_params.isParameter("verify mesh");
    if (verify_mesh_param) {
      bool verify = expert_mesh_params.get<bool>("verify mesh");
      if (verify) {
        if (rank == 0)
          std::cerr << "Verifying mesh with Mesh Audit..." << std::endl;
        if (size == 1) {
          Amanzi::MeshAudit mesh_auditor(mesh);
          int status = mesh_auditor.Verify();
          if (status == 0) {
            if (rank == 0)
              std::cerr << "Mesh Audit confirms that mesh is ok" << std::endl;
          } else {
            std::cerr << "Mesh Audit could not verify correctness of mesh" << std::endl;
            return Amanzi::Simulator::FAIL;
          }
        } else {
          std::ostringstream ofile;
          ofile << "mesh_audit_" << std::setfill('0') << std::setw(4) << rank << ".txt";
          std::ofstream ofs(ofile.str().c_str());
          if (rank == 0)
            std::cerr << "Writing Mesh Audit output to " << ofile.str() << ", etc." << std::endl;
    
          ierr = 0;
          Amanzi::MeshAudit mesh_auditor(mesh, ofs);
          int status = mesh_auditor.Verify();        // check the mesh
          if (status != 0) ierr++;
          
          comm->SumAll(&ierr, &aerr, 1);
          if (aerr == 0) {
            std::cerr << "Mesh Audit confirms that mesh is ok" << std::endl;
          } else {
            if (rank == 0)
              std::cerr << "Mesh Audit could not verify correctness of mesh" << std::endl;
            return Amanzi::Simulator::FAIL;
          }
        }
      }  // if verify
    }  // if verify_mesh_param
    if (expert_mesh_params.isSublist("mesh info")){
      Teuchos::ParameterList mesh_info_list = expert_mesh_params.sublist("mesh info");
      Teuchos::RCP<Amanzi::MeshInfo> mesh_info = Teuchos::rcp(new Amanzi::MeshInfo(mesh_info_list, comm));
      mesh_info->WriteMeshCentroids(*mesh);
    }
  }  // If expert_params_specified


  // -------------- ANALYSIS --------------------------------------------
  Amanzi::InputAnalysis analysis(mesh);
  analysis.Init(*plist_);
  analysis.RegionAnalysis();
  analysis.OutputBCs();


  Teuchos::RCP<Teuchos::ParameterList> glist = Teuchos::rcp(new Teuchos::ParameterList(*plist_));
  Amanzi::CycleDriver cycle_driver(glist, mesh, comm, observations_data);


  cycle_driver.Go();

  // Clean up
  mesh.reset();
  delete simdomain_ptr;
  delete comm;
      
  return Amanzi::Simulator::SUCCESS;
}



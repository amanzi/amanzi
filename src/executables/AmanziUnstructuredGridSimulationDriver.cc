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


#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "XMLParameterListWriter.hh"

#include "AmanziUnstructuredGridSimulationDriver.hh"
#include "CycleDriver.hh"
#include "Domain.hh"
#include "GeometricModel.hh"
#include "InputAnalysis.hh"
#include "InputConverterU.hh"
#include "MeshAudit.hh"
#include "MeshFactory.hh"
#include "MeshInfo.hh"
#include "Mesh_MSTK.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "State.hh"
#include "TimerManager.hh"

#include "bilinear_form_registration.hh"
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
#include "pks_multiphase_registration.hh"
#include "pks_shallow_water_registration.hh"
#include "wrm_flow_registration.hh"
#include "wrmmp_registration.hh"

// v1 spec constructor -- delete when we get rid of v1.2 spec.
AmanziUnstructuredGridSimulationDriver::AmanziUnstructuredGridSimulationDriver(const std::string& xmlInFileName)
{
  plist_ = Teuchos::getParametersFromXmlFile(xmlInFileName);
}


AmanziUnstructuredGridSimulationDriver::AmanziUnstructuredGridSimulationDriver(const std::string& xmlInFileName,
                                                                               xercesc::DOMDocument* input,
                                                                               const std::string& output_prefix)
{
  int rank = Teuchos::GlobalMPISession::getRank();
  int num_proc = Teuchos::GlobalMPISession::getNProc();

  Amanzi::AmanziInput::InputConverterU converter(xmlInFileName, input, output_prefix);
  plist_ = Teuchos::rcp(new Teuchos::ParameterList(converter.Translate(rank, num_proc)));
}


Amanzi::Simulator::ReturnType
AmanziUnstructuredGridSimulationDriver::Run(const MPI_Comm& mpi_comm,
                                            Amanzi::ObservationData& observations_data)
{
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  // Teuchos::OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

#ifdef HAVE_MPI
  auto comm = Teuchos::rcp(new Amanzi::MpiComm_type(mpi_comm));
#else
  auto comm = Amanzi::getCommSelf();
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
  Teuchos::ParameterList& reg_params = plist_->sublist("regions");

  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> geom_model =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(spdim, reg_params, *comm));

  // Add the geometric model to the domain
  simdomain_ptr->Add_Geometric_Model(geom_model);

  Amanzi::timer_manager.stop("Geometric Model creation");

  // If we had geometric models and free regions coexisting then we would 
  // create the free regions here and add them to the simulation domain
  // Nothing to do for now


  // ---------------- MESH -----------------------------------------------
  Amanzi::timer_manager.add("Mesh creation",Amanzi::Timer::ONCE);
  Amanzi::timer_manager.start("Mesh creation");

  // Create a mesh factory for this geometric model
  auto mesh_params = Teuchos::sublist(plist_, "mesh", true);
  auto mesh_vo = Teuchos::rcp(new Amanzi::VerboseObject("Mesh", *mesh_params));
  Amanzi::AmanziMesh::MeshFactory meshfactory(comm, geom_model, mesh_params);

  // Prepare to read/create the mesh specification
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

  // Make sure the unstructured mesh option was chosen
  ierr = 0;
  bool unstructured_option = mesh_params->isSublist("unstructured");

  if (!unstructured_option) {
    std::cerr << "Unstructured simulator invoked for structured mesh request" << std::endl;
    throw std::exception();
  }

  // Read and initialize the unstructured mesh parameters
  Teuchos::ParameterList unstr_mesh_params = mesh_params->sublist("unstructured");

  // Decide on which mesh framework to use
  bool expert_params_specified = unstr_mesh_params.isSublist("expert");

  try {
    Amanzi::AmanziMesh::Preference prefs(Amanzi::AmanziMesh::default_preference());

    if (expert_params_specified) {
      Teuchos::ParameterList expert_mesh_params = unstr_mesh_params.sublist("expert");  

      bool framework_specified = expert_mesh_params.isParameter("framework");

      // If caller has specified a particular framework to use, make
      // that the primary framework. Otherwise, use default framework
      // preferences

      if (framework_specified) {
	std::string framework = expert_mesh_params.get<std::string>("framework");

	if (framework == Amanzi::AmanziMesh::framework_names.at(Amanzi::AmanziMesh::Framework::SIMPLE)) {
	  prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::Framework::SIMPLE);
	} else if (framework == Amanzi::AmanziMesh::framework_names.at(Amanzi::AmanziMesh::Framework::MSTK)) {
	  prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::Framework::MSTK);
	} else if (framework == Amanzi::AmanziMesh::framework_names.at(Amanzi::AmanziMesh::Framework::STK)) {
	  prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::Framework::STK);
	} else if (framework == Amanzi::AmanziMesh::framework_names.at(Amanzi::AmanziMesh::Framework::MOAB)) {
	  prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::Framework::MOAB);
	// } else if (framework == "") {
	} else {
          std::string s(framework);
          s += ": specified mesh framework preference not understood";
          amanzi_throw(Errors::Message(s));
	}
      }
    }

    // Create a mesh meshfactory with default or user preferences for a
    // mesh framework
    meshfactory.set_preference(prefs);

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
	mesh = meshfactory.create(file);
        mesh->PrintMeshStatistics();
	    
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
      mesh = meshfactory.create(gen_params);
      mesh->PrintMeshStatistics();

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

  Teuchos::OSTab tab = mesh_vo->getOSTab();
  *mesh_vo->os() << "CPU time stamp: " << mesh_vo->clock() << std::endl;
  Amanzi::timer_manager.stop("Mesh creation");

  AMANZI_ASSERT(!mesh.is_null());

  // Verify mesh and geometric model compatibility
  if (geom_model->dimension() != mesh->space_dimension()) {
    amanzi_throw(Errors::Message("Geometric model and mesh have different dimensions."));
  }

  if (expert_params_specified) {
    const auto& expert_mesh_params = unstr_mesh_params.sublist("expert");  
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
  }


  // -------------- ANALYSIS --------------------------------------------
  {
    Amanzi::InputAnalysis analysis(mesh, "domain");
    analysis.Init(*plist_);
    analysis.RegionAnalysis();
    analysis.OutputBCs();
  }


  // -------------- STATE -----------------------------------------------
  Teuchos::RCP<Amanzi::State> S;
  if (plist_->isSublist("state")) {
    Teuchos::ParameterList state_plist = plist_->sublist("state");
    S = Teuchos::rcp(new Amanzi::State(state_plist));
  } else {
    Errors::Message message("AmanziUnstructuredSimuluation: xml_file does not contain 'state' sublist\n");
    Exceptions::amanzi_throw(message);
  }

  // Create meshes. This should be moved to a factory of meshes.
  S->RegisterMesh("domain", mesh); 
  
  if (unstr_mesh_params.isSublist("submesh")) {
    if (meshfactory.preference()[0] != Amanzi::AmanziMesh::Framework::MSTK) {
      std::cerr << "Cannot extract a mesh using a non-MSTK framework" << std::endl;
      return Amanzi::Simulator::FAIL;
    }
    const auto& extract_plist = unstr_mesh_params.sublist("submesh");
    std::vector<std::string> names = extract_plist.get<Teuchos::Array<std::string> >("regions").toVector();

    auto mesh_fracture = meshfactory.create(mesh, names, Amanzi::AmanziMesh::FACE);
    mesh_fracture->PrintMeshStatistics();
    S->RegisterMesh("fracture", mesh_fracture);

    {
      Amanzi::InputAnalysis analysis(mesh_fracture, "fracture");
      analysis.Init(*plist_);
      analysis.RegionAnalysis();
    }
  }


  // -------------- EXECUTION -------------------------------------------
  Amanzi::CycleDriver cycle_driver(plist_, S, comm, observations_data);

  cycle_driver.Go();

  // Clean up
  mesh.reset();
  delete simdomain_ptr;
      
  return Amanzi::Simulator::SUCCESS;
}



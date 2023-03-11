/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Brendt (brendt@lanl.gov)
*/

/*
  Simulator

*/

#include <iostream>
#include <fstream>

#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "XMLParameterListWriter.hh"

#include "AmanziComm.hh"
#include "AmanziUnstructuredGridSimulationDriver.hh"
#include "CycleDriver.hh"
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

#include "bilinear_form_registration.hh"
#include "dbc.hh"
#include "energy_tcm_registration.hh"
#include "energy_iem_registration.hh"
#include "eos_registration.hh"
#include "errors.hh"
#include "evaluators_flow_reg.hh"
#include "exceptions.hh"
#include "mdm_transport_registration.hh"
#include "models_flow_reg.hh"
#include "multiphase_evaluators_registration.hh"
#include "multiscale_flow_registration.hh"
#include "multiscale_transport_registration.hh"
#include "mpc_pks_registration.hh"
#include "numerical_flux_registration.hh"
#include "pks_chemistry_registration.hh"
#include "pks_flow_reg.hh"
#include "pks_transport_registration.hh"
#include "pks_energy_registration.hh"
#include "pks_multiphase_registration.hh"
#include "pks_shallow_water_registration.hh"
#include "wrmmp_registration.hh"

// v1 spec constructor -- delete when we get rid of v1.2 spec.
AmanziUnstructuredGridSimulationDriver::AmanziUnstructuredGridSimulationDriver(
  const std::string& xmlInFileName)
{
  plist_ = Teuchos::getParametersFromXmlFile(xmlInFileName);
}


AmanziUnstructuredGridSimulationDriver::AmanziUnstructuredGridSimulationDriver(
  const std::string& xmlInFileName,
  xercesc::DOMDocument* input,
  const std::string& output_prefix)
{
  int rank = Teuchos::GlobalMPISession::getRank();
  int num_proc = Teuchos::GlobalMPISession::getNProc();

  Amanzi::AmanziInput::InputConverterU converter(xmlInFileName, input, output_prefix);
  plist_ = Teuchos::rcp(new Teuchos::ParameterList(converter.Translate(rank, num_proc)));
}


Amanzi::Simulator::ReturnType
AmanziUnstructuredGridSimulationDriver::Run(const Amanzi::Comm_ptr_type& comm,
                                            Amanzi::ObservationData& observations_data)
{
  auto sim_timer = Teuchos::TimeMonitor::getNewCounter("Full Simulation");
  Teuchos::TimeMonitor sim_tm(*sim_timer);

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  int ierr = 0;
  int aerr = 0;
  int rank = comm->MyPID();
  int size = comm->NumProc();

  //------------ GEOMETRIC MODEL ----------------------------
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> geom_model;
  auto geom_timer = Teuchos::TimeMonitor::getNewCounter("Geometric Model Creation");

  { // context for timer
    Teuchos::TimeMonitor geom_tm(*geom_timer);

    Teuchos::ParameterList domain_params = plist_->sublist("domain");
    unsigned int spdim = domain_params.get<int>("spatial dimension");

    // For now create one geometric model from all the regions in the spec
    Teuchos::ParameterList& reg_params = plist_->sublist("regions");
    geom_model = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(spdim, reg_params, *comm));
  }

  // ---------------- MESH -----------------------------------------------
  // Prepare to read/create the mesh specification
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

  auto mesh_params = Teuchos::sublist(plist_, "mesh", true);
  bool unstructured_option = mesh_params->isSublist("unstructured");
  if (!unstructured_option) {
    Errors::Message msg("Unstructured simulator invoked for structured mesh request");
    Exceptions::amanzi_throw(msg);
  }
  Teuchos::ParameterList unstr_mesh_params = mesh_params->sublist("unstructured");
  bool expert_params_specified = unstr_mesh_params.isSublist("expert");

  // Create a mesh factory for this geometric model
  auto mesh_vo = Teuchos::rcp(new Amanzi::VerboseObject("Mesh", *mesh_params));
  Amanzi::AmanziMesh::MeshFactory meshfactory(comm, geom_model, mesh_params);

  auto mesh_timer = Teuchos::TimeMonitor::getNewCounter("Mesh Creation");
  { // context for timer
    Teuchos::TimeMonitor mesh_tm(*mesh_timer);
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

          if (framework ==
              Amanzi::AmanziMesh::framework_names.at(Amanzi::AmanziMesh::Framework::SIMPLE)) {
            prefs.clear();
            prefs.push_back(Amanzi::AmanziMesh::Framework::SIMPLE);
          } else if (framework ==
                     Amanzi::AmanziMesh::framework_names.at(Amanzi::AmanziMesh::Framework::MSTK)) {
            prefs.clear();
            prefs.push_back(Amanzi::AmanziMesh::Framework::MSTK);
          } else if (framework ==
                     Amanzi::AmanziMesh::framework_names.at(Amanzi::AmanziMesh::Framework::MOAB)) {
            prefs.clear();
            prefs.push_back(Amanzi::AmanziMesh::Framework::MOAB);
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
    if (aerr > 0) { return Amanzi::Simulator::FAIL; }

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
        if (aerr > 0) { return Amanzi::Simulator::FAIL; }
      }

    } else if (unstr_mesh_params.isSublist("generate mesh")) { // If Read parameters are specified
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
      if (aerr > 0) { return Amanzi::Simulator::FAIL; }

    } else { // If Generate parameters are specified
      std::cerr << rank << ": error: "
                << "Neither Read nor Generate options specified for mesh" << std::endl;
      throw std::exception();
    }

    Teuchos::OSTab tab = mesh_vo->getOSTab();
    *mesh_vo->os() << "CPU time stamp: " << mesh_vo->clock() << std::endl;
  }
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
        if (rank == 0) std::cerr << "Verifying mesh with Mesh Audit..." << std::endl;
        if (size == 1) {
          Amanzi::MeshAudit mesh_auditor(mesh);
          int status = mesh_auditor.Verify();
          if (status == 0) {
            if (rank == 0) std::cerr << "Mesh Audit confirms that mesh is ok" << std::endl;
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
          int status = mesh_auditor.Verify(); // check the mesh
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
      } // if verify
    }   // if verify_mesh_param
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
    Errors::Message message(
      "AmanziUnstructuredSimuluation: xml_file does not contain 'state' sublist\n");
    Exceptions::amanzi_throw(message);
  }

  // Create meshes. This should be moved to a factory of meshes.
  S->RegisterMesh("domain", mesh);

  if (unstr_mesh_params.isSublist("submesh")) {
    if (meshfactory.get_preference()[0] != Amanzi::AmanziMesh::Framework::MSTK) {
      std::cerr << "Cannot extract a mesh using a non-MSTK framework" << std::endl;
      return Amanzi::Simulator::FAIL;
    }
    const auto& extract_plist = unstr_mesh_params.sublist("submesh");
    std::string domain = extract_plist.get<std::string>("domain name");
    AMANZI_ASSERT(domain == "fracture" || domain == "surface");

    std::vector<std::string> names =
      extract_plist.get<Teuchos::Array<std::string>>("regions").toVector();

    auto submesh = meshfactory.create(mesh, names, Amanzi::AmanziMesh::FACE);
    submesh->PrintMeshStatistics();
    S->RegisterMesh(domain, submesh);

    {
      Amanzi::InputAnalysis analysis(submesh, domain);
      analysis.Init(*plist_);
      analysis.RegionAnalysis();
    }
  }


  // -------------- EXECUTION -------------------------------------------
  Amanzi::CycleDriver cycle_driver(plist_, S, comm, observations_data);
  cycle_driver.Go();

  // Clean up
  mesh.reset();
  return Amanzi::Simulator::SUCCESS;
}

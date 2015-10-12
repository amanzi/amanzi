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

#include "GlobalVerbosity.hh"
#include "VerboseObject.hh"

#include "MeshAudit.hh"
#include "MeshFactory.hh"
#include "ColumnMesh.hh"
#include "Domain.hh"
#include "GeometricModel.hh"
#include "coordinator.hh"
#include "State.hh"

#include "errors.hh"
#include "exceptions.hh"

#include "amanzi_unstructured_grid_simulation_driver.hh"
#include "InputParserIS.hh"

Amanzi::Simulator::ReturnType AmanziUnstructuredGridSimulationDriver::Run(
        const MPI_Comm& mpi_comm, Teuchos::ParameterList& input_parameter_list) {

  using Teuchos::OSTab;
  setDefaultVerbLevel(Amanzi::VerbosityLevel::level_);
  Teuchos::EVerbosityLevel verbLevel = getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = getOStream();
  OSTab tab = getOSTab(); // This sets the line prefix and adds one tab

  Epetra_MpiComm *comm = new Epetra_MpiComm(mpi_comm);

  int rank;
  MPI_Comm_rank(mpi_comm,&rank);
  int size;
  MPI_Comm_size(mpi_comm,&size);

  Teuchos::ParameterList params_copy;
  bool native = input_parameter_list.get<bool>("Native Unstructured Input",false);
  ASSERT(native);
  params_copy = input_parameter_list;

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
    std::string framework = mesh_plist.get<std::string>("Framework");
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
    std::cout << rank << ": error: " << e.what() << std::endl;
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
      file = read_params.get<std::string>("File");
    } else {
      std::cerr << "Must specify File parameter for Read option under Mesh" << std::endl;
      throw std::exception();
    }

    if (read_params.isParameter("Format")) {
      // Is the format one that we can read?
      format = read_params.get<std::string>("Format");
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

  // mesh verification
  bool expert_params_specified = mesh_plist.isSublist("Expert");
  if (expert_params_specified) {
    Teuchos::ParameterList expert_mesh_params = mesh_plist.sublist("Expert");

    bool verify_mesh_param = expert_mesh_params.isParameter("Verify Mesh");
    if (verify_mesh_param) {
      bool verify = expert_mesh_params.get<bool>("Verify Mesh");
      if (verify) {
        std::cerr << "Verifying mesh with Mesh Audit..." << std::endl;
        if (size == 1) {
          Amanzi::MeshAudit mesh_auditor(mesh);
          int status = mesh_auditor.Verify();
          if (status == 0) {
            std::cout << "Mesh Audit confirms that mesh is ok" << std::endl;
          } else {
            std::cout << "Mesh Audit could not verify correctness of mesh" << std::endl;
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
  }  // If expert_params_specified

  // Create the surface mesh if needed
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> surface_mesh = Teuchos::null;
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> surface3D_mesh = Teuchos::null;
  if (mesh_plist.isSublist("Surface Mesh")) {
    Teuchos::ParameterList surface_plist = mesh_plist.sublist("Surface Mesh");

    std::vector<std::string> setnames;
    if (surface_plist.isParameter("surface sideset name")) {
      setnames.push_back(surface_plist.get<std::string>("surface sideset name"));
    } else if (surface_plist.isParameter("surface sideset names")) {
      setnames = surface_plist.get<Teuchos::Array<std::string> >("surface sideset names").toVector();
    } else {
      Errors::Message message("Surface mesh ParameterList needs sideset names.");
      Exceptions::amanzi_throw(message);
    }

    if (mesh->cell_dimension() == 3) {
      surface3D_mesh = factory.create(&*mesh,setnames,Amanzi::AmanziMesh::FACE,false,false);
      surface_mesh = factory.create(&*mesh,setnames,Amanzi::AmanziMesh::FACE,true,false);
    } else {
      surface3D_mesh = mesh;
      surface_mesh = factory.create(&*mesh,setnames,Amanzi::AmanziMesh::CELL,true,false);
    }

    bool surf_expert_params_specified = surface_plist.isSublist("Expert");
    if (surf_expert_params_specified) {
      Teuchos::ParameterList surf_expert_mesh_params = surface_plist.sublist("Expert");
      bool verify_surf_mesh_param = surf_expert_mesh_params.isParameter("Verify Mesh");
      if (verify_surf_mesh_param) {
        bool verify = surf_expert_mesh_params.get<bool>("Verify Mesh");
        if (verify) {
          std::cerr << "Verifying surface mesh with Mesh Audit..." << std::endl;
          if (size == 1) {
            Amanzi::MeshAudit surf_mesh_auditor(surface_mesh);
            int status = surf_mesh_auditor.Verify();
            if (status == 0) {
              std::cout << "Mesh Audit confirms that surface mesh is ok" << std::endl;
            } else {
              std::cout << "Mesh Audit could not verify correctness of surface mesh" << std::endl;
              return Amanzi::Simulator::FAIL;
            }
          } else {
            std::ostringstream ofile;
            ofile << "surf_mesh_audit_" << std::setfill('0') << std::setw(4) << rank << ".txt";
            std::ofstream ofs(ofile.str().c_str());
            if (rank == 0)
              std::cerr << "Writing Surface Mesh Audit output to " << ofile.str() << ", etc." << std::endl;

            ierr = 0;
            Amanzi::MeshAudit surf_mesh_auditor(mesh, ofs);
            int status = surf_mesh_auditor.Verify();        // check the mesh
            if (status != 0) ierr++;

            comm->SumAll(&ierr, &aerr, 1);
            if (aerr == 0) {
              std::cerr << "Surface Mesh Audit confirms that mesh is ok" << std::endl;
            } else {
              if (rank == 0)
                std::cerr << "Surface Mesh Audit could not verify correctness of mesh" << std::endl;
              return Amanzi::Simulator::FAIL;
            }
          }
        }  // if verify
      }  // if verify_mesh_param

      bool export_surf_mesh = surf_expert_mesh_params.isParameter("export mesh to file");
      if (export_surf_mesh) {
        std::string export_surf_mesh_filename =
            surf_expert_mesh_params.get<std::string>("export mesh to file");
        surface3D_mesh->write_to_exodus_file(export_surf_mesh_filename);
      }
    }  // If expert_params_specified
  }

  // column meshes
  std::vector<Teuchos::RCP<Amanzi::AmanziMesh::Mesh> > col_meshes;
  if (mesh_plist.isSublist("Column Meshes")) {
    int nc = mesh->num_columns();
    col_meshes.resize(nc, Teuchos::null);
    for (int c=0; c!=nc; ++c) {
      col_meshes[c] = Teuchos::rcp(new Amanzi::AmanziMesh::ColumnMesh(*mesh, c));
    }
  }  
  
  Teuchos::TimeMonitor::summarize();
  Teuchos::TimeMonitor::zeroOutTimers();

  // Create the state.
  Teuchos::ParameterList state_plist = params_copy.sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));

  // register meshes with state
  bool deformable = mesh_plist.get<bool>("deformable mesh",false);
  S->RegisterDomainMesh(mesh, deformable);
  if (surface3D_mesh != Teuchos::null)
    S->RegisterMesh("surface_3d", surface3D_mesh, deformable);
  if (surface_mesh != Teuchos::null)
    S->RegisterMesh("surface", surface_mesh, deformable);

  if (col_meshes.size() > 0) {
    for (int c=0; c!=col_meshes.size(); ++c) {
      std::stringstream namestream;
      namestream << "column_" << c;
      S->RegisterMesh(namestream.str(), col_meshes[c], deformable);
    }
  }
  
  // create the top level Coordinator
  Amanzi::Coordinator coordinator(params_copy, S, comm);

  // run the simulation
  coordinator.cycle_driver();

  mesh.reset();
  delete comm;
  delete simdomain_ptr;

  // this is poor design
  for (int i=0; i!=geom_model_ptr->Num_Regions(); ++i) {
    delete geom_model_ptr->Region_i(i);
  }
  delete geom_model_ptr;
  return Amanzi::Simulator::SUCCESS;
}



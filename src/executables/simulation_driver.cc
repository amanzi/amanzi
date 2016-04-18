/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see ATS_DIR/COPYRIGHT
Author: Ethan Coon (ecoon@lanl.gov)

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
#include "MeshLogicalFactory.hh"
#include "MeshColumn.hh"
#include "Domain.hh"
#include "GeometricModel.hh"
#include "coordinator.hh"
#include "State.hh"

#include "errors.hh"
#include "exceptions.hh"

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

  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> geom_model_ptr =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_params, comm.get()));

  // create meshes
  Amanzi::AmanziMesh::MeshFactory factory(comm.get());
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

  // select the mesh framework
  Teuchos::ParameterList mesh_plist = plist.sublist("mesh");

  int ierr = 0;
  try {
    std::string framework = mesh_plist.get<std::string>("framework");
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
      Errors::Message msg;
      msg << "\"" << framework << "\" framework preferences not understood.";
      Exceptions::amanzi_throw(msg);
    }
    factory.preference(prefs);
  } catch (const Teuchos::Exceptions::InvalidParameterName& e) {
    // do nothing, this means that the "framework" parameter was not in the input
  } catch (const std::exception& e) {
    std::cout << rank << ": error: " << e.what() << std::endl;
    ierr++;
  }

  int aerr = 0;
  comm->SumAll(&ierr, &aerr, 1);
  if (aerr > 0) {
    return 1;
  }

  // create the base mesh
  if (mesh_plist.isSublist("read mesh file")) {
    // -- from file
    // try to read mesh from file
    Teuchos::ParameterList read_params = mesh_plist.sublist("read mesh file");

    std::string file;
    if (read_params.isParameter("file")) {
      file = read_params.get<std::string>("file");
    } else {
      Errors::Message msg("\"read mesh file\" list missing \"file\" parameter.");
      Exceptions::amanzi_throw(msg);
    }

    std::string format;
    if (read_params.isParameter("format")) {
      // Is the format one that we can read?
      format = read_params.get<std::string>("format");
      if (format != "Exodus II") {
        Errors::Message msg;
        msg << "\"read mesh file\" parameter \"format\" with value \"" << format
            << "\" not understood: valid formats are: \"Exodus II\".";
        Exceptions::amanzi_throw(msg);
      }
    } else {
      Errors::Message msg("\"read mesh file\" parameter \"format\" missing.");
      Exceptions::amanzi_throw(msg);
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
      if (aerr > 0) return 1;
    }


  } else if (mesh_plist.isSublist("generate mesh")) {
    // -- generated mesh
    // try to generate the mesh from data in plist
    Teuchos::ParameterList gen_params = mesh_plist.sublist("generate mesh");
    ierr = 0;

    try {
      mesh = factory.create(gen_params, geom_model_ptr);
    } catch (const std::exception& e) {
      std::cerr << rank << ": error: " << e.what() << std::endl;
      ierr++;
    }

    comm->SumAll(&ierr, &aerr, 1);
    if (aerr > 0) return 1;

  } else if (mesh_plist.isSublist("logical mesh")) {
    // -- from logical mesh file
    Amanzi::AmanziMesh::MeshLogicalFactory fac(comm.get(), geom_model_ptr);
    mesh = fac.Create(mesh_plist.sublist("logical mesh"));
    
  } else if (mesh_plist.isSublist("embedded logical mesh")) {
    Errors::Message msg("\"embedded logical mesh\" option not yet implemented.");
    Exceptions::amanzi_throw(msg);
    
  } else {
    Errors::Message msg("Must specify mesh sublist of type: \"read mesh file\", \"generate mesh\", or \"logical mesh\".");
    Exceptions::amanzi_throw(msg);

  }
  ASSERT(!mesh.is_null());

  // mesh verification
  bool expert_params_specified = mesh_plist.isSublist("expert");
  if (expert_params_specified) {
    Teuchos::ParameterList expert_mesh_params = mesh_plist.sublist("expert");

    bool verify_mesh_param = expert_mesh_params.isParameter("verify mesh");
    if (verify_mesh_param) {
      bool verify = expert_mesh_params.get<bool>("verify mesh");
      if (verify) {
        std::cerr << "Verifying mesh with Mesh Audit..." << std::endl;
        if (size == 1) {
          Amanzi::MeshAudit mesh_auditor(mesh);
          int status = mesh_auditor.Verify();
          if (status == 0) {
            std::cout << "Mesh Audit confirms that mesh is ok" << std::endl;
          } else {
            std::cout << "Mesh Audit could not verify correctness of mesh" << std::endl;
            return 1;
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
            if (rank == 0)
              std::cerr << "Mesh Audit confirms that mesh is ok" << std::endl;
          } else {
            if (rank == 0)
              std::cerr << "Mesh Audit could not verify correctness of mesh" << std::endl;
            return 1;
          }
        }
      }  // if verify
    }  // if verify_mesh_param
  }  // If expert_params_specified

  // Create the surface mesh if needed
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> surface_mesh = Teuchos::null;
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> surface3D_mesh = Teuchos::null;
  if (mesh_plist.isSublist("surface mesh")) {
    Teuchos::ParameterList surface_plist = mesh_plist.sublist("surface mesh");

    std::vector<std::string> setnames;
    if (surface_plist.isParameter("surface sideset name")) {
      setnames.push_back(surface_plist.get<std::string>("surface sideset name"));
    } else if (surface_plist.isParameter("surface sideset names")) {
      setnames = surface_plist.get<Teuchos::Array<std::string> >("surface sideset names").toVector();
    } else {
      Errors::Message message("Surface mesh sublist missing parameter \"surface sideset names\".");
      Exceptions::amanzi_throw(message);
    }

    if (mesh->manifold_dimension() == 3) {
      surface3D_mesh = factory.create(&*mesh,setnames,Amanzi::AmanziMesh::FACE,false,false);
      surface_mesh = factory.create(&*mesh,setnames,Amanzi::AmanziMesh::FACE,true,false);
    } else {
      surface3D_mesh = mesh;
      surface_mesh = factory.create(&*mesh,setnames,Amanzi::AmanziMesh::CELL,true,false);
    }

    bool surf_expert_params_specified = surface_plist.isSublist("expert");
    if (surf_expert_params_specified) {
      Teuchos::ParameterList surf_expert_mesh_params = surface_plist.sublist("expert");
      bool verify_surf_mesh_param = surf_expert_mesh_params.isParameter("verify mesh");
      if (verify_surf_mesh_param) {
        bool verify = surf_expert_mesh_params.get<bool>("verify mesh");
        if (verify) {
          std::cerr << "Verifying surface mesh with Mesh Audit..." << std::endl;
          if (size == 1) {
            Amanzi::MeshAudit surf_mesh_auditor(surface_mesh);
            int status = surf_mesh_auditor.Verify();
            if (status == 0) {
              std::cout << "Mesh Audit confirms that surface mesh is ok" << std::endl;
            } else {
              std::cout << "Mesh Audit could not verify correctness of surface mesh" << std::endl;
              return 1;
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
              return 1;
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
  if (mesh_plist.isSublist("column meshes")) {
    int nc = mesh->num_columns();
    col_meshes.resize(nc, Teuchos::null);
    for (int c=0; c!=nc; ++c) {
      col_meshes[c] = Teuchos::rcp(new Amanzi::AmanziMesh::MeshColumn(*mesh, c));
    }
  }  
  
  Teuchos::TimeMonitor::summarize();
  Teuchos::TimeMonitor::zeroOutTimers();

  // Create the state.
  Teuchos::ParameterList state_plist = plist.sublist("state");
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
  Amanzi::Coordinator coordinator(plist, S, comm.get());

  // run the simulation
  coordinator.cycle_driver();
  mesh.reset();
  return 0;
}



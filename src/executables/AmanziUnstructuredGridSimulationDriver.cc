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
#include "Domain.hh"
#include "GeometricModel.hh"
#include "InputTranslator.hh"
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
#include "multiscale_flow_registration.hh"
#include "multiscale_transport_registration.hh"
#include "mpc_pks_registration.hh"
#include "pks_chemistry_registration.hh"
#include "pks_flow_registration.hh"
#include "pks_transport_registration.hh"
#include "pks_energy_registration.hh"
#include "wrm_flow_registration.hh"

using namespace std;

// v1 spec constructor -- delete when we get rid of v1.2 spec.
AmanziUnstructuredGridSimulationDriver::AmanziUnstructuredGridSimulationDriver(const string& xmlInFileName)
{
  string spec;
  Teuchos::ParameterList driver_parameter_list = Amanzi::AmanziNewInput::translate(xmlInFileName, spec);
  if (driver_parameter_list.numParams() == 0) // empty list, must be native spec.
  {
    Teuchos::RCP<Teuchos::ParameterList> native = Teuchos::getParametersFromXmlFile(xmlInFileName);
    plist_ = new Teuchos::ParameterList(*native);
  }
  else
    plist_ = new Teuchos::ParameterList(Amanzi::AmanziNewInput::translate(xmlInFileName, spec));
}


AmanziUnstructuredGridSimulationDriver::AmanziUnstructuredGridSimulationDriver(const string& xmlInFileName,
                                                                               xercesc::DOMDocument* input)
{
  int rank = Teuchos::GlobalMPISession::getRank();
  int num_proc = Teuchos::GlobalMPISession::getNProc();

  Amanzi::AmanziInput::InputConverterU converter(xmlInFileName, input);
  plist_ = new Teuchos::ParameterList(converter.Translate(rank, num_proc));
}


AmanziUnstructuredGridSimulationDriver::~AmanziUnstructuredGridSimulationDriver()
{
  delete plist_;
}


Amanzi::Simulator::ReturnType
AmanziUnstructuredGridSimulationDriver::Run(const MPI_Comm& mpi_comm,
                                            Amanzi::ObservationData& output_observations)
{
  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(mpi_comm);
#else  
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  int rank, ierr, aerr, size;
  MPI_Comm_rank(mpi_comm,&rank);
  MPI_Comm_size(mpi_comm,&size);

  //------------ DOMAIN, GEOMETRIC MODEL, ETC ----------------------------
  // Create a VerboseObject to pass to the geometric model class 
  Amanzi::VerboseObject *gmverbobj = new Amanzi::VerboseObject("Geometric Model", *plist_);

  // Create the simulation domain
  Amanzi::timer_manager.add("Geometric Model creation",Amanzi::Timer::ONCE);
  Amanzi::timer_manager.start("Geometric Model creation");

  Teuchos::ParameterList domain_params = plist_->sublist("Domain");
  unsigned int spdim = domain_params.get<int>("Spatial Dimension");
  
  Amanzi::AmanziGeometry::Domain *simdomain_ptr = new Amanzi::AmanziGeometry::Domain(spdim);

  // Parse the domain description and create.

  // Under the simulation domain we have different geometric
  // models. We can also have free geometric regions not associated
  // with a geometric model.

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList reg_params = plist_->sublist("Regions");

  Amanzi::AmanziGeometry::GeometricModelPtr 
      geom_model_ptr(new Amanzi::AmanziGeometry::GeometricModel(spdim, reg_params, comm));


  // Add the geometric model to the domain
  simdomain_ptr->Add_Geometric_Model(geom_model_ptr);

  Amanzi::timer_manager.stop("Geometric Model creation");

  // If we had geometric models and free regions coexisting then we would 
  // create the free regions here and add them to the simulation domain
  // Nothing to do for now


  // ---------------- MESH -----------------------------------------------
  Amanzi::timer_manager.add("Mesh creation",Amanzi::Timer::ONCE);
  Amanzi::timer_manager.start("Mesh creation");

  // Create a Verbose object to pass to the mesh_factory and mesh
  Amanzi::VerboseObject *meshverbobj = new Amanzi::VerboseObject("Mesh", *plist_);

  // Create a mesh factory for this geometric model
  Amanzi::AmanziMesh::MeshFactory factory(comm, meshverbobj) ;

  // Prepare to read/create the mesh specification
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

  // get the Mesh sublist
  ierr = 0;
  Teuchos::ParameterList mesh_params = plist_->sublist("Mesh");

  // Make sure the unstructured mesh option was chosen
  bool unstructured_option = mesh_params.isSublist("Unstructured");

  if (!unstructured_option) {
    std::cerr << "Unstructured simulator invoked for structured mesh request" << std::endl;
    throw std::exception();
  }

  // Read and initialize the unstructured mesh parameters
  Teuchos::ParameterList unstr_mesh_params = mesh_params.sublist("Unstructured");

  // Decide on which mesh framework to use
  bool expert_params_specified = unstr_mesh_params.isSublist("Expert");

  try {
    Amanzi::AmanziMesh::FrameworkPreference prefs(Amanzi::AmanziMesh::default_preference());

    if (expert_params_specified) {
      Teuchos::ParameterList expert_mesh_params = unstr_mesh_params.sublist("Expert");  

      bool framework_specified = expert_mesh_params.isParameter("Framework");

      // If caller has specified a particular framework to use, make
      // that the primary framework. Otherwise, use default framework
      // preferences

      if (framework_specified) {
	std::string framework = expert_mesh_params.get<std::string>("Framework");

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
    }

    // Create a mesh factory with default or user preferences for a
    // mesh framework
    // prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::MSTK);
    factory.preference(prefs);

  } catch (const Teuchos::Exceptions::InvalidParameterName& e) {
    // do nothing, this means that the "Framework" parameter was 
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

  if (unstr_mesh_params.isSublist("Read Mesh File")) {
    Teuchos::ParameterList read_params = unstr_mesh_params.sublist("Read Mesh File");
    
    if (read_params.isParameter("File")) {
      file = read_params.get<std::string>("File");
    } else {
      std::cerr << "Must specify File parameter for Read option under Mesh" << std::endl;
      throw std::exception();
    }

    if (read_params.isParameter("Format")) {
      // Is the format one that we can read?
      format = read_params.get<std::string>("Format");

      if (format != "Exodus II" && format != "H5M") {	    
	std::cerr << "Can only read files in Exodus II or H5M format" << std::endl;
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
	mesh = factory.create(file, geom_model_ptr);
	    
      } catch (const std::exception& e) {
	std::cerr << rank << ": error: " << e.what() << std::endl;
	ierr++;
      }
	  
      comm->SumAll(&ierr, &aerr, 1);
      if (aerr > 0) {
	return Amanzi::Simulator::FAIL;
      }
    }

  } else if (unstr_mesh_params.isSublist("Generate Mesh")) {  // If Read parameters are specified
    Teuchos::ParameterList gen_params = unstr_mesh_params.sublist("Generate Mesh");
    ierr = 0;
    
    try {
      // create the mesh by internal generation
      mesh = factory.create(gen_params, geom_model_ptr);

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

  ASSERT(!mesh.is_null());


  if (expert_params_specified) {
    Teuchos::ParameterList expert_mesh_params = unstr_mesh_params.sublist("Expert");  
    bool verify_mesh_param = expert_mesh_params.isParameter("Verify Mesh");
    if (verify_mesh_param) {
      bool verify = expert_mesh_params.get<bool>("Verify Mesh");
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
  }  // If expert_params_specified

  // -------------- ANALYSIS --------------------------------------------
  Amanzi::InputAnalysis analysis(mesh);
  analysis.Init(*plist_);
  analysis.RegionAnalysis();
  analysis.OutputBCs();

  Teuchos::RCP<Teuchos::ParameterList> glist = Teuchos::rcp(new Teuchos::ParameterList(*plist_));
  Amanzi::CycleDriver cycle_driver(glist, mesh, comm, output_observations);

  cycle_driver.Go();

  // Clean up
  mesh.reset();
  delete meshverbobj;
  delete gmverbobj;
  delete simdomain_ptr;
  delete comm;
      
  return Amanzi::Simulator::SUCCESS;
}



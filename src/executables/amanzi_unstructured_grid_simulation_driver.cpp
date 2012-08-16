#include <iostream>
#include <fstream>

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "amanzi_unstructured_grid_simulation_driver.hpp"

#include "Teuchos_ParameterList.hpp"
// #include "Teuchos_XMLParameterListHelpers.hpp"
#include "XMLParameterListWriter.hpp"

#include "MeshFactory.hh"
#include "State.hpp"
#include "MPC.hpp"
#include "Domain.hh"
#include "GeometricModel.hh"
#include "MeshAudit.hh"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

#include "InputParserIS.hh"

#include "TimerManager.hh"


Amanzi::Simulator::ReturnType
AmanziUnstructuredGridSimulationDriver::Run(const MPI_Comm& mpi_comm,
                                            Teuchos::ParameterList& input_parameter_list,
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

  bool native = input_parameter_list.get<bool>("Native Unstructured Input",false);
  
  Teuchos::ParameterList new_list; 
  Teuchos::ParameterList sub_list;
  
  if (! native) {
    new_list = Amanzi::AmanziInput::translate( &input_parameter_list, comm->NumProc() );

    std::string verbosity = input_parameter_list.sublist("Execution Control").get<std::string>("Verbosity","Low");
    
    if ( verbosity == "None" ) {
      verbLevel = Teuchos::VERB_NONE;
    } else if ( verbosity == "Low" ) {
      verbLevel = Teuchos::VERB_LOW;
    } else if ( verbosity == "Medium" ) {
      verbLevel = Teuchos::VERB_MEDIUM;
    } else if ( verbosity == "High" ) {
      verbLevel = Teuchos::VERB_HIGH;
    } else if ( verbosity == "Extreme" ) {
      verbLevel = Teuchos::VERB_EXTREME;
    } 
      
  } else {
    new_list = input_parameter_list;
  }
  
  if ((comm->MyPID() == 0)&& (includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)))
          Amanzi::AmanziInput::output_boundary_conditions(&new_list);
  
  if (! native && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true)) { 
    std::string xmlFileName = new_list.get<std::string>("input file name");
    std::string new_extension("_native_v2.xml");
    size_t pos = xmlFileName.find(".xml");
    xmlFileName.replace(pos, (size_t)4, new_extension, (size_t)0, (size_t)14);
    if (comm->MyPID() == 0) {
        printf("Amanzi: writing the translated parameter list to file %s.\n", xmlFileName.c_str());
    
	
	Teuchos::Amanzi_XMLParameterListWriter XMLWriter;
	Teuchos::XMLObject XMLobj = XMLWriter.toXML(new_list);

	ofstream xmlfile;
	xmlfile.open(xmlFileName.c_str());
	xmlfile << XMLobj;
    }

    //Teuchos::writeParameterListToXmlFile(new_list, xmlFileName);
  }


  //------------ DOMAIN, GEOMETRIC MODEL, ETC ----------------------------

  // Create the simulation domain
  Amanzi::timer_manager.add("Geometric Model creation",Amanzi::Timer::ONCE);
  Amanzi::timer_manager.start("Geometric Model creation");

  Teuchos::ParameterList domain_params = new_list.sublist("Domain");
  unsigned int spdim = domain_params.get<int>("Spatial Dimension");
  
  Amanzi::AmanziGeometry::Domain *simdomain_ptr = new Amanzi::AmanziGeometry::Domain(spdim);


  // Parse the domain description and create.

  // Under the simulation domain we have different geometric
  // models. We can also have free geometric regions not associated
  // with a geometric model.


  // For now create one geometric model from all the regions in the spec

  Teuchos::ParameterList reg_params = new_list.sublist("Regions");

  Amanzi::AmanziGeometry::GeometricModelPtr 
      geom_model_ptr( new Amanzi::AmanziGeometry::GeometricModel(spdim, reg_params, comm) );


  // Add the geometric model to the domain

  simdomain_ptr->Add_Geometric_Model(geom_model_ptr);

  Amanzi::timer_manager.stop("Geometric Model creation");

  // If we had geometric models and free regions coexisting then we would 
  // create the free regions here and add them to the simulation domain
  // Nothing to do for now


  // ---------------- MESH -----------------------------------------------

  Amanzi::timer_manager.add("Mesh creation",Amanzi::Timer::ONCE);
  Amanzi::timer_manager.start("Mesh creation");

  // Create a mesh factory for this geometric model
  Amanzi::AmanziMesh::MeshFactory factory(comm) ;

  // Prepare to read/create the mesh specification
   Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

  // get the Mesh sublist
  ierr = 0;
  Teuchos::ParameterList mesh_params = new_list.sublist("Mesh");


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
	std::string framework = expert_mesh_params.get<string>("Framework");

	if (framework == Amanzi::AmanziMesh::framework_name(Amanzi::AmanziMesh::Simple)) {
	  prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::Simple);
	} else if (framework == Amanzi::AmanziMesh::framework_name(Amanzi::AmanziMesh::MSTK)) {
	  prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::MSTK);
	} else if (framework == Amanzi::AmanziMesh::framework_name(Amanzi::AmanziMesh::STKMESH)) {
	  prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::STKMESH);
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
        std::cerr << "Verifying mesh with Mesh Audit..." << std::endl;
        if (size == 1) {
          Amanzi::MeshAudit mesh_auditor(mesh);
          int status = mesh_auditor.Verify();
          if (status == 0) {
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

  // -------------- MULTI-PROCESS COORDINATOR------- --------------------

  Amanzi::MPC mpc(new_list, mesh, comm, output_observations);

  //--------------- DO THE SIMULATION -----------------------------------

  mpc.cycle_driver();

  //-----------------------------------------------------
  

  // Clean up
  mesh.reset();
  delete comm;
      
  return Amanzi::Simulator::SUCCESS;
}



/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $AMANZI_DIR/COPYRIGHT
Author: ??

Effectively stolen from Amanzi, with few modifications.
------------------------------------------------------------------------- */

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include "Exodus_readers.hh"
#include "Parallel_Exodus_file.hh"

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "MeshAudit.hh"
#include "MeshFactory.hh"
#include "Domain.hh"
#include "GeometricModel.hh"
#include "coordinator.hh"
#include "State.hh"

#include "errors.hh"
#include "exceptions.hh"

#include "ats_clm_driver.hh"
#include "InputParserIS.hh"

#include "global_verbosity.hh"
#include "VerboseObject_objs.hh"

#include "state_evaluators_registration.hh"
#include "constitutive_relations_eos_registration.hh"
#include "constitutive_relations_surface_subsurface_fluxes_registration.hh"
#include "flow_constitutive_relations_overland_conductivity_registration.hh"
#include "flow_constitutive_relations_porosity_registration.hh"
#include "flow_constitutive_relations_wrm_models_registration.hh"
#include "flow_constitutive_relations_wrm_registration.hh"
#include "flow_icy_overland_registration.hh"
#include "flow_overland_head_registration.hh"
#include "flow_overland_registration.hh"
#include "flow_permafrost_registration.hh"
#include "flow_richards_registration.hh"
#include "deform_constitutive_relations_porosity_registration.hh"
#include "deform_prescribed_deformation_registration.hh"
#include "energy_advection_diffusion_registration.hh"
#include "energy_constant_temperature_registration.hh"
#include "energy_constitutive_relations_internal_energy_registration.hh"
#include "energy_constitutive_relations_source_terms_registration.hh"
#include "energy_constitutive_relations_thermal_conductivity_registration.hh"
#include "energy_surface_ice_registration.hh"
#include "energy_two_phase_registration.hh"
#include "energy_three_phase_registration.hh"
#include "surface_balance_SEB_registration.hh"
#include "test_pks_divgrad_test_registration.hh"
// #include "transport_passive_tracer_registration.hh"
#include "mpc_registration.hh"


namespace Amanzi {

Teuchos::EVerbosityLevel Amanzi::VerbosityLevel::level_ = Teuchos::VERB_MEDIUM;

int32_t ATSCLMDriver::Initialize(const MPI_Comm& mpi_comm,
        int* col_types, int num_cols, int num_types) {

  using Teuchos::OSTab;
  setDefaultVerbLevel(VerbosityLevel::level_);
  Teuchos::EVerbosityLevel verbLevel = getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = getOStream();
  OSTab tab = getOSTab(); // This sets the line prefix and adds one tab

  // ======  SET UP THE COMMUNICATOR =======
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(mpi_comm);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  int rank;
  MPI_Comm_rank(mpi_comm,&rank);
  int size;
  MPI_Comm_size(mpi_comm,&size);

  // Read environment variable for XML Input File.
  char * xmlfile = getenv("ATS_XML_INPUT");
  ASSERT(xmlfile != NULL);

  std::cout << "Initing ATS with " << num_cols << " columns" << std::endl;

  // ======  SET UP THE INPUT SPEC =======
  // read the main parameter list
  Teuchos::ParameterList input_parameter_list;
  Teuchos::ParameterXMLFileReader xmlreader(xmlfile);
  input_parameter_list = xmlreader.getParameters();
  Teuchos::RCP<Teuchos::FancyOStream> fos;
  Teuchos::readVerboseObjectSublist(&input_parameter_list,&fos,&VerbosityLevel::level_);

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

  // ======  SET UP THE DOMAIN AND GEOMETRIC MODEL =======
  Teuchos::ParameterList domain_params = params_copy.sublist("Domain");
  unsigned int spdim = domain_params.get<int>("Spatial Dimension");
  simdomain_ = Teuchos::rcp(new AmanziGeometry::Domain(spdim));

  // Under the simulation domain we have different geometric
  // models. We can also have free geometric regions not associated
  // with a geometric model.

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList reg_params = params_copy.sublist("Regions");

  geom_model_ = Teuchos::rcp(new AmanziGeometry::GeometricModel(spdim, reg_params, comm) );

  // Add the geometric model to the domain
  simdomain_->Add_Geometric_Model(&*geom_model_);

  // If we had geometric models and free regions coexisting then we would
  // create the free regions here and add them to the simulation domain
  // Nothing to do for now

  // ======  SET UP THE MESH =======
  AmanziMesh::MeshFactory factory(comm);

  // select the mesh framework
  Teuchos::ParameterList mesh_plist = params_copy.sublist("Mesh");

  int ierr = 0;
  try {
    std::string framework = mesh_plist.get<string>("Framework");
    AmanziMesh::FrameworkPreference prefs(factory.preference());
    if (framework == AmanziMesh::framework_name(AmanziMesh::Simple)) {
      prefs.clear(); prefs.push_back(AmanziMesh::Simple);
    } else if (framework == AmanziMesh::framework_name(AmanziMesh::MOAB)) {
      prefs.clear(); prefs.push_back(AmanziMesh::MOAB);
    } else if (framework == AmanziMesh::framework_name(AmanziMesh::STKMESH)) {
      prefs.clear(); prefs.push_back(AmanziMesh::STKMESH);
    } else if (framework == AmanziMesh::framework_name(AmanziMesh::MSTK)) {
      prefs.clear(); prefs.push_back(AmanziMesh::MSTK);
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
    return Simulator::FAIL;
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
        mesh_ = factory.create(file, &*geom_model_);
      } catch (const std::exception& e) {
        std::cerr << rank << ": error: " << e.what() << std::endl;
        ierr++;
      }

      comm->SumAll(&ierr, &aerr, 1);
      if (aerr > 0) return Simulator::FAIL;
    }
  } else if (mesh_plist.isSublist("Generate Mesh")) {
    // try to generate the mesh from data in plist
    Teuchos::ParameterList gen_params = mesh_plist.sublist("Generate Mesh");
    ierr = 0;

    try {
      mesh_ = factory.create(gen_params, &*geom_model_);
    } catch (const std::exception& e) {
      std::cerr << rank << ": error: " << e.what() << std::endl;
      ierr++;
    }

    comm->SumAll(&ierr, &aerr, 1);
    if (aerr > 0) return Simulator::FAIL;

  } else {
    std::cerr << rank << ": error: "
              << "Neither Read nor Generate options specified for mesh" << std::endl;
    throw std::exception();
  }
  ASSERT(!mesh_.is_null());

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
          MeshAudit mesh_auditor(mesh_);
          int status = mesh_auditor.Verify();
          if (status == 0) {
            std::cout << "Mesh Audit confirms that mesh is ok" << std::endl;
          } else {
            std::cout << "Mesh Audit could not verify correctness of mesh" << std::endl;
            return Simulator::FAIL;
          }
        } else {
          std::ostringstream ofile;
          ofile << "mesh_audit_" << std::setfill('0') << std::setw(4) << rank << ".txt";
          std::ofstream ofs(ofile.str().c_str());
          if (rank == 0)
            std::cerr << "Writing Mesh Audit output to " << ofile.str() << ", etc." << std::endl;

          ierr = 0;
          MeshAudit mesh_auditor(mesh_, ofs);
          int status = mesh_auditor.Verify();        // check the mesh
          if (status != 0) ierr++;

          comm->SumAll(&ierr, &aerr, 1);
          if (aerr == 0) {
            std::cerr << "Mesh Audit confirms that mesh is ok" << std::endl;
          } else {
            if (rank == 0)
              std::cerr << "Mesh Audit could not verify correctness of mesh" << std::endl;
            return Simulator::FAIL;
          }
        }
      }  // if verify
    }  // if verify_mesh_param
  }  // If expert_params_specified

  // Create the surface mesh if needed
  Teuchos::RCP<AmanziMesh::Mesh> surface_mesh = Teuchos::null;
  Teuchos::RCP<AmanziMesh::Mesh> surface3D_mesh = Teuchos::null;
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

    if (mesh_->cell_dimension() == 3) {
      surface3D_mesh = factory.create(&*mesh_,setnames,AmanziMesh::FACE,false,false);
      surface_mesh = factory.create(&*mesh_,setnames,AmanziMesh::FACE,true,false);
    } else {
      surface3D_mesh = mesh_;
      surface_mesh = factory.create(&*mesh_,setnames,AmanziMesh::CELL,true,false);
    }
  }

  Teuchos::TimeMonitor::summarize();
  Teuchos::TimeMonitor::zeroOutTimers();

  // ======= SET UP THE STATE AND COORDINATOR =========
  // Create the state.
  Teuchos::ParameterList state_plist = params_copy.sublist("state");
  S_ = Teuchos::rcp(new State(state_plist));
  S_->RegisterDomainMesh(mesh_);
  if (surface3D_mesh != Teuchos::null) S_->RegisterMesh("surface_3d", surface3D_mesh);
  if (surface_mesh != Teuchos::null) S_->RegisterMesh("surface", surface_mesh);

  // slaved fluxes from CLM
  S_->RequireField("surface_total_energy_source","clm")->SetMesh(surface_mesh)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireField("surface_mass_source","clm")->SetMesh(surface_mesh)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireField("surface_mass_source_temperature","clm")->SetMesh(surface_mesh)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S_->FEList().sublist("surface_total_energy_source").set("field evaluator type", "primary variable");
  S_->FEList().sublist("surface_mass_source").set("field evaluator type", "primary variable");
  S_->FEList().sublist("surface_mass_source_temperature").set("field evaluator type", "primary variable");
  S_->RequireFieldEvaluator("surface_total_energy_source");
  S_->RequireFieldEvaluator("surface_mass_source");
  S_->RequireFieldEvaluator("surface_mass_source_temperature");

  // create the top level Coordinator
  coordinator_ = Teuchos::rcp(new Coordinator(params_copy, S_, comm));
  coordinator_->initialize();

  // ======= SET UP THE CLM TO ATS MAPPINGS =========

  // ASSUMES IDENTITY MAP!

  // Create maps for CLM and ATS
  // CURRENTLY ASSUMES that this can be done locally!
  // clm_type_region_names = mesh_plist.get<Teuchos::Array<std::string> >("clm region list").toVector();

  // int nindices = clm_type_region_names.size();
  // ASSERT(nindices == num_types);

  ncells_surf_ = surface_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  //  ASSERT(ncells_surf_ == num_cols); // THIS CURRENTLY FAILS!
  surf_clm_map_ = Teuchos::rcp(new Epetra_Map(-1, ncells_surf_, 0, *comm));

  ncells_sub_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  ASSERT(ncells_sub_ == ncells_surf_ * 15); // MAGIC NUMBER OF CELLS PER COL IN CLM
  sub_clm_map_ = Teuchos::rcp(new Epetra_Map(-1, ncells_sub_, 0, *comm));

  // // create list of region each cell is in
  // std::vector<int> region_ats(ncells_surf_,-1);
  // int counter = 1;
  // for (std::vector<std::string>::const_iterator region=clm_type_region_names.begin();
  //      region!=clm_type_region_names.end(); ++region) {
  //   AmanziMesh::Entity_ID_List surf_cells;
  //   surface_mesh->get_set_entities(*region, AmanziMesh::CELL,
  //           AmanziMesh::OWNED, &surf_cells);
  //   for (AmanziMesh::Entity_ID_List::const_iterator c=surf_cells.begin();
  //        c!=surf_cells.end(); ++c) {
  //     ASSERT(region_ats[*c] < 0);
  //     region_ats[*c] = counter;
  //   }
  //   counter++;
  // }

  // // loop over each index in the clm list, finding a cell of the same type in
  // // the ATS list
  // // -- create a col map for CLM
  // std::vector< std::vector<int>::const_iterator > iters(nindices, region_ats.begin());

  // std::vector<int> lids_surf(ncells_surf_,-1);
  // for (int lcv=0; lcv!=ncells_surf_; ++lcv) {
  //   int my_type = col_types[lcv];

  //   // get the iterator pointing to where we are on this type
  //   std::vector<int>::const_iterator my_iter = iters[my_type-1];

  //   // find the next instance of this type
  //   std::vector<int>::const_iterator my_loc = std::find<std::vector<int>::const_iterator, int>(my_iter, region_ats.end(), my_type);
  //   ASSERT(my_loc != region_ats.end());

  //   // grab the clm col gid at this location
  //   lids_surf[my_loc - region_ats.begin()] = lcv;

  //   // put this pointer back into the array
  //   my_iter++; // increment to point to the NEXT spot
  //   iters[my_type-1] = my_iter;
  // }

  // // Set up the subsurface gids from there
  // std::vector<int> gids_sub(ncells_sub_, -1);

  // // -- loop over surface cells
  // for (int lcv=0; lcv!=ncells_surf_; ++lcv) {
  //   // get the corresponding clm lid
  //   int ind = lids_surf[lcv];

  //   // get the face on the subsurf mesh
  //   AmanziMesh::Entity_ID f =
  //       surface_mesh->entity_get_parent(AmanziMesh::CELL, lcv);

  //   // get the cell within the face
  //   AmanziMesh::Entity_ID_List cells;
  //   mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
  //   ASSERT(cells.size() == 1);
  //   AmanziMesh::Entity_ID c = cells[0];

  //   // loop over cells_below() to find the column
  //   int cell_ind = ind * 15;
  //   while (c >= 0) {
  //     gids_sub[c] = sub_clm_map_->GID(cell_ind);
  //     cell_ind++;
  //     c = mesh_->cell_get_cell_below(c);
  //   }
  //   ASSERT(cell_ind == (ind+1)*15); // checks to make sure we got the right number of cells below
  // }

  // // lids --> gids
  // std::vector<int> gids_surf(ncells_surf_,-1);
  // for (int lcv=0; lcv!=ncells_surf_; ++lcv)
  //   gids_surf[lcv] = surf_clm_map_->GID(lids_surf[lcv]);

  // Create the imports
  const Epetra_Map& ats_cell_map = mesh_->cell_map(false);
  const Epetra_Map& ats_col_map = surface_mesh->cell_map(false);
  //  Epetra_Map ats_imp_col_map(-1, ncells_surf_, &gids_surf[0], 0, *comm);
  //  Epetra_Map ats_imp_cell_map(-1, ncells_surf_, &gids_sub[0], 0, *comm);
  Epetra_Map ats_imp_col_map(ats_col_map);
  Epetra_Map ats_imp_cell_map(ats_cell_map);

  sub_importer_ = Teuchos::rcp(new Epetra_Import(ats_cell_map, ats_imp_cell_map));
  surf_importer_ = Teuchos::rcp(new Epetra_Import(ats_col_map, ats_imp_col_map));
}

int32_t ATSCLMDriver::Finalize() {
  coordinator_->finalize();
  return 0;
}


int32_t ATSCLMDriver::SetData_(std::string key, double* data, int length) {
  Teuchos::RCP<State> S = S_next_ == Teuchos::null ? S_ : S_next_;

  // from the name, grab data from state and determine the type of the data
  Teuchos::RCP<CompositeVector> dat =
      S->GetFieldData(key, S->GetField(key)->owner());

  // Create an Epetra_Vector from the data and the appropriate map
  Epetra_Vector dat_clm(View, *sub_clm_map_, data);

  // Call the importer to map to state's vector
  Epetra_MultiVector& dat_v = *dat->ViewComponent("cell",false);
  int ierr = dat_v.Import(dat_clm, *sub_importer_, Insert);
  ASSERT(!ierr);
  return ierr;

}


int32_t ATSCLMDriver::GetData_(std::string key, double* data, int length) {
  Teuchos::RCP<State> S = S_next_ == Teuchos::null ? S_ : S_next_;

  // from the name, grab data from state and determine the type of the data
  Teuchos::RCP<CompositeVector> dat =
      S->GetFieldData(key, S->GetField(key)->owner());

  // Create an Epetra_Vector from the data and the appropriate map
  Epetra_Vector dat_clm(View, *sub_clm_map_, data);

  // Call the importer to map to state's vector
  Epetra_MultiVector& dat_v = *dat->ViewComponent("cell",false);
  int ierr = dat_v.Export(dat_clm, *sub_importer_, Insert);
  ASSERT(!ierr);
  return ierr;
}

int32_t ATSCLMDriver::SetSurfaceData_(std::string key, double* data, int length) {
  Teuchos::RCP<State> S = S_next_ == Teuchos::null ? S_ : S_next_;

  // from the name, grab data from state and determine the type of the data
  Teuchos::RCP<CompositeVector> dat =
      S->GetFieldData(key, S->GetField(key)->owner());

  // Create an Epetra_Vector from the data and the appropriate map
  Epetra_Vector dat_clm(View, *surf_clm_map_, data);

  // Call the importer to map to state's vector
  Epetra_MultiVector& dat_v = *dat->ViewComponent("cell",false);
  int ierr = dat_v.Import(dat_clm, *surf_importer_, Insert);
  ASSERT(!ierr);
  return ierr;

}

int32_t ATSCLMDriver::GetSurfaceData_(std::string key, double* data, int length) {
  Teuchos::RCP<State> S = S_next_ == Teuchos::null ? S_ : S_next_;

  // from the name, grab data from state and determine the type of the data
  Teuchos::RCP<CompositeVector> dat =
      S->GetFieldData(key, S->GetField(key)->owner());

  // Create an Epetra_Vector from the data and the appropriate map
  Epetra_Vector dat_clm(View, *surf_clm_map_, data);

  // Call the importer to map to state's vector
  Epetra_MultiVector& dat_v = *dat->ViewComponent("cell",false);
  int ierr = dat_v.Export(dat_clm, *surf_importer_, Insert);
  ASSERT(!ierr);
  return ierr;
}

int32_t ATSCLMDriver::SetInitCLMData(double* T, double* sl, double* si) {
  int ierr(0);
  // ierr |= SetData_("temperature", T, ncells_sub_);
  // S_->GetField("temperature", S_->GetField("temperature")->owner())->set_initialized();
  // ierr |= SetData_("saturation_liquid", sl, ncells_sub_);
  // S_->GetField("saturation_liquid", S_->GetField("saturation_liquid")->owner())->set_initialized();
  // ierr |= SetData_("saturation_ice", si, ncells_sub_);
  // S_->GetField("saturation_ice", S_->GetField("saturation_ice")->owner())->set_initialized();
  return ierr;
}

int32_t ATSCLMDriver::SetCLMData(double* e_flux, double* w_flux) {
  int ierr(0);
  ierr |= SetSurfaceData_("surface_total_energy_source", e_flux, ncells_surf_);
  ierr |= SetSurfaceData_("surface_mass_source", w_flux, ncells_surf_);
  return ierr;
}


int32_t ATSCLMDriver::GetCLMData(double* T, double* sl, double* si) {
  int ierr(0);
  ierr |= GetData_("temperature", T, ncells_sub_);
  ierr |= GetData_("saturation_liquid", sl, ncells_sub_);
  ierr |= GetData_("saturation_ice", si, ncells_sub_);
  return ierr;
}


int32_t ATSCLMDriver::Advance(double dt, bool force_vis) {
  int ierr(0);

  int nsteps_max = 100;
  double dt_min = 1.e-10;
  int istep = 0;
  double t0 = S_->time();
  double t1 = t0 + dt;
  double dt_local = dt;

  verbosity_ = getVerbLevel();

  bool done = false;
  Teuchos::RCP<Teuchos::FancyOStream> out = getOStream();
  while (!done) {
    if (out.get() && includesVerbLevel(verbosity_, Teuchos::VERB_MEDIUM, true)) {
      Teuchos::OSTab tab = getOSTab();
      *out << "======================================================================"
            << std::endl << std::endl;
      *out << "Cycle = " << S_->cycle();
      *out << ",  Time [days] = "<< S_->time() / (60*60*24);
      *out << ",  dt [days] = " << dt_local / (60*60*24)  << std::endl;
      *out << "----------------------------------------------------------------------"
            << std::endl;
    }

    bool fail = coordinator_->advance(dt);
    dt_local = coordinator_->get_dt();
    done = std::abs(S_next_->time() - t1) < 1.e-10 || istep > nsteps_max || dt_local < dt_min;
    istep++;
  }

  return !(std::abs(S_next_->time() - t1) < 1.e-10);
}


} // namespace Amanzi

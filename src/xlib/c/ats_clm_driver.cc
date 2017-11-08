/* -*-  mode: c++; indent-tabs-mode: nil -*- */
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
#include "test_pks_registration.hh"
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
  Epetra_MpiComm *comm = new Epetra_MpiComm(mpi_comm);

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
    std::string framework = mesh_plist.get<std::string>("Framework");
    AmanziMesh::FrameworkPreference prefs(factory.preference());
    if (framework == AmanziMesh::framework_name(AmanziMesh::MSTK)) {
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
  state_plist.print(std::cout);
  S_->RegisterDomainMesh(mesh_);
  if (surface3D_mesh != Teuchos::null) S_->RegisterMesh("surface_3d", surface3D_mesh);
  if (surface_mesh != Teuchos::null) S_->RegisterMesh("surface", surface_mesh);

  // set up primary vars
  ASSERT(!S_->FEList().isSublist("surface_total_energy_source"));
  ASSERT(!S_->FEList().isSublist("surface_mass_source"));
  S_->FEList().sublist("surface_total_energy_source").set("field evaluator type", "primary variable");
  S_->FEList().sublist("surface_mass_source").set("field evaluator type", "primary variable");

  coordinator_ = Teuchos::rcp(new Coordinator(params_copy, S_, comm));
  coord_setup_ = false;
  coord_init_ = false;

  // slaved fluxes from CLM
  S_->RequireFieldEvaluator("surface_total_energy_source");
  S_->RequireFieldEvaluator("surface_mass_source");

  S_->RequireField("surface_total_energy_source","clm")->SetMesh(surface_mesh)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireField("surface_mass_source","clm")->SetMesh(surface_mesh)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

  // ======= SET UP THE CLM TO ATS MAPPINGS =========
  // Currently assumes the identity map on the surface.
  // Currently assumes this can be done locally.

  // -- identity map for surface
  ncells_surf_ = surface_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  ASSERT(ncells_surf_ == num_cols);
  surf_clm_map_ = Teuchos::rcp(new Epetra_Map(-1, ncells_surf_, 0, *comm));

  const Epetra_Map& ats_col_map = surface_mesh->cell_map(false);
  surf_importer_ = Teuchos::rcp(new Epetra_Import(ats_col_map, *surf_clm_map_));

  // -- subsurface map -- top down
  ncells_sub_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  ASSERT(ncells_sub_ == ncells_surf_ * 15); // MAGIC NUMBER OF CELLS PER COL IN CLM
  sub_clm_map_ = Teuchos::rcp(new Epetra_Map(-1, ncells_sub_, 0, *comm));

  // Set up the subsurface gids -- GIDs in sub_clm_map
  const Epetra_Map& ats_cell_map = mesh_->cell_map(false);
  std::vector<int> gids_sub(ncells_sub_, -1);

  // -- loop over surface cells
  for (int icol=0; icol!=ncells_surf_; ++icol) {
    // get the GID of ATS's icol
    int ats_col_gid = ats_col_map.GID(icol);

    // get the LID of the corresponding GID on CLM (assumes map is local permutation, but does not assume identity map)
    int clm_col_lid = surf_clm_map_->LID(ats_col_gid);
    ASSERT(clm_col_lid >= 0);

    // get the face on the subsurf mesh
    AmanziMesh::Entity_ID f =
        surface_mesh->entity_get_parent(AmanziMesh::CELL, icol);

    // get the cell within the face
    AmanziMesh::Entity_ID_List cells;
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    ASSERT(cells.size() == 1);
    AmanziMesh::Entity_ID c = cells[0];

    // loop over cells_below() to find the column.  NOTE: CLM uses horizontal
    // cells as fastest varying, NOT vertical columns.
    int clm_cell_lid = clm_col_lid;
    int ncells_in_col = 0;
    while (c >= 0) {
      gids_sub[c] = sub_clm_map_->GID(clm_cell_lid);
      clm_cell_lid += num_cols;
      c = mesh_->cell_get_cell_below(c);
      ncells_in_col++;
    }
    ASSERT(ncells_in_col == 15);
  }

  // Create the imports
  ASSERT(*std::min_element(gids_sub.begin(), gids_sub.end()) >= 0);
  Epetra_Map ats_cell_gids(-1, ncells_sub_, &gids_sub[0], 0, *comm);
  sub_importer_ = Teuchos::rcp(new Epetra_Import(ats_cell_gids, *sub_clm_map_));

  // set up the coordinator, allocating space
  coordinator_->setup();
  coord_setup_ = true;
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

  std::cout << "SetData (ATS):" << std::endl;
  std::cout << "data[0] = " << dat_v[0][0] << std::endl;
  std::cout << "data[1] = " << dat_v[0][1] << std::endl;

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
  if (S_next_ == Teuchos::null)
    std::cout << "writing to S_old" << std::endl;
  else
    std::cout << "writing to S_next" << std::endl;

  // from the name, grab data from state and determine the type of the data
  Teuchos::RCP<CompositeVector> dat =
      S->GetFieldData(key, S->GetField(key)->owner());

  // Create an Epetra_Vector from the data and the appropriate map
  Epetra_Vector dat_clm(View, *surf_clm_map_, data);

  // Call the importer to map to state's vector
  Epetra_MultiVector& dat_v = *dat->ViewComponent("cell",false);
  int ierr = dat_v.Import(dat_clm, *surf_importer_, Insert);
  ASSERT(!ierr);

  std::cout << "Surf data in = " << std::endl;
  dat_clm.Print(std::cout);
  std::cout << "Surf data out = " << std::endl;
  dat_v.Print(std::cout);

  //  std::cout << "Surf Data In (" << key << ")[0] = " << data[0] << std::endl;
  //  std::cout << "Surf Data Out (" << key << ")[0] = " << dat_v[0][0] << std::endl;
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
  std::cout << "SetInitCLMData (ATS):" << std::endl;
  std::cout << "T_ats[0] = " << T[0] << std::endl;
  std::cout << "T_ats[374] = " << T[374] << std::endl;
  std::cout << "sl_ats[0] = " << sl[0] << std::endl;
  std::cout << "sl_ats[374] = " << sl[374] << std::endl;
  std::cout << "si_ats[0] = " << si[0] << std::endl;
  std::cout << "si_ats[374] = " << si[374] << std::endl;

  ierr |= SetData_("temperature", T, ncells_sub_);

  // ierr |= SetData_("saturation_liquid", sl, ncells_sub_);
  // S_->GetField("saturation_liquid", S_->GetField("saturation_liquid")->owner())->set_initialized();
  // ierr |= SetData_("saturation_ice", si, ncells_sub_);
  // S_->GetField("saturation_ice", S_->GetField("saturation_ice")->owner())->set_initialized();

  coordinator_->initialize();
  coord_init_ = true;
  S_next_ = coordinator_->get_next_state();

  return ierr;
}

int32_t ATSCLMDriver::SetCLMData(double* e_flux, double* w_flux) {
  int ierr(0);

  std::cout << "SetCLMData (ATS) (on " << ncells_surf_ << " cells):" << std::endl;
  std::cout << "Qe_ats[0] = " << e_flux[0] << std::endl;
  std::cout << "Qe_ats[24] = " << e_flux[24] << std::endl;
  std::cout << "Qw_ats[0] = " << w_flux[0] << std::endl;
  std::cout << "Qw_ats[24] = " << w_flux[24] << std::endl;

  ierr |= SetSurfaceData_("surface_total_energy_source", e_flux, ncells_surf_);
  ierr |= SetSurfaceData_("surface_mass_source", w_flux, ncells_surf_);
  return ierr;
}


int32_t ATSCLMDriver::GetCLMData(double* T, double* sl, double* si) {
  int ierr(0);

  ierr |= GetData_("temperature", T, ncells_sub_);
  ierr |= GetData_("saturation_liquid", sl, ncells_sub_);
  ierr |= GetData_("saturation_ice", si, ncells_sub_);

  std::cout << "GetCLMData (ATS):" << std::endl;
  std::cout << "T_ats[0] = " << T[0] << std::endl;
  std::cout << "T_ats[374] = " << T[374] << std::endl;
  std::cout << "sl_ats[0] = " << sl[0] << std::endl;
  std::cout << "sl_ats[374] = " << sl[374] << std::endl;
  std::cout << "si_ats[0] = " << si[0] << std::endl;
  std::cout << "si_ats[374] = " << si[374] << std::endl;

  return ierr;
}


int32_t ATSCLMDriver::Advance(double dt, bool force_vis) {
  ASSERT(coord_setup_);
  ASSERT(coord_init_);

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

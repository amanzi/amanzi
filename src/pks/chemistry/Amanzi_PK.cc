/*
  Chemistry PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors:

  Trilinos based process kernel for chemistry. Geochemistry
  calculations live in the chemistry library. The PK stores the 
  instance of the chemistry object and drives the chemistry 
  calculations on a cell by cell basis. It handles the movement of
  data back and forth between the State and the chemistry library 
  data structures.
*/
 
#include <string>
#include <algorithm>

// TPLs
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "Mesh.hh"
#include "Beaker.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "SimpleThermoDatabase.hh"
#include "VerboseObject.hh"

// Chemistry
#include "Amanzi_PK.hh"

namespace Amanzi {
namespace AmanziChemistry {

/* ******************************************************************
* Constructor
******************************************************************* */
Amanzi_PK::Amanzi_PK(Teuchos::ParameterList& pk_tree,
                     const Teuchos::RCP<Teuchos::ParameterList>& glist,
                     const Teuchos::RCP<State>& S,
                     const Teuchos::RCP<TreeVector>& soln) :
    soln_(soln),
    dt_max_(9.9e9),
    chem_(NULL),
    current_time_(0.0),
    saved_time_(0.0)
{
  S_ = S;
  // mesh_ = S_->GetMesh();
  glist_ = glist;

  // extract pk name
  std::string pk_name = pk_tree.name();
  auto found = pk_name.rfind("->");
  if (found != std::string::npos) pk_name.erase(0, found + 2);

  // create pointer to the chemistry parameter list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  cp_list_ = Teuchos::sublist(pk_list, pk_name, true);
  domain_ = cp_list_->get<std::string>("domain name", "domain");

  // obtain key of fields
  tcc_key_ = Keys::getKey(domain_, "total_component_concentration"); 
  poro_key_ = cp_list_->get<std::string>("porosity key", Keys::getKey(domain_, "porosity"));
  saturation_key_ = cp_list_->get<std::string>("saturation key", Keys::getKey(domain_, "saturation_liquid"));
  fluid_den_key_ = cp_list_->get<std::string>("fluid density key", Keys::getKey(domain_, "mass_density_liquid"));

  min_vol_frac_key_ = Keys::getKey(domain_,"mineral_volume_fractions");
  min_ssa_key_ = Keys::getKey(domain_,"mineral_specific_surface_area");
  sorp_sites_key_ = Keys::getKey(domain_,"sorption_sites");
  surf_cfsc_key_ = Keys::getKey(domain_,"surface_complex_free_site_conc");
  total_sorbed_key_ = Keys::getKey(domain_,"total_sorbed");
  isotherm_kd_key_ = Keys::getKey(domain_,"isotherm_kd");
  isotherm_freundlich_n_key_ = Keys::getKey(domain_,"isotherm_freundlich_n");
  isotherm_langmuir_b_key_ = Keys::getKey(domain_,"isotherm_langmuir_b");
  free_ion_species_key_ = Keys::getKey(domain_,"free_ion_species");
  primary_activity_coeff_key_ = Keys::getKey(domain_,"primary_activity_coeff");

  ion_exchange_sites_key_ = Keys::getKey(domain_,"ion_exchange_sites");
  //ion_exchange_sites_key_ = "ion_exchange_sites";

  ion_exchange_ref_cation_conc_key_ = Keys::getKey(domain_,"ion_exchange_ref_cation_conc");
  secondary_activity_coeff_key_ = Keys::getKey(domain_,"secondary_activity_coeff");
  alquimia_aux_data_key_ = Keys::getKey(domain_,"alquimia_aux_data");
  mineral_rate_constant_key_ = Keys::getKey(domain_,"mineral_rate_constant");
  first_order_decay_constant_key_ = Keys::getKey(domain_,"first_order_decay_constant");  
  
  // collect high-level information about the problem
  Teuchos::RCP<Teuchos::ParameterList> state_list = Teuchos::sublist(glist, "state", true);

  InitializeMinerals(cp_list_);
  InitializeSorptionSites(cp_list_, Teuchos::sublist(state_list, "initial conditions"));

  // grab the component names
  comp_names_.clear();
  Teuchos::RCP<Teuchos::ParameterList> cd_list = Teuchos::sublist(glist, "cycle driver", true);
  if (cd_list->isParameter("component names")) {
    comp_names_ = cd_list->get<Teuchos::Array<std::string> >("component names").toVector();
  } else{
    Errors::Message msg("Amanzi_PK: Cycle Driver has no input parameter component names.");
    Exceptions::amanzi_throw(msg);
  }
  number_aqueous_components_ = comp_names_.size();
  number_free_ion_ = number_aqueous_components_;
  number_total_sorbed_ = number_aqueous_components_;

  // verbosity object
  vo_ = Teuchos::rcp(new VerboseObject("Amanzi_PK:" + domain_, *cp_list_)); 
}


/* *******************************************************************
* Destructor
******************************************************************* */
Amanzi_PK::~Amanzi_PK() {
  delete chem_;
}


/* ******************************************************************
* Register fields and evaluators with the State
******************************************************************* */
void Amanzi_PK::Setup(const Teuchos::Ptr<State>& S)
{
  // use common registration steps
  Chemistry_PK::Setup(S);

  // no additional steps to be done.
}


/* ******************************************************************
* Can this be done during Setup phase?
******************************************************************* */
void Amanzi_PK::AllocateAdditionalChemistryStorage_(
    const Beaker::BeakerState& state)
{
  int n_secondary_comps = state.secondary_activity_coeff.size();
  if (n_secondary_comps > 0) {
    Teuchos::RCP<CompositeVectorSpace> fac = S_->RequireField(secondary_activity_coeff_key_, passwd_);
    fac->SetMesh(mesh_)->SetGhosted(false)
       ->SetComponent("cell", AmanziMesh::CELL, n_secondary_comps);

    Teuchos::RCP<CompositeVector> sac = Teuchos::rcp(new CompositeVector(*fac));
    S_->GetField(secondary_activity_coeff_key_, passwd_)->SetData(sac);
    S_->GetField(secondary_activity_coeff_key_, passwd_)->CreateData();
    S_->GetFieldData(secondary_activity_coeff_key_, passwd_)->PutScalar(1.0);
    S_->GetField(secondary_activity_coeff_key_, passwd_)->set_initialized();
  }
}


/* *******************************************************************
* Initialization
******************************************************************* */
void Amanzi_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  // initialization using base class
  Chemistry_PK::Initialize(S);

  Teuchos::RCP<Epetra_MultiVector> tcc = S->GetFieldData(tcc_key_, passwd_)->ViewComponent("cell", true);

  XMLParameters();

  // TODO: some sort of check of the state object to see if mineral_ssa,
  // CEC, site density, etc is present.

  // initial conditions for minerals etc should be handled by the
  // state/chemistry_state object before we reach this point. We just
  // resize our local memory for migrating data here.

  SizeBeakerState_();
  CopyCellStateToBeakerState(0, tcc);
  chem_->Initialize(beaker_parameters_);
  chem_->CopyStateToBeaker(beaker_state_);
  // chem_->VerifyState(beaker_state_);

  // copy the cell data into the beaker storage for initialization purposes
  // but first ensure dependencies are filled
  S->GetFieldEvaluator(poro_key_)->HasFieldChanged(S_.ptr(), name_);
  S->GetFieldEvaluator(fluid_den_key_)->HasFieldChanged(S_.ptr(), name_);
  S->GetFieldEvaluator(saturation_key_)->HasFieldChanged(S_.ptr(), name_);

  // finish setting up & testing the chemistry object
  int nprimary, ierr(0);
  std::string internal_msg;
  std::vector<double> values;

  const auto& iclist = glist_->sublist("state").sublist("initial conditions").sublist(tcc_key_);
  auto constraints = iclist.get<Teuchos::Array<std::string> >("names").toVector();

  try {
    vo_->Write(Teuchos::VERB_HIGH, "Initializing chemistry in cell 0...\n");
    chem_->Display();

    // check names of primary species
    nprimary = chem_->primary_species().size(); 
    if (nprimary == comp_names_.size()) {
      for (int i = 0; i < nprimary; ++i) {
        std::string species_name = chem_->primary_species().at(i).name();
        if (comp_names_[i] != species_name) {
          Errors::Message msg;
          msg << "Amanzi PK: mismatch of name: \"" << comp_names_[i] << "\" and \"" 
              << species_name << "\". Compare XML and BGD lists.";
          Exceptions::amanzi_throw(msg);
        }
      }
    }

    // solve for initial free-ion concentrations
    vo_->Write(Teuchos::VERB_HIGH, "Initial speciation calculations in cell 0...\n");

    for (int i = 0; i < nprimary; ++i) values.push_back((*tcc)[i][0]);
    chem_->EnforceConstraint(&beaker_state_, beaker_parameters_, constraints, values);
    // chem_->Speciate(&beaker_state_, beaker_parameters_);

    vo_->Write(Teuchos::VERB_HIGH, "\nTest solution of initial conditions in cell 0:\n");
    chem_->DisplayResults();
  } catch (Exceptions::Amanzi_exception& geochem_err) {
    ierr = 1;
    internal_msg = geochem_err.what();
  }

  ErrorAnalysis(ierr, internal_msg);

  // TODO(bandre): at this point we should know about any additional
  // storage that chemistry needs...
  AllocateAdditionalChemistryStorage_(beaker_state_);

  SetupAuxiliaryOutput();

  // solve for initial free-ion concentrations
  int num_cells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ierr = 0;
  if (fabs(initial_conditions_time_ - S->time()) <= 1e-8 * fabs(S->time())) {
    for (int c = 0; c < num_cells; ++c) {
      CopyCellStateToBeakerState(c, tcc);

      try {
        for (int i = 0; i < nprimary; ++i) values[i] = (*tcc)[i][c];
        chem_->EnforceConstraint(&beaker_state_, beaker_parameters_, constraints, values);
        // chem_->Speciate(&beaker_state_, beaker_parameters_);
        CopyBeakerStructuresToCellState(c, tcc);
      } 
      catch (Exceptions::Amanzi_exception& geochem_err) {
        ierr = 1;
        internal_msg = geochem_err.what();
      }
    }
  } else {
    if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "no data initialization due to time mismatch: " << S->time() << std::endl;
    }
  }

  // figure out if any of the processes threw an error, if so all processes will re-throw
  ErrorAnalysis(ierr, internal_msg);

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << vo_->color("green") << "Initialization of PK was successful, T=" 
        << S->time() << vo_->reset() << std::endl << std::endl;
  }
}


/* *******************************************************************
* Initialization helper functions
******************************************************************* */
void Amanzi_PK::XMLParameters()
{
  // thermo file name and format, then create the database!
  if (cp_list_->isSublist("thermodynamic database")) {
    Teuchos::RCP<Teuchos::ParameterList> tdb_list = Teuchos::sublist(glist_, "thermodynamic database", true);

    // currently we only support the simple format.
    chem_ = new SimpleThermoDatabase(tdb_list, vo_);
    beaker_parameters_.tolerance = 1e-12;
    beaker_parameters_.max_iterations = 250;
    beaker_parameters_.activity_model_name = "unit";

  } else {
    std::ostringstream msg;
    msg << "Amanzi_PK: 'thermodynamic database' sublist must be specified.\n";
    Exceptions::amanzi_throw(Errors::Message(msg.str()));    
  }

  // activity model
  beaker_parameters_.activity_model_name = cp_list_->get<std::string>("activity model", "unit");
  // -- Pitzer virial coefficients database
  if (beaker_parameters_.activity_model_name == "pitzer-hwm") {
    if (cp_list_->isParameter("Pitzer database file")) {
      beaker_parameters_.pitzer_database = cp_list_->get<std::string>("Pitzer database file");
    } else {
      std::ostringstream msg;
      msg << "Amanzi_PK: parameter 'Pitzer database file' must be specified if activity model=pitzer-hwm'.\n";
      Exceptions::amanzi_throw(Errors::Message(msg.str()));
    }
  }

  // solver parameters
  beaker_parameters_.tolerance = cp_list_->get<double>("tolerance", 1.0e-12);
  beaker_parameters_.max_iterations = cp_list_->get<int>("maximum Newton iterations", 200);

  // auxiliary data
  aux_names_.clear();
  if (cp_list_->isParameter("auxiliary data")) {
    Teuchos::Array<std::string> names = cp_list_->get<Teuchos::Array<std::string> >("auxiliary data");
    for (auto name = names.begin(); name != names.end(); ++name) {
      if (*name == "pH") {
        aux_names_.push_back(*name);
      } else {
        std::stringstream message;
        message << "XMLParameters(): unknown value in 'auxiliary data' list: " 
                << *name << std::endl;
        vo_->WriteWarning(Teuchos::VERB_LOW, message);
      }
    }
  }

  // misc other chemistry flags
  dt_control_method_ = cp_list_->get<std::string>("time step control method", "fixed");
  dt_max_ = cp_list_->get<double>("max time step (s)", 9.9e+9);
  dt_next_ = cp_list_->get<double>("initial time step (s)", dt_max_);
  dt_cut_factor_ = cp_list_->get<double>("time step cut factor", 2.0);
  dt_increase_factor_ = cp_list_->get<double>("time step increase factor", 1.2);

  dt_cut_threshold_ = cp_list_->get<int>("time step cut threshold", 8);
  dt_increase_threshold_ = cp_list_->get<int>("time step increase threshold", 4);

  num_successful_steps_ = 0;
}


/* *******************************************************************
* Requires that Beaker::Setup() has already been called!
******************************************************************* */
void Amanzi_PK::SetupAuxiliaryOutput()
{
  // TODO: this indexing scheme will not be appropriate when
  // additional types of aux data are requested, e.g. mineral SI.....
  unsigned int nvars = aux_names_.size();
  std::string name;
  aux_index_.clear();
  for (unsigned int i = 0; i < nvars; i++) {
    if (aux_names_.at(i) == "pH") {
      name = "H+";
    } else {
      name = aux_names_.at(i);
    }
    int index = chem_->GetPrimaryIndex(name);
    if (index == -1) {
      // check to make sure it is not -1, an invalid name/index
      std::stringstream message;
      message << "ChemistryPK: Output was requested for '" << aux_names_.at(i) 
              << "' (" << name 
              << ") but no chemistry varibles of this name were found.\n";
      Exceptions::amanzi_throw(Errors::Message(message.str()));
    } else {
      aux_index_.push_back(index);
    }
  }

  // create the Epetra_MultiVector that will hold the data
  if (nvars > 0) {
    aux_data_ = Teuchos::rcp(new Epetra_MultiVector(mesh_->cell_map(false), nvars));
  } else {
    aux_data_ = Teuchos::null;
  }
}


/* *******************************************************************
* Initialize the beaker component data structure
******************************************************************* */
void Amanzi_PK::SizeBeakerState_()
{
  // NOTE: The beaker already has data for site density, sorption
  // isotherms, ssa. If we want to use that single global value, then
  // we leave these arrays empty as a flag to the beaker to use its
  // own value. If we want to override the global chemistry value
  // with cell by cell data, then we resize the containers here.

  beaker_state_.total.resize(number_aqueous_components_, 0.0);
  beaker_state_.free_ion.resize(number_aqueous_components_, 1.0e-9);
  beaker_state_.mineral_volume_fraction.resize(number_minerals_, 0.0);

  if (using_sorption_) {
    beaker_state_.total_sorbed.resize(number_total_sorbed_, 0.0);
  } else {
    beaker_state_.total_sorbed.clear();
  }

  if (number_minerals_ > 0) {
    beaker_state_.mineral_specific_surface_area.resize(number_minerals_, 0.0);
  } else {
    beaker_state_.mineral_specific_surface_area.clear();
  }

  if (number_ion_exchange_sites_ > 0) {
    beaker_state_.ion_exchange_sites.resize(number_ion_exchange_sites_, 0.0);
  } else {
    beaker_state_.ion_exchange_sites.clear();
  }

  if (number_sorption_sites_ > 0) {
    beaker_state_.surface_site_density.resize(number_sorption_sites_, 0.0);
  } else {
    beaker_state_.surface_site_density.clear();
  }

  if (using_sorption_isotherms_) {
    beaker_state_.isotherm_kd.resize(number_aqueous_components_, 0.0);
    beaker_state_.isotherm_freundlich_n.resize(number_aqueous_components_, 0.0);
    beaker_state_.isotherm_langmuir_b.resize(number_aqueous_components_, 0.0);
  } else {
    beaker_state_.isotherm_kd.clear();
    beaker_state_.isotherm_freundlich_n.clear();
    beaker_state_.isotherm_langmuir_b.clear();
  }
}


/* *******************************************************************
* We must use the aqueous totals value calculated from transport
* (aqueous_components), not the value stored in state!
******************************************************************* */
void Amanzi_PK::CopyCellStateToBeakerState(
    int c, Teuchos::RCP<Epetra_MultiVector> aqueous_components)
{
  for (unsigned int i = 0; i < number_aqueous_components_; i++) {
    beaker_state_.total.at(i) = (*aqueous_components)[i][c];
  }

  const Epetra_MultiVector& free_ion = *S_->GetFieldData(free_ion_species_key_)->ViewComponent("cell");
  for (int i = 0; i < number_aqueous_components_; ++i) {
    beaker_state_.free_ion.at(i) = free_ion[i][c];
  }

  // activity coefficients
  if (beaker_state_.primary_activity_coeff.size() > 0) {
    const Epetra_MultiVector& activity = *S_->GetFieldData(primary_activity_coeff_key_)->ViewComponent("cell");
    for (int i = 0; i < beaker_state_.primary_activity_coeff.size(); ++i) {
      beaker_state_.primary_activity_coeff.at(i) = activity[i][c];
    }
  }

  if (beaker_state_.secondary_activity_coeff.size() > 0) {
    const Epetra_MultiVector& activity = *S_->GetFieldData(secondary_activity_coeff_key_)->ViewComponent("cell");
    for (int i = 0; i < beaker_state_.secondary_activity_coeff.size(); ++i) {
      beaker_state_.secondary_activity_coeff.at(i) = activity[i][c];
    }
  }

  // minerals
  if (number_minerals_ > 0) {
    const Epetra_MultiVector& mineral_vf = *S_->GetFieldData(min_vol_frac_key_)->ViewComponent("cell");
    const Epetra_MultiVector& mineral_ssa = *S_->GetFieldData(min_ssa_key_)->ViewComponent("cell");

    for (int i = 0; i < number_minerals_; ++i) {
      beaker_state_.mineral_volume_fraction[i] = mineral_vf[i][c];
      beaker_state_.mineral_specific_surface_area.at(i) = mineral_ssa[i][c];
    }
  }

  // general sorption storage
  if (using_sorption_) {
    const Epetra_MultiVector& sorbed = *S_->GetFieldData(total_sorbed_key_)->ViewComponent("cell");
    for (int i = 0; i < number_aqueous_components_; ++i) {
      beaker_state_.total_sorbed.at(i) = sorbed[i][c];
    }
  }

  // ion exchange
  // TODO: only allow one ion exchange site at the moment!
  if (number_ion_exchange_sites_ > 0) {
    const Epetra_MultiVector& ion_exchange = *S_->GetFieldData(ion_exchange_sites_key_)->ViewComponent("cell");
    for (unsigned int i = 0; i < number_ion_exchange_sites_; i++) {
      beaker_state_.ion_exchange_sites[i] = ion_exchange[i][c];
      // TODO(bandre): need to save ion exchange ref cation conc here!
    }
  }
  
  if (beaker_state_.ion_exchange_ref_cation_conc.size() > 0) {
    const Epetra_MultiVector& ion_exchange = *S_->GetFieldData(ion_exchange_ref_cation_conc_key_)->ViewComponent("cell");
    for (unsigned int i = 0; i < beaker_state_.ion_exchange_ref_cation_conc.size(); ++i) {
      beaker_state_.ion_exchange_ref_cation_conc.at(i) = ion_exchange[i][c];
    }
  }

  // surface complexation
  if (number_sorption_sites_ > 0) {
    const Epetra_MultiVector& sorption_sites = *S_->GetFieldData(sorp_sites_key_)->ViewComponent("cell");

    for (int i = 0; i < number_sorption_sites_; ++i) {
      beaker_state_.surface_site_density[i] = sorption_sites[i][c];
      // TODO(bandre): need to save surface complexation free site conc here!
    }
  }

  if (beaker_state_.surface_complex_free_site_conc.size() > 0) {
    const Epetra_MultiVector& surface_complex =
        *S_->GetFieldData(surf_cfsc_key_)->ViewComponent("cell");
    for (int i = 0; i < beaker_state_.surface_complex_free_site_conc.size(); ++i) {
      beaker_state_.surface_complex_free_site_conc.at(i) = surface_complex[i][c];
    }
  }

  // sorption isotherms
  if (using_sorption_isotherms_) {
    const Epetra_MultiVector& isotherm_kd = *S_->GetFieldData(isotherm_kd_key_)->ViewComponent("cell");
    const Epetra_MultiVector& isotherm_freundlich_n = *S_->GetFieldData(isotherm_freundlich_n_key_)->ViewComponent("cell");
    const Epetra_MultiVector& isotherm_langmuir_b = *S_->GetFieldData(isotherm_langmuir_b_key_)->ViewComponent("cell");

    for (unsigned int i = 0; i < number_aqueous_components_; ++i) {
      beaker_state_.isotherm_kd.at(i) = isotherm_kd[i][c];
      beaker_state_.isotherm_freundlich_n.at(i) = isotherm_freundlich_n[i][c];
      beaker_state_.isotherm_langmuir_b.at(i) = isotherm_langmuir_b[i][c];
    }
  }

  // copy data from state arrays into the beaker parameters
  const Epetra_MultiVector& porosity = *S_->GetFieldData(poro_key_)->ViewComponent("cell", true);
  const Epetra_MultiVector& water_saturation = *S_->GetFieldData(saturation_key_)->ViewComponent("cell", true);
  const Epetra_MultiVector& fluid_density = *S_->GetFieldData(fluid_den_key_)->ViewComponent("cell", true);

  beaker_state_.water_density = fluid_density[0][c];
  beaker_state_.porosity = porosity[0][c];
  beaker_state_.saturation = water_saturation[0][c];
  beaker_state_.volume = mesh_->cell_volume(c);
}


/* *******************************************************************
* Copy data from the beaker back into the state arrays.
******************************************************************* */
void Amanzi_PK::CopyBeakerStructuresToCellState(
    int c, Teuchos::RCP<Epetra_MultiVector> aqueous_components)
{
  for (unsigned int i = 0; i < number_aqueous_components_; ++i) {
    (*aqueous_components)[i][c] = beaker_state_.total.at(i);
  }

  const Epetra_MultiVector& free_ion = *S_->GetFieldData(free_ion_species_key_)->ViewComponent("cell");
  for (int i = 0; i < number_aqueous_components_; ++i) {
    free_ion[i][c] = beaker_state_.free_ion.at(i);
  }

  // activity coefficients
  if (beaker_state_.primary_activity_coeff.size() > 0) {
    const Epetra_MultiVector& activity = *S_->GetFieldData(primary_activity_coeff_key_)->ViewComponent("cell");
    for (int i = 0; i < beaker_state_.primary_activity_coeff.size(); ++i) {
      activity[i][c] = beaker_state_.primary_activity_coeff.at(i);
    }
  }

  if (beaker_state_.secondary_activity_coeff.size() > 0) {
    const Epetra_MultiVector& activity = *S_->GetFieldData(secondary_activity_coeff_key_)->ViewComponent("cell");
    for (int i = 0; i < beaker_state_.secondary_activity_coeff.size(); ++i) {
      activity[i][c] =  beaker_state_.secondary_activity_coeff.at(i);
    }
  }

  // minerals
  if (number_minerals_ > 0) {
    const Epetra_MultiVector& mineral_vf = *S_->GetFieldData(min_vol_frac_key_)->ViewComponent("cell");
    const Epetra_MultiVector& mineral_ssa = *S_->GetFieldData(min_ssa_key_)->ViewComponent("cell");

    for (int i = 0; i < number_minerals_; ++i) {
      mineral_vf[i][c] = beaker_state_.mineral_volume_fraction.at(i);
      mineral_ssa[i][c] = beaker_state_.mineral_specific_surface_area.at(i);
    }
  }

  // sorption
  if (using_sorption_) {
    const Epetra_MultiVector& sorbed = *S_->GetFieldData(total_sorbed_key_)->ViewComponent("cell");
    for (int i = 0; i < number_aqueous_components_; ++i) {
      sorbed[i][c] = beaker_state_.total_sorbed.at(i);
    }
  }

  // surface complexation
  if (number_sorption_sites_ > 0) {
    const Epetra_MultiVector& sorption_sites = *S_->GetFieldData(sorp_sites_key_)->ViewComponent("cell");
    for (int i = 0; i < number_sorption_sites_; i++) {
      sorption_sites[i][c] = beaker_state_.surface_site_density.at(i);
      // TODO: need to save surface complexation free site conc here!
    }
  }

  if (beaker_state_.surface_complex_free_site_conc.size() > 0) {
    const Epetra_MultiVector& surface_complex =
        *S_->GetFieldData(surf_cfsc_key_)->ViewComponent("cell");
    for (int i = 0; i < beaker_state_.surface_complex_free_site_conc.size(); ++i) {
      surface_complex[i][c] = beaker_state_.surface_complex_free_site_conc.at(i);
    }
  }

  // ion exchange
  if (number_ion_exchange_sites_ > 0) {
    const Epetra_MultiVector& ion_exchange = *S_->GetFieldData(ion_exchange_sites_key_)->ViewComponent("cell");
    for (int i = 0; i < number_ion_exchange_sites_; ++i) {
      ion_exchange[i][c] = beaker_state_.ion_exchange_sites.at(i);
      // TODO(bandre): need to save ion exchange ref cation conc here!
    }
  }

  if (beaker_state_.ion_exchange_ref_cation_conc.size() > 0) {
    const Epetra_MultiVector& ion_exchange = *S_->GetFieldData(ion_exchange_ref_cation_conc_key_)->ViewComponent("cell");
    for (int i = 0; i < beaker_state_.ion_exchange_ref_cation_conc.size(); ++i) {
      ion_exchange[i][c] = beaker_state_.ion_exchange_ref_cation_conc.at(i);
    }
  }

  // sorption isotherms
  if (using_sorption_isotherms_) {
    const Epetra_MultiVector& isotherm_kd = *S_->GetFieldData(isotherm_kd_key_)->ViewComponent("cell");
    const Epetra_MultiVector& isotherm_freundlich_n = *S_->GetFieldData(isotherm_freundlich_n_key_)->ViewComponent("cell");
    const Epetra_MultiVector& isotherm_langmuir_b = *S_->GetFieldData(isotherm_langmuir_b_key_)->ViewComponent("cell");

    for (unsigned int i = 0; i < number_aqueous_components_; ++i) {
      isotherm_kd[i][c] = beaker_state_.isotherm_kd.at(i);
      isotherm_freundlich_n[i][c] = beaker_state_.isotherm_freundlich_n.at(i);
      isotherm_langmuir_b[i][c] = beaker_state_.isotherm_langmuir_b.at(i);
    }
  }

  // TODO(bandre): if chemistry can modify the porosity or density,
  // then they should be updated here!
}


/* ******************************************************************
* This function advances concentrations in the auxialiry vector 
* aqueous_components_ (defined in the base class). This vector must 
* be set up using routine set_aqueous_components(). Tipically, it 
* contains values advected by the transport PK.
******************************************************************* */
bool Amanzi_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool failed(false);
  std::string internal_msg;

  double dt = t_new - t_old;
  current_time_ = saved_time_ + dt;

  int num_itrs, max_itrs(0), min_itrs(10000000), avg_itrs(0);
  int cmax(-1), ierr(0);

  // Ensure dependencies are filled
  S_->GetFieldEvaluator(poro_key_)->HasFieldChanged(S_.ptr(), name_);
  S_->GetFieldEvaluator(fluid_den_key_)->HasFieldChanged(S_.ptr(), name_);
  S_->GetFieldEvaluator(saturation_key_)->HasFieldChanged(S_.ptr(), name_);

  const Epetra_MultiVector& sat_vec = *S_->GetFieldData(saturation_key_)->ViewComponent("cell");
  int num_cells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < num_cells; ++c) {
    if (sat_vec[0][c] > saturation_tolerance_) {
      CopyCellStateToBeakerState(c, aqueous_components_);
      try {
        // create a backup copy of the components
        chem_->CopyState(beaker_state_, &beaker_state_copy_);

        // chemistry computations for this cell
        num_itrs = chem_->ReactionStep(&beaker_state_, dt);

        if (max_itrs < num_itrs) {
          max_itrs = num_itrs;
          cmax = c;
        }
        if (min_itrs > num_itrs) {
          min_itrs = num_itrs;
        }
        avg_itrs += num_itrs;
      } catch (Exceptions::Amanzi_exception& geochem_err) {
        ierr = 1;
        internal_msg = geochem_err.what();
        break;
      }

      if (ierr == 0) CopyBeakerStructuresToCellState(c, aqueous_components_);
      // TODO: was porosity etc changed? copy someplace
    }
  }

  ErrorAnalysis(ierr, internal_msg);
  
  int tmp(max_itrs);
  mesh_->get_comm()->MaxAll(&tmp, &max_itrs, 1);

  std::stringstream ss;
  ss << "Newton iterations: " << min_itrs << "/" << max_itrs << "/" 
     << avg_itrs / num_cells << ", maximum in gid=" << mesh_->GID(cmax, AmanziMesh::CELL) << std::endl;
  vo_->Write(Teuchos::VERB_HIGH, ss.str());

  // dumping the values of the final cell. not very helpful by itself,
  // but can be move up into the loops....
  // chem_->DisplayTotalColumnHeaders(true);
  // chem_->DisplayTotalColumns(current_time_, beaker_state_, true);

  // update time control parameters
  num_successful_steps_++;
  num_iterations_ = max_itrs;

  return failed;
}


/* ******************************************************************
* The MPC will call this function to signal to the process kernel 
* that it has accepted the state update, thus, the PK should update
* possible auxilary state variables here
******************************************************************* */
void Amanzi_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  saved_time_ = t_new;

  if (dt_control_method_ == "simple") {
    if ((num_successful_steps_ == 0) || (num_iterations_ >= dt_cut_threshold_)) {
      dt_next_ /= dt_cut_factor_;
    }
    else if (num_successful_steps_ >= dt_increase_threshold_) {
      dt_next_ *= dt_increase_factor_;
      num_successful_steps_ = 0;
    }

    dt_next_ = std::min(dt_next_, dt_max_);

    // synchronize processors since update of control parameters was 
    // made in many places.
    double tmp(dt_next_);
    mesh_->get_comm()->MinAll(&tmp, &dt_next_, 1);
  }
}


/* ******************************************************************
*
******************************************************************* */
Teuchos::RCP<Epetra_MultiVector> Amanzi_PK::extra_chemistry_output_data() {
  if (aux_data_ != Teuchos::null) {
    const Epetra_MultiVector& free_ion = *S_->GetFieldData(free_ion_species_key_)->ViewComponent("cell");
    const Epetra_MultiVector& activity = *S_->GetFieldData(primary_activity_coeff_key_)->ViewComponent("cell");
    int num_cells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

    for (int cell = 0; cell < num_cells; cell++) {
      // populate aux_data_ by copying from the appropriate internal storage
      for (unsigned int i = 0; i < aux_names_.size(); i++) {
        if (aux_names_.at(i) == "pH") {
          double* cell_aux_data = (*aux_data_)[i];
          double* cell_free_ion = free_ion[aux_index_.at(i)];
          double* activity_coeff = activity[aux_index_.at(i)];
          double activity_val = cell_free_ion[cell] * activity_coeff[cell];
          cell_aux_data[cell] = -std::log10(activity_val);
        } else {
          // don't support anything else at this time....
        }
      }
    }
  }
  return aux_data_;
}


/* ******************************************************************
*
******************************************************************* */
void Amanzi_PK::set_chemistry_output_names(std::vector<std::string>* names) {
  names->clear();
  
  for (std::vector<std::string>::const_iterator name = aux_names_.begin();
       name != aux_names_.end(); name++) {
    names->push_back(*name);
  }
}

}  // namespace AmanziChemistry
}  // namespace Amanzi

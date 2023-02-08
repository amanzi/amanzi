/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry PK

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
#include "Teuchos_XMLParameterListHelpers.hpp"

// Amanzi
#include "Mesh.hh"
#include "Beaker.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "SimpleThermoDatabase.hh"
#include "VerboseObject.hh"

// Chemistry
#include "Amanzi_PK.hh"
#include "ChemistryDefs.hh"

namespace Amanzi {
namespace AmanziChemistry {

/* ******************************************************************
* Constructor
******************************************************************* */
Amanzi_PK::Amanzi_PK(Teuchos::ParameterList& pk_tree,
                     const Teuchos::RCP<Teuchos::ParameterList>& glist,
                     const Teuchos::RCP<State>& S,
                     const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, glist, S, soln),
    Chemistry_PK(pk_tree, glist, S, soln),
    chem_(NULL),
    current_time_(0.0),
    saved_time_(0.0)
{
  // obtain key of fields
  tcc_key_ = Keys::getKey(domain_, "total_component_concentration");
  poro_key_ = plist_->get<std::string>("porosity key", Keys::getKey(domain_, "porosity"));
  saturation_key_ =
    plist_->get<std::string>("saturation key", Keys::getKey(domain_, "saturation_liquid"));
  fluid_den_key_ =
    plist_->get<std::string>("fluid density key", Keys::getKey(domain_, "mass_density_liquid"));

  temperature_key_ = Keys::getKey(domain_, "temperature");
  min_vol_frac_key_ = Keys::getKey(domain_, "mineral_volume_fractions");
  min_ssa_key_ = Keys::getKey(domain_, "mineral_specific_surface_area");
  sorp_sites_key_ = Keys::getKey(domain_, "sorption_sites");
  surf_cfsc_key_ = Keys::getKey(domain_, "surface_complex_free_site_conc");
  total_sorbed_key_ = Keys::getKey(domain_, "total_sorbed");
  isotherm_kd_key_ = Keys::getKey(domain_, "isotherm_kd");
  isotherm_freundlich_n_key_ = Keys::getKey(domain_, "isotherm_freundlich_n");
  isotherm_langmuir_b_key_ = Keys::getKey(domain_, "isotherm_langmuir_b");
  free_ion_species_key_ = Keys::getKey(domain_, "free_ion_species");
  primary_activity_coeff_key_ = Keys::getKey(domain_, "primary_activity_coeff");

  ion_exchange_sites_key_ = Keys::getKey(domain_, "ion_exchange_sites");
  ion_exchange_ref_cation_conc_key_ = Keys::getKey(domain_, "ion_exchange_ref_cation_conc");
  secondary_activity_coeff_key_ = Keys::getKey(domain_, "secondary_activity_coeff");
  alquimia_aux_data_key_ = Keys::getKey(domain_, "alquimia_aux_data");
  mineral_rate_constant_key_ = Keys::getKey(domain_, "mineral_rate_constant");
  first_order_decay_constant_key_ = Keys::getKey(domain_, "first_order_decay_constant");

  // collect high-level information about the problem
  InitializeMinerals(plist_);
  InitializeSorptionSites(plist_, S_->ICList());

  // grab the component names
  comp_names_.clear();
  Teuchos::RCP<Teuchos::ParameterList> cd_list = Teuchos::sublist(glist, "cycle driver", true);
  if (cd_list->isParameter("component names")) {
    comp_names_ = cd_list->get<Teuchos::Array<std::string>>("component names").toVector();
  } else {
    Errors::Message msg("Amanzi_PK: Cycle Driver has no input parameter component names.");
    Exceptions::amanzi_throw(msg);
  }
  number_aqueous_components_ = comp_names_.size();
  number_free_ion_ = number_aqueous_components_;
  number_total_sorbed_ = number_aqueous_components_;

  // This intentionally overrides the PK construction of vo_ to set the name to
  // what Konstantin wants it to be for Alquimia_PK.  This needs more
  // discussion -- see #672.
  //
  // overriding the vo plist for individual PKs in a collection of PKs
  Teuchos::RCP<Teuchos::ParameterList> vo_plist = plist_;
  if (plist_->isSublist(name_ + " verbose object")) {
    vo_plist = Teuchos::rcp(new Teuchos::ParameterList(*plist_));
    vo_plist->set("verbose object", plist_->sublist(name_ + " verbose object"));
  }

  //  some tests provide nullptr
  name_ = "Amanzi_PK:" + domain_;
  if (solution_.get())
    vo_ = Teuchos::rcp(new VerboseObject(solution_->Comm(), name_, *vo_plist));
  else
    vo_ = Teuchos::rcp(new VerboseObject(getDefaultComm(), name_, *vo_plist));
}


/* ******************************************************************
* Register fields and evaluators with the State
******************************************************************* */
void
Amanzi_PK::Setup()
{
  Chemistry_PK::Setup();
}


/* ******************************************************************
* Can this be done during Setup phase?
******************************************************************* */
void
Amanzi_PK::AllocateAdditionalChemistryStorage_()
{
  int n_secondary_comps = chem_->secondary_species().size();
  if (n_secondary_comps > 0) {
    S_->Require<CompositeVector, CompositeVectorSpace>(
        secondary_activity_coeff_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, n_secondary_comps);

    S_->GetRecordSetW(secondary_activity_coeff_key_).CreateData();

    S_->GetW<CompositeVector>(secondary_activity_coeff_key_, passwd_).PutScalar(1.0);
    S_->GetRecordW(secondary_activity_coeff_key_, Tags::DEFAULT, passwd_).set_initialized();
  }
}


/* *******************************************************************
* Initialization
******************************************************************* */
void
Amanzi_PK::Initialize()
{
  ncells_owned_ = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  // initialization using base class
  Chemistry_PK::Initialize();

  auto tcc =
    S_->GetPtrW<CompositeVector>(tcc_key_, Tags::DEFAULT, passwd_)->ViewComponent("cell", true);

  XMLParameters();

  // TODO: some sort of check of the state object to see if mineral_ssa,
  // CEC, site density, etc is present.

  // initial conditions for minerals etc should be handled by the
  // state/chemistry_state object before we reach this point. We just
  // resize our local memory for migrating data here.

  chem_->Initialize(beaker_state_, beaker_parameters_);

  // allocating memory for auxiliary fields
  AllocateAdditionalChemistryStorage_();
  SetupAuxiliaryOutput();

  // populate list of beaker fields
  InitializeBeakerFields_();

  if (ncells_owned_ > 0) {
    CopyCellStateToBeakerState(0, tcc);
    chem_->CopyStateToBeaker(beaker_state_);
    // chem_->VerifyState(beaker_state_);
  }

  // check names of primary species
  int nprimary = chem_->primary_species().size();
  if (nprimary == comp_names_.size()) {
    for (int i = 0; i < nprimary; ++i) {
      std::string species_name = chem_->primary_species().at(i).name();
      if (comp_names_[i] != species_name) {
        Errors::Message msg;
        msg << "Amanzi PK: mismatch of name: \"" << comp_names_[i] << "\" and \"" << species_name
            << "\". Compare XML and BGD lists.";
        Exceptions::amanzi_throw(msg);
      }
    }
  }

  // copy the cell data into the beaker storage for initialization purposes
  // but first ensure dependencies are filled
  S_->GetEvaluator(poro_key_).Update(*S_, name_);
  S_->GetEvaluator(fluid_den_key_).Update(*S_, name_);
  S_->GetEvaluator(saturation_key_).Update(*S_, name_);

  // finish setting up & testing the chemistry object
  int ierr(0);
  std::string internal_msg;
  std::vector<double> values(nprimary);

  const auto& iclist = glist_->sublist("state").sublist("initial conditions").sublist(tcc_key_);
  std::vector<std::string> constraints;
  if (iclist.isParameter("names")) {
    constraints = iclist.get<Teuchos::Array<std::string>>("names").toVector();
  } else {
    for (int i = 0; i < nprimary; ++i) constraints.push_back("total");
  }

  // print statistics
  try {
    vo_->Write(Teuchos::VERB_HIGH, "Initializing chemistry in cell 0...\n");
    chem_->Display();
    vo_->Write(Teuchos::VERB_HIGH, "Initial speciation calculations in cell 0...\n");

    if (ncells_owned_ > 0) {
      for (int i = 0; i < nprimary; ++i) values[i] = (*tcc)[i][0];
      chem_->EnforceConstraint(&beaker_state_, beaker_parameters_, constraints, values);

      vo_->Write(Teuchos::VERB_HIGH, "\nTest solution of initial conditions in cell 0:\n");
      chem_->DisplayResults();
    }
  } catch (Exceptions::Amanzi_exception& geochem_err) {
    ierr = 1;
    internal_msg = geochem_err.what();
  }

  ErrorAnalysis(ierr, internal_msg);

  // compute the equilibrium state
  int num_cells = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  ierr = 0;
  if (fabs(initial_conditions_time_ - S_->get_time()) <= 1e-8 * fabs(S_->get_time())) {
    for (int c = 0; c < num_cells; ++c) {
      CopyCellStateToBeakerState(c, tcc);

      try {
        for (int i = 0; i < nprimary; ++i) values[i] = (*tcc)[i][c];
        chem_->EnforceConstraint(&beaker_state_, beaker_parameters_, constraints, values);
        CopyBeakerStructuresToCellState(c, tcc);
      } catch (Exceptions::Amanzi_exception& geochem_err) {
        ierr = 1;
        internal_msg = geochem_err.what();
      }
    }
  } else {
    if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "no data initialization due to time mismatch: " << S_->get_time() << std::endl;
    }
  }

  // figure out if any of the processes threw an error, if so all processes will re-throw
  ErrorAnalysis(ierr, internal_msg);

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << vo_->color("green") << "Initialization of PK was successful, T=" << S_->get_time()
               << vo_->reset() << std::endl
               << std::endl;
  }
}


/* *******************************************************************
* Initialization helper functions
******************************************************************* */
void
Amanzi_PK::XMLParameters()
{
  // search for thermodynamic database
  Teuchos::RCP<Teuchos::ParameterList> tdb_list;
  if (plist_->isSublist("thermodynamic database")) {
    // -- search by the name
    if (plist_->sublist("thermodynamic database").isParameter("file")) {
      // check extension and read
      std::string filename = plist_->sublist("thermodynamic database").get<std::string>("file");
      tdb_list = Teuchos::getParametersFromXmlFile(filename);
    }
  }

  // search as a sublist of the global list
  if (tdb_list == Teuchos::null) {
    tdb_list = Teuchos::sublist(glist_, "thermodynamic database", true);
  }

  // currently we only support the simple format.
  chem_ = std::make_shared<SimpleThermoDatabase>(tdb_list, vo_);
  beaker_parameters_.tolerance = 1e-12;
  beaker_parameters_.max_iterations = 250;
  beaker_parameters_.activity_model_name = "unit";

  if (tdb_list == Teuchos::null) {
    std::ostringstream msg;
    msg << "Amanzi_PK: 'thermodynamic database' sublist must be specified.\n";
    Exceptions::amanzi_throw(Errors::Message(msg.str()));
  }

  // activity model
  beaker_parameters_.activity_model_name = plist_->get<std::string>("activity model", "unit");
  // -- Pitzer virial coefficients database
  if (beaker_parameters_.activity_model_name == "pitzer-hwm") {
    if (plist_->isParameter("Pitzer database file")) {
      beaker_parameters_.pitzer_database = plist_->get<std::string>("Pitzer database file");
    } else {
      std::ostringstream msg;
      msg << "Amanzi_PK: parameter 'Pitzer database file' must be specified if activity "
             "model=pitzer-hwm'.\n";
      Exceptions::amanzi_throw(Errors::Message(msg.str()));
    }
  }

  // solver parameters
  beaker_parameters_.tolerance = plist_->get<double>("tolerance", 1.0e-12);
  beaker_parameters_.max_iterations = plist_->get<int>("maximum Newton iterations", 200);

  // auxiliary data
  aux_names_.clear();
  if (plist_->isParameter("auxiliary data")) {
    Teuchos::Array<std::string> names = plist_->get<Teuchos::Array<std::string>>("auxiliary data");
    for (auto name = names.begin(); name != names.end(); ++name) {
      if (*name == "pH") {
        aux_names_.push_back(*name);
      } else {
        std::stringstream message;
        message << "XMLParameters(): unknown value in 'auxiliary data' list: " << *name
                << std::endl;
        vo_->WriteWarning(Teuchos::VERB_LOW, message);
      }
    }
  }

  // misc other chemistry flags
  dt_control_method_ = plist_->get<std::string>("time step control method", "fixed");
  dt_max_ = plist_->get<double>("max time step (s)", 9.9e+9);
  dt_next_ = plist_->get<double>("initial time step (s)", dt_max_);
  dt_cut_factor_ = plist_->get<double>("time step cut factor", 2.0);
  dt_increase_factor_ = plist_->get<double>("time step increase factor", 1.2);

  dt_cut_threshold_ = plist_->get<int>("time step cut threshold", 8);
  dt_increase_threshold_ = plist_->get<int>("time step increase threshold", 4);

  num_successful_steps_ = 0;
}


/* *******************************************************************
* Requires that Beaker::Setup() has already been called!
******************************************************************* */
void
Amanzi_PK::SetupAuxiliaryOutput()
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
      message << "ChemistryPK: Output was requested for '" << aux_names_.at(i) << "' (" << name
              << ") but no chemistry varibles of this name were found.\n";
      Exceptions::amanzi_throw(Errors::Message(message.str()));
    } else {
      aux_index_.push_back(index);
    }
  }

  // create the Epetra_MultiVector that will hold the data
  if (nvars > 0) {
    aux_data_ = Teuchos::rcp(new Epetra_MultiVector(mesh_->getMap(AmanziMesh::Entity_kind::CELL,false), nvars));
  } else {
    aux_data_ = Teuchos::null;
  }
}


/* *******************************************************************
* We must use the aqueous totals value calculated from transport
* (aqueous_components), not the value stored in state!
******************************************************************* */
void
Amanzi_PK::CopyCellStateToBeakerState(int c, Teuchos::RCP<Epetra_MultiVector> aqueous_components)
{
  for (unsigned int i = 0; i < number_aqueous_components_; i++) {
    beaker_state_.total.at(i) = (*aqueous_components)[i][c];
  }

  for (int i = 0; i < number_aqueous_components_; ++i) {
    beaker_state_.free_ion.at(i) = (*bf_.free_ion)[i][c];
  }

  // activity coefficients
  for (int i = 0; i < beaker_state_.primary_activity_coeff.size(); ++i) {
    beaker_state_.primary_activity_coeff.at(i) = (*bf_.activity)[i][c];
  }

  if (beaker_state_.secondary_activity_coeff.size() > 0) {
    for (int i = 0; i < beaker_state_.secondary_activity_coeff.size(); ++i) {
      beaker_state_.secondary_activity_coeff.at(i) = (*bf_.secondary_activity)[i][c];
    }
  }

  // minerals
  if (number_minerals_ > 0) {
    for (int i = 0; i < number_minerals_; ++i) {
      beaker_state_.mineral_volume_fraction[i] = (*bf_.mineral_vf)[i][c];
      beaker_state_.mineral_specific_surface_area.at(i) = (*bf_.mineral_ssa)[i][c];
    }
  }

  // general sorption storage
  if (using_sorption_) {
    for (int i = 0; i < number_aqueous_components_; ++i) {
      beaker_state_.total_sorbed.at(i) = (*bf_.sorbed)[i][c];
    }
  }

  // ion exchange
  // TODO: only allow one ion exchange site at the moment!
  if (number_ion_exchange_sites_ > 0) {
    for (unsigned int i = 0; i < number_ion_exchange_sites_; i++) {
      beaker_state_.ion_exchange_sites[i] = (*bf_.ion_exchange_sites)[i][c];
      // TODO(bandre): need to save ion exchange ref cation conc here!
      beaker_state_.ion_exchange_ref_cation_conc.at(i) = (*bf_.ion_exchange_ref_cation_conc)[i][c];
    }
  }

  // surface complexation
  if (number_sorption_sites_ > 0) {
    for (int i = 0; i < number_sorption_sites_; ++i) {
      beaker_state_.surface_site_density[i] = (*bf_.sorption_sites)[i][c];
      // TODO(bandre): need to save surface complexation free site conc here!
    }
  }

  if (beaker_state_.surface_complex_free_site_conc.size() > 0) {
    for (int i = 0; i < beaker_state_.surface_complex_free_site_conc.size(); ++i) {
      beaker_state_.surface_complex_free_site_conc.at(i) = (*bf_.surface_complex)[i][c];
    }
  }

  // sorption isotherms provided as material-based property must be copied to
  if (using_sorption_isotherms_) {
    for (unsigned int i = 0; i < number_aqueous_components_; ++i) {
      beaker_state_.isotherm_kd.at(i) = (*bf_.isotherm_kd)[i][c];
      beaker_state_.isotherm_freundlich_n.at(i) = (*bf_.isotherm_freundlich_n)[i][c];
      beaker_state_.isotherm_langmuir_b.at(i) = (*bf_.isotherm_langmuir_b)[i][c];
    }
  }

  // copy data from state arrays into the beaker parameters
  beaker_state_.water_density = (*bf_.density)[0][c];
  beaker_state_.porosity = (*bf_.porosity)[0][c];
  beaker_state_.saturation = (*bf_.saturation)[0][c];
  beaker_state_.volume = mesh_->getCellVolume(c);

  if (S_->HasRecord(temperature_key_)) { beaker_state_.temperature = (*bf_.temperature)[0][c]; }
}


/* *******************************************************************
* Copy data from the beaker back into the state arrays.
******************************************************************* */
void
Amanzi_PK::CopyBeakerStructuresToCellState(int c,
                                           Teuchos::RCP<Epetra_MultiVector> aqueous_components)
{
  for (unsigned int i = 0; i < number_aqueous_components_; ++i) {
    (*aqueous_components)[i][c] = std::max(beaker_state_.total.at(i), TCC_MIN_VALUE);
  }

  for (int i = 0; i < number_aqueous_components_; ++i) {
    (*bf_.free_ion)[i][c] = std::max(beaker_state_.free_ion.at(i), TCC_MIN_VALUE);
  }

  // activity coefficients
  for (int i = 0; i < beaker_state_.primary_activity_coeff.size(); ++i) {
    (*bf_.activity)[i][c] = beaker_state_.primary_activity_coeff.at(i);
  }

  if (beaker_state_.secondary_activity_coeff.size() > 0) {
    for (int i = 0; i < beaker_state_.secondary_activity_coeff.size(); ++i) {
      (*bf_.secondary_activity)[i][c] = beaker_state_.secondary_activity_coeff.at(i);
    }
  }

  // minerals
  if (number_minerals_ > 0) {
    for (int i = 0; i < number_minerals_; ++i) {
      (*bf_.mineral_vf)[i][c] = beaker_state_.mineral_volume_fraction.at(i);
      (*bf_.mineral_ssa)[i][c] = beaker_state_.mineral_specific_surface_area.at(i);
    }
  }

  // sorption
  if (using_sorption_) {
    for (int i = 0; i < number_aqueous_components_; ++i) {
      (*bf_.sorbed)[i][c] = beaker_state_.total_sorbed.at(i);
    }
  }

  // surface complexation
  if (number_sorption_sites_ > 0) {
    for (int i = 0; i < number_sorption_sites_; i++) {
      (*bf_.sorption_sites)[i][c] = beaker_state_.surface_site_density.at(i);
      // TODO: need to save surface complexation free site conc here!
    }
  }

  if (beaker_state_.surface_complex_free_site_conc.size() > 0) {
    for (int i = 0; i < beaker_state_.surface_complex_free_site_conc.size(); ++i) {
      (*bf_.surface_complex)[i][c] = beaker_state_.surface_complex_free_site_conc.at(i);
    }
  }

  // ion exchange
  if (number_ion_exchange_sites_ > 0) {
    for (int i = 0; i < number_ion_exchange_sites_; ++i) {
      (*bf_.ion_exchange_sites)[i][c] = beaker_state_.ion_exchange_sites.at(i);
      (*bf_.ion_exchange_ref_cation_conc)[i][c] = beaker_state_.ion_exchange_ref_cation_conc.at(i);
      // TODO(bandre): need to save ion exchange ref cation conc here!
    }
  }

  // sorption isotherms
  if (using_sorption_isotherms_) {
    for (unsigned int i = 0; i < number_aqueous_components_; ++i) {
      (*bf_.isotherm_kd)[i][c] = beaker_state_.isotherm_kd.at(i);
      (*bf_.isotherm_freundlich_n)[i][c] = beaker_state_.isotherm_freundlich_n.at(i);
      (*bf_.isotherm_langmuir_b)[i][c] = beaker_state_.isotherm_langmuir_b.at(i);
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
bool
Amanzi_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool failed(false);
  std::string internal_msg;

  dt_ = t_new - t_old;
  current_time_ = saved_time_ + dt_;

  int num_itrs, max_itrs(0), min_itrs(10000000), avg_itrs(0);
  int cmax(-1), ierr(0);

  // Ensure dependencies are filled
  S_->GetEvaluator(poro_key_).Update(*S_, name_);
  S_->GetEvaluator(fluid_den_key_).Update(*S_, name_);
  S_->GetEvaluator(saturation_key_).Update(*S_, name_);

  const auto& sat_vec = *S_->Get<CompositeVector>(saturation_key_).ViewComponent("cell");

  for (int c = 0; c < ncells_owned_; ++c) {
    if (sat_vec[0][c] > saturation_tolerance_) {
      CopyCellStateToBeakerState(c, aqueous_components_);
      try {
        // create a backup copy of the components
        chem_->CopyState(beaker_state_, &beaker_state_copy_);

        // chemistry computations for this cell
        num_itrs = chem_->ReactionStep(&beaker_state_, beaker_parameters_, dt_);

        if (max_itrs < num_itrs) {
          max_itrs = num_itrs;
          cmax = c;
        }
        if (min_itrs > num_itrs) { min_itrs = num_itrs; }
        avg_itrs += num_itrs;
      } catch (Exceptions::Amanzi_exception& geochem_err) {
        ierr = 1;
        internal_msg = geochem_err.what();
        // chem_->DisplayComponents(beaker_state_);
        break;
      } catch (...) {
        ierr = 1;
        internal_msg = "unknown error";
        break;
      }

      if (ierr == 0) CopyBeakerStructuresToCellState(c, aqueous_components_);
      // TODO: was porosity etc changed? copy some place
    }
  }

  ErrorAnalysis(ierr, internal_msg);

  int tmp0(min_itrs), tmp1(max_itrs), tmp2(avg_itrs), tmp3(ncells_owned_), ncells_total;
  mesh_->get_comm()->MinAll(&tmp0, &min_itrs, 1);
  mesh_->get_comm()->MaxAll(&tmp1, &max_itrs, 1);
  mesh_->get_comm()->SumAll(&tmp2, &avg_itrs, 1);
  mesh_->get_comm()->SumAll(&tmp3, &ncells_total, 1);

  std::stringstream ss;
  ss << "Newton iterations: " << min_itrs << "/" << max_itrs << "/"
     << avg_itrs / std::max(ncells_owned_, 1)
     << ", maximum in gid=" << mesh_->GID(cmax, AmanziMesh::CELL) << std::endl;
  vo_->Write(Teuchos::VERB_HIGH, ss.str());

  // update time control parameters
  num_successful_steps_++;
  num_iterations_ = max_itrs;

  if (!failed) EstimateNextTimeStep_(t_old, t_new);
  return failed;
}


/* ******************************************************************
* The MPC will call this function to signal to the process kernel
* that it has accepted the state update, thus, the PK should update
* possible auxilary state variables here
******************************************************************* */
void
Amanzi_PK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  saved_time_ = t_new;
  EstimateNextTimeStep_(t_old, t_new);
}


/* ******************************************************************
* Extract poiters to state fields and place them in a simple struct.
******************************************************************* */
void
Amanzi_PK::InitializeBeakerFields_()
{
  Tag tag(Tags::DEFAULT);

  // constant (tmeporarily) fields
  bf_.porosity = S_->Get<CompositeVector>(poro_key_).ViewComponent("cell");
  bf_.density = S_->Get<CompositeVector>(fluid_den_key_).ViewComponent("cell");
  bf_.saturation = S_->Get<CompositeVector>(saturation_key_).ViewComponent("cell");
  if (S_->HasRecord(temperature_key_))
    bf_.temperature = S_->Get<CompositeVector>(temperature_key_).ViewComponent("cell");

  bf_.free_ion =
    S_->GetW<CompositeVector>(free_ion_species_key_, tag, passwd_).ViewComponent("cell");
  bf_.activity =
    S_->GetW<CompositeVector>(primary_activity_coeff_key_, tag, passwd_).ViewComponent("cell");

  if (chem_->secondary_species().size() > 0) {
    bf_.secondary_activity =
      S_->GetW<CompositeVector>(secondary_activity_coeff_key_, tag, passwd_).ViewComponent("cell");
  }

  if (number_ion_exchange_sites_ > 0) {
    bf_.ion_exchange_sites =
      S_->GetW<CompositeVector>(ion_exchange_sites_key_, tag, passwd_).ViewComponent("cell");
    bf_.ion_exchange_ref_cation_conc =
      S_->GetW<CompositeVector>(ion_exchange_ref_cation_conc_key_, Tags::DEFAULT, passwd_)
        .ViewComponent("cell");
  }

  if (number_minerals_ > 0) {
    bf_.mineral_vf =
      S_->GetW<CompositeVector>(min_vol_frac_key_, tag, passwd_).ViewComponent("cell");
    bf_.mineral_ssa = S_->GetW<CompositeVector>(min_ssa_key_, tag, passwd_).ViewComponent("cell");
  }

  if (using_sorption_) {
    bf_.sorbed = S_->GetW<CompositeVector>(total_sorbed_key_, tag, passwd_).ViewComponent("cell");
  }

  if (using_sorption_isotherms_) {
    bf_.isotherm_kd =
      S_->GetW<CompositeVector>(isotherm_kd_key_, tag, passwd_).ViewComponent("cell");
    bf_.isotherm_freundlich_n =
      S_->GetW<CompositeVector>(isotherm_freundlich_n_key_, tag, passwd_).ViewComponent("cell");
    bf_.isotherm_langmuir_b =
      S_->GetW<CompositeVector>(isotherm_langmuir_b_key_, tag, passwd_).ViewComponent("cell");
  }

  if (number_sorption_sites_ > 0) {
    bf_.sorption_sites =
      S_->GetW<CompositeVector>(sorp_sites_key_, tag, passwd_).ViewComponent("cell");
  }

  if (beaker_state_.surface_complex_free_site_conc.size() > 0) {
    bf_.surface_complex =
      S_->GetW<CompositeVector>(surf_cfsc_key_, tag, passwd_).ViewComponent("cell");
  }
}


/* ******************************************************************
* Estimate next time step, including for subcycling
******************************************************************* */
void
Amanzi_PK::EstimateNextTimeStep_(double t_old, double t_new)
{
  // we keep own estimate of stable time step and report it to CD via get_dt()
  if (dt_control_method_ == "simple") {
    double dt_next_tmp(dt_next_);

    if ((num_successful_steps_ == 0) || (num_iterations_ >= dt_cut_threshold_)) {
      dt_next_ /= dt_cut_factor_;
    } else if (num_successful_steps_ >= dt_increase_threshold_) {
      dt_next_ *= dt_increase_factor_;
      num_successful_steps_ = 0;
    }

    dt_next_ = std::min(dt_next_, dt_max_);

    // synchronize processors since update of control parameters was
    // made in many places.
    double tmp(dt_next_);
    mesh_->get_comm()->MinAll(&tmp, &dt_next_, 1);

    if (dt_next_tmp != dt_next_ && vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "controller changed dt_next from " << dt_next_tmp << " to " << dt_next_
                 << std::endl;
    }
  } else if (dt_control_method_ == "fixed") {
    // dt_next could be reduced by the base PK. To avoid small time step, we reset it here.
    dt_next_ = std::max(dt_, t_new - t_old); // due to subcycling
    dt_next_ = std::min(dt_next_, dt_max_);
  }
}


/* ******************************************************************
*
******************************************************************* */
Teuchos::RCP<Epetra_MultiVector>
Amanzi_PK::extra_chemistry_output_data()
{
  if (aux_data_ != Teuchos::null) {
    const auto& free_ion = *S_->Get<CompositeVector>(free_ion_species_key_).ViewComponent("cell");
    const auto& activity =
      *S_->Get<CompositeVector>(primary_activity_coeff_key_).ViewComponent("cell");

    for (int c = 0; c < ncells_owned_; ++c) {
      // populate aux_data_ by copying from the appropriate internal storage
      for (unsigned int i = 0; i < aux_names_.size(); i++) {
        if (aux_names_.at(i) == "pH") {
          double activity_val = free_ion[aux_index_.at(i)][c] * activity[aux_index_.at(i)][c];
          (*aux_data_)[i][c] = -std::log10(activity_val);
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
void
Amanzi_PK::set_chemistry_output_names(std::vector<std::string>* names)
{
  names->clear();

  for (auto name = aux_names_.begin(); name != aux_names_.end(); name++) {
    names->push_back(*name);
  }
}

} // namespace AmanziChemistry
} // namespace Amanzi

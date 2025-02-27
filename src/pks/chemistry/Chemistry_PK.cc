/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry PK

  Base class for chemical process kernels.
*/

// Chemistry
#include "message.hh"
#include "pk_helpers.hh"
#include "TimestepControllerFactory.hh"
#include "Chemistry_PK.hh"

namespace Amanzi {
namespace AmanziChemistry {

/* ******************************************************************
* Default constructor that initializes all pointers to NULL
****************************************************************** */
Chemistry_PK::Chemistry_PK()
  : number_aqueous_components_(0),
    number_gaseous_components_(0),
    number_mineral_components_(0),
    number_aqueous_kinetics_(0),
    number_sorption_sites_(0),
    using_sorption_(false),
    using_sorption_isotherms_(false),
    number_ion_exchange_sites_(0),
    dt_next_(-1.) {};


Chemistry_PK::Chemistry_PK(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& glist,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& soln)
  : PK_Physical(pk_tree, glist, S, soln),
    PK(pk_tree, glist, S, soln),
    number_aqueous_components_(0),
    number_gaseous_components_(0),
    number_mineral_components_(0),
    number_aqueous_kinetics_(0),
    number_sorption_sites_(0),
    using_sorption_(false),
    using_sorption_isotherms_(false),
    number_ion_exchange_sites_(0),
    dt_next_(-1.)
{
  // note, we pass in null to the factory here to make sure there is no error
  // control used, which doesn't make sense for this application.
  timestep_controller_ = createTimestepController<TreeVector>(name_,
          plist_->sublist("timestep controller"), S_, Teuchos::null, Teuchos::null);
};


double
Chemistry_PK::get_dt()
{
  if (dt_next_ < 0.) dt_next_ = timestep_controller_->getInitialTimestep();
  return dt_next_;
}



/* ******************************************************************
* Parser
******************************************************************* */
void
Chemistry_PK::parseParameterList()
{
  PK_Physical::parseParameterList();

  // get all the keys used here
  // -- external fields
  poro_key_ = Keys::readKey(*plist_, domain_, "porosity", "porosity");
  temperature_key_ = Keys::readKey(*plist_, domain_, "temperature", "temperature");
  saturation_key_ = Keys::readKey(*plist_, domain_, "saturation liquid", "saturation_liquid");
  fluid_den_key_ = Keys::readKey(*plist_, domain_, "mass density liquid", "mass_density_liquid");

  // -- my fields
  mineral_vol_frac_key_ =
    Keys::readKey(*plist_, domain_, "mineral volume fractions", "mineral_volume_fractions");
  mineral_ssa_key_ = Keys::readKey(
    *plist_, domain_, "mineral specific surface area", "mineral_specific_surface_area");
  mineral_rate_constant_key_ =
    Keys::readKey(*plist_, domain_, "mineral rate constant", "mineral_rate_constant");

  total_sorbed_key_ = Keys::readKey(*plist_, domain_, "total sorbed", "total_sorbed");

  isotherm_kd_key_ = Keys::readKey(*plist_, domain_, "isotherm_kd", "isotherm_kd");
  isotherm_freundlich_n_key_ =
    Keys::readKey(*plist_, domain_, "isotherm freundlich_n", "isotherm_freundlich_n");
  isotherm_langmuir_b_key_ =
    Keys::readKey(*plist_, domain_, "isotherm langmuir_b", "isotherm_langmuir_b");

  free_ion_species_key_ = Keys::readKey(*plist_, domain_, "free ion species", "free_ion_species");
  primary_activity_coeff_key_ =
    Keys::readKey(*plist_, domain_, "primary activity coeff", "primary_activity_coeff");

  ion_exchange_sites_key_ =
    Keys::readKey(*plist_, domain_, "ion exchange sites", "ion_exchange_sites");
  ion_exchange_ref_cation_conc_key_ =
    Keys::readKey(*plist_, domain_, "ion exchange ref cation conc", "ion_exchange_ref_cation_conc");
  secondary_activity_coeff_key_ =
    Keys::readKey(*plist_, domain_, "secondary activity coeff", "secondary_activity_coeff");

  alquimia_aux_data_key_ =
    Keys::readKey(*plist_, domain_, "alquimia aux data", "alquimia_aux_data");

  ion_exchange_ref_cation_conc_key_ =
    Keys::readKey(*plist_, domain_, "ion exchange ref cation conc", "ion_exchange_ref_cation_conc");
  secondary_activity_coeff_key_ =
    Keys::readKey(*plist_, domain_, "secondary activity coeff", "secondary_activity_coeff");
  first_order_decay_constant_key_ =
    Keys::readKey(*plist_, domain_, "first order decay constant", "first_order_decay_constant");

  // other parameters
  passwd_ = plist_->get<std::string>("primary variable password", "state");
  saturation_tolerance_ = plist_->get<double>("saturation tolerance", 1e-14);

  // physical models
  parseMinerals_(*plist_, *plist_);
  parseSorptionSites_(*plist_, *plist_);
}


/* ******************************************************************
* Register fields and evaluators with the State
******************************************************************* */
void
Chemistry_PK::Setup()
{
  // require flow fields
  {
    auto keys = std::vector<Key>{ poro_key_, saturation_key_, fluid_den_key_ };
    for (const auto& key : keys) {
      requireAtCurrent(key, tag_current_, *S_)
        .SetMesh(mesh_)
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    }
  }

  // require my keys
  if (number_aqueous_components_ > 0) {
    auto keys = std::vector<Key>{ key_, free_ion_species_key_ };
    for (const auto& key : keys) {
      requireAtNext(key, tag_next_, *S_, true, passwd_)
        .SetMesh(mesh_)
        ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_aqueous_components_);
      S_->GetRecordSetW(key).set_subfieldnames(aqueous_comp_names_);
    }

    // also make sure we can recover previous value of tcc
    requireAtCurrent(key_, tag_current_, *S_, passwd_);
  }

  // require minerals
  if (number_mineral_components_ > 0) {
    auto keys = std::vector<Key>{mineral_vol_frac_key_, mineral_ssa_key_, mineral_rate_constant_key_};
    for (const auto& key : keys) {
      requireAtNext(key, tag_next_, *S_, true, passwd_)
        .SetMesh(mesh_)
        ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_mineral_components_);
      S_->GetRecordSetW(key).set_subfieldnames(mineral_comp_names_);
    }
  }

  // require aqueous kinetics
  if (number_aqueous_kinetics_ > 0) {
    // -- register field
    requireAtNext(first_order_decay_constant_key_, tag_next_, *S_, true, passwd_)
      .SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_aqueous_kinetics_);
    S_->GetRecordSetW(first_order_decay_constant_key_).set_subfieldnames(aqueous_kinetics_names_);
  }

  // require sorption sites
  if (number_sorption_sites_ > 0) {
    auto keys = std::vector<Key>{sorp_sites_key_, surf_cfsc_key_};
    for (const auto& key : keys) {
      requireAtNext(key, tag_next_, *S_, true, passwd_)
        .SetMesh(mesh_)
        ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_sorption_sites_);
      S_->GetRecordSetW(key).set_subfieldnames(sorption_site_names_);
    }
  }

  if (using_sorption_) {
    requireAtNext(total_sorbed_key_, tag_next_, *S_, true, passwd_)
      .SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_aqueous_components_);
    S_->GetRecordSetW(total_sorbed_key_).set_subfieldnames(aqueous_comp_names_);

    if (using_sorption_isotherms_) {
      auto keys = std::vector<Key>{ isotherm_kd_key_, isotherm_freundlich_n_key_, isotherm_langmuir_b_key_};
      for (const auto& key : keys) {
        requireAtNext(key, tag_next_, *S_, true, passwd_)
          .SetMesh(mesh_)
          ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_aqueous_components_);
        S_->GetRecordSetW(key).set_subfieldnames(aqueous_comp_names_);
      }
    }
  }

  // ion exchange sites
  if (number_ion_exchange_sites_ > 0) {
    auto keys = std::vector<Key>{ ion_exchange_sites_key_, ion_exchange_ref_cation_conc_key_ };
    for (const auto& key : keys) {
      requireAtNext(key, tag_next_, *S_, true, passwd_)
        .SetMesh(mesh_)
        ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_ion_exchange_sites_);
      S_->GetRecordSetW(key).set_subfieldnames(ion_exchange_site_names_);
    }
  }
}


/* ******************************************************************
* Most things are initialized through State, but State can only manage that
* if they are always initialized.  If sane defaults are available, or they
* can be derived from other initialized quantities, they are initialized
* here, where we can manage that logic.
******************************************************************* */
void
Chemistry_PK::Initialize()
{
  // Aqueous species
  if (number_aqueous_components_ > 0) {
    if (!S_->GetRecordW(key_, passwd_).initialized()) {
      initializeCVField(*S_, *vo_, key_, tag_next_, passwd_, 0.0);
    }
    initializeCVField(*S_, *vo_, primary_activity_coeff_key_, tag_next_, passwd_, 1.0);
    initializeCVField(*S_, *vo_, free_ion_species_key_, tag_next_, passwd_, 1.0e-9);

    // Sorption sites: all will have a site density, but we can default to zero
    if (using_sorption_) {
      initializeCVField(*S_, *vo_, total_sorbed_key_, tag_next_, passwd_, 0.0);
    }

    // Sorption isotherms: Kd required, Langmuir and Freundlich optional
    if (using_sorption_isotherms_) {
      initializeCVField(*S_, *vo_, isotherm_kd_key_, tag_next_, passwd_, -1.0);
      initializeCVField(*S_, *vo_, isotherm_freundlich_n_key_, tag_next_, passwd_, 1.0);
      initializeCVField(*S_, *vo_, isotherm_langmuir_b_key_, tag_next_, passwd_, 1.0);
    }
  }

  // Minerals: vol frac and surface areas
  if (number_mineral_components_ > 0) {
    initializeCVField(*S_, *vo_, mineral_vol_frac_key_, tag_next_, passwd_, 0.0);
    initializeCVField(*S_, *vo_, mineral_ssa_key_, tag_next_, passwd_, 1.0);
  }

  // Aqueous kinetics
  if (number_aqueous_kinetics_ > 0) {
    initializeCVField(*S_, *vo_, first_order_decay_constant_key_, tag_next_, passwd_, 0.0);
  }

  // Ion exchange sites: default to 1
  if (number_ion_exchange_sites_ > 0) {
    initializeCVField(*S_, *vo_, ion_exchange_sites_key_, tag_next_, passwd_, 1.0);
    initializeCVField(*S_, *vo_, ion_exchange_ref_cation_conc_key_, tag_next_, passwd_, 1.0);
  }

  if (number_sorption_sites_ > 0) {
    initializeCVField(*S_, *vo_, sorp_sites_key_, tag_next_, passwd_, 1.0);
    initializeCVField(*S_, *vo_, surf_cfsc_key_, tag_next_, passwd_, 1.0);
  }

}




/* *******************************************************************
 * Take a timestep from t_old to t_new
 ******************************************************************* */
bool
Chemistry_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  double dt = t_new - t_old;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_LOW))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << t_old
               << " t1 = " << t_new << " h = " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;

  AMANZI_ASSERT(std::abs(S_->get_time(tag_current_) - t_old) < 1.e-4);
  AMANZI_ASSERT(std::abs(S_->get_time(tag_next_) - t_new) < 1.e-4);
  State_to_Solution(Tags::NEXT, *solution_);

  // Get the number of owned (non-ghost) cells for the mesh.
  int num_cells =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  // Ensure dependencies are filled
  S_->GetEvaluator(poro_key_, tag_current_).Update(*S_, name_);
  S_->GetEvaluator(fluid_den_key_, tag_current_).Update(*S_, name_);
  S_->GetEvaluator(saturation_key_, tag_current_).Update(*S_, name_);
  Epetra_MultiVector& aqueous_components =
    *S_->GetW<CompositeVector>(key_, tag_next_, passwd_).ViewComponent("cell", false);

  // Now loop through all the cells and advance the chemistry.
  int max_itrs(0), imax(-1);
  int convergence_failure = 0;
  for (int cell = 0; cell < num_cells; ++cell) {
    int num_itrs = AdvanceSingleCell_(dt, aqueous_components, cell);
    if (num_itrs >= 0) {
      if (max_itrs < num_itrs) {
        max_itrs = num_itrs;
        imax = cell;
      }
    } else {
      // Convergence failure. Compute the next timestep size.
      convergence_failure = 1;
      break;
    }
  }

  // broadcast and check for global failed step
  bool failed = CheckForError_(convergence_failure, max_itrs, imax);

  // Compute the next global timestep.
  dt_next_ = timestep_controller_->getTimestep(dt, max_itrs, !failed);
  return failed;
}



void
Chemistry_PK::CommitStep(double t_old, double t_new, const Tag& tag_next)
{
  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;
  CopyFields_(tag_current_, tag_next);
}


/* ******************************************************************
* Process names of materials
******************************************************************* */
void
Chemistry_PK::parseMinerals_(Teuchos::ParameterList& plist,
                             Teuchos::ParameterList& ic_list)
{
  if (plist.isParameter("minerals")) {
    mineral_comp_names_ = plist.get<Teuchos::Array<std::string>>("minerals").toVector();
  }
  number_mineral_components_ = mineral_comp_names_.size();
}


/* ******************************************************************
* Process names of sorption sites
* NOTE: Do we need to worry about sorption sites?
******************************************************************* */
void
Chemistry_PK::parseSorptionSites_(Teuchos::ParameterList& plist,
        Teuchos::ParameterList& ic_list)
{
  if (plist.isParameter("sorption sites")) {
    sorption_site_names_ = plist.get<Teuchos::Array<std::string>>("sorption sites").toVector();
  }
  number_sorption_sites_ = sorption_site_names_.size();
  using_sorption_ = (number_sorption_sites_ > 0);

  // check if there is an initial condition for ion_exchange_sites
  number_ion_exchange_sites_ = 0;
  using_sorption_isotherms_ = false;

  if (ic_list.isSublist(ion_exchange_sites_key_)) {
    // there is currently only at most one site...
    using_sorption_ = true;
    number_ion_exchange_sites_ = 1;
  }

  if (ic_list.isSublist(isotherm_kd_key_)) {
    using_sorption_ = true;
    using_sorption_isotherms_ = true;
  }

  KeyTriple split;
  bool is_ds = Keys::splitDomainSet(isotherm_kd_key_, split);
  if (is_ds) {
    Key lifted_key = Keys::getKey(std::get<0>(split), "*", std::get<2>(split));

    if (ic_list.isSublist(lifted_key)) {
      using_sorption_ = true;
      using_sorption_isotherms_ = true;
    }
  }

  if (ic_list.isSublist(sorp_sites_key_)) {
    using_sorption_ = true;
  }
}


/* *******************************************************************
* I/O or error messages
******************************************************************* */
bool
Chemistry_PK::CheckForError_(int ierr, int max_itrs, int max_itrs_cell) const
{
  int ierr_l(ierr);
  mesh_->getComm()->MaxAll(&ierr_l, &ierr, 1);
  if (ierr) return true;

  ENorm_t itrs_l{ (double)max_itrs,
    mesh_->getMap(AmanziMesh::Entity_kind::CELL, false).GID(max_itrs_cell) };
  ENorm_t itrs_g;
  ierr = commMaxValLoc(*mesh_->getComm(), itrs_l, itrs_g);
  if (ierr) return true;

  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "max Newton iterations: " << (int) itrs_g.value << " in cell " << itrs_g.gid << std::endl;
  }
  return false;
}


void
Chemistry_PK::CopyFields_(const Tag& tag_dest, const Tag& tag_source) const
{
  std::vector<Key> fields{ mineral_vol_frac_key_,
    mineral_ssa_key_,
    sorp_sites_key_,
    surf_cfsc_key_,
    total_sorbed_key_,
    isotherm_kd_key_,
    isotherm_freundlich_n_key_,
    isotherm_langmuir_b_key_,
    free_ion_species_key_,
    primary_activity_coeff_key_,
    ion_exchange_sites_key_,
    ion_exchange_ref_cation_conc_key_ };

  for (const auto& field : fields) {
    assign(field, tag_dest, tag_source, *S_);
  }
}

} // namespace AmanziChemistry
} // namespace Amanzi

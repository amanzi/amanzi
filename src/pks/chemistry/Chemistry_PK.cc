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

#include "message.hh"

// Chemistry
#include "Chemistry_PK.hh"

namespace Amanzi {
namespace AmanziChemistry {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* Default constructor that initializes all pointers to NULL
****************************************************************** */
Chemistry_PK::Chemistry_PK()
  : passwd_("state"),
    number_minerals_(0),
    number_aqueous_kinetics_(0),
    number_sorption_sites_(0),
    using_sorption_(false),
    using_sorption_isotherms_(false),
    number_ion_exchange_sites_(0),
    dt_(9.9e+9),
    dt_max_(9.9e9){};


Chemistry_PK::Chemistry_PK(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& glist,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, glist, S, soln),
    PK_Physical(pk_tree, glist, S, soln),
    passwd_("state"),
    glist_(glist),
    number_minerals_(0),
    number_aqueous_kinetics_(0),
    number_sorption_sites_(0),
    using_sorption_(false),
    using_sorption_isotherms_(false),
    number_ion_exchange_sites_(0),
    dt_(9.9e+9),
    dt_max_(9.9e9){};


/* ******************************************************************
* Register fields and evaluators with the State
******************************************************************* */
void
Chemistry_PK::Setup()
{
  saturation_tolerance_ = plist_->get<double>("saturation tolerance", 1e-14);
  bool amanzi_physics = plist_->isSublist("physical models and assumptions");

  // require flow fields
  S_->Require<CV_t, CVS_t>(poro_key_, Tags::DEFAULT, poro_key_)
    .SetMesh(mesh_)
    ->SetGhosted(false)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator(poro_key_, Tags::DEFAULT);

  if (!S_->HasRecord(saturation_key_, Tags::DEFAULT)) {
    S_->Require<CV_t, CVS_t>(saturation_key_, Tags::DEFAULT, saturation_key_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S_->RequireEvaluator(saturation_key_, Tags::DEFAULT);
  }
  // REMOVE after #646 is fixed
  // meanwhile we check for the physics module via presence of a particular sublist
  if (!S_->HasRecord(saturation_key_, Tags::CURRENT) && !amanzi_physics) {
    S_->Require<CV_t, CVS_t>(saturation_key_, Tags::CURRENT, saturation_key_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S_->RequireEvaluator(saturation_key_, Tags::CURRENT);
  }

  S_->Require<CV_t, CVS_t>(fluid_den_key_, Tags::DEFAULT, fluid_den_key_)
    .SetMesh(mesh_)
    ->SetGhosted(false)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator(fluid_den_key_, Tags::DEFAULT);

  // require transport fields
  std::vector<std::string>::const_iterator it;
  if (!S_->HasRecord(tcc_key_)) {
    S_->Require<CV_t, CVS_t>(tcc_key_, tag_next_, passwd_, comp_names_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_aqueous_components_);
  }

  // require minerals
  if (number_minerals_ > 0) {
    S_->Require<CV_t, CVS_t>(min_vol_frac_key_, tag_next_, passwd_, mineral_names_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_minerals_);

    S_->Require<CV_t, CVS_t>(min_ssa_key_, tag_next_, passwd_, mineral_names_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_minerals_);

    S_->Require<CV_t, CVS_t>(mineral_rate_constant_key_, tag_next_, passwd_, mineral_names_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_minerals_);
  }

  // require aqueous kinetics
  if (number_aqueous_kinetics_ > 0) {
    std::vector<std::string> subfield_names;

    for (it = aqueous_kinetics_names_.begin(); it != aqueous_kinetics_names_.end(); ++it) {
      subfield_names.push_back(*it + std::string(" aq kin rate cnst"));
    }

    // -- register field
    S_->Require<CV_t, CVS_t>(first_order_decay_constant_key_, tag_next_, passwd_, subfield_names)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_aqueous_kinetics_);
  }

  // require sorption sites
  if (number_sorption_sites_ > 0) {
    std::vector<std::string> ss_names_cv, scfsc_names_cv;
    for (it = sorption_site_names_.begin(); it != sorption_site_names_.end(); ++it) {
      ss_names_cv.push_back(*it + std::string(" sorption site"));
      scfsc_names_cv.push_back(*it + std::string(" surface complex free site conc"));
    }

    // -- register two fields
    S_->Require<CV_t, CVS_t>(sorp_sites_key_, tag_next_, passwd_, ss_names_cv)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_sorption_sites_);
    S_->Require<CV_t, CVS_t>(surf_cfsc_key_, tag_next_, passwd_, scfsc_names_cv)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_sorption_sites_);
  }

  if (using_sorption_) {
    S_->Require<CV_t, CVS_t>(total_sorbed_key_, tag_next_, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_aqueous_components_);

    if (using_sorption_isotherms_) {
      S_->Require<CV_t, CVS_t>(isotherm_kd_key_, tag_next_, passwd_)
        .SetMesh(mesh_)
        ->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_aqueous_components_);

      S_->Require<CV_t, CVS_t>(isotherm_freundlich_n_key_, tag_next_, passwd_)
        .SetMesh(mesh_)
        ->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_aqueous_components_);

      S_->Require<CV_t, CVS_t>(isotherm_langmuir_b_key_, tag_next_, passwd_)
        .SetMesh(mesh_)
        ->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_aqueous_components_);
    }
  }

  // ion
  if (number_aqueous_components_ > 0) {
    S_->Require<CV_t, CVS_t>(free_ion_species_key_, tag_next_, passwd_, comp_names_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_aqueous_components_);

    S_->Require<CV_t, CVS_t>(primary_activity_coeff_key_, tag_next_, passwd_, comp_names_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_aqueous_components_);
  }

  if (number_ion_exchange_sites_ > 0) {
    S_->Require<CV_t, CVS_t>(ion_exchange_sites_key_, tag_next_, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_ion_exchange_sites_);

    S_->Require<CV_t, CVS_t>(ion_exchange_ref_cation_conc_key_, tag_next_, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_ion_exchange_sites_);
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
    if (!S_->GetRecordW(tcc_key_, passwd_).initialized()) {
      InitializeCVField(S_, *vo_, tcc_key_, tag_next_, passwd_, 0.0);
    }
    set_aqueous_components(
      S_->GetPtrW<CompositeVector>(tcc_key_, tag_next_, passwd_)->ViewComponent("cell", false));

    InitializeCVField(S_, *vo_, primary_activity_coeff_key_, tag_next_, passwd_, 1.0);
    InitializeCVField(S_, *vo_, free_ion_species_key_, tag_next_, passwd_, 1.0e-9);

    // Sorption sites: all will have a site density, but we can default to zero
    if (using_sorption_) {
      InitializeCVField(S_, *vo_, total_sorbed_key_, tag_next_, passwd_, 0.0);
    }

    // Sorption isotherms: Kd required, Langmuir and Freundlich optional
    if (using_sorption_isotherms_) {
      InitializeCVField(S_, *vo_, isotherm_kd_key_, tag_next_, passwd_, -1.0);
      InitializeCVField(S_, *vo_, isotherm_freundlich_n_key_, tag_next_, passwd_, 1.0);
      InitializeCVField(S_, *vo_, isotherm_langmuir_b_key_, tag_next_, passwd_, 1.0);
    }
  }

  // Minerals: vol frac and surface areas
  if (number_minerals_ > 0) {
    InitializeCVField(S_, *vo_, min_vol_frac_key_, tag_next_, passwd_, 0.0);
    InitializeCVField(S_, *vo_, min_ssa_key_, tag_next_, passwd_, 1.0);
  }

  // Aqueous kinetics
  if (number_aqueous_kinetics_ > 0) {
    InitializeCVField(S_, *vo_, first_order_decay_constant_key_, tag_next_, passwd_, 0.0);
  }

  // Ion exchange sites: default to 1
  if (number_ion_exchange_sites_ > 0) {
    InitializeCVField(S_, *vo_, ion_exchange_sites_key_, tag_next_, passwd_, 1.0);
    InitializeCVField(S_, *vo_, ion_exchange_ref_cation_conc_key_, tag_next_, passwd_, 1.0);
  }

  if (number_sorption_sites_ > 0) {
    InitializeCVField(S_, *vo_, sorp_sites_key_, tag_next_, passwd_, 1.0);
    InitializeCVField(S_, *vo_, surf_cfsc_key_, tag_next_, passwd_, 1.0);
  }

  // auxiliary fields
  InitializeCVField(S_, *vo_, alquimia_aux_data_key_, tag_next_, passwd_, 0.0);

  // miscaleneous controls
  initial_conditions_time_ = plist_->get<double>("initial conditions time", S_->get_time());
}


/* ******************************************************************
* Process names of materials
******************************************************************* */
void
Chemistry_PK::InitializeMinerals(Teuchos::RCP<Teuchos::ParameterList> plist)
{
  mineral_names_.clear();
  if (plist->isParameter("minerals")) {
    mineral_names_ = plist->get<Teuchos::Array<std::string>>("minerals").toVector();
  }

  number_minerals_ = mineral_names_.size();
}


/* ******************************************************************
* Process names of sorption sites
* NOTE: Do we need to worry about sorption sites?
******************************************************************* */
void
Chemistry_PK::InitializeSorptionSites(Teuchos::RCP<Teuchos::ParameterList> plist,
                                      Teuchos::ParameterList& ic_list)
{
  sorption_site_names_.clear();
  if (plist->isParameter("sorption sites")) {
    sorption_site_names_ = plist->get<Teuchos::Array<std::string>>("sorption sites").toVector();
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

  if (ic_list.isSublist(sorp_sites_key_)) { using_sorption_ = true; }

  // in the old version, this was only in the Block sublist... may need work?
  if (plist->isParameter("Cation Exchange Capacity")) {
    using_sorption_ = true;
    number_ion_exchange_sites_ = 1;
  }

  // surely this doesn't need to be in the global plist?
  if (glist_->isSublist("thermodynamic database")) {
    if (glist_->sublist("thermodynamic database").isSublist("isotherms")) {
      using_sorption_ = true;
    }
  }
}


/* *******************************************************************
* I/O or error messages
******************************************************************* */
void
Chemistry_PK::ErrorAnalysis(int ierr, std::string& internal_msg)
{
  int tmp_out[2], tmp_in[2] = { ierr, (int)internal_msg.size() };
  mesh_->getComm()->MaxAll(tmp_in, tmp_out, 2);

  if (tmp_out[0] != 0) {
    // update time control parameters
    num_successful_steps_ = 0;
    dt_next_ /= dt_cut_factor_;

    // get at least one error message
    int n = tmp_out[1];
    int msg_out[n + 1], msg_in[n + 1], m(mesh_->getComm()->MyPID());
    internal_msg.resize(n);

    Errors::encode_string(internal_msg, n, m, msg_in);
    mesh_->getComm()->MaxAll(msg_in, msg_out, n + 1);
    Errors::decode_string(msg_out, n, internal_msg);

    vo_->Write(Teuchos::VERB_HIGH, internal_msg);

    Errors::Message msg(internal_msg);
    Exceptions::amanzi_throw(msg);
  }
}


void
Chemistry_PK::CopyFieldstoNewState(const Teuchos::RCP<State>& S_next)
{
  std::vector<std::string> fields = { min_vol_frac_key_,
                                      min_ssa_key_,
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

  for (int i = 0; i < fields.size(); ++i) {
    if (S_->HasRecord(fields[i]) && S_next->HasRecord(fields[i])) {
      *S_next->GetW<CompositeVector>(fields[i], tag_next_, passwd_).ViewComponent("cell") =
        *S_->Get<CompositeVector>(fields[i]).ViewComponent("cell");
    }
  }
}

} // namespace AmanziChemistry
} // namespace Amanzi

/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry PK

  Base class for chemical process kernels. It should be never
  instantiated.
*/

#ifndef AMANZI_CHEMISTRY_PK_HH_
#define AMANZI_CHEMISTRY_PK_HH_

#include <vector>

// TPLs
#include "Epetra_MultiVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "VerboseObject.hh"

// Amanzi
#ifdef ALQUIMIA_ENABLED
#  include "ChemistryEngine.hh"
#endif
#include "Key.hh"
#include "Mesh.hh"
#include "PK_Physical.hh"
#include "State.hh"

namespace Amanzi {
namespace AmanziChemistry {

class Chemistry_PK : public PK_Physical {
 public:
  Chemistry_PK();
  Chemistry_PK(Teuchos::ParameterList& pk_tree,
               const Teuchos::RCP<Teuchos::ParameterList>& glist,
               const Teuchos::RCP<State>& S,
               const Teuchos::RCP<TreeVector>& soln);
  virtual ~Chemistry_PK() = default;

  // required members for PK interface
  virtual void Setup() override;
  virtual void Initialize() override;

  virtual void set_dt(double dt) override{};
  virtual double get_dt() override { return dt_next_; }

  // Required members for chemistry interface
  // -- output of auxillary cellwise data from chemistry
  virtual Teuchos::RCP<Epetra_MultiVector> extra_chemistry_output_data() = 0;

  // Basic capabilities
  // -- get/set auxiliary tcc vector that now contains only aqueous components.
  Teuchos::RCP<Epetra_MultiVector> aqueous_components() { return aqueous_components_; }
  void set_aqueous_components(Teuchos::RCP<Epetra_MultiVector> tcc) { aqueous_components_ = tcc; }

  // -- process various objects before/during setup phase
  void InitializeMinerals(Teuchos::RCP<Teuchos::ParameterList> plist);
  void InitializeSorptionSites(Teuchos::RCP<Teuchos::ParameterList> plist,
                               Teuchos::ParameterList& ic_list);

  virtual void CopyFieldstoNewState(const Teuchos::RCP<State>& S_next);
  // -- access
#ifdef ALQUIMIA_ENABLED
  Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine() { return chem_engine_; }
#endif

  // -- output of error messages.
  void ErrorAnalysis(int ierr, std::string& internal_msg);
  int num_aqueous_components() { return number_aqueous_components_; }

 protected:
  std::string passwd_;
  Teuchos::RCP<Teuchos::ParameterList> glist_;

  int number_aqueous_components_;
  std::vector<std::string> comp_names_;
  Teuchos::RCP<Epetra_MultiVector> aqueous_components_;

  int number_minerals_;
  std::vector<std::string> mineral_names_;

  int number_aqueous_kinetics_;
  std::vector<std::string> aqueous_kinetics_names_;

  int number_sorption_sites_, number_total_sorbed_;
  std::vector<std::string> sorption_site_names_;
  bool using_sorption_, using_sorption_isotherms_;

  int number_free_ion_, number_ion_exchange_sites_;
  double saturation_tolerance_;

  // names of state fields
  Key tcc_key_;
  Key poro_key_, saturation_key_, temperature_key_;
  Key fluid_den_key_, molar_fluid_den_key_;
  Key min_vol_frac_key_, min_ssa_key_;
  Key sorp_sites_key_;
  Key surf_cfsc_key_;
  Key total_sorbed_key_;
  Key isotherm_kd_key_, isotherm_freundlich_n_key_, isotherm_langmuir_b_key_;
  Key free_ion_species_key_, primary_activity_coeff_key_;
  Key ion_exchange_sites_key_, ion_exchange_ref_cation_conc_key_;
  Key secondary_activity_coeff_key_;
  Key alquimia_aux_data_key_;
  Key mineral_rate_constant_key_;
  Key first_order_decay_constant_key_;

#ifdef ALQUIMIA_ENABLED
  Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine_;
#endif

  // time controls
  int dt_cut_threshold_, dt_increase_threshold_;
  double dt_min_, dt_max_, dt_prev_, dt_next_, dt_cut_factor_, dt_increase_factor_;

  int num_iterations_, num_successful_steps_;
  double initial_conditions_time_;
};

} // namespace AmanziChemistry
} // namespace Amanzi
#endif

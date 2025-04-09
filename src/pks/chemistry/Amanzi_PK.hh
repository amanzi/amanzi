/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*!

The Amanzi chemistry process kernel uses the following parameters.

.. admonition:: chemistry_amanzi-spec

  * `"thermodynamic database`" ``[list]``

    * `"file`" ``[string]`` is the name of the chemistry database file, relative to the execution directory.

    * `"format`" ``[string]`` is the format of the database file. Actual database format is not XML and
      is the same as described for the 2010 demo with additions for the new chemical processes.
      Valid values: "simple".

  * `"minerals`" ``[Array(string)]`` is the list of mineral names.

  * `"sorption sites`" ``[Array(string)]`` is the list of sorption sites.

  * `"activity model`" ``[string]`` is the type of model used for activity corrections.
    Valid options are `"unit`", `"debye-huckel`", and `"pitzer-hwm`",

  * `"tolerance`" ``[double]`` defines tolerance in Newton solves inside the chemistry library.

  * `"maximum Newton iterations`" ``[int]`` is the maximum number of iteration the chemistry
    library can take.

  * `"auxiliary data`" ``[Array(string)]`` defines additional chemistry related data that the user
    can request be saved to vis files. Currently `"pH`" is the only variable supported.

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="_CHEMISTRY">
    <ParameterList name="thermodynamic database">
      <Parameter name="file" type="string" value="_TRITIUM.bgd"/>
      <Parameter name="format" type="string" value="simple"/>
    </ParameterList>
    <Parameter name="activity model" type="string" value="unit"/>
    <Parameter name="tolerance" type="double" value="1.5e-12"/>
    <Parameter name="maximum Newton iterations" type="int" value="25"/>
    <Parameter name="max timestep (s)" type="double" value="1.5e+07"/>
    <Parameter name="auxiliary data" type="Array(string)" value="{pH}"/>
    <Parameter name="number of component concentrations" type="int" value="1"/>
    <Parameter name="timestep control method" type="string" value="simple"/>
  </ParameterList>
  </ParameterList>

*/

#ifndef CHEMISTRY_AMANZI_PK_HH_
#define CHEMISTRY_AMANZI_PK_HH_

#include <string>
#include <vector>

// TPLs
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "Beaker.hh"
#include "BeakerFields.hh"
#include "BeakerState.hh"
#include "Chemistry_PK.hh"
#include "Key.hh"
#include "PK_Factory.hh"
#include "Mesh.hh"
#include "TreeVector.hh"

namespace Amanzi {
namespace AmanziChemistry {

// Trilinos based chemistry process kernel for the unstructured mesh
class Amanzi_PK : public Chemistry_PK {
 public:
  Amanzi_PK(Teuchos::ParameterList& pk_tree,
            const Teuchos::RCP<Teuchos::ParameterList>& glist,
            const Teuchos::RCP<State>& S,
            const Teuchos::RCP<TreeVector>& soln);

  // members required by PK interface
  virtual void parseParameterList() override;
  virtual void Setup() override final;
  virtual void Initialize() override final;

  virtual void CalculateDiagnostics(const Tag& tag) override final { extra_chemistry_output_data(); }

  // The following two routines provide the interface for
  // output of auxillary cellwise data from chemistry
  Teuchos::RCP<Epetra_MultiVector> extra_chemistry_output_data();
  void set_chemistry_output_names(std::vector<std::string>* names);

  // functions used in Transport PK
  void CopyCellStateToBeakerState(int c);

  // access
  std::shared_ptr<Beaker> get_engine() { return chem_; }
  const BeakerParameters& beaker_parameters() const { return beaker_parameters_; }
  BeakerState beaker_state() { return beaker_state_; }

  void InitializeMinerals(Teuchos::RCP<Teuchos::ParameterList> plist);
  void InitializeSorptionSites(Teuchos::RCP<Teuchos::ParameterList> plist, Teuchos::ParameterList& ic_list);

  void ErrorAnalysis(int ierr, std::string& internal_msg);

  virtual void updateSubstate() override {}

 private:
  void AllocateAdditionalChemistryStorage_();

  void XMLParameters();
  void SetupAuxiliaryOutput();

  void InitializeBeakerFields_();

  void CopyBeakerStructuresToCellState(int c);
  virtual int advanceSingleCell_(int cell, double dt) override;
  virtual void copyFields_(const Tag& tag_dest, const Tag& tag_source) override {}

 protected:
  Key tcc_key_;
  Key poro_key_;
  Key saturation_key_;
  Key fluid_den_key_;
  Key temperature_key_;
  Key min_vol_frac_key_;
  Key min_ssa_key_;
  Key surface_site_density_key_;
  Key surf_cfsc_key_;
  Key total_sorbed_key_;
  Key isotherm_kd_key_;
  Key isotherm_freundlich_n_key_;
  Key isotherm_langmuir_b_key_;
  Key free_ion_species_key_;
  Key primary_activity_coeff_key_;
  Key cation_exchange_capacity_key_;
  Key ion_exchange_ref_cation_conc_key_;
  Key secondary_activity_coeff_key_;
  Key alquimia_aux_data_key_;
  Key mineral_rate_constant_key_;
  Key first_order_decay_rate_constant_key_;

 private:
  std::shared_ptr<Beaker> chem_;
  BeakerParameters beaker_parameters_;
  BeakerState beaker_state_, beaker_state_copy_;
  BeakerFields bf_;

  Key prev_saturation_key_; // move to base class ???

  double dt_int_, dt_global_; // interpolation and global times

  bool using_sorption_;
  bool using_sorption_isotherms_;
  int number_free_ion_, number_total_sorbed_;
  int number_ion_exchange_sites_, number_sorption_sites_;
  int number_aqueous_kinetics_;
  std::vector<std::string> sorption_site_names_;
  std::vector<std::string> aqueous_kinetics_names_;

  Teuchos::RCP<Teuchos::ParameterList> glist_;
  std::vector<std::string> aux_names_;
  std::vector<int> aux_index_;
  Teuchos::RCP<Epetra_MultiVector> aux_data_;

  int ncells_owned_;

  Teuchos::RCP<Epetra_MultiVector> aqueous_components_;

 private:
  // factory registration
  static RegisteredPKFactory<Amanzi_PK> reg_;
};

} // namespace AmanziChemistry
} // namespace Amanzi
#endif

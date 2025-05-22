/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*!

The Alquimia PK leverages the Alquimia interface to implement chemical
reactions by calling out to a chemical engine.

This PK, when not coupled, solves the following set of ODEs:

.. math::
   \frac{\partial C_ij}{\partial t} = \sum_{k=1..nreactions} R_ik(C_{i=1..ncomponents,j})

for :math:`k` indexing reactions, :math:`i` indexing components, and :math:`j`
indexing grid cells.  The total component concentrations :math:`C` are in `[kg L^-1]`.

The details of how to parameterize and specify a reaction network are specific
to each geochemical engine, of which currently only PFloTran and Crunch are
supported.  Corresponding input files for these engines must be provided.
Examples are available in Amanzi's user guide tests and ATS's regression tests.

`"PK type`" = `"chemistry alquimia`"

.. _pk-chemistry-alquimia-spec
.. admonition:: pk-chemistry-alquimia-spec

   * `"timestep controller`" ``[timestep-controller-typed-spec]`` Timestep control parameters.

   * `"engine`" ``[string]`` One of:

     - `"CrunchFlow`"
     - `"PfloTran`"

   * `"engine input file`" ``[string]`` Path to the engine's input file.

   KEYS:

   - `"total component concentration`"
   - `"mass density liquid`"
   - `"porosity`"
   - `"saturation liquid`"
   - `"temperature`"
   - `"mineral volume fractions`"
   - `"mineral specific surface area`"
   - `"mineral rate constant`"
   - `"surface site density`"
   - `"total_sorbed`"
   - `"isotherm kd`"
   - `"isotherm freundlich n`"
   - `"isotherm langmuir b`"
   - `"first order decay ratae constant`"
   - `"cation exchange capacity`"
   - `"aux data`"
   - `"pH`"


See the Alquimia documentation for the needed input spec for a given engine.    

 .. code-block:: xml

   <ParameterList>  <!-- parent list -->
   <ParameterList name="alquimia PK">
     <Parameter name="engine" type="string" value="PFloTran"/>
     <Parameter name="engine input file" type="string" value="tritium.in"/>
     <ParameterList name="timestep controller">
     <Parameter name="timestep control type" type="string" value="fixed"/>
     <ParameterList name="fixed parameters">
       <Parameter name="timestep [s]" type="double" value="1.0e-02"/>
     </ParameterList>
   </ParameterList>

*/

#ifndef AMANZI_CHEMISTRY_ALQUIMIA_PK_HH_
#define AMANZI_CHEMISTRY_ALQUIMIA_PK_HH_

#include <map>
#include <string>
#include <vector>

// TPLs
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "ChemistryEngine.hh"
#include "PK_Factory.hh"
#include "TreeVector.hh"

// Chemistry PK
#include "Chemistry_PK.hh"

namespace Amanzi {
namespace AmanziChemistry {

namespace Impl {

//
// This struct is used to avoid doing State->Get<>() calls within the inner
// (cell-based) loop.
//
// NOTE: ATS at least needs both old and new values of these vectors, under the
// possiblity where chemistry succeeds but something else fails afterwards
// (e.g. transport fails due to invalid timestep size), and we therefore need
// to recover the old chemistry values.
//
// It may be worth looking into whether this can actually happen in our
// assorted coupling schemes.  In particular, I think it is possible that this
// cannot happen if transport+chemistry is subcycled relative to flow (but it
// probably can if they are not subcycled) relative to flow.
//
// If it can't happen, then we don't need to be able to back up, and so we
// don't need to store duplicate copies in memory.
//
// We'll support it for now under the assumption it is possible, and then can
// always construct both the const and non-const versions of these pointers
// with the same vector/tag.  --ETC
//
// NOTE: names here match names in Amanzi, not in Alquimia.  Translation
// between Amanzi and Alquimia happens in the Copy{To,From}Alquimia() methods.
//
struct AlquimiaSubstate {
  // alquimia state
  Epetra_MultiVector const * mass_density_liquid;  // [kg m^-3]
  Epetra_MultiVector const * porosity; // [-]
  Epetra_MultiVector const * temperature; // [K]

  Epetra_MultiVector const * tcc_old; // [mol L^-1]
  Epetra_MultiVector * tcc_new;

  Epetra_MultiVector const * total_sorbed_old; // [mol m^-3 (bulk)]
  Epetra_MultiVector * total_sorbed_new;

  Epetra_MultiVector const * mineral_volume_fraction_old; // [-]
  Epetra_MultiVector * mineral_volume_fraction_new;
  Epetra_MultiVector const * mineral_specific_surface_area_old; // [m^2 (mineral) m^-3 (bulk)]
  Epetra_MultiVector * mineral_specific_surface_area_new;

  Epetra_MultiVector const * surface_site_density_old; // [mol m^-3 (bulk)]
  Epetra_MultiVector * surface_site_density_new;

  Epetra_MultiVector const * cation_exchange_capacity_old; // [mol m^-3 (bulk)]
  Epetra_MultiVector * cation_exchange_capacity_new; // [mol m^-3 (bulk)]

  // alquimia properties -- input only, so only at old time
  Epetra_MultiVector const * saturation_liquid; // [-]

  Epetra_MultiVector const * isotherm_kd; // [kg H2O m^-3 (bulk)]
  Epetra_MultiVector const * isotherm_freundlich_n; // [-]
  Epetra_MultiVector const * isotherm_langmuir_b; // [-]

  Epetra_MultiVector const * mineral_rate_constant; // [mol m^-2 s^-1]
  Epetra_MultiVector const * first_order_decay_rate_constant; // [s^-1]

  // aux data
  Epetra_MultiVector const * aux_data_old;
  Epetra_MultiVector * aux_data_new;

  // aux output -- note, no need for "old" variants of these, output only
  Epetra_MultiVector * pH; // [-]
  Epetra_MultiVector * mineral_saturation_index; // [mol s^-1 m^-3]
  Epetra_MultiVector * mineral_reaction_rate; // [mol s^-1 m^-3 (bulk)]
  Epetra_MultiVector * aqueous_kinetic_rate; // [??]
  Epetra_MultiVector * primary_free_ion_concentration; // [mol L^-1]
  Epetra_MultiVector * primary_activity_coefficient; // [-]
  Epetra_MultiVector * secondary_free_ion_concentration; // [mol L^-1]
  Epetra_MultiVector * secondary_activity_coefficient; // [-]
  // Epetra_MultiVector * gas_partial_pressure; // [-]

  AlquimiaSubstate()
    : mass_density_liquid(nullptr),
      porosity(nullptr),
      temperature(nullptr),
      tcc_old(nullptr),
      tcc_new(nullptr),
      total_sorbed_new(nullptr),
      total_sorbed_old(nullptr),
      mineral_volume_fraction_old(nullptr),
      mineral_volume_fraction_new(nullptr),
      mineral_specific_surface_area_old(nullptr),
      mineral_specific_surface_area_new(nullptr),
      surface_site_density_old(nullptr),
      surface_site_density_new(nullptr),
      cation_exchange_capacity_old(nullptr),
      cation_exchange_capacity_new(nullptr),
      saturation_liquid(nullptr),
      isotherm_kd(nullptr),
      isotherm_freundlich_n(nullptr),
      isotherm_langmuir_b(nullptr),
      mineral_rate_constant(nullptr),
      first_order_decay_rate_constant(nullptr),
      aux_data_old(nullptr),
      aux_data_new(nullptr),
      pH(nullptr),
      mineral_saturation_index(nullptr),
      mineral_reaction_rate(nullptr),
      aqueous_kinetic_rate(nullptr),
      primary_free_ion_concentration(nullptr),
      primary_activity_coefficient(nullptr),
      secondary_free_ion_concentration(nullptr),
      secondary_activity_coefficient(nullptr) {}
};

} // namespace Impl


//
// This is just shorthand for passing all four Alquimia data structures
//
struct AlquimiaBeaker {
  AlquimiaState state;
  AlquimiaProperties properties;
  AlquimiaAuxiliaryData aux_data;
  AlquimiaAuxiliaryOutputData aux_output;
};



#ifdef ALQUIMIA_ENABLED
class Alquimia_PK : public Chemistry_PK {
 public:
  Alquimia_PK(Teuchos::ParameterList& pk_tree,
              const Teuchos::RCP<Teuchos::ParameterList>& glist,
              const Teuchos::RCP<State>& S,
              const Teuchos::RCP<TreeVector>& soln);

  ~Alquimia_PK();

  // members required by PK interface
  virtual void parseParameterList() override;
  virtual void Setup() override final;
  virtual void Initialize() override final;

  Teuchos::RCP<AmanziChemistry::ChemistryEngine> getChemEngine() {
    return chem_engine_;
  }
  void copyToAlquimia(int cell, AlquimiaBeaker& beaker);
  void updateSubstate(const Tag& water_tag) override;

 protected:
  void XMLParameters();

  // Copies the chemistry state to or from the given cell to the Alquimia containers.
  void copyFromAlquimia_(int cell);

  // customization implementations
  void copyFields_(const Tag& tag_dest, const Tag& tag_source) override;
  int advanceSingleCell_(int cell, double dt) override;

  // initialize the cell using a given condition
  int initializeSingleCell_(int cell, const std::string& condition);

  void ParseChemicalConditionRegions_(const Teuchos::ParameterList& param_list,
          std::map<std::string, std::string>& conditions);


 private:
  // Amanzi variable names:
  // -- water state
  Key dens_key_;
  Key poro_key_;
  Key temp_key_;
  Key sat_key_;

  // -- internal chemistry state
  Key total_sorbed_key_;
  Key mineral_volume_fraction_key_;
  Key mineral_specific_surface_area_key_;
  Key surface_site_density_key_;
  Key cation_exchange_capacity_key_;
  Key aux_data_key_;

  // -- chemistry parameters
  Key isotherm_kd_key_;
  Key isotherm_freundlich_n_key_;
  Key isotherm_langmuir_b_key_;
  Key mineral_rate_constant_key_;
  Key first_order_decay_rate_constant_key_;

  // -- aux output (diagnostics)
  Key pH_key_;
  Key mineral_sat_index_key_;
  Key mineral_reaction_rate_key_;
  Key aqueous_kinetic_rate_key_;
  Key primary_ion_conc_key_;
  Key primary_activity_coef_key_;
  Key secondary_ion_conc_key_;
  Key secondary_activity_coef_key_;
  // Key gas_partial_pressure_key_;

  // sizes and metadata
  int number_aqueous_kinetics_;
  std::vector<std::string> aqueous_kinetics_names_;

  int number_sorption_sites_;
  std::vector<std::string> sorption_site_names_;

  int number_ion_exchange_sites_;
  std::vector<std::string> ion_exchange_site_names_;

  int number_isotherm_species_;
  int number_aux_data_;

  std::map<std::string, Key> aux_out_names_;
  std::map<std::string, std::vector<std::string> > aux_out_subfield_names_;

  Impl::AlquimiaSubstate substate_;

  // Alquimia data structures
  AlquimiaBeaker beaker_;
  Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine_;
  bool chem_initialized_;


  // Mapping of region names to geochemical condition names. A region is identified
  // by a string, and all cells within a region will have all geochemical
  // conditions in the corresponding condition vector applied to them.
  std::map<std::string, std::string> chem_initial_conditions_;

 private:
  // factory registration
  static RegisteredPKFactory<Alquimia_PK> reg_;
};
#endif

} // namespace AmanziChemistry
} // namespace Amanzi
#endif

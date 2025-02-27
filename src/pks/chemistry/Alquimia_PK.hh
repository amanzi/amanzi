/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*!

The Alquimia chemistry process kernel only requires the *Engine* and *Engine Input File*
entries, but will also accept and respect the value given for *max timestep (s)*.
Most details are provided in the trimmed PFloTran file *1d-tritium-trim.in*.

.. admonition:: alquimia-spec

  * `"minerals`" ``[Array(string)]`` is the list of mineral names.

  * `"sorption sites`" ``[Array(string)]`` is the list of sorption sites.

  * `"auxiliary data`" ``[Array(string)]`` defines additional chemistry related data that the user
    can request be saved to vis files.

  * `"min timestep (s)`" ``[double]`` is the minimum timestep that chemistry will allow 
    the MPC to take.

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="_CHEMISTRY">
    <Parameter name="engine" type="string" value="PFloTran"/>
    <Parameter name="engine input file" type="string" value="_TRITIUM.in"/>
    <Parameter name="minerals" type="Array(string)" value="{quartz, kaolinite, goethite, opal}"/>
    <Parameter name="min timestep (s)" type="double" value="1.5778463e-07"/>
    <Parameter name="max timestep (s)" type="double" value="1.5778463e+07"/>
    <Parameter name="initial timestep (s)" type="double" value="1.0e-02"/>
    <Parameter name="timestep control method" type="string" value="simple"/>
    <Parameter name="timestep cut threshold" type="int" value="8"/>
    <Parameter name="timestep cut factor" type="double" value="2.0"/>
    <Parameter name="timestep increase threshold" type="int" value="4"/>
    <Parameter name="timestep increase factor" type="double" value="1.2"/>
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
// NOTE: ATS at least needs both old and new values of these vectors, under the
// possiblity where chemistry succeeds but transport fails (due to invalid
// timestep size), and we therefore need to recover old chemistry?
//
// It may be worth looking into whether this can actually happen in our
// assorted coupling schemes.  In particular, I think it is possible that this
// cannot happen if transport+chemistry is subcycled relative to flow (but it
// probably can if they are not subcycled) relative to flow.
//
// At any rate, we'll support it for now under the assumption it is possible,
// and then can always construct both the const and non-const versions of these
// pointers with the same vector/tag.
// --ETC
//
template<class Vector_type>
struct AlquimiaSubstate {
  Epetra_MultiVector const * porosity;
  Epetra_MultiVector const * water_saturation;
  Epetra_MultiVector const * fluid_density;
  Epetra_MultiVector const * temperature;

  Epetra_MultiVector const * aqueous_components_old;
  Epetra_MultiVector* aqueous_components_new_new;
  Epetra_MultiVector const * total_sorbed_old;
  Epetra_MultiVector* total_sorbed_new_new;

  Epetra_MultiVector const * mineral_vol_frac_old;
  Epetra_MultiVector* mineral_vol_frac_new;
  Epetra_MultiVector const * mineral_ssa_old;
  Epetra_MultiVector* mineral_ssa_new;
  Epetra_MultiVector const * mineral_rate_constant_old;
  Epetra_MultiVector* mineral_rate_constant_new;

  Epetra_MultiVector const * free_ion_species_old;
  Epetra_MultiVector* free_ion_species_new;
  Epetra_MultiVector const * ion_exchange_sites_old;
  Epetra_MultiVector* ion_exchange_sites_new;
  Epetra_MultiVector const * sorption_sites_old;
  Epetra_MultiVector* sorption_sites_new;

  Epetra_MultiVector const * isotherm_kd_old;
  Epetra_MultiVector* isotherm_kd_new;
  Epetra_MultiVector const * isotherm_freundlich_n_old;
  Epetra_MultiVector* isotherm_freundlich_n_new;
  Epetra_MultiVector const * isotherm_langmuir_b_old;
  Epetra_MultiVector* isotherm_langmuir_b_new;

  Epetra_MultiVector const * first_order_decay_constant_old;
  Epetra_MultiVector* first_order_decay_constant_new;

  // aux data is only ever copied out, not in
  Epetra_MultiVector* aux_data;
};
} // namespace Impl


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

 protected:
  // Copies the chemistry state to or from the given cell to the Alquimia containers.
  void CopyToAlquimia_(int cell_id);
  void CopyFromAlquimia_(const int cell, bool aux_out);

  // initialize the cell using a given condition
  int InitializeSingleCell_(int cell, const std::string& condition);

  void ParseChemicalConditionRegions_(const Teuchos::ParameterList& param_list,
          std::map<std::string,
          std::string>& conditions);

  // maps
  void InitializeAuxNamesMap_();

  Impl::AlquimiaSubstate createSubState();

 private:
  // Alquimia data structures for interface with Amanzi.
  AlquimiaState alq_state_;
  AlquimiaProperties alq_mat_props_;
  AlquimiaAuxiliaryData alq_aux_data_;
  AlquimiaAuxiliaryOutputData alq_aux_output_;
  Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine_;
  bool chem_initialized_;

  // Mapping of region names to geochemical condition names. A region is identified
  // by a string, and all cells within a region will have all geochemical
  // conditions in the corresponding condition vector applied to them.
  std::map<std::string, std::string> chem_initial_conditions_;

  // Auxiliary output data, requested by and stored within Amanzi.
  std::vector<std::string> aux_names_;
  std::vector<std::vector<std::string>> aux_subfield_names_;
  std::vector<std::vector<int>> aux_subfield_map_;

 private:
  // factory registration
  static RegisteredPKFactory<Alquimia_PK> reg_;
};
#endif

} // namespace AmanziChemistry
} // namespace Amanzi
#endif

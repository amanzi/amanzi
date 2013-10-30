/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi Chemistry

License: see COPYRIGHT
Author: Jeffrey Johnson

This is a point of contact for the chemistry engine exposed by Alquimia to 
the rest of Amanzi--it provides the ability to enforce geochemical conditions 
and to integrate reactions given a chemical configuration.

 ------------------------------------------------------------------------- */

#ifndef AMANZI_CHEMISTRY_ENGINE_HH_
#define AMANZI_CHEMISTRY_ENGINE_HH_

#include <string>
#include <vector>
#include <map>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "alquimia_memory.h"
#include "alquimia_util.h"
#include "alquimia_constants.h"
#include "alquimia_containers.h"
#include "alquimia_interface.h"

namespace Amanzi {
namespace AmanziChemistry {


class Chemistry_Engine {

 public:

  // Constructs a chemistry engine using the Chemistry parameter list in the input file and 
  // Amanzi's chemistry state.
  explicit Chemistry_Engine(const Teuchos::ParameterList& param_list);

  // Destructor.
  ~Chemistry_Engine();

  // Returns the number of primary chemical species in the present simulation.

  // These query methods return metadata for the chemical configuration represented 
  // within the chemistry engine.
  int NumPrimarySpecies() const;
  int NumSorbedSpecies() const;
  void GetPrimarySpeciesNames(std::vector<std::string>& speciesNames) const;
  int NumMinerals() const;
  void GetMineralNames(std::vector<std::string>& mineralNames) const;
  int NumSurfaceSites() const;
  void GetSurfaceSiteNames(std::vector<std::string>& siteNames) const;
  int NumIonExchangeSites() const;
  void GetIonExchangeNames(std::vector<std::string>& ionExchangeNames) const;
  int NumIsothermSpecies() const;
  void GetIsothermSpeciesNames(std::vector<std::string>& speciesNames) const;

  // Enforces the geochemical condition with the given name on the chemical configuration 
  // represented by the given array of concentrations at the given time. The order of the 
  // concentrations in the array matches that of the species names returned by GetSpeciesNames.
  void EnforceCondition(const std::string& conditionName,
                        double time,
                        double* concentrations);

  // Advances the species represented by the given array of concentrations, replacing old values 
  // with new values. The order of the concentrations in the array matches that of the species names 
  // returned by GetSpeciesNames.
  void Advance(double deltaTime,
               double* concentrations);

  // FIXME: How is auxiliary output handled?
 private:

  // Reads parameters from the XML input file to set up chemistry.
  void ReadXMLParameters();

  // Initializes the chemistry engine, preparing it to answer queries from other entities.
  void Initialize();

  // This helper parses geochemical conditions from the given parameter list, placing them into a
  // data structure that will be used to retrieve them by name.
  void ParseChemicalConditions(const Teuchos::ParameterList& param_list,
                               std::map<std::string, AlquimiaGeochemicalCondition*>& conditions);

  // Copies data from given buffers to Alquimia in preparation for chemistry.
  void CopyFromBuffersToAlquimia(const double* component_concentrations,
                                 const double* sorbed_components,
                                 const double* mineral_volume_fractions,
                                 const double* mineral_specific_surface_areas,
                                 const double* cation_exchange_capacity,
                                 const double* sorption_sites,
                                 double water_density,
                                 double porosity,
                                 double volume,
                                 double saturation);

  // Copies data from Alquimia to the given buffers to Alquimia after computation.
  void CopyFromAlquimiaToBuffers(double* component_concentrations,
                                 double* sorbed_components,
                                 double* mineral_volume_fractions,
                                 double* mineral_specific_surface_areas,
                                 double* cation_exchange_capacity,
                                 double* sorption_sites,
                                 double& water_density,
                                 double& porosity,
                                 double& volume,
                                 double& saturation,
                                 double* free_ion_species,
                                 double* primary_activity_coeffs,
                                 double* secondary_activity_coeffs,
                                 double* ion_exchange_ref_cation_concs,
                                 double* surface_complex_free_site_concs,
                                 double& pH) const;

  // parameter lists
  Teuchos::ParameterList main_param_list_, chem_param_list_;

  // Alquimia data structures.
  bool chem_initialized_;
  AlquimiaInterface chem_;
  AlquimiaEngineStatus chem_status_;
  AlquimiaData chem_data_;

  // Mapping of region names to geochemical conditions. A region is identified 
  // by a string, and all cells within a region will have all geochemical 
  // conditions in the corresponding condition vector applied to them. NOTE
  // that these maps do not own the geochemical conditions--they only hold 
  // pointers to the objects.
  std::map<std::string, AlquimiaGeochemicalCondition*> chem_conditions_;
  
  // Vector that takes responsibility for ownership of geochemical conditions.
  std::vector<AlquimiaGeochemicalCondition*> all_chem_conditions_;

  // Back-end engine name and input file.
  std::string chem_engine_inputfile_;
  std::string chem_engine_name_;

  // Names of auxiliary data.
  std::vector<std::string> aux_names_;

  // forbidden.
  Chemistry_Engine();
  Chemistry_Engine(const Chemistry_Engine&);
  Chemistry_Engine& operator=(const Chemistry_Engine&);

};

} // namespace
} // namespace

#endif

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
  Chemistry_Entine(const Teuchos::ParameterList& param_list,
                   Teuchos::RCP<Chemistry_State> chem_state);

  // Destructor.
  ~Chemistry_Engine();

  // Initializes the chemistry engine, preparing it to answer queries from other entities.
  void Initialize();

  // Returns the number of chemical species in the present simulation.
  int NumSpecies() const;

  // Returns an array of names for the chemical species in the present simulation.
  void GetSpeciæsNames(std::vector<std::string>& speciesNames) const;

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

  // parameter lists
  Teuchos::ParameterList main_param_list_, chem_param_list_;

  // This helper parses geochemical conditions from the given parameter list, placing them into a
  // data structure that will be used to retrieve them by name.
  void ParseChemicalConditions(const Teuchos::ParameterList& param_list,
                               std::map<std::string, AlquimiaGeochemicalCondition*>& conditions);

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
  std::map<std::string, AlquimiaGeochemicalCondition*> chem_initial_conditions_;
  std::map<std::string, AlquimiaGeochemicalCondition*> chem_boundary_conditions_;
  
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

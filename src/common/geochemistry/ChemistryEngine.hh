/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi Chemistry

License: see COPYRIGHT
Author: Jeffrey Johnson

This is a point of contact for the chemistry engine exposed by Alquimia to 
the rest of Amanzi--it provides the ability to enforce geochemical conditions 
and to integrate reactions given a chemical configuration.

 ------------------------------------------------------------------------- */

#ifndef AMANZI_ChemistryEngine_HH_
#define AMANZI_ChemistryEngine_HH_

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


// Forward declaration of GeochemicalConditionContext.
class GeochemicalConditionContext;

class ChemistryEngine {

 public:

  // Constructs a chemistry engine using the Chemistry parameter list in the input file and 
  // Amanzi's chemistry state.
  explicit ChemistryEngine(const Teuchos::ParameterList& param_list);

  // Destructor.
  ~ChemistryEngine();

  // Returns the name of the backend that does the chemistry.
  const std::string& Name() const;

  // Returns true if the chemistry engine is thread-safe, false if not.
  bool IsThreadSafe() const;

  // These query methods return metadata for the chemical configuration represented 
  // within the chemistry engine.
  int NumPrimarySpecies() const;
  int NumAqueousComplexes() const;
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
  int NumFreeIonSpecies() const;

  // Returns a reference to a "sizes" object that can be queried to find the sizes of the various 
  // arrays representing the geochemical state within the engine.
  const AlquimiaSizes& Sizes() const;

  // Initializes the data structures that hold the chemical state information.
  void InitState(AlquimiaState& chem_state, 
                 AlquimiaMaterialProperties& mat_props,
                 AlquimiaAuxiliaryData& aux_data,
                 AlquimiaAuxiliaryOutputData& aux_output);
                 
  // Frees the data structures that hold the chemical state information.
  void FreeState(AlquimiaState& chem_state,
                 AlquimiaMaterialProperties& mat_props,
                 AlquimiaAuxiliaryData& aux_data,
                 AlquimiaAuxiliaryOutputData& aux_output);

  // Enforces the geochemical condition with the given name on the chemical configuration 
  // represented by the given array of concentrations at the given time. The order of the 
  // concentrations in the array matches that of the species names returned by GetSpeciesNames.
  void EnforceCondition(const std::string& conditionName,
                        double time,
                        AlquimiaState& chem_state,
                        AlquimiaMaterialProperties& mat_props,
                        AlquimiaAuxiliaryData& aux_data,
                        AlquimiaAuxiliaryOutputData& aux_output);

  // Advances the species represented by the given array of concentrations, replacing old values 
  // with new values. The order of the concentrations in the array matches that of the species names 
  // returned by GetSpeciesNames.
  void Advance(double delta_time,
               AlquimiaState& chem_state,
               AlquimiaMaterialProperties& mat_props,
               AlquimiaAuxiliaryData& aux_data,
               AlquimiaAuxiliaryOutputData& aux_output,
               int& num_iterations);

 private:

  // Reads parameters from the XML input file to set up chemistry.
  void ReadXMLParameters();

  // Initializes the chemistry engine, preparing it to answer queries from other entities.
  void Initialize();

  // This helper parses geochemical conditions from the given parameter list, placing them into a
  // data structure that will be used to retrieve them by name.
  void ParseChemicalConditions(const Teuchos::ParameterList& param_list,
                               std::map<std::string, AlquimiaGeochemicalCondition*>& conditions);

  // parameter lists
  Teuchos::ParameterList main_param_list_, chem_param_list_;

  // Alquimia data structures.
  bool chem_initialized_;
  void* engine_state_;
  AlquimiaEngineFunctionality functionality_;
  AlquimiaSizes sizes_;
  AlquimiaInterface chem_;
  AlquimiaEngineStatus chem_status_;
  AlquimiaProblemMetaData chem_metadata_;

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
  ChemistryEngine();
  ChemistryEngine(const ChemistryEngine&);
  ChemistryEngine& operator=(const ChemistryEngine&);

};

} // namespace
} // namespace

#endif

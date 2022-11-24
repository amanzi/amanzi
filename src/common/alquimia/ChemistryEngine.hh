/*
  Alqumia

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Jeffrey Johnson

  This is a point of contact for the chemistry engine exposed by Alquimia 
  to the rest of Amanzi--it provides the ability to enforce geochemical 
  conditions and to integrate reactions given a chemical configuration.
*/

#ifndef AMANZI_CHEMISTRY_ENGINE_HH_
#define AMANZI_CHEMISTRY_ENGINE_HH_

#include <string>
#include <vector>
#include <map>

#include "alquimia/alquimia_memory.h"
#include "alquimia/alquimia_util.h"
#include "alquimia/alquimia_constants.h"
#include "alquimia/alquimia_containers.h"
#include "alquimia/alquimia_interface.h"

namespace Amanzi {
namespace AmanziChemistry {


struct GeochemicalConditionData {
  bool processed;
  AlquimiaProperties mat_props;
  AlquimiaGeochemicalCondition condition;
  AlquimiaState chem_state;
  AlquimiaAuxiliaryData aux_data;
};

typedef std::map<std::string, GeochemicalConditionData*> GeochemicalConditionMap;

class ChemistryEngine {
 public:
  // Constructs a chemistry engine using the given engine (backend) name and input file.
  ChemistryEngine(const std::string& engineName, const std::string& inputFile);

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
  void GetPrimarySpeciesNames(std::vector<std::string>& species_names) const;
  int NumMinerals() const;
  void GetMineralNames(std::vector<std::string>& mineral_names) const;
  int NumSurfaceSites() const;
  void GetSurfaceSiteNames(std::vector<std::string>& site_names) const;
  int NumIonExchangeSites() const;
  void GetIonExchangeNames(std::vector<std::string>& ion_exchange_names) const;
  int NumIsothermSpecies() const;
  void GetIsothermSpeciesNames(std::vector<std::string>& species_names) const;
  int NumFreeIonSpecies() const;
  void GetAuxiliaryOutputNames(std::vector<std::string>& aux_names,
                               std::vector<std::vector<std::string>>& subfield_names) const;
  int NumAqueousKinetics() const;
  void GetAqueousKineticNames(std::vector<std::string>& kinetics_names) const;

  // Returns a reference to a "sizes" object that can be queried to find the sizes of the various
  // arrays representing the geochemical state within the engine.
  const AlquimiaSizes& Sizes() const;

  // Initializes the data structures that hold the chemical state information.
  void InitState(AlquimiaProperties& mat_props,
                 AlquimiaState& chem_state,
                 AlquimiaAuxiliaryData& aux_data,
                 AlquimiaAuxiliaryOutputData& aux_output);

  // Frees the data structures that hold the chemical state information.
  void FreeState(AlquimiaProperties& mat_props,
                 AlquimiaState& chem_state,
                 AlquimiaAuxiliaryData& aux_data,
                 AlquimiaAuxiliaryOutputData& aux_output);

  // Creates a geochemical condition with the given name within the chemistry engine, using the
  // data in the given containers.
  void CreateCondition(const std::string& condition_name);

  /* Mineral constraints will be discontinued in Alquimia -- see Sergi
  // Adds a mineral constraint to the geochemical condition with the given name. If another 
  // constraint with the same mineral name exists for this condition, it is replaced by 
  // this one.
  void AddMineralConstraint(const std::string& condition_name,
                            const std::string& mineral_name,
                            double volume_fraction,
                            double specific_surface_area);
  Mineral constraints will be discontinued in Alquimia -- see Sergi */

  // Adds an aqueous constraint to the geochemical condition with the given name. If another
  // constraint involving the same primary species exists for this condition, it is replaced by
  // this one.
  // The constraint type may be "total_aqueous", "total_sorb", "free", "mineral", "gas",
  // "pH", or "charge".
  // The associated (mineral) species must be the name of a mineral mentioned
  // in a mineral constraint that has been previously added to this condition
  // using AddMineralConstraint().
  void AddAqueousConstraint(const std::string& condition_name,
                            const std::string& primary_species_name,
                            const std::string& constraint_type,
                            const std::string& associated_species);

  // Enforces the geochemical condition with the given name on the chemical configuration
  // represented by the given array of concentrations at the given time. The order of the
  // concentrations in the array matches that of the species names returned by GetSpeciesNames.
  void EnforceCondition(const std::string& condition_name,
                        const double time,
                        const AlquimiaProperties& mat_props,
                        AlquimiaState& chem_state,
                        AlquimiaAuxiliaryData& aux_data,
                        AlquimiaAuxiliaryOutputData& aux_output);

  // Advances the species represented by the given array of concentrations, replacing old values
  // with new values. The order of the concentrations in the array matches that of the species names
  // returned by GetSpeciesNames. Returns true if the advance is successful,
  // false if it fails.
  bool Advance(const double delta_time,
               const AlquimiaProperties& mat_props,
               AlquimiaState& chem_state,
               AlquimiaAuxiliaryData& aux_data,
               AlquimiaAuxiliaryOutputData& aux_output,
               int& num_iterations);

 private:
  // Alquimia data structures.
  bool chem_initialized_;
  void* engine_state_;
  AlquimiaEngineFunctionality functionality_;
  AlquimiaSizes sizes_;
  AlquimiaInterface chem_;
  AlquimiaEngineStatus chem_status_;
  AlquimiaProblemMetaData chem_metadata_;

  // Mapping of geochemical condition names to geochemical conditions (and flags indicating
  // whether they have been processed).
  GeochemicalConditionMap chem_conditions_;

  // Back-end engine name and input file.
  std::string chem_engine_name_;
  std::string chem_engine_inputfile_;

  // forbidden.
  ChemistryEngine();
  ChemistryEngine(const ChemistryEngine&);
  ChemistryEngine& operator=(const ChemistryEngine&);
};

} // namespace AmanziChemistry
} // namespace Amanzi

#endif

#ifndef _AlquimiaHelper_Structured_H_
#define _AlquimiaHelper_Structured_H_

#ifdef ALQUIMIA_ENABLED
#include <ChemistryEngine.hh>

#include <FArrayBox.H>

#include <vector>
#include <string>

class AlquimiaHelper_Structured
{
public:
  AlquimiaHelper_Structured(Amanzi::AmanziChemistry::ChemistryEngine* engine);

  ~AlquimiaHelper_Structured();

  void Advance(const FArrayBox& aqueous_saturation,       int sSat,
               const FArrayBox& aqueous_pressure,         int sPress,
               const FArrayBox& porosity,                 int sPhi,
               const FArrayBox& volume,                   int sVol,
               FArrayBox&       primary_species_mobile,   int sPrimMob,
               FArrayBox&       primary_species_immobile, int sPrimImmob,
               FArrayBox&       fcnCnt,                   int sFunc,
               FArrayBox&       aux_data, Real water_density, Real temperature,
               const Box& box, Real dt);

  void EnforceCondition(const FArrayBox& aqueous_saturation,       int sSat,
                        const FArrayBox& aqueous_pressure,         int sPress,
                        const FArrayBox& porosity,                 int sPhi,
                        const FArrayBox& volume,                   int sVol,
                        FArrayBox&       primary_species_mobile,   int sPrimMob,
                        FArrayBox&       primary_species_immobile, int sPrimImmob,
                        FArrayBox&       aux_data, Real water_density, Real temperature,
                        const Box& box, const std::string& condition_name, Real time);

  const std::map<std::string,int>& AuxChemVariablesMap() const {return aux_chem_variables;}


protected:
  void BL_to_Alquimia(const FArrayBox& aqueous_saturation,       int sSat,
                      const FArrayBox& aqueous_pressure,         int sPress,
                      const FArrayBox& porosity,                 int sPhi,
                      const FArrayBox& volume,                   int sVol,
                      FArrayBox&       primary_species_mobile,   int sPrimMob,
                      FArrayBox&       primary_species_immobile, int sPrimImmob,
                      FArrayBox&                   aux_data,
                      const IntVect&               iv,
                      Real                         water_density,
                      Real                         temperature,
                      AlquimiaMaterialProperties&  mat_props,
                      AlquimiaState&               chem_state,
                      AlquimiaAuxiliaryData&       aux_input,
                      AlquimiaAuxiliaryOutputData& aux_output);

  void Alquimia_to_BL(FArrayBox& primary_species_mobile,   int sPrimMob,
                      FArrayBox& primary_species_immobile, int sPrimImmob,
                      FArrayBox& aux_data,
                      const IntVect& iv,
                      AlquimiaMaterialProperties& mat_props,
                      AlquimiaState& chem_state,
                      AlquimiaAuxiliaryOutputData& aux_output);

  Amanzi::AmanziChemistry::ChemistryEngine* engine;
  int nthreads;

  std::vector<AlquimiaState> alquimia_state;
  std::vector<AlquimiaMaterialProperties> alquimia_material_properties;
  std::vector<AlquimiaAuxiliaryData> alquimia_aux_in;
  std::vector<AlquimiaAuxiliaryOutputData> alquimia_aux_out;
  std::vector<std::string> mineralNames, surfSiteNames, ionExchangeNames, isothermSpeciesNames;
  int Nmobile, Nimmobile, Nminerals;
  const AlquimiaSizes& alq_sizes;
  int num_aux_ints, ndigits_int, num_aux_doubles, ndigits_doubles, num_free_ion_species, ndigits_freeion;

  std::map<std::string,int> aux_chem_variables;
  bool using_sorption;
};

#endif

#endif

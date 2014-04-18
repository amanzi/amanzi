#ifndef _AlquimiaHelper_Structured_H_
#define _AlquimiaHelper_Structured_H_

#include <ChemistryHelper_Structured.H>

#ifdef ALQUIMIA_ENABLED
#include <ChemistryEngine.hh>

class AlquimiaHelper_Structured
  : public ChemistryHelper_Structured
{
public:
  AlquimiaHelper_Structured(Amanzi::AmanziChemistry::ChemistryEngine* engine);

  virtual ~AlquimiaHelper_Structured();

  virtual void Advance(const FArrayBox& aqueous_saturation,       int sSat,
                       const FArrayBox& aqueous_pressure,         int sPress,
                       const FArrayBox& porosity,                 int sPhi,
                       const FArrayBox& volume,                   int sVol,
                       FArrayBox&       primary_species_mobile,   int sPrimMob,
                       FArrayBox&       fcnCnt,                   int sFunc,
                       FArrayBox&       aux_data, Real water_density, Real temperature,
                       const Box& box, Real dt, bool initialize = false);

  virtual void EnforceCondition(FArrayBox& primary_species_mobile,   int sPrimMob,
                                FArrayBox& auxiliary_data, bool initialize_auxiliary_data,
                                const Box& box, const std::string& condition_name, Real time);

protected:
  void BL_to_Alquimia(const FArrayBox& aqueous_saturation,       int sSat,
                      const FArrayBox& aqueous_pressure,         int sPress,
                      const FArrayBox& porosity,                 int sPhi,
                      const FArrayBox& volume,                   int sVol,
                      FArrayBox&       primary_species_mobile,   int sPrimMob,
                      FArrayBox&                   aux_data,
                      bool                         auxiliary_data_initialized,
                      const IntVect&               iv,
                      Real                         water_density,
                      Real                         temperature,
                      AlquimiaMaterialProperties&  mat_props,
                      AlquimiaState&               chem_state,
                      AlquimiaAuxiliaryData&       aux_input,
                      AlquimiaAuxiliaryOutputData& aux_output);

  void Alquimia_to_BL(FArrayBox& primary_species_mobile,   int sPrimMob,
                      FArrayBox&                   aux_data,
                      bool                         initialize_auxiliary_data,
                      const IntVect&               iv,
                      AlquimiaMaterialProperties&  mat_props,
                      AlquimiaState&               chem_state,
                      AlquimiaAuxiliaryData&       aux_input,
                      AlquimiaAuxiliaryOutputData& aux_output);

  Amanzi::AmanziChemistry::ChemistryEngine* engine;
  int nthreads;

  std::vector<AlquimiaState> alquimia_state;
  std::vector<AlquimiaMaterialProperties> alquimia_material_properties;
  std::vector<AlquimiaAuxiliaryData> alquimia_aux_in;
  std::vector<AlquimiaAuxiliaryOutputData> alquimia_aux_out;
  std::vector<std::string> primarySpeciesNames, mineralNames, surfSiteNames, ionExchangeNames, isothermSpeciesNames;
  int Nmobile, Nimmobile, Nminerals;
  const AlquimiaSizes& alq_sizes;
  int num_ion_exchange, num_aux_ints, ndigits_int, num_aux_doubles, ndigits_doubles, num_free_ion_species, ndigits_freeion;
};

#endif

#endif

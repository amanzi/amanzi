#ifndef _AmanziChemHelper_Structured_H_
#define _AmanziChemHelper_Structured_H_

#include <ChemistryHelper_Structured.H>

#include <simple_thermo_database.hh>

class AmanziChemHelper_Structured
  : public ChemistryHelper_Structured
{
public:
  AmanziChemHelper_Structured(const std::vector<std::string>& _primarySpeciesNames,
                              const std::vector<std::string>& _sorbedPrimarySpeciesNames,
                              const std::vector<std::string>& _mineralNames,
                              const std::vector<std::string>& _surfaceComplexationSiteNames,
                              bool                            _hasCationExchangeCapacity,
                              const std::vector<std::string>& _isothermSpeciesNames,
                              const std::vector<std::string>& _freeIonSpeciesNames,
                              const std::string&              _thermo_database_filename,
                              const std::string&              _thermo_database_format,
                              const std::string&              _activity_model,
                              int                             _verbose);

  AmanziChemHelper_Structured(const std::vector<std::string>& primarySpeciesNames,
                              const std::vector<std::string>& mineralNames);

  virtual ~AmanziChemHelper_Structured();

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
  
  int nthreads;
  PArray<Amanzi::AmanziChemistry::SimpleThermoDatabase> chemSolve;
  Array<Amanzi::AmanziChemistry::Beaker::BeakerComponents> components;
  Array<Amanzi::AmanziChemistry::Beaker::BeakerParameters> parameters;


  int ntracers, nminerals;
  Array<int> guess_comp, activity_comp, sorbed_comp, kd_comp;
  std::map<int,int> freundlich_n_comp;
  std::map<int,int> langmuir_b_comp;
  Array<int> volume_fraction_comp, specific_surface_area_comp;
  int cation_exchange_capacity_comp;

  Array<std::string> primarySpeciesNames,sorbedPrimarySpeciesNames,mineralNames;
  Array<std::string> surfaceComplexationSiteNames,freeIonSpeciesNames,isothermSpeciesNames;

  bool using_sorption, hasCationExchangeCapacity;
  std::string thermo_database_file, thermo_database_format, activity_model;
  int NionExchangeSites,NprimarySpecies,NsorbedPrimarySpecies,Nminerals;
  int NsurfaceComplexationSites,NfreeIonSpecies,NisothermSpecies;
  int SmineralVolumeFraction,SmineralSpecificSurfaceArea,SsurfaceComplexationSiteDensity;
  int SsurfaceComplexationFreeSiteConcentration,ScationExchangeCapacity,SionExchangeSiteDensity;
  int SionExchangeReferenceCationConcentration,SfreeIonConcentration;
  int SactivityCoefficient,Ssorbed,SisothermKd,SisothermFreundlichN;
  int SisothermLangmuirB;

  int verbose;
};

#endif

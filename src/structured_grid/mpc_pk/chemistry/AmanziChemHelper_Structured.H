#ifndef _AmanziChemHelper_Structured_H_
#define _AmanziChemHelper_Structured_H_

#include <ChemistryHelper_Structured.H>

#include <SimpleThermoDatabase.hh>

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
                              const std::string&              _activity_model);

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
                       const Box& box, Real dt, int chem_verbose);

  virtual void EnforceCondition(FArrayBox& primary_species_mobile,   int sPrimMob,
                                FArrayBox& auxiliary_data, Real water_density, Real temperature,
                                const Box& box, const std::string& condition_name, Real time,
				int chem_verbose);

  void Initialize(const FArrayBox& aqueous_saturation,       int sSat,
                  const FArrayBox& aqueous_pressure,         int sPress,
                  const FArrayBox& porosity,                 int sPhi,
                  const FArrayBox& volume,                   int sVol,
                  FArrayBox&       primary_species_mobile,   int sPrimMob,
                  FArrayBox&       fcnCnt,                   int sFunc,
                  FArrayBox&       aux_data, Real water_density, Real temperature,
                  const Box& box);

protected:
  void DumpChemStructures(std::ostream&                                  os,
			  Amanzi::AmanziChemistry::SimpleThermoDatabase& TheChemSolve,
			  Amanzi::AmanziChemistry::BeakerState&      TheComponent,
			  Amanzi::AmanziChemistry::BeakerParameters& TheParameter);

  int nthreads;
  PArray<Amanzi::AmanziChemistry::SimpleThermoDatabase> chemSolve;
  Array<Amanzi::AmanziChemistry::BeakerState> components;
  Array<Amanzi::AmanziChemistry::BeakerParameters> parameters;


  std::string thermo_database_file, thermo_database_format, activity_model;
  int verbose;

  Teuchos::RCP<Amanzi::VerboseObject> vo_;
};

#endif

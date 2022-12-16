/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _ChemistryHelper_Structured_H_
#define _ChemistryHelper_Structured_H_

#include <FArrayBox.H>

#include <map>
#include <vector>
#include <string>

class ChemistryHelper_Structured
{
public:
  ChemistryHelper_Structured() {}

  virtual ~ChemistryHelper_Structured() {}

  virtual void Advance(const FArrayBox& aqueous_saturation,       int sSat,
                       const FArrayBox& aqueous_pressure,         int sPress,
                       const FArrayBox& porosity,                 int sPhi,
                       const FArrayBox& volume,                   int sVol,
                       FArrayBox&       primary_species_mobile,   int sPrimMob,
                       FArrayBox&       fcnCnt,                   int sFunc,
                       FArrayBox&       aux_data, Real water_density, Real temperature,
                       const Box& box, Real dt, int chem_verbose) = 0;

  virtual void EnforceCondition(FArrayBox& primary_species_mobile,   int sPrimMob,
                                FArrayBox& auxiliary_data, Real water_density, Real temperature,
                                const Box& box, const std::string& condition_name, Real time,int chem_verbose) = 0;

  const std::map<std::string,int>& AuxChemVariablesMap() const {return aux_chem_variables;}
  const std::map<std::string,Real>& AuxChemDefaultsMap() const {return aux_chem_defaults;}

  static std::string BuildPropertyParameterName(const std::string& entity,
                                                const std::string& property,
                                                const std::string& parameter);
  static std::string BuildPropertyParameterName(const std::string& entity,
                                                const std::string& parameter);

  void SetupAuxVariables();

  int NumMobile() const {return Nmobile;}

protected:
  std::map<std::string,int> aux_chem_variables;
  std::map<std::string,Real> aux_chem_defaults;

  bool using_sorption, using_isotherms;
  std::vector<std::string> primarySpeciesNames, mineralNames, surfSiteNames, ionExchangeNames;
  int Nimmobile, Nmobile, Nisotherms, NionExchange, Nminerals, NsorptionSites, NfreeIonSpecies;
};

inline std::string
ChemistryHelper_Structured::BuildPropertyParameterName(const std::string& entity,
                                                       const std::string& property,
                                                       const std::string& parameter)
{
  return entity+"_"+property+"_"+parameter;
}

inline std::string
ChemistryHelper_Structured::BuildPropertyParameterName(const std::string& entity,
                                                       const std::string& parameter)
{
  return entity+"_"+parameter;
}

#endif

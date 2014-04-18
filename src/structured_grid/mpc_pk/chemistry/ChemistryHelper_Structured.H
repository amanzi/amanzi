#ifndef _ChemistryHelper_Structured_H_
#define _ChemistryHelper_Structured_H_

#include <FArrayBox.H>

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
                       const Box& box, Real dt, bool initialize = false) = 0;

  virtual void EnforceCondition(FArrayBox& primary_species_mobile,   int sPrimMob,
                                FArrayBox& auxiliary_data, bool initialize_auxiliary_data,
                                const Box& box, const std::string& condition_name, Real time) = 0;

  const std::map<std::string,int>& AuxChemVariablesMap() const {return aux_chem_variables;}
  static std::string BuildPropertyParameterName(const std::string& entity,
                                                const std::string& property,
                                                const std::string& parameter);
  static std::string BuildPropertyParameterName(const std::string& entity,
                                                const std::string& parameter);

protected:
  std::map<std::string,int> aux_chem_variables;
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

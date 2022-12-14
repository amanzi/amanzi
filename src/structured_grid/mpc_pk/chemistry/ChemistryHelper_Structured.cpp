/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <ChemistryHelper_Structured.H>

#include <Utility.H>

void
ChemistryHelper_Structured::SetupAuxVariables()
{
  // Setup input maps, defaults
  aux_chem_variables.clear();
  aux_chem_defaults.clear();

  for (int i=0; i<primarySpeciesNames.size(); ++i) {
    const std::string label=primarySpeciesNames[i] + "_Activity_Coefficient";
    int n = aux_chem_variables.size();
    aux_chem_variables[label] = n;
    aux_chem_defaults[label] = 1;
  }

  if (NfreeIonSpecies > 0) {
    BL_ASSERT(NfreeIonSpecies == Nmobile);
    for (int i=0; i<Nmobile; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Free_Ion_Guess";
      int n = aux_chem_variables.size();
      aux_chem_variables[label] = n;
      aux_chem_defaults[label] = 1.e-20;
    }
  }

  if (using_sorption) {
    for (int i=0; i<primarySpeciesNames.size(); ++i) {
      const std::string label=primarySpeciesNames[i] + "_Sorbed_Concentration";
      int n = aux_chem_variables.size();
      aux_chem_variables[label] = n;
      aux_chem_defaults[label] = 0;
    }
  }

  if (Nminerals > 0) {
    for (int i=0; i<mineralNames.size(); ++i) {
      const std::string label=mineralNames[i] + "_Volume_Fraction";
      int n = aux_chem_variables.size();
      aux_chem_variables[label] = n;
      aux_chem_defaults[label] = 0;
    }
    for (int i=0; i<mineralNames.size(); ++i) {
      const std::string label=mineralNames[i] + "_Specific_Surface_Area";
      int n = aux_chem_variables.size();
      aux_chem_variables[label] = n;
      aux_chem_defaults[label] = 1;
    }
  }

  if (NsorptionSites > 0) {
    for (int i=0; i<surfSiteNames.size(); ++i) {
      const std::string label=surfSiteNames[i] + "_Surface_Site_Density";
      int n = aux_chem_variables.size();
      aux_chem_variables[label] = n;
      aux_chem_defaults[label] = 1;
    }
  }

  if (NionExchange > 0) {
    int ndigIES = std::log(NionExchange+1);
    for (int i=0; i<NionExchange; ++i) {
      const std::string label = BoxLib::Concatenate("Ion_Exchange_Site_Density_",i,ndigIES);
      int n = aux_chem_variables.size();
      aux_chem_variables[label] = n;
      aux_chem_defaults[label] = 1;
    }
    for (int i=0; i<NionExchange; ++i) {
      const std::string label = BoxLib::Concatenate("Ion_Exchange_Reference_Cation_Concentration_",i,ndigIES);
      int n = aux_chem_variables.size();
      aux_chem_variables[label] = n;
      aux_chem_defaults[label] = 1;
    }
  }

  if (using_isotherms) {
    for (int i=0; i<Nisotherms; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Isotherm_Kd";
      int n = aux_chem_variables.size();
      aux_chem_variables[label] = n;
      aux_chem_defaults[label] = -1;
    }
    for (int i=0; i<Nisotherms; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Isotherm_Freundlich_n";
      int n = aux_chem_variables.size();
      aux_chem_variables[label] = n;
      aux_chem_defaults[label] = 1;
    }
    for (int i=0; i<Nisotherms; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Isotherm_Langmuir_b";
      int n = aux_chem_variables.size();
      aux_chem_variables[label] = n;
      aux_chem_defaults[label] = 1;
    }
  }
}

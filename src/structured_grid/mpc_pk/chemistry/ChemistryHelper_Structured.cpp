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
      aux_chem_variables[label] = aux_chem_variables.size()-1;
      aux_chem_defaults[label] = 0;
    }
  }

  if (Nminerals > 0) {
    for (int i=0; i<mineralNames.size(); ++i) {
      const std::string label=mineralNames[i] + "_Volume_Fraction"; 
      aux_chem_variables[label] = aux_chem_variables.size()-1;
      aux_chem_defaults[label] = 0;
    }
    for (int i=0; i<mineralNames.size(); ++i) {
      const std::string label=mineralNames[i] + "_Specific_Surface_Area"; 
      aux_chem_variables[label] = aux_chem_variables.size()-1;
      aux_chem_defaults[label] = 1;
    }
  }

  if (NsorptionSites > 0) {
    for (int i=0; i<surfSiteNames.size(); ++i) {
      const std::string label=surfSiteNames[i] + "_Surface_Site_Density"; 
      aux_chem_variables[label] = aux_chem_variables.size()-1;
      aux_chem_defaults[label] = 1;
    }
  }

  if (NionExchange > 0) {
    int ndigIES = std::log(NionExchange+1);
    for (int i=0; i<NionExchange; ++i) {
      const std::string label = BoxLib::Concatenate("Ion_Exchange_Site_Density_",i,ndigIES);
      aux_chem_variables[label] = aux_chem_variables.size()-1;
      aux_chem_defaults[label] = 1;
    }
    for (int i=0; i<NionExchange; ++i) {
      const std::string label = BoxLib::Concatenate("Ion_Exchange_Reference_Cation_Concentration_",i,ndigIES);
      aux_chem_variables[label] = aux_chem_variables.size()-1;
      aux_chem_defaults[label] = 1;
    }
  }

  if (using_isotherms) {
    for (int i=0; i<Nisotherms; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Isotherm_Kd"; 
      aux_chem_variables[label] = aux_chem_variables.size()-1;
      aux_chem_defaults[label] = -1;
    }
    for (int i=0; i<Nisotherms; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Isotherm_Freundlich_n"; 
      aux_chem_variables[label] = aux_chem_variables.size()-1;
      aux_chem_defaults[label] = 1;
    }
    for (int i=0; i<Nisotherms; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Isotherm_Langmuir_b"; 
      aux_chem_variables[label] = aux_chem_variables.size()-1;
      aux_chem_defaults[label] = 1;
    }
  }
}

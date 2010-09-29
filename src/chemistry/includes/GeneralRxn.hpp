#ifndef __GeneralRxn_hpp__
#define __GeneralRxn_hpp__

#include <string>
#include <vector>
#include <cmath>

#include "Species.hpp"
#include "Block.hpp"

// Class for general forward/reverse reaction

class GeneralRxn {

 public:
  GeneralRxn();
  GeneralRxn(std::string s);
  GeneralRxn(SpeciesName name, 
             std::vector<SpeciesName>species,
             std::vector<double>stoichiometries,
             std::vector<int>species_ids,
             std::vector<double>forward_stoichiometries,
             std::vector<int>forward_species_ids,
             std::vector<double>backward_stoichiometries,
             std::vector<int>backward_species_ids,
             double kf, double kb);
  ~GeneralRxn();

  // update forward and reverse effective reaction rates
  void update_rates(const std::vector<Species> primarySpecies);
  void addContributionToResidual(std::vector<double> &residual,
                                 double por_den_sat_vol);
  void addContributionToJacobian(Block *J, 
                                 const std::vector<Species> primarySpecies,
                                 double por_den_sat_vol) ;
  void display(void) const;

 protected:

 private:

  int ncomp_; // # components in reaction
  int ncomp_forward_; // # components in forward reaction
  int ncomp_backward_; // # components in backward reaction
  std::vector<SpeciesName> species_names_;
  std::vector<int> species_ids_;       // ids of primary species in rxn
  std::vector<double> stoichiometry_;  // stoich of primary species in rxn
  std::vector<int> forward_species_ids_;       // ids species used in forward rate calc
  std::vector<double> forward_stoichiometry_;  // forward stoich of primary species in rxn
  std::vector<int> backward_species_ids_;      // ids species used in backward rate calc
  std::vector<double> backward_stoichiometry_; // backward stoich of primary species in rxn
  double kf_;     // forward rate constant
  double kb_;     // backward rate constant

  double lnQkf_;  // forward rate storage
  double lnQkb_;  // backward rate storage
  
};

#endif // __GeneralRxn_hpp__

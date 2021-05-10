/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for general forward/reverse reaction
*/

#ifndef AMANZI_CHEMISTRY_GENERAL_RXN_HH_
#define AMANZI_CHEMISTRY_GENERAL_RXN_HH_

#include <string>
#include <vector>

#include "Species.hh"

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations from chemistry
class MatrixBlock;

class GeneralRxn {
 public:
  GeneralRxn();
  GeneralRxn(const Teuchos::ParameterList& plist,
             const std::map<std::string, int>& name_to_id);
  ~GeneralRxn() {};

  // update forward and reverse effective reaction rates
  void update_rates(const std::vector<Species> primary_species);
  void AddContributionToResidual(std::vector<double> *residual,
                                 double por_den_sat_vol);
  void AddContributionToJacobian(MatrixBlock* J,
                                 const std::vector<Species> primary_species,
                                 double por_den_sat_vol);
  void Display(const Teuchos::Ptr<VerboseObject> vo) const;

 private:
  unsigned int ncomp_;  // # components in reaction
  unsigned int ncomp_forward_;  // # components in forward reaction
  unsigned int ncomp_backward_;  // # components in backward reaction
  std::vector<std::string> species_names_;
  std::vector<int> species_ids_;       // ids of primary species in rxn
  std::vector<double> stoichiometry_;  // stoich of primary species in rxn
  std::vector<int> forward_species_ids_;       // ids species used in forward rate calc
  std::vector<double> forward_stoichiometry_;  // forward stoich of primary species in rxn
  std::vector<int> backward_species_ids_;      // ids species used in backward rate calc
  std::vector<double> backward_stoichiometry_;  // backward stoich of primary species in rxn
  double kf_;  // forward rate constant
  double kb_;  // backward rate constant

  double lnQkf_;  // forward rate storage
  double lnQkb_;  // backward rate storage
};

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif

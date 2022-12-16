/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry

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
  GeneralRxn(const Teuchos::ParameterList& plist, const std::map<std::string, int>& name_to_id);
  ~GeneralRxn(){};

  // update forward and reverse effective reaction rates
  void UpdateRates(const std::vector<Species> primary_species);
  void AddContributionToResidual(std::vector<double>* residual, double por_den_sat_vol);
  void AddContributionToJacobian(MatrixBlock* J,
                                 const std::vector<Species> primary_species,
                                 double por_den_sat_vol);
  void Display(const Teuchos::Ptr<VerboseObject> vo) const;

 private:
  int ncomp_;          // # components in reaction
  int ncomp_forward_;  // # components in forward reaction
  int ncomp_backward_; // # components in backward reaction

  std::vector<std::string> species_names_;
  std::vector<int> species_ids_;      // ids of primary species in rxn
  std::vector<double> stoichiometry_; // stoich of primary species in rxn
  std::vector<int> species_ids_f_,
    species_ids_b_; // ids species used in forward and backward reactions
  std::vector<double> stoichiometry_f_,
    stoichiometry_b_; // forward and backward stoichometries of primary species

  // extension to quadratic general kinetic reactions
  bool flag_oc_;
  std::vector<double> orders_oc_;

  double kf_, lnkf_; // forward rate constant
  double kb_, lnkb_; // backward rate constant

  double lnQkf_, lnQkb_; // forward and backward rates
  double lnOcf_, lnOcb_; // fractional order in C
};

} // namespace AmanziChemistry
} // namespace Amanzi

#endif

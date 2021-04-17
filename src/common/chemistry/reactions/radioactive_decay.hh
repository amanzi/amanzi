/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for radioactive decay of aqueous and sorbed components.
  Does not deal with decay of solid phase.
*/
 
#ifndef AMANZI_CHEMISTRY_RADIOACTIVE_DECAY_HH_
#define AMANZI_CHEMISTRY_RADIOACTIVE_DECAY_HH_

#include <string>
#include <vector>

#include "Species.hh"

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations from chemistry
class MatrixBlock;

class RadioactiveDecay {
 public:
  RadioactiveDecay();
  explicit RadioactiveDecay(std::string s) {};

  RadioactiveDecay(const std::vector<std::string>& species_names,
                   const std::vector<int>& species_ids,
                   const std::vector<double>& stoichiometries,
                   const double half_life,
                   const std::string half_life_units);
  ~RadioactiveDecay() {};

  // update forward and reverse effective reaction rates
  void UpdateRate(const std::vector<double>& total,
                  const std::vector<double>& total_sorbed,
                  const double porosity,
                  const double saturation,
                  const double bulk_volume);
  void AddContributionToResidual(std::vector<double> *residual);
  void AddContributionToJacobian(const MatrixBlock& dtotal,
                                 const MatrixBlock& dtotal_sorbed,
                                 const double porosity,
                                 const double saturation,
                                 const double bulk_volume,
                                 MatrixBlock* J);
  void Display(const Teuchos::Ptr<VerboseObject> vo) const;

  int parent_id() const { return species_ids_.at(0); }
  std::string parent_name() const {return species_names_.at(0); }
  double rate() const { return rate_; }
  double rate_constant() const { return rate_constant_; }

 private:
  void ConvertHalfLifeUnits();
  void ConvertHalfLifeToRateConstant();
  std::vector<std::string> species_names_;
  std::vector<int> species_ids_;  // ids of primary species in rxn
  std::vector<double> stoichiometry_;  // stoich of primary species in rxn

  double rate_constant_;  // rate constant [1/sec]
  double half_life_user_;  // user specified units
  std::string half_life_units_;
  double half_life_seconds_;
  double rate_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif

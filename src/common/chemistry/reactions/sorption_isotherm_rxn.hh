/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Base class for sorption isotherm (linear, Langmuir, Freundlich)
  reactions
*/

#ifndef AMANZI_CHEMISTRY_SORPTION_ISOTHERM_RXN_HH_
#define AMANZI_CHEMISTRY_SORPTION_ISOTHERM_RXN_HH_

#include <memory>
#include <string>
#include <vector>

#include<species.hh>
#include<sorption_isotherm.hh>

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations from chemistry
class MatrixBlock;

class SorptionIsothermRxn {
 public:
  SorptionIsothermRxn() {};
  SorptionIsothermRxn(const SpeciesName species_name, 
                      const SpeciesId species_id,
                      std::shared_ptr<SorptionIsotherm> isotherm);
  ~SorptionIsothermRxn() {};

  const std::vector<double>& GetIsothermParameters(void) const;

  void SetIsothermParameters(const std::vector<double>& params);

  std::string IsothermName(void) const {
    return isotherm_->name();
  }

  SorptionIsotherm::SorptionIsothermType IsothermType(void) const {
    return isotherm_->isotherm_type();
  }

  SpeciesId species_id(void) const {
    return species_id_;
  }

  void Update(const std::vector<Species>& primarySpecies);
  // add sorbed concentration to sorbed total
  void AddContributionToTotal(std::vector<double> *total);
  // add sorbed concentration to sorbed total
  void AddContributionToDTotal(const std::vector<Species>& primarySpecies,
                               MatrixBlock* dtotal);
  void Display(const Teuchos::Ptr<VerboseObject> vo) const;

 private:
  SpeciesId species_id_; // ID of primary species
  SpeciesName species_name_; // Name of primary species
  double sorbed_concentration_;
  std::shared_ptr<SorptionIsotherm> isotherm_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif

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

#include "SorptionIsotherm.hh"
#include "Species.hh"

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations from chemistry
class MatrixBlock;

class SorptionIsothermRxn {
 public:
  SorptionIsothermRxn() {};
  SorptionIsothermRxn(const std::string& species_name, 
                      const int species_id,
                      std::shared_ptr<SorptionIsotherm> isotherm);
  ~SorptionIsothermRxn() {};

  const std::vector<double>& GetIsothermParameters() const;

  void SetIsothermParameters(const std::vector<double>& params);

  std::string IsothermName() const { return isotherm_->name(); }

  SorptionIsotherm::SorptionIsothermType IsothermType() const {
    return isotherm_->isotherm_type();
  }

  int species_id() const { return species_id_; }

  void Update(const std::vector<Species>& primarySpecies);
  // add sorbed concentration to sorbed total
  void AddContributionToTotal(std::vector<double> *total);
  // add sorbed concentration to sorbed total
  void AddContributionToDTotal(const std::vector<Species>& primary_species,
                               MatrixBlock* dtotal);
  void Display(const Teuchos::Ptr<VerboseObject> vo) const;

 private:
  int species_id_;  // ID of primary species
  std::string species_name_;  // name of primary species
  double sorbed_concentration_;
  std::shared_ptr<SorptionIsotherm> isotherm_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif

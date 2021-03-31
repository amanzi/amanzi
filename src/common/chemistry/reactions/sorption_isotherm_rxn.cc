/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Base class for sorption isotherm (linear, Langmuir, Freundlich)
  reactions
*/

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "chemistry_exception.hh"
#include "exceptions.hh"
#include "matrix_block.hh"
#include "sorption_isotherm_rxn.hh"

namespace Amanzi {
namespace AmanziChemistry {

SorptionIsothermRxn::SorptionIsothermRxn(const std::string& species_name, 
                                         const int species_id,
                                         std::shared_ptr<SorptionIsotherm> isotherm) 
    : species_name_(species_name), 
      species_id_(species_id), 
      isotherm_(isotherm) {
}


const std::vector<double>& SorptionIsothermRxn::GetIsothermParameters() const {
  return isotherm_->GetParameters();
}


void SorptionIsothermRxn::SetIsothermParameters(const std::vector<double>& params) {
  isotherm_->SetParameters(params);
}


void SorptionIsothermRxn::Update(const std::vector<Species>& primarySpecies) {
  sorbed_concentration_= isotherm_->Evaluate(primarySpecies.at(species_id_));
}


void SorptionIsothermRxn::AddContributionToTotal(std::vector<double> *total) {
  (*total)[species_id_] += sorbed_concentration_;
}


void SorptionIsothermRxn::AddContributionToDTotal(
    const std::vector<Species>& primarySpecies,
    MatrixBlock* dtotal) {
  dtotal->AddValue(species_id_,species_id_,
      isotherm_->EvaluateDerivative(primarySpecies.at(species_id_)));
}


void SorptionIsothermRxn::Display(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << std::setw(12) << species_name_
          << std::setw(15) << isotherm_->name()
          << std::setw(15);
  vo->Write(Teuchos::VERB_HIGH, message.str());
  isotherm_->Display(vo);
}

}  // namespace AmanziChemistry
}  // namespace Amanzi

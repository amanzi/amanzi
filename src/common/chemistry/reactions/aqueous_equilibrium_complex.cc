/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for aqueous equilibrium complexation reaction
*/

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "aqueous_equilibrium_complex.hh"
#include "matrix_block.hh"

namespace Amanzi {
namespace AmanziChemistry {

AqueousEquilibriumComplex::AqueousEquilibriumComplex(
    const SpeciesName name,
    const SpeciesId id,
    const std::vector<SpeciesName>& species,
    const std::vector<double>& stoichiometry,
    const std::vector<SpeciesId>& species_ids,
    const double h2o_stoich,
    const double charge,
    const double mol_wt,
    const double size,
    const double logK)
    : SecondarySpecies(name, id, species, stoichiometry,
                       species_ids, h2o_stoich, charge,
                       mol_wt, size, logK) {
}


void AqueousEquilibriumComplex::Update(const std::vector<Species>& primarySpecies, const Species& water_species) {
  /* This is not the true Q/K for the reaction, but is instead
  **   BC <==> cC + bB
  **   K = a_C^c * a_B^b / a_BC^1
  **   a_BC = a_C^c * a_B^b / K
  **   a_BC = QK = a_C^c * a_B^b / K
  */
  double lnQK = -lnK();
  for (int i = 0; i < ncomp(); i++) {
    lnQK += stoichiometry_.at(i) *
        primarySpecies.at(species_ids_.at(i)).ln_activity();
  }
  // Add the contribution of the water activity
  lnQK += SecondarySpecies::h2o_stoich_ * std::log(water_species.act_coef());
  lnQK_ = lnQK;
  //  molality_ = std::exp(lnQK) / act_coef_;
  update(std::exp(lnQK) / act_coef_);
}


void AqueousEquilibriumComplex::AddContributionToTotal(std::vector<double> *total) {
  for (int i = 0; i < ncomp(); i++) {
    (*total)[species_ids_.at(i)] += stoichiometry_.at(i) * molality();
  }
}


void AqueousEquilibriumComplex::AddContributionToDTotal(
    const std::vector<Species>& primarySpecies,
    MatrixBlock* dtotal) {

  // taking derivative of contribution to residual in row i with respect
  // to species in column j

  // TODO(bandre): because of memory layout, c loops should be for(i) {for(j)}...?

  // column loop
  for (int j = 0; j < ncomp(); j++) {
    int jcomp = species_ids_.at(j);
    double tempd = stoichiometry_.at(j) *
        std::exp(lnQK_ - primarySpecies.at(jcomp).ln_molality()) /
        act_coef_;  // here act_coef is from complex
    // row loop
    for (int i = 0; i < ncomp(); i++) {
      dtotal->AddValue(species_ids_.at(i), jcomp, stoichiometry_.at(i)*tempd);
    }
  }
}


void AqueousEquilibriumComplex::display(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << "    " << name() << " = ";
  // TODO(bandre): uncomment and update test output
  //   if (h2o_stoichiometry_ > 0) {
  //     message << h2o_stoichiometry_ << " " << "H2O" << " + ";
  //   }
  for (unsigned int i = 0; i < species_names_.size(); i++) {
    message << stoichiometry_.at(i) << " " << species_names_.at(i);
    if (i < species_names_.size() - 1) {
      message << " + ";
    }
  }
  if (SecondarySpecies::h2o_stoich_!=0.0) {
  	  message << " + ";
  	  message << std::setprecision(2) << h2o_stoich_ << " " << "H2O";
  }
  message << std::endl;
  message << "        logK = " << logK_ << std::endl;
  message << "        charge = " << charge() << std::endl;
  message << "        mol wt = " << gram_molecular_weight() << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void AqueousEquilibriumComplex::Display(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << "    " << name() << " = "
            << std::fixed << std::setprecision(3);
  // TODO(bandre): uncomment and update test output
  //   if (h2o_stoichiometry_ > 0) {
  //     message << h2o_stoichiometry_ << " " << "H2O" << " + ";
  //   }
  for (unsigned int i = 0; i < species_names_.size(); i++) {
    message << stoichiometry_.at(i) << " " << species_names_.at(i);
    if (i < species_names_.size() - 1) {
      message << " + ";
    }
  }
  if (SecondarySpecies::h2o_stoich_!=0.0) {
    	  message << " + ";
    	  message << std::setprecision(2) << h2o_stoich_ << " " << "H2O";
  }
  message << std::endl;
  message << std::setw(40) << " " << std::fixed
          << std::setprecision(5) << std::setw(10) << logK_
          << std::setprecision(2) << std::setw(8) << charge()
          << std::setprecision(5) << std::setw(10) << gram_molecular_weight()
          << std::setprecision(2) << std::setw(8) << ion_size_parameter()
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void AqueousEquilibriumComplex::DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << std::setw(15) << "Name"
          << std::setw(15) << "Molality"
          << std::setw(15) << "Activity Coeff"
          << std::setw(15) << "Activity"
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void AqueousEquilibriumComplex::DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << std::setw(15) << name()
          << std::scientific << std::setprecision(5)
          << std::setw(15) << molality()
          << std::setw(15) << act_coef()
          << std::setw(15) << activity()
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}

}  // namespace AmanziChemistry
}  // namespace Amanzi

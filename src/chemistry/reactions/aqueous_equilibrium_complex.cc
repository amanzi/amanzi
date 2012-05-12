/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "aqueous_equilibrium_complex.hh"

#include <cmath>

#include <iostream>
#include <iomanip>

#include "matrix_block.hh"

namespace amanzi {
namespace chemistry {

AqueousEquilibriumComplex::AqueousEquilibriumComplex()
    : SecondarySpecies() {
}  // end AqueousEquilibriumComplex() constructor


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
}  // end AqueousEquilibriumComplex() constructor

AqueousEquilibriumComplex::~AqueousEquilibriumComplex() {
}  // end AqueousEquilibriumComplex() destructor

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
}  // end update()

void AqueousEquilibriumComplex::AddContributionToTotal(std::vector<double> *total) {
  for (int i = 0; i < ncomp(); i++) {
    (*total)[species_ids_.at(i)] += stoichiometry_.at(i) * molality();
  }
}  // end addContributionToTotal()

void AqueousEquilibriumComplex::AddContributionToDTotal(
    const std::vector<Species>& primarySpecies,
    MatrixBlock* dtotal) {

  // taking derivative of contribution to residual in row i with respect
  // to species in column j

  // TODO(bandre): because of memory layout, c loops should be for(i){for(j)}...?

  // column loop
  std::cout.precision(15);
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
}  // end addContributionToDTotal()

void AqueousEquilibriumComplex::display(void) const {
  std::cout << "    " << name() << " = ";
  // TODO(bandre): uncomment and update test output
  //   if (h2o_stoichiometry_ > 0) {
  //     std::cout << h2o_stoichiometry_ << " " << "H2O" << " + ";
  //   }
  for (unsigned int i = 0; i < species_names_.size(); i++) {
    std::cout << stoichiometry_.at(i) << " " << species_names_.at(i);
    if (i < species_names_.size() - 1) {
      std::cout << " + ";
    }
  }
  if (SecondarySpecies::h2o_stoich_!=0.0) {
  	  std::cout << " + ";
  	  std::cout << std::setprecision(2) << h2o_stoich_ << " " << "H2O";
  }
  std::cout << std::endl;
  std::cout << "        logK = " << logK_ << std::endl;
  std::cout << "        charge = " << charge() << std::endl;
  std::cout << "        mol wt = " << gram_molecular_weight() << std::endl;
}  // end display()

void AqueousEquilibriumComplex::Display(void) const {
  std::cout << "    " << name() << " = "
            << std::fixed << std::setprecision(3);
  // TODO(bandre): uncomment and update test output
  //   if (h2o_stoichiometry_ > 0) {
  //     std::cout << h2o_stoichiometry_ << " " << "H2O" << " + ";
  //   }
  for (unsigned int i = 0; i < species_names_.size(); i++) {
    std::cout << stoichiometry_.at(i) << " " << species_names_.at(i);
    if (i < species_names_.size() - 1) {
      std::cout << " + ";
    }
  }
  if (SecondarySpecies::h2o_stoich_!=0.0) {
    	  std::cout << " + ";
    	  std::cout << std::setprecision(2) << h2o_stoich_ << " " << "H2O";
  }
  std::cout << std::endl;
  std::cout << std::setw(40) << " " << std::fixed
            << std::setprecision(5) << std::setw(10) << logK_
            << std::setprecision(2) << std::setw(8) << charge()
            << std::setprecision(5) << std::setw(10) << gram_molecular_weight()
            << std::setprecision(2) << std::setw(8) << ion_size_parameter()
            << std::endl;
}  // end Display()

void AqueousEquilibriumComplex::DisplayResultsHeader(void) const {
  std::cout << std::setw(15) << "Name"
            << std::setw(15) << "Molality"
            << std::setw(15) << "Activity Coeff"
            << std::setw(15) << "Activity"
            << std::endl;
}  // end DisplayResultsHeader()

void AqueousEquilibriumComplex::DisplayResults(void) const {
  std::cout << std::setw(15) << name()
            << std::scientific << std::setprecision(5)
            << std::setw(15) << molality()
            << std::setw(15) << act_coef()
            << std::setw(15) << activity()
            << std::endl;
}  // end DisplayResults()

}  // namespace chemistry
}  // namespace amanzi

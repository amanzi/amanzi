/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry

  Class for aqueous equilibrium complexation reaction
*/

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "AqueousEquilibriumComplex.hh"
#include "MatrixBlock.hh"

namespace Amanzi {
namespace AmanziChemistry {

AqueousEquilibriumComplex::AqueousEquilibriumComplex(int id,
                                                     const std::string& name,
                                                     const Teuchos::ParameterList& plist,
                                                     const std::vector<Species>& primary_species)
  : SecondarySpecies(id, name, plist, primary_species)
{}


void
AqueousEquilibriumComplex::Update(const std::vector<Species>& primary_species,
                                  const Species& water_species)
{
  /* This is not the true Q/K for the reaction, but is instead
  **   BC <==> cC + bB
  **   K = a_C^c * a_B^b / a_BC^1
  **   QK = a_BC = a_C^c * a_B^b / K
  */
  double lnQK = -lnK();
  for (int i = 0; i < ncomp(); i++) {
    lnQK += stoichiometry_.at(i) * primary_species.at(species_ids_.at(i)).ln_activity();
  }
  // Add the contribution of the water activity
  lnQK += SecondarySpecies::h2o_stoich_ * std::log(water_species.act_coef());
  lnQK_ = lnQK;
  update(std::exp(lnQK) / act_coef_); // new molality_
}


void
AqueousEquilibriumComplex::AddContributionToTotal(std::vector<double>* total)
{
  for (int i = 0; i < ncomp(); i++) {
    (*total)[species_ids_.at(i)] += stoichiometry_.at(i) * molality();
  }
}


void
AqueousEquilibriumComplex::AddContributionToDTotal(const std::vector<Species>& primary_species,
                                                   MatrixBlock* dtotal)
{
  // derivative of contribution to residual in row i with respect
  // to species in column j

  // column loop
  for (int j = 0; j < ncomp(); j++) {
    int jcomp = species_ids_.at(j);
    double tempd = stoichiometry_.at(j) *
                   std::exp(lnQK_ - primary_species.at(jcomp).ln_molality()) /
                   act_coef_; // here act_coef is from complex
    // row loop
    for (int i = 0; i < ncomp(); i++) {
      dtotal->AddValue(species_ids_.at(i), jcomp, stoichiometry_.at(i) * tempd);
    }
  }
}


void
AqueousEquilibriumComplex::Display(const Teuchos::Ptr<VerboseObject> vo) const
{
  std::stringstream message;
  message << "    " << name() << " = " << std::fixed << std::setprecision(3);
  for (int i = 0; i < species_names_.size(); i++) {
    message << stoichiometry_.at(i) << " " << species_names_.at(i);
    if (i < species_names_.size() - 1) { message << " + "; }
  }
  if (SecondarySpecies::h2o_stoich_ != 0.0) {
    message << " + ";
    message << std::setprecision(2) << h2o_stoich_ << " "
            << "H2O";
  }
  message << std::endl;
  message << std::setw(40) << " " << std::fixed << std::setprecision(5) << std::setw(10) << logK_
          << std::setprecision(2) << std::setw(8) << charge() << std::setprecision(5)
          << std::setw(10) << gram_molecular_weight() << std::setprecision(2) << std::setw(8)
          << ion_size_parameter() << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void
AqueousEquilibriumComplex::DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const
{
  std::stringstream message;
  message << std::setw(15) << "Name" << std::setw(15) << "Molality" << std::setw(15)
          << "Activity Coeff" << std::setw(15) << "Activity" << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}

} // namespace AmanziChemistry
} // namespace Amanzi

/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Base class for secondary species (aqueous equilibrium complexes,
  minerals.
*/
 
#include <iostream>
#include <iomanip>
#include <sstream>

#include "chemistry_exception.hh"
#include "chemistry_utilities.hh"
#include "exceptions.hh"
#include "matrix_block.hh"
#include "secondary_species.hh"

namespace Amanzi {
namespace AmanziChemistry {

namespace acu = Amanzi::AmanziChemistry::utilities;

SecondarySpecies::SecondarySpecies()
    : Species(),
      ncomp_(0),  // # components in reaction
      h2o_stoich_(0.0),
      lnK_(0.0),
      lnQK_(0.0),
      logK_(0.0) {
  species_names_.clear();
  species_ids_.clear();
  stoichiometry_.clear();
  logK_array_.clear();
}


SecondarySpecies::SecondarySpecies(const std::string& in_name,
                                   const int in_id,
                                   const std::vector<std::string>& in_species,
                                   const std::vector<double>& in_stoichiometries,
                                   const std::vector<int>& in_species_ids,
                                   const double in_h2o_stoich,
                                   const double in_charge,
                                   const double in_mol_wt,
                                   const double in_size,
                                   const double in_logK)
    : Species(in_id, in_name, in_charge, in_mol_wt, in_size),
      h2o_stoich_(in_h2o_stoich),
      lnK_(0.0),
      lnQK_(0.0),
      logK_(in_logK) {

  ncomp_ = in_species.size();

  // species names
  for (auto it = in_species.begin(); it != in_species.end(); ++it) {
    species_names_.push_back(*it);
  }

  // species stoichiometries
  for (auto it = in_stoichiometries.begin(); it != in_stoichiometries.end(); ++it) {
    stoichiometry_.push_back(*it);
  }

  // species ids
  for (auto it = in_species_ids.begin(); it != in_species_ids.end(); ++it) {
    species_ids_.push_back(*it);
  }

  lnK_ = acu::log_to_ln(logK());

  //
  // verify the setup
  //

  // must have ncomp > 0, or ncomp > 1?
  if (ncomp() < 1) {
    std::ostringstream error_stream;
    error_stream << "SecondarySpecies::SecondarySpecies(): \n"
                 << "invalid number of components "
                 << "(ncomp < 1), ncomp = " << ncomp() << std::endl;
    Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
  }
  // size of species names, stoichiometries and id arrays must be the same
  if (species_names_.size() != stoichiometry_.size()) {
    std::ostringstream error_stream;
    error_stream << "SecondarySpecies::SecondarySpecies(): \n"
                 << "invalid input data: \n"
                 << "species_names.size() != stoichiometries.size()" << std::endl;
    Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
  }
  if (species_names_.size() != species_ids_.size()) {
    std::ostringstream error_stream;
    error_stream << "SecondarySpecies::SecondarySpecies(): \n"
                 << "invalid input data: \n"
                 << "species_names.size() != species_ids.size()" << std::endl;
    Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
  }
}


/*
**  these functions are only needed if SecondarySpecies equilibrium is added.
*/
void SecondarySpecies::Update(const std::vector<Species>& primary_species) {
  static_cast<void>(primary_species);
}


void SecondarySpecies::AddContributionToTotal(std::vector<double> *total) {
  static_cast<void>(total);
}


void SecondarySpecies::AddContributionToDTotal(
    const std::vector<Species>& primary_species,
    MatrixBlock* dtotal) {
  static_cast<void>(primary_species);
  static_cast<void>(dtotal);
}


/*
**  Display functions
*/
void SecondarySpecies::Display(void) const {
  std::cout << "    " << name() << " = ";
  for (unsigned int i = 0; i < species_names_.size(); i++) {
    std::cout << stoichiometry_[i] << " " << species_names_[i];
    if (i < species_names_.size() - 1) {
      std::cout << " + ";
    }
  }
  std::cout << std::endl;
  std::cout << std::setw(40) << " "
            << std::setw(10) << logK_
            << std::setw(10) << gram_molecular_weight()
            << std::endl;
}

}  // namespace AmanziChemistry
}  // namespace Amanzi

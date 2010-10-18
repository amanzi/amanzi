/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <iomanip>

#include "SecondarySpecies.hpp"

SecondarySpecies::SecondarySpecies() 
    : Species(),
      ncomp_(0), // # components in reaction
      h2o_stoich_(0.0),
      lnK_(0.0),
      lnQK_(0.0),
      logK_(0.0)
{
  species_names_.clear();
  species_ids_.clear();
  stoichiometry_.clear();
  logK_array_.clear();
} // end SecondarySpecies() constructor

SecondarySpecies::SecondarySpecies(const SpeciesName in_name, 
                                   const SpeciesId in_id,
                                   std::vector<SpeciesName> in_species,
                                   std::vector<double> in_stoichiometries,
                                   std::vector<SpeciesId> in_species_ids,
                                   const double in_h2o_stoich, 
                                   const double in_charge,
                                   const double in_mol_wt, 
                                   const double in_size,
                                   const double in_logK)
    : Species(in_id, in_name, in_charge, in_mol_wt, in_size),
      h2o_stoich_(in_h2o_stoich),
      lnK_(0.0),
      lnQK_(0.0),
      logK_(in_logK)
{

  ncomp(static_cast<int>(in_species.size()));

  // species names
  for (std::vector<SpeciesName>::const_iterator i = in_species.begin(); 
       i != in_species.end(); i++) {
    species_names_.push_back(*i);
  } 
  // species stoichiometries
  for (std::vector<double>::const_iterator i = in_stoichiometries.begin();
       i != in_stoichiometries.end(); i++) {
    stoichiometry_.push_back(*i);
  }
  // species ids
  for (std::vector<int>::const_iterator i = in_species_ids.begin();
       i != in_species_ids.end(); i++) {
    species_ids_.push_back(*i);
  }

  logK_ = in_logK;
  lnK_ = log_to_ln(in_logK);

}  // end SecondarySpecies costructor


SecondarySpecies::~SecondarySpecies() 
{
} // end SecondarySpecies() destructor

/*
**
**  these functions are only needed if SecondarySpecies equilibrium is added.
**
*/
void SecondarySpecies::Update_kludge(const std::vector<Species> primary_species) 
{
  static_cast<void>(primary_species);
} // end update()

void SecondarySpecies::AddContributionToTotal(std::vector<double> &total) 
{
  static_cast<void>(total);
} // end addContributionToTotal()

void SecondarySpecies::AddContributionToDTotal(const std::vector<Species> primary_species,
                                               Block *dtotal) 
{
  static_cast<void>(primary_species);
  static_cast<void>(dtotal);
} // end addContributionToDTotal()


/*
**
**  Display functions
**
*/
void SecondarySpecies::Display(void) const
{
  std::cout << "    " << name() << " = ";
  for (int i = 0; i < (int)species_names_.size(); i++) {
    std::cout << stoichiometry_[i] << " " << species_names_[i];
    if (i < (int)species_names_.size() - 1) {
      std::cout << " + ";
    }
  }
  std::cout << std::endl;
  std::cout << std::setw(40) << " " 
            << std::setw(10) << logK_
            << std::setw(10) << gram_molecular_weight()
            << std::endl;
} // end Display()

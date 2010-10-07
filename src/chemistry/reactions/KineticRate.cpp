/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>

#include <iostream>

#include "KineticRate.hpp"
#include "StringTokenizer.hpp"
#include "Verbosity.hpp"

KineticRate::KineticRate(void)
    : verbosity_(kSilent)
{
  reactant_names.clear();
  reactant_stoichiometry.clear();
  reactant_ids.clear();

  product_names.clear();
  product_stoichiometry.clear();
  product_ids.clear();
}  // end KineticRate constructor


KineticRate::~KineticRate(void)
{
}  // end KineticRate destructor


void KineticRate::ParseReaction(const std::string rxn_string)
{
  if (0) {
    std::cout << "Reaction string: " << rxn_string << std::endl;
  }
  std::string rxn_delimiter("=");
  std::string coeff_delimiter(" ");
  StringTokenizer rxn(rxn_string, rxn_delimiter);

  StringTokenizer reactants(rxn.at(0), coeff_delimiter);

  for (std::vector<std::string>::iterator field = reactants.begin(); 
       field != reactants.end(); field++) {
    this->reactant_stoichiometry.push_back(std::atof((*field).c_str()));
    field++;
    this->reactant_names.push_back(*field);
  }

  StringTokenizer products(rxn.at(1), coeff_delimiter);

  for (std::vector<std::string>::iterator field = products.begin(); 
       field != products.end(); field++) {
    this->product_stoichiometry.push_back(std::atof((*field).c_str()));
    field++;
    this->product_names.push_back(*field);
  }
}  // end ParseReaction()


void KineticRate::SetSpeciesIds(const SpeciesArray species, 
                                const std::string species_type,
                                const std::vector<SpeciesName> in_names,
                                const std::vector<double> in_stoichiometry,
                                std::vector<SpeciesId>* out_ids,
                                std::vector<double>* out_stoichiometry)
{
  
  out_ids->clear();
  if (out_stoichiometry != NULL) {
    out_stoichiometry->clear();
  }
  for (unsigned int r = 0; r < in_names.size(); r++) {
    bool species_found = false;
    // check primary species
    SpeciesArray::const_iterator s;
    for (s = species.begin(); s != species.end(); s++) {
      if (in_names.at(r) == (*s).name()) {
        species_found = true;
        out_ids->push_back((*s).identifier());
        if (out_stoichiometry != NULL) {
          out_stoichiometry->push_back(in_stoichiometry.at(r));
        }
        if (verbosity() == kDebugMineralKinetics) {
          std::cout << "    KineticRate::SetSpeciesIds: Found " << species_type 
                    << " species " << (*s).name() << std::endl;
        }
      }
    }
    if (species_found == false && verbosity() == kDebugMineralKinetics) {
      std::cout << "    KineticRate::SetSpeciesIds: Did not find species \'" << in_names.at(r) << "\' in " << species_type << " species list! " << std::endl;
    } 
  }
}  // end SetSpeciesIds()


void KineticRate::DisplayReaction(void) const
{
  std::cout << "  Reaction: " << std::endl;
  std::cout << "    ";
  for (unsigned int species = 0; 
       species < this->reactant_names.size(); species++) {
    std::cout << this->reactant_stoichiometry.at(species) << " " 
              << this->reactant_names.at(species);
    if (species < this->reactant_names.size() - 1) {
      std::cout << " + ";
    }
  }
  std::cout << " = ";
  for (unsigned int species = 0; 
       species < this->product_names.size(); species++) {
    std::cout << this->product_stoichiometry.at(species) << " " 
              << this->product_names.at(species);
    if (species < this->product_names.size() - 1) {
      std::cout << " + ";
    }
  }
  std::cout << std::endl;
}  // end DisplayReaction

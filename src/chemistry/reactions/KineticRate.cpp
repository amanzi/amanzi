/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>

#include <iostream>

#include "KineticRate.hpp"
#include "StringTokenizer.hpp"
#include "Verbosity.hpp"

KineticRate::KineticRate(void)
    : verbosity_(kSilent)
{
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
    this->reactants_stoichiometery.push_back(std::atof((*field).c_str()));
    field++;
    this->reactants_names.push_back(*field);
  }

  StringTokenizer products(rxn.at(1), coeff_delimiter);

  for (std::vector<std::string>::iterator field = products.begin(); 
       field != products.end(); field++) {
    this->products_stoichiometery.push_back(std::atof((*field).c_str()));
    field++;
    this->products_names.push_back(*field);
  }
}  // end ParseReaction

void KineticRate::DisplayReaction(void) const
{
  std::cout << "  Reaction: " << std::endl;
  std::cout << "    ";
  for (unsigned int species = 0; 
       species < this->reactants_names.size(); species++) {
    std::cout << this->reactants_stoichiometery.at(species) << " " 
              << this->reactants_names.at(species);
    if (species < this->reactants_names.size() - 1) {
      std::cout << " + ";
    }
  }
  std::cout << " = ";
  for (unsigned int species = 0; 
       species < this->products_names.size(); species++) {
    std::cout << this->products_stoichiometery.at(species) << " " 
              << this->products_names.at(species);
    if (species < this->products_names.size() - 1) {
      std::cout << " + ";
    }
  }
  std::cout << std::endl;
}  // end DisplayReaction

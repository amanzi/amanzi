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
  if (verbosity() == kDebugMineralKinetics) {
    std::cout << "Reaction string: " << rxn_string << std::endl;
  }
  std::string rxn_delimiter("=");
  std::string coeff_delimiter(" ");
  StringTokenizer rxn(rxn_string, rxn_delimiter);

  StringTokenizer reactants(rxn.at(0), coeff_delimiter);
  std::vector<SpeciesName> reactants_names;
  std::vector<double> reactants_stoichiometery;

  for (std::vector<std::string>::iterator field = reactants.begin(); 
       field != reactants.end(); field++) {
    reactants_stoichiometery.push_back(std::atof((*field).c_str()));
    field++;
    reactants_names.push_back(*field);
  }

  StringTokenizer products(rxn.at(1), coeff_delimiter);
  std::vector<SpeciesName> products_names;
  std::vector<double> products_stoichiometery;

  for (std::vector<std::string>::iterator field = products.begin(); 
       field != products.end(); field++) {
    products_stoichiometery.push_back(std::atof((*field).c_str()));
    field++;
    products_names.push_back(*field);
  }

  if (verbosity() == kDebugMineralKinetics) {
    for (unsigned int species = 0; species < reactants_names.size(); species++) {
      std::cout << reactants_stoichiometery.at(species) << " " 
                << reactants_names.at(species) << " ";
    }
    std::cout << " = ";
    for (unsigned int species = 0; species < products_names.size(); species++) {
      std::cout << products_stoichiometery.at(species) << " " 
                << products_names.at(species) << " ";
    }
  }  
}  // end ParseReaction

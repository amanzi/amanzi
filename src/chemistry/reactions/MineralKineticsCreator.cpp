/* -*-  mode: c++; c-default-style: "google-c-style"; indent-tabs-mode: nil -*- */
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "MineralKineticsCreator.hpp"
#include "Species.hpp"
#include "StringTokenizer.hpp"
#include "Verbosity.hpp"

MineralKineticsCreator::MineralKineticsCreator(void)
  : 
  verbosity_(kSilent)
{

} // end MineralKineticsCreator constructor

MineralKineticsCreator::~MineralKineticsCreator(void)
{

} // end MineralKineticsCreator destructor

void MineralKineticsCreator::readFile(const std::string file_name)
{
  if (verbosity() == kDebugMineralKinetics) {
    std::cout << "MineralKinetics::readFile()...." << std::endl;
  }

  std::ifstream input(file_name.c_str());
  if (!input) {
    // should be some type of helpful error message and graceful exit here....
  }

  int count = 0;
  while (!input.eof() && count < 100) {
    count++;
    std::string line;
    getline(input, line);
    StringTokenizer st;
    if (line[0] != '#' && line[0] != '\0' && line[0] != ' ') {
      std::string delimiter(";");
      st.tokenize(line, delimiter);
      this->ParseReaction(st.at(0));
      this->ParseRate(st);
    } else {
      if (0) {
        std::cout << "Comment: " << line << std::endl;
      }
    }
  }

} // end readFile()

void MineralKineticsCreator::ParseReaction(const std::string rxn_string)
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
  
} // end ParseReaction

void MineralKineticsCreator::ParseRate(StringTokenizer rate)
{
  std::string space(" ");
  StringTokenizer rate_name(rate.at(1), space); // strip out spaces
  //std::cout << "rate_name[0] = \'" << rate_name.at(0) << "\'" << std::endl;

  if (!(rate_name.at(0).compare("TST"))) {
    if (verbosity() == kDebugMineralKinetics) {
      std::cout << "Rate type: \'" << rate_name.at(0) << "\' = " << "\'TST\'" << std::endl;
    }
    this->ParseTstParameters(rate);
  } else {
    std::cout << "Unknown Rate type: \'" << rate_name.at(0) << "\'" << std::endl;
    // some sort of gracefull exit here....
  }


} // end ParseRate()

void MineralKineticsCreator::ParseTstParameters(StringTokenizer rate)
{
  double area = 0.0;
  double pK = 0.0;
  double k = 0.0;
  double sat_state_exponent = 0.0;
  std::vector<SpeciesName> modifying_species_names;
  std::vector<double> modifying_exponents;
  std::vector<int> modifying_species_ids;

  std::string space(" ");
  StringTokenizer st;

  std::vector<std::string>::iterator field = rate.begin();
  field++; // field[0] is the reaction string, handle it elsewhere
  field++; // field[1] is the rate name, handled elsewhere

  for (; field != rate.end(); field++) {
    st.tokenize(*field, space);
    if (!(st.at(0).compare("area"))) {
      area = std::atof(st.at(1).c_str());
    } else if (!(st.at(0).compare("pK"))) {
      pK = std::atof(st.at(1).c_str());
    } else if (!(st.at(0).compare("k"))) {
      k = std::atof(st.at(1).c_str());
    } else if (!(st.at(0).compare("n"))) {
      sat_state_exponent = std::atof(st.at(1).c_str());
    } else {
      // assume we are dealing with the list of rate modifying speces
      for (unsigned int modifier = 0; modifier < st.size(); modifier++) {
        modifying_species_names.push_back(st.at(modifier));
        modifier++; // increment to get the exponent of this modifier
        modifying_exponents.push_back(std::atof(st.at(modifier).c_str()));
      }      
    }
    if (verbosity() == kDebugMineralKinetics) {
      std::cout << "  Field: " << *field << std::endl;
    }
  }
  if (verbosity() == kDebugMineralKinetics) {
    std::cout << "area = " << area << std::endl;
    std::cout << "pK = " << pK << std::endl;
    std::cout << "k = " << k << std::endl;
    std::cout << "n = " << sat_state_exponent << std::endl;
    std::cout << "rate modifiers: " << std::endl;
    for (unsigned int mod = 0; mod < modifying_species_names.size(); mod++) {
      std::cout << "species: " << modifying_species_names.at(mod) << " ";
      std::cout << "exponent: " << modifying_exponents.at(mod) << std::endl;
    }
  }
} // end ParseRateTST()

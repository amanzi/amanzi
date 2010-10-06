/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "MineralKineticsFactory.hpp"
#include "KineticRateTST.hpp"
#include "KineticRate.hpp"
#include "Species.hpp"
#include "StringTokenizer.hpp"
#include "Verbosity.hpp"

const std::string MineralKineticsFactory::kTST = "TST";

MineralKineticsFactory::MineralKineticsFactory(void)
    : verbosity_(kSilent)
{
}  // end MineralKineticsFactory constructor

MineralKineticsFactory::~MineralKineticsFactory(void)
{
}  // end MineralKineticsFactory destructor

std::vector<KineticRate*> MineralKineticsFactory::Create(const std::string file_name,
                                                         const std::vector<Species> primary_species)
{
  // the input file will generally contain multiple kinetic rates. so
  // this method should return an array of rate pointes....
  this->rates.clear();

  ReadFile(file_name, primary_species);

  return this->rates;
}  // end Create()


/*******************************************************************************
**
**  Mineral kinetics file format
**
**  all rate information is contained on a single semicolon line.
**
**  any line which begins with a # or a space is considered a
**  comment. blank lines are ignored
**
**  Field 0 : coeff MineralName = coeff1 Primary1 coeff2 Primary2 .....
**
**  Field 1 : rate_name ("TST", ....)
**
**  remaining fields are determined by the rate type
**
*******************************************************************************/
void MineralKineticsFactory::ReadFile(const std::string file_name,
                                      const std::vector<Species> primary_species)
{
  if (verbosity() == kDebugMineralKinetics) {
    std::cout << "MineralKineticsFactory::ReadFile()...." << std::endl;
  }

  std::ifstream input(file_name.c_str());
  if (!input) {
    // should be some type of helpful error message and graceful exit here....
  }

  std::string delimiter(";");
  int count = 0;
  while (!input.eof() && count < 100) {
    count++;
    std::string line;
    getline(input, line);
    StringTokenizer reaction_data;
    if (line[0] != '#' && line[0] != '\0' && line[0] != ' ') {
      // not a comment line or white space
      reaction_data.tokenize(line, delimiter);
      std::string reaction = reaction_data.at(0);
      std::string rate_type = reaction_data.at(1);
      reaction_data.erase(reaction_data.begin());  // erase reaction string
      reaction_data.erase(reaction_data.begin());  // erase reaction type string
      KineticRate* kinetic_rate = NULL;
      kinetic_rate = CreateRate(rate_type);
      if (kinetic_rate != NULL) {
        kinetic_rate->verbosity(verbosity());
        kinetic_rate->Setup(reaction, reaction_data, primary_species);
        this->rates.push_back(kinetic_rate);
      }
    } else {
      if (0) {
        std::cout << "Comment: " << line << std::endl;
      }
    }
  }
}  // end ReadFile()

KineticRate* MineralKineticsFactory::CreateRate(std::string rate_type)
{
  KineticRate* kinetic_rate = NULL;

  std::string space(" ");
  StringTokenizer rate_name(rate_type, space); // strip out spaces

  //std::cout << "rate_name[0] = \'" << rate_name.at(0) << "\'" << std::endl;

  if (!(rate_name.at(0).compare(this->kTST))) {
    kinetic_rate = new KineticRateTST();

  } else {
    std::cout << "Unknown Rate type: \'" << rate_name.at(0) << "\'" << std::endl;
    // some sort of gracefull exit here....
  }

  if (kinetic_rate == NULL) {
    // failed to create a rate object, error message and graceful exit here....
    std::cout << "MineralKineticsFactory::ParseRate() : could not create rate type: "
              << rate_name[0] << std::endl;
  }
  return kinetic_rate;
}  // end CreateRate()



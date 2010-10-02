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

std::vector<KineticRate*> MineralKineticsFactory::Create(std::string file_name)
{
  // the input file will generally contain multiple kinetic rates. so
  // this method should return an array of rate pointes....
  this->rates.clear();

  ReadFile(file_name);

  return this->rates;
}  // end Create()

void MineralKineticsFactory::ReadFile(const std::string file_name)
{
  if (verbosity() == kDebugMineralKinetics) {
    std::cout << "MineralKineticsFactory::ReadFile()...." << std::endl;
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
      // not a comment line or white space
      std::string delimiter(";");
      st.tokenize(line, delimiter);

      KineticRate* kinetic_rate = NULL;
      ParseRate(st, &kinetic_rate);

      kinetic_rate->verbosity(kDebugMineralKinetics);
      kinetic_rate->ParseReaction(st.at(0));
      kinetic_rate->ParseParameters(st);

      this->rates.push_back(kinetic_rate);
    } else {
      if (0) {
        std::cout << "Comment: " << line << std::endl;
      }
    }
  }
}  // end readFile()

void MineralKineticsFactory::ParseRate(StringTokenizer rate, KineticRate** kinetic_rate)
{
  std::string space(" ");
  StringTokenizer rate_name(rate.at(1), space); // strip out spaces
  //std::cout << "rate_name[0] = \'" << rate_name.at(0) << "\'" << std::endl;

  if (!(rate_name.at(0).compare(this->kTST))) {
    if (verbosity() == kDebugMineralKinetics) {
      std::cout << "Rate type: \'" << rate_name.at(0) << "\' = " 
                << "\'" << this->kTST << "\'" << std::endl;
    }

    *kinetic_rate = new KineticRateTST();

  } else {
    std::cout << "Unknown Rate type: \'" << rate_name.at(0) << "\'" << std::endl;
    // some sort of gracefull exit here....
  }

  if (*kinetic_rate == NULL) {
    // failed to create a rate object, error message and graceful exit here....
    std::cout << "MineralKineticsFactory::ParseRate() : could not create rate type: "
              << rate_name[0] << std::endl;
  }

}  // end ParseRate()



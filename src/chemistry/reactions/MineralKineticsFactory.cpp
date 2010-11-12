/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "Mineral.hpp"
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


KineticRate* MineralKineticsFactory::Create(const std::string& rate_type, 
                                            const StringTokenizer& rate_data,
                                            const Mineral& mineral,
                                            const SpeciesArray& primary_species)
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

  if (kinetic_rate != NULL) {
    kinetic_rate->verbosity(verbosity());
    kinetic_rate->Setup(dynamic_cast<const SecondarySpecies&>(mineral), rate_data, primary_species);
  } else {
    // failed to create a rate object, error message and graceful exit here....
    std::cout << "MineralKineticsFactory::ParseRate() : could not create rate type: "
              << rate_name[0] << std::endl;
  }

  return kinetic_rate;
}  // end Create()


SpeciesId MineralKineticsFactory::VerifyMineralName(const std::string mineral_name,
                                                    const std::vector<Mineral>& minerals) const
{
  bool mineral_found = false;
  int mineral_id = -1;
  for (std::vector<Mineral>::const_iterator m = minerals.begin();
       m != minerals.end(); m++) {
    if (m->name() == mineral_name) {
      mineral_found = true;
      mineral_id = m->identifier();
    }
  }
  if (!mineral_found) {
    // print helpful message and exit gracefully
    std::cout << "MineralKineticsFactory::VerifyMineralName(): Did not find mineral \'"
              << mineral_name << "\' in the mineral list." << std::endl;
  }
  return mineral_id;
}  // end VerifyMineralName()

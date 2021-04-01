/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Factory class for building a mineral kinetic rate object
*/
 
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "chemistry_exception.hh"
#include "exceptions.hh"
#include "kinetic_rate_tst.hh"
#include "kinetic_rate.hh"
#include "mineral.hh"
#include "mineral_kinetics_factory.hh"
#include "species.hh"
#include "string_tokenizer.hh"

namespace Amanzi {
namespace AmanziChemistry {

const std::string MineralKineticsFactory::kTST = "TST";

MineralKineticsFactory::MineralKineticsFactory() : debug_(false) {
} 


KineticRate* MineralKineticsFactory::Create(const std::string& rate_type,
                                            const StringTokenizer& rate_data,
                                            const Mineral& mineral,
                                            const SpeciesArray& primary_species) {
  KineticRate* kinetic_rate = NULL;

  std::string space(" ");
  StringTokenizer rate_name(rate_type, space);  // strip out spaces

  if (!(rate_name.at(0).compare(this->kTST))) {
    kinetic_rate = new KineticRateTST();
  } else {
    std::ostringstream error_stream;
    error_stream << "MineralKineticsFactory::Create(): \n";
    error_stream << "Unknown kinetic rate name: " << rate_name.at(0)
                 << "\n       valid names: " << kTST << "\n";
    Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
  }

  if (kinetic_rate != NULL) {
    kinetic_rate->set_debug(debug());
    // TODO(bandre): get rid of the dynamic cast...?
    kinetic_rate->Setup(dynamic_cast<const SecondarySpecies&>(mineral),
                        rate_data, primary_species);
  } else {
    // new failed to create a rate object, for some reason.... Do we
    // really need this check?
    std::ostringstream error_stream;
    error_stream << "MineralKineticsFactory::Create(): \n";
    error_stream << "could not create rate type: " << rate_name.at(0)
                 << "\n       new failed...?" << std::endl;
    Exceptions::amanzi_throw(ChemistryUnrecoverableError(error_stream.str()));
  }

  return kinetic_rate;
}


int MineralKineticsFactory::VerifyMineralName(const std::string& mineral_name,
                                              const std::vector<Mineral>& minerals) const {
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
    std::ostringstream error_stream;
    error_stream << "MineralKineticsFactory::VerifyMineralName(): \n";
    error_stream << "Did not find mineral: \'" << mineral_name << "\'\n"
                 << "       in the mineral list. " << std::endl;
    Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
  }

  return mineral_id;
}

}  // namespace AmanziChemistry
}  // namespace Amanzi

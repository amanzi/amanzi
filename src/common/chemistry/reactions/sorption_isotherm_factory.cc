/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/

#include <cstdlib>
#include <sstream>
#include <string>

#include "chemistry_exception.hh"
#include "exceptions.hh"
#include "sorption_isotherm_factory.hh"
#include "sorption_isotherm.hh"
#include "sorption_isotherm_linear.hh"
#include "sorption_isotherm_langmuir.hh"
#include "sorption_isotherm_freundlich.hh"
#include "string_tokenizer.hh"

namespace Amanzi {
namespace AmanziChemistry {

const std::string SorptionIsothermFactory::linear = "linear";
const std::string SorptionIsothermFactory::langmuir = "langmuir";
const std::string SorptionIsothermFactory::freundlich = "freundlich";

SorptionIsotherm* SorptionIsothermFactory::Create( 
    const std::string& isotherm_type,
    const StringTokenizer parameters) {
  SorptionIsotherm* sorption_isotherm = NULL;

  if (isotherm_type == linear) {
    SorptionIsothermLinear *linear_isotherm =
      new SorptionIsothermLinear(std::atof(parameters[0].c_str()));
    sorption_isotherm = linear_isotherm;
  } else if (isotherm_type == langmuir) {
    // require two parameters
    if (parameters.size() != 2) {
      std::ostringstream error_stream;
      error_stream << "SorptionIsothermFactory::Create(): \n"
                   << "  Langmuir isotherm requires exactly two parameters, received "
                   << parameters.size() << ".\n"
                   << "    param_1 == Kd, param_2 == b  .\n"; 
      Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
    }
    SorptionIsothermLangmuir *langmuir_isotherm =
        new SorptionIsothermLangmuir(std::atof(parameters.at(0).c_str()),
                                     std::atof(parameters.at(1).c_str()));
    sorption_isotherm = langmuir_isotherm;
  } else if (isotherm_type == freundlich) {
    // require two parameters
    if (parameters.size() != 2) {
      std::ostringstream error_stream;
      error_stream << "SorptionIsothermFactory::Create(): \n"
                   << "  Freundlich isotherm: C_sorb = Kd * C^n"
                   << "  Freundlich isotherm requires exactly two parameters, received "
                   << parameters.size() << ".\n"
                   << "    param_1 == Kd, param_2 == n  .\n"; 
      Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
    }
    SorptionIsothermFreundlich *freundlich_isotherm =
        new SorptionIsothermFreundlich(std::atof(parameters.at(0).c_str()),
                                       std::atof(parameters.at(1).c_str()));
    sorption_isotherm = freundlich_isotherm;
  } else {
    // default type, error...!
    std::ostringstream error_stream;
    error_stream << "SorptionIsothermFactory::Create(): \n"
                 << "Unknown sorption isotherm type: " << isotherm_type << "\n"
                 << "       valid names: " << linear << "\n"
                 << "                    " << langmuir << "\n"
                 << "                    " << freundlich << "\n";
    Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
  }

  if (sorption_isotherm == NULL) {
    // something went wrong, should throw an exception and exit gracefully....
    std::ostringstream error_stream;
    error_stream << "SorptionIsothermFactory::Create(): \n"
                 << "SorptionIsotherm was not created for some reason....\n";
    Exceptions::amanzi_throw(ChemistryException(error_stream.str()));
  } else {
    // finish any additional setup
    // GEH - none for now
  }

  return sorption_isotherm;
}


SpeciesId SorptionIsothermFactory::VerifySpeciesName(
    const SpeciesName species_name,
    const std::vector<Species>& species) const {
  int species_id = -1;
  for (std::vector<Species>::const_iterator s = species.begin();
       s != species.end(); s++) {
    if (s->name() == species_name) {
      species_id = s->identifier();
      break;
    }
  }
  if (species_id < 0) {
    // print helpful message and exit gracefully
    std::ostringstream error_stream;
    error_stream << "SorptionIsothermFactory::VerifySpeciesName(): \n";
    error_stream << "Did not find species: \'" << species_name << "\'\n"
                 << "       in the primary species list. " << std::endl;
    Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
  }

  return species_id;
}

}  // namespace AmanziChemistry
}  // namespace Amanzi

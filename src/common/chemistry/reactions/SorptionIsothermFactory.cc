/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry

*/

#include <cstdlib>
#include <sstream>
#include <string>

#include "exceptions.hh"
#include "errors.hh"

#include "SorptionIsotherm.hh"
#include "SorptionIsothermFactory.hh"
#include "SorptionIsothermLinear.hh"
#include "SorptionIsothermLangmuir.hh"
#include "SorptionIsothermFreundlich.hh"

namespace Amanzi {
namespace AmanziChemistry {

const std::string SorptionIsothermFactory::linear = "linear";
const std::string SorptionIsothermFactory::langmuir = "langmuir";
const std::string SorptionIsothermFactory::freundlich = "freundlich";

std::shared_ptr<SorptionIsotherm>
SorptionIsothermFactory::Create(const Teuchos::ParameterList& plist)
{
  std::string isotherm_type = plist.get<std::string>("model");
  std::vector<double> parameters = plist.get<Teuchos::Array<double>>("parameters").toVector();

  std::shared_ptr<SorptionIsotherm> sorption_isotherm = nullptr;

  if (isotherm_type == linear) {
    sorption_isotherm = std::make_shared<SorptionIsothermLinear>(parameters[0]);
  } else if (isotherm_type == langmuir) {
    // require two parameters
    if (parameters.size() != 2) {
      std::ostringstream oss;
      oss << "SorptionIsothermFactory::Create(): \n"
          << "  Langmuir isotherm requires exactly two parameters, received " << parameters.size()
          << ".\n"
          << "    param_1 == Kd, param_2 == b  .\n";
      Exceptions::amanzi_throw(Errors::Message(oss.str()));
    }
    sorption_isotherm = std::make_shared<SorptionIsothermLangmuir>(parameters[0], parameters[1]);
  } else if (isotherm_type == freundlich) {
    // require two parameters
    if (parameters.size() != 2) {
      std::ostringstream error_stream;
      error_stream << "SorptionIsothermFactory::Create(): \n"
                   << "  Freundlich isotherm: C_sorb = Kd * C^n"
                   << "  Freundlich isotherm requires exactly two parameters, received "
                   << parameters.size() << ".\n"
                   << "    param_1 == Kd, param_2 == n  .\n";
      Exceptions::amanzi_throw(Errors::Message(error_stream.str()));
    }
    sorption_isotherm = std::make_shared<SorptionIsothermFreundlich>(parameters[0], parameters[1]);
  } else {
    // default type, error...!
    std::ostringstream error_stream;
    error_stream << "SorptionIsothermFactory::Create(): \n"
                 << "Unknown sorption isotherm type: " << isotherm_type << "\n"
                 << "       valid names: " << linear << "\n"
                 << "                    " << langmuir << "\n"
                 << "                    " << freundlich << "\n";
    Exceptions::amanzi_throw(Errors::Message(error_stream.str()));
  }

  // finish any additional setup
  // GEH - none for now

  return sorption_isotherm;
}


int
SorptionIsothermFactory::VerifySpeciesName(const std::string& species_name,
                                           const std::vector<Species>& species) const
{
  int species_id = -1;
  for (auto it = species.begin(); it != species.end(); ++it) {
    if (it->name() == species_name) {
      species_id = it->identifier();
      break;
    }
  }

  if (species_id < 0) {
    // print helpful message and exit gracefully
    std::ostringstream error_stream;
    error_stream << "SorptionIsothermFactory::VerifySpeciesName(): \n";
    error_stream << "Did not find species: \'" << species_name << "\'\n"
                 << "       in the primary species list. " << std::endl;
    Exceptions::amanzi_throw(Errors::Message(error_stream.str()));
  }

  return species_id;
}

} // namespace AmanziChemistry
} // namespace Amanzi

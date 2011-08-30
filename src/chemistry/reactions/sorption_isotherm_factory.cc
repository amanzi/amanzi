/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "sorption_isotherm_factory.hh"

#include <sstream>
#include <string>

#include "sorption_isotherm.hh"
#include "sorption_isotherm_linear.hh"
#include "sorption_isotherm_langmuir.hh"
#include "sorption_isotherm_freundlich.hh"
#include "chemistry_exception.hh"
#include "verbosity.hh"
#include "exceptions.hh"

namespace amanzi {
namespace chemistry {

const std::string SorptionIsothermFactory::linear = "linear";
const std::string SorptionIsothermFactory::langmuir = "langmuir";
const std::string SorptionIsothermFactory::freundlich = "freundlich";

SorptionIsothermFactory::SorptionIsothermFactory(void) : verbosity_(kSilent) {
}  // end ActivityModelFactory constructor

SorptionIsothermFactory::~SorptionIsothermFactory(void) {
}  // end ActivityModelFactory destructor

SorptionIsotherm* SorptionIsothermFactory::Create( 
    const std::string& isotherm_type) {
  SorptionIsotherm* sorption_isotherm = NULL;

  if (isotherm_type == linear) {
    sorption_isotherm = new SorptionIsothermLinear();
  } else if (isotherm_type == langmuir) {
    sorption_isotherm = new SorptionIsothermLangmuir();
  } else if (isotherm_type == freundlich) {
    sorption_isotherm = new SorptionIsothermFreundlich();
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
    error_stream << "ActivityModelFactory::Create(): \n"
                 << "Activity model was not created for some reason....\n";
    Exceptions::amanzi_throw(ChemistryException(error_stream.str()));
  } else {
    // finish any additional setup
    // GEH - none for now
  }
  return sorption_isotherm;
}  // end Create()


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
}  // end VerifyMineralName()

}  // namespace chemistry
}  // namespace amanzi

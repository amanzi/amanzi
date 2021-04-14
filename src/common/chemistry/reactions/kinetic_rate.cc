/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Abstract base class for all kinetic rates
*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "kinetic_rate.hh"

namespace Amanzi {
namespace AmanziChemistry {

KineticRate::KineticRate()
  : debug_(false),
    name_("KineticRate"),
    identifier_(0),
    reaction_rate_(0.0)
{
  reactant_names.clear();
  reactant_stoichiometry.clear();
  reactant_ids.clear();
}


void KineticRate::SetSpeciesIds(const SpeciesArray& species,
                                const std::string& species_type,
                                const std::vector<std::string>& in_names,
                                const std::vector<double>& in_stoichiometry,
                                std::vector<int>* out_ids,
                                std::vector<double>* out_stoichiometry) {
  /*
  ** loop through a list of input names. Compare the names to the
  ** names of the species. When a match is found, set the output id to
  ** match the identifier of the speciess and possibly set the
  ** corresponding _location_ in the output stoichiometry list to the
  ** correct value.
  */
  out_ids->clear();
  if (out_stoichiometry != NULL) {
    out_stoichiometry->clear();
    out_stoichiometry->resize(species.size(), 0.0);
  }
  for (int current = 0; current < in_names.size(); current++) {
    bool species_found = false;
    // check primary species
    SpeciesArray::const_iterator s;
    for (s = species.begin(); s != species.end(); s++) {
      if (in_names.at(current) == (*s).name()) {
        species_found = true;
        out_ids->push_back((*s).identifier());
        if (out_stoichiometry != NULL) {
          // geh          (*out_stoichiometry)[current] = in_stoichiometry.at(current);
          (*out_stoichiometry)[(*s).identifier()] = in_stoichiometry.at(current);
        }
        if (debug()) {
          std::cout << "    KineticRate::SetSpeciesIds: Found " << species_type
                    << " species " << (*s).name() << std::endl;
        }
      }
    }
    if (species_found == false && debug()) {
      // TODO(bandre): is this actually a runtime error?
      std::cout << "    KineticRate::SetSpeciesIds: Did not find species \'"
                << in_names.at(current) << "\' in " << species_type
                << " species list! " << std::endl;
    }
  }
}


void KineticRate::DisplayReaction(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << "    Reaction: " << std::endl;
  message << "      ";

  message << name();
  message << " = ";
  for (int species = 0; species < this->reactant_names.size(); species++) {
    message << std::setprecision(2) << this->reactant_stoichiometry.at(species) << " "
            << this->reactant_names.at(species);
    if (species < this->reactant_names.size() - 1) {
      message << " + ";
    }
  }
  message << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}

}  // namespace AmanziChemistry
}  // namespace Amanzi

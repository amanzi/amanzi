/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Base class for secondary species (aqueous equilibrium complexes,
  minerals.
*/
 
#include <iostream>
#include <iomanip>
#include <sstream>

#include "errors.hh"

#include "ChemistryUtilities.hh"
#include "exceptions.hh"
#include "MatrixBlock.hh"
#include "SecondarySpecies.hh"

namespace Amanzi {
namespace AmanziChemistry {

namespace acu = Amanzi::AmanziChemistry::utilities;

/* *******************************************************************
* Trivial constructor
******************************************************************* */
SecondarySpecies::SecondarySpecies()
  : Species(),
    ncomp_(0),  // # components in reaction
    h2o_stoich_(0.0),
    lnK_(0.0),
    lnQK_(0.0),
    logK_(0.0)
{
  species_names_.clear();
  species_ids_.clear();
  stoichiometry_.clear();
}


/* *******************************************************************
* Actual contructor
******************************************************************* */
SecondarySpecies::SecondarySpecies(int id, const std::string& name,
                                   const Teuchos::ParameterList& plist,
                                   const std::vector<Species>& primary_species)
  : Species(id, name, plist),
    lnK_(0.0),
    lnQK_(0.0)
{
  // two ways to initialize the reaction equilibrium constant
  if (plist.isSublist("equilibrium constant")) {
    auto x = plist.sublist("equilibrium constant").get<Teuchos::Array<double> >("T").toVector();
    auto y = plist.sublist("equilibrium constant").get<Teuchos::Array<double> >("Keq").toVector();
    func_ = Teuchos::rcp(new FunctionTabular(x, y, 0));

    double T = plist.get<double>("temperature");
    logK_ = (*func_)({T});
  } else {
    logK_ = plist.get<double>("equilibrium constant");
  }
  std::string reaction = plist.get<std::string>("reaction");

  ParseReaction_(reaction, &species_names_, &species_ids_, &stoichiometry_, &h2o_stoich_, primary_species);

  ncomp_ = species_names_.size();
  lnK_ = acu::log_to_ln(logK());

  // verify the setup
  // must have ncomp > 0, or ncomp > 1?
  if (ncomp() < 1 || 
      species_names_.size() != stoichiometry_.size() || 
      species_names_.size() != species_ids_.size()) {
    std::ostringstream oss;
    oss << "Invalid data for secondary species: ncomp = " << ncomp() << std::endl
        << ", species_names.size != stoichiometries.size != species_ids.size" << std::endl;
    Exceptions::amanzi_throw(Errors::Message(oss.str()));
  }
}


/* *******************************************************************
* Recalculates equilibrium constant
******************************************************************* */
void SecondarySpecies::UpdateTemperatureDependentCoefs(double T)
{
  if (func_.get() != nullptr) {
    logK_ = (*func_)({T});
  }
std::cout << T << " " << logK_ << std::endl; 
}


/* *******************************************************************
* Parses a reaction whose products are all primary species
******************************************************************* */
void SecondarySpecies::ParseReaction_(const std::string& reaction,
                                      std::vector<std::string>* species,
                                      std::vector<int>* species_ids,
                                      std::vector<double>* stoichiometries,
                                      double* h2o_stoich,
                                      const std::vector<Species>& primary_species)
{
  double coeff;
  std::string primary_name;

  species->clear();
  species_ids->clear();
  stoichiometries->clear();
  *h2o_stoich = 0.0;

  std::istringstream iss(reaction);
  while (iss >> coeff || !iss.eof()) {
    iss >> primary_name;

    if (primary_name == "H2O") {
      *h2o_stoich = coeff;
    } else {
      int id(-1);
      for (auto it = primary_species.begin(); it != primary_species.end(); ++it) {
        if (it->name() == primary_name) {
          id = it->identifier();
          break;
        }
      }

      if (id < 0) {
        std::stringstream msg;
        msg << "Reaction primary species \'" << primary_name 
            << "\' was not found in the primary species list\n";
        Exceptions::amanzi_throw(Errors::Message(msg.str()));
      }

      species->push_back(primary_name);
      species_ids->push_back(id);
      stoichiometries->push_back(coeff);
    }
  }
}

}  // namespace AmanziChemistry
}  // namespace Amanzi

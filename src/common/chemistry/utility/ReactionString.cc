/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry

  Class for general forward/reverse reaction
*/

#include <sstream>

#include "ReactionString.hh"

namespace Amanzi {
namespace AmanziChemistry {

/* *******************************************************************
* Reads in a reaction string of format: reactants -> products
* Example:
*   30 A(aq) + 2 B(aq) <-> C(aq) + .3 D(aq) + -4 E(aq)
* returns
*   species = [5]("A(aq)","B(aq)","C(aq)","D(aq)","E(aq)")
*   stoichiometries = [5](-30.,-2.,1.,0.3,-4.)
* NOTE: Reactants and products have negative and positive
* stoichiometries, respectively.
******************************************************************* */
void
ParseReaction(const std::string& reactants,
              const std::string& products,
              std::vector<std::string>* species,
              std::vector<double>* stoichiometries)
{
  double coeff;
  std::string name;

  species->clear();
  stoichiometries->clear();

  // reactants
  std::istringstream iss1(reactants);
  while (iss1 >> coeff || !iss1.eof()) {
    iss1 >> name;

    species->push_back(name);
    stoichiometries->push_back(-coeff);
  }

  // products
  std::istringstream iss2(products);
  while (iss2 >> coeff || !iss2.eof()) {
    iss2 >> name;

    species->push_back(name);
    stoichiometries->push_back(coeff);
  }
}


/* *******************************************************************
* Parse surface complex reaction, reaction products are an array of
* primary species and a single surface site. The order of
* primary species and surface sites does not matter.
*
*  SpeciesName = coeff PrimaryName coeff SurfaceSite ...
******************************************************************* */
void
ParseReaction(const std::string& reaction,
              const std::vector<Species>& primary_species,
              const std::vector<SurfaceSite>& surface_sites,
              std::vector<std::string>* primary_names,
              std::vector<double>* primary_stoichiometries,
              std::vector<int>* primary_ids,
              std::string* surface_site_name,
              double* surface_site_stoichiometry,
              int* surface_site_id,
              double* h2o_stoich)
{
  double coeff;
  std::string search_name;

  std::istringstream iss(reaction);
  while (iss >> coeff || !iss.eof()) {
    iss >> search_name;

    if (search_name == "H2O") {
      *h2o_stoich = coeff;
    } else {
      // check to see if we have a primary species
      int id = -1;
      for (auto it = primary_species.begin(); it != primary_species.end(); ++it) {
        if (it->name() == search_name) {
          id = it->identifier();
          break;
        }
      }

      if (id >= 0) {
        primary_names->push_back(search_name);
        primary_stoichiometries->push_back(coeff);
        primary_ids->push_back(id);

      } else if (id < 0) {
        // did not match a primary. check to see if it is an surface site
        for (auto it = surface_sites.begin(); it != surface_sites.end(); ++it) {
          if (it->name() == search_name) {
            id = it->identifier();
            break;
          }
        }

        if (id >= 0) {
          *surface_site_name = search_name;
          *surface_site_stoichiometry = coeff;
          *surface_site_id = id;
        }

      } else {
        // did not match an surface site or primary species
        std::cout << "Reaction species \'" << search_name << "\' was not found in the primary"
                  << " species list or the surface site list.\n";
      }
    }
  }
}

} // namespace AmanziChemistry
} // namespace Amanzi

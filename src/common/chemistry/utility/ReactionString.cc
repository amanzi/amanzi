/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

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
void ParseReaction(const std::string& reactants,
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

}  // namespace AmanziChemistry
}  // namespace Amanzi

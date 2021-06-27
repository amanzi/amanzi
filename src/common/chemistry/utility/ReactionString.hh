/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Supporting non-member functions
*/

#ifndef AMANZI_CHEMISTRY_BEAKER_UTILITY_HH_
#define AMANZI_CHEMISTRY_BEAKER_UTILITY_HH_

#include <string>
#include <vector>

#include "Species.hh"
#include "SurfaceSite.hh"

namespace Amanzi {
namespace AmanziChemistry {

void ParseReaction(const std::string& reactants,
                   const std::string& products,
                   std::vector<std::string>* species,
                   std::vector<double>* stoichiometries);

void ParseReaction(const std::string& reaction,
                   const std::vector<Species>& primary_species,
                   const std::vector<SurfaceSite>& surface_sites,
                   std::vector<std::string>* primary_names,
                   std::vector<double>* primary_stoichiometries,
                   std::vector<int>* primary_ids,
                   std::string* surface_name,
                   double* surface_stoichiometry,
                   int* surface_id,
                   double* h2o_stoich);

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif

#ifndef __Glenn_hpp__
#define __Glenn_hpp__

#include "Species.hpp"
#include "AqueousEquilibriumComplex.hpp"
#include "Geochemistry.hpp"

#include <vector>

class Glenn {

public:

};

static void createCarbonateSystem(std::vector<double> *total, Geochemistry *g) {

  // primary species
  g->set_ncomp(2);

  total->push_back(1.e-6);
  total->push_back(1.e-3);

  int id = 1;
  SpeciesName name = "H+";
  double charge = 1.;
  double mol_wt = 1.0079;
  double size = 9.;
  g->addPrimarySpecies(Species(id,name,charge,mol_wt,size));

  id = 2;
  name = "HCO3-";
  charge = -1.;
  mol_wt = 61.0171;
  size = 4.;
  g->addPrimarySpecies(Species(id,name,charge,mol_wt,size));

  // secondary aqueous complexes
  std::vector<SpeciesName> species;
  std::vector<double> stoichiometries;
  std::vector<int> species_ids;
  double h2o_stoich;
  double logK;

  name = "OH-";
  species.push_back("H+");
  stoichiometries.push_back(-1.);
  species_ids.push_back(0);
  h2o_stoich = -1.;
  size = 3.5;
  charge = -1.;
  mol_wt = 17.0073;
  logK = 13.9951;
  g->addAqueousEquilibriumComplex(AqueousEquilibriumComplex(name,
                                 species,
                                 stoichiometries,
                                 species_ids,
                                 h2o_stoich,
                                 charge,mol_wt,size,logK));

  species.clear();
  stoichiometries.clear();
  species_ids.clear();

  name = "CO3--";
  species.push_back("H+");
  stoichiometries.push_back(-1.);
  species_ids.push_back(0);
  species.push_back("HCO3-");
  stoichiometries.push_back(1.);
  species_ids.push_back(1);
  h2o_stoich = 0.;
  size = 4.5;
  charge = -2.;
  mol_wt = 60.0092;
  logK = 10.3288;
  g->addAqueousEquilibriumComplex(AqueousEquilibriumComplex(name,
                                 species,
                                 stoichiometries,
                                 species_ids,
                                 h2o_stoich,
                                 charge,mol_wt,size,logK));


  species.clear();
  stoichiometries.clear();
  species_ids.clear();

  name = "CO2(aq)";
  species.push_back("H+");
  stoichiometries.push_back(1.);
  species_ids.push_back(0);
  species.push_back("HCO3-");
  stoichiometries.push_back(1.);
  species_ids.push_back(1);
  h2o_stoich = -1.;
  logK = -6.3447;
  g->addAqueousEquilibriumComplex(AqueousEquilibriumComplex(name,
                                 species,
                                 stoichiometries,
                                 species_ids,
                                 h2o_stoich,
                                 charge,mol_wt,size,logK));


};

#endif

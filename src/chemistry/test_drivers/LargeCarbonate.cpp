#include "LargeCarbonate.hpp"
#include "Geochemistry.hpp"

LargeCarbonate::LargeCarbonate(void)
  : Geochemistry()
{
  ncomp(4);
}

LargeCarbonate::~LargeCarbonate(void)
{

}

void LargeCarbonate::setup(std::vector<double> *total) {
  //
  // initial total component values
  //
  total->push_back(1.e-3);
  total->push_back(1.e-3);
  total->push_back(1.e-3);
  total->push_back(55.5084);

  //
  // primary species
  //
  int id = 0;
  SpeciesName name = "H+";
  double charge = 1.0;
  double mol_wt = 1.0079;
  double size = 9.;
  this->addPrimarySpecies(Species(id, name, charge, mol_wt, size));

  id++;
  name = "HCO3-";
  charge = -1.0;
  mol_wt = 61.0171;
  size = 4.;
  this->addPrimarySpecies(Species(id, name, charge, mol_wt, size));

  id++;
  name = "Ca++";
  charge = 2.0;
  mol_wt = 40.078;
  size = 0.0;
  this->addPrimarySpecies(Species(id, name, charge, mol_wt, size));  

  id++;
  name = "H2O";
  charge = 0.0;
  mol_wt = 18.0153;
  size = 0.0;
  this->addPrimarySpecies(Species(id, name, charge, mol_wt, size));  

  //
  // secondary species / aqueous complexes
  //
  std::vector<SpeciesName> species;
  std::vector<double> stoichiometries;
  std::vector<int> species_ids;
  double h2o_stoich;
  double logK;

  name = "OH-";
  species.push_back("H+");
  stoichiometries.push_back(-1.0);
  species_ids.push_back(0);

  species.push_back("H2O");
  stoichiometries.push_back(1.0);
  species_ids.push_back(3);

  h2o_stoich = -1.;
  size = 3.5;
  charge = -1.;
  mol_wt = 17.0073;
  logK = 13.9951;
  this->addAqueousEquilibriumComplex(AqueousEquilibriumComplex(name,
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
  this->addAqueousEquilibriumComplex(AqueousEquilibriumComplex(name,
							       species,
							       stoichiometries,
							       species_ids,
							       h2o_stoich,
							       charge,mol_wt,size,logK));


  species.clear();
  stoichiometries.clear();
  species_ids.clear();

  name = "H2CO3*";

  species.push_back("H+");
  stoichiometries.push_back(1.0);
  species_ids.push_back(0);

  species.push_back("HCO3-");
  stoichiometries.push_back(1.0);
  species_ids.push_back(1);

  h2o_stoich = 0.0;
  size = 0.0;
  charge = 0.0;
  mol_wt = 62.0251;
  logK = -6.3447;
  this->addAqueousEquilibriumComplex(AqueousEquilibriumComplex(name,
							       species,
							       stoichiometries,
							       species_ids,
							       h2o_stoich,
							       charge,mol_wt,size,logK));


  species.clear();
  stoichiometries.clear();
  species_ids.clear();

  name = "CaOH+";

  species.push_back("H2O");
  stoichiometries.push_back(1.0);
  species_ids.push_back(3);

  species.push_back("H+");
  stoichiometries.push_back(-1.0);
  species_ids.push_back(0);

  species.push_back("Ca++");
  stoichiometries.push_back(1.0);
  species_ids.push_back(2);

  h2o_stoich = -1.0;
  size = 0.0;
  charge = 1.0;
  mol_wt = 57.0853;
  logK = 12.78;
  this->addAqueousEquilibriumComplex(AqueousEquilibriumComplex(name,
							       species,
							       stoichiometries,
							       species_ids,
							       h2o_stoich,
							       charge,mol_wt,size,logK));


  species.clear();
  stoichiometries.clear();
  species_ids.clear();

  name = "CaHCO3+";

  species.push_back("HCO3-");
  stoichiometries.push_back(1.0);
  species_ids.push_back(1);

  species.push_back("Ca++");
  stoichiometries.push_back(1.0);
  species_ids.push_back(2);

  h2o_stoich = 1.0;
  size = 0.0;
  charge = 1.0;
  mol_wt = 101.095;
  logK = -1.043;
  this->addAqueousEquilibriumComplex(AqueousEquilibriumComplex(name,
							       species,
							       stoichiometries,
							       species_ids,
							       h2o_stoich,
							       charge,mol_wt,size,logK));


  species.clear();
  stoichiometries.clear();
  species_ids.clear();

  name = "CaCO3(aq)";

  species.push_back("HCO3-");
  stoichiometries.push_back(1.0);
  species_ids.push_back(1);

  species.push_back("H+");
  stoichiometries.push_back(-1.0);
  species_ids.push_back(0);

  species.push_back("Ca++");
  stoichiometries.push_back(1.0);
  species_ids.push_back(2);

  h2o_stoich = 1.0;
  size = 0.0;
  charge = 0.0;
  mol_wt = 100.087;
  logK = 7.0088;
  this->addAqueousEquilibriumComplex(AqueousEquilibriumComplex(name,
							       species,
							       stoichiometries,
							       species_ids,
							       h2o_stoich,
							       charge,mol_wt,size,logK));


};

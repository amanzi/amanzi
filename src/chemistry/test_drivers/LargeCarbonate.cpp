/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <vector>

#include "AqueousEquilibriumComplex.hpp"
#include "LargeCarbonate.hpp"
#include "Beaker.hpp"

LargeCarbonate::LargeCarbonate(void)
    : Beaker()
{
  ncomp(3);
  //ncomp(4);
}  // end LargeCarbonate constructor

LargeCarbonate::~LargeCarbonate(void)
{
}  // end LargeCarbonate constructor

void LargeCarbonate::setup(std::vector<double> &total, 
                           const std::string mineral_kinetics_file) {
  this->resize(ncomp());
  //
  // initial total component values
  //
  total.push_back(1.e-3);
  total.push_back(1.e-3);
  total.push_back(1.e-3);
  //total.push_back(55.5084);

  //
  // primary species
  //
  int id = 0;
  SpeciesName name = "H+";
  double size = 9.0;
  double charge = 1.0;
  double mol_wt = 1.0079;
  this->addPrimarySpecies(Species(id, name, charge, mol_wt, size));

  id++;
  name = "HCO3-";
  size = 4.0;
  charge = -1.0;
  mol_wt = 61.0171;
  this->addPrimarySpecies(Species(id, name, charge, mol_wt, size));

  id++;
  name = "Ca++";
  size = 6.0;
  charge = 2.0;
  mol_wt = 40.078;
  this->addPrimarySpecies(Species(id, name, charge, mol_wt, size));

//   id++;
//   name = "H2O";
//   size = 3.0;
//   charge = 0.0;
//   mol_wt = 18.0153;
//   this->addPrimarySpecies(Species(id, name, charge, mol_wt, size));

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

//   species.push_back("H2O");
//   stoichiometries.push_back(1.0);
//   species_ids.push_back(3);

  h2o_stoich = 1.0;
  logK = 13.9951;
  size = 3.5;
  charge = -1.0;
  mol_wt = 17.0073;
  AqueousEquilibriumComplex oh(name,
                               species,
                               stoichiometries,
                               species_ids,
                               h2o_stoich,
                               charge, mol_wt, size, logK);
  this->addAqueousEquilibriumComplex(oh);


  species.clear();
  stoichiometries.clear();
  species_ids.clear();

  name = "CO3--";

  species.push_back("H+");
  stoichiometries.push_back(-1.0);
  species_ids.push_back(0);

  species.push_back("HCO3-");
  stoichiometries.push_back(1.0);
  species_ids.push_back(1);

  h2o_stoich = 0.0;
  logK = 10.3288;
  size = 4.5;
  charge = -2.0;
  mol_wt = 60.0092;
  AqueousEquilibriumComplex co3(name,
                                species,
                                stoichiometries,
                                species_ids,
                                h2o_stoich,
                                charge, mol_wt, size, logK);
  this->addAqueousEquilibriumComplex(co3);


  species.clear();
  stoichiometries.clear();
  species_ids.clear();

  name = "CO2(aq)";

  species.push_back("H+");
  stoichiometries.push_back(1.0);
  species_ids.push_back(0);

  species.push_back("HCO3-");
  stoichiometries.push_back(1.0);
  species_ids.push_back(1);

  h2o_stoich = -1.0;
  logK = -6.3447;
  size = 3.0;
  charge = 0.0;
  mol_wt = 44.0098;
  AqueousEquilibriumComplex h2co3(name,
                                  species,
                                  stoichiometries,
                                  species_ids,
                                  h2o_stoich,
                                  charge, mol_wt, size, logK);
  this->addAqueousEquilibriumComplex(h2co3);


  species.clear();
  stoichiometries.clear();
  species_ids.clear();

  name = "CaOH+";

//   species.push_back("H2O");
//   stoichiometries.push_back(1.0);
//   species_ids.push_back(3);

  species.push_back("H+");
  stoichiometries.push_back(-1.0);
  species_ids.push_back(0);

  species.push_back("Ca++");
  stoichiometries.push_back(1.0);
  species_ids.push_back(2);

  h2o_stoich = 1.0;
  logK = 12.85;
  size = 4.0;
  charge = 1.0;
  mol_wt = 57.0853;
  AqueousEquilibriumComplex caoh(name,
                                 species,
                                 stoichiometries,
                                 species_ids,
                                 h2o_stoich,
                                 charge, mol_wt, size, logK);
  this->addAqueousEquilibriumComplex(caoh);


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

  h2o_stoich = 0.0;
  logK = -1.0467;
  size = 4.0;
  charge = 1.0;
  mol_wt = 101.0951;
  AqueousEquilibriumComplex cahco3(name,
                                   species,
                                   stoichiometries,
                                   species_ids,
                                   h2o_stoich,
                                   charge, mol_wt, size, logK);
  this->addAqueousEquilibriumComplex(cahco3);


  species.clear();
  stoichiometries.clear();
  species_ids.clear();

  name = "CaCO3(aq)";

  species.push_back("H+");
  stoichiometries.push_back(-1.0);
  species_ids.push_back(0);

  species.push_back("HCO3-");
  stoichiometries.push_back(1.0);
  species_ids.push_back(1);

  species.push_back("Ca++");
  stoichiometries.push_back(1.0);
  species_ids.push_back(2);

  h2o_stoich = 0.0;
  logK = 7.0017;
  size = 3.0;
  charge = 0.0;
  mol_wt = 100.0872;
  AqueousEquilibriumComplex caco3(name,
                                  species,
                                  stoichiometries,
                                  species_ids,
                                  h2o_stoich,
                                  charge, mol_wt, size, logK);
  this->addAqueousEquilibriumComplex(caco3);

  this->SetupMineralKinetics(mineral_kinetics_file);
}  // end setup()

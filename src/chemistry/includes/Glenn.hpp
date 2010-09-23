#ifndef __Glenn_hpp__
#define __Glenn_hpp__

#include "Species.hpp"
#include "AqueousEquilibriumComplex.hpp"
#include "Beaker.hpp"
#include "Geochemistry.hpp"
#include "FileIO.hpp"

#include <vector>

class Glenn {

public:

};

static void createCarbonateSystem(std::vector<double> *total, Beaker *g) {

  // primary species
  g->ncomp(2);

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

static void readChemistryFromFile(string filename, Beaker *g) {

  char word[32];
  int ncomp, neqcplx;
  int tempi;

  // open file with FileIO buffer
  FileIO *file = new FileIO(filename);
  // first line indicates number of primary and secondary components
  file->getLine();
  file->readInt(&ncomp);
  g->ncomp(ncomp);
  file->readInt(&neqcplx);

  int id;
  SpeciesName name;
  double charge;
  double mol_wt = 0.;
  double size;

  // for ncomp primary speces, each line lists name, Z, a0
  for (int i=0; i<ncomp; i++) {
    id = i;
    file->getLine();
    file->readWord(word);
    name = word;
    file->readDouble(&charge);
    file->readDouble(&size);
#define DEBUG
#ifdef DEBUG
    std::cout << name << " " << charge << " " << size << endl;
#endif
    g->addPrimarySpecies(Species(id,name,charge,mol_wt,size));
  }

#ifdef DEBUG
  std::cout << "-------------------------\n";
#endif
  // secondary aqueous complexes
  std::vector<SpeciesName> species;
  std::vector<double> stoichiometries;
  std::vector<int> species_ids;
  double h2o_stoich;
  double logK;
  // for neqcplx secondary speces, each line lists name, Z, a0, etc.
  for (int i=0; i<neqcplx; i++) {

    species.clear();
    stoichiometries.clear();
    species_ids.clear();

    file->getLine();
    file->readWord(word);
    name = word;
    file->readDouble(&charge);
    file->readDouble(&size);
#ifdef DEBUG
    std::cout << name << " " << charge << " " << size << endl;
#endif
    // species ids
    file->getLine();
    int ncomp_in_rxn;
    file->readInt(&ncomp_in_rxn);
    for (int j=0; j<ncomp_in_rxn; j++) {
      int tempi;
      file->readInt(&tempi);
      // decrement for zero-based indexing
      tempi--;
      species_ids.push_back(tempi);
    }
#ifdef DEBUG
    for (int j=0; j<(int)species_ids.size(); j++)
      std::cout << species_ids[j] << " " ;
    std::cout << endl;
#endif
    // species stoichiometries
    double tempd;
    file->getLine();
    // skip first value (PFLOTRAN eqcplxstroich mirrors eqcplxspecid)
    file->readDouble(&tempd);
    for (int j=0; j<ncomp_in_rxn; j++) {
      file->readDouble(&tempd);
      stoichiometries.push_back(tempd);
    }
#ifdef DEBUG
    for (int j=0; j<ncomp_in_rxn; j++)
      std::cout << stoichiometries[j] << " " ;
    std::cout << endl;
#endif
    // h2o id <- skip this
    file->getLine();
    file->readInt(&tempi);
#ifdef DEBUG
//    std::cout << eqcplxh2oid[i] << endl;
#endif
    // h2o stoichiometry
    file->getLine();
    file->readDouble(&h2o_stoich);
#ifdef DEBUG
    std::cout << h2o_stoich << endl;
#endif
    // logK
    file->getLine();
    file->readDouble(&logK);
#ifdef DEBUG
    std::cout << logK << endl;
    std::cout << "-------------------------\n";
#endif
    g->addAqueousEquilibriumComplex(AqueousEquilibriumComplex(name,
                                    species,
                                    stoichiometries,
                                    species_ids,
                                    h2o_stoich,
                                    charge,mol_wt,size,logK));
  }
  
  delete file;
};

static void readTargetTotalFromFile(string filename, int ncomp, 
                                    std::vector<double> *total) {

  total->clear();

  // open file with FileIO buffer
  FileIO *file = new FileIO(filename);
  // first line indicates number of primary and secondary components
  char word[32];
  double temp;
  for (int i=0; i<ncomp; i++) {
    file->getLine();
    file->readWord(word);
    file->readDouble(&temp);
    file->readDouble(&temp);
    total->push_back(temp);
  }
  delete file;

};

#endif

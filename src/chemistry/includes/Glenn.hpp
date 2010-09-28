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
   Glenn(Beaker *b);
   ~Glenn();

   void solve(std::vector<double> &total, double final_time, double ts_size,
              double porosity, double saturation, double water_density, 
              double volume);

 private:
   Beaker *b_;

};

static void createCarbonateSystem(std::vector<double> *total, Beaker *g) {

  // primary species
  g->resize(2);

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

}; // end createCarbonateSystem

static void readChemistryFromFile(string filename, Beaker *g) 
{

  char word[32];
  int ncomp, neqcplx, ngeneral_rxn;

  // open file with FileIO buffer
  FileIO *file = new FileIO(filename);
  // first line indicates number of primary and secondary components
  file->getLine();
  file->readInt(&ncomp);
  g->resize(ncomp);
  file->readInt(&neqcplx);
  file->readInt(&ngeneral_rxn);

  int id;
  SpeciesName name;
  double charge;
  double mol_wt = 0.;
  double size;

  // for ncomp primary speces, each line lists name, Z, a0
  for (int i = 0; i < ncomp; i++) {
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
  std::cout << "--- Aqueous Complexes----\n";
#endif
  // secondary aqueous complexes
  std::vector<SpeciesName> species;
  std::vector<double> stoichiometries;
  std::vector<int> species_ids;
  double h2o_stoich;
  double logK;
  // for neqcplx secondary speces, each line lists name, Z, a0, etc.
  for (int i = 0; i < neqcplx; i++) {

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
    for (int j = 0; j<ncomp_in_rxn; j++) {
      int tempi;
      file->readInt(&tempi);
      // decrement for zero-based indexing
      tempi--;
      species_ids.push_back(tempi);
    }
#ifdef DEBUG
    for (int j = 0; j < (int)species_ids.size(); j++)
      std::cout << species_ids[j] << " " ;
    std::cout << endl;
#endif
    // species stoichiometries
    double tempd;
    file->getLine();
    // skip first value (PFLOTRAN eqcplxstroich mirrors eqcplxspecid)
    file->readDouble(&tempd);
    for (int j = 0; j<ncomp_in_rxn; j++) {
      file->readDouble(&tempd);
      stoichiometries.push_back(tempd);
    }
#ifdef DEBUG
    for (int j = 0; j < ncomp_in_rxn; j++)
      std::cout << stoichiometries[j] << " " ;
    std::cout << endl;
#endif
    // h2o id <- skip this
    int tempi;
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
                                    species,stoichiometries,
                                    species_ids,h2o_stoich,
                                    charge,mol_wt,size,logK));
  }

#ifdef DEBUG
  std::cout << "--- General Reactions ---\n";
#endif

/*
  GeneralRxn(SpeciesName name, 
             std::vector<SpeciesName>species,
             std::vector<double>stoichiometries,
             std::vector<int>species_ids,
             std::vector<double>forward_stoichiometries,
             std::vector<int>forward_species_ids,
             std::vector<double>backward_stoichiometries,
             std::vector<int>backward_species_ids,
             double kf, double kb); */

  // secondary aqueous complexes
  name = "";
  species.clear();
  stoichiometries.clear();
  species_ids.clear();
  std::vector<double> forward_stoichiometries;
  std::vector<int> forward_species_ids;
  std::vector<double> backward_stoichiometries;
  std::vector<int> backward_species_ids;
  double kf;
  double kb;

  // for neqcplx secondary speces, each line lists name, Z, a0, etc.
  for (int i = 0; i < ngeneral_rxn; i++) {

    // name currently not used
    species.clear(); // currently not used
    stoichiometries.clear();
    species_ids.clear();
    forward_stoichiometries.clear();
    forward_species_ids.clear();
    backward_stoichiometries.clear();
    backward_species_ids.clear();

#ifdef DEBUG
    std::cout << "Species IDs" << endl;
#endif
    // species ids
    file->getLine();
    int ncomp_in_rxn;
    file->readInt(&ncomp_in_rxn);
    for (int j = 0; j<ncomp_in_rxn; j++) {
      int tempi;
      file->readInt(&tempi);
      // decrement for zero-based indexing
      tempi--;
      species_ids.push_back(tempi);
    }
#ifdef DEBUG
    for (int j = 0; j < (int)species_ids.size(); j++)
      std::cout << species_ids[j] << " " ;
    std::cout << endl;
#endif
    // species stoichiometries
    double tempd;
    file->getLine();
    for (int j = 0; j<ncomp_in_rxn; j++) {
      file->readDouble(&tempd);
      stoichiometries.push_back(tempd);
    }
#ifdef DEBUG
    for (int j = 0; j < ncomp_in_rxn; j++)
      std::cout << stoichiometries[j] << " " ;
    std::cout << endl;
#endif


#ifdef DEBUG
    std::cout << "Forward Reaction" << endl;
#endif
    // species ids
    file->getLine();
    file->readInt(&ncomp_in_rxn);
    for (int j = 0; j<ncomp_in_rxn; j++) {
      int tempi;
      file->readInt(&tempi);
      // decrement for zero-based indexing
      tempi--;
      forward_species_ids.push_back(tempi);
    }
#ifdef DEBUG
    std::cout << "IDs\n";
    for (int j = 0; j < ncomp_in_rxn; j++)
      std::cout << forward_species_ids[j] << " " ;
    std::cout << endl;
#endif
    // species stoichiometries
    file->getLine();
    for (int j = 0; j<ncomp_in_rxn; j++) {
      file->readDouble(&tempd);
      forward_stoichiometries.push_back(tempd);
    }
#ifdef DEBUG
    std::cout << "Stoichiometry\n";
    for (int j = 0; j < ncomp_in_rxn; j++)
      std::cout << forward_stoichiometries[j] << " " ;
    std::cout << endl;
#endif


#ifdef DEBUG
    std::cout << "Backward Reaction" << endl;
#endif
    // species ids
    file->getLine();
    file->readInt(&ncomp_in_rxn);
    for (int j = 0; j<ncomp_in_rxn; j++) {
      int tempi;
      file->readInt(&tempi);
      // decrement for zero-based indexing
      tempi--;
      backward_species_ids.push_back(tempi);
    }
#ifdef DEBUG
    std::cout << "IDs\n";
    for (int j = 0; j < ncomp_in_rxn; j++)
      std::cout << backward_species_ids[j] << " " ;
    std::cout << endl;
#endif
    // species stoichiometries
    file->getLine();
    for (int j = 0; j<ncomp_in_rxn; j++) {
      file->readDouble(&tempd);
      backward_stoichiometries.push_back(tempd);
    }
#ifdef DEBUG
    std::cout << "Stoichiometry\n";
    for (int j = 0; j < ncomp_in_rxn; j++)
      std::cout << backward_stoichiometries[j] << " " ;
    std::cout << endl;
#endif

    // h2o stoichiometry
    file->getLine();
    file->readDouble(&h2o_stoich);
#ifdef DEBUG
    std::cout << "H2O ID" << endl;
    std::cout << h2o_stoich << endl;
#endif

    // kf
    file->getLine();
    file->readDouble(&kf);
    // kb
    file->getLine();
    file->readDouble(&kb);
#ifdef DEBUG
    std::cout << "kf: " << kf << endl;
    std::cout << "kb: " << kb << endl;
    std::cout << "-------------------------\n";
#endif

    g->addGeneralRxn(GeneralRxn(name,species,
                                stoichiometries,species_ids,
                                forward_stoichiometries,forward_species_ids,
                                backward_stoichiometries,backward_species_ids,
//                                h2o_stoich,
                                kf,kb));
  }
  
  delete file;
}; // end readChemistryFromFile

static void readTargetTotalFromFile(string filename, int ncomp, 
                                    std::vector<double> *total) 
{
  total->clear();

  // open file with FileIO buffer
  FileIO *file = new FileIO(filename);
  // first line indicates number of primary and secondary components
  char word[32];
  double temp;
  for (int i = 0; i < ncomp; i++) {
    file->getLine();
    file->readWord(word);
    file->readDouble(&temp);
    file->readDouble(&temp);
    total->push_back(temp);
  }
  delete file;
}; // end readTargetTotalFromFile

#endif // __Glenn_hpp__

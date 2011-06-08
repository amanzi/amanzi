#ifndef __Glenn_hpp__
#define __Glenn_hpp__

#include "Species.hpp"
#include "AqueousEquilibriumComplex.hpp"
#include "Beaker.hpp"
#include "FileIO.hpp"
#include "Mineral.hpp"
#include "KineticRate.hpp"
#include "KineticRateTST.hpp"
#include "string-tokenizer.hh"
#include "MineralKineticsFactory.hpp"


#include <vector>

class Glenn {

 public:
   Glenn(Beaker *b);
   ~Glenn();

  void solve(Beaker::BeakerComponents* components, 
             double final_time, double ts_size,
             const Beaker::BeakerParameters& parameters);

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
  id = 0;
  species.push_back("H+");
  stoichiometries.push_back(-1.);
  species_ids.push_back(0);
  h2o_stoich = -1.;
  size = 3.5;
  charge = -1.;
  mol_wt = 17.0073;
  logK = 13.9951;
  g->addAqueousEquilibriumComplex(AqueousEquilibriumComplex(name,
                                 id,
                                 species,
                                 stoichiometries,
                                 species_ids,
                                 h2o_stoich,
                                 charge,mol_wt,size,logK));

  species.clear();
  stoichiometries.clear();
  species_ids.clear();

  name = "CO3--";
  id++;
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
                                 id,
                                 species,
                                 stoichiometries,
                                 species_ids,
                                 h2o_stoich,
                                 charge,mol_wt,size,logK));

  species.clear();
  stoichiometries.clear();
  species_ids.clear();

  name = "CO2(aq)";
  id++;
  species.push_back("H+");
  stoichiometries.push_back(1.);
  species_ids.push_back(0);
  species.push_back("HCO3-");
  stoichiometries.push_back(1.);
  species_ids.push_back(1);
  h2o_stoich = -1.;
  size = 3.;
  charge = 0.;
  mol_wt = 44.0098;
  logK = -6.3447;
  g->addAqueousEquilibriumComplex(AqueousEquilibriumComplex(name,
                                 id,
                                 species,
                                 stoichiometries,
                                 species_ids,
                                 h2o_stoich,
                                 charge,mol_wt,size,logK));

}; // end createCarbonateSystem

static void readChemistryFromFile(string filename, Beaker *g) 
{

  char word[32];
  int ncomp, neqcplx, ngeneral_rxn, nsurfcplx_rxn, nminerals;

  // open file with FileIO buffer
  FileIO *file = new FileIO(filename);
  // first line indicates number of primary and secondary components
  file->getLine();
  file->readInt(&ncomp);
  g->resize(ncomp);
  file->readInt(&neqcplx);
  file->readInt(&ngeneral_rxn);
  file->readInt(&nsurfcplx_rxn);
  file->readInt(&nminerals);

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

  std::vector<std::string> primary_names;
  g->GetPrimaryNames(&primary_names);

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
  for (int irxn = 0; irxn < neqcplx; irxn++) {

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
      species.push_back(primary_names[tempi]);
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
                                    irxn,
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
      species.push_back(primary_names[tempi]);
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

#ifdef DEBUG
  std::cout << "--- Surface Complexation Reactions ---\n";
#endif

  // for neqcplx secondary speces, each line lists name, Z, a0, etc.
  for (int irxn = 0; irxn < nsurfcplx_rxn; irxn++) {

    // surface complex name
    file->getLine();
    file->readWord(word);
    name = word;
    // surface site density
    double site_density;
    file->getLine();
    file->readDouble(&site_density);

    SurfaceSite *surface_site = new SurfaceSite(name,irxn,site_density);

    int num_surface_complexes;
    file->getLine();
    file->readInt(&num_surface_complexes);

    std::vector<SurfaceComplex> surface_complexes;
    for (int icplx = 0; icplx < num_surface_complexes; icplx++) {

      species.clear(); // currently not used
      stoichiometries.clear();
      species_ids.clear();

      // surface complex name
      file->getLine();
      file->readWord(word);
      name = word;

#ifdef DEBUG
      std::cout << "Surface Complex: " << name << endl;
#endif

      // charge
      file->readDouble(&charge);

#ifdef DEBUG
      std::cout << "Surface Complex Species IDs" << endl;
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
        species.push_back(primary_names[tempi]);
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

      // skip h2o id
      file->getLine();

      // h2o stoichiometry
      file->getLine();
      file->readDouble(&h2o_stoich);
#ifdef DEBUG
      std::cout << "H2O ID" << endl;
      std::cout << h2o_stoich << endl;
#endif
      
      // skip free site id
      file->getLine();

      // free site stoichiometry
      double free_site_stoich;
      file->getLine();
      file->readDouble(&free_site_stoich);
#ifdef DEBUG
      std::cout << "Free Site Stoichiometry" << endl;
      std::cout << free_site_stoich << endl;
#endif

      // logK
      file->getLine();
      file->readDouble(&logK);
#ifdef DEBUG
      std::cout << logK << endl;
      std::cout << "-------------------------\n";
#endif
     
      surface_complexes.push_back(
         SurfaceComplex(name, icplx, species, stoichiometries, species_ids,
                        h2o_stoich, free_site_stoich, charge, logK));

    }
    // have to pass in a new object here as the SurfaceSite destructor fails when
    // deleting the vector storing mineral pointers
    g->addSurfaceComplexationRxn(SurfaceComplexationRxn(surface_site,surface_complexes));
    surface_complexes.clear();
  }

#ifdef DEBUG
  std::cout << "--- Kinetic Minerals ----\n";
#endif
  // secondary aqueous complexes
  // for neqcplx secondary speces, each line lists name, Z, a0, etc.
  for (int irxn = 0; irxn < nminerals; irxn++) {

    species.clear();
    stoichiometries.clear();
    species_ids.clear();

    file->getLine();
    file->readWord(word);
    name = word;
#ifdef DEBUG
    std::cout << name << endl;
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
      species.push_back(primary_names[tempi]);
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
    // skip first value (PFLOTRAN kinmnrlstroich mirrors kinmnrlspecid)
// not longer the case - geh    file->readDouble(&tempd);
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
//    std::cout << kinmnrlh2oid[i] << endl;
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
#endif
    // molar vol
    double molar_vol;
    file->getLine();
    file->readDouble(&molar_vol);
#ifdef DEBUG
    std::cout << molar_vol << endl;
#endif
    // molar wt
    double molar_wt;
    file->getLine();
    file->readDouble(&molar_wt);
#ifdef DEBUG
    std::cout << molar_wt << endl;
#endif
    // rate
    double rate;
    file->getLine();
    file->readDouble(&rate);
#ifdef DEBUG
    std::cout << rate << endl;
#endif
    // specific surface area
    double specific_surface_area;
    file->getLine();
    file->readDouble(&specific_surface_area);
    // convert from cm^2/cm^3 to m^2/m^3
    specific_surface_area *= 100.;
#ifdef DEBUG
    std::cout << specific_surface_area << endl;
#endif


#ifdef DEBUG
    std::cout << "-------------------------\n";
#endif

    Mineral mineral(name, irxn,
                    species,
                    stoichiometries,
                    species_ids,
                    h2o_stoich,
                    molar_wt, 
                    logK, 
                    molar_vol,
                    specific_surface_area);
    g->addMineral(mineral);


    MineralKineticsFactory mkf;
//    mkf.set_verbosity(verbosity());
//    SpeciesId mineral_id = mkf.VerifyMineralName(name, minerals());
//    Mineral mineral = minerals().at(mineral_id);
    std::string rate_type = "TST";
    std::string str1("log10_rate_constant ");
    std::ostringstream oss;
    oss << rate*10000.;
    std::string str3(" moles_m2_sec");

    std::string rate_string = str1 + oss.str() + str3;
    StringTokenizer rate_data(rate_string, "");
//    KineticRate* kinetic_rate = mkf.Create(rate_type, rate_data, mineral, g->primary_species());

//    g->AddMineralKineticRate(kinetic_rate);
  }

  delete file;
}; // end readChemistryFromFile

static void readTargetTotalFromFile(string filename, int ncomp, 
                                    std::vector<double> *total) 
{
  total->clear();

  // open file with FileIO buffer
  FileIO *file = new FileIO(filename);
  char word[32];
  double temp;
  for (int i = 0; i < ncomp; i++) {
    file->getLine();
    file->readWord(word);    // name of species
    file->readDouble(&temp); // free ion concentration [m]
    file->readDouble(&temp); // total component concentration [m]
    total->push_back(temp);
  }
  delete file;
}; // end readTargetTotalFromFile

static void readTargetFreeIonFromFile(string filename, int ncomp, 
                                      std::vector<double> *free_ion) 
{
  free_ion->clear();

  // open file with FileIO buffer
  FileIO *file = new FileIO(filename);
  char word[32];
  double temp;
  for (int i = 0; i < ncomp; i++) {
    file->getLine();
    file->readWord(word);     // name of species
    file->readDouble(&temp);  // free ion concentration [m]
    free_ion->push_back(temp);
  }
  delete file;
}; // end readTargetFreeIonFromFile

static void readMineralVolFracFromFile(string filename, 
                                        const std::vector<Mineral>& minerals,
                                        std::vector<double> *mineral_vol_frac) 
{
  int nmnrl = minerals.size();
  mineral_vol_frac->clear();
  mineral_vol_frac->resize(nmnrl);
  for (int imnrl = 0; imnrl < nmnrl; imnrl++)
    (*mineral_vol_frac)[imnrl] = 0.;

  // open file with FileIO buffer
  FileIO *file = new FileIO(filename);
  char word[32];
  double temp;
  while (file->getLine() != 0) {
    file->readWord(word);     // name of mineral
    int imnrl = -1;
    for (int i = 0; i < nmnrl; i++) {
      if (minerals[i].name().compare(word) == 0) imnrl = i;
    }
    if (imnrl > -1) {
      file->readDouble(&temp);  // volume fraction [-]
      (*mineral_vol_frac)[imnrl] = temp;
      file->readDouble(&temp);  // surface area [m^2/m^3]
    }
  }
  delete file;
}; // end readMineralsVolFracFromFile

#endif // __Glenn_hpp__

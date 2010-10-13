/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __ThermoDatabase_hpp__
#define __ThermoDatabase_hpp__

#include <vector>

#include "Species.hpp"
#include "Beaker.hpp"

class ThermoDatabase : public Beaker {
 public:
  ThermoDatabase();
  virtual ~ThermoDatabase();

  void setup(std::vector<double> &total, const Beaker::BeakerParameters parameters);

 protected:
  void ReadFile(const std::string thermo_filename);
  void ParsePrimarySpecies(const std::string data);
  void ParseAqueousEquilibriumComplex(const std::string data);
  void ParseMineral(const std::string data);
  void ParseReaction(const std::string reaction, 
                     std::string *name,
                     std::vector<SpeciesName>* species, 
                     std::vector<double>* stoichiometries, 
                     std::vector<int>* species_ids, 
                     double* h2o_stoich);

 private:
  SpeciesId primary_id_;
  SpeciesId mineral_id_;
};

#endif  // __ThermoDatabase_hpp__

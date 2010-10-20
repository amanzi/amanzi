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

  void Setup(const Beaker::BeakerComponents& components,
             const Beaker::BeakerParameters& parameters);

 protected:
  void ReadFile(const std::string thermo_filename);
  void ParsePrimarySpecies(const std::string data);
  void ParseAqueousEquilibriumComplex(const std::string data);
  void ParseMineral(const std::string data);
  void ParseIonExchangeSite(const std::string data);
  void ParseIonExchangeComplex(const std::string data);
  void ParseReaction(const std::string reaction, 
                     std::string *name,
                     std::vector<SpeciesName>* species, 
                     std::vector<double>* stoichiometries, 
                     std::vector<int>* species_ids, 
                     double* h2o_stoich);
  void ParseReaction(const std::string reaction, 
                     std::string* name,
                     std::vector<SpeciesName>* primaries, 
                     std::vector<double>* primary_stoichiometries, 
                     std::vector<SpeciesId>* primary_ids, 
                     std::vector<SpeciesName>* exchange_sites, 
                     std::vector<double>* exchanger_stoichiometries, 
                     std::vector<SpeciesId>* exchanger_ids, 
                     double* h2o_stoich);

 private:
  SpeciesId primary_id_;
  SpeciesId aqueous_equilibrium_complex_id_;
  SpeciesId mineral_id_;
  SpeciesId ion_exchange_site_id_;
  SpeciesId ion_exchange_complex_id_;
};

#endif  // __ThermoDatabase_hpp__

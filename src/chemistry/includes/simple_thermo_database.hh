/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_SIMPLETHERMODATABASE_HH_
#define AMANZI_CHEMISTRY_SIMPLETHERMODATABASE_HH_

#include <vector>
#include <string>

#include "species.hh"
#include "beaker.hh"

namespace amanzi {
namespace chemistry {

class SimpleThermoDatabase : public Beaker {
 public:
  SimpleThermoDatabase();
  virtual ~SimpleThermoDatabase();

  void Setup(const Beaker::BeakerComponents& components,
             const Beaker::BeakerParameters& parameters);

 protected:
  void ReadFile(const std::string& thermo_filename);
  void ParsePrimarySpecies(const std::string& data);
  void ParseAqueousEquilibriumComplex(const std::string& data);
  void ParseGeneralKinetics(const std::string& data);
  void ParseMineral(const std::string& data);
  void ParseMineralKinetics(const std::string& data);
  void ParseIonExchangeSite(const std::string& data);
  void ParseIonExchangeComplex(const std::string& data);
  void ParseSurfaceComplexSite(const std::string& data);
  void ParseSurfaceComplex(const std::string& data);
  void ParseSorptionIsotherm(const std::string& data);
  void FinishSurfaceComplexation(void);
  void ParseReaction(const std::string& reaction,
                     std::string* name,
                     std::vector<SpeciesName>* species,
                     std::vector<double>* stoichiometries,
                     std::vector<int>* species_ids,
                     double* h2o_stoich);
  void ParseReactionString(const std::string reaction,
                           std::vector<SpeciesName>* species,
                           std::vector<double>* stoichiometries);
  int SpeciesNameToID(const std::string species_name);
  void RemoveLeadingAndTrailingSpaces(std::string* s);
  void ParseIonExchangeReaction(const std::string& reaction,
                                std::string* name,
                                SpeciesName* primary_name,
                                SpeciesId* primary_id);

  void ParseSurfaceComplexReaction(const std::string& reaction,
                                   std::string* name,
                                   std::vector<SpeciesName>* primary_name,
                                   std::vector<double>* primary_stoichiometry,
                                   std::vector<SpeciesId>* primary_id,
                                   SpeciesName* surface_name,
                                   double* surface_stoichiometry,
                                   SpeciesId* surface_id,
                                   double* h2o_stoich);
 private:
  SpeciesId primary_id_;
  SpeciesId aqueous_equilibrium_complex_id_;
  SpeciesId mineral_id_;
  SpeciesId ion_exchange_site_id_;
  SpeciesId ion_exchange_complex_id_;
  SpeciesId surface_site_id_;
  SpeciesId surface_complex_id_;
  SpeciesId surface_complexation_rxn_id_;

  std::vector<SurfaceSite> surface_sites_;
  std::vector<SurfaceComplexationRxn> surface_complexation_reactions_;
};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_SIMPLETHERMODATABASE_HH_

/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __IonExchangeComplex_hpp__
#define __IonExchangeComplex_hpp__

/*
**  Class for ion exchange complexation reaction
**
**  NaX <===> Na+ + X-
**
*/


#include <string>
#include <vector>
#include <cmath>

#include "Species.hpp"
#include "IonExchangeSite.hpp"
#include "SecondarySpecies.hpp"
#include "Block.hpp"

class IonExchangeComplex : public SecondarySpecies {

 public:
  IonExchangeComplex();
  IonExchangeComplex(const SpeciesName name,
                     SpeciesId complex_id,
                     std::vector<SpeciesName> primary_species,
                     std::vector<double> primary_stoichiometries,
                     std::vector<SpeciesId> primary_ids,
                     std::vector<SpeciesName> exchange_site_names,
                     std::vector<double> exchange_site_stoichiometries,
                     std::vector<SpeciesId> exchange_site_ids,
                     const double h2o_stoich,
                     const double log_Keq);
  ~IonExchangeComplex();

  // update molalities
  void Update(const std::vector<Species>primary_species, 
              const std::vector<IonExchangeSite>exchange_sites);
  // add stoichiometric contribution of complex to total
  void AddContributionToTotal(std::vector<double> &total);
  // add derivative of total with respect to free-ion to dtotal
  void AddContributionToDTotal(const std::vector<Species> primary_species,
                               Block *dtotal);

  void display(void) const;
  void Display(void) const;
  void DisplayReaction(void) const;
  void DisplayResultsHeader(void) const;
  void DisplayResults(void) const;

 protected:

 private:
  std::vector<SpeciesName> exchange_site_names_;
  std::vector<double> exchange_site_stoichiometries_;
  std::vector<SpeciesId> exchange_site_ids_;


};

#endif // __IonExchangeComplex_hpp__

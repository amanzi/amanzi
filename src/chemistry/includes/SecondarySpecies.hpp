/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __Secondary_Species_hpp__
#define __Secondary_Species_hpp__


#include <string>
#include <vector>
#include <cmath>

#include "Species.hpp"
#include "Block.hpp"

// Base class for secondary species (aqueous equilibrium complexes,
// minerals)

class SecondarySpecies : public Species {

 public:
  SecondarySpecies();
  SecondarySpecies(const SpeciesName name,
                   const SpeciesId secondary_id,
                   const std::vector<SpeciesName> species,
                   const std::vector<double> stoichiometries,
                   const std::vector<int> species_ids,
                   const double h2o_stoich, 
                   const double charge, 
                   const double mol_wt,
                   const double size,
                   const double logK);

  virtual ~SecondarySpecies();

  // update molalities
  virtual void Update_kludge(const std::vector<Species>primary_species);
  // add stoichiometric contribution of complex to total
  virtual void AddContributionToTotal(std::vector<double> &total) = 0;
  // add derivative of total with respect to free-ion to dtotal
  virtual void AddContributionToDTotal(const std::vector<Species> primary_species,
                                       Block *dtotal) = 0;

  void set_ncomp(const int in_ncomp) { this->ncomp_ = in_ncomp; };
  void ncomp(const int in_ncomp) { this->ncomp_ = in_ncomp; };
  int ncomp(void) const { return this->ncomp_; };

  std::vector<SpeciesName> species_names(void) const { return this->species_names_; };
  std::vector<SpeciesId> species_ids(void) const { return this->species_ids_; };
  std::vector<double> stoichiometry(void) const { return this->stoichiometry_; };

  virtual void Display(void) const;

 protected:
  
  double log_to_ln(double d) { return d*2.30258509299; }
  double ln_to_log(double d) { return d*0.434294481904; }

  int ncomp_; // # components in reaction
  std::vector<SpeciesName> species_names_;
  std::vector<SpeciesId> species_ids_;       // ids of primary species in rxn
  std::vector<double> stoichiometry_;  // stoich of primary species in rxn
  std::vector<double> logK_array_;     // for temperature dep. logK
  double h2o_stoich_;                  // stoichiometry of water in equation
  double lnK_;                         // log value of equlibrium constant
  double lnQK_;                        // store lnQK for derivatives later
  double logK_;

 private:
};

#endif

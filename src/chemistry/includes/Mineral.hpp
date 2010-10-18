/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __Mineral_hpp__
#define __Mineral_hpp__

#include <string>
#include <vector>
#include <cmath>

#include "Species.hpp"
#include "SecondarySpecies.hpp"
#include "Block.hpp"

// Class for aqueous equilibrium complexation reaction

class Mineral : public SecondarySpecies {

 public:
  Mineral();
  Mineral(std::string s);
  Mineral(const SpeciesName name,
          SpeciesId mineral_id,
          std::vector<SpeciesName> species,
          std::vector<double> stoichiometries,
          std::vector<int> species_ids,
          const double h2o_stoich, const double mol_wt,
          const double logK, const double molar_density);
  ~Mineral();

  // update molalities
  void Update(const std::vector<Species>primary_species);
  // add stoichiometric contribution of complex to total
  void AddContributionToTotal(std::vector<double> &total);
  // add derivative of total with respect to free-ion to dtotal
  void AddContributionToDTotal(const std::vector<Species> primary_species,
                               Block *dtotal);

  void Display(void) const;

  double molar_density(void) const { return this->molar_density_; }
  void molar_density(double d) { this->molar_density_ = d; }

 protected:

 private:

  double molar_density_;
  
};

#endif // __Mineral_hpp__

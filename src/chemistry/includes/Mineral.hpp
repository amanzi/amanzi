/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __Mineral_hpp__
#define __Mineral_hpp__

#include <string>
#include <vector>
#include <cmath>

#include "Species.hpp"
#include "SecondarySpecies.hpp"
#include "Block.hpp"

/* Class for mineral reaction, should be written with the mineral as
** the reactant:
**
**  Calcite = 1.0 Ca++ + 1.0 HCO3- -1.0 H+
**
*/

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
  void DisplayResultsHeader(void) const;
  void DisplayResults(void) const;

  double saturation_index(void) const { return std::log10(Q_over_K()); };  // SI = log10(Q/Keq)
  double Q_over_K(void) const { return std::exp(this->lnQK_); };

  double molar_density(void) const { return this->molar_density_; }
  void molar_density(double d) { this->molar_density_ = d; }

  double surface_area(void) const { return this->surface_area_; }
  void set_surface_area(const double d) { this->surface_area_ = d; }
  double volume_fraction(void) const { return this->volume_fraction_; }
  void set_volume_fraction(const double d) { this->volume_fraction_ = d; }

protected:

 private:
  double saturation_index_;  
  double molar_density_;     // [moles / m^3] 

  double surface_area_;      // [m^2 mineral surface / m^3 mineral]
  double volume_fraction_;   // [m^3 mineral / m^3 bulk]

};

#endif // __Mineral_hpp__

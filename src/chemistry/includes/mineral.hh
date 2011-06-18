/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_MINERAL_HH_
#define AMANZI_CHEMISTRY_MINERAL_HH_

/* Class for mineral reaction, should be written with the mineral as
** the reactant:
**
**  Calcite = 1.0 Ca++ + 1.0 HCO3- -1.0 H+
**
*/

#include <cmath>

#include <vector>

#include "species.hh"
#include "secondary_species.hh"
#include "verbosity.hh"

namespace amanzi {
namespace chemistry {

// forward declarations from chemistry
class Block;

class Mineral : public SecondarySpecies {
 public:
  Mineral();
  Mineral(const SpeciesName name,
          SpeciesId mineral_id,
          const std::vector<SpeciesName>& species,
          const std::vector<double>& stoichiometries,
          const std::vector<SpeciesId>& species_ids,
          const double h2o_stoich, const double mol_wt,
          const double logK, const double molar_volume,
          const double specific_surface_area);
  ~Mineral();

  // update molalities
  virtual void Update(const std::vector<Species>& primary_species);
  // add stoichiometric contribution of complex to total
  virtual void AddContributionToTotal(std::vector<double> *total);
  // add derivative of total with respect to free-ion to dtotal
  virtual void AddContributionToDTotal(const std::vector<Species>& primary_species,
                                       Block* dtotal);

  void Display(void) const;
  void DisplayResultsHeader(void) const;
  void DisplayResults(void) const;

  double Q_over_K(void) const {
    return std::exp(this->lnQK_);
  };
  double saturation_index(void) const {
    return std::log10(Q_over_K());
  };  // SI = log10(Q/Keq)
  double molar_volume(void) const {
    return this->molar_volume_;
  }
  double surface_area(void) const {
    return this->surface_area_;
  }
  double specific_surface_area(void) const {
    return this->specific_surface_area_;
  }
  double volume_fraction(void) const {
    return this->volume_fraction_;
  }

  void set_volume_fraction(const double d) {
    this->volume_fraction_ = d;
  }
  void UpdateSurfaceAreaFromVolumeFraction(const double total_volume);

  void set_verbosity(const Verbosity verbosity) {
    this->verbosity_ = verbosity;
  };
  Verbosity verbosity(void) const {
    return this->verbosity_;
  };

 protected:

 private:
  Verbosity verbosity_;
  double saturation_index_;
  double molar_volume_;     // [cm^3 / moles]
  double specific_surface_area_;  // [m^2/g]
  double surface_area_;      // [m^2 mineral surface / m^3 mineral]
  double volume_fraction_;   // [m^3 mineral / m^3 bulk]

  void set_surface_area(const double d) {
    this->surface_area_ = d;
  }
  // void molar_volume(double d) { this->molar_volume_ = d; }
  // void set_specific_surface_area(const double d) { this->specific_surface_area_ = d; }
};

}  // namespace chemistry
}  // namespace amanzi 
#endif  // AMANZI_CHEMISTRY_MINERAL_HH_

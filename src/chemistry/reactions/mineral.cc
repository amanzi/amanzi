/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "mineral.hh"

#include <iostream>
#include <iomanip>

#include "secondary_species.hh"
#include "verbosity.hh"

Mineral::Mineral()
    : SecondarySpecies(),
      verbosity_(kSilent),
      molar_volume_(0.0),
      specific_surface_area_(0.0),
      surface_area_(0.0),
      volume_fraction_(0.0) {
}  // end Mineral() constructor

Mineral::Mineral(const SpeciesName in_name,
                 const SpeciesId in_id,
                 const std::vector<SpeciesName>& in_species,
                 const std::vector<double>& in_stoichiometries,
                 const std::vector<int>& in_species_ids,
                 const double in_h2o_stoich,
                 const double in_mol_wt,
                 const double in_logK,
                 const double molar_volume,
                 const double specific_surface_area)
    : SecondarySpecies(in_name, in_id,
                       in_species, in_stoichiometries, in_species_ids,
                       in_h2o_stoich, 0., in_mol_wt, 0., in_logK),
      verbosity_(kSilent),
      molar_volume_(molar_volume),
      specific_surface_area_(specific_surface_area),
      surface_area_(0.0),
      volume_fraction_(0.0) {
}  // end Mineral costructor


Mineral::~Mineral() {
}  // end Mineral() destructor

/*
**
**  these functions are only needed if mineral equilibrium is added.
**
*/

void Mineral::UpdateSurfaceAreaFromVolumeFraction(const double total_volume) {
  // area = SSA * GMW * V_fraction * V_total * unit_conversion / mole volume
  // m^2 = (m^2/g * g/mol * -- * m^3 * cm^3/m^3) / (cm^3/mol)
  // double cm3_in_m3 = 10.0;
  double cm3_in_m3 = 1.0e6;
  set_surface_area(specific_surface_area() * gram_molecular_weight() *
                   volume_fraction() * total_volume * cm3_in_m3 / molar_volume());
  // hard code bulk surface area: 100 m^2/m^3
  set_surface_area(100.0 * total_volume);

  if (verbosity() == kDebugMineralKinetics) {
    std::cout << "Mineral: " << name() << "\n"
              << "   SSA: " << specific_surface_area() << " [m^2/g]\n"
              << "   GMW:" << gram_molecular_weight() << " [g/mole]\n"
              << "   mole volume^-1: " << 1.0 / molar_volume() << " [mole/cm^3]\n"
              << "   volume fraction: " << volume_fraction() << " [-]\n"
              << std::scientific << "   total volume: " << total_volume << " [m^3]\n"
              << "   cm3 in m3=" << cm3_in_m3 << "\n"
              << "   surface area=" << surface_area() << " [m^2]\n"
              << std::fixed
              << std::endl;
  }
}  // end UpdateSurfaceAreaFromVolumeFraction()

void Mineral::Update(const std::vector<Species>& primary_species) {
  double lnQK = -lnK_;
  for (int i = 0; i < ncomp(); i++) {
    lnQK += stoichiometry_.at(i) * primary_species.at(species_ids_.at(i)).ln_activity();
  }
  lnQK_ = lnQK;
}  // end update()

void Mineral::AddContributionToTotal(std::vector<double> *total) {
  static_cast<void>(total);
}  // end addContributionToTotal()

void Mineral::AddContributionToDTotal(const std::vector<Species>& primary_species,
                                      Block* dtotal) {
  static_cast<void>(primary_species);
  static_cast<void>(dtotal);
}  // end addContributionToDTotal()


/*
**
**  Display functions
**
*/
void Mineral::Display(void) const {
  std::cout << "    " << name() << " = ";
  for (unsigned int i = 0; i < species_names_.size(); i++) {
    std::cout << std::setprecision(2) << stoichiometry_[i] << " " << species_names_[i];
    if (i < species_names_.size() - 1) {
      std::cout << " + ";
    }
  }
  std::cout << std::endl;
  std::cout << std::setw(40) << " "
            << std::setw(10) << std::setprecision(5) << std::fixed << logK_
            << std::setw(13) << molar_volume()
            << std::setw(13) << gram_molecular_weight()
            << std::setw(13) << specific_surface_area()
            << std::setw(13) << std::scientific << surface_area()
            << std::setw(13) << std::fixed << volume_fraction()
            << std::endl;
}  // end Display()

void Mineral::DisplayResultsHeader(void) const {
  std::cout << std::setw(15) << "Name"
            << std::setw(15) << "Q/K"
            << std::setw(15) << "SI"
            << std::endl;
}  // end DisplayResultsHeader()

void Mineral::DisplayResults(void) const {
  std::cout << std::setw(15) << name()
            << std::scientific << std::setprecision(5)
            << std::setw(15) << Q_over_K()
            << std::fixed << std::setprecision(3)
            << std::setw(15) << saturation_index()
            << std::endl;
}  // end DisplayResults()

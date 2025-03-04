/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry

  Supporting class for Beaker.
*/

#ifndef AMANZI_CHEMISTRY_BEAKER_STATE_HH_
#define AMANZI_CHEMISTRY_BEAKER_STATE_HH_

#include <vector>

namespace Amanzi {
namespace AmanziChemistry {

struct BeakerState {
  BeakerState()
    : temperature(298.15), // 25 C
      porosity(1.0),
      saturation(1.0),
      water_density(1000.0),
      volume(1.0){};

  std::vector<double> total; // molarity
  std::vector<double> total_sorbed;
  std::vector<double> free_ion; // molality
  std::vector<double> primary_activity_coeff;
  std::vector<double> secondary_activity_coeff;

  std::vector<double> mineral_volume_fraction;       // volume fractions
  std::vector<double> mineral_specific_surface_area; // [m^2 mineral/ m^3 bulk]

  std::vector<double> ion_exchange_sites;           // CEC
  std::vector<double> ion_exchange_ref_cation_conc; // [?]

  std::vector<double> surface_site_density;
  std::vector<double> surface_complex_free_site_conc; // [moles sites / m^3 bulk]

  std::vector<double> isotherm_kd;
  std::vector<double> isotherm_freundlich_n;
  std::vector<double> isotherm_langmuir_b;

  std::vector<double> total_sorbed_colloid_mobile; // primaries on colloids
  std::vector<double> total_sorbed_colloid_immobile;

  std::vector<double> colloid_mobile_conc; // colloids fraction
  std::vector<double> colloid_immobile_conc;

  double temperature;   // [K]
  double porosity;      // [-]
  double saturation;    // [-]
  double water_density; // [kg/m^3]
  double volume;        // [m^3]
};

} // namespace AmanziChemistry
} // namespace Amanzi

#endif

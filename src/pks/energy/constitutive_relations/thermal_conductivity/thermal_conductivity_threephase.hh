/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Base class of a three-phase Thermal Conductivity relation.

UNITS: ????
------------------------------------------------------------------------- */

#ifndef PK_ENERGY_RELATIONS_TC_THREEPHASE_HH_
#define PK_ENERGY_RELATIONS_TC_THREEPHASE_HH_

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

class ThermalConductivityThreePhase {

public:
  virtual double ThermalConductivity(double porosity, double sat_liq, double sat_ice) = 0;
};

} // namespace
} // namespace
} // namespace

#endif

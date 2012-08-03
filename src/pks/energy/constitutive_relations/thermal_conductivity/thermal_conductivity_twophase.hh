/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Base class of a two-phase Thermal Conductivity relation.

UNITS: ????
------------------------------------------------------------------------- */

#ifndef PK_ENERGY_RELATIONS_TC_TWOPHASE_HH_
#define PK_ENERGY_RELATIONS_TC_TWOPHASE_HH_

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

class ThermalConductivityTwoPhase {

public:
  virtual double ThermalConductivity(double porosity, double sat_liq) = 0;
};

} // namespace
} // namespace
} // namespace

#endif

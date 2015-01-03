/*
  This is the energy component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Base class of a two-phase Thermal Conductivity relation.
*/

#ifndef PK_ENERGY_RELATIONS_TC_TWOPHASE_HH_
#define PK_ENERGY_RELATIONS_TC_TWOPHASE_HH_

namespace Amanzi {
namespace Energy {

class ThermalConductivityTwoPhase {
 public:
  virtual ~ThermalConductivityTwoPhase() {}
  virtual double ThermalConductivity(double porosity, double sat_liq) = 0;
};

}  // namespace Energy
}  // namespace Amanzi

#endif

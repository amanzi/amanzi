/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Base class of a two-phase thermal conductivity model.
*/

#ifndef PK_ENERGY_TCM_TWOPHASE_HH_
#define PK_ENERGY_TCM_TWOPHASE_HH_

namespace Amanzi {
namespace Energy {

class TCM_TwoPhase {
 public:
  virtual ~TCM_TwoPhase(){};
  virtual double ThermalConductivity(double porosity, double sat_liq) = 0;
};

} // namespace Energy
} // namespace Amanzi

#endif

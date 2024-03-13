/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Energy

  Simple model of two-phase thermal conductivity, based upon:
   - Interpolation between saturated and dry conductivities via a Kersten number.
   - Power-law Kersten number.
   - Emperical fit for dry conductivity from Peters-Lidard et al '98.

  See TCM_PetersLidard_TwoPhase.hh for more detail.
*/

#ifndef AMANZI_ENERGY_TCM_PETERSLIDARD_TWOPHASE_HH_
#define AMANZI_ENERGY_TCM_PETERSLIDARD_TWOPHASE_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "TCM_TwoPhase.hh"

namespace Amanzi {
namespace Energy {

class TCM_PetersLidard_TwoPhase : public TCM_TwoPhase {
 public:
  TCM_PetersLidard_TwoPhase(Teuchos::ParameterList& plist);

  double ThermalConductivity(double porosity, double sat_liq);

 private:
  void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double eps_;
  double alpha_;
  double k_rock_;
  double k_liquid_;
  double k_gas_;
  double d_;

 private:
  static Utils::RegisteredFactory<TCM_TwoPhase, TCM_PetersLidard_TwoPhase> reg_;
};

} // namespace Energy
} // namespace Amanzi

#endif

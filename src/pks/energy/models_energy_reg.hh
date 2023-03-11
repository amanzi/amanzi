/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Energy

  Self-registering factory for IEM and TCM models.
*/

#include "IEMFactory.hh"
#include "IEMEvaluator.hh"
#include "IEM_Linear.hh"
#include "IEM_Tabular.hh"
#include "IEM_WaterVaporEvaluator.hh"

#include "TCMFactory_TwoPhase.hh"
#include "TCM_PetersLidard_TwoPhase.hh"
#include "TCM_WetDry_TwoPhase.hh"

namespace Amanzi {
namespace Utils {

// explicity instantitate the static factory data
template <>
Amanzi::Utils::Factory<Amanzi::Energy::IEM>::map_type*
  Amanzi::Utils::Factory<Amanzi::Energy::IEM>::map_;

template <>
Factory<Energy::TCM_TwoPhase>::map_type* Factory<Energy::TCM_TwoPhase>::map_;

} // namespace Utils
} // namespace Amanzi


namespace Amanzi {
namespace Energy {

Utils::RegisteredFactory<Evaluator, IEMEvaluator> IEMEvaluator::factory_("iem");
Utils::RegisteredFactory<Evaluator, IEM_WaterVaporEvaluator>
  IEM_WaterVaporEvaluator::factory_("iem water vapor");

Utils::RegisteredFactory<IEM, IEM_Linear> IEM_Linear::factory_("linear");
Utils::RegisteredFactory<IEM, IEM_Tabular> IEM_Tabular::factory_("tabular");


// linear interpolant of thermal conductivity.
Utils::RegisteredFactory<TCM_TwoPhase, TCM_PetersLidard_TwoPhase>
  TCM_PetersLidard_TwoPhase::factory_("two-phase Peters-Lidard");

// simple model of two-phase thermal conductivity, based upon:
// - Interpolation between saturated and dry conductivities via a Kersten number.
// - Power-law Kersten number.
Utils::RegisteredFactory<TCM_TwoPhase, TCM_WetDry_TwoPhase>
  TCM_WetDry_TwoPhase::factory_("two-phase wet/dry");

} // namespace Energy
} // namespace Amanzi


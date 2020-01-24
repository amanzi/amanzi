/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Self-registering factory for IEM implementations.
*/

#include "IEMFactory.hh"
#include "IEMEvaluator.hh"
#include "IEM_Linear.hh"
#include "IEM_WaterVaporEvaluator.hh"

// explicity instantitate the static data of Factory<IEM>
template<> 
Amanzi::Utils::Factory<Amanzi::Energy::IEM>::map_type* 
   Amanzi::Utils::Factory<Amanzi::Energy::IEM>::map_;

namespace Amanzi {
namespace Energy {

Utils::RegisteredFactory<FieldEvaluator,IEMEvaluator> IEMEvaluator::factory_("iem");

Utils::RegisteredFactory<IEM,IEM_Linear> IEM_Linear::factory_("linear");
Utils::RegisteredFactory<FieldEvaluator,IEM_WaterVaporEvaluator> IEM_WaterVaporEvaluator::factory_("iem water vapor");

}  // namespace Energy
}  // namespace Amanzi

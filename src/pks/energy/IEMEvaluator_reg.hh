/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  The internal energy model evaluator simply calls the IEM with
  the correct arguments.
*/

#include "IEMEvaluator.hh"

namespace Amanzi {
namespace Energy {

Utils::RegisteredFactory<FieldEvaluator,IEMEvaluator> IEMEvaluator::factory_("iem");

}  // namespace Energy
}  // namespace Amanzi

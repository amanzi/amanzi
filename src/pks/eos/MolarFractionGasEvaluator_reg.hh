/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Determining the molar fraction of a gas component within a gas mixture.
*/

#include "MolarFractionGasEvaluator.hh"
#include "VaporPressureBaseFactory.hh"

namespace Amanzi {
namespace AmanziEOS {

// registry of method
Utils::RegisteredFactory<FieldEvaluator, MolarFractionGasEvaluator>
    MolarFractionGasEvaluator::factory_("molar fraction gas");

}  // namespace AmanziEOS
}  // namespace Amanzi


/*
  This is the energy component of the ATS and Amanzi codes. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  The WRM Evaluator simply calls the WRM with the correct arguments.
*/

#include "iem_evaluator.hh"

namespace Amanzi {
namespace Energy {

Utils::RegisteredFactory<FieldEvaluator,IEMEvaluator> IEMEvaluator::factory_("iem");


}  // namespace Energy
}  // namespace Amanzi

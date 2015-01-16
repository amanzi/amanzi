/*
  This is the EOS component of the ATS and Amanzi codes.
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  EOSFieldEvaluator is the interface between state/data and the model, an EOS.
*/

#include "eos_evaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,EOSEvaluator> EOSEvaluator::factory_("eos");

}  // namespace Relations
}  // namespace Amanzi

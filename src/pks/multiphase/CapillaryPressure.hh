/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_MULTIPHASE_CAPILLARY_PRESSURE_HH_
#define AMANZI_MULTIPHASE_CAPILLARY_PRESSURE_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "secondary_variable_field_evaluator.hh"

// Multiphase
#include "WRMmp.hh"
#include "MultiphaseTypeDefs.hh"

namespace Amanzi {
namespace Multiphase {

class CapillaryPressure : public SecondaryVariableFieldEvaluator {
 public:
  CapillaryPressure(Teuchos::ParameterList& plist, Teuchos::RCP<WRMmpPartition> wrm);

  // interface methods from FieldEvaluator
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 private:
  Teuchos::RCP<WRMmpPartition> wrm_;
};

}  // namespace Flow
}  // namespace Amanzi

#endif


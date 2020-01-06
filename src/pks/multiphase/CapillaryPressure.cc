/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "errors.hh"

#include "ModelMeshPartition.hh"
#include "MultiphaseDefs.hh"
#include "CapillaryPressure.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Simple constructor
****************************************************************** */
CapillaryPressure::CapillaryPressure(Teuchos::ParameterList& plist,
                                     Teuchos::RCP<WRMmpPartition> wrm) :
    SecondaryVariableFieldEvaluator(plist),
    wrm_(wrm)
{
  my_key_ = "capillary_pressure";
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> CapillaryPressure::Clone() const {
  return Teuchos::rcp(new CapillaryPressure(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void CapillaryPressure::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& sat_c = *S->GetFieldData("saturation_liquid")->ViewComponent("cell", false);
  auto& result_c = *result->ViewComponent("cell", false);

  int ncells = result_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = wrm_->second[(*wrm_->first)[c]]->capillaryPressure(sat_c[0][c]);
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void CapillaryPressure::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& sat_c = *S->GetFieldData("saturation_liquid")->ViewComponent("cell", false);
  auto& result_c = *result->ViewComponent("cell", false);

  int ncells = result_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    sat_c[0][c] = wrm_->second[(*wrm_->first)[c]]->dPc_dS(sat_c[0][c]);
  }
}

}  // namespace Multiphase
}  // namespace Amanzi


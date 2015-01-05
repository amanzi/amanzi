/*
  This is the energy component of the ATS and Amanzi codes. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Interface for a thermal conductivity model with two phases.
*/

#include "dbc.hh"
#include "twophase_thermal_conductivity_factory.hh"
#include "twophase_thermal_conductivity_evaluator.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Constructor.
****************************************************************** */
ThermalConductivityTwoPhaseEvaluator::ThermalConductivityTwoPhaseEvaluator(
    Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("thermal conductivity key", "thermal_conductivity");
  }

  poro_key_ = plist_.get<std::string>("porosity key", "porosity");
  dependencies_.insert(poro_key_);

  sat_key_ = plist_.get<std::string>("saturation key", "saturation_liquid");
  dependencies_.insert(sat_key_);

  ASSERT(plist_.isSublist("thermal conductivity parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("thermal conductivity parameters");
  ThermalConductivityTwoPhaseFactory fac;
  tc_ = fac.createThermalConductivityModel(sublist);
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
ThermalConductivityTwoPhaseEvaluator::ThermalConductivityTwoPhaseEvaluator(
    const ThermalConductivityTwoPhaseEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    poro_key_(other.poro_key_),
    sat_key_(other.sat_key_),
    tc_(other.tc_) {};


/* ******************************************************************
* TBW.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> ThermalConductivityTwoPhaseEvaluator::Clone() const
{
  return Teuchos::rcp(new ThermalConductivityTwoPhaseEvaluator(*this));
}


/* ******************************************************************
* Evaluator body.
****************************************************************** */
void ThermalConductivityTwoPhaseEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S, const Teuchos::Ptr<CompositeVector>& result)
{
  // pull out the dependencies
  Teuchos::RCP<const CompositeVector> poro = S->GetFieldData(poro_key_);
  Teuchos::RCP<const CompositeVector> sat = S->GetFieldData(sat_key_);

  for (CompositeVector::name_iterator comp=result->begin(); comp!=result->end(); ++comp) {
    // much more efficient to pull out vectors first
    const Epetra_MultiVector& poro_v = *poro->ViewComponent(*comp,false);
    const Epetra_MultiVector& sat_v = *sat->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = tc_->ThermalConductivity(poro_v[0][i], sat_v[0][i]);
    }
  }
}


/* ******************************************************************
* Evaluator of derivarives.
****************************************************************** */
void ThermalConductivityTwoPhaseEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S, Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  ASSERT(0); // not implemented, not yet needed
}

}  // namespace Energy
}  // namespace Amanzi

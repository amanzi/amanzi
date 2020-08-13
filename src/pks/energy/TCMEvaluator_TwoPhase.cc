/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Interface for a thermal conductivity model with two phases.
*/

#include "dbc.hh"
#include "TCMFactory_TwoPhase.hh"
#include "TCMEvaluator_TwoPhase.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Constructor.
****************************************************************** */
TCMEvaluator_TwoPhase::TCMEvaluator_TwoPhase(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  my_key_ = plist_.get<std::string>("thermal conductivity key");
  auto prefix = Keys::getDomainPrefix(my_key_);

  porosity_key_ = plist_.get<std::string>("porosity key", prefix + "porosity");
  dependencies_.insert(porosity_key_);

  saturation_key_ = plist_.get<std::string>("saturation key", prefix + "saturation_liquid");
  dependencies_.insert(saturation_key_);

  AMANZI_ASSERT(plist_.isSublist("thermal conductivity parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("thermal conductivity parameters");
  TCMFactory_TwoPhase fac;
  tc_ = fac.CreateTCM(sublist);
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
TCMEvaluator_TwoPhase::TCMEvaluator_TwoPhase(const TCMEvaluator_TwoPhase& other) :
    SecondaryVariableFieldEvaluator(other),
    porosity_key_(other.porosity_key_),
    saturation_key_(other.saturation_key_),
    tc_(other.tc_) {};


/* ******************************************************************
* TBW.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> TCMEvaluator_TwoPhase::Clone() const {
  return Teuchos::rcp(new TCMEvaluator_TwoPhase(*this));
}


/* ******************************************************************
* Evaluator body.
****************************************************************** */
void TCMEvaluator_TwoPhase::EvaluateField_(
    const Teuchos::Ptr<State>& S, const Teuchos::Ptr<CompositeVector>& result)
{
  // pull out the dependencies
  Teuchos::RCP<const CompositeVector> poro = S->GetFieldData(porosity_key_);
  Teuchos::RCP<const CompositeVector> sat = S->GetFieldData(saturation_key_);

  for (auto comp = result->begin(); comp != result->end(); ++comp) {
    // much more efficient to pull out vectors first
    const Epetra_MultiVector& poro_v = *poro->ViewComponent(*comp);
    const Epetra_MultiVector& sat_v = *sat->ViewComponent(*comp);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp);

    int ncomp = result->size(*comp);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] = tc_->ThermalConductivity(poro_v[0][i], sat_v[0][i]);
    }
  }
}


/* ******************************************************************
* Evaluator of derivarives.
****************************************************************** */
void TCMEvaluator_TwoPhase::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  AMANZI_ASSERT(0);
}

}  // namespace Energy
}  // namespace Amanzi

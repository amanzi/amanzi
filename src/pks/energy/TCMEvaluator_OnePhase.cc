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
#include "TCMEvaluator_OnePhase.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Constructor.
****************************************************************** */
TCMEvaluator_OnePhase::TCMEvaluator_OnePhase(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  my_key_ = plist_.get<std::string>("thermal conductivity key");
  auto prefix = Keys::getDomainPrefix(my_key_);

  temperature_key_ = plist_.get<std::string>("temperature key", prefix + "temperature");
  dependencies_.insert(temperature_key_);

  porosity_key_ = plist_.get<std::string>("porosity key", prefix + "porosity");
  dependencies_.insert(porosity_key_);

  AMANZI_ASSERT(plist_.isSublist("thermal conductivity parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("thermal conductivity parameters");
  tc_ = Teuchos::rcp(new AmanziEOS::ThermalConductivity_Water(sublist));

  k_rock_ = sublist.get<double>("thermal conductivity of rock");
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
TCMEvaluator_OnePhase::TCMEvaluator_OnePhase(const TCMEvaluator_OnePhase& other) :
    SecondaryVariableFieldEvaluator(other),
    temperature_key_(other.temperature_key_),
    tc_(other.tc_) {};


/* ******************************************************************
* TBW.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> TCMEvaluator_OnePhase::Clone() const {
  return Teuchos::rcp(new TCMEvaluator_OnePhase(*this));
}


/* ******************************************************************
* Evaluator body.
****************************************************************** */
void TCMEvaluator_OnePhase::EvaluateField_(
    const Teuchos::Ptr<State>& S, const Teuchos::Ptr<CompositeVector>& result)
{
  // pull out the dependencies
  const Epetra_MultiVector& temp_c = *S->GetFieldData(temperature_key_)->ViewComponent("cell");
  const Epetra_MultiVector& poro_c = *S->GetFieldData(porosity_key_)->ViewComponent("cell");
  Epetra_MultiVector& result_c = *result->ViewComponent("cell");

  int ncomp = result->size("cell", false);
  for (int i = 0; i != ncomp; ++i) {
    double phi = poro_c[0][i];
    double k_liq = tc_->ThermalConductivity(temp_c[0][i]);
    result_c[0][i] = phi * k_liq + (1.0 - phi) * k_rock_;
  }
}


/* ******************************************************************
* Evaluator of derivarives.
****************************************************************** */
void TCMEvaluator_OnePhase::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  AMANZI_ASSERT(0);
}

}  // namespace Energy
}  // namespace Amanzi

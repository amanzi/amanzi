/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Field evaluator for enthalpy, H = U + p / rho. 
*/

#include "EnthalpyEvaluator.hh"
#include "IEMEvaluator.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Constructor of a secondary evaluator.
****************************************************************** */
EnthalpyEvaluator::EnthalpyEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("enthalpy key");
  }
  auto prefix = Keys::getDomainPrefix(my_key_);  // include dash

  // Set up my dependencies.
  // -- internal energy
  ie_liquid_key_ = plist_.get<std::string>("internal energy key", prefix + "internal_energy_liquid");
  dependencies_.insert(ie_liquid_key_);

  // -- pressure work
  include_work_ = plist_.get<bool>("include work term", true);

  if (include_work_) {
    pressure_key_ = plist_.get<std::string>("pressure key", prefix + "pressure");
    dependencies_.insert(pressure_key_);

    mol_density_liquid_key_ = plist_.get<std::string>("molar density key", prefix + "molar_density_liquid");
    dependencies_.insert(mol_density_liquid_key_);
  }
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
EnthalpyEvaluator::EnthalpyEvaluator(const EnthalpyEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    pressure_key_(other.pressure_key_),
    mol_density_liquid_key_(other.mol_density_liquid_key_),
    ie_liquid_key_(other.ie_liquid_key_),
    include_work_(other.include_work_) {};


/* ******************************************************************
* TBW.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> EnthalpyEvaluator::Clone() const {
  return Teuchos::rcp(new EnthalpyEvaluator(*this));
}


/* ******************************************************************
* Evaluator's body.
****************************************************************** */
void EnthalpyEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> u_l = S->GetFieldData(ie_liquid_key_);
  *result = *u_l;

  if (include_work_) {
    Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pressure_key_);
    Teuchos::RCP<const CompositeVector> n_l = S->GetFieldData(mol_density_liquid_key_);

    for (auto comp = result->begin(); comp != result->end(); ++comp) {
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp);
      const Epetra_MultiVector& nl_v = *n_l->ViewComponent(*comp);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp);

      int ncomp = result->size(*comp);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] += pres_v[0][i] / nl_v[0][i];
      }
    }
  }
}


/* ******************************************************************
* Evaluator for derivatives.
****************************************************************** */
void EnthalpyEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  if (wrt_key == ie_liquid_key_) {
    result->PutScalar(1.0);

  } else if (wrt_key == pressure_key_) {
    AMANZI_ASSERT(include_work_);
    
    Teuchos::RCP<const CompositeVector> n_l = S->GetFieldData(mol_density_liquid_key_);

    for (CompositeVector::name_iterator comp = result->begin(); comp != result->end(); ++comp) {
      const Epetra_MultiVector& nl_v = *n_l->ViewComponent(*comp);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp);

      int ncomp = result->size(*comp);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = 1.0 / nl_v[0][i];
      }
    }

  } else if (wrt_key == mol_density_liquid_key_) {
    AMANZI_ASSERT(include_work_);
    
    Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pressure_key_);
    Teuchos::RCP<const CompositeVector> n_l = S->GetFieldData(mol_density_liquid_key_);

    for (auto comp = result->begin(); comp != result->end(); ++comp) {
      const Epetra_MultiVector& nl_v = *n_l->ViewComponent(*comp);
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp);

      int ncomp = result->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = -pres_v[0][i] / (nl_v[0][i] * nl_v[0][i]);
      }
    }
  }
}


/* ******************************************************************
* Evaluation at a point
****************************************************************** */
double EnthalpyEvaluator::EvaluateFieldSingle(
    const Teuchos::Ptr<State>& S, int c, double T)
{
  double tmp = Teuchos::rcp_dynamic_cast<IEMEvaluator>(S->GetFieldEvaluator(ie_liquid_key_))->EvaluateFieldSingle(c, T);
  if (include_work_) {
    const auto& pres_c = *S->GetFieldData(pressure_key_)->ViewComponent("cell");
    const auto& nl_c = *S->GetFieldData(mol_density_liquid_key_)->ViewComponent("cell");

    tmp += pres_c[0][c] / nl_c[0][c];
  }
  return tmp;
}

}  // namespace Energy
}  // namespace Amanzi

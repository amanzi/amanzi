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

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Constructor of a secondary evaluator.
****************************************************************** */
EnthalpyEvaluator::EnthalpyEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("enthalpy key", "enthalpy_liquid");
  }

  // Set up my dependencies.
  // -- internal energy
  ie_key_ = plist_.get<std::string>("internal energy key", "internal_energy_liquid");
  dependencies_.insert(ie_key_);

  // -- pressure work
  include_work_ = plist_.get<bool>("include work term", true);

  if (include_work_) {
    pres_key_ = plist_.get<std::string>("pressure key", "pressure");
    dependencies_.insert(pres_key_);

    dens_key_ = plist_.get<std::string>("molar density key", "molar_density_liquid");
    dependencies_.insert(dens_key_);
  }
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
EnthalpyEvaluator::EnthalpyEvaluator(const EnthalpyEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    pres_key_(other.pres_key_),
    dens_key_(other.dens_key_),
    ie_key_(other.ie_key_),
    include_work_(other.include_work_) {}


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
  Teuchos::OSTab tab = vo_->getOSTab();
  Teuchos::RCP<const CompositeVector> u_l = S->GetFieldData(ie_key_);
  *result = *u_l;

  if (include_work_) {
    Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);
    Teuchos::RCP<const CompositeVector> n_l = S->GetFieldData(dens_key_);

    for (CompositeVector::name_iterator comp=result->begin(); comp!=result->end(); ++comp) {
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *n_l->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp, false);

      int ncomp = result->size(*comp, false);
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
  if (wrt_key == ie_key_) {
    result->PutScalar(1.0);

  } else if (wrt_key == pres_key_) {
    AMANZI_ASSERT(include_work_);
    
    Teuchos::RCP<const CompositeVector> n_l = S->GetFieldData(dens_key_);

    for (CompositeVector::name_iterator comp = result->begin(); comp != result->end(); ++comp) {
      const Epetra_MultiVector& nl_v = *n_l->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp, false);

      int ncomp = result->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = 1.0 / nl_v[0][i];
      }
    }

  } else if (wrt_key == dens_key_) {
    AMANZI_ASSERT(include_work_);
    
    Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);
    Teuchos::RCP<const CompositeVector> n_l = S->GetFieldData(dens_key_);

    for (CompositeVector::name_iterator comp = result->begin(); comp != result->end(); ++comp) {
      const Epetra_MultiVector& nl_v = *n_l->ViewComponent(*comp, false);
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp, false);

      int ncomp = result->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = -pres_v[0][i] / (nl_v[0][i] * nl_v[0][i]);
      }
    }
  }
}

}  // namespace Energy
}  // namespace Amanzi

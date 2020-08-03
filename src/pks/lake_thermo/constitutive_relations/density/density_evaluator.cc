/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Svetlana Tokareva (tokareva@lanl.gov)

FieldEvaluator for water density.
----------------------------------------------------------------------------- */


#include "density_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

DensityEvaluator::DensityEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  if (my_key_.empty()) {

    my_key_ = plist_.get<std::string>("enthalpy key", "surface-enthalpy_liquid");
  }

  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_key_);
  include_work_ = plist_.get<bool>("include work term", true);

  // -- pressure
  if (include_work_) {
    pres_key_ = Keys::readKey(plist_, domain_name, "pressure", "pressure");
    dependencies_.insert(pres_key_);

    dens_key_ = Keys::readKey(plist_, domain_name, "molar density liquid", "molar_density_liquid");
    dependencies_.insert(dens_key_);
  }

  ie_key_ = Keys::readKey(plist_, domain_name, "internal energy liquid", "internal_energy_liquid");
  dependencies_.insert(ie_key_);

};

DensityEvaluator::DensityEvaluator(const DensityEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    pres_key_(other.pres_key_),
    dens_key_(other.dens_key_),
    ie_key_(other.ie_key_),
    include_work_(other.include_work_) {};

Teuchos::RCP<FieldEvaluator>
DensityEvaluator::Clone() const {
  return Teuchos::rcp(new DensityEvaluator(*this));
};


void DensityEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temperature_key_);

  double rho0 = 1000.;

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      double T = temp_v[0][i];
      result_v[0][i] = rho0 * (1.+8.0*1.e-5 + 5.88*1.e-5*T - 8.11*1.e-6*T*T + 4.77*1.e-8*T*T*T);
    }
  }
};


void DensityEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  // not implemented
  if (wrt_key == temperature_key_) {
    Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temperature_key_);

    double rho0 = 1000.;

    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        double T = temp_v[0][i];
        result_v[0][i] = rho0 * (5.88*1.e-5 - 2.*8.11*1.e-6*T + 3.*4.77*1.e-8*T*T);
      }
    }
  }
};


} //namespace
} //namespace

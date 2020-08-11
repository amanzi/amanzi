/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  The internal energy model evaluator simply calls the IEM with 
  the correct arguments.
*/

#include "IEM_WaterVaporEvaluator.hh"

namespace Amanzi {
namespace Energy {

IEM_WaterVaporEvaluator::IEM_WaterVaporEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  // defaults work fine, this sublist need not exist
  Teuchos::ParameterList sublist = plist.sublist("IEM parameters");
  iem_ = Teuchos::rcp(new IEM_WaterVapor(sublist));

  InitializeFromPlist_();
}


IEM_WaterVaporEvaluator::IEM_WaterVaporEvaluator(Teuchos::ParameterList& plist,
                                                 const Teuchos::RCP<IEM_WaterVapor>& iem) :
    SecondaryVariableFieldEvaluator(plist),
    iem_(iem) {
  InitializeFromPlist_();
}


IEM_WaterVaporEvaluator::IEM_WaterVaporEvaluator(const IEM_WaterVaporEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    iem_(other.iem_),
    temp_key_(other.temp_key_),
    mol_frac_key_(other.mol_frac_key_) {};


Teuchos::RCP<FieldEvaluator>
IEM_WaterVaporEvaluator::Clone() const {
  return Teuchos::rcp(new IEM_WaterVaporEvaluator(*this));
}


void IEM_WaterVaporEvaluator::InitializeFromPlist_()
{
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("internal energy key");
  }

  // Set up my dependencies.
  std::string domain = Keys::getDomain(my_key_);
  temp_key_ = plist_.get<std::string>("temperature key", Keys::getKey(domain, "temperature"));
  dependencies_.insert(temp_key_);

  // -- molar fraction of water vapor in the gaseous phase
  mol_frac_key_ = plist_.get<std::string>("vapor molar fraction key", Keys::getKey(domain, "molar_fraction_gas"));
  dependencies_.insert(mol_frac_key_);
}


void IEM_WaterVaporEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                                             const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const CompositeVector> mol_frac = S->GetFieldData(mol_frac_key_);

  for (CompositeVector::name_iterator comp = result->begin(); comp != result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp, false);
    const Epetra_MultiVector& molfrac_v = *mol_frac->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp, false);

    int ncomp = result->size(*comp, false);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] = iem_->InternalEnergy(temp_v[0][i], molfrac_v[0][i]);
    }
  }
}


void IEM_WaterVaporEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const CompositeVector> mol_frac = S->GetFieldData(mol_frac_key_);

  if (wrt_key == temp_key_) {
    for (CompositeVector::name_iterator comp = result->begin(); comp != result->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp, false);
      const Epetra_MultiVector& molfrac_v = *mol_frac->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp, false);

      int ncomp = result->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = iem_->DInternalEnergyDT(temp_v[0][i], molfrac_v[0][i]);
      }
    }
  } else if (wrt_key == mol_frac_key_) {
    for (CompositeVector::name_iterator comp = result->begin(); comp != result->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp, false);
      const Epetra_MultiVector& molfrac_v = *mol_frac->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp, false);

      int ncomp = result->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = iem_->DInternalEnergyDomega(temp_v[0][i], molfrac_v[0][i]);
      }
    }
  } else {
    AMANZI_ASSERT(0);
  }
}

}  // namespace Energy
}  // namespace Amanzi

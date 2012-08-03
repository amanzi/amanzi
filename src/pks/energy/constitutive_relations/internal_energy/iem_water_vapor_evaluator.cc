/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  The WRM Evaluator simply calls the WRM with the correct arguments.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "iem_water_vapor_evaluator.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

IEMWaterVaporEvaluator::IEMWaterVaporEvaluator(Teuchos::ParameterList& iem_plist) :
    SecondaryVariableFieldEvaluator(),
    iem_plist_(iem_plist) {
  ASSERT(iem_plist.isSublist("IEM parameters"));
  Teuchos::ParameterList sublist = iem_plist.sublist("IEM parameters");
  iem_ = Teuchos::rcp(new IEMWaterVapor(sublist));

  InitializeFromPlist_();
}

IEMWaterVaporEvaluator::IEMWaterVaporEvaluator(Teuchos::ParameterList& iem_plist, const Teuchos::RCP<IEMWaterVapor>& iem) :
    SecondaryVariableFieldEvaluator(),
    iem_plist_(iem_plist),
    iem_(iem) {
  InitializeFromPlist_();
}

IEMWaterVaporEvaluator::IEMWaterVaporEvaluator(const IEMWaterVaporEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    iem_plist_(other.iem_plist_),
    iem_(other.iem_),
    temp_key_(other.temp_key_) {}

Teuchos::RCP<FieldEvaluator>
IEMWaterVaporEvaluator::Clone() const {
  return Teuchos::rcp(new IEMWaterVaporEvaluator(*this));
}

void IEMEvaluator::InitializeFromPlist_() {
  my_key_ = iem_plist_.get<std::string>("internal energy key");

  // Set up my dependencies.
  std::size_t end = a_key_.find_first_of("_");
  std::string domain_name = a_key_.substr(0,end);
  if (domain_name == std::string("internal") ||
      domain_name == std::string("energy")) {
    domain_name = std::string("");
  } else {
    domain_name = domain_name+std::string("_");
  }

  // -- temperature
  temp_key_ = eos_plist_.get<std::string>("temperature key",
          domain_name+std::string("temperature"));
  dependencies_.insert(temp_key_);

  // -- molar fraction of water vapor in the gaseous phase
  mol_frac_key_ = eos_plist_.get<std::string>("vapor molar fraction key",
          domain_name+std::string("mol_frac_gas"));
  dependencies_.insert(mol_frac_key_);
}


void IEMEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const CompositeVector> mol_frac = S->GetFieldData(mol_frac_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    for (int i=0; i!=poro->size(*comp); ++i) {
      (*result)(*comp, i) =
        iem_->InternalEnergy((*temp)(*comp, i), (*mol_frac)(*comp, i));
    }
  }
}


void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const CompositeVector> mol_frac = S->GetFieldData(mol_frac_key_);

  if (wrt_key == temp_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      for (int i=0; i!=poro->size(*comp); ++i) {
        (*result)(*comp, i) =
          iem_->DInternalEnergyDT((*temp)(*comp, i), (*mol_frac)(*comp, i));
      }
    }
  } else if (wrt_key == mol_frac_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      for (int i=0; i!=poro->size(*comp); ++i) {
        (*result)(*comp, i) =
          iem_->DInternalEnergyDomega((*temp)(*comp, i), (*mol_frac)(*comp, i));
      }
    }
  } else {
    ASSERT(0);
  }
}


} //namespace
} //namespace
} //namespace

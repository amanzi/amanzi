/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Rel perm( pc ( sat ) ).
*/

#include "FlowDefs.hh"
#include "RelPermEvaluator.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Two constructors.
****************************************************************** */
RelPermEvaluator::RelPermEvaluator(Teuchos::ParameterList& plist,
                                   const Teuchos::Ptr<State>& S,
                                   const Teuchos::RCP<WRMPartition>& wrm) :
    SecondaryVariableFieldEvaluator(plist),
    wrm_(wrm) {
  InitializeFromPlist_(S);
}

RelPermEvaluator::RelPermEvaluator(const RelPermEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    wrm_(other.wrm_),
    pressure_key_(other.pressure_key_),
    patm_(other.patm_) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> RelPermEvaluator::Clone() const {
  return Teuchos::rcp(new RelPermEvaluator(*this));
}


/* ******************************************************************
* Initialization.
****************************************************************** */
void RelPermEvaluator::InitializeFromPlist_(const Teuchos::Ptr<State>& S)
{
  // my keys is for rel perm.
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("relative permeability key");
  }

  // my dependency is pressure.
  std::string domain = Keys::getDomain(my_key_);
  pressure_key_ = plist_.get<std::string>("pressure key", Keys::getKey(domain, "pressure"));
  dependencies_.insert(pressure_key_);

  // use rel perm class for calcualtion
  Teuchos::ParameterList plist;
  relperm_ = Teuchos::rcp(new RelPerm(plist, S->GetMesh(domain), patm_, wrm_));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void RelPermEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  patm_ = *S->GetScalarData("atmospheric_pressure");

  // relperm_->Compute(S->GetFieldData(pressure_key_), result);
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void RelPermEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result)
{
  patm_ = *S->GetScalarData("atmospheric_pressure");
  // relperm_->ComputeDerivative(S->GetFieldData(pressure_key_), result);
}


/* ******************************************************************
* TBW
****************************************************************** */
void RelPermEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  // Create an unowned factory to check my dependencies.
  auto dep_fac = Teuchos::rcp(new CompositeVectorSpace());
  dep_fac->SetMesh(S->GetMesh(Keys::getDomain(my_key_)));
  dep_fac->SetOwned(false);
  dep_fac->SetComponent("cell", AmanziMesh::CELL, 1);

  // Loop over my dependencies, ensuring they meet the requirements.
  for (const auto& key : dependencies_) {
    if (key == my_key_) {
      Errors::Message msg;
      msg << "Evaluator for key \"" << my_key_ << "\" depends upon itself.";
      Exceptions::amanzi_throw(msg);
    }
    Teuchos::RCP<CompositeVectorSpace> fac = S->RequireField(key);
    fac->Update(*dep_fac);
  }

  // Recurse into the tree to propagate info to leaves.
  for (const auto& key : dependencies_) {
    S->RequireFieldEvaluator(key)->EnsureCompatibility(S);
  }
}

}  // namespace Flow
}  // namespace Amanzi

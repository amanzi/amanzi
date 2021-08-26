/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Relative permeability as a function of liquid saturation.
*/

#include "MultiphaseDefs.hh"
#include "MultiphaseTypeDefs.hh"
#include "RelPermEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Two constructors.
****************************************************************** */
RelPermEvaluator::RelPermEvaluator(Teuchos::ParameterList& plist,
                                   const Teuchos::RCP<WRMmpPartition>& wrm)
  : SecondaryVariableFieldEvaluator(plist),
    wrm_(wrm)
{
  my_key_ = plist.get<std::string>("my key");
  std::string domain = Keys::getDomain(my_key_);
  saturation_liquid_key_ = plist.get<std::string>("saturation key", Keys::getKey(domain, "saturation_liquid"));

  std::string name = plist.get<std::string>("phase name");
  if (name == "liquid")
    phase_ = MULTIPHASE_PHASE_LIQUID;
  else if (name == "gas")
    phase_ = MULTIPHASE_PHASE_GAS;

  dependencies_.insert(saturation_liquid_key_);
}


RelPermEvaluator::RelPermEvaluator(const RelPermEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    wrm_(other.wrm_) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> RelPermEvaluator::Clone() const {
  return Teuchos::rcp(new RelPermEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void RelPermEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& sat_c = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  auto& result_c = *result->ViewComponent("cell");

  int ncells = result_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = wrm_->second[(*wrm_->first)[c]]->k_relative(sat_c[0][c], phase_);
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void RelPermEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& sat_c = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  auto& result_c = *result->ViewComponent("cell");

  int ncells = result_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = wrm_->second[(*wrm_->first)[c]]->dKdS(sat_c[0][c], phase_);
  }
}


/* ******************************************************************
* Ensure part of my fields exists.
****************************************************************** */
void RelPermEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  AMANZI_ASSERT(my_key_ != std::string(""));
  auto my_fac = S->RequireField(my_key_, my_key_);

  // If my requirements have not yet been set, we'll have to hope they
  // get set by someone later.  For now just defer.
  if (my_fac->Mesh() != Teuchos::null) {
    // Create an unowned factory to check my dependencies.
    auto dep_fac = Teuchos::rcp(new CompositeVectorSpace); 
    dep_fac->SetMesh(my_fac->Mesh())->SetGhosted(true)
        ->AddComponent("cell", AmanziMesh::CELL, 1)->SetOwned(false);

    // Loop over my dependencies, ensuring they meet the requirements.
    for (auto key=dependencies_.begin(); key!=dependencies_.end(); ++key) {
      Teuchos::RCP<CompositeVectorSpace> fac = S->RequireField(*key);
      fac->Update(*dep_fac);
    }

    // Recurse into the tree to propagate info to leaves.
    for (auto key=dependencies_.begin(); key!=dependencies_.end(); ++key) {
      S->RequireFieldEvaluator(*key)->EnsureCompatibility(S);
    }
  }
}

}  // namespace Flow
}  // namespace Amanzi

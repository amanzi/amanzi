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
                                   Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                                   double patm,
                                   const Teuchos::RCP<WRMPartition>& wrm) :
    SecondaryVariableFieldEvaluator(plist),
    mesh_(mesh),
    wrm_(wrm),
    patm_(patm),
    min_value_(0.0),
    max_value_(1.0) {
  InitializeFromPlist_();
}

RelPermEvaluator::RelPermEvaluator(const RelPermEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    wrm_(other.wrm_),
    pressure_key_(other.pressure_key_),
    patm_(other.patm_),
    min_value_(other.min_value_),
    max_value_(other.max_value_) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> RelPermEvaluator::Clone() const {
  return Teuchos::rcp(new RelPermEvaluator(*this));
}


/* ******************************************************************
* Initialization.
****************************************************************** */
void RelPermEvaluator::InitializeFromPlist_()
{
  // my keys is for rel perm.
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("rel perm key", "relative_permeability");
  }

  // my dependency is pressure.
  pressure_key_ = plist_.get<std::string>("pressure key", "pressure");
  dependencies_.insert(pressure_key_);

  // use rel perm class for calcualtion
  Teuchos::ParameterList plist;
  relperm_ = Teuchos::rcp(new RelPerm(plist, mesh_, patm_, wrm_));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void RelPermEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  AMANZI_ASSERT(false);
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
  AMANZI_ASSERT(wrt_key == pressure_key_);
  // relperm_->ComputeDerivative(S->GetFieldData(pressure_key_), result);
}

}  // namespace Flow
}  // namespace Amanzi

/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  The porosity model evaluator simply calls the porosity model with 
  the correct arguments.
*/

#include "FlowDefs.hh"
#include "PorosityModel_Compressible.hh"
#include "PorosityModel_Constant.hh"
#include "PorosityModelEvaluator.hh"
#include "PorosityModelPartition.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Two constructors.
****************************************************************** */
PorosityModelEvaluator::PorosityModelEvaluator(
    Teuchos::ParameterList& plist, Teuchos::RCP<PorosityModelPartition> pom) :
    SecondaryVariablesFieldEvaluator(plist),
    pom_(pom)
{
  InitializeFromPlist_();
}


PorosityModelEvaluator::PorosityModelEvaluator(const PorosityModelEvaluator& other) :
    SecondaryVariablesFieldEvaluator(other),
    pressure_key_(other.pressure_key_),
    pom_(other.pom_) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> PorosityModelEvaluator::Clone() const {
  return Teuchos::rcp(new PorosityModelEvaluator(*this));
}


/* ******************************************************************
* Initialization.
****************************************************************** */
void PorosityModelEvaluator::InitializeFromPlist_() {
  // my key is for saturation
  my_keys_.push_back(plist_.get<std::string>("porosity key"));

  // my dependency is pressure.
  pressure_key_ = plist_.get<std::string>("pressure key");
  dependencies_.insert(pressure_key_);
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void PorosityModelEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const std::vector<Teuchos::Ptr<CompositeVector> >& results)
{
  Epetra_MultiVector& phi_c = *results[0]->ViewComponent("cell", false);
  const Epetra_MultiVector& pres_c = *S->GetFieldData(pressure_key_)->ViewComponent("cell", false);

  int ncells = phi_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    phi_c[0][c] = pom_->second[(*pom_->first)[c]]->Porosity(pres_c[0][c]);
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void PorosityModelEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key,
    const std::vector<Teuchos::Ptr<CompositeVector> >& results)
{
  Epetra_MultiVector& phi_c = *results[0]->ViewComponent("cell", false);
  const Epetra_MultiVector& pres_c = *S->GetFieldData(pressure_key_)->ViewComponent("cell", false);

  int ncells = phi_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    phi_c[0][c] = pom_->second[(*pom_->first)[c]]->dPorositydPressure(pres_c[0][c]);
  }
}

}  // namespace Flow
}  // namespace Amanzi

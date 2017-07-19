/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  The fracture permeability model evaluator simply calls the 
  permeability model with the correct arguments.
*/

#include "FlowDefs.hh"
#include "FracturePermModel_CubicLaw.hh"
#include "FracturePermModel_Linear.hh"
#include "FracturePermModelEvaluator.hh"
#include "FracturePermModelPartition.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Two constructors.
****************************************************************** */
FracturePermModelEvaluator::FracturePermModelEvaluator(
    Teuchos::ParameterList& plist, Teuchos::RCP<FracturePermModelPartition> fpm) :
    SecondaryVariableFieldEvaluator(plist),
    fpm_(fpm)
{
  my_key_ = std::string("fracture_permeability");
  aperture_key_ = std::string("fracture_aperture");
  dependencies_.insert(aperture_key_);
}


FracturePermModelEvaluator::FracturePermModelEvaluator(const FracturePermModelEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    aperture_key_(other.aperture_key_),
    fpm_(other.fpm_) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> FracturePermModelEvaluator::Clone() const {
  return Teuchos::rcp(new FracturePermModelEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void FracturePermModelEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  Epetra_MultiVector& perm_c = *result->ViewComponent("cell", false);
  const Epetra_MultiVector& aperture_c = *S->GetFieldData(aperture_key_)->ViewComponent("cell", false);

  int ncells = perm_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    perm_c[0][c] = fpm_->second[(*fpm_->first)[c]]->Permeability(aperture_c[0][c]);
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void FracturePermModelEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result)
{
  result->PutScalar(0.0);
}

}  // namespace Flow
}  // namespace Amanzi

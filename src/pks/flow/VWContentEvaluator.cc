/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Field evaluator for total volumetric water content which is the 
  conserved quantity in Richards's equation.

  Wrapping this conserved quantity as a field evaluator makes it
  easier to take derivatives, keep updated, and the like.
  The equation for this is simply:

    WC = phi * (s_liquid * n_liquid + X_gas * s_gas * n_gas)

  where X_gas is the molar fraction of water in the gas phase.
*/

#include "VWContentEvaluator.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Constructor.
****************************************************************** */
VWContentEvaluator::VWContentEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {};


/* ******************************************************************
* Initialization.
****************************************************************** */
void VWContentEvaluator::Init_()
{
  my_key_ = std::string("water_content");
  saturation_key_ = plist_.get<std::string>("saturation key");
  porosity_key_ = plist_.get<std::string>("porosity key");

  dependencies_.insert(std::string(porosity_key_));
  dependencies_.insert(std::string(saturation_key_));
  dependencies_.insert(std::string("molar_density_liquid"));

  water_vapor_ = plist_.get<bool>("water vapor", false);
  if (water_vapor_) {
    dependencies_.insert(std::string("molar_density_gas"));
    dependencies_.insert(std::string("molar_fraction_gas"));
  }
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
VWContentEvaluator::VWContentEvaluator(const VWContentEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    water_vapor_(other.water_vapor_) {};


Teuchos::RCP<FieldEvaluator> VWContentEvaluator::Clone() const {
  return Teuchos::rcp(new VWContentEvaluator(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void VWContentEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                                        const Teuchos::Ptr<CompositeVector>& result)
{
  const Epetra_MultiVector& s_l = *S->GetFieldData(saturation_key_)->ViewComponent("cell");
  const Epetra_MultiVector& n_l = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell");

  const Epetra_MultiVector& phi = *S->GetFieldData(porosity_key_)->ViewComponent("cell");
  Epetra_MultiVector& result_v = *result->ViewComponent("cell");

  if (water_vapor_) {
    const Epetra_MultiVector& n_g = *S->GetFieldData("molar_density_gas")->ViewComponent("cell");
    const Epetra_MultiVector& mlf_g = *S->GetFieldData("molar_fraction_gas")->ViewComponent("cell");
    
    int ncells = result->size("cell", false);
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = phi[0][c] * (s_l[0][c] * n_l[0][c]
                                  + (1.0 - s_l[0][c]) * n_g[0][c] * mlf_g[0][c]);
    }
  } else {
    int ncells = result->size("cell", false);
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_l[0][c] * n_l[0][c];
    }
  }      
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void VWContentEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  const Epetra_MultiVector& s_l = *S->GetFieldData(saturation_key_)->ViewComponent("cell");
  const Epetra_MultiVector& n_l = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell");

  const Epetra_MultiVector& phi = *S->GetFieldData(porosity_key_)->ViewComponent("cell");
  Epetra_MultiVector& result_v = *result->ViewComponent("cell");

  int ncells = result->size("cell", false);

  if (water_vapor_) {
    const Epetra_MultiVector& n_g = *S->GetFieldData("molar_density_gas")->ViewComponent("cell");
    const Epetra_MultiVector& mlf_g = *S->GetFieldData("molar_fraction_gas")->ViewComponent("cell");

    if (wrt_key == porosity_key_) {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = (s_l[0][c]*n_l[0][c] + (1.0 - s_l[0][c]) * n_g[0][c] * mlf_g[0][c]);
      }
    } else if (wrt_key == saturation_key_) {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = phi[0][c] * n_l[0][c];
      }
    } else if (wrt_key == "molar_density_liquid") {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = phi[0][c] * s_l[0][c];
      }
    } else if (wrt_key == "molar_density_gas") {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = phi[0][c] * (1.0 - s_l[0][c]) * mlf_g[0][c];
      }
    } else if (wrt_key == "molar_fraction_gas") {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = phi[0][c] * (1.0 - s_l[0][c]) * n_g[0][c];
      }
    } else {
      AMANZI_ASSERT(0);
    }
    
  } else {
    if (wrt_key == porosity_key_) {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = s_l[0][c] * n_l[0][c];
      }
    } else if (wrt_key == "saturation_liquid") {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = phi[0][c] * n_l[0][c];
      }
    } else if (wrt_key == "molar_density_liquid") {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = phi[0][c] * s_l[0][c];
      }
    } else {
      AMANZI_ASSERT(0);
    }
  }
}

}  // namespace Flow
}  // namespace Amanzi

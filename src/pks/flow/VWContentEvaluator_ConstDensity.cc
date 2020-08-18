/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Field evaluator for total volumetric water content which is the 
  conserved quantity in Richards's equation.

  Constant water density.
*/

#include "CommonDefs.hh"
#include "VWContentEvaluator_ConstDensity.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Constructor.
****************************************************************** */
VWContentEvaluator_ConstDensity::VWContentEvaluator_ConstDensity(Teuchos::ParameterList& plist) :
    VWContentEvaluator(plist) {};


/* ******************************************************************
* Initialization.
****************************************************************** */
void VWContentEvaluator_ConstDensity::Init_()
{
  my_key_ = plist_.get<std::string>("water content key");
  pressure_key_ = plist_.get<std::string>("pressure key");
  saturation_key_ = plist_.get<std::string>("saturation key");
  porosity_key_ = plist_.get<std::string>("porosity key");

  dependencies_.insert(porosity_key_);
  dependencies_.insert(saturation_key_);

  water_vapor_ = plist_.get<bool>("water vapor", false);
  if (water_vapor_) {
    dependencies_.insert(std::string("molar_density_gas"));
    dependencies_.insert(std::string("molar_fraction_gas"));
  }
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void VWContentEvaluator_ConstDensity::EvaluateField_(
    const Teuchos::Ptr<State>& S, const Teuchos::Ptr<CompositeVector>& result)
{
  const Epetra_MultiVector& s_l = *S->GetFieldData(saturation_key_)->ViewComponent("cell");
  const Epetra_MultiVector& phi = *S->GetFieldData(porosity_key_)->ViewComponent("cell");

  double rho = *S->GetScalarData("fluid_density");
  double n_l = rho / CommonDefs::MOLAR_MASS_H2O;

  Epetra_MultiVector& result_v = *result->ViewComponent("cell");

  if (water_vapor_) {
    const Epetra_MultiVector& n_g = *S->GetFieldData("molar_density_gas")->ViewComponent("cell");
    const Epetra_MultiVector& mlf_g = *S->GetFieldData("molar_fraction_gas")->ViewComponent("cell");
    
    const Epetra_MultiVector& temp = *S->GetFieldData("temperature")->ViewComponent("cell");
    const Epetra_MultiVector& pres = *S->GetFieldData(pressure_key_)->ViewComponent("cell");
    double patm = *S->GetScalarData("atmospheric_pressure");

    double R = CommonDefs::IDEAL_GAS_CONSTANT_R;

    int ncells = result->size("cell", false);
    for (int c = 0; c != ncells; ++c) {
      double nRT = n_l * temp[0][c] * R;
      double pc = patm - pres[0][c];

      result_v[0][c] = phi[0][c] * (s_l[0][c] * n_l 
                                 + (1.0 - s_l[0][c]) * n_g[0][c] * mlf_g[0][c] * exp(-pc / nRT));
    }
  } else {
    int ncells = result->size("cell", false);
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_l[0][c] * n_l;
    }
  }      
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void VWContentEvaluator_ConstDensity::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  const Epetra_MultiVector& s_l = *S->GetFieldData(saturation_key_)->ViewComponent("cell");
  const Epetra_MultiVector& phi = *S->GetFieldData(porosity_key_)->ViewComponent("cell");

  double rho = *S->GetScalarData("fluid_density");
  double n_l = rho / CommonDefs::MOLAR_MASS_H2O;

  Epetra_MultiVector& result_v = *result->ViewComponent("cell");

  if (water_vapor_) {
    const Epetra_MultiVector& n_g = *S->GetFieldData("molar_density_gas")->ViewComponent("cell");
    const Epetra_MultiVector& mlf_g = *S->GetFieldData("molar_fraction_gas")->ViewComponent("cell");

    int ncells = result->size("cell", false);
    if (wrt_key == porosity_key_) {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = (s_l[0][c] * n_l + (1.0 - s_l[0][c]) * n_g[0][c] * mlf_g[0][c]);
      }
    } else if (wrt_key == saturation_key_) {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = phi[0][c] * n_l;
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
    int ncells = result->size("cell", false);
    if (wrt_key == porosity_key_) {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = s_l[0][c] * n_l;
      }
    } else if (wrt_key == saturation_key_) {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = phi[0][c] * n_l;
      }
    } else {
      AMANZI_ASSERT(0);
    }
  }
}

}  // namespace Flow
}  // namespace Amanzi

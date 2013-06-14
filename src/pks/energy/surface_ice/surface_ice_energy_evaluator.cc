/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Ethan Coon (ecoon@lanl.gov)

FieldEvaluator for water content.

Wrapping this conserved quantity as a field evaluator makes it easier to take
derivatives, keep updated, and the like.  The equation for this is simply:

WC = phi * (s_liquid * n_liquid + omega_gas * s_gas * n_gas)

This is simply the conserved quantity in Richards equation.
----------------------------------------------------------------------------- */


#include "surface_ice_energy_evaluator.hh"

namespace Amanzi {
namespace Energy {

SurfaceIceEnergyEvaluator::SurfaceIceEnergyEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  my_key_ = plist_.get<std::string>("energy key", "surface_energy");
  setLinePrefix(my_key_+std::string(" evaluator"));

  // densities
  dens_key_ = plist.get<std::string>("molar density liquid key",
          "surface_molar_density_liquid");
  dependencies_.insert(dens_key_);

  dens_ice_key_ = plist.get<std::string>("molar density ice key",
          "surface_molar_density_ice");
  dependencies_.insert(dens_ice_key_);

  // internal energies
  ie_key_ = plist.get<std::string>("internal energy liquid key",
          "surface_internal_energy_liquid");
  dependencies_.insert(ie_key_);

  ie_ice_key_ = plist.get<std::string>("internal energy ice key",
          "surface_internal_energy_ice");
  dependencies_.insert(ie_ice_key_);

  // unfrozen fraction
  uf_key_ = plist.get<std::string>("unfrozen fraction key", "unfrozen_fraction");
  dependencies_.insert(uf_key_);

  // ponded depth
  height_key_ = plist.get<std::string>("height key", "ponded_depth");
  dependencies_.insert(height_key_);

  //  dependencies_.insert(std::string("cell_volume"));
};

SurfaceIceEnergyEvaluator::SurfaceIceEnergyEvaluator(const SurfaceIceEnergyEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    dens_key_(other.dens_key_),
    dens_ice_key_(other.dens_ice_key_),
    ie_key_(other.ie_key_),
    ie_ice_key_(other.ie_ice_key_),
    uf_key_(other.uf_key_),
    height_key_(other.height_key_) {};

Teuchos::RCP<FieldEvaluator>
SurfaceIceEnergyEvaluator::Clone() const {
  return Teuchos::rcp(new SurfaceIceEnergyEvaluator(*this));
};


void SurfaceIceEnergyEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  const Epetra_MultiVector& n_l = *S->GetFieldData(dens_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& n_i = *S->GetFieldData(dens_ice_key_)
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& u_l = *S->GetFieldData(ie_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& u_i = *S->GetFieldData(ie_ice_key_)
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& eta = *S->GetFieldData(uf_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& height = *S->GetFieldData(height_key_)
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& cv = *S->GetFieldData("surface_cell_volume")
      ->ViewComponent("cell",false);

  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

  int ncells = result_v.MyLength();
  for (int c=0; c!=ncells; ++c) {
    result_v[0][c] = height[0][c] * ( eta[0][c] * n_l[0][c] * u_l[0][c]
            + (1. - eta[0][c]) * n_i[0][c] * u_i[0][c]);
    result_v[0][c] *= cv[0][c];
  }
};


void SurfaceIceEnergyEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S, Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result) {

  const Epetra_MultiVector& n_l = *S->GetFieldData(dens_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& n_i = *S->GetFieldData(dens_ice_key_)
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& u_l = *S->GetFieldData(ie_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& u_i = *S->GetFieldData(ie_ice_key_)
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& eta = *S->GetFieldData(uf_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& height = *S->GetFieldData(height_key_)
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& cv = *S->GetFieldData("surface_cell_volume")
      ->ViewComponent("cell",false);

  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

  int ncells = result_v.MyLength();
  if (wrt_key == height_key_) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = ( eta[0][c] * n_l[0][c] * u_l[0][c]
            + (1. - eta[0][c]) * n_i[0][c] * u_i[0][c]);
      result_v[0][c] *= cv[0][c];
    }
  } else if (wrt_key == uf_key_) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = height[0][c] * (n_l[0][c] * u_l[0][c] - n_i[0][c] * u_i[0][c]);
      result_v[0][c] *= cv[0][c];
    }
  } else if (wrt_key == dens_key_) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = height[0][c] * eta[0][c] * u_l[0][c];
      result_v[0][c] *= cv[0][c];
    }
  } else if (wrt_key == dens_ice_key_) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = height[0][c] * (1. - eta[0][c]) * u_i[0][c];
      result_v[0][c] *= cv[0][c];
    }
  } else if (wrt_key == ie_key_) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = height[0][c] * eta[0][c] * n_l[0][c];
      result_v[0][c] *= cv[0][c];
    }
  } else if (wrt_key == ie_ice_key_) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = height[0][c] * (1. - eta[0][c]) * n_i[0][c];
      result_v[0][c] *= cv[0][c];
    }
  } else {
    ASSERT(0);
  }
};


} //namespace
} //namespace

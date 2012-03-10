/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon
   ------------------------------------------------------------------------- */
#include "Epetra_MultiVector.h"
#include "passive_tracer.hh"

namespace Amanzi {
namespace Transport {

void PassiveTracer::AddAccumulation_(Teuchos::RCP<CompositeVector> f) {
  Teuchos::RCP<const CompositeVector> poro0 =
    S_inter_->GetFieldData("porosity");
  Teuchos::RCP<const CompositeVector> poro1 =
    S_next_->GetFieldData("porosity");

  Teuchos::RCP<const CompositeVector> conc0 =
    S_inter_->GetFieldData("concentration");
  Teuchos::RCP<const CompositeVector> conc1 =
    S_next_->GetFieldData("concentration");

  Teuchos::RCP<const CompositeVector> sat_liq0 =
    S_inter_->GetFieldData("saturation_liquid");
  Teuchos::RCP<const CompositeVector> sat_liq1 =
    S_next_->GetFieldData("saturation_liquid");

  Teuchos::RCP<const CompositeVector> cell_volume0 =
    S_inter_->GetFieldData("cell_volume");
  Teuchos::RCP<const CompositeVector> cell_volume1 =
    S_next_->GetFieldData("cell_volume");

  double dt = S_next_->time() - S_inter_->time();

  int c_min = S_->mesh()->cell_map(true).MinLID();
  int c_owned = S_->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=c_min; c != c_min+c_owned; ++c) {
    double mols_tracer0 = (*poro0)(c) * (*sat_liq0)(c) * (*conc0)(c);
    double mols_tracer1 = (*poro1)(c) * (*sat_liq1)(c) * (*conc1)(c);

    (*f)(c) += (mols_tracer1 - mols_tracer0)/dt;
  }
}

void PassiveTracer::AddAdvection_(Teuchos::RCP<CompositeVector> f) {
  // set darcy flux
  advection_->set_flux(S_next_->GetFieldData("darcy_flux"));

  // stuff C into the field cells
  Teuchos::RCP<CompositeVector> field = advection_->field();
  Teuchos::RCP<const CompositeVector> conc =
    S_next_->GetFieldData("concentration");
  *field->ViewComponent("cell", false) = *conc->ViewComponent(false);

  // advect and add to residual
  advection_->Apply();

  int c_min = S_->mesh()->cell_map(true).MinLID();
  int c_owned = S_->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=c_min; c != c_min+c_owned; ++c) {
    (*f)(c) += (*field)("cell",c);
  }
};

} //namespace Transport
} //namespace Amanzi

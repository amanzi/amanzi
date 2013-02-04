/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code. 
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
------------------------------------------------------------------------- */


#include "bdf1_time_integrator.hh"
#include "flow_bc_factory.hh"

#include "upwind_cell_centered.hh"
#include "upwind_arithmetic_mean.hh"
#include "upwind_total_flux.hh"
#include "upwind_gravity_flux.hh"

#include "composite_vector_function.hh"
#include "composite_vector_function_factory.hh"

#include "primary_variable_field_evaluator.hh"
#include "wrm_permafrost_evaluator.hh"
#include "rel_perm_evaluator.hh"
#include "permafrost_water_content.hh"

#include "permafrost.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory<Permafrost> Permafrost::reg_("permafrost flow");


// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void Permafrost::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  // -- Absolute permeability.
  //       For now, we assume scalar permeability.  This will change.
  S->RequireField("permeability")->SetMesh(S->GetMesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("permeability");

  // -- water content, and evaluator
  S->RequireField("water_content")->SetMesh(S->GetMesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList wc_plist = plist_.sublist("water content evaluator");
  Teuchos::RCP<PermafrostWaterContent> wc =
      Teuchos::rcp(new PermafrostWaterContent(wc_plist));
  S->SetFieldEvaluator("water_content", wc);

  // -- Water retention evaluators, for saturation and rel perm.
  Teuchos::ParameterList wrm_plist = plist_.sublist("water retention evaluator");
  Teuchos::RCP<FlowRelations::WRMPermafrostEvaluator> wrm =
      Teuchos::rcp(new FlowRelations::WRMPermafrostEvaluator(wrm_plist));
  S->SetFieldEvaluator("saturation_liquid", wrm);
  S->SetFieldEvaluator("saturation_gas", wrm);
  S->SetFieldEvaluator("saturation_ice", wrm);

  // -- the rel perm evaluator, also with the same underlying WRM.
  S->RequireField("relative_permeability")->SetMesh(S->GetMesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::RCP<FlowRelations::RelPermEvaluator> rel_perm_evaluator =
      Teuchos::rcp(new FlowRelations::RelPermEvaluator(wrm_plist, wrm->get_WRMs()));
  S->SetFieldEvaluator("relative_permeability", rel_perm_evaluator);

  // -- Liquid density and viscosity for the transmissivity.
  S->RequireField("molar_density_liquid")->SetMesh(S->GetMesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("molar_density_liquid");

  S->RequireField("viscosity_liquid")->SetMesh(S->GetMesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("viscosity_liquid");

  // -- liquid mass density for the gravity fluxes
  S->RequireField("mass_density_liquid")->SetMesh(S->GetMesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("mass_density_liquid"); // simply picks up the molar density one.
}


// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
bool Permafrost::UpdatePermeabilityData_(const Teuchos::Ptr<State>& S) {
  Teuchos::RCP<CompositeVector> uw_rel_perm = S->GetFieldData("numerical_rel_perm", name_);
  Teuchos::RCP<const CompositeVector> rel_perm = S->GetFieldData("relative_permeability");
  bool update_perm = S->GetFieldEvaluator("relative_permeability")
      ->HasFieldChanged(S, name_);

  update_perm |= S->GetFieldEvaluator("molar_density_liquid")->HasFieldChanged(S, name_);
  update_perm |= S->GetFieldEvaluator("viscosity_liquid")->HasFieldChanged(S, name_);

  // requirements due to the upwinding method
  if (Krel_method_ == FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX) {
    bool update_dir = S->GetFieldEvaluator("mass_density_liquid")
        ->HasFieldChanged(S, name_);
    update_dir |= S->GetFieldEvaluator(key_)->HasFieldChanged(S, name_);

    if (update_dir) {
      // update the direction of the flux -- note this is NOT the flux
      Teuchos::RCP<CompositeVector> flux_dir =
          S->GetFieldData("darcy_flux_direction", name_);

      // Create the stiffness matrix without a rel perm (rel perm = 1)
      matrix_->CreateMFDstiffnessMatrices(Teuchos::null);

      // Derive the pressure fluxes
      Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(key_);
      matrix_->DeriveFlux(*pres, flux_dir.ptr());

      // Add in the gravity fluxes
      Teuchos::RCP<const Epetra_Vector> gvec = S->GetConstantVectorData("gravity");
      Teuchos::RCP<const CompositeVector> rho = S->GetFieldData("mass_density_liquid");
      AddGravityFluxesToVector_(gvec.ptr(), Teuchos::null, rho.ptr(), flux_dir.ptr());
      flux_dir->ScatterMasterToGhosted();
    }

    update_perm |= update_dir;
  }

  if (update_perm) {
    // patch up the BCs
    const double& p_atm = *S->GetScalarData("atmospheric_pressure");
    const Epetra_MultiVector& pres_f = *S->GetFieldData(key_)
        ->ViewComponent("face",false);
    for (int f=0; f!=uw_rel_perm->size("face",false); ++f) {
      AmanziMesh::Entity_ID_List cells;
      uw_rel_perm->mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
      if (cells.size() == 1) {
        (*uw_rel_perm)("face",f) = wrms_->second[ (*wrms_->first)[cells[0]] ]
          ->saturation(p_atm - pres_f[0][f]);
      }
    }

    if (coupled_to_surface_via_residual_ ||
        coupled_to_surface_via_residual_new_ ||
        coupled_to_surface_via_full_ ||
        coupled_to_surface_via_head_ ||
        coupled_to_surface_via_flux_) {
      // patch up the rel perm on surface
      Teuchos::RCP<const AmanziMesh::Mesh> surface = S->GetMesh("surface");

      int ncells_surface = surface->num_entities(AmanziMesh::CELL,AmanziMesh::OWNED);

      for (int c=0; c!=ncells_surface; ++c) {
        // -- get the surface cell's equivalent subsurface face and neighboring cell
        AmanziMesh::Entity_ID f = surface->entity_get_parent(AmanziMesh::CELL, c);
        (*uw_rel_perm)("face",f) = 1.0;
      }
    }

    // upwind
    upwinding_->Update(S);
  }

  // Scale cells by n/visc if needed.
  if (update_perm) {
    const Epetra_MultiVector& n_liq = *S->GetFieldData("molar_density_liquid")
        ->ViewComponent("cell",false);
    const Epetra_MultiVector& visc = *S->GetFieldData("viscosity_liquid")
        ->ViewComponent("cell",false);

    for (int c=0; c!=uw_rel_perm->size("cell", false); ++c) {
      (*uw_rel_perm)("cell",c) *= n_liq[0][c] / visc[0][c];
    }

    // communicate
    //    uw_rel_perm->ScatterMasterToGhosted();
  }

  return update_perm;
};


} // namespace
} // namespace

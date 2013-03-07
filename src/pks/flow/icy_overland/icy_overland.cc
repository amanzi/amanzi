/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code. 
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
------------------------------------------------------------------------- */
#include "independent_variable_field_evaluator.hh"

#include "upwinding.hh"
#include "upwind_potential_difference.hh"
#include "elevation_evaluator.hh"
#include "meshed_elevation_evaluator.hh"
#include "standalone_elevation_evaluator.hh"
#include "overland_conductivity_evaluator.hh"
#include "overland_conductivity_model.hh"
#include "unfrozen_effective_depth_evaluator.hh"
#include "unfrozen_fraction_evaluator.hh"
#include "unfrozen_fraction_model.hh"
#include "icy_height_model.hh"
#include "icy_height_evaluator.hh"


#include "overland_head_icy_water_content_evaluator.hh"
#include "icy_overland.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory<IcyOverlandFlow> IcyOverlandFlow::reg_("overland flow with ice");

void IcyOverlandFlow::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  // -- evaluator for surface geometry.
  S->RequireField("elevation")->SetMesh(mesh_)->SetGhosted()
                ->SetComponents(names2, locations2, num_dofs2);
  S->RequireField("slope_magnitude")->SetMesh(mesh_)
                ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField("pres_elev")->SetMesh(mesh_)->SetGhosted()
                ->SetComponents(names2, locations2, num_dofs2);

  Teuchos::RCP<FlowRelations::ElevationEvaluator> elev_evaluator;
  if (standalone_mode_) {
    ASSERT(plist_.isSublist("elevation evaluator"));
    Teuchos::ParameterList elev_plist = plist_.sublist("elevation evaluator");
    elev_evaluator = Teuchos::rcp(new FlowRelations::StandaloneElevationEvaluator(elev_plist));
  } else {
    Teuchos::ParameterList elev_plist = plist_.sublist("elevation evaluator");
    elev_evaluator = Teuchos::rcp(new FlowRelations::MeshedElevationEvaluator(elev_plist));
  }
  S->SetFieldEvaluator("elevation", elev_evaluator);
  S->SetFieldEvaluator("slope_magnitude", elev_evaluator);
  S->SetFieldEvaluator("pres_elev", elev_evaluator);

  // -- unfrozen fraction model, eta
  S->RequireField("unfrozen_fraction")->SetMesh(mesh_)
                ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
  ASSERT(plist_.isSublist("unfrozen fraction evaluator"));
  Teuchos::ParameterList uf_plist = plist_.sublist("unfrozen fraction evaluator");
  Teuchos::RCP<FlowRelations::UnfrozenFractionEvaluator> uf_evaluator =
      Teuchos::rcp(new FlowRelations::UnfrozenFractionEvaluator(uf_plist));
  S->SetFieldEvaluator("unfrozen_fraction", uf_evaluator);
  uf_model_ = uf_evaluator->get_Model();

  // -- unfrozen effective depth, h * eta
  S->RequireField("unfrozen_effective_depth")->SetMesh(mesh_)
                ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList h_eff_plist = plist_.sublist("unfrozen effective depth");
  Teuchos::RCP<FlowRelations::UnfrozenEffectiveDepthEvaluator> h_eff_evaluator =
      Teuchos::rcp(new FlowRelations::UnfrozenEffectiveDepthEvaluator(h_eff_plist));
  S->SetFieldEvaluator("unfrozen_effective_depth", h_eff_evaluator);

  // -- conductivity evaluator
  S->RequireField("overland_conductivity")->SetMesh(mesh_)
                ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
  ASSERT(plist_.isSublist("overland conductivity evaluator"));
  Teuchos::ParameterList cond_plist = plist_.sublist("overland conductivity evaluator");
  // -- set the height key to be eta * h, not just h, for the frozen case.
  //    NOTE that this will automatically require an unfrozen_fraction
  //    field/model
  cond_plist.set("height key", "unfrozen_effective_depth");
  Teuchos::RCP<FlowRelations::OverlandConductivityEvaluator> cond_evaluator =
      Teuchos::rcp(new FlowRelations::OverlandConductivityEvaluator(cond_plist));
  S->SetFieldEvaluator("overland_conductivity", cond_evaluator);
  cond_model_ = cond_evaluator->get_Model();

  // -- source term evaluator
  if (plist_.isSublist("source evaluator")) {
    Teuchos::ParameterList source_plist = plist_.sublist("source evaluator");
    source_plist.set("evaluator name", "overland_source");
    is_source_term_ = true;
    S->RequireField("overland_source")->SetMesh(mesh_)
        ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
    Teuchos::RCP<FieldEvaluator> source_evaluator =
        Teuchos::rcp(new IndependentVariableFieldEvaluator(source_plist));
    S->SetFieldEvaluator("overland_source", source_evaluator);
  }

  // -- water content
  S->RequireField("surface_water_content")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList wc_plist = plist_.sublist("overland water content");
  Teuchos::RCP<FlowRelations::OverlandHeadIcyWaterContentEvaluator> wc_evaluator =
      Teuchos::rcp(new FlowRelations::OverlandHeadIcyWaterContentEvaluator(wc_plist));
  S->SetFieldEvaluator("surface_water_content", wc_evaluator);

  // -- ponded depth
  S->RequireField("ponded_depth")->SetMesh(mesh_)->SetGhosted()
                ->AddComponents(names2, locations2, num_dofs2);
  Teuchos::ParameterList height_plist = plist_.sublist("ponded depth");
  Teuchos::RCP<FlowRelations::IcyHeightEvaluator> height_evaluator =
      Teuchos::rcp(new FlowRelations::IcyHeightEvaluator(height_plist));
  S->SetFieldEvaluator("ponded_depth", height_evaluator);
  icy_height_model_ = height_evaluator->get_IcyModel();

  // -- overwrite the upwinding, which was created with ponded_depth, to use
  // -- unfrozen ponded depth.
  upwinding_ = Teuchos::rcp(new Operators::UpwindPotentialDifference(name_,
          "overland_conductivity", "upwind_overland_conductivity",
          "pres_elev", "unfrozen_effective_depth"));

}



// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
//
//  This version is identical to the overland.cc version except that
//    it uses eta * h instead of just h for the upwinding hack.  This
//    should be cleaned up eventually so that it doesn't need this
//    hack at all... -- etc
// -----------------------------------------------------------------------------
bool IcyOverlandFlow::UpdatePermeabilityData_(const Teuchos::Ptr<State>& S) {
  bool update_perm = S->GetFieldEvaluator("overland_conductivity")
      ->HasFieldChanged(S, name_);

  update_perm |= S->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S, name_);
  update_perm |= S->GetFieldEvaluator("pres_elev")->HasFieldChanged(S, name_);
  update_perm |= perm_update_required_;

  if (update_perm) {
    // Update the perm only if needed.
    perm_update_required_ = false;

    // This needs fixed to use the model, not assume a model.
    // Then it needs to be fixed to use a smart evaluator which picks
    // vals from cells and faces.
    // Then it needs to be fixed to work on a CompositeVector whose
    // components are cells and boundary faces. --etc
    const Epetra_MultiVector& depth = *S->GetFieldData("ponded_depth")
      ->ViewComponent("face", false);
    const Epetra_MultiVector& slope = *S->GetFieldData("slope_magnitude")
      ->ViewComponent("cell", false);
    const Epetra_MultiVector& coef = *S->GetFieldData("manning_coefficient")
      ->ViewComponent("cell", false);
    const Epetra_MultiVector& temp = *S->GetFieldData("surface_temperature")
      ->ViewComponent("face", false);

    Teuchos::RCP<CompositeVector> upwind_conductivity =
        S->GetFieldData("upwind_overland_conductivity", name_);
    Epetra_MultiVector& upwind_conductivity_f =
        *upwind_conductivity->ViewComponent("face",false);
    Epetra_MultiVector& upwind_conductivity_c =
        *upwind_conductivity->ViewComponent("cell",false);

    // initialize the face coefficients
    upwind_conductivity_f.PutScalar(0.0);
    upwind_conductivity_c.PutScalar(1.0);

    // First do the boundary
    AmanziMesh::Entity_ID_List cells;
    int nfaces = upwind_conductivity_f.MyLength();
    for (int f=0; f!=nfaces; ++f) {
      if (bc_markers_[f] != Operators::MATRIX_BC_NULL) {
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
        ASSERT(cells.size() == 1);
        int c = cells[0];
        double eff_depth = depth[0][f] * uf_model_->UnfrozenFraction(temp[0][f]);
        (*upwind_conductivity)("face",f) =
          cond_model_->Conductivity(eff_depth, slope[0][c], coef[0][c]);
      }
    }

    // Patch up zero-gradient case, which should not upwind.
    const Epetra_MultiVector& conductivity_c =
        *S->GetFieldData("overland_conductivity")->ViewComponent("cell",false);
    for (Functions::BoundaryFunction::Iterator bc=bc_zero_gradient_->begin();
         bc!=bc_zero_gradient_->end(); ++bc) {
      int f = bc->first;
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      ASSERT(cells.size() == 1);
      int c = cells[0];
      upwind_conductivity_f[0][f] = conductivity_c[0][c];
    }

    // Then upwind.  This overwrites the boundary if upwinding says so.
    upwinding_->Update(S);

    // Scale cells by n
    const Epetra_MultiVector& n_liq = *S->GetFieldData("surface_molar_density_liquid")
        ->ViewComponent("cell",false);
    int ncells = upwind_conductivity_c.MyLength();
    for (int c=0; c!=ncells; ++c) {
      upwind_conductivity_c[0][c] *= n_liq[0][c];
    }

    //    std::cout << "cond in update perm = " << upwind_conductivity_f[0][11] << std::endl;

    // Communicate.  This could be done later, but i'm not exactly sure where, so
    // we'll do it here.
    //    upwind_conductivity->ScatterMasterToGhosted();
  }

  return update_perm;
}

// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void IcyOverlandFlow::UpdateBoundaryConditions_(const Teuchos::Ptr<State>& S) {
  AmanziMesh::Entity_ID_List cells;

  const Epetra_MultiVector& elevation = *S->GetFieldData("elevation")
      ->ViewComponent("face",false);

  // initialize all as null
  for (int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::MATRIX_BC_NULL;
    bc_values_[n] = 0.0;
  }

  // Head BCs are standard Dirichlet, plus elevation
  for (Functions::BoundaryFunction::Iterator bc=bc_head_->begin();
       bc!=bc_head_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
    bc_values_[f] = bc->second + elevation[0][f];
  }

  // pressure BCs require calculating head(pressure)
  const double& p_atm = *S->GetScalarData("atmospheric_pressure");
  const Epetra_Vector& gravity = *S->GetConstantVectorData("gravity");
  double gz = -gravity[2];
  S->GetFieldEvaluator("surface_mass_density_liquid")->HasFieldChanged(S.ptr(), name_);
  const Epetra_MultiVector& rho_l = *S->GetFieldData("surface_mass_density_liquid")
      ->ViewComponent("cell",false);
  S->GetFieldEvaluator("surface_mass_density_ice")->HasFieldChanged(S.ptr(), name_);
  const Epetra_MultiVector& rho_i = *S->GetFieldData("surface_mass_density_ice")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& temp = *S->GetFieldData("surface_temperature")
      ->ViewComponent("face", false);

  for (Functions::BoundaryFunction::Iterator bc=bc_pressure_->begin();
       bc!=bc_pressure_->end(); ++bc) {
    int f = bc->first;
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    ASSERT( cells.size()==1 ) ;

    double eta = uf_model_->UnfrozenFraction(temp[0][f]);
    double height = icy_height_model_->Height(bc->second, eta, rho_l[0][cells[0]],
            rho_i[0][cells[0]], p_atm, gz);
    *out_ << "BC pressure (f" << f << "): pres = " << bc->second << ", eta = "
          << eta << " height = " << height << std::endl;

    bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
    bc_values_[f] = height + elevation[0][f];
  }

  // Standard Neumann data for flux
  for (Functions::BoundaryFunction::Iterator bc=bc_flux_->begin();
       bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MATRIX_BC_FLUX;
    bc_values_[f] = bc->second;
  }

  // zero gradient: grad h = 0 implies that q = -k grad z
  for (Functions::BoundaryFunction::Iterator bc=bc_zero_gradient_->begin();
       bc!=bc_zero_gradient_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MATRIX_BC_FLUX;
    bc_values_[f] = 0.;
  }

}



} // namespace
} // namespace

/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code. 
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
------------------------------------------------------------------------- */
#include "independent_variable_field_evaluator.hh"

#include "upwinding.hh"
#include "upwind_potential_difference.hh"
#include "pres_elev_evaluator.hh"
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
  std::vector<AmanziMesh::Entity_kind> locations2(2), locations_bf(2);
  std::vector<std::string> names2(2), names_bf(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";
  locations_bf[0] = AmanziMesh::CELL;
  locations_bf[1] = AmanziMesh::BOUNDARY_FACE;
  names_bf[0] = "cell";
  names_bf[1] = "boundary_face";

  // -- evaluator for surface geometry.
  S->RequireField("elevation")->SetMesh(mesh_)->SetGhosted()
      ->SetComponents(names2, locations2, num_dofs2);
  S->RequireField("slope_magnitude")->SetMesh(S->GetMesh("surface"))
      ->AddComponents(names_bf, locations_bf, num_dofs2);

  Teuchos::RCP<FlowRelations::ElevationEvaluator> elev_evaluator;
  if (standalone_mode_) {
    ASSERT(plist_.isSublist("elevation evaluator"));
    Teuchos::ParameterList elev_plist = plist_.sublist("elevation evaluator");
    elev_plist.set("manage communication", true);
    elev_evaluator = Teuchos::rcp(new FlowRelations::StandaloneElevationEvaluator(elev_plist));
  } else {
    Teuchos::ParameterList elev_plist = plist_.sublist("elevation evaluator");
    elev_plist.set("manage communication", true);
    elev_evaluator = Teuchos::rcp(new FlowRelations::MeshedElevationEvaluator(elev_plist));
  }
  S->SetFieldEvaluator("elevation", elev_evaluator);
  S->SetFieldEvaluator("slope_magnitude", elev_evaluator);

  // -- evaluator for potential field, h + z
  S->RequireField("pres_elev")->SetMesh(S->GetMesh("surface"))->SetGhosted()
      ->AddComponents(names2, locations2, num_dofs2);

  Teuchos::ParameterList pres_elev_plist = plist_.sublist("potential evaluator");
  pres_elev_plist.set("manage communication", true);
  Teuchos::RCP<FlowRelations::PresElevEvaluator> pres_elev_eval =
      Teuchos::rcp(new FlowRelations::PresElevEvaluator(pres_elev_plist));
  S->SetFieldEvaluator("pres_elev", pres_elev_eval);

  // -- unfrozen fraction model, eta
  S->RequireField("unfrozen_fraction")->SetMesh(mesh_)->SetGhosted()
      ->AddComponents(names_bf, locations_bf, num_dofs2);

  ASSERT(plist_.isSublist("unfrozen fraction evaluator"));
  Teuchos::ParameterList uf_plist = plist_.sublist("unfrozen fraction evaluator");
  Teuchos::RCP<FlowRelations::UnfrozenFractionEvaluator> uf_evaluator =
      Teuchos::rcp(new FlowRelations::UnfrozenFractionEvaluator(uf_plist));
  S->SetFieldEvaluator("unfrozen_fraction", uf_evaluator);
  uf_model_ = uf_evaluator->get_Model();

  // -- unfrozen effective depth, h * eta
  S->RequireField("unfrozen_effective_depth")->SetMesh(mesh_)->SetGhosted()
      ->AddComponents(names_bf, locations_bf, num_dofs2);

  Teuchos::ParameterList h_eff_plist = plist_.sublist("unfrozen effective depth");
  Teuchos::RCP<FlowRelations::UnfrozenEffectiveDepthEvaluator> h_eff_evaluator =
      Teuchos::rcp(new FlowRelations::UnfrozenEffectiveDepthEvaluator(h_eff_plist));
  S->SetFieldEvaluator("unfrozen_effective_depth", h_eff_evaluator);

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
  height_model_ = height_evaluator->get_Model();

  // -- conductivity evaluator
  locations2[1] = AmanziMesh::BOUNDARY_FACE;
  names2[1] = "boundary_face";

  S->RequireField("overland_conductivity")->SetMesh(mesh_)->SetGhosted()
      ->AddComponents(names2, locations2, num_dofs2);
  ASSERT(plist_.isSublist("overland conductivity evaluator"));
  Teuchos::ParameterList cond_plist = plist_.sublist("overland conductivity evaluator");
  // -- set the height key to be eta * h, not just h, for the frozen case.
  cond_plist.set("height key", "unfrozen_effective_depth");
  Teuchos::RCP<FlowRelations::OverlandConductivityEvaluator> cond_evaluator =
      Teuchos::rcp(new FlowRelations::OverlandConductivityEvaluator(cond_plist));
  S->SetFieldEvaluator("overland_conductivity", cond_evaluator);

  // -- overwrite the upwinding, which was created with ponded_depth, to use
  // -- unfrozen ponded depth.
  upwinding_ = Teuchos::rcp(new Operators::UpwindPotentialDifference(name_,
          "overland_conductivity", "upwind_overland_conductivity",
          "pres_elev", "unfrozen_effective_depth"));

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

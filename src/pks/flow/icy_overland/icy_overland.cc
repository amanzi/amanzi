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
#include "unfrozen_fraction_evaluator.hh"
#include "unfrozen_fraction_model.hh"

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
  S->RequireField("elevation")->SetMesh(S->GetMesh("surface"))->SetGhosted()
                ->SetComponents(names2, locations2, num_dofs2);
  S->RequireField("slope_magnitude")->SetMesh(S->GetMesh("surface"))
                ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField("pres_elev")->SetMesh(S->GetMesh("surface"))->SetGhosted()
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
  S->RequireField("unfrozen_fraction")->SetMesh(S->GetMesh("surface"))
                ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
  ASSERT(plist_.isSublist("unfrozen fraction evaluator"));
  Teuchos::ParameterList uf_plist = plist_.sublist("unfrozen fraction evaluator");
  Teuchos::RCP<FlowRelations::UnfrozenFractionEvaluator> uf_evaluator =
      Teuchos::rcp(new FlowRelations::UnfrozenFractionEvaluator(uf_plist));
  S->SetFieldEvaluator("unfrozen_fraction", uf_evaluator);
  uf_model_ = uf_evaluator->get_Model();

  // -- "rel perm" evaluator
  S->RequireField("overland_conductivity")->SetMesh(S->GetMesh("surface"))
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

  // -- coupling term evaluator
  if (plist_.isSublist("subsurface coupling evaluator")) {
    is_coupling_term_ = true;
    S->RequireField("overland_source_from_subsurface")
        ->SetMesh(mesh_)->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList source_plist = plist_.sublist("subsurface coupling evaluator");
    source_plist.set("surface mesh key", "surface");
    source_plist.set("subsurface mesh key", "domain");
    // NOTE: this should change to "overland_molar_density_liquid" or
    // whatever when such a thing exists.
    source_plist.set("source key", "overland_source_from_subsurface");

    Teuchos::RCP<FieldEvaluator> source_evaluator =
        S->RequireFieldEvaluator("overland_source_from_subsurface", source_plist);
  }


  // -- cell volume and evaluator
  S->RequireField("surface_cell_volume")->SetMesh(mesh_)->SetGhosted()
                                ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("surface_cell_volume");
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

  update_perm |= S->GetFieldEvaluator("pres_elev")->HasFieldChanged(S, name_);

  if (update_perm) {
    // Update the perm only if needed.

    // This needs fixed to use the model, not assume a model.
    // Then it needs to be fixed to use a smart evaluator which picks
    // vals from cells and faces.
    // Then it needs to be fixed to work on a CompositeVector whose
    // components are cells and boundary faces. --etc
    const Epetra_MultiVector& depth = *S->GetFieldData(key_)
      ->ViewComponent("face", false);
    const Epetra_MultiVector& slope = *S->GetFieldData("slope_magnitude")
      ->ViewComponent("cell", false);
    const Epetra_MultiVector& coef = *S->GetFieldData("manning_coefficient")
      ->ViewComponent("cell", false);
    const Epetra_MultiVector& temp = *S->GetFieldData("surface_temperature")
      ->ViewComponent("face", false);
    Teuchos::RCP<CompositeVector> upwind_conductivity =
        S->GetFieldData("upwind_overland_conductivity", name_);

    // initialize the face coefficients
    upwind_conductivity->ViewComponent("face",true)->PutScalar(0.0);
    if (upwind_conductivity->has_component("cell")) {
      upwind_conductivity->ViewComponent("cell",true)->PutScalar(1.0);
    }

    // First do the boundary
    AmanziMesh::Entity_ID_List cells;
    int nfaces = upwind_conductivity->size("face", true);
    for (int f=0; f!=nfaces; ++f) {
      if (bc_markers_[f] != Operators::MFD_BC_NULL) {
        upwind_conductivity->mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
        ASSERT(cells.size() == 1);
        int c = cells[0];
        double eff_depth = depth[0][f] * uf_model_->UnfrozenFraction(temp[0][f]);
        (*upwind_conductivity)("face",f) =
          cond_model_->Conductivity(eff_depth, slope[0][c], coef[0][c]);
      }
    }

    // Then upwind.  This overwrites the boundary if upwinding says so.
    upwinding_->Update(S);

    // Communicate.  This could be done later, but i'm not exactly sure where, so
    // we'll do it here.
    upwind_conductivity->ScatterMasterToGhosted("face");
  }

  return update_perm;
}



} // namespace
} // namespace

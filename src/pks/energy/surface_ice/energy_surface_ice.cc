/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Process kernel for energy equation for overland flow.
------------------------------------------------------------------------- */

#include "eos_evaluator.hh"
#include "iem_evaluator.hh"
#include "thermal_conductivity_surface_evaluator.hh"
#include "surface_ice_energy_evaluator.hh"
#include "enthalpy_evaluator.hh"
#include "energy_bc_factory.hh"
#include "matrix_mfd_tpfa.hh"

#include "energy_surface_ice.hh"

namespace Amanzi {
namespace Energy {

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
EnergySurfaceIce::EnergySurfaceIce(Teuchos::ParameterList& plist,
        const Teuchos::RCP<TreeVector>& solution) :
    PKDefaultBase(plist, solution),
    EnergyBase(plist, solution) {
  
  plist_.set("primary variable key", "surface_pressure");
  plist_.set("domain name", "surface");
}


// -------------------------------------------------------------
// Create the physical evaluators for energy, enthalpy, thermal
// conductivity, and any sources.
// -------------------------------------------------------------
void EnergySurfaceIce::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  standalone_mode_ = S->GetMesh() == S->GetMesh("surface");

  // Get data and evaluators needed by the PK
  // -- energy, the conserved quantity
  S->RequireField(energy_key_)->SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList ee_plist = plist_.sublist("energy evaluator");
  ee_plist.set("energy key", energy_key_);
  Teuchos::RCP<SurfaceIceEnergyEvaluator> ee =
    Teuchos::rcp(new SurfaceIceEnergyEvaluator(ee_plist));
  S->SetFieldEvaluator(energy_key_, ee);

  // -- advection of enthalpy
  S->RequireField(enthalpy_key_)->SetMesh(mesh_)
    ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList enth_plist = plist_.sublist("enthalpy evaluator");
  enth_plist.set("enthalpy key", enthalpy_key_);
  Teuchos::RCP<EnthalpyEvaluator> enth =
    Teuchos::rcp(new EnthalpyEvaluator(enth_plist));
  S->SetFieldEvaluator(enthalpy_key_, enth);

  // -- thermal conductivity
  S->RequireField(conductivity_key_)->SetMesh(mesh_)
    ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList tcm_plist =
    plist_.sublist("thermal conductivity evaluator");
  Teuchos::RCP<EnergyRelations::ThermalConductivitySurfaceEvaluator> tcm =
    Teuchos::rcp(new EnergyRelations::ThermalConductivitySurfaceEvaluator(tcm_plist));
  S->SetFieldEvaluator(conductivity_key_, tcm);

  // -- source terms
  // -- -- external source of energy
  if (plist_.isSublist("energy source evaluator")) {
    Teuchos::ParameterList source_plist = plist_.sublist("energy source evaluator");
    source_plist.set("evaluator name", "surface_energy_source");
    is_energy_source_term_ = true;
    S->RequireField("surface_energy_source")->SetMesh(mesh_)
        ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
    Teuchos::RCP<FieldEvaluator> source_evaluator =
        Teuchos::rcp(new IndependentVariableFieldEvaluator(source_plist));
    S->SetFieldEvaluator("surface_energy_source", source_evaluator);
  }

  // -- -- Air temperature
  if (plist_.isSublist("air temperature")) {
    Teuchos::ParameterList air_plist = plist_.sublist("air temperature");

    FunctionFactory fac;
    air_temp_ = Teuchos::rcp(fac.Create(air_plist.sublist("temperature function")));
    
    // If there is air temperature and a mass source, this provides enthalpy
    is_mass_source_term_ = S->IsFieldEvaluator("overland_source");

    // If there is air temperature and a conductivity, this provides heat transfer.
    if (air_plist.isParameter("air to surface conductivity")) {
      is_air_conductivity_ = true;
      K_surface_to_air_ = air_plist.get<double>("air to surface conductivity");
    }

  } else {
    // Make sure there is no mass source.
    if (S->IsFieldEvaluator("overland_source")) {
      Errors::Message m("Overland Energy equation found a mass source without an air temperature.");
      Exceptions::amanzi_throw(m);
    }
  }

  // coupling to subsurface -- kill the preconditioner and replace with a TPFA precon
  coupled_to_subsurface_via_full_ =
      plist_.get<bool>("coupled to subsurface via full coupler", false);
  if (coupled_to_subsurface_via_full_) {
    Teuchos::ParameterList mfd_pc_plist = plist_.sublist("Diffusion PC");
    Teuchos::RCP<Operators::Matrix> precon =
        Teuchos::rcp(new Operators::MatrixMFD_TPFA(mfd_pc_plist, mesh_));
    set_preconditioner(precon);
  }

  // Many quantities are based upon face areas, which are not the cell volume,
  // as the surface mesh has been flattened.
  if (!standalone_mode_) {
    S->RequireFieldEvaluator("surface_3D_cell_volume");
  }

}


// -------------------------------------------------------------
// Initialize the needed models to plug in enthalpy.
// -------------------------------------------------------------
void EnergySurfaceIce::initialize(const Teuchos::Ptr<State>& S) {
  // Call the base class's initialize.
  EnergyBase::initialize(S);

  // For the boundary conditions, we currently hack in the enthalpy to
  // the boundary faces to correctly advect in a Dirichlet temperature
  // BC.  This requires density and internal energy, which in turn
  // require a model based on p,T.
  // This will be removed once boundary faces are implemented.
  Teuchos::RCP<FieldEvaluator> eos_fe =
      S->GetFieldEvaluator("surface_molar_density_liquid");
  Teuchos::RCP<Relations::EOSEvaluator> eos_eval =
    Teuchos::rcp_dynamic_cast<Relations::EOSEvaluator>(eos_fe);
  ASSERT(eos_eval != Teuchos::null);
  eos_liquid_ = eos_eval->get_EOS();

  Teuchos::RCP<FieldEvaluator> iem_fe =
      S->GetFieldEvaluator("surface_internal_energy_liquid");
  Teuchos::RCP<EnergyRelations::IEMEvaluator> iem_eval =
    Teuchos::rcp_dynamic_cast<EnergyRelations::IEMEvaluator>(iem_fe);
  ASSERT(iem_eval != Teuchos::null);
  iem_liquid_ = iem_eval->get_IEM();
}


// -------------------------------------------------------------
// Plug enthalpy into the boundary faces manually.
// This will be removed once boundary faces exist.
// -------------------------------------------------------------
void EnergySurfaceIce::ApplyDirichletBCsToEnthalpy_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& enth) {

  // Since we don't have a surface pressure on faces, this gets a bit uglier.
  // Simply grab the internal cell for now.
  const Epetra_MultiVector& flux = *S->GetFieldData(flux_key_)
      ->ViewComponent("face",false);

  const Epetra_MultiVector& enth_c = *enth->ViewComponent("cell",false);
  Epetra_MultiVector& enth_f = *enth->ViewComponent("face",false);

  Amanzi::Entity_ID_List cells;
  for (Functions::BoundaryFunction::Iterator bc = bc_temperature_->begin();
       bc!=bc_temperature_->end(); ++bc) {
    Amanzi::Entity_ID f = bc->first;
    mesh_->face_get_cells(f, AmanziMesh::USED, cells);
    ASSERT(cells.size() == 1);

    enth_f[0][f] = enth_c[0][cells[0]] * fabs(flux[0][f]);
  }
}


// -------------------------------------------------------------
// Deal with the many source terms.
// -------------------------------------------------------------
virtual void EnergySurfaceIce::AddSources_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g) {
  Epetra_MultiVector& g_c = *g->ViewComponent("cell",false);

  // Note that energy sources are based upon face areas, which are not cell volumes.
  Teuchos::RCP<const Epetra_MultiVector> fa0;
  Teuchos::RCP<const Epetra_MultiVector> fa1;

  if (standalone_mode_) {
    fa0 = S_->GetFieldData(cell_vol_key_)->ViewComponent("cell",false);
    fa1 = S_next_->GetFieldData(cell_vol_key_)->ViewComponent("cell",false);
  } else {
    fa0 = S_->GetFieldData("surface_3D_cell_volume")
        ->ViewComponent("cell",false);
    fa1 = S_next_->GetFieldData("surface_3D_cell_volume")
        ->ViewComponent("cell",false);
  }      


  // external sources of energy
  if (is_energy_source_term_) {
    // Add in external source term.
    S_next_->GetFieldEvaluator("surface_energy_source")
        ->HasFieldChanged(S_next_.ptr(), name_);
    S_inter_->GetFieldEvaluator("surface_energy_source")
        ->HasFieldChanged(S_inter_.ptr(), name_);

    const Epetra_MultiVector& source0 =
        *S_inter_->GetFieldData("surface_energy_source")->ViewComponent("cell",false);
    const Epetra_MultiVector& source1 =
        *S_next_->GetFieldData("surface_energy_source")->ViewComponent("cell",false);

    int ncells = g_c.MyLength();
    for (int c=0; c!=ncells; ++c) {
      g_c[0][c] -= 0.5* (fa0[0][c] * source0[0][c] + fa1[0][c] * source1[0][c]);
    }
  }


  // external sources of mass, bringing enthalpy
  if (is_mass_source_term_) {
    const double& patm = *S->GetScalarData("atmospheric_pressure");

    // calculate enthalpy of incoming water
    double T_air1 = air_temp_->Compute(S_next_->time());
    double n1 = eos_liquid_->MolarDensity(T_air1, patm);
    double u1 = iem_liquid_->InternalEnergy(T_air1);
    double enth1 = u1 + patm / n1;
    
    // get enthalpy of outgoing water
    const Epetra_MultiVector& enth_surf = *S->GetFieldData(enthalpy_key_)
        ->ViewComponent("cell",false);
    
    // get the source in mols / s
    S_next_->GetFieldEvaluator("overland_source")
        ->HasFieldChanged(S_next_.ptr(), name_);
    S_inter_->GetFieldEvaluator("overland_source")
        ->HasFieldChanged(S_inter_.ptr(), name_);
    const Epetra_MultiVector& source0 =
        *S_inter_->GetFieldData("overland_source")->ViewComponent("cell",false);
    const Epetra_MultiVector& source1 =
        *S_next_->GetFieldData("overland_source")->ViewComponent("cell",false);

    // mass source done in cell_volume (rain falls straight down!)
    const Epetra_MultiVector& cv1 =
        *S_next_->GetFieldData("surface_cell_volume")->ViewComponent("cell",false);
    const Epetra_MultiVector& cv0 =
        *S_inter_->GetFieldData("surface_cell_volume")->ViewComponent("cell",false);

    int ncells = g_c.MyLength();
    for (int c=0; c!=ncells; ++c) {
      double molar_flux = 0.5 * (cv0[0][c] * source0[0][c] + cv1[0][c] * source1[0][c]);

      // upwind the enthalpy
      if (molar_flux > 0) {
        g_c[0][c] -= molar_flux * enth1;
      } else {
        g_c[0][c] -= molar_flux * enth_surf[0][c];
      }
    }
  }

  
  // conduction of heat to air
  if (is_air_conductivity_) {
    double T_air1 = air_temp_->Compute(S->time());
    const Epetra_MultiVector& temp =
        *S_inter_->GetFieldData(key_)->ViewComponent("cell",false);

    int ncells = g_c.MyLength();
    for (int c=0; c!=ncells; ++c) {
      g_c[0][c] -= K_surface_to_air_ * (T_air1 - temp[0][c]) * fa1[0][c];
    }
  }

  // coupling to subsurface
  if (coupled_to_subsurface_via_full_) {
    S->GetFieldEvaluator("overland_source_from_subsurface")
        ->HasFieldChanged(S.ptr(), name_);
    const Epetra_MultiVector& source1 =
        S->GetFieldData("overland_source_from_subsurface")->ViewComponent("cell",false);
    const Epetra_MultiVector& enth_surf =
        S->GetFieldData(enthalpy_key_)->ViewComponent("cell",false);
    const Epetra_MultiVector& enth_subsurf =
        S->GetFieldData("enthalpy")->ViewComponent("cell",false);

    AmanziMesh::Entity_ID_List cells;
    
    int ncells = g_c.MyLength();
    for (int c=0; c!=ncells; ++c) {
      double flux = source1[0][c];

      // upwind the enthalpy
      if (flux > 0.) { // exfiltration
        // get the subsurface's enthalpy
        AmanziMesh::Entity_ID f = mesh_->entity_get_parent(AmanziMesh::CELL, c);
        S->Mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
        ASSERT(cells.size() == 1);
        
        g_c[0][c] -= flux * enth_subsurf[0][cells[0]];
      } else { // infiltration
        g_c[0][c] -= flux * enth_surf[0][c];
      }
    }
  }
  
}

} // namespace
} // namespace

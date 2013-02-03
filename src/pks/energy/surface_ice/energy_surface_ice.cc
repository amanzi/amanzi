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
#include "function.hh"
#include "function-factory.hh"
#include "independent_variable_field_evaluator.hh"
#include "overland_source_from_subsurface_flux_evaluator.hh"

#include "energy_surface_ice.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory<EnergySurfaceIce> EnergySurfaceIce::reg_("surface energy");

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
EnergySurfaceIce::EnergySurfaceIce(Teuchos::ParameterList& plist,
        const Teuchos::RCP<TreeVector>& solution) :
    PKDefaultBase(plist, solution),
    EnergyBase(plist, solution),
    standalone_mode_(false),
    is_energy_source_term_(false),
    is_mass_source_term_(false),
    is_air_conductivity_(false),
    coupled_to_subsurface_via_full_(false) {

  plist_.set("primary variable key", "surface_temperature");
  plist_.set("domain name", "surface");
}


void EnergySurfaceIce::setup(const Teuchos::Ptr<State>& S) {
  // set up the meshes
  if (!S->HasMesh("surface")) {
    Teuchos::RCP<const AmanziMesh::Mesh> domain = S->GetMesh();
    ASSERT(domain->space_dimension() == 2);
    standalone_mode_ = true;
    S->RegisterMesh("surface", domain);
  } else {
    standalone_mode_ = false;
  }

  EnergyBase::setup(S);
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
  enth_plist.set("include work term", false);
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
    is_mass_source_term_ = S->HasFieldEvaluator("overland_source");

    // If there is air temperature and a conductivity, this provides heat transfer.
    if (air_plist.isParameter("air to surface conductivity")) {
      is_air_conductivity_ = true;
      K_surface_to_air_ = air_plist.get<double>("air to surface conductivity");
    }

  } else {
    // Make sure there is no mass source.
    if (S->HasFieldEvaluator("overland_source")) {
      Errors::Message m("Overland Energy equation found a mass source without an air temperature.");
      Exceptions::amanzi_throw(m);
    }
  }

  // coupling to subsurface
  coupled_to_subsurface_via_full_ =
      plist_.get<bool>("coupled to subsurface via full coupler", false);
  if (coupled_to_subsurface_via_full_) {
    // -- kill the preconditioner and replace with a TPFA precon
    Teuchos::ParameterList mfd_pc_plist = plist_.sublist("Diffusion PC");
    Teuchos::RCP<Operators::Matrix> precon =
        Teuchos::rcp(new Operators::MatrixMFD_TPFA(mfd_pc_plist, mesh_));
    set_preconditioner(precon);

    // -- ensure mass source from subsurface exists
    S->RequireFieldEvaluator("overland_source_from_subsurface");

    // -- energy source term from subsurface
    S->RequireField("overland_energy_source_from_subsurface")
        ->SetMesh(mesh_)->AddComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList source_plist =
        plist_.sublist("source from subsurface evaluator");
    if (!plist_.isParameter("source key"))
      source_plist.set("source key", "overland_energy_source_from_subsurface");
    if (!plist_.isParameter("flux key"))
      source_plist.set("flux key", "energy_flux");
    if (!plist_.isParameter("pressure key"))
      source_plist.set("pressure key", "temperature");
    source_plist.set("volume basis", false);

    Teuchos::RCP<Relations::OverlandSourceFromSubsurfaceFluxEvaluator>
        source_evaluator = Teuchos::rcp(
            new Relations::OverlandSourceFromSubsurfaceFluxEvaluator(source_plist));
    S->SetFieldEvaluator("overland_energy_source_from_subsurface", source_evaluator);
  }

  // Many quantities are based upon face areas, which are not the cell volume,
  // as the surface mesh has been flattened.
  if (!standalone_mode_) {
    S->RequireFieldEvaluator("surface_3d_cell_volume");
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
  Epetra_MultiVector& enth_f = *enth->ViewComponent("face",false);
  const Epetra_MultiVector& pres = *S->GetFieldData("surface_pressure")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& temp = *S->GetFieldData("surface_temperature")
      ->ViewComponent("face",false);

  AmanziMesh::Entity_ID_List cells;
  int nfaces = enth_f.MyLength();
  for (int f=0; f!=nfaces; ++f) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    if (cells.size() == 1) {
      double T = bc_markers_[f] == Operators::MATRIX_BC_DIRICHLET ?
          bc_values_[f] : temp[0][f];
      double p = pres[0][cells[0]];
      double dens = eos_liquid_->MolarDensity(T,p);
      double int_energy = iem_liquid_->InternalEnergy(T);
      double enthalpy = int_energy + p/dens;

      enth_f[0][f] = enthalpy * fabs(flux[0][f]);
    }
  }
}


// -------------------------------------------------------------
// Deal with the many source terms.
// -------------------------------------------------------------
void EnergySurfaceIce::AddSources_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g) {
  Epetra_MultiVector& g_c = *g->ViewComponent("cell",false);

  // Note that energy sources are based upon face areas, which are not cell volumes.
  Teuchos::RCP<const Epetra_MultiVector> fa0;
  Teuchos::RCP<const Epetra_MultiVector> fa1;
  if (standalone_mode_) {
    fa0 = S_->GetFieldData(cell_vol_key_)->ViewComponent("cell",false);
    fa1 = S_next_->GetFieldData(cell_vol_key_)->ViewComponent("cell",false);
  } else {
    fa0 = S_->GetFieldData("surface_3d_cell_volume")
        ->ViewComponent("cell",false);
    fa1 = S_next_->GetFieldData("surface_3d_cell_volume")
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
      g_c[0][c] -= 0.5* ((*fa0)[0][c] * source0[0][c] + (*fa1)[0][c] * source1[0][c]);
    }
  }

  // external sources of mass, bringing enthalpy
  if (is_mass_source_term_) {
    const double& patm = *S->GetScalarData("atmospheric_pressure");

    // calculate enthalpy of incoming water
    double time = S_next_->time();
    double T_air1 = air_temp_->operator()(&time);
    double n1 = eos_liquid_->MolarDensity(T_air1, patm);
    double u1 = iem_liquid_->InternalEnergy(T_air1);
    //    double enth1 = u1 + patm / n1;
    double enth1 = u1;

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

    // External source term is in [m water / s], not in [mols / s].  We assume
    // it comes in at the surface density, i.e. the surface temperature.  This
    // may need to be changed.
    S_next_->GetFieldEvaluator("surface_molar_density_liquid")
        ->HasFieldChanged(S_next_.ptr(), name_);
    S_inter_->GetFieldEvaluator("surface_molar_density_liquid")
        ->HasFieldChanged(S_inter_.ptr(), name_);

    const Epetra_MultiVector& nliq0 =
        *S_inter_->GetFieldData("surface_molar_density_liquid")
        ->ViewComponent("cell",false);
    const Epetra_MultiVector& nliq1 =
        *S_next_->GetFieldData("surface_molar_density_liquid")
        ->ViewComponent("cell",false);

    int ncells = g_c.MyLength();
    for (int c=0; c!=ncells; ++c) {
      double molar_flux = 0.5 * (cv0[0][c] * source0[0][c] * nliq0[0][c] + cv1[0][c] * source1[0][c] * nliq1[0][c]);

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
    double time = S->time();
    double T_air1 = air_temp_->operator()(&time);
    const Epetra_MultiVector& temp =
        *S_inter_->GetFieldData(key_)->ViewComponent("cell",false);

    int ncells = g_c.MyLength();
    for (int c=0; c!=ncells; ++c) {
      g_c[0][c] -= K_surface_to_air_ * (T_air1 - temp[0][c]) * (*fa1)[0][c];
    }
  }

  // coupling to subsurface
  // -- two parts -- conduction and advection
  if (coupled_to_subsurface_via_full_) {
    S->GetFieldEvaluator("overland_source_from_subsurface")
        ->HasFieldChanged(S.ptr(), name_);
    S->GetFieldEvaluator("enthalpy")->HasFieldChanged(S.ptr(), name_);
    S->GetFieldEvaluator(enthalpy_key_)->HasFieldChanged(S.ptr(), name_);

    // -- conduction source
    const Epetra_MultiVector& e_source1 =
        *S->GetFieldData("overland_energy_source_from_subsurface")
        ->ViewComponent("cell",false);

    // -- advection source
    const Epetra_MultiVector& source1 =
        *S->GetFieldData("overland_source_from_subsurface")->ViewComponent("cell",false);
    const Epetra_MultiVector& enth_surf =
        *S->GetFieldData(enthalpy_key_)->ViewComponent("cell",false);
    const Epetra_MultiVector& enth_subsurf =
        *S->GetFieldData("enthalpy")->ViewComponent("cell",false);

    AmanziMesh::Entity_ID_List cells;

    int ncells = g_c.MyLength();
    for (int c=0; c!=ncells; ++c) {
      double flux = source1[0][c];

      // upwind the enthalpy
      if (flux > 0.) { // exfiltration
        // get the subsurface's enthalpy
        AmanziMesh::Entity_ID f = mesh_->entity_get_parent(AmanziMesh::CELL, c);
        S->GetMesh()->face_get_cells(f, AmanziMesh::USED, &cells);
        ASSERT(cells.size() == 1);

        g_c[0][c] -= flux * enth_subsurf[0][cells[0]] + e_source1[0][cells[0]];
      } else { // infiltration
        g_c[0][c] -= flux * enth_surf[0][c] + e_source1[0][cells[0]];
      }
    }
  }

}

void EnergySurfaceIce::AddSourcesToPrecon_(const Teuchos::Ptr<State>& S, double h) {
  std::vector<double>& Acc_cells = mfd_preconditioner_->Acc_cells();

  // Note that energy sources are based upon face areas, which are not cell volumes.
  Teuchos::RCP<const Epetra_MultiVector> fa0;
  Teuchos::RCP<const Epetra_MultiVector> fa1;
  if (standalone_mode_) {
    fa0 = S_->GetFieldData(cell_vol_key_)->ViewComponent("cell",false);
    fa1 = S_next_->GetFieldData(cell_vol_key_)->ViewComponent("cell",false);
  } else {
    fa0 = S_->GetFieldData("surface_3d_cell_volume")
        ->ViewComponent("cell",false);
    fa1 = S_next_->GetFieldData("surface_3d_cell_volume")
        ->ViewComponent("cell",false);
  }

  // precon from air conductivity
  if (is_air_conductivity_) {
    int ncells = Acc_cells.size();
    for (int c=0; c!=ncells; ++c) {
      Acc_cells[c] += K_surface_to_air_ * (*fa1)[0][c];
    }
  }

}

} // namespace
} // namespace

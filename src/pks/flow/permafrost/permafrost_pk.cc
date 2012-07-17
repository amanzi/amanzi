/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code. 
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
------------------------------------------------------------------------- */


#include "bdf1_time_integrator.hh"
#include "flow_bc_factory.hh"
#include "eos_factory.hh"
#include "wrm_factory.hh"

#include "upwind_cell_centered.hh"
#include "upwind_arithmetic_mean.hh"
#include "upwind_total_flux.hh"
#include "upwind_gravity_flux.hh"

#include "permafrost.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory<Permafrost> Permafrost::reg_("permafrost flow");


// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
Permafrost::Permafrost(Teuchos::ParameterList& flow_plist,
                       const Teuchos::RCP<State>& S,
                       const Teuchos::RCP<TreeVector>& solution) {
  flow_plist_ = flow_plist;
  solution_ = solution;

  // require fields

  // -- primary variable: pressure on both cells and faces, ghosted, with 1 dof
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  S->RequireField("pressure", "flow")->SetMesh(S->Mesh())->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);

  // -- secondary variables
  CompositeVectorFactory one_cell_owned_factory;
  one_cell_owned_factory.SetMesh(S->Mesh());
  one_cell_owned_factory.SetGhosted();
  one_cell_owned_factory.SetComponent("cell", AmanziMesh::CELL, 1);

  CompositeVectorFactory one_cell_factory;
  one_cell_factory.SetMesh(S->Mesh());
  one_cell_factory.SetGhosted();
  one_cell_factory.AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireField("darcy_flux", "flow")->SetMesh(S->Mesh())->SetGhosted()
                                ->SetComponent("face", AmanziMesh::FACE, 1);
  S->RequireField("darcy_velocity", "flow")->SetMesh(S->Mesh())->SetGhosted()
                                ->SetComponent("cell", AmanziMesh::CELL, 3);

  S->RequireField("saturation_liquid", "flow")->Update(one_cell_owned_factory);
  S->RequireField("density_liquid", "flow")->Update(one_cell_owned_factory);
  S->RequireField("molar_density_liquid", "flow")->Update(one_cell_owned_factory);
  S->RequireField("viscosity_liquid", "flow")->Update(one_cell_owned_factory);

  S->RequireField("saturation_gas", "flow")->Update(one_cell_owned_factory);
  S->RequireField("density_gas", "flow")->Update(one_cell_owned_factory);
  S->RequireField("molar_density_gas", "flow")->Update(one_cell_owned_factory);
  S->RequireField("mol_frac_gas", "flow")->Update(one_cell_owned_factory);

  S->RequireField("saturation_ice", "flow")->Update(one_cell_owned_factory);
  S->RequireField("density_ice", "flow")->Update(one_cell_owned_factory);
  S->RequireField("molar_density_ice", "flow")->Update(one_cell_owned_factory);

  // -- For now, we assume scalar permeability.  This will change.
  S->RequireField("permeability", "flow")->Update(one_cell_owned_factory);
  S->RequireScalar("atmospheric_pressure", "flow");

  // -- independent variables not owned by this PK
  S->RequireField("cell_volume")->Update(one_cell_factory);
  S->RequireField("porosity")->Update(one_cell_factory);
  S->RequireField("temperature")->Update(one_cell_factory);
  S->RequireGravity();

  // -- work vectors
  S->RequireField("relative_permeability", "flow")->SetMesh(S->Mesh())->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);
  S->RequireField("numerical_rel_perm", "flow")->SetMesh(S->Mesh())->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);
  S->GetField("numerical_rel_perm","flow")->set_io_vis(false);

  // abs perm tensor
  variable_abs_perm_ = false; // currently not implemented, but may eventually want a model
  int c_owned = S->Mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  K_ = Teuchos::rcp(new std::vector<WhetStone::Tensor>(c_owned));
  for (int c=0; c!=c_owned; ++c) {
    (*K_)[c].init(S->Mesh()->space_dimension(),1);
  }

  // constitutive relations
  FlowRelations::EOSFactory eos_factory;

  // -- liquid eos
  Teuchos::ParameterList water_eos_plist = flow_plist_.sublist("Water EOS");
  eos_liquid_ = eos_factory.createEOS(water_eos_plist);

  // -- ice eos
  Teuchos::ParameterList ice_eos_plist = flow_plist_.sublist("Ice EOS");
  eos_ice_ = eos_factory.createEOS(ice_eos_plist);

  // -- gas eos
  Teuchos::ParameterList eos_gas_plist = flow_plist_.sublist("Vapor and Gas EOS");
  // EOSVaporInGas has a different interface than EOS; these sorts of
  // EOS need a factory/etc too.
  eos_gas_ = Teuchos::rcp(new FlowRelations::EOSVaporInGas(eos_gas_plist));

  // -- pc_il model.  cap pressure for ice-water interfaces are a
  // -- function of temperature.  See notes eqn 9
  Teuchos::ParameterList pc_il_plist =flow_plist_.sublist("Capillary Pressure Ice-Liquid");
  pc_ice_liq_model_ = Teuchos::rcp(new FlowRelations::PCIceWater(pc_il_plist));

  // -- water retention model
  Teuchos::ParameterList wrm_plist = flow_plist_.sublist("Water Retention Models");
  // count the number of region-model pairs
  int wrm_count = 0;
  for (Teuchos::ParameterList::ConstIterator i=wrm_plist.begin(); i!=wrm_plist.end(); ++i) {
    if (wrm_plist.isSublist(wrm_plist.name(i))) {
      wrm_count++;
    } else {
      std::string message("Permafrost: frozen water retention model list contains an entry that is not a sublist.");
      Exceptions::amanzi_throw(message);
    }
  }
  wrm_.resize(wrm_count);

  // instantiate the region-model pairs
  FlowRelations::WRMFactory wrm_factory;
  int iblock = 0;
  for (Teuchos::ParameterList::ConstIterator i=wrm_plist.begin(); i!=wrm_plist.end(); ++i) {
    Teuchos::ParameterList wrm_region_list = wrm_plist.sublist(wrm_plist.name(i));
    std::string region = wrm_region_list.get<std::string>("Region");
    wrm_[iblock] = Teuchos::rcp(new WRMRegionPair(region, wrm_factory.createWRM(wrm_region_list)));
    iblock++;
  }

  // boundary conditions
  Teuchos::ParameterList bc_plist = flow_plist_.sublist("boundary conditions", true);
  FlowBCFactory bc_factory(S->Mesh(), bc_plist);
  bc_pressure_ = bc_factory.CreatePressure();
  bc_flux_ = bc_factory.CreateMassFlux();

  // Relative permeability method
  string method_name = flow_plist_.get<string>("Relative permeability method", "Upwind with gravity");
  bool symmetric = false;
  if (method_name == "Upwind with gravity") {
    Krel_method_ = FLOW_RELATIVE_PERM_UPWIND_GRAVITY;
  } else if (method_name == "Cell centered") {
    Krel_method_ = FLOW_RELATIVE_PERM_CENTERED;
    symmetric = true;
  } else if (method_name == "Upwind with Darcy flux") {
    Krel_method_ = FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX;
  } else if (method_name == "Arithmetic mean") {
    Krel_method_ = FLOW_RELATIVE_PERM_ARITHMETIC_MEAN;
  }

  // operator for the diffusion terms
  Teuchos::ParameterList mfd_plist = flow_plist_.sublist("Diffusion");
  matrix_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_plist, S->Mesh()));
  matrix_->SetSymmetryProperty(symmetric);
  matrix_->SymbolicAssembleGlobalMatrices();

  // preconditioner for the NKA system
  Teuchos::ParameterList mfd_pc_plist = flow_plist_.sublist("Diffusion PC");
  preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_pc_plist, S->Mesh()));
  preconditioner_->SetSymmetryProperty(symmetric);
  preconditioner_->SymbolicAssembleGlobalMatrices();
  Teuchos::ParameterList mfd_pc_ml_plist = mfd_pc_plist.sublist("ML Parameters");
  preconditioner_->InitMLPreconditioner(mfd_pc_ml_plist);
};


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void Permafrost::initialize(const Teuchos::RCP<State>& S) {
  // initial timestep size
  dt_ = flow_plist_.get<double>("Initial time step", 1.);

  // initialize boundary conditions
  int nfaces = S->Mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::MFD_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  // update face pressures as a hint?
  Teuchos::RCP<CompositeVector> pres = S->GetFieldData("pressure", "flow");
  DeriveFaceValuesFromCellValues_(S, pres);

  // declare secondary variables initialized, as they will be done by
  // the commit_state call
  S->GetField("saturation_liquid","flow")->set_initialized();
  S->GetField("density_liquid","flow")->set_initialized();
  S->GetField("molar_density_liquid","flow")->set_initialized();
  S->GetField("viscosity_liquid","flow")->set_initialized();

  S->GetField("saturation_gas","flow")->set_initialized();
  S->GetField("density_gas","flow")->set_initialized();
  S->GetField("molar_density_gas","flow")->set_initialized();
  S->GetField("mol_frac_gas","flow")->set_initialized();

  S->GetField("saturation_ice","flow")->set_initialized();
  S->GetField("density_ice","flow")->set_initialized();
  S->GetField("molar_density_ice","flow")->set_initialized();

  S->GetField("relative_permeability","flow")->set_initialized();
  S->GetField("darcy_flux", "flow")->set_initialized();
  S->GetField("darcy_velocity", "flow")->set_initialized();

  // rel perm is special -- if the mode is symmetric, it needs to be
  // initialized to 1
  S->GetFieldData("numerical_rel_perm","flow")->PutScalar(1.0);
  S->GetField("numerical_rel_perm","flow")->set_initialized();

  // abs perm
  SetAbsolutePermeabilityTensor_(S);

  // operators
  matrix_->CreateMFDmassMatrices(K_.ptr());
  preconditioner_->CreateMFDmassMatrices(K_.ptr());

  // initialize the timesteppper
  solution_->set_data(pres);
  atol_ = flow_plist_.get<double>("Absolute error tolerance",1.0);
  rtol_ = flow_plist_.get<double>("Relative error tolerance",1.0);

  if (!flow_plist_.get<bool>("Strongly Coupled PK", false)) {
    // -- instantiate time stepper
    Teuchos::RCP<Teuchos::ParameterList> bdf1_plist_p =
      Teuchos::rcp(new Teuchos::ParameterList(flow_plist_.sublist("Time integrator")));
    time_stepper_ = Teuchos::rcp(new BDF1TimeIntegrator(this, bdf1_plist_p, solution_));
    time_step_reduction_factor_ = bdf1_plist_p->get<double>("time step reduction factor");

    // -- initialize time derivative
    Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp( new TreeVector(*solution_));
    solution_dot->PutScalar(0.0);

    // -- set initial state
    time_stepper_->set_initial_state(S->time(), solution_, solution_dot);
  }
};

} // namespace
} // namespace

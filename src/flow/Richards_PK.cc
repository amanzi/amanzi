/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (nnc@lanl.gov), 
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Usage: Richards_PK FPK(parameter_list, state);
         FPK->InitPK();
         FPK->Initialize();  // optional
         FPK->InitSteadyState(T, dT);
         FPK->InitTransient(T, dT);
*/

#include <vector>

#include "Epetra_IntVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "dbc.hh"
#include "exceptions.hh"

#include "Mesh.hh"
#include "Point.hh"
#include "mfd3d_diffusion.hh"

#include "Matrix_TPFA.hh"
#include "Flow_BC_Factory.hh"
#include "Richards_PK.hh"


namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* We set up only default values and call Init() routine to complete
* each variable initialization
****************************************************************** */
Richards_PK::Richards_PK(Teuchos::ParameterList& glist, Teuchos::RCP<State> S)
{
  // initialize pointers (Do we need smart pointers here? lipnikov@lanl.gov)
  bdf1_dae = NULL;

  Flow_PK::Init(glist, S);

  // extract two critical sublists
  Teuchos::ParameterList flow_list;
  if (glist.isSublist("Flow")) {
    flow_list = glist.sublist("Flow");
  } else {
    Errors::Message msg("Flow PK: input parameter list does not have <Flow> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  if (flow_list.isSublist("Richards Problem")) {
    rp_list_ = flow_list.sublist("Richards Problem");
  } else {
    Errors::Message msg("Flow PK: input parameter list does not have <Richards Problem> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  // miscalleneous default parameters
  ti_specs = NULL;
  block_picard = 1;
  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;

  mfd3d_method_ = FLOW_MFD3D_OPTIMIZED;
  mfd3d_method_preconditioner_ = FLOW_MFD3D_OPTIMIZED;

  src_sink = NULL;
  src_sink_distribution = 0;
  experimental_solver_ = FLOW_SOLVER_NKA;
}


/* ******************************************************************
* Clean memory.
****************************************************************** */
Richards_PK::~Richards_PK()
{
  delete bc_pressure;
  delete bc_flux;
  delete bc_head;
  delete bc_seepage;

  if (ti_specs != NULL) OutputTimeHistory(rp_list_, ti_specs->dT_history);
}


/* ******************************************************************
* Extract information from Richards Problem parameter list.
****************************************************************** */
void Richards_PK::InitPK()
{
  // Fundamental physical quantities
  double* gravity_data;
  S_->GetConstantVectorData("gravity")->ExtractView(&gravity_data);
  gravity_.init(dim);
  for (int k = 0; k < dim; k++) gravity_[k] = gravity_data[k];

  // Other physical quantaties
  rho_ = *(S_->GetScalarData("fluid_density"));
  mu_ = *(S_->GetScalarData("fluid_viscosity"));

  // Allocate memory for boundary data. It must go first.
  bc_tuple zero = {0.0, 0.0};
  bc_values.resize(nfaces_wghost, zero);
  bc_model.resize(nfaces_wghost, 0);
  bc_submodel.resize(nfaces_wghost, 0);

  rainfall_factor.resize(nfaces_owned, 1.0);

  // Read flow list and populate various structures. 
  ProcessParameterList(rp_list_);

  // Create water retention models
  Teuchos::ParameterList& wrm_list = rp_list_.sublist("Water retention models");
  rel_perm = Teuchos::rcp(new RelativePermeability(mesh_, wrm_list));
  rel_perm->Init(atm_pressure_, S_);
  rel_perm->ProcessParameterList();
  rel_perm->PopulateMapC2MB();
  rel_perm->set_experimental_solver(experimental_solver_);

  std::string krel_method_name = rp_list_.get<string>("relative permeability");
  rel_perm->ProcessStringRelativePermeability(krel_method_name);

  // Create the solution (pressure) vector.
  solution = Teuchos::rcp(new CompositeVector(*(S_->GetFieldData("pressure"))));
  
  const Epetra_BlockMap& cmap_owned = mesh_->cell_map(false);
  pdot_cells_prev = Teuchos::rcp(new Epetra_Vector(cmap_owned));
  pdot_cells = Teuchos::rcp(new Epetra_Vector(cmap_owned));

  // Initialize times.
  double time = S_->time();
  if (time >= 0.0) T_physics = time;

  // Initialize actions on boundary condtions. 
  ProcessShiftWaterTableList(rp_list_);

  // Process other fundamental structures.
  K.resize(ncells_wghost);
  SetAbsolutePermeabilityTensor();
  rel_perm->CalculateKVectorUnit(K, gravity_);

  /// Initialize Transmisibillities and Gravity terms contribution
  if (experimental_solver_ == FLOW_SOLVER_NEWTON) {
    const Epetra_BlockMap& fmap_wghost = mesh_->face_map(true);
    Transmis_faces = Teuchos::rcp(new Epetra_Vector(fmap_wghost));
    Grav_term_faces = Teuchos::rcp(new Epetra_Vector(fmap_wghost));
    ComputeTransmissibilities(*Transmis_faces, *Grav_term_faces);
  }

  // Allocate memory for wells
  if (src_sink_distribution & Amanzi::Functions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
    Kxy = Teuchos::rcp(new Epetra_Vector(cmap_owned));
  }


  // Select a proper matrix class
  if (experimental_solver_ == FLOW_SOLVER_NEWTON) {
    matrix_ = Teuchos::rcp(new Matrix_TPFA(S_, rel_perm->Krel_faces_ptr(), Transmis_faces, Grav_term_faces));
    preconditioner_ = Teuchos::rcp(new Matrix_TPFA(S_, rel_perm->Krel_faces_ptr(), Transmis_faces, Grav_term_faces));
  } else {
    matrix_ = Teuchos::rcp(new Matrix_MFD(mesh_));
    preconditioner_ = Teuchos::rcp(new Matrix_MFD(mesh_));
  }
  matrix_->AddActionProperty(AmanziFlow::FLOW_MATRIX_ACTION_MATRIX);
  preconditioner_->AddActionProperty(AmanziFlow::FLOW_MATRIX_ACTION_PRECONDITIONER);

  //rhs = matrix_->rhs();  // import rhs from the matrix

  is_matrix_symmetric = SetSymmetryProperty();
  matrix_->SetSymmetryProperty(is_matrix_symmetric);
  matrix_->SymbolicAssembleGlobalMatrices();


  // initialize boundary and source data 
  CompositeVector& pressure = *S_->GetFieldData("pressure", passwd_);
  UpdateSourceBoundaryData(time, pressure);

  // injected water mass
  mass_bc = 0.0;
}


/* ******************************************************************
* Initialization of auxiliary variables (lambda and two saturations).
* WARNING: Flow_PK may use complex initialization of the remaining 
* state variables.
****************************************************************** */
void Richards_PK::InitializeAuxiliaryData()
{
  // pressures
  CompositeVector& pressure = *S_->GetFieldData("pressure", passwd_);
  const Epetra_MultiVector& p_cells = *pressure.ViewComponent("cell");
  Epetra_MultiVector& p_faces = *pressure.ViewComponent("face");

  DeriveFaceValuesFromCellValues(p_cells, p_faces);

  double time = T_physics;
  UpdateSourceBoundaryData(time, pressure);

  // saturations
  Epetra_MultiVector& ws = *S_->GetFieldData("water_saturation", passwd_)->ViewComponent("cell");
  Epetra_MultiVector& ws_prev = *S_->GetFieldData("prev_water_saturation", passwd_)->ViewComponent("cell");

  DeriveSaturationFromPressure(p_cells, ws);
  ws_prev = ws;
}


/* ******************************************************************
* Initial pressure is set to the pressure for fully saturated rock.
****************************************************************** */
void Richards_PK::InitializeSteadySaturated()
{ 
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *(vo_->os()) << "initializing with a saturated steady state..." << endl;
  }
  double T = S_->time();
  SolveFullySaturatedProblem(T, *solution, ti_specs->ls_specs_ini);
}


/* ******************************************************************
* Specific initialization of the initial pressure guess calculation.
****************************************************************** */
void Richards_PK::InitPicard(double T0)
{
  ti_specs = &ti_specs_igs_;
  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;
  error_control_ |= ti_specs->error_control_options;

  InitNextTI(T0, 0.0, ti_specs_igs_);

  // calculate initial guess: cleaning is required (lipnikov@lanl.gov)
  T_physics = ti_specs_igs_.T0;
  dT = ti_specs_igs_.dT0;
  AdvanceToSteadyState_Picard(ti_specs_igs_);

  const Epetra_MultiVector& p_cells = *S_->GetFieldData("pressure")->ViewComponent("cell");
  Epetra_MultiVector& ws = *S_->GetFieldData("water_saturation", passwd_)->ViewComponent("cell");
  Epetra_MultiVector& ws_prev = *S_->GetFieldData("prev_water_saturation", passwd_)->ViewComponent("cell");

  DeriveSaturationFromPressure(p_cells, ws);
  ws_prev = ws;
}


/* ******************************************************************
* Specific initialization of a steady state time integration phase.
****************************************************************** */
void Richards_PK::InitSteadyState(double T0, double dT0)
{
  if (ti_specs != NULL) OutputTimeHistory(rp_list_, ti_specs->dT_history);
  ti_specs = &ti_specs_sss_;

  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE +  // usually 1 [Pa]
                   FLOW_TI_ERROR_CONTROL_SATURATION;  // usually 1e-4;
  error_control_ |= ti_specs->error_control_options;

  InitNextTI(T0, dT0, ti_specs_sss_);

  // DEBUG
  // AdvanceToSteadyState();
  // CommitState(FS); WriteGMVfile(FS); exit(0);
}


/* ******************************************************************
* Specific initialization of a transient time integration phase.  
* WARNING: Initialization of lambda is done in MPC and may be 
* erroneous in pure transient mode.
****************************************************************** */
void Richards_PK::InitTransient(double T0, double dT0)
{
  if (ti_specs != NULL) OutputTimeHistory(rp_list_, ti_specs->dT_history);
  ti_specs = &ti_specs_trs_;

  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE +  // usually 1 [Pa]
                   FLOW_TI_ERROR_CONTROL_SATURATION;  // usually 1e-4
  error_control_ |= ti_specs->error_control_options;

  InitNextTI(T0, dT0, ti_specs_trs_);
}


/* ******************************************************************
* Generic initialization of a next time integration phase.
****************************************************************** */
void Richards_PK::InitNextTI(double T0, double dT0, TI_Specs& ti_specs)
{
  ResetPKtimes(T0, dT0);

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    LinearSolver_Specs& ls = ti_specs.ls_specs;
    Teuchos::OSTab tab = vo_->getOSTab();
    *(vo_->os()) << endl << "****************************************" << endl
                 << "New TI phase: " << ti_specs.ti_method_name.c_str() << endl
                 << "****************************************" << endl
                 << " start T=" << T0 / FLOW_YEAR << " [y], dT=" << dT0 << " [sec]" << endl
                 << " error control id=" << error_control_ << endl
                 << " preconditioner for nonlinear solver: " << ti_specs.preconditioner_name.c_str() << endl
                 << " sources distribution id=" << src_sink_distribution << endl;

    if (ti_specs.initialize_with_darcy) {
      LinearSolver_Specs& ls_ini = ti_specs.ls_specs_ini;
      *(vo_->os()) << " initial pressure solver: " << ls_ini.solver_name << endl;
      if (ti_specs.clip_saturation > 0.0) {
        *(vo_->os()) << "  clipping saturation at " << ti_specs.clip_saturation << " [-]" << endl;
      } else if (ti_specs.clip_pressure > -5 * FLOW_PRESSURE_ATMOSPHERIC) {
        *(vo_->os()) << "  clipping pressure at " << ti_specs.clip_pressure << " [Pa]" << endl;
      }
    }
  }

  // set up new preconditioner
  Teuchos::ParameterList& tmp_list = preconditioner_list_.sublist(ti_specs.preconditioner_name);
  string mfd3d_method_name = tmp_list.get<string>("discretization method", "monotone mfd");
  ProcessStringMFD3D(mfd3d_method_name, &mfd3d_method_preconditioner_); 

  // preconditioner_->DestroyPreconditioner();
  preconditioner_->SetSymmetryProperty(is_matrix_symmetric);
  preconditioner_->SymbolicAssembleGlobalMatrices();
  preconditioner_->InitPreconditioner(ti_specs.preconditioner_name, preconditioner_list_);

  // set up new time integration or solver
  std::string ti_method_name(ti_specs.ti_method_name);

  if (ti_specs.ti_method == FLOW_TIME_INTEGRATION_BDF1) {
    Teuchos::ParameterList bdf1_list = rp_list_.sublist(ti_method_name).sublist("BDF1");
    if (! bdf1_list.isSublist("VerboseObject"))
        bdf1_list.sublist("VerboseObject") = rp_list_.sublist("VerboseObject");

    Teuchos::RCP<CompositeVector> cv = S_->GetFieldData("pressure", passwd_);
    if (bdf1_dae == NULL) bdf1_dae = new BDF1_TI<CompositeVector, CompositeVectorSpace>(*this, bdf1_list, cv);
  }

  // initialize mass matrices
  SetAbsolutePermeabilityTensor();
  for (int c = 0; c < ncells_wghost; c++) K[c] *= rho_ / mu_;

  if (experimental_solver_ != FLOW_SOLVER_NEWTON) {
    matrix_->CreateMFDmassMatrices(mfd3d_method_, K);
    preconditioner_->CreateMFDmassMatrices(mfd3d_method_preconditioner_, K);
  }

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    int missed_tmp = missed_bc_faces_;
#ifdef HAVE_MPI
    mesh_->get_comm()->SumAll(&missed_tmp, &missed_bc_faces_, 1);
#endif

    int nokay = matrix_->nokay();
    int npassed = matrix_->npassed();

    Teuchos::OSTab tab = vo_->getOSTab();
    *(vo_->os()) << " discretization method (prec): " << mfd3d_method_name.c_str() << endl
                 << "  good and repaired matrices: " << nokay << " " << npassed << endl
                 << " assigned default (no-flow) BC to " << missed_bc_faces_ << " faces" << endl << endl;
  }

  // Well modeling
  if (src_sink_distribution & Amanzi::Functions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
    CalculatePermeabilityFactorInWell();
  }

  // initialize pressure and lambda
  *solution = *S_->GetFieldData("pressure", passwd_);
  Teuchos::RCP<Epetra_MultiVector> solution_cells = solution->ViewComponent("cell", false);
  Teuchos::RCP<Epetra_MultiVector> solution_faces = solution->ViewComponent("face", false);

  if (ti_specs.initialize_with_darcy) {
    // Get a hydrostatic solution consistent with b.c.
    SolveFullySaturatedProblem(T0, *solution, ti_specs.ls_specs_ini);

    if (ti_specs.clip_saturation > 0.0) {
      double pmin = FLOW_PRESSURE_ATMOSPHERIC;
      ClipHydrostaticPressure(pmin, ti_specs.clip_saturation, *solution_cells);
      DeriveFaceValuesFromCellValues(*solution_cells, *solution_faces);
    } else if (ti_specs.clip_pressure > -5 * FLOW_PRESSURE_ATMOSPHERIC) {
      ClipHydrostaticPressure(ti_specs.clip_pressure, *solution_cells);
      DeriveFaceValuesFromCellValues(*solution_cells, *solution_faces);
    }
    pressure = *solution_cells;
  }

  // initialize saturation
  Epetra_Vector& ws = FS->ref_water_saturation();
  DeriveSaturationFromPressure(pressure, ws);
 
  // re-initialize lambda (experimental)
  if (ti_specs.pressure_lambda_constraints && experimental_solver_ == FLOW_SOLVER_NKA) {
    double Tp = T0 + dT0;
    EnforceConstraints_MFD(Tp, *solution);
  } else if (experimental_solver_ == FLOW_SOLVER_NEWTON){
    double Tp = T0 + dT0;
    Epetra_Vector* u_cells = FS->CreateCellView(*solution);
    Epetra_Vector* u_faces = FS->CreateFaceView(*solution);
    UpdateSourceBoundaryData(Tp, *u_cells, *u_faces);
    rel_perm->Compute(*solution, bc_model, bc_values);
  }

  // nonlinear solver control options
  ti_specs.num_itrs = 0;
  block_picard = 0;

  Epetra_Vector& Krel_faces = rel_perm->Krel_faces();
  Epetra_Vector& flux = FS->ref_darcy_flux();

  if (experimental_solver_ != FLOW_SOLVER_NEWTON) {
    matrix_->CreateMFDstiffnessMatrices(*rel_perm);  // We remove dT from mass matrices.
    matrix_->DeriveDarcyMassFlux(*solution, *face_importer_, bc_model, bc_values, flux);
    AddGravityFluxes_DarcyFlux(K, flux, *rel_perm);
  } else {
    matrix_->DeriveDarcyMassFlux(*solution_cells, *face_importer_, bc_model, bc_values, flux);
  }

  for (int f = 0; f < nfaces_owned; f++) flux[f] /= rho_;

  // rel_perm->Compute(*solution, bc_model, bc_values);
}


/* ******************************************************************
* this routine avoid limitations of MPC by bumping the time step.                                          
****************************************************************** */
double Richards_PK::CalculateFlowDt()
{
  if (ti_specs->ti_method == FLOW_TIME_INTEGRATION_PICARD && block_picard == 1) dT *= 1e+4;
  return dT;
}


/* ******************************************************************* 
* Performs one time step of size dT_MPC either for steady-state or 
* transient calculations.
* Warning: BDF2 and BDF1 will merge eventually.
******************************************************************* */
int Richards_PK::Advance(double dT_MPC)
{
  // update the axiliary flow state
  FS->CopyMasterFace2GhostFace(FS->ref_darcy_flux(), FS_aux->ref_darcy_flux());

  dT = dT_MPC;
  double time = S_->time();
  if (time >= 0.0) T_physics = time;

  // predict water mass change during time step
  time = T_physics;
  if (ti_specs->num_itrs == 0) {  // initialization
    Epetra_Vector udot(*super_map_);
    // I do not know how to estimate du/dt in a robust manner.
    // ComputeUDot(time, *solution, udot);  

    if (ti_specs->ti_method == FLOW_TIME_INTEGRATION_BDF1) {
      bdf1_dae->set_initial_state(time, *solution, udot);

    } else if (ti_specs->ti_method == FLOW_TIME_INTEGRATION_PICARD) {
      if (flow_status_ == FLOW_STATUS_STEADY_STATE) {
        AdvanceToSteadyState(time, dT_MPC);
        block_picard = 1;  // We will wait for transient initialization.
      }
    }

    int ierr;
    update_precon(time, *solution, dT, ierr);
    ti_specs->num_itrs++;
  }

  if (ti_specs->ti_method == FLOW_TIME_INTEGRATION_BDF1) {
    bdf1_dae->bdf1_step(dT, *solution, dTnext);
    bdf1_dae->commit_solution(dT, *solution);
    bdf1_dae->write_bdf1_stepping_statistics();

    T_physics = bdf1_dae->most_recent_time();
  }

  // Calculate time derivative and 2nd-order solution approximation.
  // Estimate of a time step multiplier overrides the above estimate.
  if (ti_specs->dT_method == FLOW_DT_ADAPTIVE) {
    Epetra_Vector& pressure = FS->ref_pressure();  // pressure at t^n

    for (int c = 0; c < ncells_owned; c++) {
      (*pdot_cells)[c] = ((*solution)[c] - pressure[c]) / dT; 
      (*solution)[c] = pressure[c] + ((*pdot_cells_prev)[c] + (*pdot_cells)[c]) * dT / 2;
    }

    double err, dTfactor;
    err = AdaptiveTimeStepEstimate(&dTfactor);
    if (err > 0.0) throw 1000;  // fix (lipnikov@lan.gov)
    dTnext = std::min(dT_MPC * dTfactor, ti_specs->dTmax);
  }

  dt_tuple times(time, dT_MPC);
  ti_specs->dT_history.push_back(times);

  ti_specs->num_itrs++;

  return 0;
}


/* ******************************************************************
* Transfer part of the internal data needed by transport to the 
* flow state FS_MPC. MPC may request to populate the original FS.
* The consistency condition is improved by adjusting saturation while
* preserving its LED property.
****************************************************************** */
void Richards_PK::CommitState(Teuchos::RCP<Flow_State> FS_MPC)
{
  // save cell-based and face-based pressures 
  Epetra_Vector& pressure = FS_MPC->ref_pressure();
  pressure = *solution_cells;
  Epetra_Vector& lambda = FS_MPC->ref_lambda();
  lambda = *solution_faces;

  Epetra_Vector& ws = FS_MPC->ref_water_saturation();
  Epetra_Vector& ws_prev = FS_MPC->ref_prev_water_saturation();
  Epetra_Vector& Krel_faces = rel_perm->Krel_faces();

  ws_prev = ws;
  DeriveSaturationFromPressure(pressure, ws);

  // calculate Darcy flux as diffusive part + advective part.
  Epetra_Vector& flux = FS_MPC->ref_darcy_flux();
  if (experimental_solver_ != FLOW_SOLVER_NEWTON) {
    matrix_->CreateMFDstiffnessMatrices(*rel_perm);  // We remove dT from mass matrices.
    matrix_->DeriveDarcyMassFlux(*solution, *face_importer_, bc_model, bc_values, flux);
    AddGravityFluxes_DarcyFlux(K, flux, *rel_perm);
  } else {
    matrix_->DeriveDarcyMassFlux(*solution_cells, *face_importer_, bc_model, bc_values, flux);
  }

  for (int f = 0; f < nfaces_owned; f++) flux[f] /= rho_;

  // update time derivative
  *pdot_cells_prev = *pdot_cells;

  // update mass balance
  // ImproveAlgebraicConsistency(flux, ws_prev, ws);
  
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH && flow_status_ >= FLOW_STATUS_TRANSIENT_STATE) {
    Epetra_Vector& phi = FS_MPC->ref_porosity();
    double mass_bc_dT = WaterVolumeChangePerSecond(bc_model, flux) * rho_ * dT;

    mass_amanzi = 0.0;
    for (int c = 0; c < ncells_owned; c++) {
      mass_amanzi += ws[c] * rho_ * phi[c] * mesh_->cell_volume(c);
    }

    double mass_amanzi_tmp = mass_amanzi, mass_bc_tmp = mass_bc_dT;
    mesh_->get_comm()->SumAll(&mass_amanzi_tmp, &mass_amanzi, 1);
    mesh_->get_comm()->SumAll(&mass_bc_tmp, &mass_bc_dT, 1);

    mass_bc += mass_bc_dT;

    Teuchos::OSTab tab = vo_->getOSTab();
    *(vo_->os()) << "water mass=" << mass_amanzi 
                 << " [kg], total boundary flux=" << mass_bc << " [kg]" << endl;
  }

  dT = dTnext;
}


/* ******************************************************************
* Estimate du/dt from the pressure equation, a du/dt = g - A*u.
* There is no meaniful way to estimate du/dt due to the abrupt change 
* of function a. The routine is not used.
****************************************************************** */
/*
double Richards_PK::ComputeUDot(double T, const Epetra_Vector& u, Epetra_Vector& udot)
{
  AssembleMatrixMFD(u, T);
  double norm_udot = matrix_->ComputeNegativeResidual(u, udot);

  Epetra_Vector* udot_cells = FS->CreateCellView(udot);
  const Epetra_Vector& phi = FS->ref_porosity();
  Epetra_Vector dSdP(mesh_->cell_map(false));
  rel_perm->DerivedSdP(u, dSdP);
 
  for (int c = 0; c < ncells_owned; c++) {
    double volume = mesh_->cell_volume(c);
    if (dSdP[c] > 0.0) (*udot_cells)[c] /= volume * dSdP[c] * phi[c] * rho_;
  }

  Epetra_Vector* udot_faces = FS->CreateFaceView(udot);
  DeriveFaceValuesFromCellValues(*udot_cells, *udot_faces);

  return norm_udot;
}
*/


/* ******************************************************************
* Temporary convertion from double to tensor.                                               
****************************************************************** */
void Richards_PK::SetAbsolutePermeabilityTensor(std::vector<WhetStone::Tensor>& K)
{
 
  if (dim == 2) {
    const Epetra_Vector& permeability_x = *(*FS->permeability())(0);
    const Epetra_Vector& permeability_y = *(*FS->permeability())(1);

    const Epetra_BlockMap& cmap_wghost = mesh_->cell_map(true);
    Epetra_Vector perm_x_gh(cmap_wghost);
    Epetra_Vector perm_y_gh(cmap_wghost);

    FS->CopyMasterCell2GhostCell(permeability_x, perm_x_gh);
    FS->CopyMasterCell2GhostCell(permeability_y, perm_y_gh);    
    
    for (int c = 0; c < ncells_wghost; c++) {
      if (perm_x_gh[c] == perm_y_gh[c]) {
	K[c].init(dim, 1);
	K[c](0, 0) = perm_x_gh[c];
      } else {
	K[c].init(dim, 2);
	K[c](0, 0) = perm_x_gh[c];
	K[c](1, 1) = perm_y_gh[c];
      }
    }    
    
  } else if (dim == 3) {
    const Epetra_Vector& permeability_x = *(*FS->permeability())(0);
    const Epetra_Vector& permeability_y = *(*FS->permeability())(1);
    const Epetra_Vector& permeability_z = *(*FS->permeability())(2);
    
    const Epetra_BlockMap& cmap_wghost = mesh_->cell_map(true);
    Epetra_Vector perm_x_gh(cmap_wghost);
    Epetra_Vector perm_y_gh(cmap_wghost);
    Epetra_Vector perm_z_gh(cmap_wghost);

    FS->CopyMasterCell2GhostCell(permeability_x, perm_x_gh);
    FS->CopyMasterCell2GhostCell(permeability_y, perm_y_gh);        
    FS->CopyMasterCell2GhostCell(permeability_z, perm_z_gh);        

    for (int c = 0; c < ncells_wghost; c++) {
      if (perm_x_gh[c] == perm_y_gh[c]  && perm_y_gh[c] == perm_z_gh[c]) {
	K[c].init(dim, 1);
	K[c](0, 0) = perm_x_gh[c];
      } else {
	K[c].init(dim, 2);
	K[c](0, 0) = perm_x_gh[c];
	K[c](1, 1) = perm_y_gh[c];
	K[c](2, 2) = perm_z_gh[c];
      }
    }        
  }
}


/* ******************************************************************
* Adds time derivative to the cell-based part of MFD algebraic system.                                               
****************************************************************** */
void Richards_PK::AddTimeDerivative_MFD(
    Epetra_Vector& pressure_cells, double dT_prec, Matrix_MFD* matrix_operator)
{
  Epetra_Vector dSdP(mesh_->cell_map(false));
  rel_perm->DerivedSdP(pressure_cells, dSdP);

  const Epetra_Vector& phi = FS->ref_porosity();
  std::vector<double>& Acc_cells = matrix_operator->Acc_cells();
  std::vector<double>& Fc_cells = matrix_operator->Fc_cells();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c = 0; c < ncells; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = rho_ * phi[c] * dSdP[c] * volume / dT_prec;
    Acc_cells[c] += factor;
    Fc_cells[c] += factor * pressure_cells[c];
  }
}


/* ******************************************************************
* Estimate dT increase factor by comparing the 1st and 2nd order
* time approximations. 
* WARNING: it is implmented for transient phase only.
****************************************************************** */
double Richards_PK::AdaptiveTimeStepEstimate(double* dTfactor)
{
  double tol, error, error_max = 0.0;
  double dTfactor_cell;

  *dTfactor = 100.0;
  for (int c = 0; c < ncells_owned; c++) {
    error = fabs((*pdot_cells)[c] - (*pdot_cells_prev)[c]) * dT / 2;
    tol = ti_specs_trs_.rtol * fabs((*solution)[c]) + ti_specs_trs_.atol;

    dTfactor_cell = sqrt(tol / std::max(error, FLOW_DT_ADAPTIVE_ERROR_TOLERANCE));
    *dTfactor = std::min(*dTfactor, dTfactor_cell);

    error_max = std::max(error_max, error - tol);
  }

  *dTfactor *= FLOW_DT_ADAPTIVE_SAFETY_FACTOR;
  *dTfactor = std::min(*dTfactor, FLOW_DT_ADAPTIVE_INCREASE);
  *dTfactor = std::max(*dTfactor, FLOW_DT_ADAPTIVE_REDUCTION);

#ifdef HAVE_MPI
    double dT_tmp = *dTfactor;
    solution->Comm().MinAll(&dT_tmp, dTfactor, 1);  // find the global minimum
 
    double error_tmp = error_max;
    solution->Comm().MaxAll(&error_tmp, &error_max, 1);  // find the global maximum
#endif
  return error_max;
}


/* ******************************************************************
* Identify type of the matrix.
****************************************************************** */
bool Richards_PK::SetSymmetryProperty()
{
  int method = rel_perm->method();
  bool sym = false;
  if ((method == FLOW_RELATIVE_PERM_CENTERED) ||
      (method == FLOW_RELATIVE_PERM_AMANZI)) sym = true;
  if (experimental_solver_ == FLOW_SOLVER_NEWTON || 
      experimental_solver_ == FLOW_SOLVER_PICARD_NEWTON) sym = false;

  return sym;
}


/* ******************************************************************
* Saturation should be in exact balance with Darcy fluxes in order to
* have local extrema dimishing (LED) property for concentrations. 
* WARNING: we can enforce it strictly only in some cells.
****************************************************************** */
void Richards_PK::ImproveAlgebraicConsistency(const Epetra_Vector& flux, 
                                              const Epetra_Vector& ws_prev, Epetra_Vector& ws)
{
  // create a disctributed flux vector
  Epetra_Vector flux_d(mesh_->face_map(true));
  for (int f = 0; f < nfaces_owned; f++) flux_d[f] = flux[f];
  FS->CopyMasterFace2GhostFace(flux_d);

  // create a distributed saturation vector
  Epetra_Vector ws_d(mesh_->cell_map(true));
  for (int c = 0; c < ncells_owned; c++) ws_d[c] = ws[c];
  FS->CopyMasterCell2GhostCell(ws_d);

  WhetStone::MFD3D_Diffusion mfd(mesh_);
  const Epetra_Vector& phi = FS->ref_porosity();
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    // calculate min/max values
    double wsmin, wsmax;
    wsmin = wsmax = ws[c];
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      int c2 = mfd.cell_get_face_adj_cell(c, f);
      if (c2 >= 0) {
        wsmin = std::min(wsmin, ws[c2]);
        wsmax = std::max(wsmax, ws[c2]);
      }
    }

    // predict new saturation
    ws[c] = ws_prev[c];
    double factor = dT / (phi[c] * mesh_->cell_volume(c));
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      ws[c] -= factor * flux_d[f] * dirs[n]; 
    }

    // limit new saturation
    ws[c] = std::max(ws[c], wsmin);
    ws[c] = std::min(ws[c], wsmax);
  }
}


}  // namespace AmanziFlow
}  // namespace Amanzi


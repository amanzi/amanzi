/*
  This is the Transport component of Amanzi. 

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>
#include <vector>

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "errors.hh"
#include "tabular-function.hh"
#include "GMVMesh.hh"

#include "Explicit_TI_RK.hh"

#include "Transport_PK.hh"
#include "Reconstruction.hh"
#include "Dispersion.hh"
#include "DispersionMatrixFactory.hh"
#include "LinearOperatorFactory.hh"


namespace Amanzi {
namespace AmanziTransport {

/* ******************************************************************
* We set up minimum default values and call Init() routine to 
* complete initialization.
****************************************************************** */
Transport_PK::Transport_PK(Teuchos::ParameterList& glist, Teuchos::RCP<State> S,
                           std::vector<std::string>& component_names)
    : S_(S), component_names_(component_names)
{
  parameter_list = glist.sublist("Transport");
  preconditioners_list = glist.sublist("Preconditioners");
  solvers_list = glist.sublist("Solvers");
  nonlin_solvers_list = glist.sublist("Nonlinear solvers");

  dT = dT_debug = T_physics = 0.0;
  double time = S->initial_time();
  if (time >= 0.0) T_physics = time;

  dispersion_method = TRANSPORT_DISPERSION_METHOD_TPFA; 
  dispersion_preconditioner = "identity";

  internal_tests = 0;
  tests_tolerance = TRANSPORT_CONCENTRATION_OVERSHOOT;

  MyPID = 0;
  mesh_ = S->GetMesh();
  dim = mesh_->space_dimension();

  bc_scaling = 0.0;
  mass_tracer_exact = 0.0;  // Tracer is defined as species #0.

  name_ = "state";  //  state password

  // require state variables
  if (!S_->HasField("darcy_flux")) {
    S_->RequireField("darcy_flux", name_)->SetMesh(mesh_)->SetGhosted(true)
        ->SetComponent("face", AmanziMesh::FACE, 1);
  }
  if (!S_->HasField("total_component_concentration")) {
    std::vector<std::vector<std::string> > subfield_names(1);
    int ncomponents = component_names_.size();
    for (int i = 0; i != ncomponents; ++i) {
      subfield_names[0].push_back(component_names_[i]);
    }

    S_->RequireField("total_component_concentration", name_, subfield_names)->SetMesh(mesh_)
      ->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, ncomponents);
  }
}


/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
Transport_PK::~Transport_PK()
{ 
  delete vo_;
  for (int i=0; i<bcs.size(); i++) delete bcs[i]; 
}


/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
int Transport_PK::InitPK()
{
  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::USED);

  MyPID = mesh_->get_comm()->MyPID();

  // extract control parameters
  ProcessParameterList();

  // state pre-prosessing
  Teuchos::RCP<CompositeVector> cv1;
  S_->GetFieldData("darcy_flux", name_)->ScatterMasterToGhosted("face");
  cv1 = S_->GetFieldData("darcy_flux", name_);
  darcy_flux = cv1->ViewComponent("face", true);

  cv1 = S_->GetFieldData("total_component_concentration", name_);
  tcc = cv1->ViewComponent("cell", false);

  Teuchos::RCP<const CompositeVector> cv2;
  cv2 = S_->GetFieldData("water_saturation");
  ws = cv2->ViewComponent("cell", false);

  cv2 = S_->GetFieldData("prev_water_saturation");
  ws_prev = cv2->ViewComponent("cell", false);

  cv2 = S_->GetFieldData("porosity");
  phi = cv2->ViewComponent("cell", false);

  // memory for new components
  tcc_tmp = Teuchos::rcp(new CompositeVector(*(S_->GetFieldData("total_component_concentration"))));
  tcc_tmp->PutScalar(0.0);

  // upwind 
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  upwind_cell_ = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));
  downwind_cell_ = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));

  IdentifyUpwindCells();

  // advection block initialization
  current_component_ = -1;
  component_local_min_.resize(ncells_owned);
  component_local_max_.resize(ncells_owned);

  const Epetra_Map& cmap_owned = mesh_->cell_map(false);
  ws_subcycle_start = Teuchos::rcp(new Epetra_Vector(cmap_owned));
  ws_subcycle_end = Teuchos::rcp(new Epetra_Vector(cmap_owned));

  // reconstruction initialization
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  limiter_ = Teuchos::rcp(new Epetra_Vector(cmap_wghost));
  lifting.ResetField(mesh_, Teuchos::null);
  lifting.Init();

  // boundary conditions initialization
  double time = T_physics;
  for (int i = 0; i < bcs.size(); i++) {
    bcs[i]->Compute(time);
  }

  CheckInfluxBC();

  // source term initialization
  if (src_sink_distribution & Functions::TransportActions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
    Kxy = Teuchos::rcp(new Epetra_Vector(mesh_->cell_map(false)));
  }

  return 0;
}


/* *******************************************************************
* Estimation of the time step based on T.Barth (Lecture Notes   
* presented at VKI Lecture Series 1994-05, Theorem 4.2.2.       
* Routine must be called every time we update a flow field.
*
* Warning: Barth calculates influx, we calculate outflux. The methods
* are equivalent for divergence-free flows and gurantee EMP. Outflux 
* takes into account sinks and sources but preserves only positivity
* of an advected mass.
* ***************************************************************** */
double Transport_PK::CalculateTransportDt()
{
  // loop over faces and accumulate upwinding fluxes
  std::vector<double> total_outflux(ncells_wghost, 0.0);

  for (int f = 0; f < nfaces_wghost; f++) {
    int c = (*upwind_cell_)[f];
    if (c >= 0) total_outflux[c] += fabs((*darcy_flux)[0][f]);
  }

  // loop over cells and calculate minimal dT
  double vol, outflux, dT_cell;
  dT = dT_cell = TRANSPORT_LARGE_TIME_STEP;
  int cmin_dT = 0;
  for (int c = 0; c < ncells_owned; c++) {
    outflux = total_outflux[c];
    if (outflux) {
      vol = mesh_->cell_volume(c);
      dT_cell = vol * (*phi)[0][c] * std::min((*ws_prev)[0][c], (*ws)[0][c]) / outflux;
    }
    if (dT_cell < dT) {
      dT = dT_cell;
      cmin_dT = c;
    }
  }
  if (spatial_disc_order == 2) dT /= 2;

  // communicate global time step
  double dT_tmp = dT;
#ifdef HAVE_MPI
  const Epetra_Comm& comm = ws_prev->Comm();
  comm.MinAll(&dT_tmp, &dT, 1);
#endif

  // incorporate developers and CFL constraints
  dT = std::min(dT, dT_debug);
  dT *= cfl_;

  // print optional diagnostics using maximum cell id as the filter
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int cmin_dT_unique = (fabs(dT_tmp * cfl_ - dT) < 1e-6 * dT) ? cmin_dT : -1;
 
#ifdef HAVE_MPI
    int cmin_dT_tmp = cmin_dT_unique;
    comm.MaxAll(&cmin_dT_tmp, &cmin_dT_unique, 1);
#endif
    if (cmin_dT == cmin_dT_unique) {
      const AmanziGeometry::Point& p = mesh_->cell_centroid(cmin_dT);

      Teuchos::OSTab tab = vo_->getOSTab();
      *(vo_->os()) << "cell id=" << cmin_dT << " has smallest dT, (" << p[0] << ", " << p[1];
      if (p.dim() == 3) *(vo_->os()) << ", ", p[2];
      *(vo_->os()) << ")" << endl;
    }
  }
  return dT;
}


/* ******************************************************************* 
* Estimate returns last time step unless it is zero.     
******************************************************************* */
double Transport_PK::EstimateTransportDt()
{
  if (dT == 0.0) CalculateTransportDt();
  return dT;
}


/* ******************************************************************* 
* MPC will call this function to advance the transport state.
* Efficient subcycling requires to calculate an intermediate state of
* saturation only once, which leads to a leap-frog-type algorithm.
* 
* WARNIN: We cannot assume that dT_MPC equals to the global time step.
******************************************************************* */
void Transport_PK::Advance(double dT_MPC)
{
  // calculate stable time step dT
  double dT_shift = 0.0, dT_global = dT_MPC;
  double time = S_->intermediate_time();
  if (time >= 0.0) { 
    T_physics = time;
    dT_shift = S_->initial_time() - time;
    dT_global = S_->final_time() - S_->initial_time();
  }

  CalculateTransportDt();
  double dT_original = dT;  // advance routines override dT
  int interpolate_ws = (dT < dT_global) ? 1 : 0;

  // start subcycling
  double dT_sum = 0.0;
  double dT_cycle;
  if (interpolate_ws) {
    dT_cycle = dT_original;
    InterpolateCellVector(*ws_prev, *ws, dT_shift, dT_global, *ws_subcycle_start);
  } else {
    dT_cycle = dT_MPC;
    ws_start = ws_prev;
    ws_end = ws;
  }

  int ncycles = 0, swap = 1;
  while (dT_sum < dT_MPC) {
    // update boundary conditions
    time = T_physics;
    for (int i = 0; i < bcs.size(); i++) bcs[i]->Compute(time);
    
    double dT_try = dT_MPC - dT_sum;
    bool final_cycle = false;
    if (dT_try >= 2 * dT_original) {
      dT_cycle = dT_original;
    } else if (dT_try > dT_original) { 
      dT_cycle = dT_try / 2; 
    } else {
      dT_cycle = dT_try;
      final_cycle = true;
    }

    T_physics += dT_cycle;
    dT_sum += dT_cycle;

    if (interpolate_ws) {
      if (swap) {  // Initial water saturation is in 'start'.
        ws_start = ws_subcycle_start;
        ws_end = ws_subcycle_end;

        double dT_int = dT_sum + dT_shift;
        InterpolateCellVector(*ws_prev, *ws, dT_int, dT_global, *ws_subcycle_end);
      } else {  // Initial water saturation is in 'end'.
        ws_start = ws_subcycle_end;
        ws_end = ws_subcycle_start;

        double dT_int = dT_sum + dT_shift;
        InterpolateCellVector(*ws_prev, *ws, dT_int, dT_global, *ws_subcycle_start);
      }
      swap = 1 - swap;
    }

    if (spatial_disc_order == 1) {  // temporary solution (lipnikov@lanl.gov)
      AdvanceDonorUpwind(dT_cycle);
    } else if (spatial_disc_order == 2 && temporal_disc_order == 1) {
      AdvanceSecondOrderUpwindRK1(dT_cycle);
    } else if (spatial_disc_order == 2 && temporal_disc_order == 2) {
      AdvanceSecondOrderUpwindRK2(dT_cycle);
      // AdvanceSecondOrderUpwindGeneric(dT_cycle);
    }

    if (! final_cycle) {  // rotate concentrations
      Teuchos::RCP<Epetra_MultiVector> tcc_next = tcc_tmp->ViewComponent("cell", false);
      *tcc = *tcc_next;
    }

    ncycles++;
  }

  dT = dT_original;  // restore the original dT (just in case)

  // We define tracer as the species #0 as calculate some statistics.
  int number_components = tcc->NumVectors();
  Teuchos::RCP<Epetra_MultiVector> tcc_next = tcc_tmp->ViewComponent("cell", false);

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *(vo_->os()) << ncycles << " sub-cycles, dT_stable=" << dT_original 
                 << " [sec]  dT_MPC=" << dT_MPC << " [sec]" << endl;

    double tccmin_vec[number_components];
    double tccmax_vec[number_components];

    tcc_next->MinValue(tccmin_vec);
    tcc_next->MaxValue(tccmax_vec);

    double tccmin, tccmax;
    tcc_next->Comm().MinAll(tccmin_vec, &tccmin, 1);  // find the global extrema
    tcc_next->Comm().MaxAll(tccmax_vec, &tccmax, 1);  // find the global extrema

    mass_tracer_exact += TracerVolumeChangePerSecond(0) * dT_MPC;
    double mass_tracer = 0.0;
    for (int c = 0; c < ncells_owned; c++) {
      double vol = mesh_->cell_volume(c);
      mass_tracer += (*ws)[0][c] * (*phi)[0][c] * (*tcc_next)[0][c] * vol;
    }

    double mass_tracer_tmp = mass_tracer, mass_exact_tmp = mass_tracer_exact, mass_exact;
    mesh_->get_comm()->SumAll(&mass_tracer_tmp, &mass_tracer, 1);
    mesh_->get_comm()->SumAll(&mass_exact_tmp, &mass_exact, 1);

    double mass_loss = mass_exact - mass_tracer;
    *(vo_->os()) << "species #0: " << tccmin << " <= concentration <= " << tccmax << endl;
    *(vo_->os()) << "species #0: mass=" << mass_tracer << " [kg], mass left domain=" << mass_loss << " [kg]" << endl;
  }

  if (dispersion_models_.size() != 0) {
    std::string scheme("tpfa");
    if (mesh_->mesh_type() == AmanziMesh::RECTANGULAR) scheme = "mfd";

    DispersionMatrixFactory mfactory;
    dispersion_matrix_ = mfactory.Create(scheme, &dispersion_models_, mesh_, S_);

    dispersion_matrix_->CalculateDispersionTensor(*darcy_flux, *phi, *ws);
    dispersion_matrix_->SymbolicAssembleMatrix();
    dispersion_matrix_->ModifySymbolicAssemble();
    dispersion_matrix_->AssembleMatrix(*phi);
    dispersion_matrix_->AddTimeDerivative(dT_MPC, *phi, *ws);

    dispersion_matrix_->InitPreconditioner(dispersion_preconditioner, preconditioners_list);
    dispersion_matrix_->UpdatePreconditioner();

    AmanziSolvers::LinearOperatorFactory<Dispersion, Epetra_Vector, Epetra_BlockMap> sfactory;
    Teuchos::RCP<AmanziSolvers::LinearOperator<Dispersion, Epetra_Vector, Epetra_BlockMap> >
       solver = sfactory.Create(dispersion_solver, solvers_list, dispersion_matrix_);

    solver->add_criteria(AmanziSolvers::LIN_SOLVER_MAKE_ONE_ITERATION);  // Make at least one iteration

    const Epetra_Map& map = dispersion_matrix_->RangeMap();
    Epetra_Vector rhs(map), sol(map);

    double residual = 0.0;
    int num_itrs = 0;
    
    for (int i = 0; i < number_components; i++) {
      for (int c = 0; c < ncells_owned; c++) {
        double factor = mesh_->cell_volume(c) * (*ws)[0][c] * (*phi)[0][c] / dT_MPC;
        rhs[c] = (*tcc_next)[i][c] * factor;
        sol[c] = (*tcc_next)[i][c];
      }

      solver->ApplyInverse(rhs, sol);

      residual += solver->residual();
      num_itrs += solver->num_itrs();

      for (int c = 0; c < ncells_owned; c++) {
        (*tcc_next)[i][c] = sol[c];
      }
    }

    if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *(vo_->os()) << "dispersion solver (" << solver->name() 
                   << ") ||r||=" << residual / number_components
                   << " itrs=" << num_itrs / number_components << endl;
    }
  }

  // DEBUG
  // WriteGMVfile(TS_nextMPC); 
}


/* ******************************************************************* 
* Updates states 
******************************************************************* */
void Transport_PK::CommitState(Teuchos::RCP<State> S) 
{
  Teuchos::RCP<CompositeVector> cv;
  cv = S->GetFieldData("total_component_concentration", name_);
  *cv = *tcc_tmp;
}


/* ******************************************************************* 
 * A simple first-order transport method 
 ****************************************************************** */
void Transport_PK::AdvanceDonorUpwind(double dT_cycle)
{
  dT = dT_cycle;  // overwrite the maximum stable transport step
  mass_tracer_source = 0.0;

  // populating next state of concentrations
  tcc_tmp->ScatterMasterToGhosted("cell");
  Teuchos::RCP<Epetra_MultiVector> tcc_next = tcc_tmp->ViewComponent("cell", true);

  // prepare conservative state in master and slave cells
  double vol_phi_ws, tcc_flux;
  int ncomponents = tcc->NumVectors();

  for (int c = 0; c < ncells_owned; c++) {
    vol_phi_ws = mesh_->cell_volume(c) * (*phi)[0][c] * (*ws_start)[0][c];

    for (int i = 0; i < ncomponents; i++)
      (*tcc_next)[i][c] = (*tcc)[i][c] * vol_phi_ws;
  }

  // advance all components at once
  for (int f = 0; f < nfaces_wghost; f++) {  // loop over master and slave faces
    int c1 = (*upwind_cell_)[f];
    int c2 = (*downwind_cell_)[f];

    double u = fabs((*darcy_flux)[0][f]);

    if (c1 >=0 && c1 < ncells_owned && c2 >= 0 && c2 < ncells_owned) {
      for (int i = 0; i < ncomponents; i++) {
        tcc_flux = dT * u * (*tcc)[i][c1];
        (*tcc_next)[i][c1] -= tcc_flux;
        (*tcc_next)[i][c2] += tcc_flux;
      }

    } else if (c1 >=0 && c1 < ncells_owned && (c2 >= ncells_owned || c2 < 0)) {
      for (int i = 0; i < ncomponents; i++) {
        tcc_flux = dT * u * (*tcc)[i][c1];
        (*tcc_next)[i][c1] -= tcc_flux;
      }

    } else if (c1 >= ncells_owned && c2 >= 0 && c2 < ncells_owned) {
      for (int i = 0; i < ncomponents; i++) {
        tcc_flux = dT * u * (*tcc_next)[i][c1];
        (*tcc_next)[i][c2] += tcc_flux;
      }
    }
  }

  // loop over exterior boundary sets
  for (int n = 0; n < bcs.size(); n++) {
    int i = bcs_tcc_index[n];

    for (Functions::TransportBoundaryFunction::Iterator bc = bcs[n]->begin(); bc != bcs[n]->end(); ++bc) {
      int f = bc->first;
      int c2 = (*downwind_cell_)[f];

      if (c2 >= 0) {
        double u = fabs((*darcy_flux)[0][f]);
        tcc_flux = dT * u * bc->second;
        (*tcc_next)[i][c2] += tcc_flux;
      }
    }
  }

  // process external sources
  if (src_sink != NULL) {
    double time = T_physics;
    ComputeAddSourceTerms(time, dT, src_sink, *tcc_next);
  }

  // recover concentration from new conservative state
  for (int c = 0; c < ncells_owned; c++) {
    vol_phi_ws = mesh_->cell_volume(c) * (*phi)[0][c] * (*ws_end)[0][c];
    for (int i = 0; i < ncomponents; i++) (*tcc_next)[i][c] /= vol_phi_ws;
  }

  // update mass balance
  mass_tracer_exact += mass_tracer_source * dT;

  if (internal_tests) {
    CheckGEDproperty(*tcc_next);
  }
}


/* ******************************************************************* 
 * We have to advance each component independently due to different
 * reconstructions. We use tcc when only owned data are needed and 
 * tcc_next when owned and ghost data. This is a special routine for 
 * transient flow and uses first-order time integrator. 
 ****************************************************************** */
void Transport_PK::AdvanceSecondOrderUpwindRK1(double dT_cycle)
{
  dT = dT_cycle;  // overwrite the maximum stable transport step
  mass_tracer_source = 0.0;

  // work memory
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  Epetra_Vector f_component(cmap_wghost);

  // distribute vector of concentrations
  S_->GetFieldData("total_component_concentration")->ScatterMasterToGhosted("cell");
  Teuchos::RCP<Epetra_MultiVector> tcc_next = tcc_tmp->ViewComponent("cell", true);

  int ncomponents = tcc->NumVectors();
  for (int i = 0; i < ncomponents; i++) {
    current_component_ = i;  // needed by BJ 

    double T = T_physics;
    Epetra_Vector*& component = (*tcc)(i);
    fun(T, *component, f_component);

    double ws_ratio;
    for (int c = 0; c < ncells_owned; c++) {
      ws_ratio = (*ws_start)[0][c] / (*ws_end)[0][c];
      (*tcc_next)[i][c] = ((*tcc)[i][c] + dT * f_component[c]) * ws_ratio;
    }
  }

  // update mass balance
  mass_tracer_exact += mass_tracer_source * dT;

  if (internal_tests) {
    CheckGEDproperty(*tcc_next);
  }
}


/* ******************************************************************* 
 * We have to advance each component independently due to different
 * reconstructions. This is a special routine for transient flow and 
 * uses second-order predictor-corrector time integrator. 
 ****************************************************************** */
void Transport_PK::AdvanceSecondOrderUpwindRK2(double dT_cycle)
{
  dT = dT_cycle;  // overwrite the maximum stable transport step
  mass_tracer_source = 0.0;

  // work memory
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  Epetra_Vector f_component(cmap_wghost);

  // distribute old vector of concentrations
  S_->GetFieldData("total_component_concentration")->ScatterMasterToGhosted("cell");
  Teuchos::RCP<Epetra_MultiVector> tcc_next = tcc_tmp->ViewComponent("cell", true);

  Epetra_Vector ws_ratio(Copy, *ws_start, 0);
  for (int c = 0; c < ncells_owned; c++) ws_ratio[c] /= (*ws_end)[0][c];

  // predictor step
  for (int i = 0; i < tcc->NumVectors(); i++) {
    current_component_ = i;  // needed by BJ 

    double T = T_physics;
    Epetra_Vector*& component = (*tcc)(i);
    fun(T, *component, f_component);

    for (int c = 0; c < ncells_owned; c++) {
      (*tcc_next)[i][c] = ((*tcc)[i][c] + dT * f_component[c]) * ws_ratio[c];
    }
  }

  tcc_tmp->ScatterMasterToGhosted("cell");

  // corrector step
  for (int i = 0; i < tcc->NumVectors(); i++) {
    current_component_ = i;  // needed by BJ 

    double T = T_physics;
    Epetra_Vector*& component = (*tcc_next)(i);
    fun(T, *component, f_component);

    for (int c = 0; c < ncells_owned; c++) {
      double value = ((*tcc)[i][c] + dT * f_component[c]) * ws_ratio[c];
      (*tcc_next)[i][c] = ((*tcc_next)[i][c] + value) / 2;
    }
  }

  // update mass balance
  mass_tracer_exact += mass_tracer_source * dT / 2;

  if (internal_tests) {
    CheckGEDproperty(*tcc_next);
  }
}


/* ******************************************************************* 
 * We have to advance each component independently due to different
 * reconstructions. We use tcc when only owned data are needed and 
 * tcc_next when owned and ghost data.
 *
 * Data flow: loop over components and apply the generic RK2 method.
 * The generic means that saturation is constant during time step. 
 ****************************************************************** */
/*
void Transport_PK::AdvanceSecondOrderUpwindGeneric(double dT_cycle)
{
  dT = dT_cycle;  // overwrite the maximum stable transport step

  Teuchos::RCP<Epetra_MultiVector> tcc = TS->total_component_concentration();
  Teuchos::RCP<Epetra_MultiVector> tcc_next = TS_nextBIG->total_component_concentration();

  // define time integration method
  Explicit_TI::RK::method_t ti_method = Explicit_TI::RK::forward_euler;
  if (temporal_disc_order == 2) {
    ti_method = Explicit_TI::RK::heun_euler;
  }
  Explicit_TI::RK TVD_RK(*this, ti_method, *component_);

  int ncomponents = tcc->NumVectors();
  for (int i = 0; i < ncomponents; i++) {
    current_component_ = i;  // it is needed in BJ called inside RK:fun

    Epetra_Vector*& tcc_component = (*tcc)(i);
    TS_nextBIG->CopyMasterCell2GhostCell(*tcc_component, *component_);

    double T = 0.0;  // requires fixes (lipnikov@lanl.gov)
    TVD_RK.step(T, dT, *component_, *component_next_);

    for (int c = 0; c < ncells_owned; c++) (*tcc_next)[i][c] = (*component_next_)[c];
  }

  if (internal_tests) {
    CheckGEDproperty(*tcc_next);
  }
}
*/


/* ******************************************************************
* Computes source and sink terms and adds them to vector tcc.
* Return mass rate for the tracer.
****************************************************************** */
void Transport_PK::ComputeAddSourceTerms(double Tp, double dTp, 
                                         Functions::TransportDomainFunction* src_sink, 
                                         Epetra_MultiVector& tcc)
{
  int ncomponents = tcc.NumVectors();
  for (int i = 0; i < ncomponents; i++) {
    std::string name(component_names_[i]);
    
    if (src_sink_distribution & 
        Functions::TransportActions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY)
      src_sink->ComputeDistributeMultiValue(Tp, name, Kxy->Values()); 
    else
      src_sink->ComputeDistributeMultiValue(Tp, name, NULL);

    Functions::TransportDomainFunction::Iterator src;
    for (src = src_sink->begin(); src != src_sink->end(); ++src) {
      int c = src->first;
      double value = mesh_->cell_volume(c) * src->second;

      tcc[i][c] += dTp * value;
      if (i == 0) mass_tracer_source += value;
    }
  }
}


/* *******************************************************************
 * Identify flux direction based on orientation of the face normal 
 * and sign of the  Darcy velocity.                               
 ****************************************************************** */
void Transport_PK::IdentifyUpwindCells()
{
  for (int f = 0; f < nfaces_wghost; f++) {
    (*upwind_cell_)[f] = -1;  // negative value indicates boundary
    (*downwind_cell_)[f] = -1;
  }

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);

    for (int i = 0; i < faces.size(); i++) {
      int f = faces[i];
      if ((*darcy_flux)[0][f] * fdirs[i] >= 0) {
        (*upwind_cell_)[f] = c;
      } else {
        (*downwind_cell_)[f] = c;
      }
    }
  }
}


/* *******************************************************************
* Interpolate linearly in time between two values v0 and v1. The time 
* is measuared relative to value v0; so that v1 is at time dT. The
* interpolated data are at time dT_int.            
******************************************************************* */
void Transport_PK::InterpolateCellVector(
    const Epetra_MultiVector& v0, const Epetra_MultiVector& v1, 
    double dT_int, double dT, Epetra_MultiVector& v_int) 
{
  double a = dT_int / dT;
  double b = 1.0 - a;
  v_int.Update(b, v0, a, v1, 0.);
}

}  // namespace AmanziTransport
}  // namespace Amanzi


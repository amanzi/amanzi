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

#include "BCs.hh"
#include "errors.hh"
#include "Explicit_TI_RK.hh"
#include "GMVMesh.hh"
#include "LinearOperatorFactory.hh"
#include "Mesh.hh"
#include "OperatorDefs.hh"
#include "OperatorDiffusionFactory.hh"
#include "OperatorDiffusion.hh"

#include "Transport_PK.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* We set up minimum default values and call Construct_() to complete 
* the initialization process.
****************************************************************** */
Transport_PK::Transport_PK(Teuchos::ParameterList& glist, Teuchos::RCP<State> S, 
			   const std::string& pk_list_name,
                           std::vector<std::string>& component_names)
{
  Construct_(glist, S, pk_list_name, component_names);
}


/* ******************************************************************
* Constructor for Alquimia. 
****************************************************************** */
#ifdef ALQUIMIA_ENABLED
Transport_PK::Transport_PK(Teuchos::ParameterList& glist, Teuchos::RCP<State> S,
			   const std::string& pk_list_name,
                           Teuchos::RCP<AmanziChemistry::Chemistry_State> chem_state,
                           Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine)
    : chem_state_(chem_state), chem_engine_(chem_engine)
{
  // Retrieve the component names from the chemistry engine.
  std::vector<std::string> component_names;
  chem_engine_->GetPrimarySpeciesNames(component_names);
  Construct_(glist, S, pk_list_name, component_names);
}
#endif


/* ******************************************************************
* High-level initialization.
****************************************************************** */
void Transport_PK::Construct_(Teuchos::ParameterList& glist, 
                              Teuchos::RCP<State> S,
			      const std::string& pk_list_name,
                              std::vector<std::string>& component_names)
{
  vo_ = NULL;
  S_ = S;
  component_names_ = component_names;

  parameter_list = glist.sublist("PKs").sublist(pk_list_name);
  preconditioners_list = glist.sublist("Preconditioners");
  solvers_list = glist.sublist("Solvers");
  nonlin_solvers_list = glist.sublist("Nonlinear solvers");

  dT = dT_debug = T_physics = 0.0;
  double time = S->time();
  if (time >= 0.0) T_physics = time;

  dispersion_models_ = 0; 
  dispersion_preconditioner = "identity";

  internal_tests = 0;
  tests_tolerance = TRANSPORT_CONCENTRATION_OVERSHOOT;

  MyPID = 0;
  mesh_ = S->GetMesh();
  dim = mesh_->space_dimension();

  dT = 0.0;
  bc_scaling = 0.0;
  passwd_ = "state";  //  state password

  // require state variables when Flow is off
  if (!S->HasField("permeability")) {
    S->RequireField("permeability", passwd_)->SetMesh(mesh_)->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, dim);
  }
  if (!S->HasField("porosity")) {
    S->RequireField("porosity", passwd_)->SetMesh(mesh_)->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S->HasField("darcy_flux")) {
    S->RequireField("darcy_flux", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);
  }
  if (!S->HasField("water_saturation")) {
    S->RequireField("water_saturation", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S->HasField("prev_water_saturation")) {
    S->RequireField("prev_water_saturation", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // require state variables when Transport is on
  if (!S->HasField("total_component_concentration")) {
    std::vector<std::vector<std::string> > subfield_names(1);
    int ncomponents = component_names_.size();
    for (int i = 0; i != ncomponents; ++i) {
      subfield_names[0].push_back(component_names_[i]);
    }

    S->RequireField("total_component_concentration", passwd_, subfield_names)->SetMesh(mesh_)
      ->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, ncomponents);
  }
}


/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
Transport_PK::~Transport_PK()
{ 
  if (vo_ != NULL) {
    delete vo_;
  }
  for (int i=0; i<bcs.size(); i++) delete bcs[i]; 
  for (int i=0; i<srcs.size(); i++) delete srcs[i]; 
}


/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
void Transport_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  // Check that stars are oriented favorably for transport PK.
  Policy(S);

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
  S->GetFieldData("darcy_flux", passwd_)->ScatterMasterToGhosted("face");
  cv1 = S->GetFieldData("darcy_flux", passwd_);
  darcy_flux = cv1->ViewComponent("face", true);

  Teuchos::RCP<const CompositeVector> cv2;
  cv2 = S->GetFieldData("water_saturation");
  ws = cv2->ViewComponent("cell", false);

  cv2 = S->GetFieldData("prev_water_saturation");
  ws_prev = cv2->ViewComponent("cell", false);

  cv2 = S->GetFieldData("porosity");
  phi = cv2->ViewComponent("cell", false);

  tcc = S->GetFieldData("total_component_concentration", passwd_);

  // memory for new components
  tcc_tmp = Teuchos::rcp(new CompositeVector(*(S->GetFieldData("total_component_concentration"))));
  *tcc_tmp = *tcc;

  // upwind 
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  upwind_cell_ = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));
  downwind_cell_ = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));

  IdentifyUpwindCells();

  // advection block initialization
  current_component_ = -1;

  const Epetra_Map& cmap_owned = mesh_->cell_map(false);
  ws_subcycle_start = Teuchos::rcp(new Epetra_Vector(cmap_owned));
  ws_subcycle_end = Teuchos::rcp(new Epetra_Vector(cmap_owned));

  // reconstruction initialization
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  lifting_ = Teuchos::rcp(new Operators::ReconstructionCell(mesh_));

  // boundary conditions initialization
  double time = T_physics;
  for (int i = 0; i < bcs.size(); i++) {
    bcs[i]->Compute(time);
  }

  VV_CheckInfluxBC();

  // source term initialization
  for (int i =0; i < srcs.size(); i++) {
    int distribution = srcs[i]->CollectActionsList();
    if (distribution & CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      Kxy = Teuchos::rcp(new Epetra_Vector(mesh_->cell_map(false)));
      Errors::Message msg;
      msg << "Transport PK: missing capability: permeability-based source distribution.\n";
      Exceptions::amanzi_throw(msg);  
      break;
    }
  }
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
  S_->GetFieldData("darcy_flux", passwd_)->ScatterMasterToGhosted("face");

  IdentifyUpwindCells();

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
      *vo_->os() << "cell " << cmin_dT << " has smallest dT, (" << p[0] << ", " << p[1];
      if (p.dim() == 3) *vo_->os() << ", " << p[2];
      *vo_->os() << ")" << std::endl;
    }
  }
  return dT;
}


/* ******************************************************************* 
* Estimate returns last time step unless it is zero.     
******************************************************************* */
double Transport_PK::get_dt()
{
  //if (dT == 0.0) 
  CalculateTransportDt();
  return dT;
}


/* ******************************************************************* 
* MPC will call this function to advance the transport state.
* Efficient subcycling requires to calculate an intermediate state of
* saturation only once, which leads to a leap-frog-type algorithm.
******************************************************************* */
int Transport_PK::Advance(double dT_MPC, double& dT_actual)
{ 
  // We use original tcc and make a copy of it later if needed.
  tcc = S_->GetFieldData("total_component_concentration", passwd_);
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell");

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
    time = T_physics + dT_cycle / 2;
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
    }

    if (! final_cycle) {  // rotate concentrations (we need new memory for tcc)
      tcc = Teuchos::RCP<CompositeVector>(new CompositeVector(*tcc_tmp));
    }

    ncycles++;
  }

  dT = dT_original;  // restore the original dT (just in case)

  // We define tracer as the species #0 as calculate some statistics.
  int num_components = tcc_prev.NumVectors();
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", false);

  bool flag_dispersion = (dispersion_models_ != TRANSPORT_DISPERSIVITY_MODEL_NULL);
  bool flag_diffusion(false);
  for (int i = 0; i < 2; i++) {
    if (diffusion_phase_[i] != Teuchos::null) {
      if (diffusion_phase_[i]->values().size() != 0) flag_diffusion = true;
    }
  }

  if (flag_dispersion || flag_diffusion) {
    Teuchos::ParameterList op_list;
    op_list.set<std::string>("discretization secondary", "mfd: two-point flux approximation");

    if (mesh_->mesh_type() == AmanziMesh::RECTANGULAR) {
      op_list.set<std::string>("discretization primary", "mfd: monotone for hex");
      Teuchos::Array<std::string> stensil(2);
      stensil[0] = "face";
      stensil[1] = "cell";
      op_list.set<Teuchos::Array<std::string> >("schema", stensil);
      stensil.remove(1);
      op_list.set<Teuchos::Array<std::string> >("preconditioner schema", stensil);
    } else {
      op_list.set<std::string>("discretization primary", "mfd: two-point flux approximation");
      Teuchos::Array<std::string> stensil(1);
      stensil[0] = "cell";
      op_list.set<Teuchos::Array<std::string> >("schema", stensil);
    }

    // default boundary conditions (none inside domain and Neumann on its boundary)
    std::vector<int> bc_model(nfaces_wghost, Operators::OPERATOR_BC_NONE);
    std::vector<double> bc_value(nfaces_wghost, 0.0);
    std::vector<double> bc_mixed;
    PopulateBoundaryData(bc_model, bc_value, -1);

    Teuchos::RCP<Operators::BCs> bc_dummy = 
        Teuchos::rcp(new Operators::BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model, bc_value, bc_mixed));
    AmanziGeometry::Point g;

    Operators::OperatorDiffusionFactory opfactory;
    Teuchos::RCP<Operators::OperatorDiffusion> op1 = opfactory.Create(mesh_, bc_dummy, op_list, g, 0);

    const CompositeVectorSpace& cvs = op1->DomainMap();
    CompositeVector sol(cvs), factor(cvs), factor0(cvs), source(cvs), zero(cvs);
    zero.PutScalar(0.0);
  
    // instantiale solver
    AmanziSolvers::LinearOperatorFactory<Operators::OperatorDiffusion, CompositeVector, CompositeVectorSpace> sfactory;
    Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::OperatorDiffusion, CompositeVector, CompositeVectorSpace> >
       solver = sfactory.Create(dispersion_solver, solvers_list, op1);

    solver->add_criteria(AmanziSolvers::LIN_SOLVER_MAKE_ONE_ITERATION);  // Make at least one iteration

    // populate the dispersion operator (if any)
    if (flag_dispersion) {
      CalculateDispersionTensor_(*darcy_flux, *phi, *ws);
    }

    int phase, num_itrs(0);
    bool flag_op1(true);
    double md_change, md_old(0.0), md_new, residual(0.0);

    // disperse and diffuse aqueous componetns
    for (int i = 0; i < num_aqueous; i++) {
      FindDiffusionValue(component_names_[i], &md_new, &phase);
      md_change = md_new - md_old;
      md_old = md_new;

      if (md_change != 0.0) {
        CalculateDiffusionTensor_(md_change, phase, *phi, *ws);
        flag_op1 = true;
      }

      // set initial guess
      Epetra_MultiVector& sol_cell = *sol.ViewComponent("cell");
      for (int c = 0; c < ncells_owned; c++) {
        sol_cell[0][c] = tcc_next[i][c];
      }
      if (sol.HasComponent("face")) {
        sol.ViewComponent("face")->PutScalar(0.0);
      }

      if (flag_op1) {
        op1->Init();
        op1->Setup(D, Teuchos::null, Teuchos::null, 1.0, 1.0);
        op1->UpdateMatrices(Teuchos::null, Teuchos::null);
        op1->ApplyBCs();

        // add accumulation term
        Epetra_MultiVector& fac = *factor.ViewComponent("cell");
        for (int c = 0; c < ncells_owned; c++) {
          fac[0][c] = (*phi)[0][c] * (*ws)[0][c];
        }
        op1->AddAccumulationTerm(sol, factor, dT_MPC, "cell");
 
        int schema_prec_dofs = op1->schema_prec_dofs();
        op1->SymbolicAssembleMatrix(schema_prec_dofs);
        op1->AssembleMatrix(schema_prec_dofs);
        op1->InitPreconditioner(dispersion_preconditioner, preconditioners_list);
      } else {
        Epetra_MultiVector& rhs_cell = *op1->rhs()->ViewComponent("cell");
        for (int c = 0; c < ncells_owned; c++) {
          double tmp = mesh_->cell_volume(c) * (*ws)[0][c] * (*phi)[0][c] / dT_MPC;
          rhs_cell[0][c] = tcc_next[i][c] * tmp;
        }
      }
  
      CompositeVector& rhs = *op1->rhs();
      int ierr = solver->ApplyInverse(rhs, sol);

      if (ierr != 0) {
        Errors::Message msg;
        msg << "\nLinear solver returned an unrecoverable error code.\n";
        Exceptions::amanzi_throw(msg);
      }

      residual += solver->residual();
      num_itrs += solver->num_itrs();

      for (int c = 0; c < ncells_owned; c++) {
        tcc_next[i][c] = sol_cell[0][c];
      }
    }

    // Disperse and diffuse gaseous components. We ignore dispersion 
    // tensor (D is reset). Inactive cells (s[c] = 1 and D[c] = 0) 
    // are treated with a hack of the accumulation term.
    D.clear();
    md_old = 0.0;
    for (int i = num_aqueous; i < num_components; i++) {
      FindDiffusionValue(component_names_[i], &md_new, &phase);
      md_change = md_new - md_old;
      md_old = md_new;

      if (md_change != 0.0) {
        CalculateDiffusionTensor_(md_change, phase, *phi, *ws);
      }

      // set initial guess
      Epetra_MultiVector& sol_cell = *sol.ViewComponent("cell");
      for (int c = 0; c < ncells_owned; c++) {
        sol_cell[0][c] = tcc_next[i][c];
      }
      if (sol.HasComponent("face")) {
        sol.ViewComponent("face")->PutScalar(0.0);
      }

      op1->Init();
      op1->Setup(D, Teuchos::null, Teuchos::null, 1.0, 1.0);
      op1->UpdateMatrices(Teuchos::null, Teuchos::null);

      // add boundary conditions and sources for gaseous components
      PopulateBoundaryData(bc_model, bc_value, i);

      Epetra_MultiVector& rhs_cell = *op1->rhs()->ViewComponent("cell");
      double time = T_physics + dT_MPC;
      ComputeAddSourceTerms(time, 1.0, srcs, rhs_cell, i, i);
      op1->ApplyBCs();

      // add accumulation term
      Epetra_MultiVector& fac1 = *factor.ViewComponent("cell");
      Epetra_MultiVector& fac0 = *factor0.ViewComponent("cell");

      for (int c = 0; c < ncells_owned; c++) {
        fac1[0][c] = (*phi)[0][c] * (1.0 - (*ws)[0][c]);
        fac0[0][c] = (*phi)[0][c] * (1.0 - (*ws_prev)[0][c]);
        if ((*ws)[0][c] == 1.0) fac1[0][c] = 1.0;  // hack so far
      }
      op1->AddAccumulationTerm(sol, factor0, factor, dT_MPC, "cell");
 
      int schema_prec_dofs = op1->schema_prec_dofs();
      op1->SymbolicAssembleMatrix(schema_prec_dofs);
      op1->AssembleMatrix(schema_prec_dofs);
      op1->InitPreconditioner(dispersion_preconditioner, preconditioners_list);
  
      CompositeVector& rhs = *op1->rhs();
      int ierr = solver->ApplyInverse(rhs, sol);

      if (ierr != 0) {
        Errors::Message msg;
        msg << "\nLinear solver returned an unrecoverable error code.\n";
        Exceptions::amanzi_throw(msg);
      }

      residual += solver->residual();
      num_itrs += solver->num_itrs();

      for (int c = 0; c < ncells_owned; c++) {
        tcc_next[i][c] = sol_cell[0][c];
      }
    }

    if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "dispersion solver (" << solver->name() 
                 << ") ||r||=" << residual / num_components
                 << " itrs=" << num_itrs / num_components << std::endl;
    }
  }

  // statistics output
  nsubcycles = ncycles;
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << ncycles << " sub-cycles, dT_stable=" << dT_original 
               << " [sec]  dT_MPC=" << dT_MPC << " [sec]" << std::endl;

    VV_PrintSoluteExtrema(tcc_next, dT_MPC);
  }

  return 0;
}


/* ******************************************************************* 
* Copy theadvected tcc to the provided state.
******************************************************************* */
void Transport_PK::CommitState(double dummy_dT, const Teuchos::Ptr<State>& S) 
{
  Teuchos::RCP<CompositeVector> cv;
  cv = S->GetFieldData("total_component_concentration", passwd_);
  *cv = *tcc_tmp;
}


/* ******************************************************************* 
 * A simple first-order transport method 
 ****************************************************************** */
void Transport_PK::AdvanceDonorUpwind(double dT_cycle)
{
  dT = dT_cycle;  // overwrite the maximum stable transport step
  mass_solutes_source_.assign(num_aqueous + num_gaseous, 0.0);

  // populating next state of concentrations
  tcc->ScatterMasterToGhosted("cell");
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell", true);
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", true);

  // prepare conservative state in master and slave cells
  double vol_phi_ws, tcc_flux;

  // We advect only aqueous components.
  int num_advect = num_aqueous;

  for (int c = 0; c < ncells_owned; c++) {
    vol_phi_ws = mesh_->cell_volume(c) * (*phi)[0][c] * (*ws_start)[0][c];

    for (int i = 0; i < num_advect; i++)
      tcc_next[i][c] = tcc_prev[i][c] * vol_phi_ws;
  }

  // advance all components at once
  for (int f = 0; f < nfaces_wghost; f++) {  // loop over master and slave faces
    int c1 = (*upwind_cell_)[f];
    int c2 = (*downwind_cell_)[f];

    double u = fabs((*darcy_flux)[0][f]);

    if (c1 >=0 && c1 < ncells_owned && c2 >= 0 && c2 < ncells_owned) {
      for (int i = 0; i < num_advect; i++) {
        tcc_flux = dT * u * tcc_prev[i][c1];
        tcc_next[i][c1] -= tcc_flux;
        tcc_next[i][c2] += tcc_flux;
      }

    } else if (c1 >=0 && c1 < ncells_owned && (c2 >= ncells_owned || c2 < 0)) {
      for (int i = 0; i < num_advect; i++) {
        tcc_flux = dT * u * tcc_prev[i][c1];
        tcc_next[i][c1] -= tcc_flux;
      }

    } else if (c1 >= ncells_owned && c2 >= 0 && c2 < ncells_owned) {
      for (int i = 0; i < num_advect; i++) {
        tcc_flux = dT * u * tcc_prev[i][c1];
        tcc_next[i][c2] += tcc_flux;
      }
    }
  }

  // loop over exterior boundary sets
  for (int m = 0; m < bcs.size(); m++) {
    std::vector<int>& tcc_index = bcs[m]->tcc_index();
    std::vector<int>& faces = bcs[m]->faces();
    std::vector<std::vector<double> >& values = bcs[m]->values();

    int ncomp = tcc_index.size();
    int nbfaces = faces.size();
    for (int n = 0; n < nbfaces; n++) {
      int f = faces[n];
      int c2 = (*downwind_cell_)[f];

      if (c2 >= 0) {
        double u = fabs((*darcy_flux)[0][f]);
        for (int i = 0; i < ncomp; i++) {
          int k = tcc_index[i];
          if (k < num_advect) {
            tcc_flux = dT * u * values[n][i];
            tcc_next[k][c2] += tcc_flux;
          }
        }
      }
    }
  }

  // process external sources
  if (srcs.size() != 0) {
    double time = T_physics;
    ComputeAddSourceTerms(time, dT, srcs, tcc_next, 0, num_advect - 1);
  }

  // recover concentration from new conservative state
  for (int c = 0; c < ncells_owned; c++) {
    vol_phi_ws = mesh_->cell_volume(c) * (*phi)[0][c] * (*ws_end)[0][c];
    for (int i = 0; i < num_advect; i++) tcc_next[i][c] /= vol_phi_ws;
  }

  // update mass balance
  for (int i = 0; i < mass_solutes_exact_.size(); i++) {
    mass_solutes_exact_[i] += mass_solutes_source_[i] * dT;
  }

  if (internal_tests) {
    VV_CheckGEDproperty(*tcc_tmp->ViewComponent("cell"));
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
  mass_solutes_source_.assign(num_aqueous + num_gaseous, 0.0);

  // work memory
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  Epetra_Vector f_component(cmap_wghost);

  // distribute vector of concentrations
  S_->GetFieldData("total_component_concentration")->ScatterMasterToGhosted("cell");
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell", true);
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", true);

  // We advect only aqueous components.
  int num_advect = num_aqueous;

  for (int i = 0; i < num_advect; i++) {
    current_component_ = i;  // needed by BJ 

    double T = T_physics;
    Epetra_Vector*& component = tcc_prev(i);
    Functional(T, *component, f_component);

    double ws_ratio;
    for (int c = 0; c < ncells_owned; c++) {
      ws_ratio = (*ws_start)[0][c] / (*ws_end)[0][c];
      tcc_next[i][c] = (tcc_prev[i][c] + dT * f_component[c]) * ws_ratio;
    }
  }

  // update mass balance
  for (int i = 0; i < num_aqueous + num_gaseous; i++) {
    mass_solutes_exact_[i] += mass_solutes_source_[i] * dT;
  }

  if (internal_tests) {
    VV_CheckGEDproperty(*tcc_tmp->ViewComponent("cell"));
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
  mass_solutes_source_.assign(num_aqueous + num_gaseous, 0.0);

  // work memory
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  Epetra_Vector f_component(cmap_wghost);

  // distribute old vector of concentrations
  S_->GetFieldData("total_component_concentration")->ScatterMasterToGhosted("cell");
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell", true);
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", true);

  Epetra_Vector ws_ratio(Copy, *ws_start, 0);
  for (int c = 0; c < ncells_owned; c++) ws_ratio[c] /= (*ws_end)[0][c];

  // We advect only aqueous components.
  int num_advect = num_aqueous;

  // predictor step
  for (int i = 0; i < num_advect; i++) {
    current_component_ = i;  // needed by BJ 

    double T = T_physics;
    Epetra_Vector*& component = tcc_prev(i);
    Functional(T, *component, f_component);

    for (int c = 0; c < ncells_owned; c++) {
      tcc_next[i][c] = (tcc_prev[i][c] + dT * f_component[c]) * ws_ratio[c];
    }
  }

  tcc_tmp->ScatterMasterToGhosted("cell");

  // corrector step
  for (int i = 0; i < num_advect; i++) {
    current_component_ = i;  // needed by BJ 

    double T = T_physics;
    Epetra_Vector*& component = tcc_next(i);
    Functional(T, *component, f_component);

    for (int c = 0; c < ncells_owned; c++) {
      double value = (tcc_prev[i][c] + dT * f_component[c]) * ws_ratio[c];
      tcc_next[i][c] = (tcc_next[i][c] + value) / 2;
    }
  }

  // update mass balance
  for (int i = 0; i < num_aqueous + num_gaseous; i++) {
    mass_solutes_exact_[i] += mass_solutes_source_[i] * dT / 2;
  }

  if (internal_tests) {
    VV_CheckGEDproperty(*tcc_tmp->ViewComponent("cell"));
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

  // We advect only aqueous components.
  // int ncomponents = tcc_next.NumVectors();
  int ncomponents = num_aqueous;

  for (int i = 0; i < ncomponents; i++) {
    current_component_ = i;  // it is needed in BJ called inside RK:fun

    Epetra_Vector*& tcc_component = (*tcc)(i);
    TS_nextBIG->CopyMasterCell2GhostCell(*tcc_component, *component_);

    double T = 0.0;  // requires fixes (lipnikov@lanl.gov)
    TVD_RK.step(T, dT, *component_, *component_next_);

    for (int c = 0; c < ncells_owned; c++) (*tcc_next)[i][c] = (*component_next_)[c];
  }
}
*/


/* ******************************************************************
* Computes source and sink terms and adds them to vector tcc.
* Return mass rate for the tracer.
* The routine treats two cases of tcc with one and all components.
****************************************************************** */
void Transport_PK::ComputeAddSourceTerms(double Tp, double dTp, 
                                         std::vector<TransportDomainFunction*>& srcs, 
                                         Epetra_MultiVector& tcc, int n0, int n1)
{
  int num_vectors = tcc.NumVectors();
  int nsrcs = srcs.size();

  for (int m = 0; m < nsrcs; m++) {
    int i = srcs[m]->tcc_index();
    if (i < n0 || i > n1) continue;

    int imap = i;
    if (num_vectors == 1) imap = 0;

    double T0 = Tp - dTp;
    int distribution = srcs[m]->CollectActionsList();
    if (distribution & CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      srcs[m]->ComputeDistributeMultiValue(T0, Tp, Kxy->Values()); 
    } else {
      srcs[m]->ComputeDistributeMultiValue(T0, Tp, NULL);
    }

    TransportDomainFunction::Iterator it;

    for (it = srcs[m]->begin(); it != srcs[m]->end(); ++it) {
      int c = it->first;
      double value = mesh_->cell_volume(c) * it->second;

      tcc[imap][c] += dTp * value;
      mass_solutes_source_[i] += value;
    }
  }
}


/* *******************************************************************
* Populates operators' boundary data for given component.
* Returns true if at least one face was populated.
******************************************************************* */
bool Transport_PK::PopulateBoundaryData(
    std::vector<int>& bc_model, std::vector<double>& bc_value, int component)
{
  bool flag = false;

  for (int i = 0; i < bc_model.size(); i++) {
    bc_model[i] = Operators::OPERATOR_BC_NONE;
    bc_value[i] = 0.0;
  }

  AmanziMesh::Entity_ID_List cells;
  for (int f = 0; f < nfaces_wghost; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    if (cells.size() == 1) bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
  }

  for (int m = 0; m < bcs.size(); m++) {
    std::vector<int>& tcc_index = bcs[m]->tcc_index();
    std::vector<int>& faces = bcs[m]->faces();
    std::vector<std::vector<double> >& values = bcs[m]->values();

    int ncomp = tcc_index.size();
    int nbfaces = faces.size();
    for (int n = 0; n < nbfaces; n++) {
      int f = faces[n];

      for (int i = 0; i < ncomp; i++) {
        int k = tcc_index[i];
        if (k == component) {
          bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_value[f] = values[n][i];
          flag = true;
        }
      }
    }
  }

  return flag;
}


/* *******************************************************************
* Identify flux direction based on orientation of the face normal 
* and sign of the  Darcy velocity.                               
******************************************************************* */
void Transport_PK::IdentifyUpwindCells()
{
  for (int f = 0; f < nfaces_wghost; f++) {
    (*upwind_cell_)[f] = -1;  // negative value indicates boundary
    (*downwind_cell_)[f] = -1;

    }
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

    for (int i = 0; i < faces.size(); i++) {
      int f = faces[i];
      double tmp = (*darcy_flux)[0][f] * dirs[i];
      if (tmp > 0.0) {
        (*upwind_cell_)[f] = c;
      } else if (tmp < 0.0) {
        (*downwind_cell_)[f] = c;
      } else if (dirs[i] > 0) {
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

}  // namespace Transport
}  // namespace Amanzi


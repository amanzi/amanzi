/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "Teuchos_RCP.hpp"

#include "errors.hh"
#include "tabular-function.hh"

#include "Mesh.hh"
#include "Flow_PK.hh"


namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
void Flow_PK::ProcessParameterList(Teuchos::ParameterList& plist)
{
  double rho = *(S_->GetScalarData("fluid_density"));

  // Process main one-line options (not sublists)
  atm_pressure_ = plist.get<double>("atmospheric pressure", FLOW_PRESSURE_ATMOSPHERIC);
 
  std::string mfd3d_method_name = plist.get<std::string>("discretization method", "optimized mfd");
  ProcessStringMFD3D(mfd3d_method_name, &mfd3d_method_); 

  // Create the BC objects.
  Teuchos::RCP<Teuchos::ParameterList>
      bc_list = Teuchos::rcp(new Teuchos::ParameterList(plist.sublist("boundary conditions", true)));
  FlowBCFactory bc_factory(mesh_, bc_list);

  bc_pressure = bc_factory.CreatePressure(bc_submodel);
  bc_head = bc_factory.CreateStaticHead(atm_pressure_, rho, gravity_, bc_submodel);
  bc_flux = bc_factory.CreateMassFlux(bc_submodel);
  bc_seepage = bc_factory.CreateSeepageFace(atm_pressure_, bc_submodel);

  ValidateBCs();
  ProcessBCs();

  // Create the source object if any
  if (plist.isSublist("source terms")) {
    std::string distribution_method_name = plist.get<std::string>("source and sink distribution method", "none");
    ProcessStringSourceDistribution(distribution_method_name, &src_sink_distribution); 

    Teuchos::RCP<Teuchos::ParameterList> src_list = Teuchos::rcpFromRef(plist.sublist("source terms", true));
    FlowSourceFactory src_factory(mesh_, src_list);
    src_sink = src_factory.createSource();
    src_sink_distribution = src_sink->CollectActionsList();
  }

  // Time integrator for period I, temporary called initial guess initialization
  if (plist.isSublist("initial guess pseudo time integrator")) {
    Teuchos::ParameterList& igs_list = plist.sublist("initial guess pseudo time integrator");

    std::string ti_method_name = igs_list.get<std::string>("time integration method", "none");
    ProcessStringTimeIntegration(ti_method_name, &ti_specs_igs_.ti_method);
    ProcessSublistTimeIntegration(igs_list, ti_method_name, ti_specs_igs_);
    ti_specs_igs_.ti_method_name = "initial guess pseudo time integrator";

    ti_specs_igs_.preconditioner_name = FindStringPreconditioner(igs_list);
    ProcessStringPreconditioner(ti_specs_igs_.preconditioner_name, &ti_specs_igs_.preconditioner_method);

    std::string linear_solver_name = FindStringLinearSolver(igs_list);
    ProcessStringLinearSolver(linear_solver_name, &ti_specs_igs_.ls_specs);

    ProcessStringPreconditioner(ti_specs_igs_.preconditioner_name, &ti_specs_igs_.preconditioner_method);
    ProcessStringErrorOptions(igs_list, &ti_specs_igs_.error_control_options);
  }

  // Time integrator for period II, temporary called steady-state time integrator
  if (plist.isSublist("steady state time integrator")) {
    Teuchos::ParameterList& sss_list = plist.sublist("steady state time integrator");

    std::string ti_method_name = sss_list.get<std::string>("time integration method", "none");
    ProcessStringTimeIntegration(ti_method_name, &ti_specs_sss_.ti_method);
    ProcessSublistTimeIntegration(sss_list, ti_method_name, ti_specs_sss_);
    ti_specs_sss_.ti_method_name = "steady state time integrator";

    ti_specs_sss_.preconditioner_name = FindStringPreconditioner(sss_list);
    ProcessStringPreconditioner(ti_specs_sss_.preconditioner_name, &ti_specs_sss_.preconditioner_method);

    ti_specs_sss_.ls_specs.solver_name = FindStringLinearSolver(sss_list);
    ProcessStringLinearSolver(ti_specs_sss_.ls_specs.solver_name, &ti_specs_sss_.ls_specs);

    ProcessStringPreconditioner(ti_specs_sss_.preconditioner_name, &ti_specs_sss_.preconditioner_method);
    ProcessStringErrorOptions(sss_list, &ti_specs_sss_.error_control_options);

  } else if (vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << vo_->color("yellow") << "steady-state calculation was not requested." 
               << vo_->reset() << std::endl;
  }

  // Time integrator for period III, called transient time integrator
  if (plist.isSublist("transient time integrator")) {
    Teuchos::ParameterList& trs_list = plist.sublist("transient time integrator");

    std::string ti_method_name = trs_list.get<std::string>("time integration method", "none");
    ProcessStringTimeIntegration(ti_method_name, &ti_specs_trs_.ti_method);
    ProcessSublistTimeIntegration(trs_list, ti_method_name, ti_specs_trs_);
    ti_specs_trs_.ti_method_name = "transient time integrator";

    ti_specs_trs_.preconditioner_name = FindStringPreconditioner(trs_list);
    ProcessStringPreconditioner(ti_specs_trs_.preconditioner_name, &ti_specs_trs_.preconditioner_method);

    std::string linear_solver_name = FindStringLinearSolver(trs_list);
    ProcessStringLinearSolver(linear_solver_name, &ti_specs_trs_.ls_specs);

    ProcessStringPreconditioner(ti_specs_trs_.preconditioner_name, &ti_specs_trs_.preconditioner_method);
    ProcessStringErrorOptions(trs_list, &ti_specs_trs_.error_control_options);

  } else if (vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << vo_->color("yellow") << "no sublist '\"transient time integrator'\"" 
               << vo_->reset() << std::endl;
  }
}


/* ******************************************************************
* Process time integration sublist.
**************************************************************** */
void Flow_PK::ProcessSublistTimeIntegration(
    Teuchos::ParameterList& list, const std::string name, TI_Specs& ti_specs)
{
  Errors::Message msg;

  if (list.isSublist(name)) {
    Teuchos::ParameterList bdf1_list = list.sublist(name);
    // obsolete way to define parameters
    if (list.isSublist("obsolete parameters")) {
      Teuchos::ParameterList olist = list.sublist("obsolete parameters");
      ti_specs.initialize_with_darcy = olist.get<bool>("initialize with darcy", false);
      ti_specs.clip_saturation = olist.get<double>("clipping saturation value", -1.0);
      ti_specs.clip_pressure = olist.get<double>("clipping pressure value", -1e+10);

      ti_specs.pressure_lambda_constraints = olist.get<bool>("enforce pressure-lambda constraints", true);

      ti_specs.dT_method = 0;
      std::string dT_name = olist.get<std::string>("time stepping strategy", "standard");
      if (dT_name == "adaptive") ti_specs.dT_method = FLOW_DT_ADAPTIVE;

      ti_specs.dTfactor = olist.get<double>("time step increase factor", 1.0);

      ti_specs.T0 = olist.get<double>("start time", 0.0);
      ti_specs.T1 = olist.get<double>("end time", 100 * AmanziFlow::FLOW_YEAR);
      ti_specs.dT0 = olist.get<double>("initial time step", AmanziFlow::FLOW_INITIAL_DT);
      ti_specs.dTmax = olist.get<double>("max time step", AmanziFlow::FLOW_MAXIMUM_DT);

      ti_specs.residual_tol = olist.get<double>("convergence tolerance", FLOW_TI_NONLINEAR_RESIDUAL_TOLERANCE);
      ti_specs.max_itrs = olist.get<int>("maximum number of iterations", FLOW_TI_MAX_ITERATIONS);
    }

    // new way to define parameters overrides the above values.
    if (bdf1_list.isParameter("timestep controller type")) {
      std::string dT_method_name = bdf1_list.get<std::string>("timestep controller type");

      ti_specs.dT_method = 0;
      if (dT_method_name == "standard") {
        Teuchos::ParameterList dtlist = bdf1_list.sublist("timestep controller standard parameters");
        ti_specs.dTfactor = dtlist.get<double>("time step increase factor");
      } else if (dT_method_name == "fixed") {
        Teuchos::ParameterList dtlist = bdf1_list.sublist("timestep controller fixed parameters");
        ti_specs.dTfactor = dtlist.get<double>("time step increase factor");
      } else if (dT_method_name == "adaptive") {
        ti_specs.dT_method = FLOW_DT_ADAPTIVE;
      }
    }

    if (list.isSublist("initialization")) {
      Teuchos::ParameterList& ini_list = list.sublist("initialization");
      std::string name = ini_list.get<std::string>("method", "none");
      ti_specs.initialize_with_darcy = (name == "saturated solver");
      ti_specs.clip_saturation = ini_list.get<double>("clipping saturation value", -1.0);
      ti_specs.clip_pressure = ini_list.get<double>("clipping pressure value", -1e+10);

      std::string linear_solver_name = FindStringLinearSolver(ini_list);
      ProcessStringLinearSolver(linear_solver_name, &ti_specs.ls_specs_ini);
    }

    if (list.isSublist("pressure-lambda constraints")) {
      Teuchos::ParameterList& pl_list = list.sublist("pressure-lambda constraints");
      ti_specs.pressure_lambda_constraints = true;

      std::string linear_solver_name = FindStringLinearSolver(pl_list);
      ProcessStringLinearSolver(linear_solver_name, &ti_specs.ls_specs_constraints);
    }

  } else if (name != "none") {
    msg << "\nFlow PK: specified time integration sublist does not exist.";
    Exceptions::amanzi_throw(msg);
  }
}


/* ****************************************************************
* Process string for the time integration method.
**************************************************************** */
void Flow_PK::ProcessStringTimeIntegration(const std::string name, int* method)
{
  Errors::Message msg;
  if (name == "Picard") {
    *method = AmanziFlow::FLOW_TIME_INTEGRATION_PICARD;
  } else if (name == "backward Euler") {
    *method = AmanziFlow::FLOW_TIME_INTEGRATION_BACKWARD_EULER;
  } else if (name == "BDF1") {
    *method = AmanziFlow::FLOW_TIME_INTEGRATION_BDF1;
  } else if (name == "BDF2") {
    *method = AmanziFlow::FLOW_TIME_INTEGRATION_BDF2;
  } else {
    msg << "Flow PK: time integration method \"" << name.c_str() << "\" is not known.";
    Exceptions::amanzi_throw(msg);
  }
}


/* ****************************************************************
* Process string for the linear solver.
**************************************************************** */
void Flow_PK::ProcessStringSourceDistribution(const std::string name, int* method)
{
  if (name != "none") {
    Errors::Message msg;
    msg << "\nFlow_PK: \"source and sink distribution method\" is obsolete.\n"
        << "         see desription of sublist \"source terms\" in the native spec.\n";
    Exceptions::amanzi_throw(msg);
  }
}


/* ****************************************************************
* Process string for the discretization method.
**************************************************************** */
void Flow_PK::ProcessStringMFD3D(const std::string name, int* method)
{
  if (name == "monotone mfd hex") {  // two monotone methods
    *method = FLOW_MFD3D_HEXAHEDRA_MONOTONE;
  } else if (name == "monotone mfd") {
    *method = FLOW_MFD3D_POLYHEDRA_MONOTONE;
  } else if (name == "support operator") {
    *method = FLOW_MFD3D_SUPPORT_OPERATOR;
  } else if (name == "two point flux approximation") {
    *method = FLOW_MFD3D_TPFA;
  } else if (name == "finite volume") {
    *method = FLOW_FV_TPFA;
  } else if (name == "optimized mfd") {
    *method = FLOW_MFD3D_OPTIMIZED;
  } else if (name == "optimized mfd scaled") {
    *method = FLOW_MFD3D_OPTIMIZED_SCALED;
  } else if (name == "mfd") {  // first basic mfd
    *method = FLOW_MFD3D_POLYHEDRA;
  } else if (name == "mfd scaled") {  // second basic mfd
    *method = FLOW_MFD3D_POLYHEDRA_SCALED;
  } else {
    *method = FLOW_MFD3D_POLYHEDRA;
  }
}


/* ****************************************************************
* Process string for the linear solver.
**************************************************************** */
void Flow_PK::ProcessStringLinearSolver(const std::string& name, LinearSolver_Specs* ls_specs)
{
  Errors::Message msg;

  if (! linear_operator_list_.isSublist(name)) {
    msg << "Flow PK: linear solver does not exist for a time integrator.";
    Exceptions::amanzi_throw(msg);
  }

  Teuchos::ParameterList& tmp_list = linear_operator_list_.sublist(name);
  ls_specs->max_itrs = tmp_list.get<int>("maximum number of iterations", 100);
  ls_specs->convergence_tol = tmp_list.get<double>("error tolerance", 1e-14);

  ls_specs->preconditioner_name = FindStringPreconditioner(tmp_list);
  ProcessStringPreconditioner(ls_specs->preconditioner_name, &ls_specs->preconditioner_method);
}


/* ****************************************************************
* Find string for the preconditoner.
**************************************************************** */
std::string Flow_PK::FindStringPreconditioner(const Teuchos::ParameterList& list)
{   
  Errors::Message msg;
  std::string name; 

  if (list.isParameter("preconditioner")) {
    name = list.get<std::string>("preconditioner");
  } else {
    msg << "Flow PK: parameter <preconditioner> is missing either in TI or LS list.";
    Exceptions::amanzi_throw(msg);
  }

  if (! preconditioner_list_.isSublist(name)) {
    msg << "Flow PK: preconditioner \"" << name.c_str() << "\" does not exist.";
    Exceptions::amanzi_throw(msg);
  }
  return name;
}


/* ****************************************************************
* Process string for the preconitioner.
**************************************************************** */
void Flow_PK::ProcessStringPreconditioner(const std::string& name, int* preconditioner)
{
  Errors::Message msg;

  if (name == "Trilinos ML") {
    *preconditioner = FLOW_PRECONDITIONER_TRILINOS_ML;
  } else if (name == "Hypre AMG") {
    *preconditioner = FLOW_PRECONDITIONER_HYPRE_AMG;
  } else if (name == "Block ILU") {
    *preconditioner = FLOW_PRECONDITIONER_TRILINOS_BLOCK_ILU;
  } else {
    msg << "\nFlow PK: specified preconditioner does not exist.";
    Exceptions::amanzi_throw(msg);
  }
}


/* ****************************************************************
* Find string for the preconditoner.
**************************************************************** */
std::string Flow_PK::FindStringLinearSolver(const Teuchos::ParameterList& list)
{   
  Errors::Message msg;
  std::string name; 

  if (list.isParameter("linear solver")) {
    name = list.get<std::string>("linear solver");
  } else {
    msg << "Flow PK: time integrator does not define <linear solver>.";
    Exceptions::amanzi_throw(msg);
  }

  if (! linear_operator_list_.isSublist(name)) {
    msg << "Flow PK: linear solver \"" << name.c_str() << "\" does not exist.";
    Exceptions::amanzi_throw(msg);
  }
  return name;
}


/* ****************************************************************
* Process string for error control options
**************************************************************** */
void Flow_PK::ProcessStringErrorOptions(Teuchos::ParameterList& list, int* control)
{
  *control = 0;
  if (list.isParameter("error control options")){
    std::vector<std::string> options;
    options = list.get<Teuchos::Array<std::string> >("error control options").toVector();

    for (int i=0; i < options.size(); i++) {
      if (options[i] == "pressure") {
        *control += FLOW_TI_ERROR_CONTROL_PRESSURE;
      } else if (options[i] == "saturation") {
        *control += FLOW_TI_ERROR_CONTROL_SATURATION;
      } else if (options[i] == "residual") {
        *control += FLOW_TI_ERROR_CONTROL_RESIDUAL;
      } else {
        Errors::Message msg;
        msg << "Flow PK: unknown error control option has been specified.";
        Exceptions::amanzi_throw(msg);
      }
    }
  }
}


/* ****************************************************************
* Find string for the preconditoner.
**************************************************************** */
void Flow_PK::OutputTimeHistory(
    const Teuchos::ParameterList& plist, std::vector<dt_tuple>& dT_history)
{
  if (plist.isParameter("plot time history") && 
      vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "saving time history in file flow_dt_history.txt..." << std::endl;

    char file_name[30];
    sprintf(file_name, "flow_dt_history_%d.txt", ti_phase_counter++);

    std::ofstream ofile;
    ofile.open(file_name);

    for (double n = 0; n < dT_history.size(); n++) {
      ofile << std::setprecision(10) << dT_history[n].first / FLOW_YEAR << " " << dT_history[n].second << std::endl;
    }
    ofile.close();
  }
}


}  // namespace AmanziFlow
}  // namespace Amanzi


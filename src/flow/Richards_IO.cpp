/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Neil Carlson (nnc@lanl.gov), 
         Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <set>
#include <vector>
#include <string>

#include "Teuchos_ParameterList.hpp"

#include "WRM_vanGenuchten.hpp"
#include "WRM_BrooksCorey.hpp"
#include "WRM_fake.hpp"
#include "Richards_PK.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
void Richards_PK::ProcessParameterList()
{
  Errors::Message msg;

  // create verbosity list if it does not exist
  if (! rp_list_.isSublist("VerboseObject")) {
    Teuchos::ParameterList verbosity_list;
    verbosity_list.set<std::string>("Verbosity Level", "none");
    rp_list_.set("VerboseObject", verbosity_list);
  }

  // extract verbosity level
  Teuchos::ParameterList verbosity_list = rp_list_.get<Teuchos::ParameterList>("VerboseObject");
  std::string verbosity_name = verbosity_list.get<std::string>("Verbosity Level");
  ProcessStringVerbosity(verbosity_name, &verbosity);

  // Process main one-line options (not sublists)
  std::string krel_method_name = rp_list_.get<string>("relative permeability");
  ProcessStringRelativePermeability(krel_method_name, &Krel_method);
 
  atm_pressure = rp_list_.get<double>("atmospheric pressure", 101325.0);
 
  string mfd3d_method_name = rp_list_.get<string>("discretization method", "optimized mfd");
  ProcessStringMFD3D(mfd3d_method_name, &mfd3d_method_); 

  // Create the BC objects.
  Teuchos::RCP<Teuchos::ParameterList>
      bc_list = Teuchos::rcp(new Teuchos::ParameterList(rp_list_.sublist("boundary conditions", true)));
  FlowBCFactory bc_factory(mesh_, bc_list);

  bc_pressure = bc_factory.CreatePressure(bc_submodel);
  bc_head = bc_factory.CreateStaticHead(atm_pressure, rho, gravity_, bc_submodel);
  bc_flux = bc_factory.CreateMassFlux(bc_submodel);
  bc_seepage = bc_factory.CreateSeepageFace(bc_submodel);

  ValidateBoundaryConditions(bc_pressure, bc_head, bc_flux);
  ProcessStaticBCsubmodels(bc_submodel, rainfall_factor);

  /*
  double time = T_physics;
  bc_pressure->Compute(time);
  bc_head->Compute(time);
  bc_flux->Compute(time);
  bc_seepage->Compute(time);
  */

  // Create the source object if any
  if (rp_list_.isSublist("source terms")) {
    string distribution_method_name = rp_list_.get<string>("source and sink distribution method", "none");
    ProcessStringSourceDistribution(distribution_method_name, &src_sink_distribution); 

    Teuchos::RCP<Teuchos::ParameterList> src_list = Teuchos::rcpFromRef(rp_list_.sublist("source terms", true));
    FlowSourceFactory src_factory(mesh_, src_list);
    src_sink = src_factory.createSource();
  } else {
    src_sink = NULL;
  }

  // Create water retention models
  if (! rp_list_.isSublist("Water retention models")) {
    msg << "There is no Water retention models list";
    Exceptions::amanzi_throw(msg);
  }
  Teuchos::ParameterList& vG_list = rp_list_.sublist("Water retention models");

  int nblocks = 0;  // Find out how many WRM entries there are.
  for (Teuchos::ParameterList::ConstIterator i = vG_list.begin(); i != vG_list.end(); i++) {
    if (vG_list.isSublist(vG_list.name(i))) {
      nblocks++;
    } else {
      msg << "Richards Problem: water retention model contains an entry that is not a sublist.";
      Exceptions::amanzi_throw(msg);
    }
  }

  WRM.resize(nblocks);

  int iblock = 0;
  for (Teuchos::ParameterList::ConstIterator i = vG_list.begin(); i != vG_list.end(); i++) {
    if (vG_list.isSublist(vG_list.name(i))) {
      Teuchos::ParameterList& wrm_list = vG_list.sublist(vG_list.name(i));

      if (wrm_list.get<string>("Water retention model") == "van Genuchten") {
        std::string region = wrm_list.get<std::string>("Region");  // associated mesh block

        double m = wrm_list.get<double>("van Genuchten m");
        double alpha = wrm_list.get<double>("van Genuchten alpha");
        double l = wrm_list.get<double>("van Genuchten l", 0.5);
        double sr = wrm_list.get<double>("residual saturation");
        double pc0 = wrm_list.get<double>("regularization interval", 0.0);
        std::string krel_function = wrm_list.get<std::string>("relative permeability model", "Mualem");
        VerifyWRMparameters(m, alpha, sr, pc0);
        VerifyStringMualemBurdine(krel_function);

        WRM[iblock] = Teuchos::rcp(new WRM_vanGenuchten(region, m, l, alpha, sr, krel_function, pc0));

      } else if (wrm_list.get<string>("Water retention model") == "Brooks Corey") {
        std::string region = wrm_list.get<std::string>("Region");  // associated mesh block

        double lambda = wrm_list.get<double>("Brooks Corey lambda");
        double alpha = wrm_list.get<double>("Brooks Corey alpha");
        double l = wrm_list.get<double>("Brooks Corey l", 0.5);
        double sr = wrm_list.get<double>("residual saturation");
        double pc0 = wrm_list.get<double>("regularization interval", 0.0);
        std::string krel_function = wrm_list.get<std::string>("relative permeability model", "Mualem");
        VerifyWRMparameters(lambda, alpha, sr, pc0);
        VerifyStringMualemBurdine(krel_function);

        WRM[iblock] = Teuchos::rcp(new WRM_BrooksCorey(region, lambda, l, alpha, sr, krel_function, pc0));

      } else if (wrm_list.get<string>("Water retention model") == "fake") {
        std::string region = wrm_list.get<std::string>("Region");  // associated mesh block
        WRM[iblock] = Teuchos::rcp(new WRM_fake(region));

      } else {
        msg << "Richards Problem: unknown water retention model.";
        Exceptions::amanzi_throw(msg);
      }
      iblock++;
    }
  }

  // Time integrator for period I, temporary called initial guess initialization
  if (rp_list_.isSublist("initial guess pseudo time integrator")) {
    Teuchos::ParameterList& igs_list = rp_list_.sublist("initial guess pseudo time integrator");

    std::string ti_method_name = igs_list.get<string>("time integration method", "none");
    ProcessStringTimeIntegration(ti_method_name, &ti_specs_igs_.ti_method);
    ProcessSublistTimeIntegration(igs_list, ti_method_name, ti_specs_igs_);
    ti_specs_igs_.ti_method_name = "initial guess pseudo time integrator";

    std::string linear_solver_name = FindStringLinearSolver(igs_list, solver_list_);
    LinearSolver_Specs& ls_specs = ti_specs_igs_.ls_specs;
    ProcessStringLinearSolver(linear_solver_name, &ls_specs.max_itrs, &ls_specs.convergence_tol, 
                              &ti_specs_igs_.preconditioner_name);

    ProcessStringPreconditioner(ti_specs_igs_.preconditioner_name, &ti_specs_igs_.preconditioner_method);
    ProcessStringErrorOptions(igs_list, &ti_specs_igs_.error_control_options);
  }

  // Time integrator for period II, temporary called steady-state time integrator
  if (rp_list_.isSublist("steady state time integrator")) {
    Teuchos::ParameterList& sss_list = rp_list_.sublist("steady state time integrator");

    std::string ti_method_name = sss_list.get<string>("time integration method", "none");
    ProcessStringTimeIntegration(ti_method_name, &ti_specs_sss_.ti_method);
    ProcessSublistTimeIntegration(sss_list, ti_method_name, ti_specs_sss_);
    ti_specs_sss_.ti_method_name = "steady state time integrator";

    std::string linear_solver_name = FindStringLinearSolver(sss_list, solver_list_);
    LinearSolver_Specs& ls_specs = ti_specs_sss_.ls_specs;
    ProcessStringLinearSolver(linear_solver_name, &ls_specs.max_itrs, &ls_specs.convergence_tol,
                              &ti_specs_sss_.preconditioner_name);

    ProcessStringPreconditioner(ti_specs_sss_.preconditioner_name, &ti_specs_sss_.preconditioner_method);
    ProcessStringErrorOptions(sss_list, &ti_specs_sss_.error_control_options);

  } else if (verbosity >= FLOW_VERBOSITY_LOW) {
    printf("Flow PK: mandatory sublist for steady-state calculations is missing.\n");
  }

  // Time integrator for period III, called transient time integrator
  if (rp_list_.isSublist("transient time integrator")) {
    Teuchos::ParameterList& trs_list = rp_list_.sublist("transient time integrator");

    string ti_method_name = trs_list.get<string>("time integration method", "none");
    ProcessStringTimeIntegration(ti_method_name, &ti_specs_trs_.ti_method);
    ProcessSublistTimeIntegration(trs_list, ti_method_name, ti_specs_trs_);
    ti_specs_trs_.ti_method_name = "transient time integrator";

    std::string linear_solver_name = FindStringLinearSolver(trs_list, solver_list_);
    LinearSolver_Specs& ls_specs = ti_specs_trs_.ls_specs;
    ProcessStringLinearSolver(linear_solver_name, &ls_specs.max_itrs, &ls_specs.convergence_tol,
                              &ti_specs_trs_.preconditioner_name);

    ProcessStringPreconditioner(ti_specs_trs_.preconditioner_name, &ti_specs_trs_.preconditioner_method);
    ProcessStringErrorOptions(trs_list, &ti_specs_trs_.error_control_options);

  } else if (verbosity >= FLOW_VERBOSITY_LOW) {
    printf("Flow PK: missing sublist \"transient time integrator\".\n");
  }

  // allowing developer to use non-standard simulation modes
  if (! rp_list_.isParameter("developer access granted")) AnalysisTI_Specs();

  // experimental solver
  string experimental_solver_name = rp_list_.get<string>("experimental solver", "nka");
  ProcessStringExperimentalSolver(experimental_solver_name, &experimental_solver_);

  // optional debug output
  CalculateWRMcurves(rp_list_);
}


/* ****************************************************************
* Process string for error control options
**************************************************************** */
void Richards_PK::ProcessStringErrorOptions(Teuchos::ParameterList& list, int* control)
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
        msg << "Richards Problem: unknown error control option has been specified.";
        Exceptions::amanzi_throw(msg);
      }
    }
  }
}


/* ****************************************************************
* Process string for the relative permeability
**************************************************************** */
void Richards_PK::ProcessStringRelativePermeability(const std::string name, int* method)
{
  Errors::Message msg;
  if (name == "upwind with gravity") {
    *method = AmanziFlow::FLOW_RELATIVE_PERM_UPWIND_GRAVITY;
  } else if (name == "cell centered") {
    *method = AmanziFlow::FLOW_RELATIVE_PERM_CENTERED;
  } else if (name == "upwind with Darcy flux") {
    *method = AmanziFlow::FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX;
  } else if (name == "arithmetic mean") {
    *method = AmanziFlow::FLOW_RELATIVE_PERM_ARITHMETIC_MEAN;
  } else if (name == "upwind experimental") {
    *method = AmanziFlow::FLOW_RELATIVE_PERM_EXPERIMENTAL;
  } else {
    msg << "Richards Problem: unknown relative permeability method has been specified.";
    Exceptions::amanzi_throw(msg);
  }
}


/* ****************************************************************
* Verify string for the relative permeability model.
**************************************************************** */
void Richards_PK::VerifyStringMualemBurdine(const std::string name)
{
  Errors::Message msg;
  if (name != "Mualem" && name != "Burdine") {
    msg << "Flow PK: supported relative permeability models are Mualem and Burdine.";
    Exceptions::amanzi_throw(msg);
  }
}


/* ****************************************************************
* Verify string for the relative permeability model.
**************************************************************** */
void Richards_PK::VerifyWRMparameters(double m, double alpha, double sr, double pc0)
{
  Errors::Message msg;
  if (m < 0.0 || alpha < 0.0 || sr < 0.0 || pc0 < 0.0) {
    msg << "Flow PK: Negative parameter in one of the water retention models.";
    Exceptions::amanzi_throw(msg);    
  }
  if (sr > 1.0) {
    msg << "Flow PK: residual saturation is greater than 1.";
    Exceptions::amanzi_throw(msg);    
  }
}


/* ****************************************************************
* Process string for the time integration method.
**************************************************************** */
void Richards_PK::ProcessStringTimeIntegration(const std::string name, int* method)
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
    msg << "Flow PK: unknown time integration method has been specified.";
    Exceptions::amanzi_throw(msg);
  }
}


/* ****************************************************************
* Process string for the linear solver.
**************************************************************** */
void Richards_PK::ProcessStringLinearSolver(
    const std::string name, int* max_itrs, double* convergence_tol, std::string* preconditioner_name)
{
  Errors::Message msg;

  if (! solver_list_.isSublist(name)) {
    msg << "Flow PK: linear solver does not exist for a time integrator.";
    Exceptions::amanzi_throw(msg);
  }

  Teuchos::ParameterList& tmp_list = solver_list_.sublist(name);
  *max_itrs = tmp_list.get<int>("maximum number of iterations", 100);
  *convergence_tol = tmp_list.get<double>("error tolerance", 1e-10);

  *preconditioner_name = FindStringPreconditioner(tmp_list);
}


/* ****************************************************************
* Process string for the relative permeability
**************************************************************** */
void Richards_PK::ProcessStringExperimentalSolver(const std::string name, int* method)
{
  if (name == "newton") {
    *method = AmanziFlow::FLOW_SOLVER_NEWTON;
  } else if (name == "picard-newton") {
    *method = AmanziFlow::FLOW_SOLVER_PICARD_NEWTON;
  } else {
    *method = AmanziFlow::FLOW_SOLVER_NKA;
  }
}


/* ****************************************************************
* Find string for the preconditoner.
**************************************************************** */
std::string Richards_PK::FindStringPreconditioner(const Teuchos::ParameterList& list)
{   
  Errors::Message msg;
  std::string name; 

  if (list.isParameter("preconditioner")) {
    name = list.get<string>("preconditioner");
  } else {
    msg << "Flow PK: specified linear solver does not define <preconditioner>.";
    Exceptions::amanzi_throw(msg);
  }

  if (! preconditioner_list_.isSublist(name)) {
    msg << "Flow PK: preconditioner \"" << name.c_str() << "\" does not exist.";
    Exceptions::amanzi_throw(msg);
  }
  return name;
}


/* ****************************************************************
* Find string for the preconditoner.
**************************************************************** */
void Richards_PK::CalculateWRMcurves(Teuchos::ParameterList& list)
{
  if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_MEDIUM) {
    if (list.isParameter("calculate krel-pc curves")) {
      std::printf("Flow PK: saving krel-pc curves in file flow_krel_pc.txt...\n");
      ofstream ofile;
      ofile.open("flow_krel_pc.txt");

      std::vector<double> spe;
      spe = list.get<Teuchos::Array<double> >("calculate krel-pc curves").toVector();

      for (double pc = spe[0]; pc < spe[2]; pc += spe[1]) {
        ofile << pc << " ";
        for (int mb = 0; mb < WRM.size(); mb++) {
          double krel = WRM[mb]->k_relative(pc);
          ofile << krel << " ";
        }
        ofile << endl;
      }
      ofile.close();
    }

    if (list.isParameter("calculate krel-sat curves")) {
      std::printf("Flow PK: saving krel-sat curves in file flow_krel_sat.txt...\n");
      ofstream ofile;
      ofile.open("flow_krel_sat.txt");

      std::vector<double> spe;
      spe = list.get<Teuchos::Array<double> >("calculate krel-sat curves").toVector();

      for (double s = spe[0]; s < spe[2]; s += spe[1]) {
        ofile << s << " ";
        for (int mb = 0; mb < WRM.size(); mb++) {
          double ss = std::max<double>(s, WRM[mb]->residualSaturation());
          double pc = WRM[mb]->capillaryPressure(ss);
          double krel = WRM[mb]->k_relative(pc);
          ofile << krel << " ";
        }
        ofile << endl;
      }
      ofile.close();
    }
  }
}


/* ****************************************************************
* Analyzin time integration specs for logical consistency
**************************************************************** */
void Richards_PK::AnalysisTI_Specs()
{
  Errors::Message msg;
  if (ti_specs_igs_.initialize_with_darcy) {
    if (ti_specs_sss_.initialize_with_darcy || ti_specs_trs_.initialize_with_darcy) { 
      msg << "Flow PK: cannot re-initialize pressure without developer password.";
      Exceptions::amanzi_throw(msg);
    }
  }

  if (ti_specs_igs_.dT_method == FLOW_DT_ADAPTIVE || 
      ti_specs_sss_.dT_method == FLOW_DT_ADAPTIVE) {
    msg << "Flow PK: adaptive time stepping is allowed only for transient TI phase.";
    Exceptions::amanzi_throw(msg);
  }
}


/* ****************************************************************
* Printing information about Flow status 
**************************************************************** */
void Richards_PK::PrintStatistics() const
{
  if (! MyPID && verbosity > 0) {
    cout << "Flow PK:" << endl;
    cout << "  Verbosity level = " << verbosity << endl;
    cout << "  Enable internal tests = " << (internal_tests ? "yes" : "no")  << endl;
    cout << "  Upwind = " << ((Krel_method == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) ? "gravity" : "other") << endl;
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi



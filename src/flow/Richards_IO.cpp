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
#include <string>

#include "WRM_vanGenuchten.hpp"
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
  if (!rp_list.isSublist("VerboseObject")) {
    Teuchos::ParameterList verbosity_list;
    verbosity_list.set<std::string>("Verbosity Level", "none");
    rp_list.set("VerboseObject", verbosity_list);
  }

  // extract verbosity level
  Teuchos::ParameterList verbosity_list = rp_list.get<Teuchos::ParameterList>("VerboseObject");
  std::string verbosity_name = verbosity_list.get<std::string>("Verbosity Level");
  if (verbosity_name == "none") {
    verbosity = FLOW_VERBOSITY_NONE;
  } else if (verbosity_name == "low") {
    verbosity = FLOW_VERBOSITY_LOW;
  } else if (verbosity_name == "medium") {
    verbosity = FLOW_VERBOSITY_MEDIUM;
  } else if (verbosity_name == "high") {
    verbosity = FLOW_VERBOSITY_HIGH;
  } else if (verbosity_name == "extreme") {
    verbosity = FLOW_VERBOSITY_EXTREME;
  }

  Teuchos::ParameterList preconditioner_list;
  preconditioner_list = rp_list.get<Teuchos::ParameterList>("Diffusion Preconditioner");

  // Relative permeability method
  std::string method_name = rp_list.get<string>("Relative permeability method", "Upwind with gravity");
  if (method_name == "Upwind with gravity") {
    Krel_method = AmanziFlow::FLOW_RELATIVE_PERM_UPWIND_GRAVITY;
  } else if (method_name == "Cell centered") {
    Krel_method = AmanziFlow::FLOW_RELATIVE_PERM_CENTERED;
  } else if (method_name == "Upwind with Darcy flux") {
    Krel_method = AmanziFlow::FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX;
  } else if (method_name == "Arithmetic mean") {
    Krel_method = AmanziFlow::FLOW_RELATIVE_PERM_ARITHMETIC_MEAN;
  }

  // Miscaleneous
  if (rp_list.isParameter("atmospheric pressure")) {
    atm_pressure = rp_list.get<double>("atmospheric pressure");
  } else {
    msg << "Richards Problem: no <atmospheric pressure> entry.";
    Exceptions::amanzi_throw(msg);
  }

  // Create the BC objects.
  Teuchos::RCP<Teuchos::ParameterList>
      bc_list = Teuchos::rcp(new Teuchos::ParameterList(rp_list.sublist("boundary conditions", true)));
  FlowBCFactory bc_factory(mesh_, bc_list);

  bc_pressure = bc_factory.createPressure();
  bc_head = bc_factory.createStaticHead(atm_pressure, rho, gravity_);
  bc_flux = bc_factory.createMassFlux();
  bc_seepage = bc_factory.createSeepageFace();

  validate_boundary_conditions(bc_pressure, bc_head, bc_flux);

  double T_physical = FS->get_time();  // set-up internal clock
  T_internal = (standalone_mode) ? T_internal : T_physical;

  double time = T_internal;
  bc_pressure->Compute(time);
  bc_head->Compute(time);
  bc_flux->Compute(time);
  bc_seepage->Compute(time);

  // Create water retention models
  if (!rp_list.isSublist("Water retention models")) {
    msg << "There is no Water retention models list";
    Exceptions::amanzi_throw(msg);
  }
  Teuchos::ParameterList& vG_list = rp_list.sublist("Water retention models");

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

        double vG_m = wrm_list.get<double>("van Genuchten m");
        double vG_alpha = wrm_list.get<double>("van Genuchten alpha");
        double vG_sr = wrm_list.get<double>("van Genuchten residual saturation");
        double vG_pc0 = wrm_list.get<double>("regularization interval", 0.0);

        WRM[iblock] = Teuchos::rcp(new WRM_vanGenuchten(region, vG_m, vG_alpha, vG_sr, vG_pc0));

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

  // Time integrator for period I, temporary called steady state time integrator
  if (rp_list.isSublist("steady state time integrator")) {
    Teuchos::ParameterList& sss_list = rp_list.sublist("steady state time integrator");

    string ti_method_name = sss_list.get<string>("method", "Picard");
    if (ti_method_name == "Picard") {
      ti_method_sss = AmanziFlow::FLOW_TIME_INTEGRATION_PICARD;
    } else if (ti_method_name == "backward Euler") {
      ti_method_sss = AmanziFlow::FLOW_TIME_INTEGRATION_BACKWARD_EULER;
    } else if (ti_method_name == "BDF1") {
      ti_method_sss = AmanziFlow::FLOW_TIME_INTEGRATION_BDF1;
    } else if (ti_method_name == "BDF2") {
      ti_method_sss = AmanziFlow::FLOW_TIME_INTEGRATION_BDF2;
    } else {
      msg << "Richards Problem: steady state defines an unknown time integration method.";
      Exceptions::amanzi_throw(msg);
    }

    initialize_with_darcy = (sss_list.get<std::string>("initialize with darcy", "no") == "yes");

    if (sss_list.isSublist("error control")) {
      Teuchos::ParameterList& err_list = sss_list.sublist("error control");
      absolute_tol_sss = err_list.get<double>("absolute error tolerance", 1.0);
      relative_tol_sss = err_list.get<double>("relative error tolerance", 1e-5);
      convergence_tol_sss = err_list.get<double>("convergence tolerance", AmanziFlow::FLOW_TIME_INTEGRATION_TOLERANCE);
      max_itrs_sss = err_list.get<int>("maximal number of iterations", AmanziFlow::FLOW_TIME_INTEGRATION_MAX_ITERATIONS);
    } else {
      msg << "Richards Problem: steady state sublist has no <error control> sublist.";
      Exceptions::amanzi_throw(msg);
    }

    if (sss_list.isSublist("time control")) {
      Teuchos::ParameterList& time_list = sss_list.sublist("time control");
      T0_sss = time_list.get<double>("start time", -1e+9);
      T1_sss = time_list.get<double>("end time", 0.0);
      dT0_sss = time_list.get<double>("initial time step", AmanziFlow::FLOW_INITIAL_DT);
      dTmax_sss = time_list.get<double>("maximal time step", dT0_sss);
    } else {
      msg << "Richards Problem: steady state time integrator has no <time control> sublist.";
      Exceptions::amanzi_throw(msg);
    }

    if (sss_list.isSublist("linear solver")) {
      Teuchos::ParameterList& solver_list = sss_list.sublist("linear solver");
      max_itrs = solver_list.get<int>("maximal number of iterations", 100);
      convergence_tol = solver_list.get<double>("error tolerance", 1e-12);
    } else {
      msg << "Richards Problem: steady state time integrator has no <linear solver> sublist.";
      Exceptions::amanzi_throw(msg);
    }
  } else if (verbosity >= FLOW_VERBOSITY_LOW) {
    printf("Richards Problem: there is no sublist for steady-state calculations.\n");
  }

  // Time integrator for period II, called transient time integrator
  if (rp_list.isSublist("transient time integrator")) {
    Teuchos::ParameterList& trs_list = rp_list.sublist("transient time integrator");

    string ti_method_name = trs_list.get<string>("method", "BDF2");
    if (ti_method_name == "backward Euler") {
      ti_method_trs = AmanziFlow::FLOW_TIME_INTEGRATION_BACKWARD_EULER;
    } else if (ti_method_name == "BDF1") {
      ti_method_trs = AmanziFlow::FLOW_TIME_INTEGRATION_BDF1;
    } else if (ti_method_name == "BDF2") {
      ti_method_trs = AmanziFlow::FLOW_TIME_INTEGRATION_BDF2;
    } else if (ti_method_name == "Picard") {
      ti_method_trs = AmanziFlow::FLOW_TIME_INTEGRATION_PICARD;
    } else {
      msg << "Richards Problem: transient sublist defines an unknown time integration method.";
      Exceptions::amanzi_throw(msg);
    }

    if (trs_list.isSublist("error control")) {
      Teuchos::ParameterList& err_list = trs_list.sublist("error control");
      absolute_tol_trs = err_list.get<double>("absolute error tolerance", 1.0);
      relative_tol_trs = err_list.get<double>("relative error tolerance", 1e-5);
      convergence_tol_trs = err_list.get<double>("convergence tolerance", AmanziFlow::FLOW_TIME_INTEGRATION_TOLERANCE);
      max_itrs_trs = err_list.get<int>("maximal number of iterations", AmanziFlow::FLOW_TIME_INTEGRATION_MAX_ITERATIONS);
    } else {
      msg << "Richards Problem: transient sublist has no <error control> sublist.";
      Exceptions::amanzi_throw(msg);
    }

    if (trs_list.isSublist("time control")) {
      Teuchos::ParameterList& time_list = trs_list.sublist("time control");
      T0_trs = time_list.get<double>("start time", 0.0);
      T1_trs = time_list.get<double>("end time", 0.0);
      dT0_trs = time_list.get<double>("initial time step", AmanziFlow::FLOW_INITIAL_DT);
      dTmax_trs = time_list.get<double>("maximal time step", dT0_trs);
    } else {
      msg << "Richards Problem: transient sublist has no <time control> sublist.";
      Exceptions::amanzi_throw(msg);
    }
  } else if (verbosity >= FLOW_VERBOSITY_LOW) {
    printf("Warning: Richards Problem has no sublist <transient time integration>.\n");
  }
}


/* ****************************************************************
* Printing information about Flow status 
**************************************************************** */
void Richards_PK::PrintStatistics() const
{
  if (!MyPID && verbosity > 0) {
    cout << "Flow PK:" << endl;
    cout << "  Execution mode = " << (standalone_mode ? "standalone" : "MPC") << endl;
    cout << "  Verbosity level = " << verbosity << endl;
    cout << "  Enable internal tests = " << (internal_tests ? "yes" : "no")  << endl;
    cout << "  Upwind = " << ((Krel_method == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) ? "gravity" : "other") << endl;
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi



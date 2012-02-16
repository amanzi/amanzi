/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (nnc@lanl.gov), 
         Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "WRM_vanGenuchten.hpp"
#include "WRM_fake.hpp"
#include "Richards_PK.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
void Richards_PK::processParameterList()
{
  Teuchos::ParameterList preconditioner_list;
  preconditioner_list = rp_list.get<Teuchos::ParameterList>("Diffusion Preconditioner");

  // Relative permeability method
  string method_name = rp_list.get<string>("Relative permeability method", "Upwind with gravity");
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
  atm_pressure = rp_list.get<double>("Atmospheric pressure");

  // Create the BC objects.
  Teuchos::RCP<Teuchos::ParameterList> 
      bc_list = Teuchos::rcp(new Teuchos::ParameterList(rp_list.sublist("boundary conditions", true)));
  FlowBCFactory bc_factory(mesh_, bc_list);

  bc_pressure = bc_factory.createPressure();
  bc_head = bc_factory.createStaticHead(atm_pressure, rho, gravity_);
  bc_flux = bc_factory.createMassFlux();

  validate_boundary_conditions(bc_pressure, bc_head, bc_flux);  

  double time = (standalone_mode) ? T_internal : T_physical;
  bc_pressure->Compute(time);
  bc_head->Compute(time);
  bc_flux->Compute(time);

  // Create water retention models
  if (!rp_list.isSublist("Water retention models")) {
    Errors::Message m("There is no Water retention models list");
    Exceptions::amanzi_throw(m);
  }
  Teuchos::ParameterList& vG_list = rp_list.sublist("Water retention models");

  int nblocks = 0;  // Find out how many WRM entries there are.
  for (Teuchos::ParameterList::ConstIterator i = vG_list.begin(); i != vG_list.end(); i++) {
    if (vG_list.isSublist(vG_list.name(i))) {
      nblocks++;
    } else {
      Errors::Message msg("Water retention models sublist contains an entry that is not a sublist.");
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
	      
        WRM[iblock] = Teuchos::rcp(new WRM_vanGenuchten(region, vG_m, vG_alpha, vG_sr));
      } 
      else if (wrm_list.get<string>("Water retention model") == "fake") {
        std::string region = wrm_list.get<std::string>("Region");  // associated mesh block
        WRM[iblock] = Teuchos::rcp(new WRM_fake(region));
      }
      iblock++;
    }
  }

  // Steady state solution
  if (rp_list.isSublist("Steady state solution")) {
    Teuchos::ParameterList& sss_list = rp_list.sublist("Steady state solution");

    method_name = sss_list.get<string>("method", "Picard");
    if (method_name == "Picard") {
      method_sss = AmanziFlow::FLOW_STEADY_STATE_PICARD;
    } else if (method_name == "backward Euler") {
      method_sss = AmanziFlow::FLOW_STEADY_STATE_BACKWARD_EULER;
    } else if (method_name == "BDF2") {
      method_sss = AmanziFlow::FLOW_STEADY_STATE_BDF2; 
    }

    Teuchos::ParameterList& err_list = sss_list.sublist("Error control");
    absolute_tol_sss = err_list.get<double>("Absolute error tolerance", 1.0);
    relative_tol_sss = err_list.get<double>("Relative error tolerance", 1e-5); 
    convergence_tol_sss = err_list.get<double>("Convergence tolerance", AmanziFlow::FLOW_STEADY_STATE_TOLERANCE);
    max_itrs_sss = err_list.get<int>("Maximal number of iterations", AmanziFlow::FLOW_STEADY_STATE_MAX_ITERATIONS);

    Teuchos::ParameterList& time_list = sss_list.sublist("Time control");
    T0_sss = time_list.get<double>("Start pseudo time", -1e+9);
    T1_sss = time_list.get<double>("End pseudo time", 0.0);
    dT0_sss = time_list.get<double>("Initial time step", AmanziFlow::FLOW_STEADY_STATE_INITIAL_DT);
    dTmax_sss = time_list.get<double>("Maximal time step", dT0_sss);

    Teuchos::ParameterList& solver_list = sss_list.sublist("Linear solvers");
    max_itrs = solver_list.get<int>("Maximal number of iterations", 100);
    convergence_tol = solver_list.get<double>("Error tolerance", 1e-12);
  } else if (verbosity >= FLOW_VERBOSITY_LOW) {
    printf("FLOW: no sublist for steady-state calculation was provided\n");
  }

  // Transient solution
  if (rp_list.isSublist("Transient solution")) {
    Teuchos::ParameterList& trs_list = rp_list.sublist("Transient solution");

    method_name = trs_list.get<string>("method", "BDF2");
    if (method_name == "backward Euler") {
      method_trs = AmanziFlow::FLOW_STEADY_STATE_BACKWARD_EULER;
    } else if (method_name == "BDF2") {
      method_trs = AmanziFlow::FLOW_STEADY_STATE_BDF2; 
    }

    Teuchos::ParameterList& err_list = trs_list.sublist("Error control");
    absolute_tol_trs = err_list.get<double>("Absolute error tolerance", 1.0);
    relative_tol_trs = err_list.get<double>("Relative error tolerance", 1e-5); 
    convergence_tol_trs = err_list.get<double>("Convergence tolerance", AmanziFlow::FLOW_STEADY_STATE_TOLERANCE);
    max_itrs_trs = err_list.get<int>("Maximal number of iterations", AmanziFlow::FLOW_STEADY_STATE_MAX_ITERATIONS);

    Teuchos::ParameterList& time_list = trs_list.sublist("Time control");
    T0_trs = time_list.get<double>("Start time", 0.0);
    T1_trs = time_list.get<double>("End time", 0.0);
    dT0_trs = time_list.get<double>("Initial time step", AmanziFlow::FLOW_STEADY_STATE_INITIAL_DT);
    dTmax_trs = time_list.get<double>("Maximal time step", dT0_trs);
  } else if (verbosity >= FLOW_VERBOSITY_LOW) {
    printf("FLOW: no sublist for transient calculation was provided\n");
  }
}


/* ****************************************************************
* Printing information about Flow status 
**************************************************************** */
void Richards_PK::print_statistics() const
{
  if (!MyPID && verbosity > 0) {
    cout << "Flow PK:" << endl;
    cout << "    Execution mode = " << (standalone_mode ? "standalone" : "MPC") << endl;
    cout << "    Verbosity level = " << verbosity << endl;
    cout << "    Enable internal tests = " << (internal_tests ? "yes" : "no")  << endl;
    cout << "    Upwind = " << ((Krel_method == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) ? "gravity" : "other") << endl;
  }
}
 

}  // namespace AmanziFlow
}  // namespace Amanzi



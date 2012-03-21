/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include "Flow_BC_Factory.hpp"
#include "Darcy_PK.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
void Darcy_PK::processParameterList()
{
  Teuchos::ParameterList preconditioner_list;
  preconditioner_list = dp_list.get<Teuchos::ParameterList>("Diffusion Preconditioner");

  atm_pressure = dp_list.get<double>("Atmospheric pressure");

  // Create the BC objects.
  Teuchos::RCP<Teuchos::ParameterList> bc_list = Teuchos::rcpFromRef(dp_list.sublist("boundary conditions", true));
  FlowBCFactory bc_factory(mesh_, bc_list);

  bc_pressure = bc_factory.createPressure();
  bc_head = bc_factory.createStaticHead(atm_pressure, rho_, gravity_);
  bc_flux = bc_factory.createMassFlux();
  bc_seepage = bc_factory.createSeepageFace();

  validate_boundary_conditions(bc_pressure, bc_head, bc_flux);  

  double T_physical = FS->get_time();
  T_internal = (standalone_mode) ? T_internal : T_physical;

  double time = T_internal;
  bc_pressure->Compute(time);
  bc_head->Compute(time);
  bc_flux->Compute(time);
  bc_seepage->Compute(time);

  // Steady state solution
  Teuchos::ParameterList& sss_list = dp_list.sublist("Steady state solution");
  Teuchos::ParameterList& solver_list = sss_list.sublist("Linear solvers");

  max_itrs_sss = solver_list.get<int>("Maximal number of iterations", 100);
  convergence_tol_sss = solver_list.get<double>("Error tolerance", 1e-12);

  string method_name = sss_list.get<string>("Discretization method hint", "monotone");
  if (method_name == "monotone") {
    mfd3d_method = FLOW_MFD3D_HEXAHEDRA_MONOTONE;
  } else if (method_name == "none") {
    mfd3d_method = FLOW_MFD3D_POLYHEDRA;
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi



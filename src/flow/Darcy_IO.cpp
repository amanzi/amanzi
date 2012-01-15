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
  bc_head = bc_factory.createStaticHead(atm_pressure, rho, -gravity[dim - 1]);  // We need |g| here.
  bc_flux = bc_factory.createMassFlux();

  validate_boundary_conditions(bc_pressure, bc_head, bc_flux);  

  double time = (standalone_mode) ? T_internal : T_physical;
  bc_pressure->Compute(time);
  bc_head->Compute(time);
  bc_flux->Compute(time);

  // Steady state solution
  Teuchos::ParameterList& sss_list = dp_list.sublist("Steady state solution");

  max_itrs_sss = preconditioner_list.get<int>("Linear solver maximal iterations", 100);
  convergence_tol_sss = preconditioner_list.get<double>("Linear solver tolerance", 1e-12);
}

}  // namespace AmanziFlow
}  // namespace Amanzi



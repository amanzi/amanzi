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
void Darcy_PK::process_parameter_list()
{
  Teuchos::ParameterList preconditioner_list;
  preconditioner_list = dp_list.get<Teuchos::ParameterList>("Diffusion Preconditioner");

  max_itrs = preconditioner_list.get<int>("Max Iterations");
  err_tol = preconditioner_list.get<double>("Error Tolerance");

  // Create the BC objects.
  Teuchos::RCP<Teuchos::ParameterList> bc_list = Teuchos::rcpFromRef(dp_list.sublist("boundary conditions", true));
  FlowBCFactory bc_factory(mesh_, bc_list);

  bc_pressure = bc_factory.CreatePressure();
  bc_head = bc_factory.CreateStaticHead(0.0, rho, gravity[dim - 1]);
  bc_flux = bc_factory.CreateMassFlux();

  validate_boundary_conditions(bc_pressure, bc_head, bc_flux);  

  double time = (standalone_mode) ? T_internal : T_physical;
  bc_pressure->Compute(time);
  bc_head->Compute(time);
  bc_flux->Compute(time);
}

}  // namespace AmanziFlow
}  // namespace Amanzi



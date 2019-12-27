/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy

  Process kernel for coupling Transport and Chemistry PKs in the 
  matrix and fracture network.
*/

#include "Alquimia_PK.hh"
#include "ReactiveTransportMatrixFracture_PK.hh"
#include "Transport_PK.hh"
#include "TransportMatrixFracture_PK.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Standard constructor
// -----------------------------------------------------------------------------
ReactiveTransportMatrixFracture_PK::ReactiveTransportMatrixFracture_PK(
    Teuchos::ParameterList& pk_tree,
    const Teuchos::RCP<Teuchos::ParameterList>& global_list,
    const Teuchos::RCP<State>& S,
    const Teuchos::RCP<TreeVector>& soln) :
    PK_MPCAdditive<PK>(pk_tree, global_list, S, soln)
{
  coupled_chemistry_pk_ = Teuchos::rcp_dynamic_cast<PK_MPCWeak>(sub_pks_[0]);
  coupled_transport_pk_ = Teuchos::rcp_dynamic_cast<PK_MPCWeak>(sub_pks_[1]);
}


// -----------------------------------------------------------------------------
// Setup delegates work to base PK
// -----------------------------------------------------------------------------
void ReactiveTransportMatrixFracture_PK::Setup(const Teuchos::Ptr<State>& S)
{
  Amanzi::PK_MPCAdditive<PK>::Setup(S);
  
  // communicate chemistry engine to transport.
#ifdef ALQUIMIA_ENABLED
  auto ic = coupled_chemistry_pk_->begin();
  if ((*ic)->name() == "chemistry alquimia") {
    for (auto it = coupled_transport_pk_->begin(); it != coupled_transport_pk_->end(); ++it, ++ic) {
      auto it1 = Teuchos::rcp_dynamic_cast<Transport::Transport_PK>(*it);
      auto ic1 = Teuchos::rcp_dynamic_cast<AmanziChemistry::Alquimia_PK>(*ic);
      it1->SetupAlquimia(ic1, ic1->chem_engine());
    }
  }
#endif
}


// -----------------------------------------------------------------------------
// Calculate minimum of sub PKs timestep sizes
// -----------------------------------------------------------------------------
double ReactiveTransportMatrixFracture_PK::get_dt()
{
  double dTtran = coupled_transport_pk_->get_dt();
  double dTchem = coupled_chemistry_pk_->get_dt();

  if (!chem_step_succeeded_ && (dTchem / dTtran > 0.99)) {
    dTchem *= 0.5;
  } 

  if (dTtran > dTchem) dTtran = dTchem; 
  
  return dTchem;
}


// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
void ReactiveTransportMatrixFracture_PK::set_dt(double dt)
{
  coupled_chemistry_pk_->set_dt(dt);
  coupled_transport_pk_->set_dt(dt);
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible
// -----------------------------------------------------------------------------
bool ReactiveTransportMatrixFracture_PK::AdvanceStep(
    double t_old, double t_new, bool reinit)
{
  bool fail = false;
  chem_step_succeeded_ = false;


  /*
  // First we do a transport step.
  bool pk_fail = tranport_pk_->AdvanceStep(t_old, t_new, reinit);

  // Right now transport step is always succeeded.
  if (!pk_fail) {
    *total_component_concentration_stor = *tranport_pk_->total_component_concentration()->ViewComponent("cell", true);
  } else {
    Errors::Message message("MPC: Transport PK returned an unexpected error.");
    Exceptions::amanzi_throw(message);
  }

  // Second, we do a chemistry step.
  try {
    chemistry_pk_->set_aqueous_components(total_component_concentration_stor);

    pk_fail = chemistry_pk_->AdvanceStep(t_old, t_new, reinit);
    chem_step_succeeded_ = true;
 
    *S_->GetFieldData("total_component_concentration", "state")
       ->ViewComponent("cell", true) = *chemistry_pk_->aqueous_components();
  }
  catch (const Errors::Message& chem_error) {
    fail = true;
  }
  */


  return fail;
};

}  // namespace Amanzi


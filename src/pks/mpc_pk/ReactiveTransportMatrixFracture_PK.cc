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
  coupled_chemistry_pk_ = Teuchos::rcp_dynamic_cast<ChemistryMatrixFracture_PK>(sub_pks_[0]);
  coupled_transport_pk_ = Teuchos::rcp_dynamic_cast<PK_MPC<PK> >(sub_pks_[1]);
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
  double dTchem = sub_pks_[0]->get_dt();
  double dTtran = sub_pks_[1]->get_dt();

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
  sub_pks_[0]->set_dt(dt);
  sub_pks_[1]->set_dt(dt);
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible
// -----------------------------------------------------------------------------
bool ReactiveTransportMatrixFracture_PK::AdvanceStep(
    double t_old, double t_new, bool reinit)
{
  bool fail = sub_pks_[1]->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;

  // save copy of fields (FIXME)
  Teuchos::RCP<Epetra_MultiVector> tcc_m_copy, tcc_f_copy;
  tcc_m_copy = Teuchos::rcp(new Epetra_MultiVector(
      *S_->GetFieldData("total_component_concentration")->ViewComponent("cell", true)));
  tcc_f_copy = Teuchos::rcp(new Epetra_MultiVector(
      *S_->GetFieldData("fracture-total_component_concentration")->ViewComponent("cell", true)));

  try {
    std::vector<Teuchos::RCP<AmanziChemistry::Chemistry_PK> > subpks;
    for (auto ic = coupled_chemistry_pk_->begin(); ic != coupled_chemistry_pk_->end(); ++ic) { 
      auto ic1 = Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(*ic);
      subpks.push_back(ic1);
    }
    
    subpks[0]->set_aqueous_components(tcc_m_copy);
    subpks[1]->set_aqueous_components(tcc_f_copy);

    fail = coupled_chemistry_pk_->AdvanceStep(t_old, t_new, reinit);
 
    *S_->GetFieldData("total_component_concentration", "state")
      ->ViewComponent("cell", true) = *subpks[0]->aqueous_components();

    *S_->GetFieldData("fracture-total_component_concentration", "state")
      ->ViewComponent("cell", true) = *subpks[1]->aqueous_components();
  }
  catch (const Errors::Message& chem_error) {
    fail = true;
  }

  return fail;
};

}  // namespace Amanzi


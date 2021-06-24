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

#include "Chemistry_PK.hh"

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
}


// -----------------------------------------------------------------------------
// Setup delegates work to base PK
// -----------------------------------------------------------------------------
void ReactiveTransportMatrixFracture_PK::Setup(const Teuchos::Ptr<State>& S)
{
  mesh_domain_ = S->GetMesh();
  mesh_fracture_ = S->GetMesh("fracture");

  // darcy fluxes use non-uniform distribution of DOFs
  // -- darcy flux for matrix
  if (!S->HasField("darcy_flux")) {
    auto cvs = Operators::CreateFracturedMatrixCVS(mesh_domain_, mesh_fracture_);
    auto mmap = cvs->Map("face", false);
    auto gmap = cvs->Map("face", true);
    S->RequireField("darcy_flux", "transport")->SetMesh(mesh_domain_)->SetGhosted(true) 
      ->SetComponent("face", AmanziMesh::FACE, mmap, gmap, 1);
  }

  // -- darcy flux for fracture
  if (!S->HasField("fracture-darcy_flux")) {
    auto cvs = Operators::CreateNonManifoldCVS(mesh_fracture_);
    *S->RequireField("fracture-darcy_flux", "transport")->SetMesh(mesh_fracture_)->SetGhosted(true) = *cvs;
  }

  tcc_matrix_key_="total_component_concentration";
  tcc_fracture_key_="fracture-total_component_concentration";
  
  // evaluators in fracture
  Teuchos::ParameterList& elist = S->FEList();
  Teuchos::ParameterList& ilist = S->ICList();

  double rho = ilist.sublist("const_fluid_density").get<double>("value");
  elist.sublist("fracture-mass_density_liquid").sublist("function").sublist("All")
      .set<std::string>("region", "All")
      .set<std::string>("component", "cell")
      .sublist("function").sublist("function-constant")
      .set<double>("value", rho);
  elist.sublist("fracture-mass_density_liquid")
        .set<std::string>("field evaluator type", "independent variable");

  // communicate chemistry engine to transport.
  auto ic = coupled_chemistry_pk_->begin();

  auto mpc_implicit = Teuchos::rcp_dynamic_cast<PK_MPC<PK_BDF> >(sub_pks_[1]);
  if (mpc_implicit.get() != NULL) {
    for (auto it = mpc_implicit->begin(); it != mpc_implicit->end(); ++it, ++ic) {
      auto it1 = Teuchos::rcp_dynamic_cast<Transport::Transport_PK>(*it);
      auto ic1 = Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(*ic);
      it1->SetupChemistry(ic1);
    }
  } else {
    auto mpc_explicit = Teuchos::rcp_dynamic_cast<PK_MPC<PK> >(sub_pks_[1]);
    for (auto it = mpc_explicit->begin(); it != mpc_explicit->end(); ++it, ++ic) {
      auto it1 = Teuchos::rcp_dynamic_cast<Transport::Transport_PK>(*it);
      auto ic1 = Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(*ic);
      it1->SetupChemistry(ic1);
    }
  }

  Amanzi::PK_MPCAdditive<PK>::Setup(S);
}

  
// -----------------------------------------------------------------------------
// Initialization of copies requires fileds to exists
// -----------------------------------------------------------------------------
void ReactiveTransportMatrixFracture_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  S->RequireFieldCopy(tcc_matrix_key_, "matrix_copy", "state");
  S->RequireFieldCopy(tcc_fracture_key_, "fracture_copy", "state");

  Amanzi::PK_MPCAdditive<PK>::Initialize(S);
}
  

// -----------------------------------------------------------------------------
// Calculate minimum of sub PKs timestep sizes
// -----------------------------------------------------------------------------
double ReactiveTransportMatrixFracture_PK::get_dt()
{
  double dTchem = sub_pks_[0]->get_dt();
  double dTtran = sub_pks_[1]->get_dt();

  if (dTchem / dTtran > 0.99) {
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
  S_->CopyField(tcc_matrix_key_, "matrix_copy", "state");
  S_->CopyField(tcc_fracture_key_, "fracture_copy", "state");  
  
  try {
    std::vector<Teuchos::RCP<AmanziChemistry::Chemistry_PK> > subpks;
    for (auto ic = coupled_chemistry_pk_->begin(); ic != coupled_chemistry_pk_->end(); ++ic) { 
      auto ic1 = Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(*ic);
      subpks.push_back(ic1);
    }
    
    // tell chemistry PKs to work with copies
    auto tcc_m_copy = S_->GetFieldCopyData(tcc_matrix_key_, "matrix_copy", "state")->ViewComponent("cell", true);
    auto tcc_f_copy = S_->GetFieldCopyData(tcc_fracture_key_, "fracture_copy", "state")->ViewComponent("cell", true);  

    subpks[0]->set_aqueous_components(tcc_m_copy);
    subpks[1]->set_aqueous_components(tcc_f_copy);

    fail = coupled_chemistry_pk_->AdvanceStep(t_old, t_new, reinit);
 
    *S_->GetFieldData(tcc_matrix_key_, "state")
      ->ViewComponent("cell", true) = *subpks[0]->aqueous_components();

    *S_->GetFieldData(tcc_fracture_key_, "state")
      ->ViewComponent("cell", true) = *subpks[1]->aqueous_components();
  }
  catch (const Errors::Message& chem_error) {
    fail = true;
  }

  return fail;
};

}  // namespace Amanzi


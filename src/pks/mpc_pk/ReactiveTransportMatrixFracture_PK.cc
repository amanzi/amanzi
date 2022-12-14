/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
*/

/*
  This is the mpc_pk component of the Amanzi code.

  Process kernel for coupling Transport and Chemistry PKs in the
  matrix and fracture network.
*/

#include "Chemistry_PK.hh"

#include "ReactiveTransportMatrixFracture_PK.hh"
#include "Transport_PK.hh"
#include "TransportMatrixFracture_PK.hh"
#include "TransportMatrixFractureImplicit_PK.hh"

namespace Amanzi {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

// -----------------------------------------------------------------------------
// Standard constructor
// -----------------------------------------------------------------------------
ReactiveTransportMatrixFracture_PK::ReactiveTransportMatrixFracture_PK(
  Teuchos::ParameterList& pk_tree,
  const Teuchos::RCP<Teuchos::ParameterList>& global_list,
  const Teuchos::RCP<State>& S,
  const Teuchos::RCP<TreeVector>& soln)
  : PK_MPCSubcycled(pk_tree, global_list, S, soln)
{
  coupled_chemistry_pk_ = Teuchos::rcp_dynamic_cast<ChemistryMatrixFracture_PK>(sub_pks_[0]);

  subcycling_ = my_list_->get<bool>("subcycle chemistry", true);

  AMANZI_ASSERT(master_ == 1);
}


// -----------------------------------------------------------------------------
// Setup delegates work to base PK
// -----------------------------------------------------------------------------
void
ReactiveTransportMatrixFracture_PK::Setup()
{
  mesh_domain_ = S_->GetMesh();
  mesh_fracture_ = S_->GetMesh("fracture");

  // darcy fluxes use non-uniform distribution of DOFs
  // -- darcy flux for matrix
  if (!S_->HasRecord("volumetric_flow_rate")) {
    auto cvs = Operators::CreateFracturedMatrixCVS(mesh_domain_, mesh_fracture_);
    auto mmap = cvs->Map("face", false);
    auto gmap = cvs->Map("face", true);
    S_->Require<CV_t, CVS_t>("volumetric_flow_rate", Tags::DEFAULT, "transport")
      .SetMesh(mesh_domain_)
      ->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, mmap, gmap, 1);
  }

  // -- darcy flux for fracture
  if (!S_->HasRecord("fracture-volumetric_flow_rate")) {
    auto cvs = Operators::CreateManifoldCVS(mesh_fracture_);
    *S_->Require<CV_t, CVS_t>("fracture-volumetric_flow_rate", Tags::DEFAULT, "transport")
       .SetMesh(mesh_fracture_)
       ->SetGhosted(true) = *cvs;
  }

  tcc_matrix_key_ = "total_component_concentration";
  tcc_fracture_key_ = "fracture-total_component_concentration";

  // evaluators in fracture
  Teuchos::ParameterList& elist = S_->FEList();
  Teuchos::ParameterList& ilist = S_->ICList();

  double rho = ilist.sublist("const_fluid_density").get<double>("value");
  elist.sublist("fracture-mass_density_liquid")
    .sublist("function")
    .sublist("All")
    .set<std::string>("region", "All")
    .set<std::string>("component", "cell")
    .sublist("function")
    .sublist("function-constant")
    .set<double>("value", rho);
  elist.sublist("fracture-mass_density_liquid")
    .set<std::string>("evaluator type", "independent variable");

  // copies
  S_->Require<CV_t, CVS_t>(tcc_matrix_key_, Tags::COPY, "state");
  S_->Require<CV_t, CVS_t>(tcc_fracture_key_, Tags::COPY, "state");

  // communicate chemistry engine to transport.
  auto ic = coupled_chemistry_pk_->begin();

  auto mpc_implicit = Teuchos::rcp_dynamic_cast<PK_MPC<PK_BDF>>(sub_pks_[1]);
  if (mpc_implicit.get() != NULL) {
    for (auto it = mpc_implicit->begin(); it != mpc_implicit->end(); ++it, ++ic) {
      auto it1 = Teuchos::rcp_dynamic_cast<Transport::Transport_PK>(*it);
      auto ic1 = Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(*ic);
      it1->SetupChemistry(ic1);
    }
  } else {
    auto mpc_explicit = Teuchos::rcp_dynamic_cast<PK_MPC<PK>>(sub_pks_[1]);
    for (auto it = mpc_explicit->begin(); it != mpc_explicit->end(); ++it, ++ic) {
      auto it1 = Teuchos::rcp_dynamic_cast<Transport::Transport_PK>(*it);
      auto ic1 = Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(*ic);
      it1->SetupChemistry(ic1);
    }
  }

  Amanzi::PK_MPCSubcycled::Setup();
}


// -----------------------------------------------------------------------------
// Initialization of copies requires fileds to exists
// -----------------------------------------------------------------------------
void
ReactiveTransportMatrixFracture_PK::Initialize()
{
  Amanzi::PK_MPCSubcycled::Initialize();
}


// -----------------------------------------------------------------------------
// Calculate minimum of sub PKs timestep sizes
// -----------------------------------------------------------------------------
double
ReactiveTransportMatrixFracture_PK::get_dt()
{
  double dTchem = sub_pks_[0]->get_dt();
  double dTtran = sub_pks_[1]->get_dt();

  if (subcycling_) return dTtran;

  if (dTchem / dTtran > 0.99) dTchem *= 0.5;
  if (dTtran > dTchem) dTtran = dTchem;

  return dTchem;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void
ReactiveTransportMatrixFracture_PK::set_dt(double dt)
{
  sub_pks_[0]->set_dt(dt);
  sub_pks_[1]->set_dt(dt);
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible
// -----------------------------------------------------------------------------
bool
ReactiveTransportMatrixFracture_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = sub_pks_[1]->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;

  std::vector<Teuchos::RCP<AmanziChemistry::Chemistry_PK>> subpks_chem;
  for (auto ic = coupled_chemistry_pk_->begin(); ic != coupled_chemistry_pk_->end(); ++ic) {
    auto ic1 = Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(*ic);
    subpks_chem.push_back(ic1);
  }

  std::vector<Teuchos::RCP<Transport::Transport_PK>> subpks_tran;
  if (sub_pks_[1]->name() == "coupled transport") {
    auto tpk = Teuchos::rcp_dynamic_cast<TransportMatrixFracture_PK>(sub_pks_[1]);

    for (auto ic = tpk->begin(); ic != tpk->end(); ++ic) {
      auto ic1 = Teuchos::rcp_dynamic_cast<Transport::Transport_PK>(*ic);
      subpks_tran.push_back(ic1);
    }
  } else if (sub_pks_[1]->name() == "coupled transport implicit") {
    auto tpk = Teuchos::rcp_dynamic_cast<TransportMatrixFractureImplicit_PK>(sub_pks_[1]);

    for (auto ic = tpk->begin(); ic != tpk->end(); ++ic) {
      auto ic1 = Teuchos::rcp_dynamic_cast<Transport::Transport_PK>(*ic);
      subpks_tran.push_back(ic1);
    }
  }

  // at this moment, we have old tcc in State and new tcc in a temporary vector
  // tell chemistry PKs to work with the temporary vector
  for (int i = 0; i < 2; ++i) {
    auto tcc_copy = subpks_tran[i]->total_component_concentration()->ViewComponent("cell", true);
    subpks_chem[i]->set_aqueous_components(tcc_copy);
  }

  // advance the slave, subcycling if needed
  double dTchem = sub_pks_[0]->get_dt();
  double dt_next(dTchem), dt_done(0.0);

  try {
    S_->set_intermediate_time(t_old);
    bool done = false;

    while (!done) {
      // do not overstep
      if (t_old + dt_done + dt_next > t_new) { dt_next = t_new - t_old - dt_done; }
      S_->set_intermediate_time(t_old + dt_done + dt_next);

      fail = coupled_chemistry_pk_->AdvanceStep(t_old + dt_done, t_old + dt_done + dt_next, reinit);

      if (fail) {
        dt_next /= 2;
        for (int i = 0; i < 2; ++i) {
          auto tcc_copy =
            subpks_tran[i]->total_component_concentration()->ViewComponent("cell", true);
          subpks_chem[i]->set_aqueous_components(tcc_copy);
        }
      } else {
        dt_done += dt_next;
        dt_next = sub_pks_[0]->get_dt();
      }

      // no state recovery is made, so the only option is to fail
      if (dt_next < min_dt_)
        Exceptions::amanzi_throw("Failure in coupled reactive transport PK: small time step.");

      // check for subcycling condition
      done = std::abs(t_old + dt_done - t_new) / (t_new - t_old) < 0.1 * min_dt_;
    }

    *S_->GetW<CV_t>(tcc_matrix_key_, Tags::DEFAULT, "state").ViewComponent("cell", true) =
      *subpks_chem[0]->aqueous_components();

    *S_->GetW<CV_t>(tcc_fracture_key_, Tags::DEFAULT, "state").ViewComponent("cell", true) =
      *subpks_chem[1]->aqueous_components();
  } catch (const Errors::Message& chem_error) {
    fail = true;
  } catch (...) {
    fail = true;
  }

  return fail;
};

} // namespace Amanzi
